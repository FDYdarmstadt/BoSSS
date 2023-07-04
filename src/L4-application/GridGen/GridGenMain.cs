/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System.Linq;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution;
using System.Collections.Generic;
using System.Collections;

namespace BoSSS.Application.GridGen {

    /// <summary>
    /// A simple app for grid creation, mainly for Cartesian meshes.
    /// Typically, one could create such meshes using a Jupyter notebook.
    /// Instead, this is most useful when really large meshes should be created, e.g. for benchmarking, 
    /// which exceed the capabilities of a workstation running Jupyter.
    /// In such scenarios, this app can be used to create a mesh in parallel.
    /// </summary>
    public class GridGenMain : BoSSS.Solution.Application<GridGenControl> {

        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

            
            _Main(
                args,
                false,
                () => new GridGenMain());
        }



        /// <summary>
        /// 
        /// </summary>
        protected override IGrid CreateOrLoadGrid() {


            Console.Write("Grid instantiation...");
            GridCommons[] grds = new GridCommons[this.Control.GridBlocks.Length];
            for (int iBlock = 0; iBlock < grds.Length; iBlock++) {
                Console.Write($" Block{iBlock+1}of{grds.Length}");
                grds[iBlock] = this.Control.GridBlocks[iBlock].CreateGrid();
                Console.Write(".");
            }
            Console.WriteLine(" done.");

            GridCommons grd;
            if(grds.Length <= 1) {
                grd = grds[0];
            } else {
                Console.Write("Sealing Blocks...");
                grd = GridCommons.MergeLogically(grds);
                grd = GridCommons.Seal(grd);
                Console.WriteLine(" done.");
            }

            if (Control.BoundaryRegions.Count > 0) {
                Console.WriteLine("Setting Boundary Regions...");
                string EdgeTagNameFunc(Vector X) {
                    foreach (var pair in Control.BoundaryRegions) {
                        if (pair.Region == null)
                            return pair.EdgeTagName;

                        if (pair.Region.Contains(X))
                            return pair.EdgeTagName;

                    }

                    return "undefined";
                }


                grd.DefineEdgeTags(EdgeTagNameFunc, Control.EdgeTagNamesToEnsure);

                Console.WriteLine("done.");
            }

            Console.WriteLine("Grid created.");

            grd.Name = Control.GridName;
            grd.Description = Control.Description;
            IGrid _grd = grd;

            Console.WriteLine("Saving to database...");
            this.DatabaseDriver.SaveGridIfUnique(ref _grd, out bool EquivFound, this.GetDatabase());
            if (EquivFound) {
                Console.WriteLine($"Note: Equivalent grid {_grd.ID} already in database; Using the grid from database.");
            } else {
                Console.WriteLine($"Saveg grid: {_grd.ID}");
            }

            bool bkup = ilPSP.Environment.StdoutOnlyOnRank0;
            ilPSP.Environment.StdoutOnlyOnRank0 = false;
            Console.WriteLine($" cells on MPI rank {MPIRank}: {_grd.CellPartitioning.LocalLength}");
            ilPSP.Environment.StdoutOnlyOnRank0 = bkup;

            return _grd;
        }


        public override void Init(BoSSS.Solution.Control.AppControl control) {
            control.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;
            base.Init(control);
        }


        SinglePhaseField m_MPIRank;

        SinglePhaseField[] bndyMarkers;


        protected override void CreateFields() {
            m_MPIRank = new SinglePhaseField(new Basis(this.Grid, 0), "MPIrank");
            this.m_MPIRank.AccConstant(base.MPIRank);
            base.RegisterField(m_MPIRank, IOListOption.Always);


            BndyMarkerSet();
        }


        List<(string EdgeTagName, int Count)> BoundaryStat = new List<(string EdgeTagName, int Count)>();

        void BndyMarkerSet() {

            var grd = this.GridData;
            
            string SanitizeName(string s) {
                char[] ot = s.ToCharArray();
                for (int k = 0; k < ot.Length; k++) {
                    if (char.IsWhiteSpace(ot[k])) {
                        ot[k] = '_';
                    }

                    if (ot[k] == '(')
                        ot[k] = 'L';
                    if (ot[k] == ')')
                        ot[k] = 'R';
                }
                return new string(ot);
            }



            var et2Name = grd.EdgeTagNames;
            Console.WriteLine($"Grid containing {et2Name.Count} EdgeTag names: ");
            int i = 0;
            foreach (var t in et2Name) { // loop over all different edge tag names...
                string name = t.Value;
                byte tag2color = t.Key;

                string sname = SanitizeName(name);

                if (name.Equals(sname)) {
                    Console.WriteLine($"   {i}: {name} -- tag = {tag2color}");
                } else {
                    Console.WriteLine($"   {i}: {name} -- tag = {tag2color}   (marked as '{sname}') in output file.");
                }
                i++;
            }


            var B0 = new Basis(grd, 0);
            bndyMarkers = new SinglePhaseField[et2Name.Count + 1];

            int[,] Edge2GeomCell = grd.iGeomEdges.CellIndices;
            int[] G2L = grd.iGeomCells.GeomCell2LogicalCell;
            byte[] EdgeTags = grd.iGeomEdges.EdgeTags;
            int Jup = grd.iLogicalCells.NoOfLocalUpdatedCells;
            

            i = 0;
            foreach (var t in et2Name) { // loop over all different edge tag names...
                string name = t.Value;
                byte tag2color = t.Key;
                string sname = SanitizeName(name);

                int count = 0;

                var FI = new SinglePhaseField(B0, "Marker-" + sname);
                bndyMarkers[i] = FI;
                i++;

                var logMarker = new BitArray(Jup);

                for (int e = 0; e < EdgeTags.Length; e++) { // loop over edges...
                    byte tag_e = EdgeTags[e];

                    if (tag_e == tag2color) {
                        // mar cells next to edge e

                        

                        foreach (int jG in Edge2GeomCell.GetRow(e)) {
                            if (jG < 0)
                                continue;

                            // convert geometrical cell index to logical cell index
                            int jL;
                            if (G2L == null)
                                jL = jG;
                            else
                                jL = G2L[jG];

                            // color respective cell
                            if (jL < Jup) {
                                if (logMarker[jL] == false)
                                    count++;
                                logMarker[jL] = true;
                                FI.SetMeanValue(jL, tag2color);
                            }
                        }

                    }
                }

                this.BoundaryStat.Add((name, count));

            }

            var dummy = new SinglePhaseField(B0, "DummyData");
            bndyMarkers[bndyMarkers.Length - 1] = dummy;
            
            foreach (var f in bndyMarkers) {
                base.RegisterField(f, IOListOption.Always);
            }
            
        }




        /// <summary>
        /// 
        /// </summary>
        protected override void SetInitial(double t) {
            

        }


        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
           
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            Console.WriteLine("Mesh created:");

            Console.WriteLine("   Number of cells: " + this.Grid.NumberOfCells);
            Console.WriteLine("   Boundary conditions: ");
            foreach(var t in BoundaryStat) {
                Console.WriteLine("     " + t.EdgeTagName + ": " + t.Count.MPISum() + " edges.");
            }
            Console.WriteLine("   Grid ID    : " + this.Grid.ID);
            Console.WriteLine("   Grid Name  : " + this.Grid.Name);
            Console.WriteLine("   Description: " + this.Grid.Description);

            dt = 1e100;
            base.TerminationKey = true;

            Console.WriteLine("terminating app.");
            return dt;
        }
        
    

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "GridGen." + timestepNo;
            Tecplot.PlotFields(ArrayTools.Cat(bndyMarkers, MPIRank), filename, 0, superSampling);
        }
    }
}
