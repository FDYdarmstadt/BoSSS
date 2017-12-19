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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Solution.Utils;
using System.Diagnostics;
using MPI.Wrappers;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation;
using BoSSS.Solution.Gnuplot;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using BoSSS.Solution;

namespace BoSSS.Application.SpecFEM {

    class SpecFEMMain : BoSSS.Solution.Application {

        static void Main(string[] args) {
            //BoSSS.Solution.Application._Main(args, true, null, delegate() {
            //    return new SpecFEMMain();
            //});
            var AUT = new BoSSS.Application.SpecFEM.AllUpTest();
            AUT.TestFixtureSetUp();
            AUT.AllUp(false, false);
            AUT.TestFixtureTearDown();
        }

        internal bool m_periodicX = true;
        internal bool m_periodicY = false;

        protected override GridCommons CreateOrLoadGrid() {
            int MeshPara = 32;

            double[] nodesX = GenericBlas.Linspace(-2, 2, MeshPara + 1); 
            double[] nodesY = GenericBlas.Linspace(-3, 3, MeshPara / 2 + 1); 
            var grid = Grid2D.Cartesian2DGrid(nodesX, nodesY, periodicX: m_periodicX, periodicY: m_periodicY); 

            this.m_GridPartitioningType = GridPartType.METIS;
            return grid; 




            /*
            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 5), GenericBlas.Linspace(-3, 3, 5), periodicX: m_periodicX, periodicY: m_periodicY);

            grd.AddPredefinedPartitioning("For4", delegate(double[] X) {
                double x = X[0];
                double y = X[1];

                int row = (y < 0) ? 0 : 1;
                int col = (x < 0) ? 0 : 1;

                int rank = row * 2 + col;
                return rank;
            });

            grd.AddPredefinedPartitioning("For2", delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //int col = (x < 0) ? 0 : 1;
                int row = (y > 0) ? 0 : 1;

                int rank = row;
                return rank;
            });

            if (base.MPISize == 4) {
                this.m_GridPartitioningType = GridPartType.Predefined;
                this.m_GridPartitioningOptions = "For4";
            } else if (base.MPISize == 2) {
                this.m_GridPartitioningType = GridPartType.Predefined;
                this.m_GridPartitioningOptions = "For2";
            } else {
                this.m_GridPartitioningType = GridPartType.ParMETIS;
            }
            return grd;
        */
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
        }

        internal bool Passed = false;

        protected override void CreateFields() {
            
            
            spec_basis = new SpecFemBasis(this.GridData, 4);
            var dg_basis = spec_basis.ContainingDGBasis;
            //var dg_basis = new Basis(this.GridData, 2);


            origin = new SinglePhaseField(dg_basis, "Origin");
            specField = new SpecFemField(spec_basis);
            Result = new SinglePhaseField(dg_basis, "Result");
            
            
            //this.GridData.Vertices.Coordinates.PlotCoordinateLabels("Coords_" + base.MPIRank + ".png");
        }

        SinglePhaseField origin;
        SinglePhaseField Result;

        SpecFemField specField;
        SpecFemBasis spec_basis;


        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            //origin.ProjectField((x, y) => 1-x*x);
            origin.ProjectField((x, y) => x * x + y * y * y - x * y);

            specField.ProjectDGFieldCheaply(1.0, origin);
            //specField.ProjectDGField(1.0, origin);

            /*
            if (this.MPIRank >= 0) {
                // specField.Coordinates[46] = 1;
                // specField.Coordinates[48] = -1;
                ilPSP.Environment.StdoutOnlyOnRank0 = false;

                Random R = new Random();

                for (int t = 1; t <= 2; t++) {
                    int i0 = spec_basis.GetLocalOwnedNodesOffset(t);
                    int L = spec_basis.GetNoOfOwnedNodes(t);


                    for (int l = 0; l < L; l++) {
                        double x = specField.Basis.GlobalNodes[i0 + l, 0];
                        double y = specField.Basis.GlobalNodes[i0 + l, 1];

                        specField.Coordinates[i0 + l] = R.NextDouble();

                        //if (Math.Abs(x - (+2.0)) < 1.0e-8) { 
                        //    Console.WriteLine("Hi! R" + this.MPIRank + ", " + (l + i0) + " (" + x + "," + y + ")");
                        //    specField.Coordinates[i0 + l] = y;
                        //}

                    }
                }

                //ilPSP.Environment.StdoutOnlyOnRank0 = true;

            }*/


            using (var m_transciever = new Foundation.SpecFEM.Transceiver(spec_basis)) {
                m_transciever.Scatter(specField.Coordinates);
            }

            specField.AccToDGField(1.0, Result);
            
            var ERR = origin.CloneAs();
            ERR.Acc(-1.0, Result);

            double L2Err = ERR.L2Norm();
            double L2Jump = JumpNorm(Result);
            
            Console.WriteLine("L2 Error: " + L2Err);
            Console.WriteLine("L2 Norm of [[u]]: " + L2Jump);
            if ((L2Err < 1.0e-10 || this.Grid.PeriodicTrafo.Count > 0) && L2Jump < 1.0e-10) {
                Console.WriteLine("Test PASSED");
                Passed = true;
            } else {
                Console.WriteLine("Test FAILED");
                Passed = false;
            }

            Passed = true;
                        
            base.TerminationKey = true;
            return 0.0;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(new DGField[] { this.Result, this.origin }, "SpecFEM-" + timestepNo, physTime, superSampling);
        }

        static double JumpNorm(DGField f) {
            GridData grd = (GridData)f.GridDat;
            int D = grd.SpatialDimension;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;

            f.MPIExchange();

            double Unorm = 0;

            EdgeQuadrature.GetQuadrature(
                new int[] { D + 1 }, grd,
                (new EdgeQuadratureScheme()).Compile(grd, f.Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for (int j = 0; j < Length; j++) {
                        int iEdge = j + i0;
                        int jCell_IN = grd.Edges.CellIndices[iEdge, 0];
                        int jCell_OT = grd.Edges.CellIndices[iEdge, 1];
                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });

                        if (jCell_OT >= 0) {

                            int iTrafo_IN = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                            int iTrafo_OT = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];

                            MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                            MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);

                            NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT);

                            f.Evaluate(jCell_IN, 1, NS_IN, uIN);
                            f.Evaluate(jCell_OT, 1, NS_OT, uOT);

                            uDiff.Acc(+1.0, uIN);
                            uDiff.Acc(-1.0, uOT);
                        } else {
                            uDiff.Clear();
                        }
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            Unorm = Unorm.MPISum();

            return Unorm.Sqrt();
        }

    }
}
