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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Solution.Timestepping;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Reinit.FastMarch {

    /// <summary>
    /// Holds various fast-marching tools
    /// </summary>
    public class FastMarchReinit {

        public FastMarchReinit(Basis __LevelSetBasis) {
            GridDat = (GridData)(__LevelSetBasis.GridDat);
            LevelSetBasis = __LevelSetBasis;
            LevelSetMapping = new UnsetteledCoordinateMapping(__LevelSetBasis);

            lsEllipt = new LocalSolver_Elliptic(__LevelSetBasis);
            lsGeom = new LocalSolver_Geometric(__LevelSetBasis);
            gradModule = new GradientModule();

            //This plotter can plot each single update of the resp. field. 
            plotter = new PlotStepper();
        }
        PlotStepper plotter;

        GridData GridDat;

        LocalSolver_Geometric lsGeom;
        LocalSolver_Elliptic lsEllipt;
        GradientModule gradModule;

        Basis LevelSetBasis;
        UnsetteledCoordinateMapping LevelSetMapping;

        /// <summary>
        /// Solves the Reinitialisation-problem locally, i.e. cell by cell.
        /// </summary>
        /// <param name="Accepted">
        /// </param>
        /// <param name="Recalc"></param>
        /// <param name="phi_accepted">
        /// Input: level-set field values for the accepted cells/vallues in these cells are unchanged on exit.
        /// Output: level-set-field values for cells in <see cref="Recalc"/>.
        /// </param>
        /// <param name="_sign"></param>
        /// <param name="gradPhi"></param>
        void LocalSolve(CellMask Accepted, CellMask Recalc, SinglePhaseField phi_accepted, VectorField<SinglePhaseField> gradPhi, double _sign, SinglePhaseField DiffusionCoeff) {

            // loop over 'Recalc' - cells...
            foreach (int jCell in Recalc.ItemEnum) {

                double mini, maxi;
                this.Stpw_geomSolver.Start();
                this.lsGeom.LocalSolve(jCell, Accepted.GetBitMaskWithExternal(), phi_accepted, _sign, out mini, out maxi);
                this.Stpw_geomSolver.Stop();
                this.Stpw_reinitSolver.Start();
                this.lsEllipt.LocalSolve(jCell, Accepted.GetBitMaskWithExternal(), phi_accepted, gradPhi);
                this.Stpw_reinitSolver.Stop();

                //bool q = LocalSolve_Iterative(jCell, Accepted.GetBitMask(), phi_accepted, gradPhi);
                //if(!q)
                //    Console.WriteLine("Local solver not converged.");

            }

        }

        Stopwatch Stpw_gradientEval = new Stopwatch();
        Stopwatch Stpw_total = new Stopwatch();
        Stopwatch Stpw_geomSolver = new Stopwatch();
        Stopwatch Stpw_extVelSolver = new Stopwatch();
        Stopwatch Stpw_reinitSolver = new Stopwatch();

        public void PrintInstrumentation() {
            double tot = Stpw_total.Elapsed.TotalSeconds;
            double grd = Stpw_gradientEval.Elapsed.TotalSeconds;
            double geo = Stpw_geomSolver.Elapsed.TotalSeconds;
            double ext = Stpw_extVelSolver.Elapsed.TotalSeconds;
            double rin = Stpw_reinitSolver.Elapsed.TotalSeconds;

            Console.WriteLine(" Total runtime: " + tot);
            Console.WriteLine("    gradient:   {0:0.##E-00} {1:0.##E-00}%", grd, 100 * grd / tot);
            Console.WriteLine("    geometric:  {0:0.##E-00} {1:0.##E-00}%", geo, 100 * geo / tot);
            Console.WriteLine("    reinit:     {0:0.##E-00} {1:0.##E-00}%", rin, 100 * rin / tot);
            Console.WriteLine("    extension:  {0:0.##E-00} {1:0.##E-00}%", ext, 100 * ext / tot);

            this.lsEllipt.PrintInstrumentation();
        }

        public double[] m_PhiAvg;

        public void GradientUpdate(SubGrid NEAr, SinglePhaseField Phi, VectorField<SinglePhaseField> GradPhi) {

            this.gradModule.GradientUpdate(NEAr, this.m_PhiAvg, Phi, GradPhi);

            //GradPhi.Clear(NEAr.VolumeMask);
            //GradPhi.Gradient(1, Phi, NEAr.VolumeMask);
        }

        public void AvgInit(SinglePhaseField Phi, CellMask _Accepted) {
            int J = this.GridDat.Cells.Count;
            double[] PhiAvg;
            if (m_PhiAvg == null) {
                PhiAvg = new double[J];
                m_PhiAvg = PhiAvg;
            } else {
                PhiAvg = m_PhiAvg;
            }

            foreach (int jCell in _Accepted.ItemEnum)
                PhiAvg[jCell] = Phi.GetMeanValue(jCell);
        }



        /// <summary>
        /// Reinit on un-cut cells.
        /// </summary>
        /// <param name="Phi">The level set</param>
        /// <param name="ReInitSpecies">Cell mask wich is to be reinitialized</param>
        /// <param name="sign">Sign of the level set for this <paramref name="ReInitSpecies"/></param>
        /// <param name="_Accepted">CellMask which is taken as boundray values</param>
        /// <param name="GradPhi">LEvel Set gradient</param>
        /// <param name="callBack">A delegate, which might be called after the execution of the reinitialization</param>
        public void Reinitialize(SinglePhaseField Phi, CellMask ReInitSpecies, double sign,
            CellMask _Accepted,
            //ConventionalDGField[] ExtProperty, double[][] ExtPropertyMin, double[][] ExtPropertyMax,
            VectorField<SinglePhaseField> GradPhi, Action<int> callBack) //
        {
            using (new FuncTrace()) {
                Tracer.InstrumentationSwitch = false; // lots of tracing on calls acting on singe cells causes massive overhead (up to 5x slower).
                Stpw_total.Start();

                SinglePhaseField DiffusionCoeff = new SinglePhaseField(new Basis(this.GridDat, 1), "DiffusionCoeff");


                // check args and init
                // ===================

                /*
                ExtVelSolver extVelSlv = null;
                if(ExtProperty != null) {
                    if(ExtProperty.Length != ExtPropertyMin.Length)
                        throw new ArgumentException();
                    if(ExtProperty.Length != ExtPropertyMax.Length)
                        throw new ArgumentException();

                    extVelSlv = new ExtVelSolver(ExtProperty[0].Basis);
                }
                 */

                BitArray Acceped_Mutuable = _Accepted.GetBitMask().CloneAs();
                BitArray Trial_Mutuable = ((_Accepted.AllNeighbourCells().Intersect(ReInitSpecies)).Except(_Accepted)).GetBitMask().CloneAs();
                BitArray Recalc_Mutuable = Trial_Mutuable.CloneAs();
                BitArray PosSpecies_Bitmask = ReInitSpecies.GetBitMask();

                int J = this.GridDat.Cells.Count;
                int D = this.GridDat.SpatialDimension;
                int N = this.LevelSetBasis.Length;

                double _sign = sign >= 0 ? 1.0 : -1.0;
                double[] PhiAvg = m_PhiAvg;
                if (PhiAvg == null)
                    throw new ApplicationException();

                foreach (int jCell in _Accepted.ItemEnum)
                    PhiAvg[jCell] = Phi.GetMeanValue(jCell);

                int NoOfNew;
                {
                    var Neu = ReInitSpecies.Except(_Accepted);
                    NoOfNew = Neu.NoOfItemsLocally;
                    Phi.Clear(Neu);
                    Phi.AccConstant(_sign, Neu);

                    foreach (int jCell in Neu.ItemEnum)
                        PhiAvg[jCell] = 1.0e10;
                }

                if (this.GridDat.MpiSize > 1)
                    throw new NotSupportedException("Currently not MPI parallel.");


                for (int d = 0; d < this.GridDat.SpatialDimension; d++) {
                    if (!GradPhi[d].Basis.Equals(Phi.Basis))
                        throw new ArgumentException("Level-set and level-set gradient field should have the same DG basis."); // ein grad niedriger wrürde auch genügen...
                }



                // perform marching...
                // ===================

                // update gradient for cut-cells
                GradPhi.Clear(_Accepted);
                GradPhi.Gradient(1.0, Phi, _Accepted);

                // marching loop../
                int cnt = 0;
                while (true) {
                    cnt++;



                    CellMask Recalc = new CellMask(this.GridDat, Recalc_Mutuable);
                    CellMask Accepted = new CellMask(this.GridDat, Acceped_Mutuable);
                    CellMask Trial = new CellMask(this.GridDat, Trial_Mutuable);

                    int NoOfTrial = Trial.NoOfItemsLocally;
                    int NoOfAccpt = Accepted.NoOfItemsLocally;
                    int NoOfRcalc = Recalc.NoOfItemsLocally;

                    if (Trial.NoOfItemsLocally <= 0) {
                        //Ploti(Recalc, Accepted, Trial, Phi, Phi_gradient, optEikonalOut, cnt);
                        break;
                    }

                    // Local solver for all 'Recalc'-cells
                    // --------------------------------------


                    if (Recalc.NoOfItemsLocally > 0) {
                        this.LocalSolve(Accepted, Recalc, Phi, GradPhi, _sign, DiffusionCoeff);
                    }

                    // find the next cell to accept
                    // ----------------------------

                    // get mean value in all cells
                    foreach (int jCell in Recalc.ItemEnum) {
                        PhiAvg[jCell] = Phi.GetMeanValue(jCell);
                        Recalc_Mutuable[jCell] = false;
                    }

                    //Ploti(Recalc, Accepted, Trial, Phi, Phi_gradient, optEikonalOut, cnt);

                    // find trial-cell with minimum average value
                    // this should be done with heap-sort (see fast-marching algorithm)
                    int jCellAccpt = int.MaxValue;
                    double TrialMin = double.MaxValue;
                    foreach (int jCell in Trial.ItemEnum) {
                        if (PhiAvg[jCell] * _sign < TrialMin) {
                            TrialMin = PhiAvg[jCell] * _sign;
                            jCellAccpt = jCell;
                        }
                    }

                    if (callBack != null)
                        callBack(cnt);

                    /*
                    // update the gradient
                    // -------------------

                    this.Stpw_gradientEval.Start();
                    gradModule.GradientUpdate(jCellAccpt, Acceped_Mutuable, Phi, GradPhi);
                    this.Stpw_gradientEval.Stop();

                    /*
                    // solve for the extension properties
                    // ----------------------------------

                    if(ExtProperty != null) {
                        int[] Neight, dummy33;
                        GridDat.Cells.GetCellNeighbours(jCellAccpt, GridData.CellData.GetCellNeighbours_Mode.ViaEdges, out Neight, out dummy33);

                        for(int iComp = 0; iComp < ExtProperty.Length; iComp++) {

                            ExtPropertyMax[iComp][jCellAccpt] = -double.MaxValue;
                            ExtPropertyMin[iComp][jCellAccpt] = double.MaxValue;

                            foreach(int jNeig in Neight) {
                                if(Acceped_Mutuable[jNeig]) {
                                    ExtPropertyMax[iComp][jCellAccpt] = Math.Max(ExtPropertyMax[iComp][jCellAccpt], ExtPropertyMax[iComp][jNeig]);
                                    ExtPropertyMin[iComp][jCellAccpt] = Math.Min(ExtPropertyMin[iComp][jCellAccpt], ExtPropertyMin[iComp][jNeig]);
                                }
                            }

                            this.Stpw_extVelSolver.Start();
                            extVelSlv.ExtVelSolve_Far(Phi, GradPhi, ExtProperty[iComp], ref ExtPropertyMin[iComp][jCellAccpt], ref ExtPropertyMax[iComp][jCellAccpt], jCellAccpt, Accepted, _sign);
                            this.Stpw_extVelSolver.Stop();
                        }
                    }
                    /*
                    {
                        int[] Neight, dummy33;
                        GridDat.Cells.GetCellNeighbours(jCellAccpt, GridData.CellData.GetCellNeighbours_Mode.ViaEdges, out Neight, out dummy33);
                        foreach(int jNeig in Neight) {
                            if(Acceped_Mutuable[jNeig]) {
                                plotDependencyArrow(cnt, jCellAccpt, jNeig);
                            }
                        }
                    }
                    */

                    // the mimium is moved to accepted
                    // -------------------------------
                    Acceped_Mutuable[jCellAccpt] = true;
                    Trial_Mutuable[jCellAccpt] = false;
                    Recalc_Mutuable[jCellAccpt] = false;
                    NoOfNew--;

                    // recalc on all neighbours
                    // ------------------------
                    int[] Neighs, dummy;
                    this.GridDat.GetCellNeighbours(jCellAccpt, GetCellNeighbours_Mode.ViaEdges, out Neighs, out dummy);
                    foreach (int jNeig in Neighs) {
                        if (!Acceped_Mutuable[jNeig] && PosSpecies_Bitmask[jNeig]) {
                            Trial_Mutuable[jNeig] = true;
                            Recalc_Mutuable[jNeig] = true;
                        }
                    }

                }

                if (NoOfNew > 0)
                    throw new ArithmeticException("Unable to perform reinitialization for all requested cells - maybe they are not reachable from the initialy 'accepted' domain?");

                //PlottAlot("dependencies.csv");

                Tracer.InstrumentationSwitch = true;
                Stpw_total.Stop();
            }
        }



        /// <summary>
        /// Extension velocity on un-cut cells.
        /// </summary>
        /// <param name="Phi">Input; a signed-distance field.</param>
        /// <param name="Domain">Cells in which the extension velocity should be computed.</param>
        /// <param name="cut">Cells in which a value for the extension property is given, must be
        /// disjoint from and share a boundary with <paramref name="Domain"/>.
        /// </param>
        /// <param name="ExtProperty">
        /// Input/Output: on cells in <paramref name="cut"/>, valid values for the extension property,
        /// i.e. boundary values for the extension problem.
        /// On exit, the result of the marching algorithm is stored in the <paramref name="Domain"/>-cells.
        /// Multiple extension properties can be constructed at once.
        /// </param>
        /// <param name="ExtPropertyMin">
        /// Input/Output: Helper array to check the minimum/maximum principle of the extension velocity problem
        ///   - 1st index: extension property index
        ///   - 2nd index: cell index
        /// On entry, the minimum values for cells in <paramref name="cut"/> must contain valid entries.
        /// </param>
        /// <param name="ExtPropertyMax">
        /// Analogous to <paramref name="ExtPropertyMin"/>.
        /// </param>
        /// <param name="GradPhi">
        /// Helper variable to compute the gradient values of <paramref name="Phi"/>.
        /// </param>
        public void ConstructExtension(SinglePhaseField Phi, CellMask Domain, CellMask cut,
            ConventionalDGField[] ExtProperty, double[][] ExtPropertyMin, double[][] ExtPropertyMax,
            VectorField<SinglePhaseField> GradPhi, int TimestepNo, bool plotMarchingSteps = false) //
        {
            using (new FuncTrace()) {
                Tracer.InstrumentationSwitch = false; // lots of tracing on calls acting on singe cells causes massive overhead (up to 5x slower).
                Stpw_total.Start();

                // check args and init
                // ===================


                //ExtVelSolver extVelSlv = null;
                ExtVelSolver_Geometric extVelSlv = null;
                if (ExtProperty != null) {
                    if (ExtProperty.Length != ExtPropertyMin.Length)
                        throw new ArgumentException();
                    if (ExtProperty.Length != ExtPropertyMax.Length)
                        throw new ArgumentException();


                    //extVelSlv = new ExtVelSolver(ExtProperty[0].Basis);
                    extVelSlv = new ExtVelSolver_Geometric(ExtProperty[0].Basis);
                }

                BitArray Accepted_Mutuable = cut.GetBitMask().CloneAs();

                int J = this.GridDat.Cells.Count;
                int D = this.GridDat.SpatialDimension;
                int N = this.LevelSetBasis.Length;

                int[] DomainCellIndices = Domain.ItemEnum.ToArray();
                double[] PhiAvg = new double[DomainCellIndices.Length];
                {
                    int L = DomainCellIndices.Length;
                    for (int jSub = 0; jSub < L; jSub++) {
                        int jCell = DomainCellIndices[jSub];
                        PhiAvg[jSub] = Phi.GetMeanValue(jCell);

                    }
                    Array.Sort(PhiAvg, DomainCellIndices);
                }

                if (this.GridDat.MpiSize > 1)
                    throw new NotSupportedException("Currently not MPI parallel.");


                int[] PosDomain, NegDomain;
                {
                    int median = 0;
                    for (; median < DomainCellIndices.Length; median++) {
                        if (PhiAvg[median] >= 0)
                            break;
                    }

                    NegDomain = new int[median];
                    for (int i = 0; i < median; i++) {
                        NegDomain[i] = DomainCellIndices[median - i - 1];
                    }

                    PosDomain = new int[DomainCellIndices.Length - median];
                    Array.Copy(DomainCellIndices, median, PosDomain, 0, PosDomain.Length);

                    Debug.Assert(PosDomain.Length + NegDomain.Length == DomainCellIndices.Length);
                }
                if (plotMarchingSteps) {
                    this.plotter.setup(new DGField[] { Phi, ExtProperty[0], ExtProperty[1] }, TimestepNo);
                    this.plotter.plotstep(Accepted_Mutuable);
                }


                // perform marching...
                // ===================

                // marching loop..
                for (int iMinusPlus = -1; iMinusPlus <= 1; iMinusPlus += 2) {
                    double _sign = iMinusPlus;
                    int[] _Domain;
                    switch (iMinusPlus) {
                        case -1:
                            _Domain = NegDomain;
                            break;
                        case +1:
                            _Domain = PosDomain;
                            break;
                        default:
                            throw new Exception();
                    }

                    for (int iSub = 0; iSub < _Domain.Length; iSub++) {

                        CellMask Accepted = new CellMask(this.GridDat, Accepted_Mutuable);

                        int jCellAccpt = _Domain[iSub];
                        //this.Stpw_gradientEval.Start();
                        //gradModule.GradientUpdate(jCellAccpt, Acceped_Mutuable, Phi, GradPhi);
                        //this.Stpw_gradientEval.Stop();


                        // solve for the extension properties
                        // ----------------------------------

                        if (ExtProperty != null) {
                            int[] Neighb, dummy33;
                            GridDat.GetCellNeighbours(jCellAccpt, GetCellNeighbours_Mode.ViaEdges, out Neighb, out dummy33);

                            // solve for each component seperately
                            for (int iComp = 0; iComp < ExtProperty.Length; iComp++) {

                                ExtPropertyMax[iComp][jCellAccpt] = -double.MaxValue;
                                ExtPropertyMin[iComp][jCellAccpt] = double.MaxValue;

                                foreach (int jNeig in Neighb) {
                                    if (Accepted_Mutuable[jNeig]) {
                                        ExtPropertyMax[iComp][jCellAccpt] = Math.Max(ExtPropertyMax[iComp][jCellAccpt], ExtPropertyMax[iComp][jNeig]);
                                        ExtPropertyMin[iComp][jCellAccpt] = Math.Min(ExtPropertyMin[iComp][jCellAccpt], ExtPropertyMin[iComp][jNeig]);
                                    }
                                }

                                this.Stpw_extVelSolver.Start();
                                //extVelSlv.ExtVelSolve_Far(Phi, GradPhi, ExtProperty[iComp], ref ExtPropertyMin[iComp][jCellAccpt], ref ExtPropertyMax[iComp][jCellAccpt], jCellAccpt, Accepted, _sign);
                                extVelSlv.ExtVelSolve_Geometric(Phi, ExtProperty[iComp], Accepted_Mutuable, jCellAccpt, _sign);
                                this.Stpw_extVelSolver.Start();
                            }
                        }
                        Accepted_Mutuable[jCellAccpt] = true;
                        if (plotMarchingSteps) {
                            this.plotter.plotstep(Accepted_Mutuable);
                        }
                    }
                }
                Tracer.InstrumentationSwitch = true;
                Stpw_total.Stop();
            }
        }

        /// <summary>
        /// Reinitializes <paramref name="Phi"/> on <paramref name="reinitField"/> using the <paramref name="Accepted"/> cells as start value. 
        /// </summary>
        /// <param name="Phi">Level Set</param>
        /// <param name="Accepted">Given Field</param>
        /// <param name="NegativeField">Field, on which the level set function is negative</param>
        /// <param name="reinitField">Field , on which the level set function should be initialized</param>
        public void FirstOrderReinit(SinglePhaseField Phi, CellMask Accepted, CellMask NegativeField, CellMask reinitField) {

            //Idea: Add options : FirstOrder, Elliptic, Geometric, Iterative. That is, include the existing Local Solvers.
            CellMask ReinitField;
            if (reinitField == null) {
                ReinitField = CellMask.GetFullMask(Phi.GridDat);
            } else {
                ReinitField = reinitField;
            }
            //Build Local Solver that solves the Eikonal in each cell.
            FastMarching.ILocalSolver localSolver = new FastMarching.LocalMarcher.LocalMarcher_2DStructured(Phi.Basis);

            //Build Global Solver that marches through all cells
            FastMarching.GlobalMarcher.CellMarcher fastMarcher = new FastMarching.GlobalMarcher.CellMarcher(Phi.Basis, localSolver);
            
            //Solve
            fastMarcher.Reinit(Phi, Accepted, ReinitField);

            //Invert Negative Domain that is part of ReinitField
            Phi.Scale(-1, NegativeField.Intersect(ReinitField).Except(Accepted));
        }

        
        /// <summary>
        /// Reinitializes levelset in <paramref name="_LevelSetTracker"/> on a near field with width <paramref name="nearFieldWidth"/>
        /// </summary>
        /// <param name="_LevelSetTracker"></param>
        /// <param name="NameNegativeSpecies">Name of the negative species of Levelset</param>
        /// <param name="nearFieldWidth">Width of band arround cutcells</param>
        public void FirstOrderReinit(LevelSetTracker _LevelSetTracker, string NameNegativeSpecies, int nearFieldWidth) {
            //Extract Data from LevelSetTracker
            CellMask accepted = _LevelSetTracker.Regions.GetCutCellMask();
            CellMask negativeDomain = _LevelSetTracker.Regions.GetSpeciesMask(NameNegativeSpecies);
            CellMask nearField = _LevelSetTracker.Regions.GetNearFieldMask(nearFieldWidth);

            //Solve
            FirstOrderReinit((SinglePhaseField)_LevelSetTracker.LevelSets[0], accepted, negativeDomain, nearField);
        }


        /// <summary>
        /// Reinitializes levelset in <paramref name="LevelSetTracker"/> on all cells except the cut cells
        /// </summary>
        /// <param name="LevelSetTracker"></param>
        /// <param name="NameNegativeSpecies">Name of the negative species of Levelset</param>
        public void FirstOrderReinit(LevelSetTracker LevelSetTracker, string NameNegativeSpecies) {

            //Extract Data from LevelSetTracker
            CellMask accepted = LevelSetTracker.Regions.GetCutCellMask();
            CellMask negativeDomain = LevelSetTracker.Regions.GetSpeciesMask(NameNegativeSpecies);
            CellMask ReInitField = CellMask.GetFullMask(LevelSetTracker.GridDat).Except(accepted);

            //Solve
            FirstOrderReinit((SinglePhaseField)LevelSetTracker.LevelSets[0], accepted, negativeDomain, ReInitField);
        }
    }
}
