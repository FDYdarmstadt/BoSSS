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
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.LevelSetTools.Smoothing;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution.LevelSetTools.Advection {
    
    /// <summary>
    /// Extension Velocity by a particle method.
    /// </summary>
    public static class NarrowMarchingBand {
        
        const int iLevSet = 0;

        /*
        /// <summary>
        /// Level-Set evolution using a scalar extension velocity.
        /// </summary>
		/// <param name="Velocity">Input; velocity field to compute the extension velocity from</param>
        /// <param name="OldLevSet">Input, level set fied at current timestep </param>
        /// <param name="NewLevelSet">New level set.</param>
        /// <param name="ExtVel">Output, the computed extension velocity.</param>
        /// <param name="LevelSetGrad"></param>
        /// <param name="Tracker">Tracker.</param>
        /// <param name="dt">Dt.</param>
        public static void Evolve_ScalarExtvel(double dt, LevelSetTracker Tracker,
            SinglePhaseField OldLevSet, SinglePhaseField NewLevelSet, VectorField<SinglePhaseField> LevelSetGrad,
            VectorField<SinglePhaseField> Velocity, SinglePhaseField ExtVel) //
        {
            GridData gdat = Tracker.GridDat;


            // find (remember old) cut cell, positive and negative domain
            // ----------------------------------------------------------
            SubGrid CCgrid = Tracker.GetCutCellSubgrid4LevSet(iLevSet);
            CellMask CC = CCgrid.VolumeMask;
            CellMask Pos = Tracker.GetLevelSetWing(iLevSet, +1.0).VolumeMask.Except(CC);
            CellMask Neg = Tracker.GetLevelSetWing(iLevSet, -1.0).VolumeMask.Except(CC);
            
            SubGrid PosGrid = new SubGrid(Pos);
            SubGrid NegGrid = new SubGrid(Neg);

            EdgeMask touching = PosGrid.BoundaryEdgesMask.Intersect(NegGrid.BoundaryEdgesMask);
            if(touching.NoOfItemsLocally > 0)
                throw new ArithmeticException("Error in level-set topology.");


            // set values in positive and negative region to +1 and -1
            // -------------------------------------------------------

            //SinglePhaseField OldLevSet = (SinglePhaseField)(Tracker.LevelSets[iLevSet]);
            if(!NewLevelSet.Basis.Equals(OldLevSet.Basis))
                throw new ArgumentException("Basis of new and old level-set must match.");

            if(object.ReferenceEquals(OldLevSet, NewLevelSet)) {
                OldLevSet = OldLevSet.CloneAs();
            } else {
                NewLevelSet.Clear(CC);
                NewLevelSet.Acc(1.0, OldLevSet, CC);
            }
            
            NewLevelSet.Clear(Pos);
            NewLevelSet.AccConstant(1.0, Pos);

            NewLevelSet.Clear(Neg);
            NewLevelSet.AccConstant(-1.0, Neg);
                        
            // evolution on cut-cells
            // ----------------------
            ExtVel.Clear();
            int J = gdat.Cells.NoOfLocalUpdatedCells;
            double[] ExtVelMin = new double[J], ExtVelMax = new double[J];
            //ConstructScalarExtVel_CC(Tracker, ExtVel, Velocity, OldLevSet, LevelSetGrad, ExtVelMin, ExtVelMax);
                                   
            LevelSetGrad.Clear(CC);
            LevelSetGrad.Gradient(1.0, NewLevelSet, CC);
                        
            NewLevelSet.ProjectFunction(-dt,
                LevelSetEvoTerm_Scalar,
                (new CellQuadratureScheme(true, CC)),
                ArrayTools.Cat<DGField>(new DGField[] { ExtVel }, LevelSetGrad));
            
            LevelSetGrad.Clear(CC.Complement());

            VectorField<SinglePhaseField> LevelSetMeanGrad = new VectorField<SinglePhaseField>(gdat.SpatialDimension.ForLoop(d => new SinglePhaseField(new Basis(gdat, 1), "MeanGrad_" + d)));
            for(int d = 0; d < gdat.SpatialDimension; d++) {
                LevelSetMeanGrad[d].AccLaidBack(1.0, LevelSetGrad[d]);
            }

            ExtVel.Clear(CCgrid.VolumeMask.Complement());
            
            // detect whether the levset has moved to another cell, and 
			// perform Reinit & Extension on those cells.
            // ------------------------------------------------------------------------------------------

            Reinit.FastMarch.FastMarchReinit hack = null;
            CellMask AllExtensionMask = null;
                        
            
            BitArray PosExtensionBitMask = FindReinitDomain(NewLevelSet, CCgrid, PosGrid, +1);
            if(PosExtensionBitMask != null) {
                if(hack == null)
                    hack = new Reinit.FastMarch.FastMarchReinit(NewLevelSet.Basis);

                CellMask PosExtensionMask = new CellMask(gdat, PosExtensionBitMask);
                hack.Reinitialize(OldLevSet, PosExtensionMask, CC,
                    new ConventionalDGField[] { ExtVel }, new double[][] { ExtVelMin }, new double[][] { ExtVelMax }, 
                    LevelSetGrad, +1, null);

                //PosExtensionMask.ToTxtFile("posExtension.csv", WriteHeader:false);

                AllExtensionMask = PosExtensionMask;
            }
                        
            BitArray NegExtensionBitMask = FindReinitDomain(NewLevelSet, CCgrid, NegGrid, -1);
            if(NegExtensionBitMask != null) {
                if(hack == null)
                    hack = new Reinit.FastMarch.FastMarchReinit(NewLevelSet.Basis);

                CellMask NegExtensionMask = new CellMask(gdat, NegExtensionBitMask);
                hack.Reinitialize(OldLevSet, NegExtensionMask, CC, 
                    new ConventionalDGField[] { ExtVel }, new double[][] { ExtVelMin }, new double[][] { ExtVelMax }, 
                    LevelSetGrad, -1, null);

                //NegExtensionMask.ToTxtFile("negExtension.csv", WriteHeader: false);
                
                if(AllExtensionMask == null)
                    AllExtensionMask = NegExtensionMask;
                else
                    AllExtensionMask = AllExtensionMask.Union(NegExtensionMask);
            }

           

            // perform extension velocity on extended cells (if necessary)
            // -----------------------------------------------------------

            if(AllExtensionMask != null) {
                NewLevelSet.Clear(AllExtensionMask);
                NewLevelSet.Acc(1.0, OldLevSet, AllExtensionMask);
                NewLevelSet.ProjectFunction(-dt, 
                    LevelSetEvoTerm_Scalar,
                    (new CellQuadratureScheme(true, AllExtensionMask)),
                    ArrayTools.Cat<DGField>(new DGField[] { ExtVel }, LevelSetGrad));

            }

            
            
            // finally, apply the stabilization
            // --------------------------------

            SubGrid stabiRegion;
            if(AllExtensionMask == null)
                stabiRegion = CCgrid;
            else
                stabiRegion = new SubGrid(CC.Union(AllExtensionMask));

            var jp = new JumpPenalization();
            jp.ImplicitEuler(dt*1, stabiRegion, NewLevelSet);

        }


        /// <summary>
        /// The term $ s | \nabla \varphi |$, where $s$ denotes the scalar extension velocity.
        /// </summary>
        static double LevelSetEvoTerm_Scalar(double[] X, double[] U, int jCell) {
            double s = U[0];
            double normGradPhi = 0;
            for(int i = 1; i < U.Length; i++) {
                normGradPhi += U[i].Pow2();
            }
            normGradPhi = Math.Sqrt(normGradPhi);
            return s * normGradPhi;
        }
        */

        /// <summary>
        /// Level-Set evolution using a scalar or vector extension velocity (mark 2).
        /// </summary>
        /// <param name="Velocity">Input: velocity field to compute the extension velocity.</param>
        /// <param name="OldLevSet">Input: level set fied at current time-step.</param>
        /// <param name="NewLevelSet">New level set.</param>
        /// <param name="ExtVel">
        /// Output, the computed extension velocity;
        /// For one entry, the algorithm uses a scalar extension velocity.
        /// For D entries, the algorithm uses a vector extension velocity.
        /// </param>
        /// <param name="LevelSetGrad">
        /// Output/Work: velocity gradient on near- and cut-cells;
        /// </param>
        /// <param name="Tracker"></param>
        /// <param name="dt">Size of time-step.</param>
        /// <param name="HMForder"></param>
        public static void Evolve_Mk2(
            double dt, LevelSetTracker Tracker,
            SinglePhaseField OldLevSet, SinglePhaseField NewLevelSet, VectorField<SinglePhaseField> LevelSetGrad,
            ConventionalDGField[] Velocity, SinglePhaseField[] ExtVel,
            int HMForder, int TimestepNo = 0, bool plotMarchingSteps = false) //
        {
            GridData gdat = Tracker.GridDat;
            int D = gdat.SpatialDimension;
            int J = gdat.Cells.NoOfLocalUpdatedCells;
            
            if(Velocity.Length != D)
                throw new ArgumentException("wrong spatial dimension of velocity");

            if(Tracker.NearRegionWidth < 1)
                throw new NotSupportedException("At least a width of 1 is required.");
            if(!NewLevelSet.Basis.Equals(OldLevSet.Basis))
                throw new ArgumentException("Basis of new and old level-set must match.");

            bool ComponentMode; // whether scalar or vector extension velocity is used
            if(ExtVel.Length == D) {
                // 
                ComponentMode = true;

            } else {
                ComponentMode = false;
                if(ExtVel.Length != 1)
                    throw new ArgumentException();
            }


            SubGrid NEARgrid = Tracker.Regions.GetNearFieldSubgrid4LevSet(iLevSet, 1);
            CellMask NEAR = NEARgrid.VolumeMask;

            SubGrid CCgrid = Tracker.Regions.GetCutCellSubgrid4LevSet(iLevSet);
            CellMask CC = CCgrid.VolumeMask;
            BitArray CCBitMask = CC.GetBitMask();
            CellMask Pos = Tracker.Regions.GetLevelSetWing(iLevSet, +1.0).VolumeMask.Except(CC);
            CellMask Neg = Tracker.Regions.GetLevelSetWing(iLevSet, -1.0).VolumeMask.Except(CC);
            BitArray PosBitMask = Pos.GetBitMask();
            BitArray NegBitMask = Neg.GetBitMask();

            SubGrid PosGrid = new SubGrid(Pos);
            SubGrid NegGrid = new SubGrid(Neg);

            EdgeMask touching = PosGrid.BoundaryEdgesMask.Intersect(NegGrid.BoundaryEdgesMask);
            if(touching.NoOfItemsLocally > 0)
                throw new ArithmeticException("Error in level-set topology.");

            if(object.ReferenceEquals(OldLevSet, NewLevelSet)) {
                NewLevelSet = OldLevSet.CloneAs();
            } else {
                NewLevelSet.Clear(NEAR);
                NewLevelSet.Acc(1.0, OldLevSet, NEAR);
            }

            
            
            // set values in positive and negative FAR region to +1 and -1
            // -----------------------------------------------------------

            CellMask PosFAR = Pos.Except(NEAR);
            CellMask NegFAR = Neg.Except(NEAR);

            NewLevelSet.Clear(PosFAR);
            NewLevelSet.AccConstant(1.0, PosFAR);

            NewLevelSet.Clear(NegFAR);
            NewLevelSet.AccConstant(-1.0, NegFAR);


            // recalculate Level-Set Gradient
            //-------------------------------
            LevelSetGrad.Clear(NEAR);
            LevelSetGrad.Gradient(1.0, NewLevelSet, NEAR);


            // identify cells for Reinit
            // -------------------------
            
            BitArray ReinitPosBitmask = new BitArray(J);
            BitArray ReinitNegBitmask = new BitArray(J);
            BitArray KnownBitmask = new BitArray(J);
            int NoOfPosReinit = 0, NoOfNegReinit = 0;
            foreach(var j in NEAR.ItemEnum) {
                double GradNorm = 0;
                for(int d = 0; d < D; d++) {
                    GradNorm += LevelSetGrad[d].Coordinates.GetRow(j).L2NormPow2();
                }

                if (GradNorm < 1.0e-12) {
                    //if (!CCBitMask[j] && (GradNorm < 0.9 || GradNorm > 1.1)) { 

                    // the level-set field is flat (+1 or -1) in cell j
                    // it must be initialized

                    if (CCBitMask[j])
                        throw new ArithmeticException("Found vanishing level-set gradient on cut-cells.");
                    if(PosBitMask[j] == NegBitMask[j])
                        throw new ArithmeticException("Level set in un-cut cell seems to be positive and negative at the same time -- something wrong here.");

                    if(PosBitMask[j]) {
                        ReinitPosBitmask[j] = true;
                        NoOfPosReinit++;
                    }
                    if(NegBitMask[j]) {
                        ReinitNegBitmask[j] = true;
                        NoOfNegReinit++;
                    }
                } else {
                    KnownBitmask[j] = true;
                }
            }

            //Console.WriteLine("No of pos/neg reinit: {0}, {1}", NoOfPosReinit, NoOfNegReinit);

            CellMask ReintPos = new CellMask(gdat, ReinitPosBitmask);
            CellMask ReintNeg = new CellMask(gdat, ReinitNegBitmask);
            CellMask Known = new CellMask(gdat, KnownBitmask);


            // perform Reinitialization
            // ------------------------
                        
            var marcher = new Reinit.FastMarch.FastMarchReinit(NewLevelSet.Basis);

            marcher.AvgInit(NewLevelSet, Known);

            if(NoOfPosReinit > 0)
                marcher.Reinitialize(NewLevelSet, ReintPos, +1, Known, LevelSetGrad, null);
            if(NoOfNegReinit > 0)
                marcher.Reinitialize(NewLevelSet, ReintNeg, -1, Known, LevelSetGrad, null);

            // save reinit at OLD levelset, to prevent doing it multiple times
            OldLevSet.Clear(ReintPos);
            OldLevSet.Clear(ReintNeg);
            OldLevSet.Acc(1.0, NewLevelSet, ReintPos);
            OldLevSet.Acc(1.0, NewLevelSet, ReintNeg);



            // bring gradient up to date
            marcher.GradientUpdate(NEARgrid, NewLevelSet, LevelSetGrad);


            // compute velocity extension
            // --------------------------

            double[][] ExtVelMin = new double[ExtVel.Length][];
            double[][] ExtVelMax = new double[ExtVel.Length][];
            for(int i = 0; i < ExtVel.Length; i++) {
                ExtVelMin[i] = new double[J];
                ExtVelMax[i] = new double[J];
            }


            // marching
            // ========

            // first, on cut cells
            //MultidimensionalArray x0;
            //x0 = ConstructExtVel_GOM_CC(Tracker, ExtVel, Velocity.ToArray(), NewLevelSet, LevelSetGrad, ExtVelMin, ExtVelMax);

            ConstructExtVel_PDE(Tracker, CCgrid, ExtVel, Velocity.ToArray(), NewLevelSet, LevelSetGrad, ExtVelMin, ExtVelMax, HMForder);

            // then, on the rest of the domain
            marcher.ConstructExtension(NewLevelSet, NEAR.Except(CC), CC, ExtVel, ExtVelMin, ExtVelMax, LevelSetGrad, TimestepNo, plotMarchingSteps);


            // construct via PDE all cells at once
            // ===================================

            // for the elliptic solve far field need to be Double.MaxValue or Double.MinValue
            //NewLevelSet.Clear(PosFAR);
            //NewLevelSet.AccConstant(Double.MaxValue, PosFAR);

            //NewLevelSet.Clear(NegFAR);
            //NewLevelSet.AccConstant(Double.MinValue, NegFAR);

            //ConstructExtVel_PDE(Tracker, NEARgrid, ExtVel, Velocity.ToArray(), NewLevelSet, LevelSetGrad, ExtVelMin, ExtVelMax, HMFversion, HMForder);

            //// set far field back to +1.0/-1.0
            //NewLevelSet.Clear(PosFAR);
            //NewLevelSet.AccConstant(1.0, PosFAR);

            //NewLevelSet.Clear(NegFAR);
            //NewLevelSet.AccConstant(-1.0, NegFAR);


            for (int d = 0; d < ExtVel.Length; d++)
                ExtVel[0].CheckForNanOrInf(true, true, true);

            //marcher.PrintInstrumentation();




            // time evolution
            // --------------

            if (ComponentMode) {
                MoveLevelSet(dt, NewLevelSet, LevelSetGrad, ExtVel, NEARgrid, marcher);
            }
            else {
                throw new NotImplementedException("todo");
            }


            // finally, apply the stabilization
            // --------------------------------
            //if (dt > 0) {

            //    if (ComponentMode) {

            //        // upwind stabilization can only be performed for vector extension velocity

            //        Debug.Assert(!object.ReferenceEquals(OldLevSet, NewLevelSet));

            //        var UpwdStabi = (new UpwindStabiForm()).Operator(2);
            //        UpwdStabi.Evaluate(-1.0 * dt, 1.0,
            //            OldLevSet.Mapping, ExtVel, NewLevelSet.Mapping,
            //            qInsEdge: new EdgeQuadratureScheme(domain: NEARgrid.InnerEdgesMask),
            //            qInsVol: new CellQuadratureScheme(domain: CellMask.GetEmptyMask(gdat)),
            //            bndMode: SpatialOperator.SubGridBoundaryModes.InnerEdge);



            //    } else {
            //        // lets interpret (\/Phi / |\/Phi|)*s as velocity?

            //        throw new NotImplementedException("todo");
            //    }

            //    var jp = new JumpPenalization();
            //    jp.ImplicitEuler(dt * 0.001, NEARgrid, NewLevelSet);
            //}



            //return x0;
        }

        private static void MoveLevelSet(double dt, SinglePhaseField LevelSet, VectorField<SinglePhaseField> LevelSetGrad, SinglePhaseField[] ExtVel, SubGrid NEARgrid, Reinit.FastMarch.FastMarchReinit marcher) {

            GridData gdat = (GridData)(LevelSet.GridDat);
            int D = gdat.SpatialDimension;

            //var TimeEvoOp = (new LevelSetEvoTerm_Vector()).Operator(2);


            SpatialOperator TimeEvoOp = new SpatialOperator(1, 2 * D, 1, QuadOrderFunc.NonLinear(2), "Phi", "dPhi_dx", "dPhi_dy", "Sx", "Sy", "c1");
            TimeEvoOp.EquationComponents["c1"].Add(new LevelSetEvoTerm_Vector());
            TimeEvoOp.EquationComponents["c1"].Add(new UpwindStabiForm());
            TimeEvoOp.Commit();

            RungeKutta RunschCjuda = new RungeKutta(RungeKuttaScheme.TVD3, TimeEvoOp,
                LevelSet.Mapping, new CoordinateMapping(ArrayTools.Cat<DGField>(LevelSetGrad, ExtVel)), sgrd: NEARgrid);

            RunschCjuda.OnBeforeComputeChangeRate += delegate (double AbsTime, double RelTime) {
                marcher.GradientUpdate(NEARgrid, LevelSet, LevelSetGrad);
            };

            double dt_CFL = gdat.ComputeCFLTime(ExtVel.ToArray(), 10000, NEARgrid.VolumeMask);
            dt_CFL *= 1.0 / (((double)(LevelSet.Basis.Degree)).Pow2());
            if (dt / dt_CFL >= 1.0) {
                Console.WriteLine(" Warning: exceeding Level-Set CFL: dt = {0:e4}, dt_CFL = {1:e4}, frac = {2:e4}", dt, dt_CFL, dt / dt_CFL);
                throw new ArithmeticException("Levelset CFL exceeded");
            }

            RunschCjuda.Perform(dt);

            LevelSet.CheckForNanOrInf(true, true, true);

            // stabilization term:
            // -------------------

            //var jp = new JumpPenalization();
            //jp.ImplicitEuler(dt * 0.001, NEARgrid, LevelSet); // how large should the parameter be?
        }

        class UpwindStabiForm : IEdgeForm {
            
            public TermActivationFlags BoundaryEdgeTerms {
                get { return TermActivationFlags.None; }
            }

            public TermActivationFlags InnerEdgeTerms {
                get { return TermActivationFlags.UxV | TermActivationFlags.V; }
            }

            public double InnerEdgeForm(ref CommonParams inp, double[] PhiIn, double[] PhiOt, double[,] GradPhiIn, double[,] GradPhiOt, double vIn, double vOt, double[] Grad_vIn, double[] Grad_vOt) {
                int D = inp.D;
                double Disc = 0, FluxIn = 0, FluxOt = 0;
                 
                for(int d = 0; d < D; d++) {
                    Disc = (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normale[d];
                    FluxIn += (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * 0.5 * inp.Normale[d] * PhiIn[0];
                    FluxOt += (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * 0.5 * inp.Normale[d] * PhiOt[0];
                }

                double Flux = 0;
                if(Disc >= 0) {
                    // use inner values
                    Flux = FluxIn;
                } else {
                    // use outer values
                    Flux = FluxOt;
                }

                return (Flux - FluxIn) * vIn - (Flux - FluxOt) * vOt;
            }

            public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
                return 0;
            }

            public IList<string> ArgumentOrdering {
                get { return new string[] { "Phi" }; }
            }

            public IList<string> ParameterOrdering {
                get { 
                    return new string[] { "Sx", "Sy" }; 
                }
            }
        }

        class LevelSetEvoTerm_Vector : IVolumeForm {

            public TermActivationFlags VolTerms {
                get { return TermActivationFlags.V | TermActivationFlags.UxV; }
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
                int D = cpv.D;
                Debug.Assert(U.Length % 2 == 0);
                double acc = 0;
                for(int i = 0; i < D; i++) {
                    acc += cpv.Parameters[i] * cpv.Parameters[i + D];
                }
                return acc*V;
            }

            public IList<string> ArgumentOrdering {
                get { 
                    return new string[0]; 
                }
            }

            public IList<string> ParameterOrdering {
                get { 
                    return new string[] { "dPhi_dx", "dPhi_dy", "Sx", "Sy" }; 
                }
            }
        }


        /*
        /// <summary>
        /// The term $ \vec{S} \cdot \nabla \varphi $, where $\vec{S}$ denotes the vector extension velocity.
        /// </summary>
        static double LevelSetEvoTerm_Vector(double[] X, double[] U, int jCell) {
            //ouble s = U[0];
            int D = U.Length / 2;
            Debug.Assert(U.Length % 2 == 0);
            double acc = 0;
            for(int i = 0; i < D; i++) {
                acc += U[i] * U[i + D];
            }
            return acc;
        }
        */

        const double RelNearThreshold = 0.4;

        private static BitArray FindReinitDomain(SinglePhaseField NewLevelSet, SubGrid CCgrid, SubGrid signGrid, int sign) {
            GridData gdat = (GridData)(NewLevelSet.GridDat);
            
            // find the level-set extremals
            // ============================

            //var signGrid = new SubGrid(signDomain);
            var signDomain = signGrid.VolumeMask;
            
            var signCritEdges = signGrid.BoundaryEdgesMask.Intersect(CCgrid.BoundaryEdgesMask);
            double[] CritEdges_minis, CritEdges_maxis;
            ExtremalsOnEdge(NewLevelSet, CCgrid.VolumeMask, signCritEdges, out CritEdges_minis, out CritEdges_maxis);

            if(sign < 0) {
                CritEdges_minis = CritEdges_maxis;
                CritEdges_minis.ScaleV(-1);
                CritEdges_maxis = null;
            } else {
                CritEdges_maxis = null;
            }
            

            // find if we have to do some reinit
            // =================================

            BitArray PosBitMask = signDomain.GetBitMask();
            BitArray CCbitmask = null;
            BitArray ExtensionBitMask = null;

            int J = gdat.Cells.NoOfLocalUpdatedCells;
            int[,] Edge2Cell = gdat.Edges.CellIndices;

            int[] Edges = signCritEdges.ItemEnum.ToArray();
            int I = Edges.Length;

            int NoOfFoundEdges = 0;

            var h_min = gdat.Cells.h_min;

            for(int i = 0; i < I; i++) {
                int iEdge = Edges[i];
                int jCell_IN = Edge2Cell[iEdge, 0];
                int jCell_OT = Edge2Cell[iEdge, 1];
                        
                double h; // threshold which triggers reinit if the levset is this close to the edge
                {
                    if(jCell_OT < 0)
                        h = h_min[jCell_IN];
                    else
                        h = Math.Min(h_min[jCell_IN], h_min[jCell_OT]);

                    h *= RelNearThreshold;
                }


                if(CritEdges_minis[i] <= h) {
                    NoOfFoundEdges++;

                    Console.WriteLine("LS moved: " + sign);

                    // we need the masks, so initialize them
                    if(ExtensionBitMask == null)
                        ExtensionBitMask = new BitArray(J);
                    if(CCbitmask == null)
                        CCbitmask = CCgrid.VolumeMask.GetBitMask();

                    // find the cell that has been left by the level-set:
                    // The neighbours of this cell must be re-initialized.
                    int CutCell = -1;
                    {
                        if(jCell_OT < 0) {
                            // level-set has moved out ouf bounds: we are done
                            //  
                            Debug.Assert(CCbitmask[jCell_IN] == true);
                            CutCell = -1;
                        } else if(CCbitmask[jCell_IN] == true) {
                            Debug.Assert(CCbitmask[jCell_OT] == false);
                            CutCell = jCell_IN;
                        } else if(CCbitmask[jCell_OT] == true) {
                            Debug.Assert(CCbitmask[jCell_IN] == false);
                            CutCell = jCell_OT;
                        } else {
                            Debug.Assert(false, "should never occur");
                        }
                    }

                    if(CutCell >= 0) {
                        int[] NeighCells, ConnectingEdges;
                        gdat.GetCellNeighbours(CutCell, GetCellNeighbours_Mode.ViaVertices, out NeighCells, out ConnectingEdges);

                        foreach(int jNeigh in NeighCells) {
                            // the neighbour cell 'jNeigh' must be reinitialized, if it is part of the (old) pos-domain (or neg-dom.)
                            ExtensionBitMask[jNeigh] = PosBitMask[jNeigh];
                        }
                    }
                }
            }

            //if(NoOfFoundEdges != 15 && NoOfFoundEdges != 0)
            //    Debugger.Break();

            return ExtensionBitMask;
        }

        /// <summary>
        /// Sets the level-set field values to +1 or -1 for all un-cut cells
        /// </summary>
        public static void NormalizeFar(LevelSetTracker Tracker, double val = 1) {
            if(val <= 0)
                throw new ArgumentOutOfRangeException();
            var CC = Tracker.Regions.GetCutCellMask4LevSet(iLevSet);
            var Pos = Tracker.Regions.GetLevelSetWing(iLevSet, +1).VolumeMask.Except(CC);
            var Neg = Tracker.Regions.GetLevelSetWing(iLevSet, -1).VolumeMask.Except(CC);
            
            LevelSet LevSet = (LevelSet)(Tracker.LevelSets[iLevSet]);
            
            LevSet.Clear(Pos);
            LevSet.AccConstant(val, Pos);
            
            LevSet.Clear(Neg);
            LevSet.AccConstant(-val, Neg);
        }


        /// <summary>
        /// Construction of extension velocity (vector and scalar) on the cut cells by geometric means.
        /// </summary>
        /// <param name="Velocity">Input: Velocity Vector field</param>
        /// <param name="ExtVel">output: constructed extension velocity</param>
        /// <param name="ExtVelMax">output: maximum values for extension velocity for each cut cell</param>
        /// <param name="ExtVelMin">output: minimum values for extension velocity for each cut cell</param>
        /// <param name="Tracker"></param>
        /// <param name="Levset">Input</param>
        /// <param name="LevSetGrad">Output/Work</param>
        public static MultidimensionalArray ConstructExtVel_GOM_CC(LevelSetTracker Tracker, SinglePhaseField[] ExtVel, ConventionalDGField[] Velocity,
            SinglePhaseField Levset, VectorField<SinglePhaseField> LevSetGrad,
            double[][] ExtVelMin, double[][] ExtVelMax) {
            GridData grdDat = Tracker.GridDat;
            int D = grdDat.SpatialDimension;
            int pDG = Levset.Basis.Degree;
            RefElement[] Krefs = grdDat.Grid.RefElements;
            
            bool ComponentMode; // whether scalar or vector extension velocity is used
            if(ExtVel.Length == Velocity.Length) {
                // 
                ComponentMode = true;

            } else {
                ComponentMode = false;
                if(ExtVel.Length != 1)
                    throw new ArgumentException();
                if(Velocity.Length != D)
                    throw new ArgumentException();
            }

            if(ExtVelMin.GetLength(0) != ExtVel.Length)
                throw new ArgumentException();
            if(ExtVelMax.GetLength(0) != ExtVel.Length)
                throw new ArgumentException();


            // ================================
            // find closest points on level-set
            // ================================

            QuadRule[] qrs = Krefs.Select(Kref => Kref.GetQuadratureRule(pDG * 2)).ToArray();
            NodeSet[] qrs_nodes = qrs.Select(qr => qr.Nodes).ToArray();

            SubGrid cutCellSgrd = Tracker.Regions.GetNearFieldSubgrid4LevSet(iLevSet, 0);
            
            ClosestPointFinder cp = new ClosestPointFinder(Tracker, 0, cutCellSgrd, _Nodes: qrs_nodes);
            
            
            // ==============================
            // evaluate velocity and gradient
            // ==============================

            MultidimensionalArray[] VelocityAtLevSet = new MultidimensionalArray[Velocity.Length];
            MultidimensionalArray[] LevelSetGradient = new MultidimensionalArray[D];
            if(!ComponentMode) {
                for(int d = 0; d < D; d++) {
                    LevSetGrad[d].Clear(Tracker.Regions.GetCutCellMask4LevSet(iLevSet));
                    LevSetGrad[d].Derivative(1.0, Levset, d, Tracker.Regions.GetCutCellMask4LevSet(iLevSet));

                    LevelSetGradient[d] = cp.EvaluateAtCp(LevSetGrad[d]);
                }
            }
            for(int d = 0; d < Velocity.Length; d++) {
                VelocityAtLevSet[d] = cp.EvaluateAtCp(Velocity[d]);
            }
            
            // =============================================
            // construct extension velocity on the cut cells
            // =============================================

            CellMask cutCellMask = cutCellSgrd.VolumeMask;
            foreach(int jCell in cutCellMask.ItemEnum) {
                for(int d = 0; d < ExtVel.Length; d++) {
                    ExtVelMin[d][jCell] = double.MaxValue;
                    ExtVelMax[d][jCell] = -double.MaxValue;
                }
            }
            
            CellQuadratureScheme cqs = (new CellQuadratureScheme(false, cutCellMask));
            for(int iKref = 0; iKref < Krefs.Length; iKref++) {
                cqs = cqs.AddFixedRule(qrs[iKref], grdDat.Cells.GetCells4Refelement(iKref));
            }

            int[] jCell2jSub = cutCellSgrd.LocalCellIndex2SubgridIndex;

            if(ComponentMode) {
                for(int d = 0; d < ExtVel.Length; d++) {
                    ExtVel[d].Clear(cutCellMask);
                    ExtVel[d].ProjectField(1.0,
                        (ScalarFunctionEx)delegate(int j0, int Len, NodeSet nodes, MultidimensionalArray result) {
                            for(int j = 0; j < Len; j++) {
                                int jCell = j + j0;
                                int jSub = jCell2jSub[jCell];

                                for(int k = 0; k < nodes.NoOfNodes; k++) {
                                    double S = VelocityAtLevSet[d][jSub, k];
                                    result[j, k] = S;
                                }
                            }
                        },
                        cqs);
                }
            } else {

                ExtVel[0].Clear(cutCellMask);
                ExtVel[0].ProjectField(1.0,
                    (ScalarFunctionEx)delegate(int j0, int Len, NodeSet nodes, MultidimensionalArray result) {
                        for(int j = 0; j < Len; j++) {
                            int jCell = j + j0;
                            int jSub = jCell2jSub[jCell];

                            for(int k = 0; k < nodes.NoOfNodes; k++) {

                                double norm = 0;
                                double S = 0;
                                for(int d = 0; d < D; d++) {
                                    double Nd = LevelSetGradient[d][jSub, k];
                                    double Ud = VelocityAtLevSet[d][jSub, k];

                                    norm += Nd * Nd;
                                    S += Ud * Nd;
                                }

                                norm = 1.0 / Math.Sqrt(norm);
                                S = S * norm;

                                ExtVelMin[0][jCell] = Math.Min(ExtVelMin[0][jCell], S);
                                ExtVelMax[0][jCell] = Math.Max(ExtVelMax[0][jCell], S);

                                result[j, k] = S;
                            }
                        }
                    },
                    cqs);
            }

            // ========================
            // minium/maximum detection
            // ========================
            
            // Minima und Maxima für jede Zelle sollten für die exakte Lösung niemals
            // die Werte am Level-Set überschreiten.
            // Wir akzeptieren das trotzdem.
            
            
            foreach(int j in cutCellMask.ItemEnum) {
                for(int d = 0; d < ExtVel.Length; d++) {
                    double miniCell, maxiCell;
                    ExtVel[d].GetExtremalValuesInCell(out miniCell, out maxiCell, j);

                    ExtVelMin[d][j] = Math.Min(ExtVelMin[d][j], miniCell);
                    ExtVelMax[d][j] = Math.Max(ExtVelMax[d][j], maxiCell);
                }

                
            }


            // return the closest points
            // -------------------------

            return cp.X0_global.CloneAs();
        }

        /// <summary>
        /// Construction of extension velocity (vector and scalar) on the cut cells by means of PDE.
        /// </summary>
        /// <param name="Velocity">Input: Velocity Vector field</param>
        /// <param name="ExtVel">output: constructed extension velocity</param>
        /// <param name="ExtVelMax">output: maximum values for extension velocity for each cut cell
        /// 1st index: correlates with index into <paramref name="ExtVel"/>
        /// 2nd index: local cell index.
        /// </param>
        /// <param name="ExtVelMin">output: minimum values for extension velocity for each cut cell</param>
        /// <param name="Tracker"></param>
        /// <param name="Levset"></param>
        /// <param name="LevSetGrad"></param>
        /// <param name="HMForder"></param>
        public static void ConstructExtVel_PDE(LevelSetTracker Tracker, SubGrid subGrid, SinglePhaseField[] ExtVel, ConventionalDGField[] Velocity,
            SinglePhaseField Levset, VectorField<SinglePhaseField> LevSetGrad,
            double[][] ExtVelMin, double[][] ExtVelMax,
            int HMForder) //
        {
            GridData grdDat = Tracker.GridDat;
            int D = grdDat.SpatialDimension;
            LevelSet LevSet0 = ((LevelSet)(Tracker.LevelSets[iLevSet])); // Level-Set at initial time
            int pDG = LevSet0.Basis.Degree;
            RefElement[] Krefs = grdDat.Grid.RefElements;


            bool ComponentMode;
            if(ExtVel.Length == Velocity.Length) {
                // 
                ComponentMode = true;

            } else {
                ComponentMode = false;
                if(ExtVel.Length != 1)
                    throw new ArgumentException();
                if(Velocity.Length != D)
                    throw new ArgumentException();
            }

            if(ExtVelMin.GetLength(0) != ExtVel.Length)
                throw new ArgumentException();
            if(ExtVelMax.GetLength(0) != ExtVel.Length)
                throw new ArgumentException();


            // ===========================================================
            // solve extension velocity problem on cut-cells (all at once)
            // ===========================================================

            var subMask = subGrid.VolumeMask;

            double penaltyBase = ((double)(ExtVel[0].Basis.Degree)).Pow2();

            // matrix assembly
            // ---------------

            MsrMatrix ExtVelMatrix;
            double[][] ExtVelRHS = new double[ComponentMode ? ExtVel.Length : 1][];
            List<int> UsedRows = new List<int>();
            {
                var map = ExtVel[0].Mapping;

                // mem alloc
                ExtVelMatrix = new MsrMatrix(map, map);
                for(int d = 0; d < ExtVelRHS.Length; d++) {
                    ExtVelRHS[d] = new double[map.LocalLength];
                }

                // generate operators
                ILevelSetComponent InterfaceFlux;
                if (ComponentMode)
                    InterfaceFlux = new EllipticExtension.SingleComponentInterfaceForm(+penaltyBase, Tracker);
                else
                    InterfaceFlux = new EllipticExtension.ScalarVelocityInterfaceForm(+penaltyBase, Tracker);

                XSpatialOperator InterfaceOperator = InterfaceFlux.XOperator((int[] A, int[] B, int[] C) => HMForder);

                var BulkForm = new EllipticExtension.ExtVelForm_bulk(penaltyBase, 0.0 ,InterfaceFlux,Tracker, subMask.GetBitMaskWithExternal());

                var BulkOperator = BulkForm.Operator();

                VectorField<SinglePhaseField> meanLevSetGrad = new VectorField<SinglePhaseField>(D, new Basis(grdDat, 0), SinglePhaseField.Factory);
                meanLevSetGrad.AccLaidBack(1.0, LevSetGrad, subMask);
                

                // quadrature schemes
                //var SchHelper = new XQuadSchemeHelper(Tracker, HMFversion, Tracker.SpeciesIdS.ToArray());
                //CellQuadratureScheme surfScheme = SchHelper.GetLevelSetquadScheme(0, subMask);
                var SchHelper = Tracker.GetXDGSpaceMetrics(Tracker.SpeciesIdS.ToArray(), HMForder, 1).XQuadSchemeHelper;
                
                EdgeQuadratureScheme emptyEdgeScheme = new EdgeQuadratureScheme(domain: EdgeMask.GetEmptyMask(grdDat));

                EdgeQuadratureScheme posEdgesScheme, negEdgesScheme;
                posEdgesScheme = SchHelper.GetEdgeQuadScheme(Tracker.GetSpeciesId("B"), UseDefaultFactories: false, IntegrationDomain: subGrid.InnerEdgesMask, fixedOrder: HMForder);
                negEdgesScheme = SchHelper.GetEdgeQuadScheme(Tracker.GetSpeciesId("A"), UseDefaultFactories: false, IntegrationDomain: subGrid.InnerEdgesMask, fixedOrder: HMForder);

                // integrate interface
                // (Only the interface part contributes to the RHS; mutst be integrated fore EACH component.)
                DGField[] IfParams;
                if (ComponentMode) {
                    IfParams = new DGField[] { Velocity[0] };
                }
                else {
                    IfParams = ArrayTools.Cat<DGField>(Velocity);
                }
                for (int d = 0; d < ExtVelRHS.Length; d++) {
                    if (ComponentMode) {
                        IfParams = new DGField[] { Velocity[d] };
                    }

                    InterfaceOperator.ComputeMatrixEx(
                        lsTrk: Tracker,
                        DomainMap: map,
                        Parameters: IfParams,
                        CodomainMap: map,
                        Matrix: (d == 0 ? ExtVelMatrix : default(MsrMatrix)),
                        AffineOffset: ExtVelRHS[d],
                        OnlyAffine: d != 0,
                        time: 0,
                        MPIParameterExchange:false,
                        whichSpc:Tracker.GetSpeciesId("A"));

                    ExtVelRHS[d].ScaleV(-1.0);

                }

                // integrate bulk, volume part
                // (The linear part is equal for ALL components, integration is required only once.)
                
                
                double[] dummy_affine = new double[map.LocalLength];
                //BulkForm.phi_sign = double.NaN; // sollte keine Auswirkungen haben, da im volumenteil nicht vorkommend!
                DGField[] BulkParams = ArrayTools.Cat<DGField>(LevSetGrad, Levset, meanLevSetGrad, IfParams.ToArray());
                BulkOperator.ComputeMatrixEx(
                    map, BulkParams, map,
                    ExtVelMatrix, dummy_affine,
                    edgeQuadScheme: emptyEdgeScheme,
                    volQuadScheme: new CellQuadratureScheme(UseDefaultFactories: true, domain: subMask));
                Debug.Assert(dummy_affine.L2NormPow2() == 0);

                // integrate bulk, edges
                // The information flow direction on edges depends on the sign of the level-set.
                // Therefore, to be "more consistent", whatever that means, we use the HMF edge rules.
                // (In the worst case, this is just costly, but it should not harm.)
                for(int i = 0; i < 2; i++) {
                    EdgeQuadratureScheme edgSch;
                    switch(i) {
                        case 0:
                        //BulkForm.phi_sign = +1;
                        edgSch = posEdgesScheme;
                        break;
                        case 1:
                        //BulkForm.phi_sign = -1;
                        edgSch = negEdgesScheme;
                        break;
                        default: throw new Exception("Should never occur.");
                    }

                    BulkOperator.ComputeMatrixEx(
                        map, BulkParams, map,
                        ExtVelMatrix, dummy_affine,
                        edgeQuadScheme: edgSch,
                        volQuadScheme: new CellQuadratureScheme(domain: subMask));
                    Debug.Assert(dummy_affine.L2NormPow2() == 0);
                }


                // test matrix and affine vectors
#if DEBUG
                int J = grdDat.Cells.NoOfLocalUpdatedCells;
                BitArray CCbitmask = subMask.GetBitMaskWithExternal();
                for(int j = 0; j < J; j++) {
                    int N = map.BasisS[0].GetLength(j);
                    for(int n = 0; n < N; n++) {
                        int iL = map.LocalUniqueCoordinateIndex(0, j, n);
                        int iG = map.GlobalUniqueCoordinateIndex(0, j, n);

                        if(CCbitmask[j]) {
                            //var Row = ExtVelMatrix.GetRow(iG);
                            //foreach (var entry in Row) {
                            int Lr;
                            int[] row_cols = null;
                            double[] row_vals = null;
                            Lr = ExtVelMatrix.GetRow(iG, ref row_cols, ref row_vals);
                            for (int lr = 0; lr < Lr; lr++) {
                                int ColIndex = row_cols[lr];
                                double Value = row_vals[lr];
                                if (Value != 0.0) {
                                    int jL = map.Global2LocalIndex(ColIndex);
                                    int iFld, jCellCol, m;
                                    map.LocalFieldCoordinateIndex(jL, out iFld, out jCellCol, out m);
                                    Debug.Assert(CCbitmask[jCellCol] == true);
                                    Debug.Assert(iFld == 0);
                                }
                            }

                        } else {
                            Debug.Assert(ExtVelMatrix.GetNoOfNonZerosPerRow(iG) == 0);
                            for(int d = 0; d < D; d++) {
                                Debug.Assert(ExtVelRHS[d][iL] == 0);
                            }
                        }
                    }
                }
#endif

                int i0 = map.i0;
                int iE = map.iE;
                for(int i = i0; i < iE; i++) {
                    if(ExtVelMatrix.GetNoOfNonZerosPerRow(i) == 0) {
                        ExtVelMatrix[i, i] = 1.0;
                    } else {
                        UsedRows.Add(i);
                    }
                }

            }

            // solving
            // -------
            MsrMatrix essExtVelMatrix;
            {
                essExtVelMatrix = new MsrMatrix(new Partitioning(UsedRows.Count));
                ExtVelMatrix.WriteSubMatrixTo(essExtVelMatrix, UsedRows, default(int[]), UsedRows, default(int[]));
            }


            // using a sparse direct solver
            //////////////////////////////////////////
            

            using(var slv = new ilPSP.LinSolvers.PARDISO.PARDISOSolver()) {
                slv.CacheFactorization = true;
                slv.DefineMatrix(essExtVelMatrix);

                Debug.Assert(ExtVelRHS.Length == ExtVel.Length);

                for(int d = 0; d < ExtVelRHS.Length; d++) {
                    // extract rhs
                    double[] rhs = new double[UsedRows.Count];
                    double[] x = new double[UsedRows.Count];
                    rhs.AccV(1.0, ExtVelRHS[d], default(int[]), UsedRows);

                    double[] Resi_d = rhs.CloneAs();
                    double NORM_rhs = rhs.L2NormPow2().MPISum().Sqrt();
                    
                    // solve
                    slv.Solve(x, rhs);
                    
                    // check the solution; never trust PARDISO!
                    essExtVelMatrix.SpMVpara(-1.0, x, 1.0, Resi_d);
                    double NORM_Resi_d = Resi_d.L2NormPow2().MPISum().Sqrt();
                    double RelResiNorm = NORM_Resi_d/NORM_rhs;

                    if(RelResiNorm > 1.0e-8) {
                        Console.WriteLine("WARNING: high solver residual in extension velocity solver: relative residual norm: {0:E2}", RelResiNorm);
                    }

                    // write back solution
                    ExtVel[d].Clear();
                    ExtVel[d].CoordinateVector.AccV(1.0, x, UsedRows, default(int[]));
                }
            }
            
            /*
              
            // using a dense solver from LAPACK
            //////////////////////////////////////////
            
            if(UsedRows.Count > 4000) {
                Console.WriteLine("Warning: extension velocity, dense direct solver on {0}x{0} -- matrix.", UsedRows.Count);
            }
            
            MultidimensionalArray MM = essExtVelMatrix.ToFullMatrixOnProc0();
            for(int d = 0; d < ExtVelRHS.Length; d++) {
                double[] rhs = new double[UsedRows.Count];
                double[] x = new double[UsedRows.Count];
                rhs.AccV(1.0, ExtVelRHS[d], default(int[]), UsedRows);
                
                MM.Solve(x, rhs);

                ExtVel[d].CoordinateVector.ClearEntries();
                ExtVel[d].CoordinateVector.AccV(1.0, x, UsedRows, default(int[]));
            }
            */

            // ========================
            // minium/maximum detection
            // ========================

            // Eigentlich bräuchten wir Min/Max am Interface, aber dafür bräuchten wir wieder den ClosestPointfinder

            foreach(int j in subMask.ItemEnum) {
                for(int d = 0; d < ExtVel.Length; d++) {
                    double miniCell, maxiCell;
                    ExtVel[d].GetExtremalValuesInCell(out miniCell, out maxiCell, j);

                    ExtVelMin[d][j] = miniCell;
                    ExtVelMax[d][j] = maxiCell;
                }
            }

            
        }

        
 
        /// <summary>
        /// Finds the extremals of the level-set field on the boundary between cut- and near cells:
        /// we need this to detect wether the level-set has left the cell.
        /// </summary>
        public static void ExtremalsOnEdge(ConventionalDGField LevSet, CellMask cutCellsMask, EdgeMask em, out double[] minis, out double[] maxis) {
            GridData gdat = (GridData)(LevSet.GridDat);
            minis = new double[em.NoOfItemsLocally];
            maxis = new double[em.NoOfItemsLocally];
            BitArray cutCellsBitmask = cutCellsMask.GetBitMask();

#if DEBUG
            {
                SubGrid cc = new SubGrid(cutCellsMask);
                EdgeMask defect = em.Except(cc.BoundaryEdgesMask);
                if(defect.NoOfItemsLocally > 0)
                    throw new ArgumentException("The given edge mask must be on the boundary of the cut cells.");
            }
#endif

            int[,] Edge2Cell = gdat.Edges.CellIndices;
            int[,] TrafoIdx = gdat.Edges.Edge2CellTrafoIndex;

            NodeSet[] TestNodes = new NodeSet[gdat.Edges.EdgeRefElements.Length];

            int i = -1;
            foreach(int iEdge in em.ItemEnum) {
                i++;

                // find cell and edge-to-cell trafo index
                // --------------------------------------

                int jCell, iTrafo;
                {
                    int jCellIn = Edge2Cell[iEdge, 0];
                    int jCellOt = Edge2Cell[iEdge, 1];
                    if(cutCellsBitmask[jCellIn] == true) {
                        jCell = jCellIn;
                        iTrafo = TrafoIdx[iEdge, 0];
                        if(jCellOt >= 0 && cutCellsBitmask[jCellOt] == true)
                            throw new ApplicationException("The given edge mask must be on the boundary of the cut cells.");

                    } else if(jCellOt >= 0 && cutCellsBitmask[jCellOt] == true) {
                        jCell = jCellOt;
                        iTrafo = TrafoIdx[iEdge, 1];
                        if(cutCellsBitmask[jCellIn] == true)
                            throw new ArgumentException("The given edge mask must be on the boundary of the cut cells.");
                    } else {
                        throw new ArgumentException("The given edge mask must be on the boundary of the cut cells.");
                    }
                }

                // get the test nodes at the edge
                // ------------------------------

                int iKrefEdge = gdat.Edges.GetRefElementIndex(iEdge);
                if(TestNodes[iKrefEdge] == null) {
                    Debug.Assert(gdat.Edges.EdgeRefElements[iKrefEdge].SpatialDimension == 1, "Number of thest vertices maybe is high.");
                    TestNodes[iKrefEdge] = gdat.Edges.EdgeRefElements[iKrefEdge].GetSubdivisionTree((LevSet.Basis.Degree + 1) * 2).GlobalVertice;
                }

                // evaluate the level-set field at the edge, find extremals
                // --------------------------------------------------------

                int NoOfNodes = TestNodes[iKrefEdge].NoOfNodes;
                MultidimensionalArray LevelSetValAtEdge = MultidimensionalArray.Create(1, NoOfNodes);
                LevSet.Evaluate(jCell, 1, TestNodes[iKrefEdge].GetVolumeNodeSet(gdat, iTrafo), LevelSetValAtEdge);

                minis[i] = LevelSetValAtEdge.Min();
                maxis[i] = LevelSetValAtEdge.Max();
            }

        }

        /// <summary>
        /// Performs Level-Set ReInitialization on the Cut-Cells
        /// by searching the closest point to each quadrature point
        /// </summary>
        /// <param name="lsTrk"></param>
        /// <param name="NewLevelSet">Return Level-Set</param>
        public static void CutCellReinit(LevelSetTracker lsTrk, SinglePhaseField NewLevelSet) {
            GridData grdDat = lsTrk.GridDat;
            int D = grdDat.SpatialDimension;

            LevelSet LevSet0 = ((LevelSet)(lsTrk.LevelSets[iLevSet])); // Level-Set at initial time
            RefElement[] Krefs = grdDat.Grid.RefElements;

            Basis LevSetBasis = NewLevelSet.Basis;
            int pDG = LevSetBasis.Degree;


            // ================================
            // find closest points on level-set
            // ================================

            QuadRule[] qrs = Krefs.Select(Kref => Kref.GetQuadratureRule(pDG * 2)).ToArray();
            NodeSet[] qrs_nodes = qrs.Select(qr => qr.Nodes).ToArray();

            SubGrid cutCellSgrd = lsTrk.Regions.GetNearFieldSubgrid4LevSet(iLevSet, 0);

            ClosestPointFinder cp = new ClosestPointFinder(lsTrk, iLevSet, cutCellSgrd, _Nodes: qrs_nodes);

            // ==================================
            // get level-set sign at quad. nodes
            // ==================================

            CellMask cutFieldMask = cutCellSgrd.VolumeMask;
            int[] jCell2jSub = cutCellSgrd.LocalCellIndex2SubgridIndex;

            MultidimensionalArray LevSet0_Vals = MultidimensionalArray.Create(cutCellSgrd.LocalNoOfCells, qrs_nodes.Max(nds => nds.NoOfNodes));
            int jSub = 0;
            foreach(Chunk chk in cutCellSgrd.VolumeMask) {
                NodeSet nds = qrs_nodes[grdDat.Cells.GetRefElementIndex(chk.i0)];

                LevSet0.Evaluate(chk.i0, chk.Len, nds, LevSet0_Vals.ExtractSubArrayShallow(new int[] { jSub, 0 }, new int[] { jSub + chk.Len - 1, nds.NoOfNodes - 1 }));

                jSub += chk.Len;
            }


            // ==========================================================
            // project closest-point-distance onto new level set
            // ==========================================================

            
            CellQuadratureScheme cqs = (new CellQuadratureScheme(false, cutFieldMask));
            for(int iKref = 0; iKref < Krefs.Length; iKref++) {
                cqs = cqs.AddFixedRule(qrs[iKref], grdDat.Cells.GetCells4Refelement(iKref));
            }

            MultidimensionalArray X = cp.X_global;
            MultidimensionalArray X0 = cp.X0_global_Resorted;


            NewLevelSet.Clear(null);
            NewLevelSet.ProjectField(1.0,
                (ScalarFunctionEx)delegate(int j0, int Len, NodeSet nodes, MultidimensionalArray result) //
                {
                    for(int j = 0; j < Len; j++) {
                        int jCell = j + j0;
                        int _jSub = jCell2jSub[jCell];

                        for(int k = 0; k < nodes.NoOfNodes; k++) {

                            double acc = 0;
                            for(int d = 0; d < D; d++) {
                                acc += (X0[_jSub, k, d] - X[_jSub, k, d]).Pow2();
                            }
                            acc = Math.Sqrt(acc);

                            double sign = Math.Sign(LevSet0_Vals[_jSub, k]);

                            result[j, k] = acc*sign;
                        }
                    }
                },
                cqs);

            // ================================
            // set positive and negative region
            // ================================

            var CC = lsTrk.Regions.GetCutCellMask4LevSet(iLevSet);
            var Pos = lsTrk.Regions.GetLevelSetWing(iLevSet, +1).VolumeMask.Except(CC);
            var Neg = lsTrk.Regions.GetLevelSetWing(iLevSet, -1).VolumeMask.Except(CC);


            NewLevelSet.Clear(Pos);
            NewLevelSet.AccConstant(1, Pos);

            NewLevelSet.Clear(Neg);
            NewLevelSet.AccConstant(-1, Neg);

        }

      /*
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Velocity">
        /// Input; velocity field to compute the extension velocity from
        /// </param>
        /// <param name="NewLevelSet">
        /// output;
        /// </param>
        /// <param name="Tracker"></param>
        public static void Evolve_Lagrange(VectorField<SinglePhaseField> Velocity, SinglePhaseField NewLevelSet, LevelSetTracker Tracker, double dt) {
            GridData grdDat = Tracker.GridDat;
            int D = grdDat.SpatialDimension;
            LevelSet LevSet0 = ((LevelSet)(Tracker.LevelSets[iLevSet])); // Level-Set at initial time
            int pDG = LevSet0.Basis.Degree;
            RefElement[] Krefs = grdDat.Grid.RefElements;

            if(!NewLevelSet.Basis.Equals(LevSet0.Basis))
                throw new ArgumentException();
            Basis LevSetBasis = NewLevelSet.Basis;


            // ================================
            // find closest points on level-set
            // ================================

            NodeSet[] UnitedNodes = new NodeSet[Krefs.Length];
            {
                NodeSet[] VolumeNodes = Krefs.Select(Kref => Kref.GetQuadratureRule(pDG * 4).Nodes).ToArray();
                //NodeSet[] Vertices = Krefs.Select(Kref => Kref.Vertices).ToArray();

                for(int iKref = 0; iKref < Krefs.Length; iKref++) {

                    int N_vol = VolumeNodes[iKref].GetLength(0);
                    //int N_vtx = Vertices[iKref].GetLength(0);
                    int N_vtx = 0;

                    UnitedNodes[iKref] = new NodeSet(grdDat.Cells.RefElements[iKref], N_vol + N_vtx, D);
                    UnitedNodes[iKref].SetSubArray(VolumeNodes[iKref], new int[] { 0, 0 }, new int[] { N_vol - 1, D - 1 });
                    //UnitedNodes[iKref].SetSubArray(Vertices[iKref], new int[] { N_vol, 0 }, new int[] { N_vol + N_vtx - 1, D - 1 });
                    UnitedNodes[iKref].LockForever();
                }
            }

            SubGrid oldCutCellSgrd = Tracker.GetNearFieldSubgrid4LevSet(0, 0);

            ClosestPointFinder cp = new ClosestPointFinder(Tracker, 0, oldCutCellSgrd, _Nodes: UnitedNodes);

            
            // ==============================
            // evaluate velocity and gradient
            // ==============================

            MultidimensionalArray[] VelocityAtLevSet = new MultidimensionalArray[D];  // velocity at closest points
            MultidimensionalArray[] LevelSetGradient = new MultidimensionalArray[D];  // level-set gradient at closest points
            SinglePhaseField lsGrad = new SinglePhaseField(LevSet0.Basis);
            for(int d = 0; d < D; d++) {
                lsGrad.Clear(Tracker.GetCutCellMask4LevSet(iLevSet));
                lsGrad.Derivative(1.0, LevSet0, d, Tracker.GetCutCellMask4LevSet(iLevSet));

                VelocityAtLevSet[d] = cp.EvaluateAtCp(Velocity[d]);
                LevelSetGradient[d] = cp.EvaluateAtCp(lsGrad);
            }

            // ==============================================
            // compute extension velocity and translate nodes
            // ==============================================
            int JSUB_oldCut = oldCutCellSgrd.LocalNoOfCells;
            int K = VelocityAtLevSet[0].GetLength(1); // number of nodes
            Debug.Assert(JSUB_oldCut == VelocityAtLevSet[0].GetLength(0));
            Debug.Assert(JSUB_oldCut == VelocityAtLevSet[0].GetLength(0));
            Debug.Assert(K == UnitedNodes[0].NoOfNodes);

            MultidimensionalArray extVel = MultidimensionalArray.Create(JSUB_oldCut, K);
            Debug.Assert(cp.X_global.GetLength(0) == JSUB_oldCut);
            Debug.Assert(cp.X_global.GetLength(1) == K);
            
            MultidimensionalArray NewGlobalNodes = cp.X_global.ResizeShallow(JSUB_oldCut * K, D).CloneAs();
            ClosestPointFinder.PlottiPunkti(NewGlobalNodes, "OldNodes.csv");



            double[] VectorVel = new double[D];
            for(int iPoint = 0; iPoint < NewGlobalNodes.GetLength(0); iPoint++) {
                int jSub = cp.CellIndexMap[iPoint];
                int k = cp.NodeIndexMap[iPoint];

                double acc = 0;
                double GradNorm = 0;

                for(int d = 0; d < D; d++) {
                    VectorVel[d] = LevelSetGradient[d][jSub, k];
                    acc += VectorVel[d] * VelocityAtLevSet[d][jSub, k];
                    GradNorm += VectorVel[d].Pow2();
                }

                GradNorm = 1.0 / Math.Sqrt(GradNorm);
                double _extVel = acc * GradNorm;
                extVel[jSub, k] = _extVel;

                for(int d = 0; d < D; d++) {
                    VectorVel[d] *= GradNorm;
                    VectorVel[d] *= _extVel;

                    NewGlobalNodes[iPoint, d] += VectorVel[d] * dt;
                }
            }

            ClosestPointFinder.PlottiPunkti(NewGlobalNodes, "NewNodes.csv");

            // assign moved nodes to new cells
            //int[] CellIndexMap = cp.CellIndexMap.CloneAs(); // mapping: node index of X0 among all nodes ----> cell index of original point X    (X0 is closest point for X)
            //int[] NodeIndexMap = cp.NodeIndexMap.CloneAs(); // mapping: node index of X0 among all nodes ----> node index of original point X
            int[] I = cp.x0I2SgrdCell.CloneAs();
            int[] Perm, jSub2jCell_newCut;
            AssignPoints(oldCutCellSgrd, ref NewGlobalNodes, ref I, out Perm, out jSub2jCell_newCut);


            // ==================================
            // get level-set values at old points
            // ==================================

            MultidimensionalArray LevelSetValues = MultidimensionalArray.Create(JSUB_oldCut, UnitedNodes.Max(ns => ns.NoOfNodes));
            int[] jSub2jCell_oldCut = oldCutCellSgrd.SubgridIndex2LocalCellIndex;
            for(int jSub = 0; jSub < JSUB_oldCut; jSub++) {
                int jCell = jSub2jCell_oldCut[jSub];
                int iKref = grdDat.Cells.GetRefElementIndex(jCell);

                LevSet0.Evaluate(jCell, 1, UnitedNodes[iKref], LevelSetValues.ExtractSubArrayShallow(new int[] { jSub, 0}, new int[] {jSub, UnitedNodes[iKref].NoOfNodes - 1}));
            }

            
            // ==========================
            // re-construct the level-set
            // ==========================

            int[] jCell2jSub_oldCut = oldCutCellSgrd.LocalCellIndex2SubgridIndex;

            int[] InvPerm = new int[Perm.Length];
            for(int i = 0; i < Perm.Length; i++) {
                InvPerm[Perm[i]] = i;
            }


            NewLevelSet.Clear();

            for(int jsub = 0; jsub < (I.Length - 1); jsub++) {
                int I0 = I[jsub];
                int IE = I[jsub + 1];
                int NoOfNodes = IE - I0;
                int jCell = jSub2jCell_newCut[jsub];

                // transform points to local coordinates in 'jCell'
                // ------------------------------------------------
                NodeSet X0_loc = new NodeSet(grdDat.Cells.GetRefElement(jCell), NoOfNodes, D);
                MultidimensionalArray X0_glb = NewGlobalNodes.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });
                grdDat.TransformGlobal2Local(X0_glb, X0_loc, jCell); // transform to local
                X0_loc.LockForever();

                // ansatz for the polynomial reconstruction
                // ----------------------------------------
                MultidimensionalArray basisVals = LevSetBasis.CellEval(X0_loc, jCell, 1).ExtractSubArrayShallow(0, -1, -1);
                MultidimensionalArray basisGrad = LevSetBasis.CellEvalGradient(X0_loc, jCell, 1).ExtractSubArrayShallow(0, -1, -1, -1);
                int N = LevSetBasis.GetLength(jCell);
                MultidimensionalArray Matrix = MultidimensionalArray.Create((1) * NoOfNodes, N);
                //MultidimensionalArray RHS = MultidimensionalArray.Create((1 + D) * NoOfNodes);
                double[] RHS = new double[((1) * NoOfNodes)];

                for(int k = 0; k < NoOfNodes; k++) { // loop over matrix rows/points ...
                    int iPoint = I0 + k;

                    double ptNeu_x0 = X0_glb[k, 0];
                    double ptNeu_x1 = X0_glb[k, 1];

                    

                    int iOldPoint = Perm[iPoint];
                    int jSub_old = cp.CellIndexMap[iOldPoint];
                    int iPt_oldPoint = cp.NodeIndexMap[iOldPoint];


                    double ptAlt_x0 = cp.X_global[jSub_old, iPt_oldPoint, 0];
                    double ptAlt_x1 = cp.X_global[jSub_old, iPt_oldPoint, 1];


                    for(int n = 0; n < N; n++) { // loop over matrix columns/test functions ..
                        Matrix[k, n] = basisVals[k, n];
                    }
                    RHS[k] = cp.X_global[jSub_old, iPt_oldPoint, 0];
                    //RHS[k] = LevelSetValues[jSub_old, iPt_oldPoint];
                }

                
                
                //for(int k = 0; k < NoOfNodes; k++) { // loop over matrix rows/points ...
                //    //double[] point = X0_glb.GetRow(k);
                //    //double[] Normale = point.Normalize();
                //    double[] Normale = new double[] { 1, 0 };
                //    Normale = Normale.Normalize();

                //    for(int d = 0; d < D; d++) {
                //        int iPoint = I0 + k;
                //        int iOrgPoint = Perm[iPoint];
                //        int jCell_orgPoint = cp.CellIndexMap[iOrgPoint];
                //        int iPt_orgpoint = cp.NodeIndexMap[iOrgPoint];

                //        int iRow = k + (d + 1) * NoOfNodes;

                //        for(int n = 0; n < N; n++) { // loop over matrix columns/test functions ..
                //            Matrix[iRow, n] = basisGrad[k, n, d];
                //        }

                //        RHS[iRow] = Normale[d];// LevelSetGradient[d][jCell_orgPoint, iPt_orgpoint];
                //    }
                //}
                 

                // solve system and store data
                // ---------------------------
                double[] phiCoord = new double[N];
                Matrix.LeastSquareSolve(phiCoord, RHS);
                for(int n = 0; n < N; n++) {
                    NewLevelSet.Coordinates[jCell, n] = phiCoord[n];
                }
            }

        }
        */


        static void AssignPoints(SubGrid sgrd, ref MultidimensionalArray X0_global, ref int[] ref_I, 
            //ref int[] _CellIndexMap, ref int[] _NodeIndexMap, 
            out int[] Permutation, out int[] out_jSub2jCell) {
            int[] I = ref_I;
            int JSUB = I.Length - 1;
            if(JSUB != sgrd.LocalNoOfCells)
                throw new ArgumentException();
            int NoOfPts = X0_global.GetLength(0);
            GridData _GridData = (GridData)(sgrd.GridData);
            int D = _GridData.SpatialDimension;
            //int[] CellIndexMap = _CellIndexMap;
            //int[] NodeIndexMap = _NodeIndexMap;
            
            //if(CellIndexMap.Length != NoOfPts)
            //    throw new ArgumentException();
            //if(NodeIndexMap.Length != NoOfPts)
            //    throw new ArgumentException();
            if(X0_global.Dimension != 2)
                throw new ArgumentException();
            if(X0_global.GetLength(0) != NoOfPts)
                throw new ArgumentException();
            if(X0_global.GetLength(1) != D)
                throw new ArgumentException();
            if(I[0] != 0)
                throw new ArgumentException();
            if(I[JSUB] != NoOfPts)
                throw new ArgumentException();

            int[] jSub_2_jCell = sgrd.SubgridIndex2LocalCellIndex;
            
            double[] _X0 = new double[D];
            double[] _X0_G = new double[D];
            double[] _X0_L = new double[D];
            double[] _X0_L2 = new double[D];
            double[] dummyPt1 = new double[D];
            double[] dummyPt2 = new double[D];

            var h_min = _GridData.iGeomCells.h_min;

#if DEBUG
            BitArray sgrdMask = sgrd.VolumeMask.GetBitMask().CloneAs();
#endif

            BitArray newSgrdMask = new BitArray(_GridData.Cells.NoOfLocalUpdatedCells);
            int[] Point2CellIndex = new int[NoOfPts];
            int newNoOfPts = 0;


            // ==================================
            // find the cells for all points
            // ==================================

            // loop over all cells in the subgrid...
            bool doResort = false; // if this flag stays false, all points stay whithin their cells and no resorting is necessary
            for(int jsub = 0; jsub < JSUB; jsub++) {
                int I0 = I[jsub]; // ............................................. index of first point in cell
                int IE = I[jsub + 1]; // ......................................... index of last point in next cell + 1
                int K = IE - I0; // .............................................. number of points in cell
                int jCell = jSub_2_jCell[jsub]; // ............................... cell index
                RefElement Kref = _GridData.Cells.GetRefElement(jCell); // ....... reference element

#if DEBUG
                Debug.Assert(sgrdMask[jCell] == true);
#endif

                if(K <= 0)
                    // no points in cell
                    continue;

                // transform points to local coordinates in 'jCell'
                // ------------------------------------------------
                MultidimensionalArray X0_loc = MultidimensionalArray.Create(K, D);
                MultidimensionalArray X0_glb = X0_global.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });
                _GridData.TransformGlobal2Local(X0_glb, X0_loc, jCell, null); // transform to local


                // test if any nodes have left their cell, and find new cells
                // ----------------------------------------------------------

                int[] Neighs = null, dummy = null;
                //double h_min_Neighs = 0.0;   // minimum cell diameter

                for(int nn = 0; nn < K; nn++) { // loop over all nodes in sub-cell 'jsub'...
                   
                    for(int d = 0; d < D; d++) {
                        _X0[d] = X0_loc[nn, d];
                        _X0_G[d] = X0_glb[nn, d];
                    }

#if DEBUG
                    Debug.Assert(Kref.IsWithin(_X0) == _GridData.Cells.IsInCell(_X0_G, jCell));
#endif

                    if(!Kref.IsWithin(_X0)) {

                        // the point may has jumped to another cell...
                        // +++++++++++++++++++++++++++++++++++++++++++


                        doResort = true;  // it is necessary to sort the X0-vertices again.

                        if(Neighs == null) {
                            _GridData.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaVertices, out Neighs, out dummy);
                        }


                        // first search attempt ...
                        // ----------------------------

                        int new_jCell = int.MinValue;
                        double dist = double.MaxValue;
                        for(int nc = 0; nc < Neighs.Length; nc++) { // search all cell neighbours, whether they contain the point...
                            int trial_jCell = Neighs[nc];

                            if(_GridData.Cells.IsInCell(_X0_G, trial_jCell, _X0_L)) {
                                new_jCell = trial_jCell;
                                dist = 0.0;
                                break;
                            }
                        }

                        // more sophisticated second search attempt, if the first one failed...
                        // --------------------------------------------------------------------
                        
                        if(dist != 0.0) {
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // search in neigbour cells failed, maybe due to 
                            // round-of errors in the 'IsInCell'-routine
                            // Therefore, we have to do a more sophisticated search:
                            // we are now seaching for the minimum distance
                            //                             ----------------
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                            for(int nc = -1; nc < Neighs.Length; nc++) {
                                int trial_jCell;
                                if(nc < 0)
                                    trial_jCell = jCell; // we also search in the current cell, because the point might just be at the edge
                                //                          and 'IsInCell(...)' just didn't recognized it.
                                else
                                    trial_jCell = Neighs[nc];


                                double dist_to_trialCell = _GridData.Cells.ClosestPointInCell(_X0_G, trial_jCell, _X0_L2, dummyPt1, dummyPt2);

                                if(dist_to_trialCell < dist) {
                                    new_jCell = trial_jCell;
                                    dist = dist_to_trialCell;
                                    Array.Copy(_X0_L2, _X0_L, D);
                                }
                            }
                        }

                        double cell_h = h_min[new_jCell];
                        if(dist <= 0.01 * cell_h) {
                            newSgrdMask[new_jCell] = true;
                            Point2CellIndex[I0 + nn] = new_jCell;
                            newNoOfPts++;
                            
                        } else {
                            // ++++++++++++++++++++++++++
                            // unable to locate the point
                            // ++++++++++++++++++++++++++

                            // just forget about it....

                            Point2CellIndex[I0 + nn] = int.MaxValue;
                        }

                    } else {

                        newSgrdMask[jCell] = true;
                        Point2CellIndex[I0 + nn] = jCell;
                        newNoOfPts++;
                    }
                }

            }

            


            // ==================
            // re-sort the points
            // ==================


            SubGrid newSgrd;


            // resort nodes, if necessary
            if(doResort) {

                // create the subgrid with new cut-cells
                // -------------------------------------
#if DEBUG
                newSgrd = new SubGrid(new CellMask(_GridData, newSgrdMask));
                int[] new_jSub2jCell = newSgrd.SubgridIndex2LocalCellIndex;
#endif
                // determine order: sort according to cell-index
                // ---------------------------------------------

                int[] NewIndex = new int[NoOfPts];
                for(int i =  NewIndex.Length - 1; i >= 0; i--)
                    NewIndex[i] = i;
                Permutation = NewIndex;

                Array.Sort(NewIndex, delegate(int a, int b) {
                    int diff;
                    diff = Point2CellIndex[a] - Point2CellIndex[b];
                    return diff;
                });

#if DEBUG
                for(int i = newNoOfPts; i < NoOfPts; i++) {
                    Debug.Assert(NewIndex[Point2CellIndex[i]] == int.MaxValue); // int.MaxValue marks the points which o to nirvana
                }

                for(int i = 1; i < NoOfPts; i++) {
                    Debug.Assert(Point2CellIndex[NewIndex[i - 1]] <= Point2CellIndex[NewIndex[i]]); // just check that i do not mix up the pointers...
                }
#endif

                
                // sort data
                // ---------
                
                //int[] new_R_NodeIndexMap = new int[newNoOfPts];
                //int[] new_R_CellIndexMap = new int[newNoOfPts];
                //var new_R_X0_local = MultidimensionalArray.Create(newNoOfPts, D);
                var new_R_X0_global = MultidimensionalArray.Create(newNoOfPts, D);

                for(int i = 0; i < newNoOfPts; i++) {
                    for(int d = 0; d < D; d++) {
                        new_R_X0_global[i, d] = X0_global[Permutation[i], d];
                    }
                }
                X0_global = new_R_X0_global;

                // -----------

                List<int> newI = new List<int>();
                List<int> CellIdx = new List<int>();
                newI.Add(0);
                int jSub = 0;
                int NodesPerSubCell;
                int jCell;
                if(newNoOfPts > 0) {
                    jCell = Point2CellIndex[NewIndex[0]];
                    CellIdx.Add(jCell);
                    NodesPerSubCell = 1;
                } else {
                    jCell = -1;
                    NodesPerSubCell = -1;
                }

                for(int hh = 1; hh < newNoOfPts; hh++) {
                    int tt = NewIndex[hh];

                    int jCellPt = Point2CellIndex[tt];
                    if(jCell != jCellPt) {
                        Debug.Assert(NodesPerSubCell > 0);
                        jSub++;
                        newI.Add(newI[jSub - 1] + NodesPerSubCell);
                        CellIdx.Add(jCellPt);
                        jCell = jCellPt;
                        NodesPerSubCell = 1;
                    } else {
                        NodesPerSubCell++;
#if DEBUG
                        Debug.Assert(new_jSub2jCell[jSub] >= jCellPt);
#endif
                    }
#if DEBUG
                    Debug.Assert(newSgrdMask[jCell] == true);
#endif

                }

                Debug.Assert(NodesPerSubCell > 0);
                jSub++;
                newI.Add(newI[jSub - 1] + NodesPerSubCell);

                ref_I = newI.ToArray();
                out_jSub2jCell = CellIdx.ToArray();
#if DEBUG
                Debug.Assert(ref_I.Length == newSgrd.LocalNoOfCells + 1);
#endif
            } else {
                Permutation = NoOfPts.ForLoop(i => i);
                out_jSub2jCell = jSub_2_jCell;
            }





        }


    }
}
