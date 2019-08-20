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
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Numerics;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.Tracing;
using ilPSP.LinSolvers;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.XDG;

using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;

namespace BoSSS.Application.XNSE_Solver  {

    public abstract class XBase_Solver<T> : BoSSS.Solution.Application<T>
        where T : XNSE_Control, new() {

        // comment in for start project
        // ============================
        //static void Main(string[] args)
        //{
        //    _Main(args, false, delegate () {
        //        var p = new XNSE_SolverMain();
        //        return p;
        //    });
        //}


        /// <summary>
        /// the DG representation of the level set.
        /// This one is used for level-set evolution in time; it is in general discontinuous.
        /// </summary>
        //[InstantiateFromControlFile("PhiDG", "Phi", IOListOption.ControlFileDetermined)]
        protected ScalarFieldHistory<SinglePhaseField> DGLevSet;

        /// <summary>
        /// The continuous level set field which defines the XDG space; 
        /// it is obtained from the projection of the discontinuous <see cref="DGLevSet"/> onto the 
        /// continuous element space.
        /// </summary>
        //[InstantiateFromControlFile("Phi", "Phi", IOListOption.ControlFileDetermined)]
        protected LevelSet LevSet;

        /// <summary>
        /// Guess what?
        /// </summary>
        protected VectorField<SinglePhaseField> DGLevSetGradient;

        protected VectorField<SinglePhaseField> LevSetGradient;

        /// <summary>
        /// Curvature; DG-polynomial degree should be 2 times the polynomial degree of <see cref="LevSet"/>.
        /// </summary>
        [InstantiateFromControlFile(VariableNames.Curvature, VariableNames.Curvature, IOListOption.ControlFileDetermined)]
        protected SinglePhaseField Curvature;



        /// <summary>
        /// If requested, performs the projection of the level-set on a continuous field
        /// </summary>
        protected ContinuityProjection ContinuityEnforcer;

        /// <summary>
        /// Lauritz' Fast Marching Solver
        /// !!! Caution !!! Only works in Single-Core
        /// </summary>
        protected FastMarchReinit FastMarchReinitSolver;


        /// <summary>
        /// Bundling of variables which are either DG or XDG (see <see cref="XNSE_Control.UseXDG4Velocity"/>);
        /// </summary>
        protected class VelocityRelatedVars<TX> where TX : DGField
        {
            /// <summary>
            /// velocity
            /// </summary>
            [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                null,
                true, true,
                IOListOption.ControlFileDetermined)]
            public VectorField<TX> Velocity;


            /// <summary>
            /// Volume Force, dimension is acceleration, i.e. length per time-square.
            /// </summary>
            [InstantiateFromControlFile(
                new string[] { VariableNames.GravityX, VariableNames.GravityY, VariableNames.GravityZ },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true,
                IOListOption.ControlFileDetermined)]
            public VectorField<TX> Gravity;

            /// <summary>
            /// Residual in the momentum equation.
            /// </summary>
            [InstantiateFromControlFile(new string[] { "ResidualMomentumX", "ResidualMomentumY", "ResidualMomentumZ" },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true,
                IOListOption.ControlFileDetermined)]
            public VectorField<TX> ResidualMomentum;
        }


        /// <summary>
        /// Velocity and related variables for the extended case, <see cref="XNSE_Control.UseXDG4Velocity"/> == false.
        /// </summary>
        VelocityRelatedVars<XDGField> XDGvelocity;


        /// <summary>
        /// The velocity for the level-set evolution; 
        /// since the velocity representation (<see cref="XDGvelocity"/>) is in the XDG space, int cannot be used directly for the level-set evolution.
        /// </summary>
        [InstantiateFromControlFile(
            new string[] { "Extension" + VariableNames.VelocityX, "Extension" + VariableNames.VelocityY, "Extension" + VariableNames.VelocityZ },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        VectorFieldHistory<SinglePhaseField> ExtensionVelocity;


        /// <summary>
        /// Motion Algorithm for a Extension Velocity based on the density averaged velocity directly at the interface;
        /// </summary>
        ExtensionVelocityBDFMover ExtVelMover;




        IncompressibleMultiphaseBoundaryCondMap m_BcMap;

        /// <summary>
        /// Boundary conditions.
        /// </summary>
        IncompressibleMultiphaseBoundaryCondMap BcMap {
            get {
                if (m_BcMap == null)
                {
                    m_BcMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
                }
                return m_BcMap;
            }
        }


        /// <summary>
        /// HMF order/degree which is used globally in this solver.
        /// </summary>
        int m_HMForder;


        int bdfOrder = -1000;


        // =========
        // level-set
        // =========
        #region level-set

        /// <summary>
        /// Information of the current Fourier Level-Set
        /// DFT_coeff
        /// </summary>
        FourierLevSetBase Fourier_LevSet;

        FourierLevSetTimestepper Fourier_Timestepper;

        /// <summary>
        /// saves the vector Guid for the sample points 
        /// </summary>
        TextWriter Log_FourierLS;

        /// <summary>
        /// init routine for the specialized Fourier level-set
        /// </summary>
        private void InitFourier()
        {
            if (this.Control.FourierLevSetControl == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of FourierLevSetControl!");

            Fourier_LevSet = FourierLevelSetFactory.Build(this.Control.FourierLevSetControl);
            if (this.Control.EnforceLevelSetConservation)
            {
                throw new NotSupportedException("mass conservation correction currently not supported");
            }
            Fourier_LevSet.ProjectToDGLevelSet(this.DGLevSet.Current, this.LsTrk);

            if (base.MPIRank == 0 && this.CurrentSessionInfo.ID != Guid.Empty)
            {
                // restart information for Fourier LS
                Log_FourierLS = base.DatabaseDriver.FsDriver.GetNewLog("Log_FourierLS", this.CurrentSessionInfo.ID);
                Guid vecSamplP_id = this.DatabaseDriver.SaveVector<double>(Fourier_LevSet.getRestartInfo());
                Log_FourierLS.WriteLine(vecSamplP_id);
                Log_FourierLS.Flush();

                //if(this.Control.FourierLevSetControl.WriteFLSdata)
                //    Fourier_LevSet.setUpLogFiles(base.DatabaseDriver, this.CurrentSessionInfo, TimestepNo, PhysTime);

            }
            //create specialized fourier timestepper
            Fourier_Timestepper = FourierLevelSetFactory.Build_Timestepper(this.Control.FourierLevSetControl, Fourier_LevSet.GetFLSproperty(),
                                                            Fourier_LevSet.ComputeChangerate, Fourier_LevSet.EvolveFourierLS);
        }


        /// <summary>
        /// setUp for the Level set initialization (Level-set algorithm, continuity, conservation)
        /// </summary>
        private void InitLevelSet()
        {
            using (new FuncTrace())
            {

                // check level-set
                if (this.LevSet.L2Norm() == 0)
                {
                    throw new NotSupportedException("Level set is not initialized - norm is 0.0 - ALL cells will be cut, no gradient can be defined!");
                }

                // tracker needs to be updated to get access to the cut-cell mask
                this.LsTrk.UpdateTracker();

                // ==============================
                // level-set initialization
                // ==============================

                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 0 }), 3);

                #region Initialize Level Set Evolution Algorithm
                switch (this.Control.Option_LevelSetEvolution)
                {
                    case LevelSetEvolution.Fourier:
                        InitFourier();
                        break;
                    case LevelSetEvolution.None:
                        if (this.Control.AdvancedDiscretizationOptions.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Fourier)
                        {
                            Fourier_LevSet = FourierLevelSetFactory.Build(this.Control.FourierLevSetControl);
                            Fourier_LevSet.ProjectToDGLevelSet(this.DGLevSet.Current, this.LsTrk);
                        }
                        else
                        {
                            goto default;
                        }
                        break;
                    case LevelSetEvolution.ExtensionVelocity:
                        {
                            // Create ExtensionVelocity Motion Algorithm
                            this.DGLevSet.Current.Clear();
                            this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);
                            DGLevSetGradient.Gradient(1.0, DGLevSet.Current);
                            //VectorField<SinglePhaseField> VectorExtVel = ExtensionVelocity.Current;
                            base.RegisterField(ExtensionVelocity.Current);

                            //ReInitPDE = new EllipticReInit(this.LsTrk, this.Control.ReInitControl, DGLevSet.Current);
                            FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);

                            // full initial reinitialization
                            //ReInitPDE.ReInitialize(Restriction: LsTrk.Regions.GetNearFieldSubgrid(1));

                            CellMask Accepted = LsTrk.Regions.GetNearFieldMask(1);
                            CellMask ActiveField = Accepted.Complement();
                            CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                            FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            //ReInitPDE.ReInitialize();

                            // setup extension velocity mover
                            switch (this.Control.Timestepper_Scheme)
                            {
                                case XNSE_Control.TimesteppingScheme.RK_CrankNicolson:
                                case XNSE_Control.TimesteppingScheme.CrankNicolson:
                                    {
                                        //do not instantiate rksch, use bdf instead
                                        bdfOrder = -1;
                                        break;
                                    }
                                case XNSE_Control.TimesteppingScheme.RK_ImplicitEuler:
                                case XNSE_Control.TimesteppingScheme.ImplicitEuler:
                                    {
                                        //do not instantiate rksch, use bdf instead
                                        bdfOrder = 1;
                                        break;
                                    }
                                default:
                                    {
                                        if (this.Control.Timestepper_Scheme.ToString().StartsWith("BDF"))
                                        {
                                            //do not instantiate rksch, use bdf instead
                                            bdfOrder = Convert.ToInt32(this.Control.Timestepper_Scheme.ToString().Substring(3));
                                            break;
                                        }
                                        else
                                            throw new NotImplementedException();
                                    }
                            }

                            ExtVelMover = new ExtensionVelocityBDFMover(LsTrk, DGLevSet.Current, DGLevSetGradient, new VectorField<DGField>(XDGvelocity.Velocity.ToArray()),
                                Control.EllipticExtVelAlgoControl, BcMap, bdfOrder, ExtensionVelocity.Current, new double[2] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B });


                            break;
                        }
                    case LevelSetEvolution.FastMarching:
                    case LevelSetEvolution.Prescribed:
                    case LevelSetEvolution.ScalarConvection:
                    default:
                        // evolution algorithms need a signed-distance level-set:
                        // do some reinit at startup
                        //BoSSS.Solution.LevelSetTools.Advection.NarrowMarchingBand.CutCellReinit(this.LsTrk, this.DGLevSet.Current);
                        // apply only the minimal necessary change
                        this.DGLevSet.Current.Clear();
                        this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);

                        //FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);

                        break;
                }
                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 1 }), 3);
                #endregion

                // =========================================
                // Enforcing the continuity of the level-set
                // =========================================

                ContinuityEnforcer = new ContinuityProjection(
                    ContBasis: this.LevSet.Basis,
                    DGBasis: this.DGLevSet.Current.Basis,
                    gridData: GridData,
                    Option: Control.LSContiProjectionMethod
                    );

                //var CC = this.LsTrk.Regions.GetCutCellMask4LevSet(0);
                var Near1 = this.LsTrk.Regions.GetNearMask4LevSet(0, 1);
                var Near = this.LsTrk.Regions.GetNearMask4LevSet(0, this.Control.LS_TrackerWidth);
                var PosFF = this.LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;

                if (this.Control.Option_LevelSetEvolution != LevelSetEvolution.ExtensionVelocity)
                    ContinuityEnforcer.SetFarField(this.DGLevSet.Current, Near1, PosFF);

                ContinuityEnforcer.MakeContinuous(this.DGLevSet.Current, this.LevSet, Near, PosFF);

                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 2 }), 3);

                this.LsTrk.UpdateTracker();

            }

        }

        /// <summary>
        /// 
        /// </summary>
        public void PushLevelSetAndRelatedStuff()
        {

            if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.Fourier)
            {
                Fourier_Timestepper.updateFourierLevSet();
            }

            this.ExtensionVelocity.IncreaseHistoryLength(1);
            this.ExtensionVelocity.Push();

            this.DGLevSet.IncreaseHistoryLength(1);
            this.DGLevSet.Push();
        }


        int hack_TimestepIndex;


        /// <summary>
        /// Computes the new level set field at time <paramref name="Phystime"/> + <paramref name="dt"/>.
        /// This is a 'driver function' which provides a universal interface to the various level set evolution algorithms.
        /// It also acts as a callback to the time stepper (see <see cref="m_BDF_Timestepper"/> resp. <see cref="m_RK_Timestepper"/>),
        /// i.e. it matches the signature of 
        /// <see cref="BoSSS.Solution.XdgTimestepping.DelUpdateLevelset"/>.
        /// </summary>
        /// <param name="Phystime"></param>
        /// <param name="dt"></param>
        /// <param name="CurrentState">
        /// The current solution (velocity and pressure), since the complete solution is provided by the time stepper,
        /// only the velocity components(supposed to be at the beginning) are used.
        /// </param>
        /// <param name="underrelax">
        /// </param>
        double DelUpdateLevelSet(DGField[] CurrentState, double Phystime, double dt, double underrelax, bool incremental)
        {
            using (new FuncTrace())
            {

                //dt *= underrelax;
                int D = base.Grid.SpatialDimension;
                int iTimestep = hack_TimestepIndex;
                DGField[] EvoVelocity = CurrentState.GetSubVector(0, D);


                // ========================================================
                // Backup old level-set, in order to compute the residual
                // ========================================================

                SinglePhaseField LsBkUp = new SinglePhaseField(this.LevSet.Basis);
                LsBkUp.Acc(1.0, this.LevSet);
                CellMask oldCC = LsTrk.Regions.GetCutCellMask();

                // ====================================================
                // set evolution velocity, but only on the CUT-cells
                // ====================================================

                #region Calculate density averaged Velocity for each cell

                ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(EvoVelocity);

                #endregion


                // ===================================================================
                // backup interface properties (mass conservation, surface changerate)
                // ===================================================================

                #region backup interface props

                double oldSurfVolume = 0.0;
                double oldSurfLength = 0.0;
                double SurfChangerate = 0.0;
                if (this.Control.CheckInterfaceProps)
                {
                    oldSurfVolume = XNSEUtils.GetSpeciesArea(this.LsTrk, LsTrk.GetSpeciesId("A"));
                    oldSurfLength = XNSEUtils.GetInterfaceLength(this.LsTrk);
                    SurfChangerate = EnergyUtils.GetSurfaceChangerate(this.LsTrk, meanVelocity, this.m_HMForder);
                }

                #endregion

                // ====================================================
                // perform level-set evolution
                // ====================================================

                #region level-set evolution

                // set up for Strang splitting
                SinglePhaseField DGLevSet_old;
                if (incremental)
                    DGLevSet_old = this.DGLevSet.Current.CloneAs();
                else
                    DGLevSet_old = this.DGLevSet[0].CloneAs();


                // set up for underrelaxation
                SinglePhaseField DGLevSet_oldIter = this.DGLevSet.Current.CloneAs();

                //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 0 }), 2);

                // actual evolution
                switch (this.Control.Option_LevelSetEvolution)
                {
                    case LevelSetEvolution.None:
                        throw new ArgumentException("illegal call");

                    case LevelSetEvolution.FastMarching:
                        {

                            NarrowMarchingBand.Evolve_Mk2(
                             dt, this.LsTrk, DGLevSet_old, this.DGLevSet.Current, this.DGLevSetGradient,
                             meanVelocity, this.ExtensionVelocity.Current.ToArray(), //new DGField[] { LevSetSrc },
                             this.m_HMForder, iTimestep);

                            //FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                            //CellMask Accepted = LsTrk.Regions.GetCutCellMask();
                            //CellMask ActiveField = LsTrk.Regions.GetNearFieldMask(1);
                            //CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                            //FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            break;
                        }

                    case LevelSetEvolution.Fourier:
                        {
                            Fourier_Timestepper.moveLevelSet(dt, meanVelocity);
                            if (incremental)
                                Fourier_Timestepper.updateFourierLevSet();
                            Fourier_LevSet.ProjectToDGLevelSet(this.DGLevSet.Current, this.LsTrk);
                            break;
                        }

                    case LevelSetEvolution.Prescribed:
                        {
                            this.DGLevSet.Current.Clear();
                            this.DGLevSet.Current.ProjectField(1.0, Control.Phi.Vectorize(Phystime + dt));
                            break;
                        }

                    case LevelSetEvolution.ScalarConvection:
                        {
                            var LSM = new LevelSetMover(EvoVelocity,
                                this.ExtensionVelocity,
                                this.LsTrk,
                                XVelocityProjection.CutCellVelocityProjectiontype.L2_plain,
                                this.DGLevSet,
                                this.BcMap);

                            int check1 = this.ExtensionVelocity.PushCount;
                            int check2 = this.DGLevSet.PushCount;

                            this.DGLevSet[1].Clear();
                            this.DGLevSet[1].Acc(1.0, DGLevSet_old);
                            LSM.Advect(dt);

                            if (check1 != this.ExtensionVelocity.PushCount)
                                throw new ApplicationException();
                            if (check2 != this.DGLevSet.PushCount)
                                throw new ApplicationException();

                            break;
                        }
                    case LevelSetEvolution.ExtensionVelocity:
                        {

                            DGLevSetGradient.Clear();
                            DGLevSetGradient.Gradient(1.0, DGLevSet.Current);

                            ExtVelMover.Advect(dt);

                            // Fast Marching: Specify the Domains first
                            // Perform Fast Marching only on the Far Field
                            if (this.Control.AdaptiveMeshRefinement)
                            {
                                int NoCells = ((GridData)this.GridData).Cells.Count;
                                BitArray Refined = new BitArray(NoCells);
                                for (int j = 0; j < NoCells; j++)
                                {
                                    if (((GridData)this.GridData).Cells.GetCell(j).RefinementLevel > 0)
                                        Refined[j] = true;
                                }
                                CellMask Accepted = new CellMask(this.GridData, Refined);
                                CellMask AcceptedNeigh = Accepted.AllNeighbourCells();

                                Accepted = Accepted.Union(AcceptedNeigh);
                                CellMask ActiveField = Accepted.Complement();
                                CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                                FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            }
                            else
                            {
                                CellMask Accepted = LsTrk.Regions.GetNearFieldMask(1);
                                CellMask ActiveField = Accepted.Complement();
                                CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                                FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            }
                            //SubGrid AcceptedGrid = new SubGrid(Accepted);
                            //ReInitPDE.ReInitialize(Restriction: AcceptedGrid);

                            //CellMask ActiveField = Accepted.Complement();
                            //CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                            //FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            //ReInitPDE.ReInitialize();

                            break;
                        }
                    default:
                        throw new ApplicationException();
                }


                // performing underrelaxation
                if (underrelax < 1.0)
                {
                    this.DGLevSet.Current.Scale(underrelax);
                    this.DGLevSet.Current.Acc((1.0 - underrelax), DGLevSet_oldIter);
                }

                //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 1 }), 2);


                #endregion


                // ======================
                // postprocessing  
                // =======================

                if (this.Control.ReInitPeriod > 0 && hack_TimestepIndex % this.Control.ReInitPeriod == 0)
                {
                    Console.WriteLine("Filtering DG-LevSet");
                    SinglePhaseField FiltLevSet = new SinglePhaseField(DGLevSet.Current.Basis);
                    FiltLevSet.AccLaidBack(1.0, DGLevSet.Current);
                    Filter(FiltLevSet, 2, oldCC);
                    DGLevSet.Current.Clear();
                    DGLevSet.Current.Acc(1.0, FiltLevSet);

                    Console.WriteLine("FastMarchReInit performing FirstOrderReInit");
                    FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                    CellMask Accepted = LsTrk.Regions.GetCutCellMask();
                    CellMask ActiveField = LsTrk.Regions.GetNearFieldMask(1);
                    CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                    FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);
                }

                #region ensure continuity

                // make level set continuous
                CellMask CC = LsTrk.Regions.GetCutCellMask4LevSet(0);
                CellMask Near1 = LsTrk.Regions.GetNearMask4LevSet(0, 1);
                CellMask PosFF = LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
                ContinuityEnforcer.MakeContinuous(this.DGLevSet.Current, this.LevSet, Near1, PosFF);

                if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.FastMarching)
                {
                    CellMask Nearband = Near1.Union(CC);
                    this.DGLevSet.Current.Clear(Nearband);
                    this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet, Nearband);
                    //ContinuityEnforcer.SetFarField(this.DGLevSet.Current, Near1, PosFF);
                }

                //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 2 }), 2);

                #endregion


                for (int d = 0; d < D; d++)
                    this.XDGvelocity.Velocity[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;



                // ===============
                // tracker update
                // ===============

                this.LsTrk.UpdateTracker(incremental: true);

                // update near field (in case of adaptive mesh refinement)
                if (this.Control.AdaptiveMeshRefinement && this.Control.Option_LevelSetEvolution == LevelSetEvolution.FastMarching)
                {
                    Near1 = LsTrk.Regions.GetNearMask4LevSet(0, 1);
                    PosFF = LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
                    ContinuityEnforcer.SetFarField(this.DGLevSet.Current, Near1, PosFF);
                    ContinuityEnforcer.SetFarField(this.LevSet, Near1, PosFF);
                }


                // ==================================================================
                // check interface properties (mass conservation, surface changerate)
                // ==================================================================

                if (this.Control.CheckInterfaceProps)
                {

                    double currentSurfVolume = XNSEUtils.GetSpeciesArea(this.LsTrk, LsTrk.GetSpeciesId("A"));
                    double massChange = ((currentSurfVolume - oldSurfVolume) / oldSurfVolume) * 100;
                    Console.WriteLine("Change of mass = {0}%", massChange);

                    double currentSurfLength = XNSEUtils.GetInterfaceLength(this.LsTrk);
                    double actualSurfChangerate = (currentSurfLength - oldSurfLength) / dt;
                    Console.WriteLine("Interface divergence = {0}", SurfChangerate);
                    Console.WriteLine("actual surface changerate = {0}", actualSurfChangerate);

                }


                // ==================
                // compute residual
                // ==================

                var newCC = LsTrk.Regions.GetCutCellMask();
                LsBkUp.Acc(-1.0, this.LevSet);
                double LevSetResidual = LsBkUp.L2Norm(newCC.Union(oldCC));

                return LevSetResidual;
            }
        }


        private void Filter(SinglePhaseField FiltrdField, int NoOfSweeps, CellMask CC)
        {

            Basis patchRecoveryBasis = FiltrdField.Basis;

            L2PatchRecovery l2pr = new L2PatchRecovery(patchRecoveryBasis, patchRecoveryBasis, CC, true);

            SinglePhaseField F_org = FiltrdField.CloneAs();

            for (int pass = 0; pass < NoOfSweeps; pass++)
            {
                F_org.Clear();
                F_org.Acc(1.0, FiltrdField);
                FiltrdField.Clear();
                l2pr.Perform(FiltrdField, F_org);
            }
        }


        /// <summary>
        ///  Take density-weighted mean value in cut-cells
        /// </summary>
        /// <param name="EvoVelocity"></param>
        /// <returns></returns>
        private ConventionalDGField[] GetMeanVelocityFromXDGField(DGField[] EvoVelocity)
        {
            int D = EvoVelocity.Length;
            ConventionalDGField[] meanVelocity;

            Debug.Assert(this.XDGvelocity != null);

            meanVelocity = new ConventionalDGField[D];

            double rho_A = this.Control.PhysicalParameters.rho_A, rho_B = this.Control.PhysicalParameters.rho_B;
            double mu_A = this.Control.PhysicalParameters.mu_A, mu_B = this.Control.PhysicalParameters.mu_B;
            CellMask CC = this.LsTrk.Regions.GetCutCellMask4LevSet(0);
            CellMask Neg = this.LsTrk.Regions.GetLevelSetWing(0, -1).VolumeMask;
            CellMask Pos = this.LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
            CellMask posNear = this.LsTrk.Regions.GetNearMask4LevSet(0, 1).Except(Neg);
            CellMask negNear = this.LsTrk.Regions.GetNearMask4LevSet(0, 1).Except(Pos);

            for (int d = 0; d < D; d++)
            {
                Basis b = this.XDGvelocity.Velocity[d].Basis.NonX_Basis;
                meanVelocity[d] = new SinglePhaseField(b);


                foreach (string spc in this.LsTrk.SpeciesNames)
                {
                    double rhoSpc;
                    double muSpc;
                    switch (spc)
                    {
                        case "A": rhoSpc = rho_A; muSpc = mu_A; break;
                        case "B": rhoSpc = rho_B; muSpc = mu_B; break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }

                    double scale = 1.0;
                    switch (this.Control.InterAverage)
                    {
                        case XNSE_Control.InterfaceAveraging.mean:
                            {
                                scale = 0.5;
                                break;
                            }
                        case XNSE_Control.InterfaceAveraging.density:
                            {
                                scale = rhoSpc / (rho_A + rho_B);
                                break;
                            }
                        case XNSE_Control.InterfaceAveraging.viscosity:
                            {
                                scale = muSpc / (mu_A + mu_B);
                                break;
                            }
                    }

                    meanVelocity[d].Acc(scale, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), CC);
                    switch (spc)
                    {
                        //case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), Neg.Except(CC)); break;
                        case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), negNear); break;
                        case "B": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), posNear); break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }
                }

            }

            return meanVelocity;
        }

        #endregion


    }
}
