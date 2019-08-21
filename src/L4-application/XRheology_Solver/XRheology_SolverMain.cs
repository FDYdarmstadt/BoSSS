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
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;

using MPI.Wrappers;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;

using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Control;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.AdvancedSolvers;

using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.RheologyCommon;

using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;

using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;

using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XRheology_Solver {
    public class XRheology_SolverMain : BoSSS.Application.XNSE_Solver.XBase_Solver<XRheology_Control> {
        static void Main(string[] args) {
            XRheology_SolverMain._Main(args, false,
            delegate () {
                var app = new XRheology_SolverMain();
                return app;
            });
        }

        #region instantiation

        // Attributes for fields (Names), initialization of DG fields
        //==============================================================

        /// <summary>
        /// Velocity and related variables for the extended case, <see cref="XNSE_Control.UseXDG4Velocity"/> == false.
        /// </summary>
        VelocityRelatedVars<XDGField> XDGvelocity;

        /// <summary>
        /// Pressure domain
        /// </summary>
        XDGField Pressure;

        /// <summary>
        /// Pressure codomain: Residuum in continuity equation
        /// </summary>
        XDGField ResidualContinuity;

        /// <summary>
        /// Extra stress domain (2D): StressXX
        /// </summary>
        XDGField StressXX;

        /// <summary>
        /// Extra stress domain (2D): StressXY
        /// </summary>
        XDGField StressXY;

        /// <summary>
        /// Extra stress domain (2D): StressYY
        /// </summary>
        XDGField StressYY;

        /// <summary>
        /// Extra stress codomain (2D): StressXX
        /// </summary>
        XDGField ResidualStressXX;

        /// <summary>
        /// Extra stresses codomain (2D): StressXY
        /// </summary>
        XDGField ResidualStressXY;

        /// <summary>
        /// Extra stresses codomain (2D): StressYY
        /// </summary>
        XDGField ResidualStressYY;

        /// <summary>
        /// Extra stresses parameter (2D): StressXX
        /// </summary>
        XDGField StressXXP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        XDGField StressXYP;

        /// <summary>
        /// Extra stresses parameter (2D): StressXY
        /// </summary>
        XDGField StressYYP;

        /// <summary>
        /// Extra source (e.g. gravity)
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.GravityX, VariableNames.GravityY, VariableNames.GravityZ },
                    new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                    true, true,
                    IOListOption.ControlFileDetermined)]
        public VectorField<XDGField> Gravity;

        //// Gravity source constitutive
        //[InstantiateFromControlFile("GravityXX", "StressXX", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityXX;

        //[InstantiateFromControlFile("GravityXY", "StressXY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityXY;

        //[InstantiateFromControlFile("GravityYY", "StressYY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityYY;

        ////Gravity source for divergence of u
        //[InstantiateFromControlFile("GravityDiv", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        //public SinglePhaseField GravityDiv;

        // Parameters: Velocity Gradient
        public VectorField<XDGField> VelocityXGradient;
        public VectorField<XDGField> VelocityYGradient;


        //Parameters: external analytical velocity
        XDGField U;
        XDGField V;


        // necessary XDG Stuff
        // ===========================
        /// <summary>
        /// PDE based elliptic reInitialization by Thomas
        /// </summary>
        EllipticReInit ReInitPDE;

        /// <summary>
        /// Lauritz' Fast Marching Solver
        /// !!! Caution !!! Only works in Single-Core
        /// </summary>
        FastMarchReinit FastMarchReinitSolver;

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

        /// <summary>
        /// Artificial force term at the fluid interface, usually only to support manufactured solutions.
        /// </summary>
        [InstantiateFromControlFile(
            new string[] { VariableNames.SurfaceForceX, VariableNames.SurfaceForceY, VariableNames.SurfaceForceZ },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        VectorField<SinglePhaseField> SurfaceForce;

        // Some initialisation of variables
        //============================================

        /// <summary>
        /// the spatial operator (momentum, continuity and constitutive equation)
        /// </summary>
        XRheology_OperatorFactory XRheology_Operator;

        /// <summary>
        /// OperatorConfiguration for the <see cref="XRheology_Operator"/>
        /// </summary>
        XRheology_OperatorConfiguration XOpConfig;

        //not sure if needed....
        // ===================================================
        /// <summary>
        /// Current Velocity
        /// </summary>
        XDGField[] CurrentVel {
            get {
                return this.XDGvelocity.Velocity.ToArray();
            }
        }

        XDGField[] prevVel;

        /// <summary>
        /// output of <see cref="AssembleMatrix"/>;
        /// </summary>
        MassMatrixFactory MassFact;
        // ======================================================

        int D; // Spatial Dimension
        /// <summary>
        /// current Weissenberg number
        /// </summary>
        public double currentWeissenberg;
        bool ChangeMesh = true;

        /// <summary>
        /// initialisation of BDF Timestepper
        /// </summary>
        protected XdgBDFTimestepping m_BDF_Timestepper;


        // Persson sensor and artificial viscosity
        //=============================================
        /// <summary>
        /// initialisation of Persson sensor
        /// </summary>
        protected PerssonSensor perssonsensor;

        /// <summary>
        /// initialisation of artificial viscosity
        /// </summary>
        protected SinglePhaseField artificalViscosity;

        /// <summary>
        /// initialisation of max value of artificial viscosity
        /// </summary>
        protected double artificialMaxViscosity;


        // Settings for calculation
        //===============================================
        /// <summary>
        /// Set true if Navier Stokes is solved, then the mean velocities as parameters for calculation of convective terms are needed
        /// </summary>
        protected bool U0MeanRequired {
            get {
                return (this.Control.PhysicalParameters.IncludeConvection);
            }
        }

        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double reynolds_A = this.Control.PhysicalParameters.reynolds_A,
                    reynolds_B = this.Control.PhysicalParameters.reynolds_B;

                int D = this.GridData.SpatialDimension;

                double[] _reynolds_A = new double[D + 1];
                _reynolds_A.SetAll(reynolds_A); // mass matrix in momentum equation
                _reynolds_A[D] = 0; // no  mass matrix for continuity equation
                double[] _reynolds_B = new double[D + 1];
                _reynolds_B.SetAll(reynolds_B); // mass matrix in momentum equation
                _reynolds_B[D] = 0; // no  mass matrix for continuity equation


                //double[] _rho = new double[D + 4];
                //_rho.SetAll(rho);
                ////No MassMatrix for the pressure
                //_rho[D] = 0;

                //_rho[D + 1] = 1;
                //_rho[D + 2] = 1;
                //_rho[D + 3] = 1;

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _reynolds_A);
                R.Add(this.LsTrk.GetSpeciesId("B"), _reynolds_B);

                return R;
            }
        }

        CoordinateVector m_CurrentSolution = null;
        CoordinateVector m_CurrentResidual = null;

        /// <summary>
        /// Current solution vector
        /// </summary>
        public CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(XDGvelocity.Velocity, this.Pressure, this.StressXX, this.StressXY, this.StressYY));
                }
                return m_CurrentSolution;
            }
        }

        /// <summary>
        /// Current residual vector
        /// </summary>
        public CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(XDGvelocity.ResidualMomentum, this.ResidualContinuity, this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY));
                }
                return m_CurrentResidual;
            }
        }

        /// <summary>
        /// DG Field instantiation.
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();
            int D = this.GridData.SpatialDimension;

            // LEVEL SET FIELDS
            DGLevSet = new ScalarFieldHistory<SinglePhaseField>(
                   new SinglePhaseField(new Basis(this.GridData, this.Control.FieldOptions["Phi"].Degree), "PhiDG"));

            if (this.Control.FieldOptions["PhiDG"].Degree >= 0 && this.Control.FieldOptions["PhiDG"].Degree != this.DGLevSet.Current.Basis.Degree) {
                throw new ApplicationException("Specification of polynomial degree for 'PhiDG' is not supportet, since it is induced by polynomial degree of 'Phi'.");
            }

            this.LsTrk = new LevelSetTracker((GridData)this.GridData, base.Control.CutCellQuadratureType, base.Control.LS_TrackerWidth, new string[] { "A", "B" }, this.LevSet);
            base.RegisterField(this.LevSet);
            this.LevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(this.LevSet.Basis, "dPhi_dx[" + d + "]")));
            base.RegisterField(this.LevSetGradient);

            base.RegisterField(this.DGLevSet.Current);
            this.DGLevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(this.DGLevSet.Current.Basis, "dPhiDG_dx[" + d + "]")));
            base.RegisterField(this.DGLevSetGradient);


            //PRESSURE FIELD
            this.Pressure = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree), VariableNames.Pressure);
            base.RegisterField(this.Pressure);


            // CONTI FIELD
            this.ResidualContinuity = new XDGField(this.Pressure.Basis, "ResidualConti");
            base.RegisterField(this.ResidualContinuity);


            // ALL VELOCITY RELATED FIELDS
            this.XDGvelocity = new VelocityRelatedVars<XDGField>();
            InitFromAttributes.CreateFieldsAuto(this.XDGvelocity, this.GridData, base.Control.FieldOptions, base.Control.CutCellQuadratureType, base.IOFields, base.m_RegisteredFields);


            // ALL STRESS RELATED FIELDS
            this.StressXX = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), VariableNames.StressXX);
            base.RegisterField(this.StressXX);
            this.StressXY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXY].Degree), VariableNames.StressXY);
            base.RegisterField(this.StressXY);
            this.StressYY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressYY].Degree), VariableNames.StressYY);
            base.RegisterField(this.StressYY);

            this.ResidualStressXX = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), "ResidualStressXX");
            base.RegisterField(this.ResidualStressXX);
            this.ResidualStressXY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXY].Degree), "ResidualStressXY");
            base.RegisterField(this.ResidualStressXY);
            this.ResidualStressYY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressYY].Degree), "ResidualStressYY");
            base.RegisterField(this.ResidualStressYY);

            this.StressXXP = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), VariableNames.StressXXP);
            base.RegisterField(this.StressXXP);
            this.StressXYP = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXY].Degree), VariableNames.StressXYP);
            base.RegisterField(this.StressXYP);
            this.StressYYP = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressYY].Degree), VariableNames.StressYYP);
            base.RegisterField(this.StressYYP);


            //PERSSON SENSOR FIELD
            if (Control.UsePerssonSensor == true) {
                perssonsensor = new PerssonSensor(StressXX);
                this.IOFields.Add(perssonsensor.GetField());
                base.RegisterField(this.perssonsensor.GetField());
            }


            //ARTIFICIAL VISCOSITY FIELD
            if (Control.UseArtificialDiffusion == true) {
                artificalViscosity = new SinglePhaseField(new Basis(GridData, 1), "artificalViscosity");
                this.IOFields.Add(artificalViscosity);
                base.RegisterField(this.artificalViscosity);

            }
        }


        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
        }

        #endregion

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            #region Checks
            // CreateEquationsAndSolvers might be called multiple times
            // exit if so, and no LoadBalancing
            if (XRheology_Operator != null && L == null)
                return;


            if (Control.CompMode == AppControl._CompMode.Steady) {
                if (Control.Timestepper_LevelSetHandling != LevelSetHandling.None)
                    throw new ApplicationException(string.Format("Illegal control file: for a steady computation ({0}), the level set handling must be {1}.", AppControl._CompMode.Steady, LevelSetHandling.None));
            }

            int degU = this.CurrentVel[0].Basis.Degree;
            int stressDegree = this.StressXX.Basis.Degree;

            #endregion

            //Quadrature Order
            //----------------

            m_HMForder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);

            // Create Spatial Operator
            // ======================= 

            XOpConfig = new XRheology_OperatorConfiguration(this.Control);

            XRheology_Operator = new XRheology_OperatorFactory(XOpConfig, this.LsTrk, this.m_HMForder, this.BcMap, degU);

            #region Create Timestepper
            // =======================
            if (L == null) {

                switch (this.Control.Timestepper_Scheme) {
                    case XRheology_Control.TimesteppingScheme.RK_ImplicitEuler: {
                            rksch = RungeKuttaScheme.ImplicitEuler;
                            break;
                        }
                    case XRheology_Control.TimesteppingScheme.RK_CrankNicolson: {
                            rksch = RungeKuttaScheme.CrankNicolson;
                            break;
                        }
                    case XRheology_Control.TimesteppingScheme.CrankNicolson: {
                            //do not instantiate rksch, use bdf instead
                            bdfOrder = -1;
                            break;
                        }
                    case XRheology_Control.TimesteppingScheme.ImplicitEuler: {
                            //do not instantiate rksch, use bdf instead
                            bdfOrder = 1;
                            break;
                        }
                    default: {
                            if (this.Control.Timestepper_Scheme.ToString().StartsWith("BDF")) {
                                //do not instantiate rksch, use bdf instead
                                bdfOrder = Convert.ToInt32(this.Control.Timestepper_Scheme.ToString().Substring(3));
                                break;
                            } else
                                throw new NotImplementedException("The chosen timestepper is not implemented!");
                        }

                }

                if (rksch == null) {
                    m_BDF_Timestepper = new XdgBDFTimestepping(
                        this.CurrentSolution.Mapping.Fields,
                        this.CurrentResidual.Mapping.Fields,
                        LsTrk,
                        true,
                        DelComputeOperatorMatrix, null, DelUpdateLevelSet,
                        (this.Control.CompMode == AppControl._CompMode.Transient) ? bdfOrder : 1,
                        this.Control.Timestepper_LevelSetHandling,
                        this.XOpConfig.mmsd,
                        (this.Control.PhysicalParameters.IncludeConvection) ? SpatialOperatorType.Nonlinear : SpatialOperatorType.LinearTimeDependent,
                        MassScale,
                        this.MultigridOperatorConfig, base.MultigridSequence,
                        this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder,
                        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                        true,
                        this.Control.NonLinearSolver,
                        this.Control.LinearSolver
                        );
                    m_BDF_Timestepper.m_ResLogger = base.ResLogger;
                    m_BDF_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                    m_BDF_Timestepper.Timestepper_Init = (this.Control.CompMode == AppControl._CompMode.Transient) ? this.Control.Timestepper_BDFinit : TimeStepperInit.SingleInit;
                    m_BDF_Timestepper.incrementTimesteps = this.Control.incrementTimesteps;
                    m_BDF_Timestepper.PushLevelSet = this.PushLevelSetAndRelatedStuff;
                    m_BDF_Timestepper.IterUnderrelax = this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? this.Control.LSunderrelax : 1.0;

                    m_BDF_Timestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                    //m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;

                    // solver 
                    this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations; //m_BDF_Timestepper.config_NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;

                    if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverConfig.Code.NewtonGMRES) {
                        m_BDF_Timestepper.XdgSolverFactory.Selfmade_precond =
                                            new Schwarz() {
                                                m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                                    NoOfPartsPerProcess = this.CurrentSolution.Count / 10000,
                                                },
                                                Overlap = 1,
                                                CoarseSolver = new SparseSolver() { WhichSolver = SparseSolver._whichSolver.MUMPS }
                                            };
                    } else {
                        //m_BDF_Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = this.Control.LinearSolver };
                    }

                    m_BDF_Timestepper.XdgSolverFactory.Update(this.Control.NonLinearSolver, this.Control.LinearSolver); //Changes made to configs need to be updated afterwards

                    //Console.WriteLine("noofpartsperprocess = {0}", this.CurrentSolution.Count / 10000);
                } else {

                    throw new NotSupportedException();
                }
            } else {

                //PlotCurrentState(hack_Phystime, new TimestepNumber(hack_TimestepIndex, 12), 2);

                Debug.Assert(object.ReferenceEquals(this.MultigridSequence[0].ParentGrid, this.GridData));


                if (this.Control.AdaptiveMeshRefinement && hack_TimestepIndex == 0) {
                    base.SetInitial();
                    this.InitLevelSet();
                }


                m_BDF_Timestepper.DataRestoreAfterBalancing(L,
                    ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure),
                    ArrayTools.Cat<DGField>(this.XDGvelocity.ResidualMomentum.ToArray(), this.ResidualContinuity),
                    this.LsTrk, this.MultigridSequence);

                //PlotCurrentState(hack_Phystime, new TimestepNumber(hack_TimestepIndex, 13), 2);

                ContinuityEnforcer = new ContinuityProjection(
                    ContBasis: this.LevSet.Basis,
                    DGBasis: this.DGLevSet.Current.Basis,
                    gridData: GridData,
                    Option: Control.LSContiProjectionMethod);

                if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
                    ReInitPDE = new EllipticReInit(this.LsTrk, this.Control.ReInitControl, DGLevSet.Current);
                    FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                    ExtVelMover = new ExtensionVelocityBDFMover(LsTrk, DGLevSet.Current, DGLevSetGradient, new VectorField<DGField>(XDGvelocity.Velocity.ToArray()),
                    Control.EllipticExtVelAlgoControl, BcMap, bdfOrder, ExtensionVelocity.Current, new double[2] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B });
                }

            }
        }
        #endregion

        void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {

            int D = this.GridData.SpatialDimension;

            // ============================
            // treatment of surface tension
            // ============================

            VectorField<SinglePhaseField> filtLevSetGradient;
            switch (this.Control.AdvancedDiscretizationOptions.SST_isotropicMode) {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine: {
                        CurvatureAlgorithms.LaplaceBeltramiDriver(
                            this.Control.AdvancedDiscretizationOptions.SST_isotropicMode,
                            this.Control.AdvancedDiscretizationOptions.FilterConfiguration,
                            out filtLevSetGradient, this.LsTrk,
                            this.DGLevSet.Current);
                        if (this.Control.AdvancedDiscretizationOptions.CurvatureNeeded) {
                            VectorField<SinglePhaseField> filtLevSetGradient_dummy;
                            CurvatureAlgorithms.CurvatureDriver(
                                SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                                CurvatureAlgorithms.FilterConfiguration.Default,
                                this.Curvature, out filtLevSetGradient_dummy, this.LsTrk,
                                this.m_HMForder,
                                this.DGLevSet.Current);
                        }
                        break;
                    }
                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                    CurvatureAlgorithms.CurvatureDriver(
                        this.Control.AdvancedDiscretizationOptions.SST_isotropicMode,
                        this.Control.AdvancedDiscretizationOptions.FilterConfiguration,
                        this.Curvature, out filtLevSetGradient, this.LsTrk,
                        this.m_HMForder,
                        this.DGLevSet.Current);
                    //CurvatureAlgorithms.MakeItConservative(LsTrk, this.Curvature, this.Control.PhysicalParameters.Sigma, this.SurfaceForce, filtLevSetGradient, MomentFittingVariant, this.m_HMForder);
                    break;

                case SurfaceStressTensor_IsotropicMode.Curvature_Fourier:
                    if (Fourier_LevSet != null) {
                        Fourier_LevSet.ProjectToDGCurvature(this.Curvature, out filtLevSetGradient, this.LsTrk.Regions.GetCutCellMask());
                    } else {
                        throw new NotImplementedException("Curvature_Fourier needs an instance of Fourier_LevSet");
                    }
                    break;

                default: throw new NotImplementedException("Unknown SurfaceTensionMode");
            }

            // ============================
            // matrix assembly
            // ============================

            var codMap = Mapping;
            var domMap = Mapping;

            this.XRheology_Operator.AssembleMatrix(
                OpMtx, OpAffine, codMap, domMap,
                CurrentState, AgglomeratedCellLengthScales, phystime,
                this.m_HMForder, SurfaceForce, filtLevSetGradient, Curvature);


            if (filtLevSetGradient != null) {
                if (this.Control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource == CurvatureAlgorithms.LevelSetSource.fromC0) {
                    this.LevSetGradient.Clear();
                    this.LevSetGradient.Acc(1.0, filtLevSetGradient);
                } else {
                    this.DGLevSetGradient.Clear();
                    this.DGLevSetGradient.Acc(1.0, filtLevSetGradient);
                }
            }



            // ====================================
            // something with surface tension ?????
            // ====================================

            {
                if (this.Control.PhysicalParameters.useArtificialSurfaceForce == true)
                    throw new NotSupportedException("Not supported for this hack.");


                var TmpRhs = new CoordinateVector(CurrentState.Select(f => (DGField)f.Clone()).ToArray());
                TmpRhs.Clear();

                var VelA = new CoordinateVector(TmpRhs.Mapping.Fields.Take(D).Select(f => (DGField)(((XDGField)f).GetSpeciesShadowField("A"))).ToArray());
                var VelB = new CoordinateVector(TmpRhs.Mapping.Fields.Take(D).Select(f => (DGField)(((XDGField)f).GetSpeciesShadowField("B"))).ToArray());

                int N = ((XDGBasis)(CurrentState[0].Basis)).NonX_Basis.Length;

                foreach (int jCell in this.LsTrk.Regions.GetCutCellMask4LevSet(0).ItemEnum) {
                    for (int d = 0; d < D; d++) {
                        for (int n = 0; n < N; n++) {
                            ((XDGField)(TmpRhs.Mapping.Fields[d])).GetSpeciesShadowField("A").Coordinates[jCell, n] = 0.5 * SurfaceForce[d].Coordinates[jCell, n];
                            ((XDGField)(TmpRhs.Mapping.Fields[d])).GetSpeciesShadowField("B").Coordinates[jCell, n] = 0.5 * SurfaceForce[d].Coordinates[jCell, n];
                        }
                    }
                }

                OpAffine.AccV(1.0, TmpRhs);
            }

            // so far, 'SaddlePointRHS' is on the left-hand-side, since it is the output of ComputeMatrix
            // multiply by -1 to make it RHS
            OpAffine.ScaleV(-1.0);


            // ============================
            // Generate MassMatrix
            // ============================

            // mass matrix factory
            MassFact = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), m_HMForder, 1).MassMatrixFactory;// new MassMatrixFactory(maxB, CurrentAgg);
            var WholeMassMatrix = MassFact.GetMassMatrix(Mapping, MassScale); // mass matrix scaled with density rho


            // ============================
            //  Add Gravity
            // ============================
            // Dimension: [ rho * G ] = mass / time^2 / len^2 == [ d/dt rho U ]
            var WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity.Gravity.ToArray<DGField>(), new XDGField(this.Pressure.Basis)));
            WholeMassMatrix.SpMV(1.0, WholeGravity, 1.0, OpAffine);

            // ============================
            // Gravity Source (default should be zero!) for manufactured solution
            // ============================
            if (Control.GravitySource == true) {
                bool test = false;

                if (this.Control.GravityX != null && this.Control.GravityY != null) {
                    Gravity[0].ProjectField(this.Control.GravityX.Vectorize(0.0));
                    Gravity[1].ProjectField(this.Control.GravityY.Vectorize(0.0));
                    int[] MomEqIdx = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 0, 1);
                    OpAffine.AccV(-1.0, this.Gravity.CoordinateVector, MomEqIdx, default(int[]));
                    test = true;
                }

                //if (this.Control.GravityXX != null && this.Control.GravityXY != null && this.Control.GravityYY != null)
                //{
                //    GravityXX.ProjectField(this.Control.GravityXX.Vectorize(0.0));
                //    int[] ConstEqIdx1 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 3);
                //    OpAffine.AccV(-1.0, this.GravityXX.CoordinateVector, ConstEqIdx1, default(int[]));

                //    GravityXY.ProjectField(this.Control.GravityXY.Vectorize(0.0));
                //    int[] ConstEqIdx2 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 4);
                //    OpAffine.AccV(-1.0, this.GravityXY.CoordinateVector, ConstEqIdx2, default(int[]));

                //    GravityYY.ProjectField(this.Control.GravityYY.Vectorize(0.0));
                //    int[] ConstEqIdx3 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 5);
                //    OpAffine.AccV(-1.0, this.GravityYY.CoordinateVector, ConstEqIdx3, default(int[]));
                //    test = true;
                //}

                //if (this.Control.GravityDiv != null)
                //{
                //    GravityDiv.ProjectField(this.Control.GravityDiv.Vectorize(0.0));
                //    int[] ContiEqIdx = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 2);
                //    OpAffine.AccV(-1.0, this.GravityDiv.CoordinateVector, ContiEqIdx, default(int[]));
                //    test = true;
                //}

                if (!test) {
                    throw new ApplicationException("Gravity is true, but no values set!");
                }
            }

            // ============================
            // Set Pressure Reference Point
            // ============================

            if (OpMtx != null) {
                if (!this.BcMap.DirichletPressureBoundary) {
                    XNSEUtils.SetPressureReferencePoint(
                        Mapping,
                        this.GridData.SpatialDimension,
                        this.LsTrk, OpMtx, OpAffine);
                }
            } else {
                if (!this.BcMap.DirichletPressureBoundary) {
                    XNSEUtils.SetPressureReferencePointResidual(
                        new CoordinateVector(CurrentState),
                        this.GridData.SpatialDimension,
                        this.LsTrk, OpAffine);
                }
            }

            // transform from RHS to Affine
            OpAffine.ScaleV(-1.0);

        }

        double hack_Phystime;

        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (var tr = new FuncTrace()) {

                if (this.Control.OperatorMatrixAnalysis == true) {

                    OpAnalysisBase myAnalysis = new OpAnalysisBase(DelComputeOperatorMatrix, CurrentSolution.Mapping, CurrentSolution.Mapping.Fields.ToArray(), null, phystime);
                    //myAnalysis.VarGroup = new int[] { 0};
                    myAnalysis.Analyse();
                }

                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;
                base.ResLogger.TimeStep = TimestepInt;

                if (TimestepNo[0] > 1) {
                    this.Control.RaiseWeissenberg = false;
                }

                hack_TimestepIndex = TimestepInt;
                hack_Phystime = phystime;
                int NoIncrementTimestep;

                // NOT IMPLEMENTED!!!
                //Preprocessing(TimestepInt, phystime, dt, TimestepNo);


                if (Control.SkipSolveAndEvaluateResidual) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    // setup: project exact solution -- for consistency tests
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    foreach (string spc in LsTrk.SpeciesNames) {
                        for (int d = 0; d < this.GridData.SpatialDimension; d++) {
                            ConventionalDGField Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);
                            Vel_d.ProjectField(Control.ExactSolutionVelocity[spc][d].Convert_Xt2X(phystime + dt));
                        }
                        Pressure.GetSpeciesShadowField(spc).ProjectField(Control.ExactSolutionPressure[spc].Convert_Xt2X(phystime + dt));
                    }
                }


                // =====================================================
                // setup stationary 
                // =====================================================


                if (base.Control.CompMode == AppControl._CompMode.Steady) {
                    dt = 1.0e100;
                    Console.WriteLine("Steady-state solve ...", TimestepNo, dt);

                    if (this.Control.Option_LevelSetEvolution != LevelSetEvolution.None) {
                        throw new ApplicationException("For steady-state solutions, the only allowed level-set-evolution option is '" + LevelSetEvolution.None + "'.");
                    }



                    // =====================================================
                    // setup transient 
                    // =====================================================
                } else if (base.Control.CompMode == AppControl._CompMode.Transient) {

                    // push stacks
                    // -----------

                    PushLevelSetAndRelatedStuff();


                    // backup old velocity for energy checks
                    // -------------------------------------
                    if (this.Control.ComputeEnergy && this.Control.CompMode == AppControl._CompMode.Transient) {
                        for (int d = 0; d < D; d++) {
                            this.prevVel[d].Clear();
                            this.prevVel[d].Acc(1.0, this.CurrentVel[d]);
                        }
                    }


                    // fields setup
                    // ------------
                    for (int d = 0; d < D; d++) {
                        // Gravity must be set up like this to avoid regions of zero gravity when updating the level-set
                        this.XDGvelocity.Gravity[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
                    }


                    // +++++++++++++++++++++++++++++++++++++
                    // compute/check time step restrictions
                    // +++++++++++++++++++++++++++++++++++++

                    dt = base.Control.dtFixed;

                    // Level-Set motion-CFL
                    double LevSet_Deg2 = this.DGLevSet.Current.Basis.Degree;
                    LevSet_Deg2 = LevSet_Deg2 * LevSet_Deg2;
                    double dt_LevSetCFL = base.GridData.ComputeCFLTime(this.ExtensionVelocity.Current, dt * LevSet_Deg2);
                    dt_LevSetCFL = dt_LevSetCFL / LevSet_Deg2;
                    if (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative && this.Control.LSunderrelax == 1.0) {
                        if (dt / dt_LevSetCFL > 1.0) {
                            double underrelax = Math.Round(dt_LevSetCFL / dt, 1);
                            m_BDF_Timestepper.IterUnderrelax = underrelax;
                            Console.WriteLine("Exceeding Level-Set CFL: Setting underrelaxation factor to {0}", underrelax);
                        } else {
                            m_BDF_Timestepper.IterUnderrelax = 1.0;
                        }
                    }


                    //dt = Math.Min(dt, dt_LevSetCFL);

                    // Capillary Timestep restriction
                    if (this.Control.PhysicalParameters.Sigma != 0.0) {
                        MultidimensionalArray h_mins = ((GridData)this.GridData).Cells.h_min;
                        double h = h_mins.Min();
                        double LevSet_Deg = this.LevSet.Basis.Degree + 1;
                        h /= LevSet_Deg;
                        double dt_sigma = Math.Sqrt((this.Control.PhysicalParameters.rho_A + this.Control.PhysicalParameters.rho_B)
                            * Math.Pow(h, 3) / (2 * Math.PI * this.Control.PhysicalParameters.Sigma));
                        if (dt > dt_sigma)
                            Console.WriteLine("Warning: exceeding Capillary timestep: dt = {0}, dt_sigma = {1}, frac = {2}", dt, dt_sigma, dt / dt_sigma);
                    }


                    // elo
                    // ---

                    Console.WriteLine("Instationary solve, timestep #{0}, dt = {1} ...", TimestepNo, dt);

                } else {
                    throw new NotImplementedException("Option " + base.Control.CompMode + " not supported yet.");
                }

                // =======================================================================
                // call timestepper
                // =======================================================================

                //if ((m_BDF_Timestepper == null) == (m_RK_Timestepper == null))
                //    throw new ApplicationException();

                //CurvatureAlgorithms.CurvatureDriver(
                //    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                //    CurvatureAlgorithms.FilterConfiguration.NoFilter,
                //    this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                //    this.m_HMForder, this.DGLevSet.Current);

                //double[] momBal_Norm = XNSEUtils.MomentumBalanceNormAtInterface(this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                //    this.Control.PhysicalParameters, this.Control.AdvancedDiscretizationOptions.SurfStressTensor, this.m_HMForder);

                //Console.WriteLine("x-momentum balance norm = {0}", momBal_Norm[0]);
                //Console.WriteLine("y-momentum balance norm = {0}", momBal_Norm[1]);


                // ++++++++++++++++++++++++++++++++++++++++++
                // The actual solution of the System
                // ++++++++++++++++++++++++++++++++++++++++++

                using (new BlockTrace("Solve", tr)) {

                    if (m_BDF_Timestepper != null) {
                        if (Control.RaiseWeissenberg == true) {

                            currentWeissenberg = 0.0;

                            if (Control.PhysicalParameters.Weissenberg_a != 0.0) {

                                if (Control.WeissenbergIncrement != 0.0) {
                                    NoIncrementTimestep = (int)(Control.PhysicalParameters.Weissenberg_a / Control.WeissenbergIncrement);
                                } else {
                                    throw new ArgumentException("Raise Weissenberg is turned on, but WeissenbergIncrement is zero!");
                                }

                            } else if (Control.PhysicalParameters.Weissenberg_b != 0.0) {

                                if (Control.WeissenbergIncrement != 0.0) {
                                    NoIncrementTimestep = (int)(Control.PhysicalParameters.Weissenberg_b / Control.WeissenbergIncrement);
                                } else {
                                    throw new ArgumentException("Raise Weissenberg is turned on, but WeissenbergIncrement is zero!");
                                }

                            } else {
                                throw new ArgumentException("Raise Weissenberg is turned on, but aim Weissenberg is 0.0 (Newtonian)!");
                            }


                            for (int i = 0; i <= NoIncrementTimestep; i++) {

                                if (Control.UseArtificialDiffusion == true) {


                                    artificialMaxViscosity = 1.0;

                                    for (int j = 0; j < 3; j++) {

                                        if (Control.UsePerssonSensor == true) {
                                            perssonsensor.Update(StressXX);
                                        } else {
                                            throw new ArgumentException("artificial viscosity is turned on, but Persson sensor is turned off!");
                                        }

                                        m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

                                        //this.ResLogger.NextTimestep(false);

                                        // this evaluation must later out of this loop. now here for comparing resluts with  
                                        PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                                        SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                                        artificialMaxViscosity = artificialMaxViscosity - 0.5;
                                    }
                                } else {
                                    m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

                                    //this.ResLogger.NextTimestep(false);

                                    // this evaluation must later out of this loop. now here for comparing resluts with  
                                    PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                                    SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                                }

                                ChangeMesh = Control.AdaptiveMeshRefinement;
                                while (ChangeMesh == true) {
                                    this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                                    perssonsensor.Update(StressXX);
                                    PlotCurrentState(phystime, TimestepNo);
                                }

                                if (currentWeissenberg < Control.PhysicalParameters.Weissenberg_a || currentWeissenberg < Control.PhysicalParameters.Weissenberg_b) {
                                    currentWeissenberg = currentWeissenberg + Control.WeissenbergIncrement;
                                    Console.WriteLine();
                                    Console.WriteLine("Raise Weissenberg number to " + currentWeissenberg);
                                    Console.WriteLine();
                                }

                            }
                        } else {
                            //current Weissenberg is set to the HIGHER value... DIRTY HACK AT THE MOMENT!

                            if (Control.PhysicalParameters.Weissenberg_b < Control.PhysicalParameters.Weissenberg_a) {
                                currentWeissenberg = Control.PhysicalParameters.Weissenberg_a;
                            } else {
                                currentWeissenberg = Control.PhysicalParameters.Weissenberg_b;
                            }

                            if (Control.UseArtificialDiffusion == true) {
                                artificialMaxViscosity = 1.0;

                                for (int j = 0; j < 3; j++) {

                                    if (Control.UsePerssonSensor == true) {
                                        perssonsensor.Update(StressXX);
                                    } else {
                                        throw new ArgumentException("artificial viscosity is turned on, but Persson sensor is turned off!");
                                    }

                                    m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

                                    // this evaluation must later out of this loop. now here for comparing resluts with  
                                    //PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                                    //SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                                    artificialMaxViscosity = artificialMaxViscosity - 0.5;
                                }

                                ChangeMesh = Control.AdaptiveMeshRefinement;
                                while (ChangeMesh == true) {
                                    this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                                }
                            } else {

                                m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
                            }
                        }
                    }

                }

                // NOT IMPLEMENTED!!!
                //Postprocessing(TimestepInt, phystime, dt, TimestepNo);

                //NOT YET IMPLEMENTED!!!
                //if (Control.ComputeL2Error == true) {
                //    this.ComputeL2Error();
                //}

                // ================
                // Good bye
                // ================
                if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
                    ExtVelMover.FinishTimeStep();
                }

#if DEBUG
            // in case of Debugging Save first Timesteps
            //if(TimestepNo[1] <= 2) {
            //    this.SaveToDatabase(TimestepNo, phystime);
            //}
#endif

                Console.WriteLine("done.");
                return dt;
            }
        }

        protected override void SetInitial() {
            base.SetInitial();

            this.InitLevelSet();

            this.CreateEquationsAndSolvers(null);

            m_BDF_Timestepper.SingleInit();
            VelocityXGradient = new VectorField<XDGField>(D, this.CurrentVel[0].Basis, "VelocityX_Gradient", XDGField.Factory);
            VelocityYGradient = new VectorField<XDGField>(D, this.CurrentVel[1].Basis, "VelocityY_Gradient", XDGField.Factory);

            if (this.Control.SetParamsAnalyticalSol == true) {
                U = new XDGField(new XDGBasis(LsTrk, this.CurrentVel[0].Basis.Degree), "UAnalytical");
                V = new XDGField(new XDGBasis(LsTrk, this.CurrentVel[1].Basis.Degree), "VAnalytical");
                U.ProjectField(this.Control.VelFunctionU);
                V.ProjectField(this.Control.VelFunctionV);

                VelocityXGradient.Clear();
                VelocityXGradient.Gradient(1.0, U);
                VelocityYGradient.Clear();
                VelocityYGradient.Gradient(1.0, V);
            }

            // =========================================
            // XDG BDF Timestepper initialization
            // =========================================

            if (m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DelayedTimestepperInit(0.0, 0, this.Control.GetFixedTimestep(),
                    // delegate for the initialization of previous timesteps from an analytic solution
                    BDFDelayedInitSetIntial);
            }
        }

        /// <summary>
        /// performs restart
        /// </summary>
        /// <param name="Time">
        /// on exit, the physical time associated with the field state
        /// </param>
        /// <param name="TimestepNo">
        /// on exit, the physical time associated with the field state
        /// </param>
        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            base.LoadRestart(out Time, out TimestepNo);

            //this.InitLevelSet();
            if (this.Control.ReInitPeriod > 0) {
                Console.WriteLine("FastMarchReInit performing FirstOrderReInit");
                FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                CellMask Accepted = LsTrk.Regions.GetCutCellMask();
                CellMask ActiveField = LsTrk.Regions.GetNearFieldMask(1);
                CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);
            }

            ContinuityEnforcer = new ContinuityProjection(
                    ContBasis: this.LevSet.Basis,
                    DGBasis: this.DGLevSet.Current.Basis,
                    gridData: GridData,
                    Option: Control.LSContiProjectionMethod
                    );

            if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
                ExtVelMover = new ExtensionVelocityBDFMover(LsTrk, DGLevSet.Current, DGLevSetGradient, new VectorField<DGField>(XDGvelocity.Velocity.ToArray()),
                    Control.EllipticExtVelAlgoControl, BcMap, bdfOrder, ExtensionVelocity.Current, new double[2] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B });
            }

            //this.LsTrk.UpdateTracker();

            this.CreateEquationsAndSolvers(null);

            // =========================================
            // XDG BDF Timestepper initialization
            // =========================================

            if (m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DelayedTimestepperInit(Time, TimestepNo.MajorNumber, this.Control.GetFixedTimestep(),
                    // delegate for the initialization of previous timesteps from restart session
                    BDFDelayedInitLoadRestart);
            }

        }

        /// <summary>
        /// overriding the method to implement any user-specific tasks which
        /// should be carried out after a restart file has been loaded (e.g.,
        /// setting the correct time for a time-stepper)
        /// </summary>
        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            VelocityXGradient = new VectorField<XDGField>(this.GridData.SpatialDimension, this.CurrentVel[0].Basis, "VelocityX_Gradient", XDGField.Factory);
            VelocityYGradient = new VectorField<XDGField>(this.GridData.SpatialDimension, this.CurrentVel[1].Basis, "VelocityY_Gradient", XDGField.Factory);

            if (this.Control.FourierLevSetControl != null) {

                Guid sessionGuid = this.Control.RestartInfo.Item1;
                string FourierLog_path = this.Control.DbPath;
                FourierLog_path = FourierLog_path + "/sessions/" + sessionGuid + "/Log_FourierLS.txt";

                IList<Guid> spUids = new List<Guid>();
                try {
                    using (StreamReader samplPLogReader = new StreamReader(FourierLog_path)) {

                        while (!samplPLogReader.EndOfStream) {
                            spUids.Add(Guid.Parse(samplPLogReader.ReadLine()));
                        }
                    }
                } catch (FileNotFoundException) {
                    spUids = new Guid[0];
                }
                Guid[] samplPUids = spUids.ToArray();

                TimestepNumber tsNmbr = this.Control.RestartInfo.Item2;
                Guid samplPUid_restart = Guid.Empty;
                if (tsNmbr == null) {
                    samplPUid_restart = samplPUids[samplPUids.Length - 1];
                } else {
                    samplPUid_restart = samplPUids[tsNmbr.MajorNumber];
                }
                Partitioning p = null;
                double[] samplP = base.DatabaseDriver.LoadVector<double>(samplPUid_restart, ref p).ToArray();

                if (this.Control.FourierLevSetControl.FType == FourierType.Polar) {
                    double[] center = samplP.GetSubVector(0, 2);
                    samplP = samplP.GetSubVector(2, samplP.Length - 2);
                    this.Control.FourierLevSetControl.center = center;
                    this.Control.FourierLevSetControl.samplP = samplP;
                } else {
                    this.Control.FourierLevSetControl.samplP = samplP;
                }
            }
        }

        private void BDFDelayedInitSetIntial(int TimestepIndex, double Time, DGField[] St) {
            using (new FuncTrace()) {
                Console.WriteLine("Timestep index {0}, time {1} ", TimestepIndex, Time);

                // level-set
                // ---------
                this.DGLevSet.Current.ProjectField(X => this.Control.Phi(X, Time));
                this.LevSet.ProjectField(X => this.Control.Phi(X, Time));

                this.LsTrk.UpdateTracker(incremental: true);

                // solution
                // --------
                int D = this.LsTrk.GridDat.SpatialDimension;

                for (int d = 0; d < D; d++) {
                    XDGField _u = (XDGField)St[d];
                    _u.Clear();
                    _u.GetSpeciesShadowField("A").ProjectField(X => this.Control.ExactSolutionVelocity["A"][d](X, Time));
                    _u.GetSpeciesShadowField("B").ProjectField((X => this.Control.ExactSolutionVelocity["B"][d](X, Time)));
                }
                XDGField _p = (XDGField)St[D];
                _p.Clear();
                _p.GetSpeciesShadowField("A").ProjectField(X => this.Control.ExactSolutionPressure["A"](X, Time));
                _p.GetSpeciesShadowField("B").ProjectField((X => this.Control.ExactSolutionPressure["B"](X, Time)));
            }
        }

        private void BDFDelayedInitLoadRestart(int TimestepIndex, double time, DGField[] St) {

            Console.WriteLine("Timestep index {0}, time {1} ", TimestepIndex, time);

            ITimestepInfo tsi_toLoad;
            if (TimestepIndex < 0) {
                throw new ArgumentOutOfRangeException("Not enough Timesteps to restart with desired Timestepper");
            } else {
                ISessionInfo reloadSession = GetDatabase().Controller.GetSessionInfo(this.CurrentSessionInfo.RestartedFrom);
                tsi_toLoad = reloadSession.Timesteps.Single(t => t.TimeStepNumber.Equals(new TimestepNumber(TimestepIndex)));
            }
            DatabaseDriver.LoadFieldData(tsi_toLoad, this.GridData, this.IOFields);

            // level-set
            // ---------
            this.DGLevSet.Current.Clear();
            this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);

            //this.LsTrk.UpdateTracker(incremental: true);

            // solution
            // --------
            int D = this.LsTrk.GridDat.SpatialDimension;

            for (int d = 0; d < D; d++) {
                St[d] = this.XDGvelocity.Velocity[d].CloneAs();
            }
            St[D] = this.Pressure.CloneAs();
        }

        /// <summary>
        /// 
        /// </summary>
        protected override ITimestepInfo SaveToDatabase(TimestepNumber timestepno, double t) {
            var tsi = base.SaveToDatabase(timestepno, t);

            if (tsi != null && m_BDF_Timestepper != null) {
                int S = m_BDF_Timestepper.GetNumberOfStages;

                SinglePhaseField LsBkUp = new SinglePhaseField(this.LevSet.Basis);
                LsBkUp.Acc(1.0, this.LevSet);

                ICollection<DGField>[] restartFields = m_BDF_Timestepper.GetRestartInfos();

                if (S > 1 && this.Control.saveperiod >= S && restartFields != null) {

                    // save additional timesteps/information for restart
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    for (int ti = 1; ti < S; ti++) {

                        //SinglePhaseField LsBkUp = new SinglePhaseField(this.LevSet.Basis);
                        //LsBkUp.Acc(1.0, this.LevSet);

                        ICollection<DGField> restartIOFields = new List<DGField>();
                        foreach (DGField f in this.IOFields) {

                            int rfidx = restartFields[ti - 1].IndexWhere(rf => rf.Identification == f.Identification);
                            if (rfidx > -1) {
                                DGField rf = restartFields[ti - 1].ElementAt(rfidx);
                                if (f.Identification == "Phi") {
                                    this.LevSet.Clear();
                                    this.LevSet.Acc(1.0, rf);
                                    restartIOFields.Add(this.LevSet);
                                } else {
                                    restartIOFields.Add(rf);
                                }
                            } else {
                                DGField rf = f.CloneAs();
                                rf.Clear();
                                restartIOFields.Add(rf);
                            }
                        }

                        //this.LevSet.Clear();
                        //this.LevSet.Acc(1.0, LsBkUp);

                        ITimestepInfo rtsi;
                        TimestepNumber tsn = new TimestepNumber(timestepno.MajorNumber - ti);

                        //Console.WriteLine("saving to Database ...");

                        //Exception e = null;
                        try {
                            rtsi = new TimestepInfo(
                                t - ti * this.Control.GetFixedTimestep(),
                                this.CurrentSessionInfo,
                                tsn,
                                restartIOFields);
                        } catch (Exception ee) {
                            Console.Error.WriteLine(ee.GetType().Name + " on rank " + this.MPIRank + " saving time-step " + tsn + ": " + ee.Message);
                            Console.Error.WriteLine(ee.StackTrace);
                            //tsi = null;
                            //e = ee;

                            if (ContinueOnIOError) {
                                Console.WriteLine("Ignoring IO error: " + DateTime.Now);
                            } else {
                                throw ee;
                            }

                            tsi = null;
                        }
                        // e.ExceptionBcast();
                        csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                    }

                }

                this.LevSet.Clear();
                this.LevSet.Acc(1.0, LsBkUp);
            }
            return tsi;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {
            Tecplot.PlotFields(base.m_RegisteredFields, "XRheology_Solver" + timestepNo, physTime, superSampling);
            //Tecplot.PlotFields(new DGField[] { this.LevSet }, "grid" + timestepNo, physTime, 0);
        }

        protected void PlotOnIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, iterIndex }), 2);
        }

        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int pVel = this.CurrentVel[0].Basis.Degree;
                int pPrs = this.Pressure.Basis.Degree;
                int D = this.GridData.SpatialDimension;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 1];

                    // configurations for velocity
                    for (int d = 0; d < D; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            Degree = Math.Max(1, pVel - iLevel),
                            mode = this.Control.VelocityBlockPrecondMode,
                            VarIndex = new int[] { d }
                        };
                    }
                    // configuration for pressure
                    configs[iLevel][D] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, pPrs - iLevel),
                        mode = this.Control.PressureBlockPrecondMode,
                        VarIndex = new int[] { D }
                    };
                }


                return configs;
            }
        }

        internal double[] ComputeL2Error(double time) {
            int D = this.GridData.SpatialDimension;
            double[] Ret = new double[D + 1];

            if (this.Control.ExactSolutionVelocity == null && this.Control.ExactSolutionPressure == null)
                // nothing to do
                return Ret;

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;


            // Velocity error
            // ==============
            if (this.Control.ExactSolutionVelocity != null) {
                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[D];

                foreach (var spc in this.LsTrk.SpeciesNames) {
                    L2Error_Species.Add(spc, new double[D]);

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                    for (int d = 0; d < D; d++) {
                        ConventionalDGField Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);

                        L2Error_Species[spc][d] = Vel_d.L2Error(this.Control.ExactSolutionVelocity[spc][d].Vectorize(time), order, scheme);
                        L2Error[d] += L2Error_Species[spc][d].Pow2();

                        base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);
                    }
                }
                L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

                for (int d = 0; d < D; d++) {
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                    Ret[d] = L2Error[d];
                }

            }


            // pressure error
            // ==============
            if (this.Control.ExactSolutionPressure != null) {

                // pass 1: mean value of pressure difference
                double DiffInt = 0;
                foreach (var spc in this.LsTrk.SpeciesNames) {

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.GridData, order);

                    DiffInt += this.Pressure.GetSpeciesShadowField(spc).LxError(this.Control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b), rule);
                    //Volume +=  this.Pressure.GetSpeciesShadowField(spc).LxError(null, (a, b) => (1.0), rule);
                }
                double Volume2 = (new SubGrid(CellMask.GetFullMask(this.GridData))).Volume;
                double PressureDiffMean = DiffInt / Volume2;


                double L2Error = 0;
                Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                foreach (var spc in this.LsTrk.SpeciesNames) {

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.GridData, order);

                    double IdV = this.Pressure.GetSpeciesShadowField(spc).LxError(this.Control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b - PressureDiffMean).Pow2(), rule);
                    L2Error += IdV;
                    L2Error_Species.Add(spc, IdV.Sqrt());

                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure + "#" + spc, L2Error_Species[spc], true);
                }


                L2Error = L2Error.Sqrt();
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
                Ret[D] = L2Error;

            }

            // Stress error
            // =============================================================
            if (this.Control.ExSol_Stress != null) {

                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[3];


                //must be changed for multiphase!!!
                foreach (var spc in this.LsTrk.SpeciesNames) {
                    L2Error_Species.Add(spc, new double[3]);

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                    //stressXX
                    //ConventionalDGField Stress_XX = ((XDGField)this.StressXX).GetSpeciesShadowField(spc);
                    //L2Error_Species[spc] = Stress_XX.L2Error(this.Control.ExSol_Stress[spc].Vectorize(0.0), order, scheme);
                    //L2Error[0] += L2Error_Species[spc].Pow2();

                    //base.QueryHandler.ValueQuery("L2err_" + VariableNames.Stress_XX + "#" + spc, L2Error_Species[spc], true);

                    //______________________________________________________________________________________________________________________
                    for (int d = 0; d < D; d++) {
                        ConventionalDGField Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);

                        L2Error_Species[spc][d] = Vel_d.L2Error(this.Control.ExactSolutionVelocity[spc][d].Vectorize(time), order, scheme);
                        L2Error[d] += L2Error_Species[spc][d].Pow2();

                        base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);
                    }
                }
                L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

                for (int d = 0; d < D; d++) {
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                    Ret[d] = L2Error[d];
                }
                //______________________________________________________________________________________________________________________



                L2Error[0] = this.StressXX.L2Error(this.Control.ExSol_Stress[0].Vectorize(0.0), order);
                L2Error[1] = this.StressXY.L2Error(this.Control.ExSol_Stress[1].Vectorize(0.0), order);
                L2Error[2] = this.StressYY.L2Error(this.Control.ExSol_Stress[2].Vectorize(0.0), order);

                base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressXX, L2Error[0], true);
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressXY, L2Error[1], true);
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressYY, L2Error[2], true);

                Console.WriteLine("L2err " + VariableNames.StressXX + " is " + L2Error[0]);
                Console.WriteLine("L2err " + VariableNames.StressXY + " is " + L2Error[1]);
                Console.WriteLine("L2err " + VariableNames.StressYY + " is " + L2Error[2]);
            }

            return Ret;
        }

        // =========================
        // adaptive mesh refinement
        // ========================
        #region AMR

        CellMask NScm;

        CellMask NSbuffer;

        /// <summary>
        /// refinement indicator for a constant near band refinement
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel) {

            if (this.Control.BaseRefinementLevel == 0)
                return 0;

            CellMask ccm = this.LsTrk.Regions.GetCutCellMask();
            CellMask near = this.LsTrk.Regions.GetNearFieldMask(1);
            CellMask nearBnd = near.AllNeighbourCells();
            CellMask buffer = nearBnd.AllNeighbourCells().Union(nearBnd).Except(near);


            int DesiredLevel_j = CurrentLevel;

            if (near.Contains(j)) {
                if (CurrentLevel < this.Control.BaseRefinementLevel) {
                    DesiredLevel_j++;
                } else {
                    // additional refinement
                    switch (this.Control.RefineStrategy) {
                        case XNSE_Control.RefinementStrategy.CurvatureRefined: {
                                double curv_max = 1.0 / (2.0 * ((GridData)this.GridData).Cells.h_min[j]);
                                double mean_curv = Math.Abs(this.Curvature.GetMeanValue(j));
                                double minCurv, maxCurv;
                                this.Curvature.GetExtremalValuesInCell(out minCurv, out maxCurv, j);
                                double max_AbsCurv = Math.Max(Math.Abs(minCurv), Math.Abs(maxCurv));

                                double curv_thrshld = mean_curv;
                                if (curv_thrshld > curv_max && CurrentLevel == this.Control.RefinementLevel) {
                                    DesiredLevel_j++;
                                } else if (curv_thrshld < (curv_max / 2) && CurrentLevel == this.Control.RefinementLevel + 1) {
                                    DesiredLevel_j--;
                                }
                                break;
                            }
                        case XNSE_Control.RefinementStrategy.ContactLineRefined: {
                                CellMask BCells = ((GridData)this.GridData).BoundaryCells.VolumeMask;
                                if (ccm.Contains(j) && BCells.Contains(j) && CurrentLevel < this.Control.RefinementLevel) {
                                    DesiredLevel_j++;
                                } else if (!BCells.Contains(j)) { // && CurrentLevel == this.Control.RefinementLevel + 1) {
                                    DesiredLevel_j--;
                                }
                                break;
                            }
                        case XNSE_Control.RefinementStrategy.constantInterface:
                        default:
                            break;
                    }
                }

            } else if (NScm.Contains(j)) {
                if (CurrentLevel < this.Control.BaseRefinementLevel)
                    DesiredLevel_j++;

            } else if (buffer.Contains(j) || NSbuffer.Contains(j)) {
                if (CurrentLevel < this.Control.BaseRefinementLevel - 1)
                    DesiredLevel_j++;
            } else {
                DesiredLevel_j = 0;
            }

            return DesiredLevel_j;

        }

        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            using (new FuncTrace()) {

                if (this.Control.AdaptiveMeshRefinement) {

                    //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 0), 2);

                    // Check grid changes
                    // ==================

                    CellMask BlockedCells;
                    if (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Once
                        || this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                        int prevInd = LsTrk.PopulatedHistoryLength;
                        CellMask prevNear = LsTrk.RegionsHistory[-prevInd + 1].GetNearFieldMask(1);
                        BlockedCells = (TimestepNo > 1) ? prevNear : null; // CellMask.Union(currNear, prevNear);
                    } else {
                        CellMask currNear = LsTrk.Regions.GetNearFieldMask(1);
                        BlockedCells = currNear;
                    }

                    // compute curvature for levelindicator 
                    if (this.Control.RefineStrategy == XNSE_Control.RefinementStrategy.CurvatureRefined) {
                        CurvatureAlgorithms.CurvatureDriver(
                            SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint,
                            CurvatureAlgorithms.FilterConfiguration.Default,
                            this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                            this.m_HMForder, this.DGLevSet.Current);
                    }


                    // navier slip boundary cells
                    NScm = new CellMask(this.GridData);
                    NSbuffer = new CellMask(this.GridData);
                    if (this.Control.RefineNavierSlipBoundary) {
                        BitArray NSc = new BitArray(((GridData)this.GridData).Cells.Count);
                        CellMask bnd = ((GridData)this.GridData).BoundaryCells.VolumeMask;
                        int[][] c2e = ((GridData)this.GridData).Cells.Cells2Edges;
                        foreach (Chunk cnk in bnd) {
                            for (int i = cnk.i0; i < cnk.JE; i++) {
                                foreach (int e in c2e[i]) {
                                    int eId = (e < 0) ? -e - 1 : e - 1;
                                    byte et = ((GridData)this.GridData).Edges.EdgeTags[eId];
                                    if (this.GridData.EdgeTagNames[et].Contains("navierslip_linear"))
                                        NSc[i] = true;
                                }
                            }
                        }
                        NScm = new CellMask(this.GridData, NSc);
                        CellMask bndNScm = NScm.AllNeighbourCells();
                        int bndLvl = 2;
                        for (int lvl = 1; lvl < bndLvl; lvl++) {
                            NScm = NScm.Union(bndNScm);
                            bndNScm = NScm.AllNeighbourCells();
                            NSbuffer = NScm.Union(bndNScm);
                        }
                        NSbuffer = NSbuffer.AllNeighbourCells();
                    }


                    //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 1), 2);


                    bool AnyChange = GridRefinementController.ComputeGridChange((BoSSS.Foundation.Grid.Classic.GridData)this.GridData, BlockedCells, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                    int NoOfCellsToRefine = 0;
                    int NoOfCellsToCoarsen = 0;
                    if (AnyChange) {
                        int[] glb = (new int[] {

                    CellsToRefineList.Count,
                    Coarsening.Sum(L => L.Length),
                }).MPISum();
                        NoOfCellsToRefine = glb[0];
                        NoOfCellsToCoarsen = glb[1];
                    }
                    int oldJ = this.GridData.CellPartitioning.TotalLength;

                    // Update Grid
                    // ===========

                    if (AnyChange) {

                        //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 1 }), 2);

                        Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                        Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                        newGrid = ((GridData)this.GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                        //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 2 }), 2);


                    } else {

                        newGrid = null;
                        old2NewGrid = null;
                    }

                } else {

                    newGrid = null;
                    old2NewGrid = null;
                }

            }
        }

        #endregion
    }
}