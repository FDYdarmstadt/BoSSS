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

using BoSSS.Platform;

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
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.Grid.Aggregation;
using NUnit.Framework;
using MPI.Wrappers;
using BoSSS.Solution.LevelSetTools.Advection;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Solver for Incompressible Multiphase flows
    /// </summary>
    public class XNSE_SolverMain : BoSSS.Solution.Application<XNSE_Control> {



        static void Main(string[] args) {
        
            _Main(args, false, delegate () {
                var p = new XNSE_SolverMain();
                return p;
            });
        }

        // Instantiate Fields from Control File

        #region instantiation
#pragma warning disable 649
        /// <summary>
        /// Pressure
        /// </summary>
        //[InstantiateFromControlFile(VariableNames.Pressure, null, IOListOption.ControlFileDetermined)]
        XDGField Pressure;

        /// <summary>
        /// Artificial force term at the fluid interface, usually only to support manufactured solutions.
        /// </summary>
        [InstantiateFromControlFile(
            new string[] { "SurfaceForceX", "SurfaceForceY", "SurfaceForceZ" },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        VectorField<SinglePhaseField> SurfaceForce;

        /// <summary>
        /// the DG representation of the level set.
        /// This one is used for level-set evolution in time; it is in general discontinuous.
        /// </summary>
        //[InstantiateFromControlFile("PhiDG", "Phi", IOListOption.ControlFileDetermined)]
        ScalarFieldHistory<SinglePhaseField> DGLevSet;

        /// <summary>
        /// The continuous level set field which defines the XDG space; 
        /// it is obtained from the projection of the discontinuous <see cref="DGLevSet"/> onto the 
        /// continuous element space.
        /// </summary>
        //[InstantiateFromControlFile("Phi", "Phi", IOListOption.ControlFileDetermined)]
        LevelSet LevSet;

        /// <summary>
        /// Guess what?
        /// </summary>
        VectorField<SinglePhaseField> DGLevSetGradient;

        VectorField<SinglePhaseField> LevSetGradient;

        /// <summary>
        /// Curvature; DG-polynomial degree should be 2 times the polynomial degree of <see cref="LevSet"/>.
        /// </summary>
        [InstantiateFromControlFile("Curvature", "Curvature", IOListOption.ControlFileDetermined)]
        SinglePhaseField Curvature;

        /// <summary>
        /// Residual of the continuity equation
        /// </summary>
        //[InstantiateFromControlFile("ResidualConti", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        XDGField ResidualContinuity;

        /// <summary>
        /// Divergence of velocity ->
        /// Conservation of mass for incompressibility
        /// </summary>
        //XDGField divVelocity;


        SinglePhaseField MassBalanceAtInterface;

        VectorField<SinglePhaseField> MomentumBalanceAtInterface;

        SinglePhaseField EnergyBalanceAtInterface;



        /// <summary>
        /// If requested, performs the projection of the level-set on a continuous field
        /// </summary>
        ContinuityProjection ContinuityEnforcer;

        /// <summary>
        /// Lauritz' Fast Marching Solver
        /// !!! Caution !!! Only works in Single-Core
        /// </summary>
        FastMarchReinit FastMarchReinitSolver;

        /// <summary>
        /// PDE based elliptic reInitialization by Thomas
        /// </summary>
        EllipticReInit ReInitPDE;

        /// <summary>
        /// Bundling of variables which are either DG or XDG (see <see cref="XNSE_Control.UseXDG4Velocity"/>);
        /// </summary>
        class VelocityRelatedVars<TX> where TX : DGField {
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
        /// Velocity and related variables for the non-extended case, <see cref="XNSE_Control.UseXDG4Velocity"/> == false.
        /// </summary>
        VelocityRelatedVars<SinglePhaseField> DGvelocity;

        /// <summary>
        /// Velocity and related variables for the extended case, <see cref="XNSE_Control.UseXDG4Velocity"/> == false.
        /// </summary>
        VelocityRelatedVars<XDGField> XDGvelocity;

        //// <summary>
        //// Continuous high-order finite element basis for the level-set in toder to ensure continuity,
        //// see <see cref="XNSE_Control.EnforceLevelSetContinuity"/>.
        //// </summary>
        //SpecFemBasis ContinuousLevelSetBasis;

        //Basis ContinuousLevelSetDGBasis;

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

#pragma warning restore 649

        protected override void CreateFields() {
            using (new FuncTrace()) {
                base.CreateFields();
                int D = this.GridData.SpatialDimension;


                this.DGLevSet = new ScalarFieldHistory<SinglePhaseField>(
                       new SinglePhaseField(new Basis(this.GridData, this.Control.FieldOptions["Phi"].Degree), "PhiDG"));

                if (this.Control.FieldOptions["PhiDG"].Degree >= 0 && this.Control.FieldOptions["PhiDG"].Degree != this.DGLevSet.Current.Basis.Degree) {
                    throw new ApplicationException("Specification of polynomial degree for 'PhiDG' is not supportet, since it is induced by polynomial degree of 'Phi'.");
                }

                // ==============================
                // Initialize ContinuityProjection
                // if needed, if not , Option: None
                // ==============================
                this.LevSet = ContinuityProjection.CreateField(
                    DGLevelSet: this.DGLevSet.Current,
                    gridData: GridData,
                    Option: Control.LSContiProjectionMethod
                    );

                this.LsTrk = new LevelSetTracker(this.GridData, base.Control.CutCellQuadratureType, base.Control.LS_TrackerWidth, new string[] { "A", "B" }, this.LevSet);
                base.RegisterField(this.LevSet);
                this.LevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(this.LevSet.Basis, "dPhi_dx[" + d + "]")));
                base.RegisterField(this.LevSetGradient);

                base.RegisterField(this.DGLevSet.Current);
                this.DGLevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(this.DGLevSet.Current.Basis, "dPhiDG_dx[" + d + "]")));
                base.RegisterField(this.DGLevSetGradient);

                this.Pressure = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree), VariableNames.Pressure);
                base.RegisterField(this.Pressure);
                this.ResidualContinuity = new XDGField(this.Pressure.Basis, "ResidualConti");
                base.RegisterField(this.ResidualContinuity);
        

                if (base.Control.UseXDG4Velocity) {

                    XDGvelocity = new VelocityRelatedVars<XDGField>();
                    InitFromAttributes.CreateFieldsAuto(XDGvelocity, this.GridData, base.Control.FieldOptions, base.Control.CutCellQuadratureType, base.IOFields, base.m_RegisteredFields);

                } else {

                    DGvelocity = new VelocityRelatedVars<SinglePhaseField>();
                    InitFromAttributes.CreateFieldsAuto(DGvelocity, this.GridData, base.Control.FieldOptions, base.Control.CutCellQuadratureType, base.IOFields, base.m_RegisteredFields);

                    //if(base.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold > 0.0)
                    //    throw new NotSupportedException("Agglomearion currently not supported for non-extended velocity.");
                }

                //XDGBasis b = new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree);
                //this.divVelocity = new XDGField(b, "DivergenceVelocity");
                //base.RegisterField(this.divVelocity);

                if (this.Control.ComputeEnergy && this.Control.CompMode == AppControl._CompMode.Transient) {
                    prevVel = new XDGField[D];
                    for (int d = 0; d < D; d++) {
                        prevVel[d] = new XDGField(this.XDGvelocity.Velocity[d].Basis);
                    }
                }

                if (this.Control.CheckJumpConditions) {
                    Basis basis = new Basis(this.GridData, 0);

                    //Basis basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.VelocityX].Degree);
                    this.MassBalanceAtInterface = new SinglePhaseField(basis, "MassBalanceAtInterface");
                    base.RegisterField(this.MassBalanceAtInterface);

                    //basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree + this.Control.FieldOptions["Phi"].Degree);
                    this.MomentumBalanceAtInterface = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(basis, d + "-MomentumBalanceAtInterface")));
                    base.RegisterField(this.MomentumBalanceAtInterface);

                    //basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree + this.Control.FieldOptions[VariableNames.VelocityX].Degree + this.Control.FieldOptions["Phi"].Degree);
                    this.EnergyBalanceAtInterface = new SinglePhaseField(basis, "EnergyBalanceAtInterface");
                    base.RegisterField(this.EnergyBalanceAtInterface);
                }

            }
        }

        #endregion
        

        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double rho_A = this.Control.PhysicalParameters.rho_A,
                    rho_B = this.Control.PhysicalParameters.rho_B;

                int D = this.GridData.SpatialDimension;

                double[] _rho_A = new double[D + 1];
                _rho_A.SetAll(rho_A); // mass matrix in momentum equation
                _rho_A[D] = 0; // no  mass matrix for continuity equation
                double[] _rho_B = new double[D + 1];
                _rho_B.SetAll(rho_B); // mass matrix in momentum equation
                _rho_B[D] = 0; // no  mass matrix for continuity equation

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _rho_A);
                R.Add(this.LsTrk.GetSpeciesId("B"), _rho_B);

                return R;
            }
        }

        IncompressibleMultiphaseBoundaryCondMap m_BcMap;

        /// <summary>
        /// Boundary conditions.
        /// </summary>
        IncompressibleMultiphaseBoundaryCondMap BcMap {
            get {
                if (m_BcMap == null) {
                    m_BcMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
                }
                return m_BcMap;
            }
        }

        /// <summary>
        /// the spatial operator (momentum and continuity equation)
        /// </summary>
        OperatorFactory XNSE_Operator;

        /// <summary>
        /// OperatorConfiguration for the <see cref="XNSE_Operator"/>
        /// </summary>
        OperatorConfiguration XOpConfig;

        /// <summary>
        /// Current Velocity: either extended or non-extended DG.
        /// </summary>
        DGField[] CurrentVel {
            get {
                Debug.Assert((this.XDGvelocity != null) != (this.DGvelocity != null));
                Debug.Assert((this.XDGvelocity != null) == base.Control.UseXDG4Velocity);
                Debug.Assert((this.DGvelocity == null) == base.Control.UseXDG4Velocity);
                if (base.Control.UseXDG4Velocity)
                    return this.XDGvelocity.Velocity.ToArray();
                else
                    return this.DGvelocity.Velocity.ToArray();
            }
        }

        DGField[] prevVel;

        CoordinateVector m_CurrentSolution;

        /// <summary>
        /// Current velocity and pressure;
        /// </summary>
        internal CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(this.CurrentVel, this.Pressure));
                } else {
                    for (int d = 0; d < base.GridData.SpatialDimension; d++) {
                        Debug.Assert(object.ReferenceEquals(m_CurrentSolution.Mapping.Fields[d], this.CurrentVel[d]));
                    }
                }

                return m_CurrentSolution;
            }
        }

        CoordinateVector m_CurrentResidual;

        /// <summary>
        /// Current residual for momentum and continuity equation.
        /// </summary>
        internal CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    if (base.Control.UseXDG4Velocity)
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat<DGField>(XDGvelocity.ResidualMomentum, ResidualContinuity));
                    else
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat<DGField>(DGvelocity.ResidualMomentum, ResidualContinuity));
                }
                return m_CurrentResidual;
            }
        }

        /*
        CoordinateVector m_PreviousSolution;

        /// <summary>
        /// In instationary calculations, the previous timestep, i.e. 
        /// the initial condition for the current timestep.
        /// </summary>
        CoordinateVector PreviousSolution {
            get {
                //if (m_PreviousSolution == null) {
                if (base.Control.UseXDG4Velocity)
                    m_PreviousSolution = new CoordinateVector(ArrayTools.Cat(this.XDGvelocity.Velocity, this.Pressure));
                else
                    m_PreviousSolution = new CoordinateVector(ArrayTools.Cat(this.DGvelocity.Velocity, this.Pressure));
                //} else {
                //    for (int d = 0; d < base.GridData.SpatialDimension; d++) {
                //        Debug.Assert(object.ReferenceEquals(m_PreviousSolution.Mapping.Fields[d], this.XDGvelocity.Velocity[0][d]));
                //    }
                //}
                return m_PreviousSolution;
            }
        }
        */

        /// <summary>
        /// output of <see cref="AssembleMatrix"/>;
        /// </summary>
        MassMatrixFactory MassFact;

        /// <summary>
        /// HMF order/degree which is used globally in this solver.
        /// </summary>
        int m_HMForder;

        /// <summary>
        /// Implicit timestepping using Backward-Differentiation-Formulas (BDF),
        /// specialized for XDG applications.
        /// </summary>
        XdgBDFTimestepping m_BDF_Timestepper;

        
        ///// <summary>
        ///// Explicit or implicit timestepping using Runge-Kutta formulas,
        ///// specialized for XDG applications.
        ///// </summary>
        //XdgRKTimestepping m_RK_Timestepper;

        RungeKuttaScheme rksch = null;
        int bdfOrder = -1000;




        /// <summary>
        /// Create XOperator and Timestepper
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            #region Checks
            // CreateEquationsAndSolvers might be called multiple times
            // exit if so, and no LoadBalancing
            if (XNSE_Operator != null && L == null)
                return;


            if (Control.CompMode == AppControl._CompMode.Steady) {
                if (Control.Timestepper_LevelSetHandling != LevelSetHandling.None)
                    throw new ApplicationException(string.Format("Illegal control file: for a steady computation ({0}), the level set handling must be {1}.", AppControl._CompMode.Steady, LevelSetHandling.None));
            }

            int degU = this.CurrentVel[0].Basis.Degree;

            if (base.Control.FakePoisson) {
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("ACHTUNG: Fake-Poisson aktiviert!");
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            }
            #endregion

            #region Config and Generate XOperator

            XOpConfig = new OperatorConfiguration() {
                continuity = true,
                Viscous = !base.Control.FakePoisson,
                PressureGradient = true,
                Transport = !base.Control.FakePoisson,
                CodBlocks = new bool[] { true, true },
                DomBlocks = new bool[] { true, true },
                dntParams = this.Control.AdvancedDiscretizationOptions,
                physParams = this.Control.PhysicalParameters,
                UseXDG4Velocity = this.Control.UseXDG4Velocity
            };

            //Quadrature Order
            //----------------

            m_HMForder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);


            // Is Moving Mesh required?
            //-------------------------

            bool movingmesh;
            MassMatrixShapeandDependence mmsd;
            switch (this.Control.Timestepper_LevelSetHandling) {
                case LevelSetHandling.Coupled_Once:
                movingmesh = true;
                mmsd = MassMatrixShapeandDependence.IsTimeDependent;
                break;

                case LevelSetHandling.Coupled_Iterative:
                movingmesh = true;
                mmsd = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                break;

                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.StrangSplitting:
                movingmesh = false;
                mmsd = MassMatrixShapeandDependence.IsTimeDependent;
                break;

                case LevelSetHandling.None:
                movingmesh = false;
                mmsd = MassMatrixShapeandDependence.IsNonIdentity;
                break;

                default:
                throw new NotImplementedException();
            }


            // Create Spatial Operator
            // ======================= 

            XNSE_Operator = new OperatorFactory(
               XOpConfig,
               this.LsTrk,
               this.m_HMForder,
               degU,
               this.BcMap,
               movingmesh);
            #endregion

            #region Create Timestepper
            // ==================
            if (L == null) {

                switch (this.Control.Timestepper_Scheme) {
                    case XNSE_Control.TimesteppingScheme.RK_ImplicitEuler: {
                        rksch = RungeKuttaScheme.ImplicitEuler;
                        break;
                    }
                    case XNSE_Control.TimesteppingScheme.RK_CrankNicolson: {
                        rksch = RungeKuttaScheme.CrankNicolson;
                        break;
                    }
                    case XNSE_Control.TimesteppingScheme.CrankNicolson: {
                        //do not instantiate rksch, use bdf instead
                        bdfOrder = -1;
                        break;
                    }
                    case XNSE_Control.TimesteppingScheme.ImplicitEuler: {
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
                            throw new NotImplementedException();
                    }

                }


                //if no Runge Kutta Timesteper is initialized
                // we are using bdf or screwed things up
                if (rksch == null) {
                    m_BDF_Timestepper = new XdgBDFTimestepping(
                        this.CurrentSolution.Mapping.Fields,
                        this.CurrentResidual.Mapping.Fields,
                        LsTrk,
                        true,
                        DelComputeOperatorMatrix, DelUpdateLevelSet,
                        (this.Control.CompMode == AppControl._CompMode.Transient) ? bdfOrder : 1,
                        this.Control.Timestepper_LevelSetHandling,
                        mmsd,
                        (this.Control.PhysicalParameters.IncludeConvection) ? SpatialOperatorType.Nonlinear : SpatialOperatorType.LinearTimeDependent,
                        MassScale,
                        this.MultigridOperatorConfig, base.MultigridSequence,
                        this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder,
                        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                        true
                        );
                    m_BDF_Timestepper.m_ResLogger = base.ResLogger;
                    m_BDF_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                    m_BDF_Timestepper.Config_SolverConvergenceCriterion = this.Control.Solver_ConvergenceCriterion;
                    m_BDF_Timestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                    m_BDF_Timestepper.Config_MaxIterations = this.Control.Solver_MaxIterations;
                    m_BDF_Timestepper.Config_MinIterations = this.Control.Solver_MinIterations;
                    m_BDF_Timestepper.Timestepper_Init = (this.Control.CompMode == AppControl._CompMode.Transient) ? this.Control.Timestepper_BDFinit : TimeStepperInit.SingleInit;
                    m_BDF_Timestepper.incrementTimesteps = this.Control.incrementTimesteps;
                    m_BDF_Timestepper.PushLevelSet = this.PushLevelSetAndRelatedStuff;
                    m_BDF_Timestepper.IterUnderrelax = this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? this.Control.LSunderrelax : 1.0;
                    m_BDF_Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = this.Control.LinearSolver, TestSolution = true };
                    //m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;

                } else {

                    throw new NotSupportedException();

                    //m_RK_Timestepper = new XdgRKTimestepping(
                    //    this.CurrentSolution.Mapping.Fields.ToArray(),
                    //    this.CurrentResidual.Mapping.Fields.ToArray(),
                    //    LsTrk,
                    //    DelComputeOperatorMatrix, DelUpdateLevelSet, DelUpdateCutCellMetrics,
                    //    rksch,
                    //    this.Control.Timestepper_LevelSetHandling,
                    //    mmsd,
                    //    (this.Control.PhysicalParameters.IncludeConvection) ? SpatialOperatorType.Nonlinear : SpatialOperatorType.LinearTimeDependent,
                    //    MassScale,
                    //    this.MultigridOperatorConfig, base.MultigridSequence,
                    //    this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold, 
                    //    true);
                    //m_RK_Timestepper.m_ResLogger = base.ResLogger;
                    //m_RK_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                }

            } else {

                Debug.Assert(object.ReferenceEquals(this.MultigridSequence[0].ParentGrid, this.GridData));

                //this.LsTrk.UpdateTracker();

                //FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                //CellMask Accepted = LsTrk.Regions.GetCutCellMask();
                //CellMask ActiveField = LsTrk.Regions.GetNearFieldMask(1);
                //CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                //FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                ContinuityEnforcer = new ContinuityProjection(DGLevelSet: this.DGLevSet.Current, gridData: GridData, Option: Control.LSContiProjectionMethod);
                //var Near = this.LsTrk.Regions.GetNearMask4LevSet(0, 1);
                //var PosFF = this.LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
                //ContinuityEnforcer.SetFarField(this.DGLevSet.Current, Near, PosFF);

                m_BDF_Timestepper.DataRestoreAfterBalancing(L,
                    ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure),
                    ArrayTools.Cat<DGField>(this.XDGvelocity.ResidualMomentum.ToArray(), this.ResidualContinuity),
                    this.LsTrk, this.MultigridSequence);

                //Console.WriteLine("number of cells {0}", this.Grid.NumberOfCells);
                //PlotCurrentState(hack_Phystime, new TimestepNumber(hack_TimestepIndex, 0), 2);


            }
            #endregion

        }


        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
        }



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

            this.XNSE_Operator.AssembleMatrix_Timestepper(
                this.m_HMForder,
                OpMtx, OpAffine,
                AgglomeratedCellLengthScales,
                CurrentState,
                SurfaceForce,
                filtLevSetGradient,
                this.Curvature,
                codMap,
                domMap,
                phystime);

            if (filtLevSetGradient != null) {
                if (this.Control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource == CurvatureAlgorithms.LevelSetSource.fromC0) {
                    this.LevSetGradient.Clear();
                    this.LevSetGradient.Acc(1.0, filtLevSetGradient);
                } else {
                    this.DGLevSetGradient.Clear();
                    this.DGLevSetGradient.Acc(1.0, filtLevSetGradient);
                }
            }
                       


            // ============================
            // something with surface tension ?????
            // ============================

            {
                if (this.Control.PhysicalParameters.useArtificialSurfaceForce == true)
                    throw new NotSupportedException("Not supported for this hack.");
                if (this.Control.UseXDG4Velocity != true)
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
            //  'FakePoisson'
            // ============================

            //if (base.Control.FakePoisson)
            //    OpMtx.Acc(1.0, WholeMassMatrix);

            // ============================
            //  Add Gravity
            // ============================
            // Dimension: [ rho * G ] = mass / time^2 / len^2 == [ d/dt rho U ]
            var WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity != null ? this.XDGvelocity.Gravity.ToArray<DGField>() : this.DGvelocity.Gravity.ToArray<DGField>(), new XDGField(this.Pressure.Basis)));
            WholeMassMatrix.SpMV(1.0, WholeGravity, 1.0, OpAffine);

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

        int hack_TimestepIndex;
        double hack_Phystime;


        /// <summary>
        /// 
        /// </summary>
        protected override ITimestepInfo SaveToDatabase(TimestepNumber timestepno, double t) {
            var tsi = base.SaveToDatabase(timestepno, t);

            if(tsi != null && m_BDF_Timestepper != null) {
                int S = m_BDF_Timestepper.GetNumberOfStages;

                SinglePhaseField LsBkUp = new SinglePhaseField(this.LevSet.Basis);
                LsBkUp.Acc(1.0, this.LevSet);

                ICollection<DGField>[] restartFields = m_BDF_Timestepper.GetRestartInfos();

                if(S > 1 && this.Control.saveperiod >= S && restartFields != null) {

                    // save additional timesteps/information for restart
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    for(int ti = 1; ti < S; ti++) {

                        //SinglePhaseField LsBkUp = new SinglePhaseField(this.LevSet.Basis);
                        //LsBkUp.Acc(1.0, this.LevSet);

                        ICollection<DGField> restartIOFields = new List<DGField>();
                        foreach(DGField f in this.IOFields) {

                            int rfidx = restartFields[ti - 1].IndexWhere(rf => rf.Identification == f.Identification);
                            if(rfidx > -1) {
                                DGField rf = restartFields[ti - 1].ElementAt(rfidx);
                                if(f.Identification == "Phi") {
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
                            rtsi = this.DatabaseDriver.SaveTimestep(
                                t - ti * this.Control.GetFixedTimestep(),
                                tsn,
                                this.CurrentSessionInfo,
                                this.GridData,
                                restartIOFields);
                        } catch(Exception ee) {
                            Console.Error.WriteLine(ee.GetType().Name + " on rank " + this.MPIRank + " saving time-step " + tsn + ": " + ee.Message);
                            Console.Error.WriteLine(ee.StackTrace);
                            //tsi = null;
                            //e = ee;

                            if(ContinueOnIOError) {
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

#if DEBUG
            //Debug/Test code for XDG database interaction

            if(tsi != null) {
                // checking some neccessary reference-equalities BEFORE serialisation
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                LevelSet.LevelSetInitializer lsi_1 = (LevelSet.LevelSetInitializer)(tsi.FieldInitializers.Single(fi => fi.Identification == this.LevSet.Identification));
                XDGField.FieldInitializer pri = (XDGField.FieldInitializer)(tsi.FieldInitializers.Single(fi => fi.Identification == this.Pressure.Identification));

                LevelSetTracker.LevelSetTrackerInitializer trki = ((XDGBasis.XDGBasisInitializer)(pri.BasisInfo)).TrackerInitializer;

                LevelSet.LevelSetInitializer lsi_2 = trki.LevelSets[0];

                Debug.Assert(object.ReferenceEquals(lsi_1, lsi_2));

                foreach(XDGField.FieldInitializer fi in tsi.FieldInitializers.Where(fii => fii is XDGField.XDGFieldInitializer)) {
                    LevelSetTracker.LevelSetTrackerInitializer trki_alt = ((XDGBasis.XDGBasisInitializer)(fi.BasisInfo)).TrackerInitializer;
                    Debug.Assert(object.ReferenceEquals(trki, trki_alt));
                }
            }


            if(tsi != null) {
                // checking some neccessary equalities AFTER serialisation
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                

                var tsi_alt = this.DatabaseDriver.LoadTimestepInfo(tsi.ID, base.CurrentSessionInfo, base.GetDatabase());

                Debug.Assert(!object.ReferenceEquals(tsi, tsi_alt));


                LevelSet.LevelSetInitializer lsi_1 = (LevelSet.LevelSetInitializer)(tsi_alt.FieldInitializers.Single(fi => fi.Identification == "Phi"));
                XDGField.FieldInitializer pri = (XDGField.FieldInitializer)(tsi_alt.FieldInitializers.Single(fi => fi.Identification == this.Pressure.Identification));

                LevelSetTracker.LevelSetTrackerInitializer trki = ((XDGBasis.XDGBasisInitializer)(pri.BasisInfo)).TrackerInitializer;

                LevelSet.LevelSetInitializer lsi_2 = trki.LevelSets[0];

                Debug.Assert(lsi_1.Equals(lsi_2));

                foreach(XDGField.FieldInitializer fi in tsi_alt.FieldInitializers.Where(fii => fii is XDGField.XDGFieldInitializer)) {
                    LevelSetTracker.LevelSetTrackerInitializer trki_alt = ((XDGBasis.XDGBasisInitializer)(fi.BasisInfo)).TrackerInitializer;
                    Debug.Assert(trki.Equals(trki_alt));
                }


                var Fields = this.DatabaseDriver.LoadFields(tsi_alt, this.GridData);
                LevelSet Rphi_1 = (LevelSet)(Fields.Single(f => f.Identification.Equals(this.LevSet.Identification)));

                XDGField Rpressure = (XDGField)(Fields.Single(f => f.Identification.Equals(this.Pressure.Identification)));

                LevelSetTracker Rtracker = Rpressure.Basis.Tracker;
                Debug.Assert(!object.ReferenceEquals(this.LsTrk, Rtracker));
                Debug.Assert(object.ReferenceEquals(Rtracker.LevelSets[0], Rphi_1));

                foreach(XDGField xf in Fields.Where(fii => fii is XDGField)) {
                    Debug.Assert(object.ReferenceEquals(xf.Basis.Tracker, Rtracker));
                }
            }

#endif

            return tsi;
        }



        /// <summary>
        /// Depending on settings <see cref="AppControl.CompMode"/>, computes either one timestep or a steady-state solution.
        /// </summary>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {

            TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
            int D = this.GridData.SpatialDimension;
            base.ResLogger.TimeStep = TimestepInt;
            hack_TimestepIndex = TimestepInt;
            hack_Phystime = phystime;


            if (Control.SkipSolveAndEvaluateResidual) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++
                // setup: project exact solution -- for consistency tests
                // +++++++++++++++++++++++++++++++++++++++++++++++++

                foreach (string spc in LsTrk.SpeciesNames) {
                    for (int d = 0; d < this.GridData.SpatialDimension; d++) {
                        ConventionalDGField Vel_d;
                        if (this.CurrentVel[d] is XDGField)
                            Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);
                        else
                            Vel_d = (ConventionalDGField)this.CurrentVel[d];
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
                if (this.XDGvelocity != null) {
                    for (int d = 0; d < D; d++) {
                        // Gravity must be set up like this to avoid regions of zero gravity when updating the level-set
                        this.XDGvelocity.Gravity[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
                    }
                } else {
                    throw new NotSupportedException();
                }

                /// +++++++++++++++++++++++++++++++++++++
                /// compute/check time step restrictions
                /// +++++++++++++++++++++++++++++++++++++

                dt = base.Control.dtFixed;

                // Level-Set motion-CFL
                double LevSet_Deg2 = this.DGLevSet.Current.Basis.Degree;
                LevSet_Deg2 = LevSet_Deg2 * LevSet_Deg2;
                double dt_LevSetCFL = base.GridData.ComputeCFLTime(this.ExtensionVelocity.Current, dt * LevSet_Deg2);
                dt_LevSetCFL = dt_LevSetCFL / LevSet_Deg2;

                //dt = Math.Min(dt, dt_LevSetCFL);

                // Capillary Timestep restriction
                if (this.Control.PhysicalParameters.Sigma != 0.0) {
                    MultidimensionalArray h_mins = this.GridData.Cells.h_min;
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

            if (m_BDF_Timestepper != null) {
                m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            } else {
                //m_RK_Timestepper.Solve(phystime, dt);
            }



            Postprocessing(TimestepInt, phystime, dt, TimestepNo);
            


            // ================
            // Good bye
            // ================
            if(this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
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



        /// <summary>
        /// 
        /// </summary>
        public void PushLevelSetAndRelatedStuff() {

            if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.Fourier) {
                Fourier_Timestepper.updateFourierLevSet();
            }

            this.ExtensionVelocity.IncreaseHistoryLength(1);
            this.ExtensionVelocity.Push();

            this.DGLevSet.IncreaseHistoryLength(1);
            this.DGLevSet.Push();
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="TimestepInt"></param>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <param name="TimestepNo"></param>
        private void Postprocessing(int TimestepInt, double phystime, double dt, TimestepNumber TimestepNo) {


            // ======================================
            // Check jump conditions at the interface 
            // ======================================

            #region check jump conditions

            if (this.Control.CheckJumpConditions) {

                // mass balance
                double velJump_Norm = XNSEUtils.VelocityJumpNorm(this.XDGvelocity.Velocity, false, this.m_HMForder);

                Console.WriteLine("velocity jump norm = {0}", velJump_Norm);

                this.MassBalanceAtInterface.Clear();
                XNSEUtils.ProjectMassBalanceNorm(this.MassBalanceAtInterface, 1.0, this.XDGvelocity.Velocity, this.m_HMForder);


                // momentum balance
                CurvatureAlgorithms.CurvatureDriver(
                    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                    CurvatureAlgorithms.FilterConfiguration.NoFilter,
                    this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                    this.m_HMForder, this.DGLevSet.Current);

                ConventionalDGField[] meanVelocity = XNSEUtils.GetMeanVelocity(this.XDGvelocity.Velocity, this.LsTrk,
                    this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B);


                double[] momBal_Norm = XNSEUtils.MomentumBalanceNormAtInterface(this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                    this.Control.PhysicalParameters, this.Control.AdvancedDiscretizationOptions.SurfStressTensor, this.m_HMForder);

                Console.WriteLine("x-momentum balance norm = {0}", momBal_Norm[0]);
                Console.WriteLine("y-momentum balance norm = {0}", momBal_Norm[1]);

                for (int d = 0; d < this.Grid.SpatialDimension; d++) {
                    this.MomentumBalanceAtInterface[d].Clear();
                    XNSEUtils.ProjectMomentumBalanceNorm(this.MomentumBalanceAtInterface[d], 1.0, this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                        this.Control.PhysicalParameters, this.Control.AdvancedDiscretizationOptions.SurfStressTensor, d, this.m_HMForder);
                }


                // energy balance
                double energyBal_Norm = XNSEUtils.EnergyBalanceNormAtInterface(this.Pressure, this.XDGvelocity.Velocity, meanVelocity, this.Curvature,
                    this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, this.Control.PhysicalParameters.Sigma, this.m_HMForder);

                Console.WriteLine("energy balance norm = {0}", energyBal_Norm);

                this.EnergyBalanceAtInterface.Clear();
                XNSEUtils.ProjectEnergyBalanceNorm(this.EnergyBalanceAtInterface, 1.0, this.Pressure, this.XDGvelocity.Velocity, meanVelocity, this.Curvature,
                    this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, this.Control.PhysicalParameters.Sigma, this.m_HMForder);

            }

            #endregion


            // ====================================================
            // Compute energy of system
            // ====================================================

            #region energy computation

            if (this.Control.ComputeEnergy) {

                // compute current energies
                double[] RhoS = new double[] { this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B };
                double currentKinEnergy = XNSEUtils.GetKineticEnergy(this.LsTrk, this.XDGvelocity.Velocity.ToArray(), RhoS, this.m_HMForder);
                double currentSurfEnergy = XNSEUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma, this.m_HMForder);

                // compute changerates (kinetic, surface)
                double CR_KinEnergy = 0.0;
                double CR_SurfEnergy = 0.0;
                if (this.Control.CompMode == AppControl._CompMode.Transient) {
                    double prevKinEnergy = XNSEUtils.GetKineticEnergy(this.LsTrk, this.prevVel, RhoS, this.m_HMForder, 0);
                    CR_KinEnergy = (currentKinEnergy - prevKinEnergy) / dt;

                    double prevSurfEnergy = XNSEUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma, this.m_HMForder, 0);
                    CR_SurfEnergy = (currentSurfEnergy - prevSurfEnergy) / dt;

                    Console.WriteLine("current kinetic energy = {0}; actual changerate = {1}", currentKinEnergy, CR_KinEnergy);
                    Console.WriteLine("current surface energy = {0}; actual changerate = {1}", currentSurfEnergy, CR_SurfEnergy);
                }


                ConventionalDGField[] meanVelocity = XNSEUtils.GetMeanVelocity(this.XDGvelocity.Velocity, this.LsTrk,
                        this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B);
                double SurfDivergence = XNSEUtils.GetSurfaceChangerate(this.LsTrk, meanVelocity, this.m_HMForder);


                // changerate of the discretization
                double[] muS = new double[] { this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B };
                //double prevCR_discEnergy = XNSEUtils.GetEnergyChangerate(this.LsTrk, prevVel, muS, this.m_HMForder, 0);
                double currentCR_discEnergy = XNSEUtils.GetEnergyChangerate(this.LsTrk, this.XDGvelocity.Velocity.ToArray(), muS, this.m_HMForder);


                this.EnergyLogger.TimeStep = TimestepInt;
                this.EnergyLogger.CustomValue(phystime + dt, "PhysicalTime");
                this.EnergyLogger.CustomValue(currentKinEnergy, "KineticEnergy");
                this.EnergyLogger.CustomValue(currentSurfEnergy, "SurfaceEnergy");
                this.EnergyLogger.CustomValue(CR_KinEnergy, "ChangerateKineticEnergy");
                this.EnergyLogger.CustomValue(CR_SurfEnergy, "ChangerateSurfaceEnergy");
                this.EnergyLogger.CustomValue(SurfDivergence, "SurfaceDivergence");
                this.EnergyLogger.CustomValue(currentCR_discEnergy, "BulkDissipationrate");


                // surface viscosity parts
                if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor != SurfaceSressTensor.Isotropic) {

                    double shearViscEnergyCR = 0.0;
                    double dilViscEnergyCR = 0.0;

                    // surface shear viscosity energy
                    if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation
                        || this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        shearViscEnergyCR = XNSEUtils.GetInterfaceShearViscosityEnergyCR(this.LsTrk, meanVelocity, this.Control.PhysicalParameters.mu_I, this.m_HMForder);
                    }

                    // surface dilatational viscosity energy
                    if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation
                        || this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        dilViscEnergyCR = XNSEUtils.GetInterfaceDilatationalViscosityEnergyCR(this.LsTrk, meanVelocity, this.Control.PhysicalParameters.lambda_I, this.m_HMForder);
                    }


                    this.EnergyLogger.CustomValue(shearViscEnergyCR, "ShearViscosityDR");
                    this.EnergyLogger.CustomValue(dilViscEnergyCR, "DilatationalViscosityDR");

                    Console.WriteLine("current discretization Energy changerate = {0} / {1}", currentCR_discEnergy - shearViscEnergyCR - dilViscEnergyCR, CR_KinEnergy + CR_SurfEnergy);

                } else {

                    Console.WriteLine("current discretization Energy changerate = {0} / {1}", currentCR_discEnergy, CR_KinEnergy + CR_SurfEnergy);
                }


                // logging
                // =======

                this.EnergyLogger.NextTimestep(true);



                //double[] RhoS = new double[] { this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B };
                //double newKinEnergy = XNSEUtils.GetKineticEnergy(this.LsTrk, this.CurrentVel, RhoS, this.m_HMForder);
                //double oldKinEnergy;
                //if (base.Control.CompMode == AppControl._CompMode.Transient) {
                //    DGField[] prevVel;
                //    if (this.XDGvelocity != null)
                //        prevVel = this.XDGvelocity.Velocity.ToArray();//<DGField>();
                //    else
                //        prevVel = this.DGvelocity.Velocity.ToArray();
                //    oldKinEnergy = XNSEUtils.GetKineticEnergy(this.LsTrk, prevVel, RhoS, this.m_HMForder);
                //} else if (base.Control.CompMode == AppControl._CompMode.Steady) {
                //    oldKinEnergy = newKinEnergy;
                //} else {
                //    throw new NotSupportedException();
                //}
                //double surfEnergy = XNSEUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma.Abs(), this.m_HMForder);


                //// Logging and Console Output
                //// ===========================

                //this.EnergyLogger.TimeStep = TimestepInt;
                //this.EnergyLogger.CustomValue(phystime + dt, "PhysicalTime");
                //this.EnergyLogger.CustomValue(oldKinEnergy, "OldKineticEnergy");
                //this.EnergyLogger.CustomValue(newKinEnergy, "NewKineticEnergy");
                //this.EnergyLogger.CustomValue(surfEnergy, "SurfaceEnergy");

                //this.EnergyLogger.NextTimestep(true);

            }

            #endregion


            // ====================================
            // divergence of velocity
            // ====================================

            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            //this.divVelocity.Clear();
            //this.divVelocity.Divergence(1.0, this.XDGvelocity.Velocity);
            //ilPSP.Environment.StdoutOnlyOnRank0 = true;


            // ====================================
            // L2 error against exact solution
            // ====================================

            this.ComputeL2Error(phystime + dt);

            // =========== 
            // check area
            // ===========

            //SpeciesId spcId = LsTrk.SpeciesIdS[0];
            //double area = XNSEUtils.GetSpeciesArea(LsTrk, spcId, MomentFittingVariant);
            //Console.WriteLine("Area of species 'A' = {0}", area);


            // ====================================
            // IO related to Fourier level set
            // ====================================

            if (base.MPIRank == 0) {
                // save restart infos for FLS
                if (Log_FourierLS != null) {
                    Guid vecSamplP_id = this.DatabaseDriver.SaveVector<double>(Fourier_LevSet.getRestartInfo());
                    Log_FourierLS.WriteLine(vecSamplP_id);
                    Log_FourierLS.Flush();
                }
                // Log_files for FLS
                //if (this.Control.FourierLevSetControl.WriteFLSdata) {
                //    Fourier_LevSet.saveToLogFiles(TimestepNo.MajorNumber, phystime + dt);
                //}
            }


            // ====================================================================== 
            // IO for further external postprocessing/ Query handling for Testprogram
            // ======================================================================

            if (this.Control.TestMode == true) {
                LogQueryValue(phystime + dt);
            } else {
                if (Log != null && this.Control.LogValues != XNSE_Control.LoggingValues.none && base.MPIRank == 0 && (TimestepNo.MajorNumber % this.Control.LogPeriod == 0))
                    WriteLogLine(TimestepNo, phystime + dt);
            }

            //Console.WriteLine("Pause");

            //=======================
            //var jmpNorm = XNSEUtils.VelocityJumpNorm(this.XDGvelocity.Velocity, true, MomentFittingVariant, -1);
            //Console.WriteLine("Velocity Jump Norm: " + jmpNorm);
            //var jmpStressNorm = XNSEUtils.MomentumJumpNorm(this.XDGvelocity.Velocity, this.Pressure, this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, MomentFittingVariant, -1);
            //Console.WriteLine("Stress Jump Norm [0]: " + jmpStressNorm[0]);
            //Console.WriteLine("Stress Jump Norm [1]: " + jmpStressNorm[1]);
            //PrintVelocityAtLevSet(TimestepNo.MajorNumber);

        }


        /// <summary>
        /// testcase specific LogFile
        /// </summary>
        TextWriter Log;

        /// <summary>
        /// saves interface points
        /// </summary>
        TextWriter LogInterfaceP;


        /// <summary>
        /// initializes the format of the Log File
        /// </summary>
        /// <param name="sessionID"></param>
        public void InitLogFile(Guid sessionID) {

            if (this.Control.WriteInterfaceP) {
                LogInterfaceP = base.DatabaseDriver.FsDriver.GetNewLog("InterfaceP", sessionID);
                string header = String.Format("{0}\t{1}\t{2}", "#timestep", "#time", "interfacePoints");
                LogInterfaceP.WriteLine(header);
                LogInterfaceP.Flush(); 
            }

            switch (this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        // File for physical data
                        TextWriter setUpData = base.DatabaseDriver.FsDriver.GetNewLog("SetUpData", sessionID);
                        string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", "lambda", "H0", "rho1", "rho2", "mu1", "mu2", "sigma", "g");
                        setUpData.WriteLine(header);
                        setUpData.Flush();
                        string data = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", this.Control.AdditionalParameters[1], this.Control.AdditionalParameters[2], this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B,
                            this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, this.Control.PhysicalParameters.Sigma, this.Control.AdditionalParameters[3]);
                        setUpData.WriteLine(data);
                        setUpData.Flush();


                        // Log file for the interface height
                        Log = base.DatabaseDriver.FsDriver.GetNewLog("Amplitude", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#timestep", "#time", "magnitude", "real", "imaginary");
                        Log.WriteLine(header);
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.RisingBubble: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("BenchmarkQuantities_RisingBubble", sessionID);
                        string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", "#timestep", "#time", "center of mass - x", "center of mass - y", "circularity", "rise velocity");
                        Log.WriteLine(header);
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.MovingContactLine: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("ContactAngle", sessionID);
                        string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "#time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-angle");
                        Log.WriteLine(header);
                        Log.Flush();

                        return;
                    }
                default:
                    throw new ArgumentException("No specified LogFormat");
            }
        }

        /// <summary>
        /// writes one line to the Log File
        /// </summary>
        public void WriteLogLine(TimestepNumber TimestepNo, double phystime) {

            if (this.Control.WriteInterfaceP) {
                double[] interfaceP;
                if (Fourier_LevSet != null) {
                    interfaceP = Fourier_LevSet.current_interfaceP.To1DArray();
                } else {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                    interfaceP = interP.ResizeShallow(interP.Length).To1DArray();
                }
                string logline = String.Format("{0}\t{1}", TimestepNo, phystime);
                //for (int ip = 0; ip < interfaceP.Length; ip++) {
                //    logline = logline + "\t" + interfaceP[ip].ToString();
                //}
                logline = logline + "\t" + String.Join("\t", interfaceP.Select(ip => ip.ToString()).ToArray());
                LogInterfaceP.WriteLine(logline);
                LogInterfaceP.Flush();
            }

            switch (this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        Complex DFT_k;
                        int numP;
                        if (Fourier_LevSet != null) {
                            //amplitude = 2.0 * (Fourier_LevSet.DFT_coeff[1].Magnitude / Fourier_LevSet.current_samplP.Length);
                            DFT_k = Fourier_LevSet.DFT_coeff[(int)this.Control.AdditionalParameters[0]];
                            numP = Fourier_LevSet.current_samplP.Length;
                        } else {
                            MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                            Complex[] DFT_coeff = DiscreteFourierTransformation.TransformForward_nonequidistant(interP, this.Control.AdditionalParameters[1]);
                            DFT_k = DFT_coeff[(int)this.Control.AdditionalParameters[0]];
                            numP = interP.Lengths[0];
                            //amplitude = -2.0 * DFT_coeff[1].Imaginary / (double)interP.Lengths[0];
                            //amplitude = DiscreteFourierTransformation.SingleSidedPowerSpectrum(DFT_coeff)[1];
                        }
                        string logline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, 2.0 * DFT_k.Magnitude / numP, 2.0 * DFT_k.Real / numP, -2.0 * DFT_k.Imaginary / numP);
                        Log.WriteLine(logline);
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.RisingBubble: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", TimestepNo, phystime, BmQ_RB[0], BmQ_RB[1], BmQ_RB[2], BmQ_RB[4]);
                        Log.WriteLine(line);
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.MovingContactLine: {

                        // contact angles at contact points
                        //=================================

                        List<double[]> contactPoints = new List<double[]>();
                        List<double[]> contactVelocities = new List<double[]>();
                        List<double> contactAngles = new List<double>();

                        ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(this.CurrentVel);

                        var Phi = (LevelSet)LsTrk.LevelSets[0];
                        var LevelSetGradient = new VectorField<SinglePhaseField>(Grid.SpatialDimension, Phi.Basis, SinglePhaseField.Factory);
                        LevelSetGradient.Gradient(1.0, (SinglePhaseField)LsTrk.LevelSets[0]);
                        SinglePhaseField[] Normals = LevelSetGradient.ToArray();

                        XQuadSchemeHelper SchemeHelper = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadSchemeHelper;
                        EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"));

                        var QuadDom = SurfaceElement_Edge.Domain;
                        var boundaryCutEdge = QuadDom.Intersect(this.GridData.GetBoundaryEdgeMask());

                        var factory = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, LsTrk.GridDat.Grid.RefElements[0]);
                        SurfaceElement_Edge = new EdgeQuadratureScheme(factory, boundaryCutEdge);

                        EdgeQuadrature.GetQuadrature(new int[] { 5 }, LsTrk.GridDat,
                            SurfaceElement_Edge.Compile(LsTrk.GridDat, 0),
                            delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {

                                // contact point
                                NodeSet Enode_l = QR.Nodes;
                                int trf = LsTrk.GridDat.Edges.Edge2CellTrafoIndex[i0, 0];
                                NodeSet Vnode_l = Enode_l.GetVolumeNodeSet(LsTrk.GridDat, trf);
                                NodeSet Vnode_g = Vnode_l.CloneAs();
                                int cell = LsTrk.GridDat.Edges.CellIndices[i0, 0];
                                LsTrk.GridDat.TransformLocal2Global(Vnode_l, Vnode_g, cell);
                                //Console.WriteLine("contact point: ({0},{1})", Vnode_g[0, 0], Vnode_g[0, 1]);

                                int D = Grid.SpatialDimension;
                                for (int d = 0; d < D; d++) {
                                    EvalResult[0, 0, d] = Vnode_g[0, d];
                                }

                                // contact line velocity
                                MultidimensionalArray U_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                MultidimensionalArray U_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                for (int d = 0; d < D; d++) {
                                    (meanVelocity[d] as SinglePhaseField).EvaluateEdge(i0, length, QR.Nodes, U_IN.ExtractSubArrayShallow(-1, -1, d), U_OUT.ExtractSubArrayShallow(-1, -1, d));
                                }

                                for (int d = 0; d < D; d++) {
                                    EvalResult[0, 0, 2 + d] = U_IN[0, 0, d];
                                }

                                // contact angle
                                MultidimensionalArray normal_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                MultidimensionalArray normal_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                for (int d = 0; d < D; d++) {
                                    Normals[d].EvaluateEdge(i0, length, QR.Nodes, normal_IN.ExtractSubArrayShallow(-1, -1, d), normal_OUT.ExtractSubArrayShallow(-1, -1, d));
                                }

                                double theta_surf = Math.Atan2(normal_IN[0, 0, 1], normal_IN[0, 0, 0]);
                                double theta_edge = Math.Atan2(LsTrk.GridDat.Edges.NormalsForAffine[i0, 1], LsTrk.GridDat.Edges.NormalsForAffine[i0, 0]);
                                double theta = (theta_surf - theta_edge) * (180 / Math.PI);

                                EvalResult[0, 0, 2 * D] = (theta > 180) ? theta - 180 : theta;
                                //Console.WriteLine("contact angle = {0}", EvalResult[0, 0, 2]);

                            },
                            delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                                int D = Grid.SpatialDimension;
                                for (int i = 0; i < length; i++) {
                                    if (ResultsOfIntegration[i, 2*D] != 0.0) {
                                        contactAngles.Add(Math.Abs(ResultsOfIntegration[i, 2*D]));
                                        double[] cp = new double[D];
                                        double[] cpV = new double[D];
                                        for (int d = 0; d < D; d++) {
                                            cp[d] = ResultsOfIntegration[i, d];
                                            cpV[d] = ResultsOfIntegration[i, 2 + d];
                                        }
                                        contactPoints.Add(cp);
                                        contactVelocities.Add(cpV);
                                    }
                                }
                            }
                        ).Execute();


                        for (int p = 0; p < contactAngles.Count; p++) {
                            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime, contactPoints.ElementAt(p)[0], contactPoints.ElementAt(p)[1], contactVelocities.ElementAt(p)[0], contactVelocities.ElementAt(p)[1], contactAngles.ElementAt(p));
                            Log.WriteLine(line);
                        }
                        Log.Flush();

                        return;
                    }
                default:
                    throw new ArgumentException("No specified LogFormat");
            }

        }


        /// <summary>
        /// encapsulated handling of query values
        /// </summary>
        public void LogQueryValue(double phystime) {

            base.QueryResultTable.LogValue("time", phystime);

            if (this.Control.WriteInterfaceP) {
                double[] interfaceP;
                if (Fourier_LevSet != null) {
                    interfaceP = Fourier_LevSet.current_interfaceP.To1DArray();
                } else {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                    interfaceP = interP.ResizeShallow(interP.Length).To1DArray();
                }

                base.QueryResultTable.LogValue("interfaceP", interfaceP);
                
            }

            switch (this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        double amplitude;
                        if (Fourier_LevSet != null) {
                            amplitude = 2.0 * (Fourier_LevSet.DFT_coeff[1].Magnitude / Fourier_LevSet.current_samplP.Length);
                        } else {
                            MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                            Complex[] DFT_coeff = DiscreteFourierTransformation.TransformForward_nonequidistant(interP, this.Control.AdditionalParameters[1]);
                            amplitude = -2.0 * DFT_coeff[1].Imaginary / (double)interP.Lengths[0];
                            //amplitude = DiscreteFourierTransformation.SingleSidedPowerSpectrum(DFT_coeff)[1];
                        }

                        base.QueryResultTable.LogValue("amplitude", amplitude);
                        return;
                    }
                case XNSE_Control.LoggingValues.RisingBubble: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        base.QueryResultTable.LogValue("yCM", BmQ_RB[1]);
                        base.QueryResultTable.LogValue("circ", BmQ_RB[2]);
                        base.QueryResultTable.LogValue("riseV", BmQ_RB[4]);
                        
                        return;
                    }
                case XNSE_Control.LoggingValues.LinelikeLS: {
                        break;
                    }
                case XNSE_Control.LoggingValues.CirclelikeLS: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        base.QueryResultTable.LogValue("xM", BmQ_RB[0]);
                        base.QueryResultTable.LogValue("yM", BmQ_RB[1]);
                        base.QueryResultTable.LogValue("circ", BmQ_RB[2]);
                        base.QueryResultTable.LogValue("vM_x", BmQ_RB[3]);
                        base.QueryResultTable.LogValue("vM_y", BmQ_RB[4]);

                        break;
                    }
                default:
                    return;
            }

        }



        public double[] ComputeBenchmarkQuantities_RisingBubble() {

            int order = 0;
            if(LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;
            
            // area of bubble
            double area = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        area += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            // center of mass/geometric center (for incommpressible fluid)
            int D = this.Grid.SpatialDimension;
            MultidimensionalArray center = MultidimensionalArray.Create(1, D);
            CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    NodeSet nodes_global = QR.Nodes.CloneAs();
                    for (int i = i0; i < i0 + Length; i++) {
                        LsTrk.GridDat.TransformLocal2Global(QR.Nodes, nodes_global, i);
                        EvalResult.AccSubArray(1.0, nodes_global, new int[] { i - i0, -1, -1 });
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            center[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();

            center.Scale(1.0 / area);

            // rise velocity
            MultidimensionalArray VelocityAtCenter = MultidimensionalArray.Create(1, D);

            // integral computation
            CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    for (int d = 0; d < D; d++) {
                        this.CurrentVel[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            VelocityAtCenter[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();
            VelocityAtCenter.Scale(1.0 / area);

            double v_rise = VelocityAtCenter[0, 1];

            //Console.WriteLine("rise velocity = " + v_rise);


            // circularity
            double diamtr_c = Math.Sqrt(4 * area / Math.PI);
            double perimtr_b = 0.0;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        perimtr_b += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            double circ = Math.PI * diamtr_c / perimtr_b;

            return new double[] { center[0, 0], center[0, 1], circ, VelocityAtCenter[0, 0], VelocityAtCenter[0, 1] };
        }


        UnsetteledCoordinateMapping SaddlePointProblemMapping {
            get {
                return this.CurrentSolution.Mapping;
            }
        }

        
       
        int RepairZeroRows(MsrMatrix Mtx) {
            int NoOfZeroRows = 0;
            for (int iRow = Mtx.RowPartitioning.i0; iRow < Mtx.RowPartitioning.iE; iRow++) {
                if (Mtx.GetNoOfNonZerosPerRow(iRow) == 0) {
                    Mtx[iRow, iRow] = +1.0;
                    NoOfZeroRows++;
                }
            }
            return NoOfZeroRows;
        }

        /// <summary>
        /// Computes condition number, etc. of the current system matrix.
        /// </summary>
        /// <param name="CheckAssertions"></param>
        /// <param name="AnalysisLevel">
        /// - equal 0: check that pressure gradient and velocity divergence are transpose
        /// - equal 1: in addition, positive definiteness test.
        /// - equal 2: in addition, check condition number and eigenvalues using MATLAB
        /// </param>
        public void SpatialOperatorMatrixAnalysis(bool CheckAssertions, int AnalysisLevel) {
            using (new FuncTrace()) {
                int D = this.Grid.SpatialDimension;

                if (AnalysisLevel < 0 || AnalysisLevel > 2)
                    throw new ArgumentException();

                // ========================================================
                // compute agglomeration & operator (saddle-point) matrix
                // ========================================================

                //var CurrentAgglomeration = new MultiphaseCellAgglomerator(
                //        this.DelUpdateCutCellMetrics(),
                //        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                //        AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);
                var CurrentAgglomeration = this.LsTrk.GetAgglomerator(
                    this.LsTrk.SpeciesIdS.ToArray(), m_HMForder, 
                    this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold, 
                    AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);

                BlockMsrMatrix SaddlePointMatrix = new BlockMsrMatrix(this.SaddlePointProblemMapping);
                double[] AffineDummy = new double[this.SaddlePointProblemMapping.LocalLength];

                

                DelComputeOperatorMatrix(SaddlePointMatrix, AffineDummy, this.SaddlePointProblemMapping,
                    this.CurrentSolution.Mapping.Fields.ToArray(), CurrentAgglomeration.CellLengthScales, 0.0);

                // =============================
                // AnalysisLevel 0
                // =============================
                {
                    var SaddlePointMatrixT = SaddlePointMatrix.Transpose();

                    CoordinateVector TestVec = new CoordinateVector(this.CurrentSolution.Mapping.Fields.Select(f => f.CloneAs()).ToArray());

                    double testsumPos = 0.0;
                    double testsumNeg = 0.0;
                    for (int rnd_seed = 0; rnd_seed < 20; rnd_seed++) {

                        // fill the pressure components of the test vector
                        TestVec.Clear();
                        Random rnd = new Random(rnd_seed);
                        XDGField Pressack = TestVec.Mapping.Fields[D] as XDGField;
                        int J = this.GridData.Cells.NoOfLocalUpdatedCells;
                        for (int j = 0; j < J; j++) {
                            int N = Pressack.Basis.GetLength(j);

                            for (int n = 0; n < N; n++)
                                Pressack.Coordinates[j, n] = rnd.NextDouble();
                        }

                        //Pressack.Clear(H.Complement());
                        //Pressack.GetSpeciesShadowField("A").Clear();

                        // Gradient times P:
                        double[] R1 = new double[TestVec.Count];
                        SaddlePointMatrix.SpMV(1.0, TestVec, 0.0, R1);       // R1 = Grad * P
                                                                     //Console.WriteLine("L2 of 'Grad * P': " + R1.L2Norm());

                        // transpose of Divergence times P: 
                        double[] R2 = new double[TestVec.Count];
                        SaddlePointMatrixT.SpMV(1.0, TestVec, 0.0, R2);      // R2 = divT * P
                                                                     //Console.WriteLine("L2 of 'divT * P': " + R2.L2Norm());


                        TestVec.Clear();
                        TestVec.Acc(1.0, R1);
                        TestVec.Acc(1.0, R2);


                        // analyze!
                        testsumNeg += GenericBlas.L2Dist(R1, R2);

                        R2.ScaleV(-1.0);
                        testsumPos += GenericBlas.L2Dist(R1, R2);

                    }

                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (+): " + testsumPos);
                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (-): " + testsumNeg);

                    if(CheckAssertions)
                        Assert.LessOrEqual(Math.Abs(testsumNeg), testsumPos*1.0e-13);
                }


                // =============================
                // AnalysisLevel 1 and 2
                // =============================

                if (AnalysisLevel > 0) {
                    AggregationGridBasis[][] MgBasis = AggregationGridBasis.CreateSequence(this.MultigridSequence, this.CurrentSolution.Mapping.BasisS);
                    //todo: AsyncCallback update
                    MgBasis.UpdateXdgAggregationBasis(CurrentAgglomeration);
                    MultigridOperator mgOp = new MultigridOperator(MgBasis, this.SaddlePointProblemMapping,
                        SaddlePointMatrix, this.MassFact.GetMassMatrix(this.SaddlePointProblemMapping, false),
                        this.MultigridOperatorConfig);

                    // extract
                    ////////////

                    MsrMatrix FullMatrix = mgOp.OperatorMatrix.ToMsrMatrix();

                    MsrMatrix DiffMatrix;
                    {
                        int[] VelVarIdx = D.ForLoop(d => d);

                        int[] USubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int[] USubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int L = USubMatrixIdx_Row.Length;

                        DiffMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(int[]), USubMatrixIdx_Col, default(int[]));
                    }

                    int Zeros_FullMatrix = RepairZeroRows(FullMatrix);
                    int Zeros_DiffMatrix = RepairZeroRows(DiffMatrix);

                    Console.WriteLine("Indefinite Basis elements in Diffusion matrix:\t" + Zeros_DiffMatrix);
                    Console.WriteLine("Indefinite Basis elements in Saddle-point matrix:\t" + Zeros_FullMatrix);

                    base.QueryHandler.ValueQuery("NoOfIndef_DiffMtx", Zeros_DiffMatrix, false);
                    base.QueryHandler.ValueQuery("NoOfIndef_FullMtx", Zeros_FullMatrix, false);

                    // operator analysis
                    //////////////////////

                    bool posDef;
                    if (AnalysisLevel > 1) {
                        // +++++++++++++++++++++++++++++++
                        // check condition number, etc
                        // +++++++++++++++++++++++++++++++

                        //DiffMatrix.SaveToTextFileSparse("C:\\tmp\\DiffMatrix.txt");

                        MultidimensionalArray ret = MultidimensionalArray.Create(1, 5);
                        Console.WriteLine("Calling MATLAB/Octave...");
                        using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                            bmc.PutSparseMatrix(FullMatrix, "FullMatrix");
                            bmc.PutSparseMatrix(DiffMatrix, "DiffMatrix");
                            bmc.Cmd("DiffMatrix = 0.5*(DiffMatrix + DiffMatrix');");
                            bmc.Cmd("condNoDiffMatrix = condest(DiffMatrix);");
                            bmc.Cmd("condNoFullMatrix = condest(FullMatrix);");
                            bmc.Cmd("eigiMaxi = eigs(DiffMatrix,1,'lm')");
                            bmc.Cmd("eigiMini = eigs(DiffMatrix,1,'sm')");
                            bmc.Cmd("lasterr");
                            bmc.Cmd("[V,r]=chol(DiffMatrix);");
                            bmc.Cmd("ret = [condNoDiffMatrix, condNoFullMatrix, eigiMaxi, eigiMini, r]");
                            bmc.GetMatrix(ret, "ret");

                            bmc.Execute(false);
                        }

                        double condNoDiffMatrix = ret[0, 0];
                        double condNoFullMatrix = ret[0, 1];
                        double eigiMaxi = ret[0, 2];
                        double eigiMini = ret[0, 3];
                        posDef = ret[0, 4] == 0;

                        Console.WriteLine("Eigenvalue range of diffusion matrix: {0} to {1}", eigiMini, eigiMaxi);

                        Console.WriteLine("Condition number diffusion operator: {0:0.####E-00}", condNoDiffMatrix);
                        Console.WriteLine("Condition number full operator: {0:0.####E-00}", condNoFullMatrix);

                    } else {
                        // +++++++++++++++++++++++++++++++++++++++
                        // test only for positive definiteness
                        // +++++++++++++++++++++++++++++++++++++++

                        var DiffMatrixFull = DiffMatrix.ToFullMatrixOnProc0();

                        posDef = true;
                        try {
                            DiffMatrixFull.Cholesky();
                        } catch (ArithmeticException) {
                            posDef = false;
                        }
                    }

                    double DiffSymm = DiffMatrix.SymmetryDeviation();
                    Console.WriteLine("Symmetry deviation of diffusion matrix: " + DiffSymm);

                    if (posDef)
                        Console.WriteLine("Good news: Diffusion operator matrix seems to be positive definite.");
                    else
                        Console.WriteLine("WARNING: Diffusion operator matrix is not positive definite.");

                    if (CheckAssertions) {
                        if (Control.AdvancedDiscretizationOptions.ViscosityMode == ViscosityMode.FullySymmetric && Control.PhysicalParameters.IncludeConvection == false) {
                            double compVal = DiffMatrix.InfNorm() * 1e-13;
                            Assert.LessOrEqual(DiffSymm, compVal, "Diffusion matrix seems to be non-symmetric.");
                            Assert.IsTrue(posDef, "Positive definiteness test failed.");
                        }
                    }
                }
            }
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {
            Tecplot.PlotFields(base.m_RegisteredFields, "XNSE_Solver" + timestepNo, physTime, superSampling);
            //Tecplot.PlotFields(new DGField[] { this.LevSet }, "grid" + timestepNo, physTime, 0);
        }


        protected void PlotOnIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, iterIndex }), 2);
        }


        protected override void SetInitial() {
            base.SetInitial();

            this.InitLevelSet();

            this.CreateEquationsAndSolvers(null);

            // =========================================
            // XDG BDF Timestepper initialization
            // =========================================

            if (m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DelayedTimestepperInit(0.0, 0, this.Control.GetFixedTimestep(),
                    // delegate for the initialization of previous timesteps from an analytic solution
                    BDFDelayedInitSetIntial);
            }

            After_SetInitialOrLoadRestart(0.0, 0);

        }

        /// <summary>
        /// delegate for the initialization of previous timesteps from an analytic solution
        /// </summary>
        /// <param name="TimestepIndex"></param>
        /// <param name="Time"></param>
        /// <param name="St"></param>
        private void BDFDelayedInitSetIntial(int TimestepIndex, double Time, DGField[] St) {
            using(new FuncTrace()) {
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


        private void After_SetInitialOrLoadRestart(double PhysTime, int TimestepNo)
        {

            // =============================================
            // LogFile initialization
            // =============================================  

            if (this.Control.TestMode == true)
            {
                LogQueryValue(PhysTime);
            }
            else
            {
                if (this.Control.LogValues != XNSE_Control.LoggingValues.none && this.CurrentSessionInfo.ID != Guid.Empty && base.MPIRank == 0)
                {
                    InitLogFile(this.CurrentSessionInfo.ID);
                    WriteLogLine(TimestepNo, PhysTime);
                }
            }
        }


        /// <summary>
        /// setUp for the Level set initialization (Level-set algorithm, continuity, conservation)
        /// </summary>
        private void InitLevelSet() {
            using (new FuncTrace())
            {

                // check level-set
                if (this.LevSet.L2Norm() == 0) {
                    throw new NotSupportedException("Level set is not initialized - norm is 0.0 - ALL cells will be cut, no gradient can be defined!");
                }

                // tracker needs to be updated to get access to the cut-cell mask
                this.LsTrk.UpdateTracker();

                // ==============================
                // level-set initialization
                // ==============================

                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 0 }), 3);

                #region Initialize Level Set Evolution Algorithm
                switch (this.Control.Option_LevelSetEvolution) {
                    case LevelSetEvolution.Fourier:
                        InitFourier();
                        break;
                    case LevelSetEvolution.None:
                        if (this.Control.AdvancedDiscretizationOptions.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Fourier) {
                            Fourier_LevSet = FourierLevelSetFactory.Build(this.Control.FourierLevSetControl);
                            Fourier_LevSet.ProjectToDGLevelSet(this.DGLevSet.Current, this.LsTrk);
                        } else {
                            goto default;
                        }
                        break;
                    case LevelSetEvolution.ExtensionVelocity: {
                            // Create ExtensionVelocity Motion Algorithm
                            this.DGLevSet.Current.Clear();
                            this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);
                            DGLevSetGradient.Gradient(1.0, DGLevSet.Current);
                            //VectorField<SinglePhaseField> VectorExtVel = ExtensionVelocity.Current;
                            base.RegisterField(ExtensionVelocity.Current);

                            ReInitPDE = new EllipticReInit(this.LsTrk, this.Control.ReInitControl, DGLevSet.Current);
                            FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);

                            // full initial reinitialization
                            ReInitPDE.ReInitialize(Restriction: LsTrk.Regions.GetNearFieldSubgrid(1));

                            CellMask Accepted = LsTrk.Regions.GetNearFieldMask(1);
                            CellMask ActiveField = Accepted.Complement();
                            CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                            FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            ReInitPDE.ReInitialize();

                            // setup extension velocity mover
                            switch (this.Control.Timestepper_Scheme) {
                                case XNSE_Control.TimesteppingScheme.RK_CrankNicolson:
                                case XNSE_Control.TimesteppingScheme.CrankNicolson: {
                                        //do not instantiate rksch, use bdf instead
                                        bdfOrder = -1;
                                        break;
                                    }
                                case XNSE_Control.TimesteppingScheme.RK_ImplicitEuler:
                                case XNSE_Control.TimesteppingScheme.ImplicitEuler: {
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
                                            throw new NotImplementedException();
                                    }
                            }

                            ExtVelMover = new ExtensionVelocityBDFMover(LsTrk, DGLevSet.Current, DGLevSetGradient, new VectorField<DGField>(XDGvelocity.Velocity.ToArray()), Control.EllipticExtVelAlgoControl, BcMap, bdfOrder, ExtensionVelocity.Current, new double[2] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B });


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

                        FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);

                        break;
                }
                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 1 }), 3);
                #endregion

                // =========================================
                // Enforcing the continuity of the level-set
                // =========================================

                ContinuityEnforcer = new ContinuityProjection(
                    DGLevelSet: this.DGLevSet.Current,
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

                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 4 }), 3);

                this.LsTrk.UpdateTracker();

            }

        }

        /// <summary>
        /// init routine for the specialized Fourier level-set
        /// </summary>
        private void InitFourier() {
            if (this.Control.FourierLevSetControl == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of FourierLevSetControl!");

            Fourier_LevSet = FourierLevelSetFactory.Build(this.Control.FourierLevSetControl);
            if (this.Control.EnforceLevelSetConservation) {
                throw new NotSupportedException("mass conservation correction currently not supported");
            }
            Fourier_LevSet.ProjectToDGLevelSet(this.DGLevSet.Current, this.LsTrk);

            if (base.MPIRank == 0 && this.CurrentSessionInfo.ID != Guid.Empty) {
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


        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            base.LoadRestart(out Time, out TimestepNo);

            //this.InitLevelSet();

            ContinuityEnforcer = new ContinuityProjection(
                    DGLevelSet: this.DGLevSet.Current,
                    gridData: GridData,
                    Option: Control.LSContiProjectionMethod
                    );

            //this.LsTrk.UpdateTracker();

            this.CreateEquationsAndSolvers(null);

            // =========================================
            // XDG BDF Timestepper initialization
            // =========================================

            if (m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DelayedTimestepperInit(Time, TimestepNo.MajorNumber, this.Control.GetFixedTimestep(),
                    // delegate for the initialization of previous timesteps from restart session
                    BDFDelayedInitLoadRestart );
            }

            After_SetInitialOrLoadRestart(Time, TimestepNo.MajorNumber);

        }

        /// <summary>
        /// delegate for the initialization of previous timesteps from restart session
        /// </summary>
        /// <param name="TimestepIndex"></param>
        /// <param name="time"></param>
        /// <param name="St"></param>
        private void BDFDelayedInitLoadRestart(int TimestepIndex, double time, DGField[] St) {

            Console.WriteLine("Timestep index {0}, time {1} ", TimestepIndex, time);

            ITimestepInfo tsi_toLoad;
            if(TimestepIndex < 0) {
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


        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 20 }), 2);

            // Load the sample Points for the restart of the Fourier LevelSet
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

        ResidualLogger m_EnergyLogger;

        /// <summary>
        /// Logger for kinetic and surface energy.
        /// </summary>
        ResidualLogger EnergyLogger {
            get {
                if (!this.Control.ComputeEnergy)
                    return null;

                if (m_EnergyLogger == null) {
                    m_EnergyLogger = new ResidualLogger(base.MPIRank, base.DatabaseDriver, base.CurrentSessionInfo.ID);
                    m_EnergyLogger.WriteResidualsToConsole = false;
                    m_EnergyLogger.WriteResidualsToTextFile = true;
                    m_EnergyLogger.TextFileFileName = "Energy";
                }

                return m_EnergyLogger;
            }
        }


        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionVelocity"/> and <see cref="XNSE_Control.ExactSolutionPressure"/>).
        /// </summary>
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
                        ConventionalDGField Vel_d;
                        if(this.CurrentVel[d] is XDGField) {
                            Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);
                        } else {
                            Vel_d = ((ConventionalDGField)this.CurrentVel[d]);
                        }
                        
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
                
            } //*/

            return Ret;
        }

        protected override void Bye() {
            base.Bye();
            if (EnergyLogger != null)
                EnergyLogger.Close();
        }



        /// <summary>
        /// refinement indicator
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel) {

            //int minRefineLevelLS = 0;
            //int maxRefineLevelLS = 1;

            CellMask ccm = this.LsTrk.Regions.GetCutCellMask();
            CellMask near = this.LsTrk.Regions.GetNearFieldMask(1);

            //double curv_max = 2.0 / this.GridData.Cells.h_min[j];

            int DesiredLevel_j = CurrentLevel;

            if (near.Contains(j)) {

                DesiredLevel_j = this.Control.RefinementLevel;

            } else {
                DesiredLevel_j = 0;
            }


            //if (ccm.Contains(j)) {

                //if (DesiredLevel_j < minRefineLevelLS) {
                //    // set minimum refinement level for the interface
                //    DesiredLevel_j = minRefineLevelLS;

                //} else {
                //    // further localized refinement

                //    // check for high curvature
                //    int DesiredLevelj_highCurv = DesiredLevel_j;
                //    this.Curvature.GetExtremalValuesInCell(out double curv_jMin, out double curv_jMax, j);
                //    if ((curv_jMax >= curv_max || Math.Abs(curv_jMin) >= curv_max) && DesiredLevel_j < maxRefineLevelLS) {
                //        DesiredLevelj_highCurv++;
                //    } else if ((curv_jMax < curv_max / 2) || (Math.Abs(curv_jMin) < curv_max / 2)) {
                //        DesiredLevelj_highCurv--;
                //    }

                //    //double mean_curv = Math.Abs(this.Curvature.GetMeanValue(j));
                //    //if ((mean_curv >= curv_max) && CurrentLevel < maxRefineLevelLS)
                //    //    DesiredLevel_j = CurrentLevel + 1;

                //    // check for small cut cells
                //    int DesiredLevelj_agglom = DesiredLevel_j;
                //    double cellVol = this.GridData.Cells.GetCellVolume(j);
                //    var spcIds = this.LsTrk.SpeciesIdS.ToArray();
                //    double ratioVolSpcMin = 1.0;
                //    foreach (SpeciesId spc in this.LsTrk.SpeciesIdS) {
                //        double cellVolSpc = this.LsTrk.GetXDGSpaceMetrics(spcIds, m_HMForder, 1).CutCellMetrics.CutCellVolumes[spc][j];
                //        double ratioVolSpc = cellVolSpc / cellVol;
                //        if (ratioVolSpc < ratioVolSpcMin)
                //            ratioVolSpcMin = ratioVolSpc;
                //    }
                //    double thrshld = this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold;
                //    if (ratioVolSpcMin < thrshld && DesiredLevel_j < maxRefineLevelLS) {
                //        DesiredLevelj_agglom++;
                //    } else if (ratioVolSpcMin > 4 * thrshld) {
                //        DesiredLevelj_agglom--;
                //    }

                //    // check for a change of sign in the curvature
                //    int DesiredLevelj_inflection = DesiredLevel_j;
                //    //this.Curvature.GetExtremalValuesInCell(out double curv_jMin, out double curv_jMax, j);
                //    if (Math.Sign(curv_jMin) != Math.Sign(curv_jMax) && DesiredLevel_j < maxRefineLevelLS) 
                //        DesiredLevelj_inflection++;

                //    DesiredLevel_j = (new int[] { DesiredLevelj_highCurv, DesiredLevelj_agglom, DesiredLevelj_inflection }).Max();

                //}

            //} else {
            //    // non cut cells don't need to be refined
            //    DesiredLevel_j = 0;
            //}

            return DesiredLevel_j;

        }

        //CellMask refinedInterfaceCells;

        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {

            if (this.Control.AdaptiveMeshRefinement) {

                //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 0), 2);

                // Check grid changes
                // ==================

                CellMask BlockedCells = LsTrk.Regions.GetNearFieldMask(1);

                // compute curvature for levelindicator 
                //CurvatureAlgorithms.CurvatureDriver(
                //    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                //    CurvatureAlgorithms.FilterConfiguration.NoFilter,
                //    this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                //    this.m_HMForder, this.DGLevSet.Current);

                //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 1), 2);

                bool AnyChange = GridRefinementController.ComputeGridChange(this.GridData, BlockedCells, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
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

                    newGrid = this.GridData.Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                    //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 2 }), 2);#

                } else {

                    newGrid = null;
                    old2NewGrid = null;
                }

            } else {

                newGrid = null;
                old2NewGrid = null;
            }
        }


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
        /// only the velocity components(supposed to be at the beginning)  are used.
        /// </param>
        /// <param name="underrelax">
        /// </param>
        double DelUpdateLevelSet(DGField[] CurrentState, double Phystime, double dt, double underrelax, bool incremental) {
            using (new FuncTrace()) {

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

                ConventionalDGField[] meanVelocity;

                if (EvoVelocity[0] is XDGField) {
                    meanVelocity = GetMeanVelocityFromXDGField(EvoVelocity);
                }
                else if (EvoVelocity[0] is ConventionalDGField) {
                    // +++++++++++++++++
                    // plain DG velocity 
                    // +++++++++++++++++
                    Debug.Assert(this.XDGvelocity == null);
                    Debug.Assert(this.DGvelocity != null);  // make sure that we *do not* run velocity enrichment

                    meanVelocity = EvoVelocity.Select(v => ((ConventionalDGField)v)).ToArray();
                } else {
                    throw new ApplicationException("Should not happen.");
                }

                #endregion


                // ==================================================================
                // backup interface properties (mass conservation, surface changerate)
                // ==================================================================

                #region backup interface props

                double oldSurfVolume = 0.0;
                double oldSurfLength = 0.0;
                double SurfChangerate = 0.0;
                if (this.Control.CheckInterfaceProps) {
                    oldSurfVolume = XNSEUtils.GetSpeciesArea(this.LsTrk, LsTrk.GetSpeciesId("A"));
                    oldSurfLength = XNSEUtils.GetInterfaceLength(this.LsTrk);
                    SurfChangerate = XNSEUtils.GetSurfaceChangerate(this.LsTrk, meanVelocity, this.m_HMForder);
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
                switch (this.Control.Option_LevelSetEvolution) {
                    case LevelSetEvolution.None:
                    throw new ArgumentException("illegal call");

                    case LevelSetEvolution.FastMarching: {

                        NarrowMarchingBand.Evolve_Mk2(
                             dt, this.LsTrk, DGLevSet_old, this.DGLevSet.Current, this.DGLevSetGradient,
                             meanVelocity, this.ExtensionVelocity.Current.ToArray(),
                             this.m_HMForder, iTimestep);

                            //CellMask Accepted = LsTrk.Regions.GetCutCellMask();
                            //CellMask ActiveField = LsTrk.Regions.GetNearFieldMask(1);
                            //CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                            //FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            break;
                    }

                    case LevelSetEvolution.Fourier: {
                        Fourier_Timestepper.moveLevelSet(dt, meanVelocity);
                        if (incremental)
                            Fourier_Timestepper.updateFourierLevSet();
                        Fourier_LevSet.ProjectToDGLevelSet(this.DGLevSet.Current, this.LsTrk);
                        break;
                    }

                    case LevelSetEvolution.Prescribed: {
                        this.DGLevSet.Current.Clear();
                        this.DGLevSet.Current.ProjectField(1.0, Control.Phi.Vectorize(Phystime + dt));
                        break;
                    }

                    case LevelSetEvolution.ScalarConvection: {
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
                    case LevelSetEvolution.ExtensionVelocity: {

                            DGLevSetGradient.Clear();
                            DGLevSetGradient.Gradient(1.0, DGLevSet.Current);

                            ExtVelMover.Advect(dt);

                            //ReInitPDE.ReInitialize(Restriction: LsTrk.Regions.GetNearFieldSubgrid(1));

                            // Fast Marching: Specify the Domains first
                            // Perform Fast Marching only on the Far Field
                            CellMask Accepted =  LsTrk.Regions.GetNearFieldMask(1);
                            CellMask ActiveField = Accepted.Complement();
                            CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                            FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            //ReInitPDE.ReInitialize();
                            
                            break;
                        }
                    default:
                        throw new ApplicationException();
                }


                // performing underrelaxation
                if (underrelax < 1.0) {
                    this.DGLevSet.Current.Scale(underrelax);
                    this.DGLevSet.Current.Acc((1.0 - underrelax), DGLevSet_oldIter);
                }

                //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 1 }), 2);

                #endregion


                // ======================
                // postprocessing  
                // =======================

                #region ensure continuity

                // make level set continuous
                //CellMask CC = LsTrk.Regions.GetCutCellMask4LevSet(0);
                CellMask Near1 = LsTrk.Regions.GetNearMask4LevSet(0, 1);
                CellMask PosFF = LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
                ContinuityEnforcer.MakeContinuous(this.DGLevSet.Current, this.LevSet, Near1, PosFF);

                if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.FastMarching) {
                    this.DGLevSet.Current.Clear(Near1);
                    this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet, Near1);
                }

                //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 2 }), 2);

                #endregion


                for (int d = 0; d < D; d++)
                    this.XDGvelocity.Velocity[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;



                // ===============
                // tracker update
                // ===============

                this.LsTrk.UpdateTracker(incremental: true);


                // ==================================================================
                // check interface properties (mass conservation, surface changerate)
                // ==================================================================

                if (this.Control.CheckInterfaceProps) {

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


        private void EnforceVolumeConservation() {
            double spcArea = XNSEUtils.GetSpeciesArea(LsTrk, LsTrk.SpeciesIdS[0]);
            Console.WriteLine("area = {0}", spcArea);
            double InterLength = XNSEUtils.GetInterfaceLength(LsTrk);

            //double cmc = (consvRefArea - spcArea) / InterLength;
            //Console.WriteLine("add constant: {0}", -cmc);
            //this.DGLevSet.Current.AccConstant(-cmc);
            //this.LevSet.AccConstant(-cmc);
        }


        /// <summary>
        ///  Take density-weighted mean value in cut-cells
        /// </summary>
        /// <param name="EvoVelocity"></param>
        /// <returns></returns>
        private ConventionalDGField[] GetMeanVelocityFromXDGField(DGField[] EvoVelocity) {
            int D = EvoVelocity.Length;
            ConventionalDGField[] meanVelocity;

            Debug.Assert(this.XDGvelocity != null);
            Debug.Assert(this.DGvelocity == null);  // make sure that we run velocity enrichment

            meanVelocity = new ConventionalDGField[D];

            double rho_A = this.Control.PhysicalParameters.rho_A, rho_B = this.Control.PhysicalParameters.rho_B;
            CellMask CC = this.LsTrk.Regions.GetCutCellMask4LevSet(0);
            CellMask Neg = this.LsTrk.Regions.GetLevelSetWing(0, -1).VolumeMask;
            CellMask Pos = this.LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
            CellMask posNear = this.LsTrk.Regions.GetNearMask4LevSet(0, 1).Except(Neg);
            CellMask negNear = this.LsTrk.Regions.GetNearMask4LevSet(0, 1).Except(Pos);

            for (int d = 0; d < D; d++) {
                Basis b = this.XDGvelocity.Velocity[d].Basis.NonX_Basis;
                meanVelocity[d] = new SinglePhaseField(b);


                foreach (string spc in this.LsTrk.SpeciesNames) {
                    double rhoSpc;
                    switch (spc) {
                        case "A": rhoSpc = rho_A; break;
                        case "B": rhoSpc = rho_B; break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }

                    double scale = rhoSpc / (rho_A + rho_B);

                    meanVelocity[d].Acc(scale, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), CC);
                    switch (spc) {
                        //case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), Neg.Except(CC)); break;
                        case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), negNear); break;
                        case "B": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), posNear); break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }
                }

            }

            return meanVelocity;
        }

    }
}
