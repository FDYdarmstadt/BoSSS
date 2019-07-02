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
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator.Energy;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.Grid.Aggregation;
using NUnit.Framework;
using MPI.Wrappers;
using System.Collections;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Solver for Incompressible Multiphase flows
    /// </summary>
    public class XNSE_SolverMain : BoSSS.Solution.Application<XNSE_Control> {



        static void Main(string[] args) {

            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.TestFixtureSetUp();
            ////BoSSS.Application.XNSE_Solver.Tests.UnitTest.PolynomialTestForConvectionTest(3, 0, false);
            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.TestCapillaryWave();
            ////BoSSS.Application.XNSE_Solver.Tests.ElementalTestProgramm.LineMovementTest(LevelSetEvolution.ScalarConvection, LevelSetHandling.Coupled_Once, XNSE_Control.TimesteppingScheme.ImplicitEuler, 0.5);
            //Assert.IsFalse(true, "remove me");


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
        XDGField divVelocity;

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
                    gridData: (GridData)GridData,
                    Option: Control.LSContiProjectionMethod
                    );

                this.LsTrk = new LevelSetTracker((GridData) this.GridData, base.Control.CutCellQuadratureType, base.Control.LS_TrackerWidth, new string[] { "A", "B" }, this.LevSet);
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


                if(this.Control.solveCoupledHeatEquation) {
                    this.Temperature = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Temperature].Degree), VariableNames.Temperature);
                    base.RegisterField(this.Temperature);
                    this.ResidualHeat = new XDGField(this.Temperature.Basis, "ResidualHeat");
                    base.RegisterField(this.ResidualHeat);

                    this.Heatflux = new VectorField<XDGField>(D.ForLoop(d => new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Temperature].Degree), "Heatflux_[" + d + "]")));
                    base.RegisterField(this.Heatflux);

                    this.DisjoiningPressure = new SinglePhaseField(new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree), "DisjoiningPressure");
                    if(this.Control.DisjoiningPressureFunc != null) {
                        DisjoiningPressure.ProjectField(this.Control.DisjoiningPressureFunc);
                    }
                    base.RegisterField(this.DisjoiningPressure);
                }


                XDGBasis b = new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree);
                this.divVelocity = new XDGField(b, "DivergenceVelocity");
                base.RegisterField(this.divVelocity);

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

                if(this.Control.ComputeEnergy) {

                    this.DerivedKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions["KineticEnergy"].Degree)), "DerivedKineticEnergy");
                    base.RegisterField(this.DerivedKineticEnergy);

                    this.GeneratedKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions["KineticEnergy"].Degree)), "KineticEnergyProduction");
                    base.RegisterField(this.GeneratedKineticEnergy);

                    this.KineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions["KineticEnergy"].Degree)), "KineticEnergy");
                    base.RegisterField(this.KineticEnergy);
                    this.ResidualKineticEnergy = new XDGField(this.KineticEnergy.Basis, "ResidualKineticEnergy");
                    base.RegisterField(this.ResidualKineticEnergy);

                    this.prevKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions["KineticEnergy"].Degree)));

                    this.KineticDissipation = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions["KineticEnergy"].Degree)), "KineticDissipation");
                    base.RegisterField(this.KineticDissipation);

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
        //XNSE_OperatorFactory XNSE_Operator;
        XNSFE_OperatorFactory XNSFE_Operator;

        /// <summary>
        /// OperatorConfiguration for the <see cref="XNSE_Operator"/>
        /// </summary>
        XNSFE_OperatorConfiguration XOpConfig;

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
            if (XNSFE_Operator != null && L == null)
                return;


            if (Control.CompMode == AppControl._CompMode.Steady) {
                if (Control.Timestepper_LevelSetHandling != LevelSetHandling.None)
                    throw new ApplicationException(string.Format("Illegal control file: for a steady computation ({0}), the level set handling must be {1}.", AppControl._CompMode.Steady, LevelSetHandling.None));
            }

            int degU = this.CurrentVel[0].Basis.Degree;

            //if(base.Control.FakePoisson) {
            //    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //    Console.WriteLine("ACHTUNG: Fake-Poisson aktiviert!");
            //    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //    Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //}

            //if(this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SemiImplicit)
            //    this.Control.PhysicalParameters.mu_I = this.Control.dtFixed * this.Control.PhysicalParameters.Sigma; //--> added to XNSE-operator config

            #endregion

            #region Config and Generate XOperator


            //Quadrature Order
            //----------------

            m_HMForder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);


            // Create Spatial Operator
            // ======================= 

            XOpConfig = new XNSFE_OperatorConfiguration(this.Control);

            XNSFE_Operator = new XNSFE_OperatorFactory(XOpConfig, this.LsTrk, this.m_HMForder, this.BcMap, degU);


            // kinetic energy balance Operator
            // ===============================

            if(this.Control.ComputeEnergy) {
                this.generateKinEnergyOperator();
            }


            // coupled heat Operator
            // =====================

            if(this.Control.solveCoupledHeatEquation) {
                this.generateCoupledOperator();
            }

            #endregion

            #region Create Timestepper
            // =======================
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

                //Herr Smuda bitte anpassen ...


                //if no Runge Kutta Timesteper is initialized
                // we are using bdf or screwed things up
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
                    //m_BDF_Timestepper.config_LinearSolver.SolverCode = LinearSolverConfig.Code.classic_mumps;
                    //m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;


                    // solver 
                    //m_BDF_Timestepper.config_NonLinearSolver.ConvergenceCriterion= this.Control.Solver_ConvergenceCriterion;
                    m_BDF_Timestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                    //m_BDF_Timestepper.config_NonLinearSolver.MaxSolverIterations = this.Control.Solver_MaxIterations;
                    
                    this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations; //m_BDF_Timestepper.config_NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;

                    //m_BDF_Timestepper.config_NonLinearSolver.MaxKrylovDim = this.Control.Solver_MaxKrylovDim;
                    //m_BDF_Timestepper.config_NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Picard;
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

                    if(this.Control.ComputeEnergy) {
                        m_BDF_energyTimestepper = new XdgBDFTimestepping(
                        this.CurrentEnergySolution.Mapping.Fields,
                        this.CurrentEnergyResidual.Mapping.Fields,
                        LsTrk,
                        false,
                        DelComputeEnergyOperatorMatrix, null, DelUpdateLevelSet_EnergyOperator,
                        (this.Control.CompMode == AppControl._CompMode.Transient) ? bdfOrder : 1,
                        this.Control.Timestepper_LevelSetHandling,
                        MassMatrixShapeandDependence.IsTimeDependent,   // only for Lie-Splitting and coupled_Once
                        SpatialOperatorType.LinearTimeDependent,
                        MassScaleForEnergy,
                        this.MultigridEnergyOperatorConfig, base.MultigridSequence,
                        this.LsTrk.SpeciesIdS.ToArray(), this.KineticEnergy.Basis.Degree * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2),
                        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                        true,
                        this.Control.NonLinearSolver,
                        this.Control.LinearSolver
                        );
                        m_BDF_energyTimestepper.m_ResLogger = this.EnergyResLogger;
                        m_BDF_energyTimestepper.m_ResidualNames = this.CurrentEnergyResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                        //m_BDF_coupledTimestepper.Config_SolverConvergenceCriterion = this.Control.Solver_ConvergenceCriterion;
                        m_BDF_energyTimestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                        //m_BDF_coupledTimestepper.Config_MaxIterations = this.Control.Solver_MaxIterations;
                        this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations;
                        //m_BDF_energyTimestepper.Config_MinIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;
                        m_BDF_energyTimestepper.Timestepper_Init = TimeStepperInit.SingleInit;
                        m_BDF_energyTimestepper.PushLevelSet = delegate () { };    // dummy push
                        m_BDF_energyTimestepper.coupledOperator = true;
                    }

                    if(this.Control.solveCoupledHeatEquation) {
                        m_BDF_coupledTimestepper = new XdgBDFTimestepping(
                        this.CurrentCoupledSolution.Mapping.Fields,
                        this.CurrentCoupledResidual.Mapping.Fields,
                        LsTrk,
                        false,
                        DelComputeCoupledOperatorMatrix, null, DelUpdateLevelSet_CoupledOperator,
                        (this.Control.CompMode == AppControl._CompMode.Transient) ? bdfOrder : 1,
                        this.Control.Timestepper_LevelSetHandling,
                        MassMatrixShapeandDependence.IsTimeDependent,   // only for Lie-Splitting and coupled_Once
                        SpatialOperatorType.LinearTimeDependent,
                        HeatScale,
                        this.MultigridCoupledOperatorConfig, base.MultigridSequence,
                        this.LsTrk.SpeciesIdS.ToArray(), m_HMForder,
                        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                        true,
                        this.Control.NonLinearSolver,
                        this.Control.LinearSolver
                        );           
                        m_BDF_coupledTimestepper.m_ResLogger = this.CouplededResLogger;
                        m_BDF_coupledTimestepper.m_ResidualNames = this.CurrentCoupledResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                        //m_BDF_coupledTimestepper.Config_SolverConvergenceCriterion = this.Control.Solver_ConvergenceCriterion;
                        m_BDF_coupledTimestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                        //m_BDF_coupledTimestepper.Config_MaxIterations = this.Control.Solver_MaxIterations;

                        this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations; 
                        //m_BDF_coupledTimestepper.Config_MinIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;

                        m_BDF_coupledTimestepper.Timestepper_Init = TimeStepperInit.SingleInit;
                        m_BDF_coupledTimestepper.PushLevelSet = delegate() { };    // dummy push
                        m_BDF_coupledTimestepper.coupledOperator = true;

                        m_BDF_coupledTimestepper.XdgSolverFactory.Update(this.Control.NonLinearSolver, this.Control.LinearSolver); //do not forget to update your changes to Solver configurations!
                    }

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

                //PlotCurrentState(hack_Phystime, new TimestepNumber(hack_TimestepIndex, 12), 2);

                Debug.Assert(object.ReferenceEquals(this.MultigridSequence[0].ParentGrid, this.GridData));

                m_BDF_Timestepper.DataRestoreAfterBalancing(L, 
                    ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure), 
                    ArrayTools.Cat<DGField>(this.XDGvelocity.ResidualMomentum.ToArray(), this.ResidualContinuity), 
                    this.LsTrk, this.MultigridSequence);

                if(this.Control.solveCoupledHeatEquation)
                    m_BDF_coupledTimestepper.DataRestoreAfterBalancing(L,
                          this.Temperature.ToEnumerable(),
                          this.ResidualHeat.ToEnumerable(),
                          this.LsTrk, this.MultigridSequence);

                //PlotCurrentState(hack_Phystime, new TimestepNumber(hack_TimestepIndex, 13), 2);

                ContinuityEnforcer = new ContinuityProjection(
                    ContBasis: this.LevSet.Basis, 
                    DGBasis: this.DGLevSet.Current.Basis, 
                    gridData: GridData, 
                    Option: Control.LSContiProjectionMethod);

                if(this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
                    ReInitPDE = new EllipticReInit(this.LsTrk, this.Control.ReInitControl, DGLevSet.Current);
                    FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                    ExtVelMover = new ExtensionVelocityBDFMover(LsTrk, DGLevSet.Current, DGLevSetGradient, new VectorField<DGField>(XDGvelocity.Velocity.ToArray()),
                    Control.EllipticExtVelAlgoControl, BcMap, bdfOrder, ExtensionVelocity.Current, new double[2] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B });
                }

            }
            #endregion

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
                        if(this.Control.AdvancedDiscretizationOptions.CurvatureNeeded || XOpConfig.isEvaporation) {
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

            this.XNSFE_Operator.AssembleMatrix(
                OpMtx, OpAffine, codMap, domMap,
                CurrentState, AgglomeratedCellLengthScales, phystime,
                this.m_HMForder, SurfaceForce, filtLevSetGradient, Curvature,
                (this.Control.solveCoupledHeatEquation ? this.Temperature.ToEnumerable() : null),
                (this.Control.solveCoupledHeatEquation ? this.DisjoiningPressure.ToEnumerable() : null));


            if(filtLevSetGradient != null) {
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
                            rtsi = new TimestepInfo(
                                t - ti * this.Control.GetFixedTimestep(),
                                this.CurrentSessionInfo,
                                tsn,
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
            using(var tr = new FuncTrace()) {

                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;
                base.ResLogger.TimeStep = TimestepInt;
                hack_TimestepIndex = TimestepInt;
                hack_Phystime = phystime;

  
                Preprocessing(TimestepInt, phystime, dt, TimestepNo);


                if(Control.SkipSolveAndEvaluateResidual) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    // setup: project exact solution -- for consistency tests
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    foreach(string spc in LsTrk.SpeciesNames) {
                        for(int d = 0; d < this.GridData.SpatialDimension; d++) {
                            ConventionalDGField Vel_d;
                            if(this.CurrentVel[d] is XDGField)
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


                if(base.Control.CompMode == AppControl._CompMode.Steady) {
                    dt = 1.0e100;
                    Console.WriteLine("Steady-state solve ...", TimestepNo, dt);

                    if(this.Control.Option_LevelSetEvolution != LevelSetEvolution.None) {
                        throw new ApplicationException("For steady-state solutions, the only allowed level-set-evolution option is '" + LevelSetEvolution.None + "'.");
                    }



                    // =====================================================
                    // setup transient 
                    // =====================================================
                } else if(base.Control.CompMode == AppControl._CompMode.Transient) {

                    // push stacks
                    // -----------

                    PushLevelSetAndRelatedStuff();


                    // backup old velocity for energy checks
                    // -------------------------------------
                    if(this.Control.ComputeEnergy && this.Control.CompMode == AppControl._CompMode.Transient) {
                        for(int d = 0; d < D; d++) {
                            this.prevVel[d].Clear();
                            this.prevVel[d].Acc(1.0, this.CurrentVel[d]);
                        }
                    }


                    // fields setup
                    // ------------
                    if(this.XDGvelocity != null) {
                        for(int d = 0; d < D; d++) {
                            // Gravity must be set up like this to avoid regions of zero gravity when updating the level-set
                            this.XDGvelocity.Gravity[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
                        }
                    } else {
                        throw new NotSupportedException();
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
                    if(this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative && this.Control.LSunderrelax == 1.0) {
                        if(dt / dt_LevSetCFL > 1.0) {
                            double underrelax = Math.Round(dt_LevSetCFL / dt, 1);
                            m_BDF_Timestepper.IterUnderrelax = underrelax;
                            Console.WriteLine("Exceeding Level-Set CFL: Setting underrelaxation factor to {0}", underrelax);
                        } else {
                            m_BDF_Timestepper.IterUnderrelax = 1.0;
                        }
                    }


                    //dt = Math.Min(dt, dt_LevSetCFL);

                    // Capillary Timestep restriction
                    if(this.Control.PhysicalParameters.Sigma != 0.0) {
                        MultidimensionalArray h_mins = ((GridData)this.GridData).Cells.h_min;
                        double h = h_mins.Min();
                        double LevSet_Deg = this.LevSet.Basis.Degree + 1;
                        h /= LevSet_Deg;
                        double dt_sigma = Math.Sqrt((this.Control.PhysicalParameters.rho_A + this.Control.PhysicalParameters.rho_B)
                            * Math.Pow(h, 3) / (2 * Math.PI * this.Control.PhysicalParameters.Sigma));
                        if(dt > dt_sigma)
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

                using(new BlockTrace("Solve", tr)) {

                    if(m_BDF_Timestepper != null) {
                        m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

                        if(this.Control.solveCoupledHeatEquation && m_BDF_coupledTimestepper != null) {
                            m_BDF_coupledTimestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
                            ComputeHeatflux();
                        }
                        
                    } else {
                        //m_RK_Timestepper.Solve(phystime, dt);
                    }

                }


                if(this.Control.ComputeEnergy && m_BDF_energyTimestepper != null) {

                    this.prevKineticEnergy.Clear();
                    this.prevKineticEnergy.Acc(1.0, this.KineticEnergy);

                    // solve kinetic energy balance
                    m_BDF_energyTimestepper.Solve(phystime, dt);

                    // derive kinetic Energy from flow solution
                    double[] rhoS = new double[] { this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B };
                    EnergyUtils.ProjectKineticEnergy(this.DerivedKineticEnergy, this.LsTrk, this.XDGvelocity.Velocity.ToArray(), rhoS, this.m_HMForder);

                    // compute generated kinetic energy
                    GeneratedKineticEnergy.Clear();
                    GeneratedKineticEnergy.Acc(1.0, this.KineticEnergy);
                    GeneratedKineticEnergy.Acc(-1.0, this.DerivedKineticEnergy);

                    // changerate of kinetic energy from discretization 
                    double[] muS = new double[] { this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B };
                    EnergyUtils.ProjectKineticDissipation(this.KineticDissipation, this.LsTrk, this.XDGvelocity.Velocity.ToArray(), muS, this.m_HMForder);

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
                        int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
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
                if (this.Control.LogValues != XNSE_Control.LoggingValues.None && this.CurrentSessionInfo.ID != Guid.Empty && base.MPIRank == 0)
                {
                    InitLogFile(this.CurrentSessionInfo.ID);
                    WriteLogLine(TimestepNo, PhysTime);
                }
            }
        }


        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            base.LoadRestart(out Time, out TimestepNo);

            //this.InitLevelSet();
            if(this.Control.ReInitPeriod > 0) {
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

            if(this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
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

            // hack for restart
            //Console.WriteLine("Warning: hack for restart!");
            //this.XDGvelocity.Gravity[1].GetSpeciesShadowField("A").AccConstant(-9.81e2);
            //this.XDGvelocity.Gravity[1].GetSpeciesShadowField("B").AccConstant(-9.81e2);

            if(this.Control.ClearVelocitiesOnRestart) {
                Console.WriteLine("clearing all velocities");
                this.XDGvelocity.Velocity.Clear();
            }                
             

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


        protected override void Bye() {
            base.Bye();
            if (EnergyLogger != null)
                EnergyLogger.Close();
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

            if(this.Control.BaseRefinementLevel == 0)
                return 0;

            CellMask ccm = this.LsTrk.Regions.GetCutCellMask();
            CellMask near = this.LsTrk.Regions.GetNearFieldMask(1);
            CellMask nearBnd = near.AllNeighbourCells();
            CellMask buffer = nearBnd.AllNeighbourCells().Union(nearBnd).Except(near);


            int DesiredLevel_j = CurrentLevel;

            if(near.Contains(j)) {
                if(CurrentLevel < this.Control.BaseRefinementLevel) {
                    DesiredLevel_j++;
                } else {
                    // additional refinement
                    switch(this.Control.RefineStrategy) {
                        case XNSE_Control.RefinementStrategy.CurvatureRefined: {
                                double curv_max = 1.0 / (2.0 * ((GridData)this.GridData).Cells.h_min[j]);
                                double mean_curv = Math.Abs(this.Curvature.GetMeanValue(j));
                                double minCurv, maxCurv;
                                this.Curvature.GetExtremalValuesInCell(out minCurv, out maxCurv, j);
                                double max_AbsCurv = Math.Max(Math.Abs(minCurv), Math.Abs(maxCurv));

                                double curv_thrshld = mean_curv;
                                if(curv_thrshld > curv_max && CurrentLevel == this.Control.RefinementLevel) {
                                    DesiredLevel_j++;
                                } else if(curv_thrshld < (curv_max / 2) && CurrentLevel == this.Control.RefinementLevel + 1) {
                                    DesiredLevel_j--;
                                }
                                break;
                            }
                        case XNSE_Control.RefinementStrategy.ContactLineRefined: {
                                CellMask BCells = ((GridData)this.GridData).BoundaryCells.VolumeMask;
                                if(ccm.Contains(j) && BCells.Contains(j) && CurrentLevel < this.Control.RefinementLevel) {
                                    DesiredLevel_j++;
                                } else if(!BCells.Contains(j)) { // && CurrentLevel == this.Control.RefinementLevel + 1) {
                                    DesiredLevel_j--;
                                }
                                break;
                            }
                        case XNSE_Control.RefinementStrategy.constantInterface:
                        default:
                            break;
                    }
                }

            } else if(NScm.Contains(j)) {
                if(CurrentLevel < this.Control.BaseRefinementLevel)
                    DesiredLevel_j++;

            } else if(buffer.Contains(j) || NSbuffer.Contains(j)) {
                if(CurrentLevel < this.Control.BaseRefinementLevel - 1)
                    DesiredLevel_j++;
            } else {
                DesiredLevel_j = 0;
            }

            return DesiredLevel_j;

        }


        //int LevelIndicator(int j, int CurrentLevel) {

        //    CellMask spc = this.LsTrk.Regions.GetSpeciesMask("A");

        //    if(spc.Contains(j)) {
        //        return 2;
        //    } else {
        //        return 0;
        //    }

        //}


        /// <summary>
        /// refinement indicator
        /// </summary>
        //int LevelIndicator(int j, int CurrentLevel) {

        //    int minRefineLevelLS = 1;
        //    int maxRefineLevelLS = 2;

        //    CellMask ccm = this.LsTrk.Regions.GetCutCellMask();
        //    CellMask near = this.LsTrk.Regions.GetNearFieldMask(1);

        //    double curv_max = 1.0 / this.GridData.Cells.h_min[j];

        //    int DesiredLevel_j = CurrentLevel;

        //    if(near.Contains(j)) {

        //        if(DesiredLevel_j < minRefineLevelLS) {
        //            // set minimum refinement level for the interface
        //            DesiredLevel_j = minRefineLevelLS;

        //        } else if (ccm.Contains(j)) {
        //            // further localized refinement

        //            // check for high curvature
        //            //int DesiredLevelj_highCurv = DesiredLevel_j;
        //            //this.Curvature.GetExtremalValuesInCell(out double curv_jMin, out double curv_jMax, j);
        //            //if((curv_jMax >= curv_max || Math.Abs(curv_jMin) >= curv_max) && DesiredLevel_j < maxRefineLevelLS) {
        //            //    DesiredLevelj_highCurv++;
        //            //} else if((curv_jMax < curv_max / 2) || (Math.Abs(curv_jMin) < curv_max / 2)) {
        //            //    DesiredLevelj_highCurv--;
        //            //}

        //            double mean_curv = Math.Abs(this.Curvature.GetMeanValue(j));
        //            if((mean_curv >= curv_max) && CurrentLevel < maxRefineLevelLS) {
        //                DesiredLevel_j++;
        //            } else if(mean_curv < curv_max / 2 && CurrentLevel > minRefineLevelLS) {
        //                DesiredLevel_j--;
        //            }

        //            //// check for small cut cells
        //            //int DesiredLevelj_agglom = DesiredLevel_j;
        //            //double cellVol = this.GridData.Cells.GetCellVolume(j);
        //            //var spcIds = this.LsTrk.SpeciesIdS.ToArray();
        //            //double ratioVolSpcMin = 1.0;
        //            //foreach(SpeciesId spc in this.LsTrk.SpeciesIdS) {
        //            //    double cellVolSpc = this.LsTrk.GetXDGSpaceMetrics(spcIds, m_HMForder, 1).CutCellMetrics.CutCellVolumes[spc][j];
        //            //    double ratioVolSpc = cellVolSpc / cellVol;
        //            //    if(ratioVolSpc < ratioVolSpcMin)
        //            //        ratioVolSpcMin = ratioVolSpc;
        //            //}
        //            //double thrshld = this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold;
        //            //if(ratioVolSpcMin < thrshld && DesiredLevel_j < maxRefineLevelLS) {
        //            //    DesiredLevelj_agglom++;
        //            //} else if(ratioVolSpcMin > 4 * thrshld) {
        //            //    DesiredLevelj_agglom--;
        //            //}

        //            //// check for a change of sign in the curvature
        //            //int DesiredLevelj_inflection = DesiredLevel_j;
        //            ////this.Curvature.GetExtremalValuesInCell(out double curv_jMin, out double curv_jMax, j);
        //            //if(Math.Sign(curv_jMin) != Math.Sign(curv_jMax) && DesiredLevel_j < maxRefineLevelLS)
        //            //    DesiredLevelj_inflection++;

        //            //DesiredLevel_j = (new int[] { DesiredLevelj_highCurv, DesiredLevelj_agglom, DesiredLevelj_inflection }).Max();

        //        }

        //    } else {
        //        // non cut cells don't need to be refined
        //        DesiredLevel_j = 0;
        //    }

        //    return DesiredLevel_j;

        //}

        //CellMask refinedInterfaceCells;

        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            using(new FuncTrace()) {

                if(this.Control.AdaptiveMeshRefinement) {

                    //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 0), 2);

                    // Check grid changes
                    // ==================

                    CellMask BlockedCells;
                    if(this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Once
                        || this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                        int prevInd = LsTrk.PopulatedHistoryLength;
                        CellMask prevNear = LsTrk.RegionsHistory[-prevInd + 1].GetNearFieldMask(1);
                        BlockedCells = (TimestepNo > 1) ? prevNear : null; // CellMask.Union(currNear, prevNear);
                    } else {
                        CellMask currNear = LsTrk.Regions.GetNearFieldMask(1);
                        BlockedCells = currNear;
                    }

                    // compute curvature for levelindicator 
                    if(this.Control.RefineStrategy == XNSE_Control.RefinementStrategy.CurvatureRefined) {
                        CurvatureAlgorithms.CurvatureDriver(
                            SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint,
                            CurvatureAlgorithms.FilterConfiguration.Default,
                            this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                            this.m_HMForder, this.DGLevSet.Current);
                    }


                    // navier slip boundary cells
                    NScm = new CellMask(this.GridData);
                    NSbuffer = new CellMask(this.GridData);
                    if(this.Control.RefineNavierSlipBoundary) {
                        BitArray NSc = new BitArray(((GridData)this.GridData).Cells.Count);
                        CellMask bnd = ((GridData)this.GridData).BoundaryCells.VolumeMask;
                        int[][] c2e = ((GridData)this.GridData).Cells.Cells2Edges;
                        foreach(Chunk cnk in bnd) {
                            for(int i = cnk.i0; i < cnk.JE; i++) {
                                foreach(int e in c2e[i]) {
                                    int eId = (e < 0) ? -e - 1 : e - 1;
                                    byte et = ((GridData)this.GridData).Edges.EdgeTags[eId];
                                    if(this.GridData.EdgeTagNames[et].Contains("navierslip_linear"))
                                        NSc[i] = true;
                                }
                            }
                        }
                        NScm = new CellMask(this.GridData, NSc);
                        CellMask bndNScm = NScm.AllNeighbourCells();
                        int bndLvl = 2;
                        for(int lvl = 1; lvl < bndLvl; lvl++) {
                            NScm = NScm.Union(bndNScm);
                            bndNScm = NScm.AllNeighbourCells();
                            NSbuffer = NScm.Union(bndNScm);
                        }
                        NSbuffer = NSbuffer.AllNeighbourCells();
                    }


                    //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 1), 2);


                bool AnyChange = GridRefinementController.ComputeGridChange((BoSSS.Foundation.Grid.Classic.GridData) this.GridData, BlockedCells, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
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

                    if(AnyChange) {

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



        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
            if(this.Control.solveCoupledHeatEquation)
                m_BDF_coupledTimestepper.DataBackupBeforeBalancing(L);
        }




        #endregion


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
        /// init routine for the specialized Fourier level-set
        /// </summary>
        private void InitFourier() {
            if(this.Control.FourierLevSetControl == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of FourierLevSetControl!");

            Fourier_LevSet = FourierLevelSetFactory.Build(this.Control.FourierLevSetControl);
            if(this.Control.EnforceLevelSetConservation) {
                throw new NotSupportedException("mass conservation correction currently not supported");
            }
            Fourier_LevSet.ProjectToDGLevelSet(this.DGLevSet.Current, this.LsTrk);

            if(base.MPIRank == 0 && this.CurrentSessionInfo.ID != Guid.Empty) {
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
        private void InitLevelSet() {
            using(new FuncTrace()) {

                // check level-set
                if(this.LevSet.L2Norm() == 0) {
                    throw new NotSupportedException("Level set is not initialized - norm is 0.0 - ALL cells will be cut, no gradient can be defined!");
                }

                // tracker needs to be updated to get access to the cut-cell mask
                this.LsTrk.UpdateTracker();

                // ==============================
                // level-set initialization
                // ==============================

                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 0 }), 3);

                #region Initialize Level Set Evolution Algorithm
                switch(this.Control.Option_LevelSetEvolution) {
                    case LevelSetEvolution.Fourier:
                        InitFourier();
                        break;
                    case LevelSetEvolution.None:
                        if(this.Control.AdvancedDiscretizationOptions.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Fourier) {
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
                            switch(this.Control.Timestepper_Scheme) {
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
                                        if(this.Control.Timestepper_Scheme.ToString().StartsWith("BDF")) {
                                            //do not instantiate rksch, use bdf instead
                                            bdfOrder = Convert.ToInt32(this.Control.Timestepper_Scheme.ToString().Substring(3));
                                            break;
                                        } else
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

                if(this.Control.Option_LevelSetEvolution != LevelSetEvolution.ExtensionVelocity)
                    ContinuityEnforcer.SetFarField(this.DGLevSet.Current, Near1, PosFF);

                ContinuityEnforcer.MakeContinuous(this.DGLevSet.Current, this.LevSet, Near, PosFF);

                //PlotCurrentState(0.0, new TimestepNumber(new int[] { 0, 2 }), 3);

                this.LsTrk.UpdateTracker();

            }

        }

        /// <summary>
        /// 
        /// </summary>
        public void PushLevelSetAndRelatedStuff() {

            if(this.Control.Option_LevelSetEvolution == LevelSetEvolution.Fourier) {
                Fourier_Timestepper.updateFourierLevSet();
            }

            this.ExtensionVelocity.IncreaseHistoryLength(1);
            this.ExtensionVelocity.Push();

            this.DGLevSet.IncreaseHistoryLength(1);
            this.DGLevSet.Push();
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
        /// only the velocity components(supposed to be at the beginning) are used.
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


                // =============================================
                // compute interface velocity due to evaporation
                // =============================================

                #region Compute evaporative velocity

                //SinglePhaseField LevSetSrc = new SinglePhaseField(meanVelocity[0].Basis, "LevelSetSource");

                if(this.Control.solveCoupledHeatEquation &&
                    this.Control.ThermalParameters.hVap_A != 0.0 && this.Control.ThermalParameters.hVap_B != 0.0) {

                    SinglePhaseField[] evapVelocity = new SinglePhaseField[D];
                    BitArray EvapMicroRegion = new BitArray(this.LsTrk.GridDat.Cells.Count);  //this.LsTrk.GridDat.GetBoundaryCells().GetBitMask();


                    double kA = this.Control.ThermalParameters.k_A;
                    double kB = this.Control.ThermalParameters.k_B;

                    for(int d = 0; d < D; d++) {
                        evapVelocity[d] = new SinglePhaseField(meanVelocity[0].Basis, "evapVelocity_d" + d);

                        evapVelocity[d].ProjectField(1.0,
                           delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                               int K = result.GetLength(1); // No nof Nodes

                               MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                               MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                               this.Temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                               this.Temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

                               MultidimensionalArray TempA_Res = MultidimensionalArray.Create(Len, K);
                               MultidimensionalArray TempB_Res = MultidimensionalArray.Create(Len, K);
                               MultidimensionalArray Curv_Res = MultidimensionalArray.Create(Len, K);
                               MultidimensionalArray Pdisp_Res = MultidimensionalArray.Create(Len, K);

                               this.Temperature.GetSpeciesShadowField("A").Evaluate(j0, Len, NS, TempA_Res);
                               this.Temperature.GetSpeciesShadowField("B").Evaluate(j0, Len, NS, TempB_Res);
                               this.Curvature.Evaluate(j0, Len, NS, Curv_Res);
                               this.DisjoiningPressure.Evaluate(j0, Len, NS, Pdisp_Res);

                               var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                               for(int j = 0; j < Len; j++) {
                                   for(int k = 0; k < K; k++) {

                                       double hVap = 0.0;
                                       double rho_l = 0.0;
                                       double rho_v = 0.0;
                                       double qEvap = 0.0;
                                       if(EvapMicroRegion[j]) {
                                           // micro region
                                           double Tsat = this.Control.ThermalParameters.T_sat;
                                           double pc = this.Control.ThermalParameters.pc;
                                           double pc0 = (pc < 0.0) ? this.Control.PhysicalParameters.Sigma * Curv_Res[j, k] + Pdisp_Res[j, k] : pc;
                                           double f = this.Control.ThermalParameters.fc;
                                           double R = this.Control.ThermalParameters.Rc;
                                           if(this.Control.ThermalParameters.hVap_A > 0) {
                                               hVap = this.Control.ThermalParameters.hVap_A;
                                               rho_l = this.Control.PhysicalParameters.rho_A;
                                               rho_v = this.Control.PhysicalParameters.rho_B;
                                               double TintMin = Tsat * (1 + (pc0 / (hVap * rho_l)));
                                               double Rint = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rho_v * hVap.Pow2());
                                               if(TempA_Res[j, k] > TintMin)
                                                   qEvap = -(TempA_Res[j, k] - TintMin) / Rint;
                                           } else {
                                               hVap = -this.Control.ThermalParameters.hVap_A;
                                               rho_l = this.Control.PhysicalParameters.rho_B;
                                               rho_v = this.Control.PhysicalParameters.rho_A;
                                               double TintMin = Tsat * (1 + (pc0 / (hVap * rho_l)));
                                               double Rint = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rho_v * hVap.Pow2());
                                               if(TempB_Res[j, k] > TintMin)
                                                   qEvap = (TempB_Res[j, k] - TintMin) / Rint;
                                           }

                                       } else {
                                           //macro region
                                           if(this.Control.ThermalParameters.hVap_A > 0) {
                                               hVap = this.Control.ThermalParameters.hVap_A;
                                               rho_l = this.Control.PhysicalParameters.rho_A;
                                               rho_v = this.Control.PhysicalParameters.rho_B;
                                               for(int dd = 0; dd < D; dd++)
                                                   qEvap += (kA * GradTempA_Res[j, k, dd] - kB * GradTempB_Res[j, k, dd]) * Normals[j, k, dd];
                                           } else {
                                               hVap = -this.Control.ThermalParameters.hVap_A;
                                               rho_l = this.Control.PhysicalParameters.rho_B;
                                               rho_v = this.Control.PhysicalParameters.rho_A;
                                               for(int dd = 0; dd < D; dd++)
                                                   qEvap += (kB * GradTempB_Res[j, k, dd] - kA * GradTempA_Res[j, k, dd]) * Normals[j, k, dd];
                                           }
                                       }


                                       double mEvap = -0.1; // qEvap / hVap; // mass flux
                                       //result[j, k] = mEvap * ((1 / rho_v) - (1 / rho_l)) * Normals[j, k, d];   //
                                       result[j, k] = mEvap * (1 / rho_v) * Normals[j, k, d];   //
                                       //result[j, k] = - Normals[j, k, d];   //
                                   }
                               }
                           }, new Foundation.Quadrature.CellQuadratureScheme(true, LsTrk.Regions.GetCutCellMask()));

                    }

                    SinglePhaseField[] Mevap = new SinglePhaseField[D];
                    for(int d = 0; d < D; d++) {
                        Mevap[d] = new SinglePhaseField(meanVelocity[0].Basis, "Mevap_d" + d);
                        double rho_v = 0.0;
                        if(this.Control.ThermalParameters.hVap_A > 0) {
                            rho_v = this.Control.PhysicalParameters.rho_B;
                        } else {
                            rho_v = this.Control.PhysicalParameters.rho_A;
                        }
                        Mevap[d].Acc(rho_v, evapVelocity[d]);
                    }


                    // evaporation for micro region 
                    #region micro evaporation

                    //double f = this.Control.ThermalParameters.fc;
                    //double Tsat = this.Control.ThermalParameters.T_sat;
                    //double R = this.Control.ThermalParameters.Rc;

                    //double rho_l = 0.0;
                    //double h_Vap = 0.0;
                    //double R_int = 0.0;
                    //double dir = 0.0;   // direction of volume flow
                    //DGField Temp_Vap = new SinglePhaseField(meanVelocity[0].Basis);
                    //if(this.Control.ThermalParameters.hVap_A > 0 && this.Control.ThermalParameters.hVap_B < 0) {
                    //    rho_l = this.Control.PhysicalParameters.rho_A;
                    //    h_Vap = this.Control.ThermalParameters.hVap_A;
                    //    Temp_Vap = Temperature.GetSpeciesShadowField("A");
                    //    R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (this.Control.PhysicalParameters.rho_B * h_Vap.Pow2());
                    //    dir = 1.0;
                    //} else if(this.Control.ThermalParameters.hVap_A < 0 && this.Control.ThermalParameters.hVap_B > 0) {
                    //    rho_l = this.Control.PhysicalParameters.rho_B;
                    //    h_Vap = this.Control.ThermalParameters.hVap_B;
                    //    Temp_Vap = Temperature.GetSpeciesShadowField("B");
                    //    R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (this.Control.PhysicalParameters.rho_A * h_Vap.Pow2());
                    //    dir = -1.0;
                    //}
                    //double p_c = this.Control.ThermalParameters.pc;

                    //LevSetSrc.ProjectField(1.0,
                    //delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    //    int K = result.GetLength(1); // No nof Nodes

                    //    //BitArray sbba = new BitArray(this.Grid.NoOfUpdateCells);
                    //    //for(int j = j0; j < j0 + Len; j++)
                    //    //    sbba[j] = true;
                    //    //SubGrid sbgrd = new SubGrid(new CellMask(this.GridData, sbba));
                    //    //var cp = new BoSSS.Solution.LevelSetTools.ClosestPointFinder(this.LsTrk, 0, sbgrd, new NodeSet[] { NS });

                    //    //MultidimensionalArray CurvIntRes = cp.EvaluateAtCp(this.Curvature);
                    //    //MultidimensionalArray pDisIntRes = cp.EvaluateAtCp(this.DisjoiningPressure);
                    //    //MultidimensionalArray TempIntRes = cp.EvaluateAtCp(Temp_Vap);

                    //    MultidimensionalArray CurvIntRes = MultidimensionalArray.Create(Len, K);
                    //    MultidimensionalArray pDisIntRes = MultidimensionalArray.Create(Len, K);
                    //    MultidimensionalArray TempIntRes = MultidimensionalArray.Create(Len, K);

                    //    this.Curvature.Evaluate(j0, Len, NS, CurvIntRes);
                    //    this.DisjoiningPressure.Evaluate(j0, Len, NS, pDisIntRes);
                    //    Temp_Vap.Evaluate(j0, Len, NS, TempIntRes);

                    //    for(int j = 0; j < Len; j++) {
                    //        for(int k = 0; k < K; k++) {

                    //            double pc0 = (p_c < 0.0) ? this.Control.PhysicalParameters.Sigma * CurvIntRes[j, k] + pDisIntRes[j, k] : p_c;      // augmented capillary pressure (without nonlinear evaporative masss part)

                    //            double T_intMin = Tsat * (1 + (pc0 / (rho_l * h_Vap)));

                    //            double qEvap = 0.0;
                    //            double T_int = TempIntRes[j, k];
                    //            if(T_int > T_intMin)
                    //                qEvap = dir * (T_int - T_intMin) / R_int;

                    //            result[j, k] = qEvap * (h_Vap / rho_l); // volume flux
                    //        }
                    //    }
                    //}, new Foundation.Quadrature.CellQuadratureScheme(true, LsTrk.Regions.GetNearFieldMask(1)));

                    #endregion


                    // check interface velocity
                    int p = evapVelocity[0].Basis.Degree;
                    SubGrid sgrd = LsTrk.Regions.GetCutCellSubgrid4LevSet(0);
                    NodeSet[] Nodes = LsTrk.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(p * 2).Nodes).ToArray();

                    var cp = new ClosestPointFinder(LsTrk, 0, sgrd, Nodes);

                    MultidimensionalArray[] VelocityEval = evapVelocity.Select(sf => cp.EvaluateAtCp(sf)).ToArray();

                    double nNodes = VelocityEval[0].Length;
                    double evapVelX = VelocityEval[0].Sum() / nNodes;
                    double evapVelY = VelocityEval[1].Sum() / nNodes;
                    Console.WriteLine("EvapVelocity: ({0},{1})", evapVelX, evapVelY);


                    // construct evolution velocity
                    for(int d = 0; d < D; d++) {
                        //SinglePhaseField FiltEvapVeloc = new SinglePhaseField(evapVelocity[d].Basis);
                        //FiltEvapVeloc.AccLaidBack(1.0, evapVelocity[d]);
                        //Filter(FiltEvapVeloc, 2, LsTrk.Regions.GetCutCellMask());
                        //evapVelocity[d].Clear();
                        //evapVelocity[d].Acc(1.0, FiltEvapVeloc);

                        meanVelocity[d].Clear();
                        if(this.Control.ThermalParameters.hVap_A > 0.0)
                            meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField("B"), this.LsTrk.Regions.GetCutCellMask());
                        else
                            meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField("A"), this.LsTrk.Regions.GetCutCellMask());

                        meanVelocity[d].Acc(1.0, evapVelocity[d]);
                    }

                    // plot
                    //Tecplot.PlotFields(Mevap.ToArray(), "Mevap" + hack_TimestepIndex, hack_Phystime, 2);
                    //Tecplot.PlotFields(evapVelocity.ToArray(), "EvapVelocity" + hack_TimestepIndex, hack_Phystime, 2);
                    //Tecplot.PlotFields(meanVelocity.ToArray(), "meanVelocity" + hack_TimestepIndex, hack_Phystime, 2);
                }

                #endregion

                // ===================================================================
                // backup interface properties (mass conservation, surface changerate)
                // ===================================================================

                #region backup interface props

                double oldSurfVolume = 0.0;
                double oldSurfLength = 0.0;
                double SurfChangerate = 0.0;
                if (this.Control.CheckInterfaceProps) {
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
                switch (this.Control.Option_LevelSetEvolution) {
                    case LevelSetEvolution.None:
                    throw new ArgumentException("illegal call");

                    case LevelSetEvolution.FastMarching: {

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

                            // Fast Marching: Specify the Domains first
                            // Perform Fast Marching only on the Far Field
                            if(this.Control.AdaptiveMeshRefinement) {
                                int NoCells = ((GridData)this.GridData).Cells.Count;
                                BitArray Refined = new BitArray(NoCells);
                                for(int j = 0; j < NoCells; j++) {
                                    if(((GridData)this.GridData).Cells.GetCell(j).RefinementLevel > 0)
                                        Refined[j] = true;
                                }
                                CellMask Accepted = new CellMask(this.GridData, Refined);
                                CellMask AcceptedNeigh = Accepted.AllNeighbourCells();

                                Accepted = Accepted.Union(AcceptedNeigh);
                                CellMask ActiveField = Accepted.Complement();
                                CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                                FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);

                            } else {
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
                if (underrelax < 1.0) {
                    this.DGLevSet.Current.Scale(underrelax);
                    this.DGLevSet.Current.Acc((1.0 - underrelax), DGLevSet_oldIter);
                }

                //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 1 }), 2);


                #endregion


                // ======================
                // postprocessing  
                // =======================

                if(this.Control.ReInitPeriod > 0 && hack_TimestepIndex % this.Control.ReInitPeriod == 0) {
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

                if(this.Control.Option_LevelSetEvolution == LevelSetEvolution.FastMarching) {
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
                if(this.Control.AdaptiveMeshRefinement && this.Control.Option_LevelSetEvolution == LevelSetEvolution.FastMarching) {
                    Near1 = LsTrk.Regions.GetNearMask4LevSet(0, 1);
                    PosFF = LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
                    ContinuityEnforcer.SetFarField(this.DGLevSet.Current, Near1, PosFF);
                    ContinuityEnforcer.SetFarField(this.LevSet, Near1, PosFF);
                }


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


        private void Filter(SinglePhaseField FiltrdField, int NoOfSweeps, CellMask CC) {

            Basis patchRecoveryBasis = FiltrdField.Basis;

            L2PatchRecovery l2pr = new L2PatchRecovery(patchRecoveryBasis, patchRecoveryBasis, CC, true);

            SinglePhaseField F_org = FiltrdField.CloneAs();

            for(int pass = 0; pass < NoOfSweeps; pass++) {
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
        private ConventionalDGField[] GetMeanVelocityFromXDGField(DGField[] EvoVelocity) {
            int D = EvoVelocity.Length;
            ConventionalDGField[] meanVelocity;

            Debug.Assert(this.XDGvelocity != null);
            Debug.Assert(this.DGvelocity == null);  // make sure that we run velocity enrichment

            meanVelocity = new ConventionalDGField[D];

            double rho_A = this.Control.PhysicalParameters.rho_A, rho_B = this.Control.PhysicalParameters.rho_B;
            double mu_A = this.Control.PhysicalParameters.mu_A, mu_B = this.Control.PhysicalParameters.mu_B;
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
                    double muSpc;
                    switch (spc) {
                        case "A": rhoSpc = rho_A; muSpc = mu_A; break;
                        case "B": rhoSpc = rho_B; muSpc = mu_B; break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }

                    double scale = 1.0;
                    switch(this.Control.InterAverage) {
                        case XNSE_Control.InterfaceAveraging.mean: {
                                scale = 0.5;
                                break;
                            }
                        case XNSE_Control.InterfaceAveraging.density: {
                                scale = rhoSpc / (rho_A + rho_B);
                                break;
                            }
                        case XNSE_Control.InterfaceAveraging.viscosity: {
                                scale = muSpc / (mu_A + mu_B);
                                break;
                            }
                    }
                     
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

        #endregion


        // ===================================================
        // pre-/postprocessing, compute properteis and logging
        // ===================================================

        private void Preprocessing(int TimestepInt, double phystime, double dt, TimestepNumber TimestepNo) {

            if(this.Control.CheckInterfaceProps) {
                double CL_length = this.GetContactLineLength();
                Console.WriteLine("contact line length = {0}", CL_length);

                double[] props = this.ComputeSphericalPorperties();
                Console.WriteLine("volume = {0}", props[0]);
                Console.WriteLine("surface = {0}", props[1]);
            }

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

            if(this.Control.CheckJumpConditions) {

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

                for(int d = 0; d < this.Grid.SpatialDimension; d++) {
                    this.MomentumBalanceAtInterface[d].Clear();
                    XNSEUtils.ProjectMomentumBalanceNorm(this.MomentumBalanceAtInterface[d], 1.0, this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                        this.Control.PhysicalParameters, this.Control.AdvancedDiscretizationOptions.SurfStressTensor, d, this.m_HMForder);
                }


                // energy balance
                //double energyBal_Norm = XNSEUtils.EnergyBalanceNormAtInterface(this.Pressure, this.XDGvelocity.Velocity, meanVelocity, this.Curvature,
                //    this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, this.Control.PhysicalParameters.Sigma, this.m_HMForder);

                //Console.WriteLine("energy balance norm = {0}", energyBal_Norm);

                //this.EnergyBalanceAtInterface.Clear();
                //XNSEUtils.ProjectEnergyBalanceNorm(this.EnergyBalanceAtInterface, 1.0, this.Pressure, this.XDGvelocity.Velocity, meanVelocity, this.Curvature,
                //    this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, this.Control.PhysicalParameters.Sigma, this.m_HMForder);

            }

            #endregion


            // ====================================================
            // Compute integral energy of the system
            // ====================================================

            #region integral energy computation 

            if (this.Control.ComputeEnergy) {

                // compute current energies
                double[] rhoS = new double[] { this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B };
                double currentKinEnergy = EnergyUtils.GetKineticEnergy(this.LsTrk, this.XDGvelocity.Velocity.ToArray(), rhoS, this.m_HMForder);
                double currentSurfEnergy = EnergyUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma, this.m_HMForder);

                // compute changerates (kinetic, surface)
                double CR_KinEnergy = 0.0;
                double CR_SurfEnergy = 0.0;
                if (this.Control.CompMode == AppControl._CompMode.Transient) {
                    double prevKinEnergy = EnergyUtils.GetKineticEnergy(this.LsTrk, this.prevVel, rhoS, this.m_HMForder, 0);
                    CR_KinEnergy = (currentKinEnergy - prevKinEnergy) / dt;

                    double prevSurfEnergy = EnergyUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma, this.m_HMForder, 0);
                    CR_SurfEnergy = (currentSurfEnergy - prevSurfEnergy) / dt;

                    Console.WriteLine("current kinetic energy = {0}; actual changerate = {1}", currentKinEnergy, CR_KinEnergy);
                    Console.WriteLine("current surface energy = {0}; actual changerate = {1}", currentSurfEnergy, CR_SurfEnergy);
                }

                // changerate of kinetic energy from discretization
                double[] muS = new double[] { this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B };
                double kineticDissipationBulk = EnergyUtils.GetKineticDissipation(this.LsTrk, this.XDGvelocity.Velocity.ToArray(), muS, this.m_HMForder);
                EnergyUtils.ProjectKineticDissipation(this.KineticDissipation, this.LsTrk, this.XDGvelocity.Velocity.ToArray(), muS, this.m_HMForder);

                // changerate of surface energy form discretization
                ConventionalDGField[] meanVelocity = XNSEUtils.GetMeanVelocity(this.XDGvelocity.Velocity, this.LsTrk,
                    this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B);
                double SurfDivergence = EnergyUtils.GetSurfaceChangerate(this.LsTrk, meanVelocity, this.m_HMForder);


                // logging
                this.EnergyLogger.TimeStep = TimestepInt;
                this.EnergyLogger.CustomValue(phystime + dt, "PhysicalTime");
                this.EnergyLogger.CustomValue(currentKinEnergy, "KineticEnergy");
                this.EnergyLogger.CustomValue(currentSurfEnergy, "SurfaceEnergy");
                this.EnergyLogger.CustomValue(CR_KinEnergy, "ChangerateKineticEnergy");
                this.EnergyLogger.CustomValue(CR_SurfEnergy, "ChangerateSurfaceEnergy");
                this.EnergyLogger.CustomValue(SurfDivergence, "SurfaceDivergence");
                this.EnergyLogger.CustomValue(kineticDissipationBulk, "KineticDissipationBulk");


                // surface viscosity parts
                if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor != SurfaceSressTensor.Isotropic) {

                    double shearViscEnergyCR = 0.0;
                    double dilViscEnergyCR = 0.0;

                    // surface shear viscosity energy
                    if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation
                        || this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        shearViscEnergyCR = EnergyUtils.GetInterfaceShearViscosityEnergyCR(this.LsTrk, meanVelocity, this.Control.PhysicalParameters.mu_I, this.m_HMForder);
                    }

                    // surface dilatational viscosity energy
                    if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation
                        || this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        dilViscEnergyCR = EnergyUtils.GetInterfaceDilatationalViscosityEnergyCR(this.LsTrk, meanVelocity, this.Control.PhysicalParameters.lambda_I, this.m_HMForder);
                    }


                    this.EnergyLogger.CustomValue(shearViscEnergyCR, "ShearViscosityDR");
                    this.EnergyLogger.CustomValue(dilViscEnergyCR, "DilatationalViscosityDR");

                    Console.WriteLine("current kinetic energy dissipation from discretization = {0}", kineticDissipationBulk + shearViscEnergyCR + dilViscEnergyCR);

                } else {

                    Console.WriteLine("current kinetic energy dissipation from discretization = {0}", kineticDissipationBulk);
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
            this.divVelocity.Clear();
            this.divVelocity.Divergence(1.0, this.XDGvelocity.Velocity);
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


            //double[] props = this.ComputeSphericalPorperties();
            //Console.WriteLine("volume = {0}", props[0]);
            //Console.WriteLine("surface = {0}", props[1]);

            //double CL_length = this.GetContactLineLength();
            //Console.WriteLine("contact line length = {0}", CL_length);

            //double CapHeight = GetCapillaryHeight();
            //Console.WriteLine("Capillary height = {0}", CapHeight);

            //ContinuityEnforcer = new ContinuityProjection(
            //        ContBasis: this.LevSet.Basis,
            //        DGBasis: this.DGLevSet.Current.Basis,
            //        gridData: GridData,
            //        Option: Control.LSContiProjectionMethod
            //        );


            // ====================================
            // IO related to Fourier level set
            // ====================================

            if (base.MPIRank == 0) {
                // save restart infos for FLS
                if(Log_FourierLS != null) {
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

            if(this.Control.TestMode == true) {
                LogQueryValue(phystime + dt);
            } else {
                if(Log != null && this.Control.LogValues != XNSE_Control.LoggingValues.None && base.MPIRank == 0 && (TimestepNo.MajorNumber % this.Control.LogPeriod == 0))
                    try {
                        WriteLogLine(TimestepNo, phystime + dt);
                    } catch(Exception e) {
                        Console.WriteLine("An error occured during WriteLogLine: '{0}'", e);
                    }

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

        

        #region property computation


        public double[] ComputeSphericalPorperties() {

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), this.m_HMForder, 1).XQuadSchemeHelper;

            // area/volume
            double volume = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        volume += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            // surface
            double surface = 0.0;
            //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            var surfElemVol = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                surfElemVol.Compile(LsTrk.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return new double[] { volume, surface };

        }


        public double GetContactLineLength() {

            double CL_length = 0.0;

            if(this.LsTrk.GridDat.SpatialDimension == 3) {

                var metrics = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder);

                XQuadSchemeHelper SchemeHelper = metrics.XQuadSchemeHelper;
                EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"));

                var QuadDom = SurfaceElement_Edge.Domain;
                var boundaryCutEdge = QuadDom.Intersect(this.GridData.GetBoundaryEdgeMask());

                var innerDom = QuadDom.Except(this.GridData.GetBoundaryEdgeMask());

                System.Collections.BitArray lowerBits = new System.Collections.BitArray(((GridData)this.GridData).Edges.Count);
                foreach(Chunk cnk in boundaryCutEdge) {
                    for(int iE = cnk.i0; iE < cnk.JE; iE++) {
                        if(((GridData)this.GridData).Edges.EdgeTags[iE] == 1) {
                            lowerBits[iE] = true;
                        }
                    }
                }
                EdgeMask lowerDom = new EdgeMask(this.GridData, lowerBits);

                EdgeMask dom = lowerDom;

                var factory = metrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, LsTrk.GridDat.Grid.RefElements[0]);
                SurfaceElement_Edge = new EdgeQuadratureScheme(factory, dom);

                EdgeQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    SurfaceElement_Edge.Compile(LsTrk.GridDat, this.m_HMForder),
                    delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {
                        EvalResult.SetAll(1.0);
                    },
                    delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < length; i++)
                            CL_length += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
            }

            return CL_length;

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
                    for(int i = 0; i < Length; i++)
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
                    for(int i = i0; i < i0 + Length; i++) {
                        LsTrk.GridDat.TransformLocal2Global(QR.Nodes, nodes_global, i);
                        EvalResult.AccSubArray(1.0, nodes_global, new int[] { i - i0, -1, -1 });
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++) {
                        for(int d = 0; d < D; d++) {
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
                    for(int d = 0; d < D; d++) {
                        this.CurrentVel[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++) {
                        for(int d = 0; d < D; d++) {
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
                    for(int i = 0; i < Length; i++)
                        perimtr_b += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            double circ = Math.PI * diamtr_c / perimtr_b;

            return new double[] { area, center[0, 0], center[0, 1], circ, VelocityAtCenter[0, 0], VelocityAtCenter[0, 1] };
        }


        public double[] ComputeBenchmarkQuantities_LineInterface() {

            // interface length
            double length = 0.0;
            length = XNSEUtils.GetInterfaceLength(LsTrk);

            // species area
            double area = 0.0;
            area = XNSEUtils.GetSpeciesArea(LsTrk, LsTrk.SpeciesIdS[0]);

            // interface mean angle

            return new double[] { length, area };
        }


        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionVelocity"/> and <see cref="XNSE_Control.ExactSolutionPressure"/>).
        /// </summary>
        internal double[] ComputeL2Error(double time) {
            int D = this.GridData.SpatialDimension;
            double[] Ret = new double[D + 1];

            if(this.Control.ExactSolutionVelocity == null && this.Control.ExactSolutionPressure == null)
                // nothing to do
                return Ret;

            int order = 0;
            if(LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;


            // Velocity error
            // ==============
            if(this.Control.ExactSolutionVelocity != null) {
                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[D];

                foreach(var spc in this.LsTrk.SpeciesNames) {
                    L2Error_Species.Add(spc, new double[D]);

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                    for(int d = 0; d < D; d++) {
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

                for(int d = 0; d < D; d++) {
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                    Ret[d] = L2Error[d];
                }

            }


            // pressure error
            // ==============
            if(this.Control.ExactSolutionPressure != null) {

                // pass 1: mean value of pressure difference
                double DiffInt = 0;
                foreach(var spc in this.LsTrk.SpeciesNames) {

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

                foreach(var spc in this.LsTrk.SpeciesNames) {

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


        #endregion


        #region logging

        /// <summary>
        /// saves the vector Guid for the sample points 
        /// </summary>
        TextWriter Log_FourierLS;


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

            if(this.Control.WriteInterfaceP) {
                LogInterfaceP = base.DatabaseDriver.FsDriver.GetNewLog("InterfaceP", sessionID);
                string header = String.Format("{0}\t{1}\t{2}", "#timestep", "#time", "interfacePoints");
                LogInterfaceP.WriteLine(header);
                LogInterfaceP.Flush();
            }

            switch(this.Control.LogValues) {
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
                        string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "#time", "area", "center of mass - x", "center of mass - y", "circularity", "rise velocity");
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
                case XNSE_Control.LoggingValues.CapillaryHeight: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("CapillaryHeight", sessionID);
                        string header = String.Format("{0}\t{1}\t{2}\t{3}", "#timestep", "#time", "capillary-height", "at-PositionX");
                        Log.WriteLine(header);
                        Log.Flush();

                        break;
                    }
                default:
                    throw new ArgumentException("No specified LogFormat");
            }
        }

        /// <summary>
        /// writes one line to the Log File
        /// </summary>
        public void WriteLogLine(TimestepNumber TimestepNo, double phystime) {

            if(this.Control.WriteInterfaceP) {
                double[] interfaceP;
                if(Fourier_LevSet != null) {
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

            switch(this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        Complex DFT_k;
                        int numP;
                        if(Fourier_LevSet != null) {
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

                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime, BmQ_RB[0], BmQ_RB[1], BmQ_RB[2], BmQ_RB[3], BmQ_RB[5]);
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
                        var boundaryEdge = ((GridData)this.GridData).GetBoundaryEdgeMask().GetBitMask();
                        var boundaryCutEdge = QuadDom.Intersect(new EdgeMask((GridData)this.GridData, boundaryEdge, MaskType.Geometrical));

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
                                for(int d = 0; d < D; d++) {
                                    EvalResult[0, 0, d] = Vnode_g[0, d];
                                }

                                // contact line velocity
                                MultidimensionalArray U_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                MultidimensionalArray U_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                for(int d = 0; d < D; d++) {
                                    (meanVelocity[d] as SinglePhaseField).EvaluateEdge(i0, length, QR.Nodes, U_IN.ExtractSubArrayShallow(-1, -1, d), U_OUT.ExtractSubArrayShallow(-1, -1, d));
                                }

                                for(int d = 0; d < D; d++) {
                                    EvalResult[0, 0, 2 + d] = U_IN[0, 0, d];
                                }

                                // contact angle
                                MultidimensionalArray normal_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                MultidimensionalArray normal_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                                for(int d = 0; d < D; d++) {
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
                                for(int i = 0; i < length; i++) {
                                    if(ResultsOfIntegration[i, 2 * D] != 0.0) {
                                        contactAngles.Add(Math.Abs(ResultsOfIntegration[i, 2 * D]));
                                        double[] cp = new double[D];
                                        double[] cpV = new double[D];
                                        for(int d = 0; d < D; d++) {
                                            cp[d] = ResultsOfIntegration[i, d];
                                            cpV[d] = ResultsOfIntegration[i, 2 + d];
                                        }
                                        contactPoints.Add(cp);
                                        contactVelocities.Add(cpV);
                                    }
                                }
                            }
                        ).Execute();


                        for(int p = 0; p < contactAngles.Count; p++) {
                            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime, contactPoints.ElementAt(p)[0], contactPoints.ElementAt(p)[1], contactVelocities.ElementAt(p)[0], contactVelocities.ElementAt(p)[1], contactAngles.ElementAt(p));
                            Log.WriteLine(line);
                        }
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.CapillaryHeight: {

                        MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);

                        double h_min = double.MaxValue, x_pos = 0.0;
                        for(int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                            if(InterfacePoints[i, 1] < h_min) {
                                h_min = InterfacePoints[i, 1];
                                x_pos = InterfacePoints[i, 0];
                            }
                        }

                        string line = String.Format("{0}\t{1}\t{2}\t{3}", TimestepNo, phystime, h_min, x_pos);
                        Log.WriteLine(line);
                        Log.Flush();

                        break;
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

            if(this.Control.WriteInterfaceP) {
                double[] interfaceP;
                if(Fourier_LevSet != null) {
                    interfaceP = Fourier_LevSet.current_interfaceP.To1DArray();
                } else {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                    interfaceP = interP.ResizeShallow(interP.Length).To1DArray();
                }

                base.QueryResultTable.LogValue("interfaceP", interfaceP);

            }

            switch(this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        double amplitude;
                        if(Fourier_LevSet != null) {
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

                        base.QueryResultTable.LogValue("area", BmQ_RB[0]);
                        base.QueryResultTable.LogValue("yCM", BmQ_RB[2]);
                        base.QueryResultTable.LogValue("circ", BmQ_RB[3]);
                        base.QueryResultTable.LogValue("riseV", BmQ_RB[5]);

                        return;
                    }
                case XNSE_Control.LoggingValues.LinelikeLS: {
                        break;
                    }
                case XNSE_Control.LoggingValues.CirclelikeLS: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        base.QueryResultTable.LogValue("area", BmQ_RB[0]);
                        base.QueryResultTable.LogValue("xM", BmQ_RB[1]);
                        base.QueryResultTable.LogValue("yM", BmQ_RB[2]);
                        base.QueryResultTable.LogValue("circ", BmQ_RB[3]);
                        base.QueryResultTable.LogValue("vM_x", BmQ_RB[4]);
                        base.QueryResultTable.LogValue("vM_y", BmQ_RB[5]);

                        break;
                    }
                default:
                    return;
            }

        }


        #endregion


        // =========================================================
        // related stuff for property tracking (e.g. kinetic energy)
        // =========================================================
        #region tracking

        ResidualLogger m_EnergyLogger;

        /// <summary>
        /// Logger for kinetic and surface energy.
        /// </summary>
        ResidualLogger EnergyLogger {
            get {
                if(!this.Control.ComputeEnergy)
                    return null;

                if(m_EnergyLogger == null) {
                    m_EnergyLogger = new ResidualLogger(base.MPIRank, base.DatabaseDriver, base.CurrentSessionInfo.ID);
                    m_EnergyLogger.WriteResidualsToConsole = false;
                    m_EnergyLogger.WriteResidualsToTextFile = true;
                    m_EnergyLogger.TextFileFileName = "Energy";
                }

                return m_EnergyLogger;
            }
        }


        SinglePhaseField MassBalanceAtInterface;

        VectorField<SinglePhaseField> MomentumBalanceAtInterface;

        SinglePhaseField EnergyBalanceAtInterface;

        /// <summary>
        /// kinetic energy derived via \f$ \rho \frac{vec{u} \cdot \vec{u}}{ 2 } \f$
        /// </summary>
        XDGField DerivedKineticEnergy;

        XDGField GeneratedKineticEnergy;

        /// <summary>
        /// kinetic energy computed via <see cref="KineticEnergyBalanceOperator"/>
        /// </summary>
        XDGField KineticEnergy;

        XDGField prevKineticEnergy;

        /// <summary>
        /// Residual of the kinetic energy balance
        /// </summary>
        XDGField ResidualKineticEnergy;

        /// <summary>
        /// source term for the kinetic energy
        /// </summary>
        XDGField KineticDissipation;

        /// <summary>
        /// spatial Operator for the kinetic energy balance
        /// </summary>
        XSpatialOperatorMk2 KineticEnergyBalanceOperator;

        IDictionary<SpeciesId, IEnumerable<double>> MassScaleForEnergy {
            get {
                double rho_A = this.Control.PhysicalParameters.rho_A,
                    rho_B = this.Control.PhysicalParameters.rho_B;

                double[] _rho_A = new double[1];
                _rho_A[0] = rho_A;
                double[] _rho_B = new double[1];
                _rho_B[0] = rho_B;

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _rho_A);
                R.Add(this.LsTrk.GetSpeciesId("B"), _rho_B);

                return R;
            }
        }

        MultigridOperator.ChangeOfBasisConfig[][] MultigridEnergyOperatorConfig {
            get {
                int pEnergy = this.KineticEnergy.Basis.Degree;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[1][];
                for(int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[1];

                    // configuration for Temperature
                    configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, pEnergy - iLevel),
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { 0 }
                    };
                }

                return configs;
            }
        }


        EnergyMultiphaseBoundaryCondMap m_energyBcMap;

        /// <summary>
        /// Boundary conditions.
        /// </summary>
        EnergyMultiphaseBoundaryCondMap energyBcMap {
            get {
                if(m_energyBcMap == null) {
                    m_energyBcMap = new EnergyMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
                }
                return m_energyBcMap;
            }
        }

        CoordinateVector m_CurrentEnergySolution;

        /// <summary>
        /// Current temperature;
        /// </summary>
        internal CoordinateVector CurrentEnergySolution {
            get {
                if(m_CurrentEnergySolution == null) {
                    m_CurrentEnergySolution = new CoordinateVector(this.KineticEnergy);
                }
                return m_CurrentEnergySolution;
            }
        }

        CoordinateVector m_CurrentEnergyResidual;

        /// <summary>
        /// Current residual for coupled heat equation.
        /// </summary>
        internal CoordinateVector CurrentEnergyResidual {
            get {
                if(m_CurrentEnergyResidual == null) {
                    m_CurrentEnergyResidual = new CoordinateVector(this.ResidualKineticEnergy);
                }
                return m_CurrentEnergyResidual;
            }
        }


        /// <summary>
        /// Implicit timestepping using Backward-Differentiation-Formulas (BDF),
        /// specialized for XDG applications.
        /// </summary>
        XdgBDFTimestepping m_BDF_energyTimestepper;


        public void generateKinEnergyOperator() {

            int degK = this.KineticEnergy.Basis.Degree;
            
            int D = this.GridData.SpatialDimension;

            string[] CodName = new string[] { "kinBalance" };
            string[] Params = ArrayTools.Cat(
                 VariableNames.VelocityVector(D),
                 (new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0, D),
                 VariableNames.VelocityX_GradientVector(),
                 VariableNames.VelocityY_GradientVector(),
                 VariableNames.Pressure,
                 (new string[] { "PressureGradX", "PressureGradY", "PressureGradZ" }).GetSubVector(0, D),
                 (new string[] { "GravityX", "GravityY", "GravityZ" }).GetSubVector(0, D),
                 (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D),
                 "Curvature");
            string[] DomName = new string[] { "KineticEnergy" };

            double rhoA = this.Control.PhysicalParameters.rho_A;
            double rhoB = this.Control.PhysicalParameters.rho_B;
            double muA = this.Control.PhysicalParameters.mu_A;
            double muB = this.Control.PhysicalParameters.mu_B;
            double sigma = this.Control.PhysicalParameters.Sigma;

            double LFFA = this.Control.AdvancedDiscretizationOptions.LFFA;
            double LFFB = this.Control.AdvancedDiscretizationOptions.LFFB;


            var dntParams = this.Control.AdvancedDiscretizationOptions;

            // create operator
            // ===============
            KineticEnergyBalanceOperator = new XSpatialOperatorMk2(DomName, Params, CodName, (A, B, C) => degK * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2), null);


            // build the operator
            // ==================
            {

                // convective part
                // ================
                {
                    if(this.Control.PhysicalParameters.IncludeConvection) {

                        var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];

                        // kinetic energy
                        var convK = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyConvectionInBulk(D, energyBcMap, rhoA, rhoB, LFFA, LFFB, LsTrk);
                        comps.Add(convK); // Bulk component


                        bool movingmesh;
                        switch(this.Control.Timestepper_LevelSetHandling) {
                            case LevelSetHandling.Coupled_Once:
                                movingmesh = true;
                                break;
                            case LevelSetHandling.LieSplitting:
                            case LevelSetHandling.StrangSplitting:
                            case LevelSetHandling.None:
                                movingmesh = false;
                                break;
                            case LevelSetHandling.Coupled_Iterative:
                            default:
                                throw new NotImplementedException();
                        }

                        //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyConvectionAtLevelSet(D, LsTrk, rhoA, rhoB, LFFA, LFFB, this.Control.PhysicalParameters.Material, energyBcMap, movingmesh));       // LevelSet component
                    }
                }

                // Laplace of kinetic energy
                // =========================
                {
                    var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];

                    double penalty = dntParams.PenaltySafety;

                    var Visc = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyLaplace(
                        dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                        energyBcMap, D, muA, muB);

                    comps.Add(Visc);

                    if(dntParams.UseGhostPenalties) {
                        var ViscPenalty = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyLaplace(penalty * 1.0, 0.0, energyBcMap, D, muA, muB);
                        Xheat_Operator.GhostEdgesOperator.EquationComponents[CodName[0]].Add(ViscPenalty);
                    }

                    // Level-Set operator:
                    //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergylaplceAtLevelSet(LsTrk, muA, muB, penalty * 1.0));
                }

                // Divergence of stress tensor
                // ===========================
                {
                    var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
                    comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.StressDivergence(D, energyBcMap, muA, muB));

                    // Level-Set operator:
                    //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.StressDivergenceAtLevelSet(LsTrk, muA, muB));
                }

                // surface energy (surface tension)
                // ================================
                {
                    //var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
                    //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.SurfaceEnergy(D, LsTrk, sigma));
                }

                // pressure term
                // =============
                {
                    var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
                    comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.DivergencePressureEnergy(D, energyBcMap));
                    //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.PressureConvectionInBulk(D, energyBcMap, LFFA, LFFB, LsTrk));
                    //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.PressureGradientConvection(D));

                    // Level-Set operator:
                    //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.DivergencePressureEnergyAtLevelSet(LsTrk));
                }

                // dissipation
                // ===========
                {
                    var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
                    comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.Dissipation(D, muA, muB));
                }

                // gravity (volume forces)
                // =======================
                {
                    var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
                    comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.PowerofGravity(D, rhoA, rhoB));
                }


                // finalize
                // ========

                KineticEnergyBalanceOperator.Commit();

            }

        }


        void DelComputeEnergyOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {

            int D = this.GridData.SpatialDimension;

            SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();

            // parameter assembly
            // ==================    

            // velocity
            var VelMap = new CoordinateMapping(this.XDGvelocity.Velocity.ToArray());
            DGField[] VelParam = VelMap.Fields.ToArray();

            // velocity mean
            VectorField<XDGField> VelMeanParam = new VectorField<XDGField>(D, new XDGBasis(LsTrk, 0), "VelMean_", XDGField.Factory);
            XheatUtils.ComputeAverageU(VelParam, VelMeanParam, m_HMForder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, m_HMForder, 1).XQuadSchemeHelper, this.LsTrk);

            // velocity gradient vectors
            VectorField<DGField> GradVelX = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradient", XDGField.Factory);
            for(int d = 0; d < D; d++) {
                foreach(var Spc in this.LsTrk.SpeciesIdS) { 
                    DGField f_Spc = ((VelParam[0] as XDGField).GetSpeciesShadowField(Spc));
                    (GradVelX[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
                }
            }
            GradVelX.ForEach(F => F.CheckForNanOrInf(true, true, true));

            VectorField<DGField> GradVelY = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradient", XDGField.Factory);
            for(int d = 0; d < D; d++) {
                foreach(var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((VelParam[1] as XDGField).GetSpeciesShadowField(Spc));
                    (GradVelY[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
                }
            }
            GradVelY.ForEach(F => F.CheckForNanOrInf(true, true, true));

            // pressure and gradient
            var PressMap = new CoordinateMapping(this.Pressure);
            DGField[] PressParam = PressMap.Fields.ToArray();

            VectorField<DGField> PressGrad = new VectorField<DGField>(D, PressParam[0].Basis, "PressureGrad", XDGField.Factory);
            for(int d = 0; d < D; d++) {
                foreach(var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((PressParam[0] as XDGField).GetSpeciesShadowField(Spc));
                    (PressGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
                }
            }
            PressGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

            // gravity
            var GravMap = new CoordinateMapping(this.XDGvelocity.Gravity.ToArray());
            DGField[] GravParam = GravMap.Fields.ToArray();

            // normals:
            SinglePhaseField[] Normals; // Normal vectors: length not normalized - will be normalized at each quad node within the flux functions.
            var LevelSetGradient = new VectorField<SinglePhaseField>(D, LevSet.Basis, SinglePhaseField.Factory);
            LevelSetGradient.Gradient(1.0, LevSet);
            Normals = LevelSetGradient.ToArray();

            // Curvature
            CurvatureAlgorithms.CurvatureDriver(
                SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                CurvatureAlgorithms.FilterConfiguration.Default,
                this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                this.m_HMForder, this.DGLevSet.Current);


            // concatenate everything
            var Params = ArrayTools.Cat<DGField>(
                VelParam,
                VelMeanParam,
                GradVelX,
                GradVelY,
                PressParam,
                PressGrad,
                GravMap,
                Normals, 
                this.Curvature);



            // assemble the matrix & affine vector
            // ===================================


            // compute matrix
            if(OpMtx != null) {

                var mtxBuilder = KineticEnergyBalanceOperator.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping, SpcToCompute);

                mtxBuilder.time = phystime;

                foreach(var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                }

                mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin eval = KineticEnergyBalanceOperator.GetEvaluatorEx(LsTrk,
                    CurrentState.ToArray(), Params, Mapping,
                    SpcToCompute);

                foreach(var kv in AgglomeratedCellLengthScales) {
                    eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                }

                eval.time = phystime;

                eval.Evaluate(1.0, 1.0, OpAffine);

            }


            //OpAffine.ScaleV(-1.0);

            //// mass matrix factory
            //MassFact = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), m_HMForder, 1).MassMatrixFactory;// new MassMatrixFactory(maxB, CurrentAgg);
            //var WholeMassMatrix = MassFact.GetMassMatrix(Mapping, MassScale); // mass matrix scaled with density rho

            //// add power of gravity forces
            //var WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity.Gravity.ToArray<DGField>()));
            //var WholeVelocity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray<DGField>()));
            //WholeMassMatrix.SpMV(1.0, WholeVelocity, 1.0, OpAffine);

            //// transform from RHS to Affine
            //OpAffine.ScaleV(-1.0);

        }

        /// <summary>
        /// dummy delegate for coupled operators
        /// </summary>
        /// <param name="CurrentState"></param>
        /// <param name="Phystime"></param>
        /// <param name="dt"></param>
        /// <param name="underrelax"></param>
        /// <param name="incremental"></param>
        /// <returns></returns>
        double DelUpdateLevelSet_EnergyOperator(DGField[] CurrentState, double Phystime, double dt, double underrelax, bool incremental) {
            // do nothing
            return 0.0;
        }


        /// <summary>
        /// The residual logger for this application.
        /// </summary>
        public ResidualLogger EnergyResLogger {
            get {
                return m_EnergyResLogger;
            }
        }

        ResidualLogger m_EnergyResLogger;


        #endregion


        // =====================================
        // related stuff for coupled heat solver
        // =====================================
        #region coupled heat solver

        /// <summary>
        /// prescribed volume flux for testing. If volume flux > 0, mass flux for domain A is > 0.
        /// </summary>
        //double m_prescribedVolumeFlux;


        /// <summary>
        /// prescribed disjoining pressure field for evaporation near wall 
        /// </summary>
        SinglePhaseField DisjoiningPressure;


        /// <summary>
        /// Temperature
        /// </summary>
        XDGField Temperature;

        /// <summary>
        /// Residual of the heat equation
        /// </summary>
        XDGField ResidualHeat;


        VectorField<XDGField> Heatflux;


        /// <summary>
        /// the spatial operator (heat equation)
        /// </summary>
        XSpatialOperatorMk2 Xheat_Operator;


        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        IDictionary<SpeciesId, IEnumerable<double>> HeatScale {
            get {
                double rho_A = this.Control.ThermalParameters.rho_A,
                    rho_B = this.Control.ThermalParameters.rho_B;

                double c_A = this.Control.ThermalParameters.c_A,
                    c_B = this.Control.ThermalParameters.c_B;

                double[] scale_A = new double[1];
                scale_A[0] = rho_A * c_A;
                double[] scale_B = new double[1];
                scale_B[0] = rho_B * c_B;

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), scale_A);
                R.Add(this.LsTrk.GetSpeciesId("B"), scale_B);

                return R;
            }
        }

        MultigridOperator.ChangeOfBasisConfig[][] MultigridCoupledOperatorConfig {
            get {
                int pTemp = this.Temperature.Basis.Degree;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[1][];
                for(int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[1];

                    // configuration for Temperature
                    configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, pTemp - iLevel),
                        mode = this.Control.TemperatureBlockPrecondMode,
                        VarIndex = new int[] { 0 }
                    };
                }


                return configs;
            }
        }


        ThermalMultiphaseBoundaryCondMap m_coupledBcMap;

        /// <summary>
        /// Boundary conditions.
        /// </summary>
        ThermalMultiphaseBoundaryCondMap coupledBcMap {
            get {
                if(m_coupledBcMap == null) {
                    m_coupledBcMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
                }
                return m_coupledBcMap;
            }
        }

        CoordinateVector m_CurrentCoupledSolution;

        /// <summary>
        /// Current temperature;
        /// </summary>
        internal CoordinateVector CurrentCoupledSolution {
            get {
                if(m_CurrentCoupledSolution == null) {
                    m_CurrentCoupledSolution = new CoordinateVector(this.Temperature);
                }
                return m_CurrentCoupledSolution;
            }
        }

        CoordinateVector m_CurrentCoupledResidual;

        /// <summary>
        /// Current residual for coupled heat equation.
        /// </summary>
        internal CoordinateVector CurrentCoupledResidual {
            get {
                if(m_CurrentCoupledResidual == null) {
                    m_CurrentCoupledResidual = new CoordinateVector(this.ResidualHeat);
                }
                return m_CurrentCoupledResidual;
            }
        }


        /// <summary>
        /// Implicit timestepping using Backward-Differentiation-Formulas (BDF),
        /// specialized for XDG applications.
        /// </summary>
        XdgBDFTimestepping m_BDF_coupledTimestepper;


        public void generateCoupledOperator() {

            int degT = this.Temperature.Basis.Degree;

            int D = this.GridData.SpatialDimension;

            string[] CodName = new string[] { "heat" };
            string[] Params = ArrayTools.Cat(
                 VariableNames.VelocityVector(D),
                 (new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0, D),
                 (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D),
                 (new string[] { "GradTemp0_X", "GradTemp0_Y", "GradTemp0_Z" }.GetSubVector(0, D)),
                 "Temperature0", "Curvature", "DisjoiningPressure");
            string[] DomName = new string[] { VariableNames.Temperature };


            // create operator
            // ===============
            Xheat_Operator = new XSpatialOperatorMk2(DomName, Params, CodName, (A, B, C) => m_HMForder, this.LsTrk.SpeciesIdS.ToArray());


            // build the operator
            // ==================
            {

                // species bulk components
                for (int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
                    // heat equation
                    Solution.XheatCommon.XOperatorComponentsFactory.AddSpeciesHeatEq(Xheat_Operator,
                        CodName[0], D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], coupledBcMap, XOpConfig, LsTrk);
                }

                // interface components
                Solution.XheatCommon.XOperatorComponentsFactory.AddInterfaceHeatEq(Xheat_Operator,
                        CodName[0], D, coupledBcMap, XOpConfig, LsTrk);


                if (XOpConfig.isEvaporation)
                    XOperatorComponentsFactory.AddInterfaceHeatEq_withEvaporation(Xheat_Operator, CodName[0], D, XOpConfig, LsTrk);


                // finalize
                // ========

                Xheat_Operator.Commit();

            }

        }


        void DelComputeCoupledOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {

            int D = this.GridData.SpatialDimension;

            SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();

            // parameter assembly
            // ==================    

            // velocity
            var VelMap = new CoordinateMapping(this.XDGvelocity.Velocity.ToArray());
            DGField[] VelParam = VelMap.Fields.ToArray();

            // velocity mean
            VectorField<XDGField> VelMeanParam = new VectorField<XDGField>(D, new XDGBasis(LsTrk, 0), "VelMean_", XDGField.Factory);
            XheatUtils.ComputeAverageU(VelParam, VelMeanParam, m_HMForder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, m_HMForder, 1).XQuadSchemeHelper, this.LsTrk);

            // normals:
            SinglePhaseField[] Normals; // Normal vectors: length not normalized - will be normalized at each quad node within the flux functions.
            var LevelSetGradient = new VectorField<SinglePhaseField>(D, LevSet.Basis, SinglePhaseField.Factory);
            LevelSetGradient.Gradient(1.0, LevSet);
            Normals = LevelSetGradient.ToArray();

            // Temperature0
            var TempMap = new CoordinateMapping(this.Temperature);
            DGField[] TempParam = TempMap.Fields.ToArray();

            // Temperature gradient for evaporation
            VectorField<DGField> GradTempParam = new VectorField<DGField>(D, TempParam[0].Basis, XDGField.Factory);
            GradTempParam = new VectorField<DGField>(D, TempParam[0].Basis, "GradTemp0_", XDGField.Factory);
            XNSEUtils.ComputeGradientForParam(TempParam[0], GradTempParam, this.LsTrk);


            // concatenate everything
            var Params = ArrayTools.Cat<DGField>(
                VelParam,
                VelMeanParam,
                Normals,
                GradTempParam,
                TempParam,
                this.Curvature,
                this.DisjoiningPressure);



            BitArray EvapMicroRegion = this.LsTrk.GridDat.GetBoundaryCells().GetBitMask();
            EvapMicroRegion.SetAll(false);


            // assemble the matrix & affine vector
            // ===================================

            // compute matrix
            if (OpMtx != null) {

                var mtxBuilder = Xheat_Operator.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping, SpcToCompute);

                mtxBuilder.time = phystime;

               
                foreach(var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                    //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("SlipLengths", kv.Value);
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("EvapMicroRegion", EvapMicroRegion);
                }

                mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin eval = Xheat_Operator.GetEvaluatorEx(LsTrk,
                    CurrentState.ToArray(), Params, Mapping,
                    SpcToCompute);

                foreach(var kv in AgglomeratedCellLengthScales) {
                    eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                    //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("SlipLengths", kv.Value);
                    eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("EvapMicroRegion", EvapMicroRegion);
                }
                

                //if(Op.SurfaceElementOperator.TotalNoOfComponents > 0) {
                //    foreach(var kv in InterfaceLengths)
                //        eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                //}

                eval.time = phystime;

                eval.Evaluate(1.0, 1.0, OpAffine);

            }

        }

        /// <summary>
        /// dummy delegate for coupled operators
        /// </summary>
        /// <param name="CurrentState"></param>
        /// <param name="Phystime"></param>
        /// <param name="dt"></param>
        /// <param name="underrelax"></param>
        /// <param name="incremental"></param>
        /// <returns></returns>
        double DelUpdateLevelSet_CoupledOperator(DGField[] CurrentState, double Phystime, double dt, double underrelax, bool incremental) {
            // do nothing
            return 0.0;
        }


        /// <summary>
        /// The residual logger for this application.
        /// </summary>
        public ResidualLogger CouplededResLogger {
            get {
                return m_CoupledResLogger;
            }
        }

        ResidualLogger m_CoupledResLogger;


        /// <summary>
        /// 
        /// </summary>
        public void ComputeHeatflux() {
            using(FuncTrace ft = new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                for(int d = 0; d < D; d++) {

                    foreach(var Spc in LsTrk.SpeciesNames) { // loop over species...
                        // shadow fields
                        DGField Temp_Spc = (this.Temperature.GetSpeciesShadowField(Spc));

                        double kSpc = 0.0;
                        switch(Spc) {
                            case "A": kSpc = this.Control.ThermalParameters.k_A; break;
                            case "B": kSpc = this.Control.ThermalParameters.k_B; break;
                            default: throw new NotSupportedException("Unknown species name '" + Spc + "'");
                        }

                        (this.Heatflux[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(kSpc, Temp_Spc, d);
                    }
                }

                this.Heatflux.ForEach(F => F.CheckForNanOrInf(true, true, true));

            }
        }


        #endregion
    }
}
