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
using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution.RheologyCommon;
using BoSSS.Solution.Statistic;

namespace BoSSS.Application.XRheology_Solver {

    /// <summary>
    /// Solver for Incompressible Multiphase flows
    /// </summary>
    public class XRheology_SolverMain : BoSSS.Application.XNSE_Solver.XBase_Solver<XRheology_Control> {



        static void Main(string[] args) {

            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.TestFixtureSetUp();
            ////BoSSS.Application.XNSE_Solver.Tests.UnitTest.PolynomialTestForConvectionTest(3, 0, false);
            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.TestCapillaryWave();
            ////BoSSS.Application.XNSE_Solver.Tests.ElementalTestProgramm.LineMovementTest(LevelSetEvolution.ScalarConvection, LevelSetHandling.Coupled_Once, XNSE_Control.TimesteppingScheme.ImplicitEuler, 0.5);
            //Assert.IsFalse(true, "remove me");


            _Main(args, false, delegate () {
                var p = new XRheology_SolverMain();
                return p;
            });
        }

        // Instantiate Fields from Control File

        #region instantiation
#pragma warning disable 649
        /// <summary>
        /// Pressure
        /// </summary>
        XDGField Pressure;

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

        // Parameters: Velocity Gradient
        public VectorField<XDGField> VelocityXGradient;
        public VectorField<XDGField> VelocityYGradient;


        //Parameters: external analytical velocity
        XDGField U;
        XDGField V;



        /// <summary>
        /// Extra source (e.g. gravity)
        /// </summary>
        //[InstantiateFromControlFile(new string[] { VariableNames.GravityX, VariableNames.GravityY, VariableNames.GravityZ },
        //            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
        //            true, true,
        //            IOListOption.ControlFileDetermined)]
        //public VectorField<XDGField> Gravity;
        //XDGField GravityX;
        //XDGField GravityY;

        // Gravity source constitutive
        //[InstantiateFromControlFile("GravityXX", "StressXX", IOListOption.ControlFileDetermined)]
        //XDGField GravityXX;

        //[InstantiateFromControlFile("GravityXY", "StressXY", IOListOption.ControlFileDetermined)]
        //XDGField GravityXY;

        //[InstantiateFromControlFile("GravityYY", "StressYY", IOListOption.ControlFileDetermined)]
        //XDGField GravityYY;

        //Gravity source for divergence of u
        //[InstantiateFromControlFile("GravityDiv", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        //public XDGField GravityDiv;


        /// <summary>
        /// Artificial force term at the fluid interface, usually only to support manufactured solutions.
        /// </summary>
        [InstantiateFromControlFile(
                new string[] { VariableNames.SurfaceForceX, VariableNames.SurfaceForceY, VariableNames.SurfaceForceZ },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true,
                IOListOption.ControlFileDetermined)]
        VectorField<SinglePhaseField> SurfaceForce;

        /// <summary>
        /// Curvature; DG-polynomial degree should be 2 times the polynomial degree of <see cref="LevSet"/>.
        /// </summary>
        [InstantiateFromControlFile(VariableNames.Curvature, VariableNames.Curvature, IOListOption.ControlFileDetermined)]
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
        /// PDE based elliptic reInitialization by Thomas
        /// </summary>
        EllipticReInit ReInitPDE;

        /// <summary>
        /// Set true if Navier Stokes is solved, then the mean velocities as parameters for calculation of convective terms are needed
        /// </summary>
        protected bool U0MeanRequired {
            get {
                return (this.Control.PhysicalParameters.IncludeConvection);
            }
        }
#pragma warning restore 649

        protected override void CreateFields() {
            using (new FuncTrace()) {
                base.CreateFields();
                int D = this.GridData.SpatialDimension;

                if (D > 2)
                    throw new NotImplementedException("The viscoelastic solver is only implemented for 2D cases!");

                //PRESSURE FIELD
                this.Pressure = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree), VariableNames.Pressure);
                base.RegisterField(this.Pressure);
                // CONTI FIELD
                this.ResidualContinuity = new XDGField(this.Pressure.Basis, "ResidualConti");
                base.RegisterField(this.ResidualContinuity);

                // ALL VELOCITY RELATED FIELDS
                this.XDGvelocity = new VelocityRelatedVars<XDGField>();
                InitFromAttributes.CreateFieldsAuto(this.XDGvelocity, this.GridData, base.Control.FieldOptions, base.Control.CutCellQuadratureType, base.IOFields, base.m_RegisteredFields);

                //this.GravityX = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), "GravityX");
                //base.RegisterField(this.GravityX);
                //this.GravityY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), "GravityY");
                //base.RegisterField(this.GravityY);

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

                //this.GravityXX = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXX].Degree), "GravityXX");
                //base.RegisterField(this.GravityXX);
                //this.GravityXY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressXY].Degree), "GravityXY");
                //base.RegisterField(this.GravityXY);
                //this.GravityYY = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.StressYY].Degree), "GravityYY");
                //base.RegisterField(this.GravityYY);


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

                XDGBasis b = new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree);
                this.divVelocity = new XDGField(b, "DivergenceVelocity");
                base.RegisterField(this.divVelocity);

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

                if (this.Control.ComputeEnergy) {

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
                double reynolds_A = this.Control.PhysicalParameters.reynolds_A,
                    reynolds_B = this.Control.PhysicalParameters.reynolds_B;

                int D = this.GridData.SpatialDimension;

                double[] _reynolds_A = new double[D + 4];
                _reynolds_A.SetAll(reynolds_A); // mass matrix in momentum equation
                _reynolds_A[D] = 0; // no  mass matrix for continuity equation
                _reynolds_A[D + 1] = 0; // no  mass matrix for constitutive equation
                _reynolds_A[D + 2] = 0; // no  mass matrix for constitutive equation
                _reynolds_A[D + 3] = 0; // no  mass matrix for constitutive equation
                double[] _reynolds_B = new double[D + 4];
                _reynolds_B.SetAll(reynolds_B); // mass matrix in momentum equation
                _reynolds_B[D] = 0; // no  mass matrix for continuity equation
                _reynolds_B[D + 1] = 0; // no  mass matrix for constitutive equation
                _reynolds_B[D + 2] = 0; // no  mass matrix for constitutive equation
                _reynolds_B[D + 3] = 0; // no  mass matrix for constitutive equation


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

        //IncompressibleMultiphaseBoundaryCondMap m_BcMap;

        ///// <summary>
        ///// Boundary conditions.
        ///// </summary>
        //IncompressibleMultiphaseBoundaryCondMap BcMap {
        //    get {
        //        if (m_BcMap == null) {
        //            m_BcMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
        //        }
        //        return m_BcMap;
        //    }
        //}

        /// <summary>
        /// the spatial operator (momentum and continuity equation)
        /// </summary>
        XRheology_OperatorFactory XRheology_Operator;

        /// <summary>
        /// OperatorConfiguration for the <see cref="XNSE_Operator"/>
        /// </summary>
        XRheology_OperatorConfiguration XOpConfig;

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
        /// current Weissenberg number
        /// </summary>
        public double[] currentWeissenberg = new double[] { 0.0, 0.0 };
        bool ChangeMesh = true;

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

        CoordinateVector m_CurrentSolution;

        /// <summary>
        /// Current velocity and pressure;
        /// </summary>
        internal CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(this.CurrentVel, this.Pressure, this.StressXX, this.StressXY, this.StressYY));
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
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat<DGField>(XDGvelocity.ResidualMomentum, ResidualContinuity, this.ResidualStressXX, this.ResidualStressXY, this.ResidualStressYY));
                }
                return m_CurrentResidual;
            }
        }


        /// <summary>
        /// output of <see cref="AssembleMatrix"/>;
        /// </summary>
        MassMatrixFactory MassFact;

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



        /// <summary>
        /// Create XOperator and Timestepper
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            #region Checks
            // CreateEquationsAndSolvers might be called multiple times
            // exit if so, and no LoadBalancing
            if (XRheology_Operator != null && L == null)
                return;


            if (Control.TimesteppingMode == AppControl._TimesteppingMode.Steady) {
                if (Control.Timestepper_LevelSetHandling != LevelSetHandling.None)
                    throw new ApplicationException(string.Format("Illegal control file: for a steady computation ({0}), the level set handling must be {1}.", AppControl._TimesteppingMode.Steady, LevelSetHandling.None));
            }

            int degU = this.CurrentVel[0].Basis.Degree;
            int stressDegree = this.StressXX.Basis.Degree;

            #endregion

            #region Config and Generate XOperator


            //Quadrature Order
            //----------------

            m_HMForder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 4 : 3);


            // Create Spatial Operator
            // ======================= 

            XOpConfig = new XRheology_OperatorConfiguration(this.Control);

            XRheology_Operator = new XRheology_OperatorFactory(XOpConfig, this.LsTrk, this.m_HMForder, this.BcMap, stressDegree, degU);
            //XRheology_Operator = new XNSE_OperatorFactory(XOpConfig, this.LsTrk, this.m_HMForder, this.BcMap, degU);

            #endregion

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
                                throw new NotImplementedException();
                        }

                }


                if (rksch == null) {
                    m_BDF_Timestepper = new XdgBDFTimestepping(
                        this.CurrentSolution.Mapping.Fields,
                        this.CurrentResidual.Mapping.Fields,
                        LsTrk,
                        true,
                        DelComputeOperatorMatrix, null, DelUpdateLevelSet,
                        (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? bdfOrder : 1,
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
                    m_BDF_Timestepper.Timestepper_Init = (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? this.Control.Timestepper_BDFinit : TimeStepperInit.SingleInit;
                    m_BDF_Timestepper.incrementTimesteps = this.Control.incrementTimesteps;
                    m_BDF_Timestepper.PushLevelSet = this.PushLevelSetAndRelatedStuff;
                    m_BDF_Timestepper.IterUnderrelax = this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? this.Control.LSunderrelax : 1.0;

                    m_BDF_Timestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                    //m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;


                    // solver 
                    this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations; //m_BDF_Timestepper.config_NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;

                    if (m_BDF_Timestepper.XdgSolverFactory.IsNewtonGmresType) {
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

                    //Console.WriteLine("noofpartsperprocess = {0}", this.CurrentSolution.Count / 10000);               

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


                if (this.Control.AdaptiveMeshRefinement && hack_TimestepIndex == 0) {
                    base.SetInitial();
                    this.InitLevelSet();
                }


                m_BDF_Timestepper.DataRestoreAfterBalancing(L,
                    ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure, this.StressXX, this.StressXY, this.StressYY),
                    ArrayTools.Cat<DGField>(this.XDGvelocity.ResidualMomentum.ToArray(), this.ResidualContinuity, this.ResidualStressXX, ResidualStressXY, this.ResidualStressYY),
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
                this.m_HMForder, SurfaceForce, filtLevSetGradient, Curvature, currentWeissenberg);

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

                //foreach (int jCell in this.LsTrk.Regions.GetCutCellMask4LevSet(0).ItemEnum) {
                //    for (int d = 0; d < D; d++) {
                //        for (int n = 0; n < N; n++) {
                //            ((XDGField)(TmpRhs.Mapping.Fields[d])).GetSpeciesShadowField("A").Coordinates[jCell, n] = 0.5 * SurfaceForce[d].Coordinates[jCell, n];
                //            ((XDGField)(TmpRhs.Mapping.Fields[d])).GetSpeciesShadowField("B").Coordinates[jCell, n] = 0.5 * SurfaceForce[d].Coordinates[jCell, n];
                //        }
                //    }
                //}

                OpAffine.AccV(1.0, TmpRhs);
            }

            // so far, 'SaddlePointRHS' is on the left-hand-side, since it is the output of ComputeMatrix
            // multiply by -1 to make it RHS
            OpAffine.ScaleV(-1.0);

            //foreach (string spc in LsTrk.SpeciesNames) {
            //    GravityX.GetSpeciesShadowField(spc).ProjectField(Control.GravityX[spc].Convert_Xt2X(phystime));
            //    int[] MomEqIdx1 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 0);
            //    OpAffine.AccV(-1.0, this.XDGvelocity.Gravity.CoordinateVector, MomEqIdx1, default(int[]));

            //    GravityY.GetSpeciesShadowField(spc).ProjectField(Control.GravityY[spc].Convert_Xt2X(phystime));
            //    int[] MomEqIdx2 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 1);
            //    OpAffine.AccV(-1.0, this.XDGvelocity.Gravity.CoordinateVector, MomEqIdx2, default(int[]));

                //GravityXX.GetSpeciesShadowField(spc).ProjectField(Control.GravityXX[spc].Convert_Xt2X(phystime));
                //int[] ConstEqIdx1 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 3);
                //OpAffine.AccV(-1.0, this.GravityXX.CoordinateVector, ConstEqIdx1, default(int[]));

                //GravityXY.GetSpeciesShadowField(spc).ProjectField(Control.GravityXY[spc].Convert_Xt2X(phystime));
                //int[] ConstEqIdx2 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 3);
                //OpAffine.AccV(-1.0, this.GravityXY.CoordinateVector, ConstEqIdx2, default(int[]));

                //GravityYY.GetSpeciesShadowField(spc).ProjectField(Control.GravityYY[spc].Convert_Xt2X(phystime));
                //int[] ConstEqIdx3 = this.CurrentSolution.Mapping.GetSubvectorIndices(true, 3);
                //OpAffine.AccV(-1.0, this.GravityXX.CoordinateVector, ConstEqIdx3, default(int[]));
            //}

                // ============================
                // Generate MassMatrix
                // ============================

                // mass matrix factory
                MassFact = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), m_HMForder, 1).MassMatrixFactory;// new MassMatrixFactory(maxB, CurrentAgg);
            var WholeMassMatrix = MassFact.GetMassMatrix(Mapping, MassScale); // mass matrix scaled with density rho

            // For ResidualTest from Markus only
            //======================================
            //BlockMsrMatrix inverseMassMatrix = MassFact.GetMassMatrix(CurrentResidual.Mapping, true);
            //inverseMassMatrix.SpMV(1.0, OpAffine, 0.0, CurrentResidual);


            // ============================
            //  Add Gravity
            // ============================
            // Dimension: [ rho * G ] = mass / time^2 / len^2 == [ d/dt rho U ]
            var WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity.Gravity.ToArray<DGField>(), new XDGField(this.Pressure.Basis), new XDGField(this.StressXX.Basis), new XDGField(this.StressXY.Basis), new XDGField(this.StressYY.Basis)));
            WholeMassMatrix.SpMV(1.0, WholeGravity, 1.0, OpAffine);


            // ============================
            // Set Pressure Reference Point
            // ============================

            if (Control.NonLinearSolver.UsePresRefPoint == true) {
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
            }

            // transform from RHS to Affine
            OpAffine.ScaleV(-1.0);

            //OpMtx.SaveToTextFile("OpMatrix");
        }

        int hack_TimestepIndex;
        double hack_Phystime;



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

#if DEBUG
            //Debug/Test code for XDG database interaction

            if (tsi != null) {
                // checking some neccessary reference-equalities BEFORE serialisation
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                LevelSet.LevelSetInitializer lsi_1 = (LevelSet.LevelSetInitializer)(tsi.FieldInitializers.Single(fi => fi.Identification == this.LevSet.Identification));
                XDGField.FieldInitializer pri = (XDGField.FieldInitializer)(tsi.FieldInitializers.Single(fi => fi.Identification == this.Pressure.Identification));

                LevelSetTracker.LevelSetTrackerInitializer trki = ((XDGBasis.XDGBasisInitializer)(pri.BasisInfo)).TrackerInitializer;

                LevelSet.LevelSetInitializer lsi_2 = trki.LevelSets[0];

                Debug.Assert(object.ReferenceEquals(lsi_1, lsi_2));

                foreach (XDGField.FieldInitializer fi in tsi.FieldInitializers.Where(fii => fii is XDGField.XDGFieldInitializer)) {
                    LevelSetTracker.LevelSetTrackerInitializer trki_alt = ((XDGBasis.XDGBasisInitializer)(fi.BasisInfo)).TrackerInitializer;
                    Debug.Assert(object.ReferenceEquals(trki, trki_alt));
                }
            }


            if (tsi != null) {
                // checking some neccessary equalities AFTER serialisation
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++


                var tsi_alt = this.DatabaseDriver.LoadTimestepInfo(tsi.ID, base.CurrentSessionInfo, base.GetDatabase());

                Debug.Assert(!object.ReferenceEquals(tsi, tsi_alt));


                LevelSet.LevelSetInitializer lsi_1 = (LevelSet.LevelSetInitializer)(tsi_alt.FieldInitializers.Single(fi => fi.Identification == "Phi"));
                XDGField.FieldInitializer pri = (XDGField.FieldInitializer)(tsi_alt.FieldInitializers.Single(fi => fi.Identification == this.Pressure.Identification));

                LevelSetTracker.LevelSetTrackerInitializer trki = ((XDGBasis.XDGBasisInitializer)(pri.BasisInfo)).TrackerInitializer;

                LevelSet.LevelSetInitializer lsi_2 = trki.LevelSets[0];

                Debug.Assert(lsi_1.Equals(lsi_2));

                foreach (XDGField.FieldInitializer fi in tsi_alt.FieldInitializers.Where(fii => fii is XDGField.XDGFieldInitializer)) {
                    LevelSetTracker.LevelSetTrackerInitializer trki_alt = ((XDGBasis.XDGBasisInitializer)(fi.BasisInfo)).TrackerInitializer;
                    Debug.Assert(trki.Equals(trki_alt));
                }


                var Fields = this.DatabaseDriver.LoadFields(tsi_alt, this.GridData);
                LevelSet Rphi_1 = (LevelSet)(Fields.Single(f => f.Identification.Equals(this.LevSet.Identification)));

                XDGField Rpressure = (XDGField)(Fields.Single(f => f.Identification.Equals(this.Pressure.Identification)));

                LevelSetTracker Rtracker = Rpressure.Basis.Tracker;
                Debug.Assert(!object.ReferenceEquals(this.LsTrk, Rtracker));
                Debug.Assert(object.ReferenceEquals(Rtracker.LevelSets[0], Rphi_1));

                foreach (XDGField xf in Fields.Where(fii => fii is XDGField)) {
                    Debug.Assert(object.ReferenceEquals(xf.Basis.Tracker, Rtracker));
                }
            }

#endif

            return tsi;
        }



        /// <summary>
        /// Depending on settings <see cref="AppControl.TimesteppingMode"/>, computes either one timestep or a steady-state solution.
        /// </summary>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (var tr = new FuncTrace()) {

                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;
                base.ResLogger.TimeStep = TimestepInt;

                if (TimestepNo[0] > 1) {
                    this.Control.RaiseWeissenberg = false;
                }

                hack_TimestepIndex = TimestepInt;
                hack_Phystime = phystime;
                int NoIncrementTimestep;


                Preprocessing(TimestepInt, phystime, dt, TimestepNo);


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
                        StressXX.GetSpeciesShadowField(spc).ProjectField(Control.ExactSolutionStressXX[spc].Convert_Xt2X(phystime + dt));
                        StressXY.GetSpeciesShadowField(spc).ProjectField(Control.ExactSolutionStressXY[spc].Convert_Xt2X(phystime + dt));
                        StressYY.GetSpeciesShadowField(spc).ProjectField(Control.ExactSolutionStressYY[spc].Convert_Xt2X(phystime + dt));
                    }
                }


                // =====================================================
                // setup stationary 
                // =====================================================


                if (base.Control.TimesteppingMode == AppControl._TimesteppingMode.Steady) {
                    dt = 1.0e100;
                    Console.WriteLine("Steady-state solve ...", TimestepNo, dt);

                    if (this.Control.Option_LevelSetEvolution != LevelSetEvolution.None) {
                        throw new ApplicationException("For steady-state solutions, the only allowed level-set-evolution option is '" + LevelSetEvolution.None + "'.");
                    }



                    // =====================================================
                    // setup transient 
                    // =====================================================
                } else if (base.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {

                    // push stacks
                    // -----------

                    PushLevelSetAndRelatedStuff();


                    // backup old velocity for energy checks
                    // -------------------------------------
                    if (this.Control.ComputeEnergy && this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
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
                    throw new NotImplementedException("Option " + base.Control.TimesteppingMode + " not supported yet.");
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

                           // currentWeissenberg = new double[] { 0.0, 0.0 };

                            if (Control.PhysicalParameters.Weissenberg_a != 0.0 || Control.PhysicalParameters.Weissenberg_b != 0.0) {

                                if (Control.WeissenbergIncrement != 0.0) {
                                    NoIncrementTimestep = 1;
                                    if(Control.PhysicalParameters.Weissenberg_a > Control.PhysicalParameters.Weissenberg_b)
                                        NoIncrementTimestep = (int)(Control.PhysicalParameters.Weissenberg_a / Control.WeissenbergIncrement);
                                    else if(Control.PhysicalParameters.Weissenberg_b > Control.PhysicalParameters.Weissenberg_a)
                                        NoIncrementTimestep = (int)(Control.PhysicalParameters.Weissenberg_b / Control.WeissenbergIncrement);
                                    else if (Control.PhysicalParameters.Weissenberg_b == Control.PhysicalParameters.Weissenberg_a)
                                        NoIncrementTimestep = (int)(Control.PhysicalParameters.Weissenberg_a / Control.WeissenbergIncrement);
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

                                    // this evaluation must later out of this loop. now here for comparing results with  
                                    PlotCurrentState(phystime, new TimestepNumber(TimestepNo.MajorNumber, i));
                                    SaveToDatabase(new TimestepNumber(TimestepNo.MajorNumber, i), phystime);

                                }

                                //ChangeMesh = Control.AdaptiveMeshRefinement;
                                //while (ChangeMesh == true) {
                                //    this.MpiRedistributeAndMeshAdapt(TimestepNo.MajorNumber, phystime);
                                //    perssonsensor.Update(StressXX);
                                //    PlotCurrentState(phystime, TimestepNo);
                                //}

                                if (currentWeissenberg[0] < Control.PhysicalParameters.Weissenberg_a) {
                                    currentWeissenberg[0] = currentWeissenberg[0] + Control.WeissenbergIncrement;
                                    Console.WriteLine();
                                    Console.WriteLine("Raise Weissenberg number A to " + currentWeissenberg[0]);
                                    Console.WriteLine();
                                }

                                if (currentWeissenberg[1] < Control.PhysicalParameters.Weissenberg_b) {
                                    currentWeissenberg[1] = currentWeissenberg[1] + Control.WeissenbergIncrement;
                                    Console.WriteLine();
                                    Console.WriteLine("Raise Weissenberg number B to " + currentWeissenberg[1]);
                                    Console.WriteLine();
                                }

                            }
                        } else {
                            //current Weissenberg is set to the HIGHER value... DIRTY HACK AT THE MOMENT!

                                currentWeissenberg[0] = Control.PhysicalParameters.Weissenberg_a;
                                currentWeissenberg[1] = Control.PhysicalParameters.Weissenberg_b;


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

                                if (this.Control.OperatorMatrixAnalysis == true) {

                                    BlockMsrMatrix SaddlePointMatrix = new BlockMsrMatrix(this.CurrentSolution.Mapping);
                                    double[] AffineDummy = new double[this.CurrentSolution.Mapping.LocalLength];

                                    var agg = LsTrk.GetAgglomerator(LsTrk.SpeciesIdS.ToArray(), m_HMForder, this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold);

                                    DelComputeOperatorMatrix(SaddlePointMatrix, AffineDummy, this.CurrentSolution.Mapping,
                                    this.CurrentSolution.Mapping.Fields.ToArray(), agg.CellLengthScales, 0.0);

                                    AggregationGridBasis[][] MgBasis = AggregationGridBasis.CreateSequence(this.MultigridSequence, this.CurrentSolution.Mapping.BasisS);
                                    //todo: AsyncCallback update
                                    MgBasis.UpdateXdgAggregationBasis(agg);
                                    MultigridOperator mgOp = new MultigridOperator(MgBasis, CurrentSolution.Mapping,
                                        SaddlePointMatrix, this.MassFact.GetMassMatrix(CurrentSolution.Mapping, false),
                                        this.MultigridOperatorConfig);

                                    MsrMatrix FullMatrix = mgOp.OperatorMatrix.ToMsrMatrix();

                                    MsrMatrix DiffMatrix;
                                    {
                                        int[] VelVarIdx = new int[] { 3, 4, 5 };

                                        int[] USubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                                        int[] USubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                                        int L = USubMatrixIdx_Row.Length;

                                        DiffMatrix = new MsrMatrix(L, L, 1, 1);
                                        FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(int[]), USubMatrixIdx_Col, default(int[]));
                                    }

                                    MultidimensionalArray ret = MultidimensionalArray.Create(1, 2);
                                    Console.WriteLine("Calling MATLAB/Octave...");
                                    using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                                        bmc.PutSparseMatrix(FullMatrix, "FullMatrix");
                                        bmc.PutSparseMatrix(DiffMatrix, "DiffMatrix");
                                        bmc.Cmd("DiffMatrix = 0.5*(DiffMatrix + DiffMatrix');");
                                        bmc.Cmd("condNoDiffMatrix = condest(DiffMatrix);");
                                        bmc.Cmd("condNoFullMatrix = condest(FullMatrix);");
                                        //bmc.Cmd("eigiMaxi = eigs(DiffMatrix,1,'lm')");
                                        //bmc.Cmd("eigiMini = eigs(DiffMatrix,1,'sm')");
                                        //bmc.Cmd("lasterr");
                                        //bmc.Cmd("[V,r]=chol(DiffMatrix);");
                                        bmc.Cmd("ret = [condNoFullMatrix, condNoDiffMatrix]");
                                        bmc.GetMatrix(ret, "ret");

                                        bmc.Execute(false);
                                    }

                                    double condNoFullMatrix = ret[0, 0];
                                    double condNoDiffMatrix = ret[0, 1];
                                    //double eigiMaxi = ret[0, 2];
                                    //double eigiMini = ret[0, 3];
                                    //posDef = ret[0, 4] == 0;

                                    //Console.WriteLine("Eigenvalue range of diffusion matrix: {0} to {1}", eigiMini, eigiMaxi);

                                    Console.WriteLine("Condition number diffusion operator: {0:0.####E-00}", condNoDiffMatrix);
                                    Console.WriteLine("Condition number full operator: {0:0.####E-00}", condNoFullMatrix);
                                    base.QueryHandler.ValueQuery("condFull", condNoFullMatrix, true);
                                    base.QueryHandler.ValueQuery("condDiff", condNoDiffMatrix, true);

                                    //OpAnalysisBase myAnalysis = new OpAnalysisBase(DelComputeOperatorMatrix, CurrentSolution.Mapping, CurrentSolution.Mapping.Fields.ToArray(), agg.CellLengthScales, phystime);
                                    //myAnalysis.VarGroup = new int[] { 3, 4, 5};
                                    ////myAnalysis.Analyse();
                                    //double[] condest = myAnalysis.CondNum();
                                    //Console.WriteLine("Condition number full system, full matrix: " + condest[0] + "full system inner matrix (excl. BC): " + condest[1]);

                                }


                                m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);

                                //===============================================================================

                                #region Write residuals to text file
                                if (Control.InterfaceTest == true) {
                                    // Sample points

                                    int noOfPoints = 1000;

                                    double[] nodes = GenericBlas.Linspace(-2, 2, noOfPoints);

                                    MultidimensionalArray points = MultidimensionalArray.Create(noOfPoints, 2);

                                    for (int i = 0; i < noOfPoints; i++) {

                                        points[i, 0] = nodes[i];

                                        points[i, 1] = 0.5;

                                    }



                                    // FieldEvaluation
                                    MultidimensionalArray results = MultidimensionalArray.Create(noOfPoints, CurrentResidual.Mapping.Count);

                                    for (int i = 0; i < CurrentResidual.Length; i++) {

                                        FieldEvaluation fieldEvaluator = new FieldEvaluation((GridData)this.GridData);

                                        fieldEvaluator.Evaluate(1.0, CurrentResidual.Mapping, points, 0.0, results);

                                    }



                                    // StreamWriter

                                    using (System.IO.StreamWriter sw = new System.IO.StreamWriter(String.Format("Residuals{0}.txt", dt))) {

                                        //Console.WriteLine("x \t y \t result");

                                        sw.WriteLine("x \t y \t momX \t momY \t conti \t constXX \t constXY \t constYY");

                                        string resultLine;

                                        for (int i = 0; i < noOfPoints; i++) {

                                            resultLine = points[i, 0] + "\t" + points[i, 1] + "\t" + results[i, 0] + "\t" + results[i, 1] + "\t" + results[i, 2] + "\t" + results[i, 3] + "\t" + results[i, 4] + "\t" + results[i, 5] + "\t";

                                            //Console.WriteLine(resultLine);

                                            sw.WriteLine(resultLine);

                                        }

                                        sw.Flush();

                                    }
                                }
                                #endregion
                                //=============================================================================================
                            }
                        }
                    }


                    if (this.Control.ComputeEnergy && m_BDF_energyTimestepper != null) {

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
        }


        //UnsetteledCoordinateMapping SaddlePointProblemMapping {
        //    get {
        //        return this.CurrentSolution.Mapping;
        //    }
        //}

        //int RepairZeroRows(MsrMatrix Mtx) {
        //    int NoOfZeroRows = 0;
        //    for (int iRow = Mtx.RowPartitioning.i0; iRow < Mtx.RowPartitioning.iE; iRow++) {
        //        if (Mtx.GetNoOfNonZerosPerRow(iRow) == 0) {
        //            Mtx[iRow, iRow] = +1.0;
        //            NoOfZeroRows++;
        //        }
        //    }
        //    return NoOfZeroRows;
        //}


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {
            Tecplot.PlotFields(base.m_RegisteredFields, "XRheology_Solver" + timestepNo, physTime, superSampling);
            //Tecplot.PlotFields(new DGField[] { this.LevSet }, "grid" + timestepNo, physTime, 0);
        }


        protected void PlotOnIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, iterIndex }), 2);
        }


        protected override void SetInitial() {
            base.SetInitial();

            this.InitLevelSet();

            this.CreateEquationsAndSolvers(null);


            m_BDF_Timestepper.SingleInit();

            int D = LsTrk.GridDat.SpatialDimension;

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

            After_SetInitialOrLoadRestart(0.0, 0);

        }

        /// <summary>
        /// delegate for the initialization of previous timesteps from an analytic solution
        /// </summary>
        /// <param name="TimestepIndex"></param>
        /// <param name="Time"></param>
        /// <param name="St"></param>
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


        private void After_SetInitialOrLoadRestart(double PhysTime, int TimestepNo) {

            // =============================================
            // LogFile initialization
            // =============================================  

            if (this.Control.TestMode == true) {
                LogQueryValue(PhysTime);
            } else {
                if (this.Control.LogValues != XRheology_Control.LoggingValues.None && this.CurrentSessionInfo.ID != Guid.Empty && base.MPIRank == 0) {
                    InitLogFile(this.CurrentSessionInfo.ID);
                    WriteLogLine(TimestepNo, PhysTime);
                }
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
        /// overriding the method to implement any user-specific tasks which
        /// should be carried out after a restart file has been loaded (e.g.,
        /// setting the correct time for a time-stepper)
        /// </summary>
        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            VelocityXGradient = new VectorField<XDGField>(this.GridData.SpatialDimension, this.CurrentVel[0].Basis, "VelocityX_Gradient", XDGField.Factory);
            VelocityYGradient = new VectorField<XDGField>(this.GridData.SpatialDimension, this.CurrentVel[1].Basis, "VelocityY_Gradient", XDGField.Factory);

            //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 20 }), 2);


            //if (this.Control.ClearVelocitiesOnRestart) {
            //    Console.WriteLine("clearing all velocities");
            //    this.XDGvelocity.Velocity.Clear();
            //}


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
                int pStr = this.StressXX.Basis.Degree;
                int D = this.GridData.SpatialDimension;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 4];

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
                    //configurations for stresses
                    for (int d = 3; d < 6; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            Degree = Math.Max(1, pStr - iLevel),
                            mode = this.Control.StressBlockPrecondMode,
                            VarIndex = new int[] { d }
                        };
                    }
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
                        case XRheology_Control.RefinementStrategy.CurvatureRefined: {
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
                        case XRheology_Control.RefinementStrategy.ContactLineRefined: {
                                CellMask BCells = ((GridData)this.GridData).BoundaryCells.VolumeMask;
                                if (ccm.Contains(j) && BCells.Contains(j) && CurrentLevel < this.Control.RefinementLevel) {
                                    DesiredLevel_j++;
                                } else if (!BCells.Contains(j)) { // && CurrentLevel == this.Control.RefinementLevel + 1) {
                                    DesiredLevel_j--;
                                }
                                break;
                            }
                        case XRheology_Control.RefinementStrategy.constantInterface:
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
                    if (this.Control.RefineStrategy == XRheology_Control.RefinementStrategy.CurvatureRefined) {
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


        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
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


        private void EnforceVolumeConservation() {
            double spcArea = XNSEUtils.GetSpeciesArea(LsTrk, LsTrk.SpeciesIdS[0]);
            Console.WriteLine("area = {0}", spcArea);
            double InterLength = XNSEUtils.GetInterfaceLength(LsTrk);

            //double cmc = (consvRefArea - spcArea) / InterLength;
            //Console.WriteLine("add constant: {0}", -cmc);
            //this.DGLevSet.Current.AccConstant(-cmc);
            //this.LevSet.AccConstant(-cmc);
        }


        //private void Filter(SinglePhaseField FiltrdField, int NoOfSweeps, CellMask CC) {

        //    Basis patchRecoveryBasis = FiltrdField.Basis;

        //    L2PatchRecovery l2pr = new L2PatchRecovery(patchRecoveryBasis, patchRecoveryBasis, CC, true);

        //    SinglePhaseField F_org = FiltrdField.CloneAs();

        //    for (int pass = 0; pass < NoOfSweeps; pass++) {
        //        F_org.Clear();
        //        F_org.Acc(1.0, FiltrdField);
        //        FiltrdField.Clear();
        //        l2pr.Perform(FiltrdField, F_org);
        //    }
        //}


        #endregion


        // ===================================================
        // pre-/postprocessing, compute properteis and logging
        // ===================================================

        private void Preprocessing(int TimestepInt, double phystime, double dt, TimestepNumber TimestepNo) {

            if (this.Control.CheckInterfaceProps) {
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
                if (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
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
                if (Log != null && this.Control.LogValues != XRheology_Control.LoggingValues.None && base.MPIRank == 0 && (TimestepNo.MajorNumber % this.Control.LogPeriod == 0))
                    try {
                        WriteLogLine(TimestepNo, phystime + dt);
                    } catch (Exception e) {
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
                    for (int i = 0; i < Length; i++)
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
                    for (int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return new double[] { volume, surface };

        }


        public double GetContactLineLength() {

            double CL_length = 0.0;

            if (this.LsTrk.GridDat.SpatialDimension == 3) {

                var metrics = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder);

                XQuadSchemeHelper SchemeHelper = metrics.XQuadSchemeHelper;
                EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"));

                var QuadDom = SurfaceElement_Edge.Domain;
                var boundaryCutEdge = QuadDom.Intersect(this.GridData.GetBoundaryEdgeMask());

                var innerDom = QuadDom.Except(this.GridData.GetBoundaryEdgeMask());

                System.Collections.BitArray lowerBits = new System.Collections.BitArray(((GridData)this.GridData).Edges.Count);
                foreach (Chunk cnk in boundaryCutEdge) {
                    for (int iE = cnk.i0; iE < cnk.JE; iE++) {
                        if (((GridData)this.GridData).Edges.EdgeTags[iE] == 1) {
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
                        for (int i = 0; i < length; i++)
                            CL_length += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
            }

            return CL_length;

        }


        public double[] ComputeBenchmarkQuantities_RisingBubble() {

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
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

            } //*/
              // Stress error
              // =============================================================
              //if (this.Control.ExSol_Stress != null) {

            //    Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
            //    double[] L2Error = new double[3];


            //    //must be changed for multiphase!!!
            //    foreach (var spc in this.LsTrk.SpeciesNames) {
            //        L2Error_Species.Add(spc, new double[3]);

            //        SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
            //        var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

            //        //stressXX
            //        //ConventionalDGField Stress_XX = ((XDGField)this.StressXX).GetSpeciesShadowField(spc);
            //        //L2Error_Species[spc] = Stress_XX.L2Error(this.Control.ExSol_Stress[spc].Vectorize(0.0), order, scheme);
            //        //L2Error[0] += L2Error_Species[spc].Pow2();

            //        //base.QueryHandler.ValueQuery("L2err_" + VariableNames.Stress_XX + "#" + spc, L2Error_Species[spc], true);

            //        //______________________________________________________________________________________________________________________
            //        for (int d = 0; d < D; d++) {
            //            ConventionalDGField Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);

            //            L2Error_Species[spc][d] = Vel_d.L2Error(this.Control.ExactSolutionVelocity[spc][d].Vectorize(time), order, scheme);
            //            L2Error[d] += L2Error_Species[spc][d].Pow2();

            //            base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);
            //        }
            //    }
            //    L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

            //    for (int d = 0; d < D; d++) {
            //        base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
            //        Ret[d] = L2Error[d];
            //    }
            //    //______________________________________________________________________________________________________________________



            //    L2Error[0] = this.StressXX.L2Error(this.Control.ExSol_Stress[0].Vectorize(0.0), order);
            //    L2Error[1] = this.StressXY.L2Error(this.Control.ExSol_Stress[1].Vectorize(0.0), order);
            //    L2Error[2] = this.StressYY.L2Error(this.Control.ExSol_Stress[2].Vectorize(0.0), order);

            //    base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressXX, L2Error[0], true);
            //    base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressXY, L2Error[1], true);
            //    base.QueryHandler.ValueQuery("L2err_" + VariableNames.StressYY, L2Error[2], true);

            //    Console.WriteLine("L2err " + VariableNames.StressXX + " is " + L2Error[0]);
            //    Console.WriteLine("L2err " + VariableNames.StressXY + " is " + L2Error[1]);
            //    Console.WriteLine("L2err " + VariableNames.StressYY + " is " + L2Error[2]);
            //}


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

            if (this.Control.WriteInterfaceP) {
                LogInterfaceP = base.DatabaseDriver.FsDriver.GetNewLog("InterfaceP", sessionID);
                string header = String.Format("{0}\t{1}\t{2}", "#timestep", "#time", "interfacePoints");
                LogInterfaceP.WriteLine(header);
                LogInterfaceP.Flush();
            }

            switch (this.Control.LogValues) {
                case XRheology_Control.LoggingValues.Wavelike: {

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
                case XRheology_Control.LoggingValues.RisingBubble: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("BenchmarkQuantities_RisingBubble", sessionID);
                        string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "#time", "area", "center of mass - x", "center of mass - y", "circularity", "rise velocity");
                        Log.WriteLine(header);
                        Log.Flush();

                        return;
                    }
                case XRheology_Control.LoggingValues.MovingContactLine: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("ContactAngle", sessionID);
                        string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "#time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-angle");
                        Log.WriteLine(header);
                        Log.Flush();

                        return;
                    }
                case XRheology_Control.LoggingValues.CapillaryHeight: {

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
                case XRheology_Control.LoggingValues.Wavelike: {

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
                case XRheology_Control.LoggingValues.RisingBubble: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime, BmQ_RB[0], BmQ_RB[1], BmQ_RB[2], BmQ_RB[3], BmQ_RB[5]);
                        Log.WriteLine(line);
                        Log.Flush();

                        return;
                    }
                case XRheology_Control.LoggingValues.MovingContactLine: {

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
                                    if (ResultsOfIntegration[i, 2 * D] != 0.0) {
                                        contactAngles.Add(Math.Abs(ResultsOfIntegration[i, 2 * D]));
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
                case XRheology_Control.LoggingValues.CapillaryHeight: {

                        MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);

                        double h_min = double.MaxValue, x_pos = 0.0;
                        for (int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                            if (InterfacePoints[i, 1] < h_min) {
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
                case XRheology_Control.LoggingValues.Wavelike: {

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
                case XRheology_Control.LoggingValues.RisingBubble: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        base.QueryResultTable.LogValue("area", BmQ_RB[0]);
                        base.QueryResultTable.LogValue("yCM", BmQ_RB[2]);
                        base.QueryResultTable.LogValue("circ", BmQ_RB[3]);
                        base.QueryResultTable.LogValue("riseV", BmQ_RB[5]);

                        return;
                    }
                case XRheology_Control.LoggingValues.LinelikeLS: {
                        break;
                    }
                case XRheology_Control.LoggingValues.CirclelikeLS: {

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
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
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
                if (m_energyBcMap == null) {
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
                if (m_CurrentEnergySolution == null) {
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
                if (m_CurrentEnergyResidual == null) {
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
                    if (this.Control.PhysicalParameters.IncludeConvection) {

                        var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];

                        // kinetic energy
                        var convK = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyConvectionInBulk(D, energyBcMap, rhoA, rhoB, LFFA, LFFB, LsTrk);
                        comps.Add(convK); // Bulk component


                        bool movingmesh;
                        switch (this.Control.Timestepper_LevelSetHandling) {
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

                    if (dntParams.UseGhostPenalties) {
                        var ViscPenalty = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyLaplace(penalty * 1.0, 0.0, energyBcMap, D, muA, muB);
                        KineticEnergyBalanceOperator.GhostEdgesOperator.EquationComponents[CodName[0]].Add(ViscPenalty);
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
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((VelParam[0] as XDGField).GetSpeciesShadowField(Spc));
                    (GradVelX[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
                }
            }
            GradVelX.ForEach(F => F.CheckForNanOrInf(true, true, true));

            VectorField<DGField> GradVelY = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((VelParam[1] as XDGField).GetSpeciesShadowField(Spc));
                    (GradVelY[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
                }
            }
            GradVelY.ForEach(F => F.CheckForNanOrInf(true, true, true));

            // pressure and gradient
            var PressMap = new CoordinateMapping(this.Pressure);
            DGField[] PressParam = PressMap.Fields.ToArray();

            VectorField<DGField> PressGrad = new VectorField<DGField>(D, PressParam[0].Basis, "PressureGrad", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
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
            if (OpMtx != null) {

                var mtxBuilder = KineticEnergyBalanceOperator.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping, SpcToCompute);

                mtxBuilder.time = phystime;

                foreach (var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                }

                mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin eval = KineticEnergyBalanceOperator.GetEvaluatorEx(LsTrk,
                    CurrentState.ToArray(), Params, Mapping,
                    SpcToCompute);

                foreach (var kv in AgglomeratedCellLengthScales) {
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

    }
}




















