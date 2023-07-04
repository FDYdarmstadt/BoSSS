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

            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.OneTimeSetUp();
            ////BoSSS.Application.XNSE_Solver.Tests.UnitTest.PolynomialTestForConvectionTest(3, 0, false);
            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.TestCapillaryWave();
            ////BoSSS.Application.XNSE_Solver.Tests.ElementalTestProgramm.LineMovementTest(LevelSetEvolution.ScalarConvection, LevelSetHandling.Coupled_Once, TimeSteppingScheme.ImplicitEuler, 0.5);
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
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {

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

            if (Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                m_HMForder = 2 * degU * (this.Control.PhysicalParameters.IncludeConvection ? 4 : 3) + 1;
            } else {
                m_HMForder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 4 : 3);
            }


            // Create Spatial Operator
            // ======================= 

            XOpConfig = new XRheology_OperatorConfiguration(this.Control);

            XRheology_Operator = new XRheology_OperatorFactory(XOpConfig, this.LsTrk, this.m_HMForder, this.BcMap, stressDegree, degU);
            //XRheology_Operator = new XNSE_OperatorFactory(XOpConfig, this.LsTrk, this.m_HMForder, this.BcMap, degU);

            {
                var tempOp = new ConstantXTemporalOperator(XRheology_Operator.Xop, 0.0);
                foreach (var kv in this.MassScale) {
                    tempOp.DiagonalScale[LsTrk.GetSpeciesName(kv.Key)].SetV(kv.Value.ToArray());
                }
                XRheology_Operator.Xop.TemporalOperator = tempOp;

            }


            #endregion

            #region Create Timestepper
            // =======================
            if (L == null) {

                switch (this.Control.TimeSteppingScheme) {
                    case TimeSteppingScheme.RK_ImplicitEuler: {
                            rksch = RungeKuttaScheme.ImplicitEuler;
                            break;
                        }
                    case TimeSteppingScheme.RK_CrankNic: {
                            rksch = RungeKuttaScheme.CrankNicolson;
                            break;
                        }
                    case TimeSteppingScheme.CrankNicolson: {
                            //do not instantiate rksch, use bdf instead
                            bdfOrder = -1;
                            break;
                        }
                    case TimeSteppingScheme.ImplicitEuler: {
                            //do not instantiate rksch, use bdf instead
                            bdfOrder = 1;
                            break;
                        }
                    default: {
                            if (this.Control.TimesteppingMode.ToString().StartsWith("BDF")) {
                                //do not instantiate rksch, use bdf instead
                                bdfOrder = Convert.ToInt32(this.Control.TimesteppingMode.ToString().Substring(3));
                                break;
                            } else
                                throw new NotImplementedException();
                        }

                }


                if (rksch == null) {
                    m_BDF_Timestepper = new XdgBDFTimestepping(
                        this.CurrentSolution.Mapping.Fields,
                        this.XRheology_Operator.Xop.InvokeParameterFactory(this.CurrentSolution.Mapping.Fields),
                        this.CurrentResidual.Mapping.Fields,
                        LsTrk,
                        true,
                        DelComputeOperatorMatrix, this.XRheology_Operator.Xop, () => new LevelSetTimeIntegratorWrapper(this),
                        (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? bdfOrder : 1,
                        this.Control.Timestepper_LevelSetHandling,
                        this.XOpConfig.mmsd,
                        (this.Control.PhysicalParameters.IncludeConvection) ? SpatialOperatorType.Nonlinear : SpatialOperatorType.LinearTimeDependent,
                        this.MultigridOperatorConfig, base.MultigridSequence,
                        this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder,
                        this.Control.AgglomerationThreshold,
                        true,
                        this.Control.NonLinearSolver,
                        this.Control.LinearSolver
                        );

                    m_BDF_Timestepper.m_ResLogger = base.ResLogger;
                    m_BDF_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                    m_BDF_Timestepper.Timestepper_Init = (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? this.Control.Timestepper_BDFinit : TimeStepperInit.SingleInit;
                    m_BDF_Timestepper.IterUnderrelax = this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? this.Control.LSunderrelax : 1.0;

                    m_BDF_Timestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                    //m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;


                    // solver 
                    this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations; //m_BDF_Timestepper.config_NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;

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
                    base.SetInitial(0);
                    this.InitLevelSet();
                }


                m_BDF_Timestepper.DataRestoreAfterBalancing(L,
                    ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure, this.StressXX, this.StressXY, this.StressYY),
                    this.XRheology_Operator.Xop.InvokeParameterFactory(this.CurrentSolution.Mapping.Fields),
                    ArrayTools.Cat<DGField>(this.XDGvelocity.ResidualMomentum.ToArray(), this.ResidualContinuity, this.ResidualStressXX, ResidualStressXY, this.ResidualStressYY),
                    this.LsTrk, this.MultigridSequence, this.XRheology_Operator.Xop);

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



        void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime, int LsTrkHistIdx) {
            if(LsTrkHistIdx != 1)
                throw new NotSupportedException("No supported for anything but the current tracker time level.");

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


            /*  not required anymore; 
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
            */


            // transform from RHS to Affine
            OpAffine.ScaleV(-1.0);

            //OpMtx.SaveToTextFile("OpMatrix");
        }

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

                            if (Control.PhysicalParametersRheology.Weissenberg_a != 0.0 || Control.PhysicalParametersRheology.Weissenberg_b != 0.0) {

                                if (Control.WeissenbergIncrement != 0.0) {
                                    NoIncrementTimestep = 1;
                                    if(Control.PhysicalParametersRheology.Weissenberg_a > Control.PhysicalParametersRheology.Weissenberg_b)
                                        NoIncrementTimestep = (int)(Control.PhysicalParametersRheology.Weissenberg_a / Control.WeissenbergIncrement);
                                    else if(Control.PhysicalParametersRheology.Weissenberg_b > Control.PhysicalParametersRheology.Weissenberg_a)
                                        NoIncrementTimestep = (int)(Control.PhysicalParametersRheology.Weissenberg_b / Control.WeissenbergIncrement);
                                    else if (Control.PhysicalParametersRheology.Weissenberg_b == Control.PhysicalParametersRheology.Weissenberg_a)
                                        NoIncrementTimestep = (int)(Control.PhysicalParametersRheology.Weissenberg_a / Control.WeissenbergIncrement);
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

                                if (currentWeissenberg[0] < Control.PhysicalParametersRheology.Weissenberg_a) {
                                    currentWeissenberg[0] = currentWeissenberg[0] + Control.WeissenbergIncrement;
                                    Console.WriteLine();
                                    Console.WriteLine("Raise Weissenberg number A to " + currentWeissenberg[0]);
                                    Console.WriteLine();
                                }

                                if (currentWeissenberg[1] < Control.PhysicalParametersRheology.Weissenberg_b) {
                                    currentWeissenberg[1] = currentWeissenberg[1] + Control.WeissenbergIncrement;
                                    Console.WriteLine();
                                    Console.WriteLine("Raise Weissenberg number B to " + currentWeissenberg[1]);
                                    Console.WriteLine();
                                }

                            }
                        } else {
                            //current Weissenberg is set to the HIGHER value... DIRTY HACK AT THE MOMENT!

                                currentWeissenberg[0] = Control.PhysicalParametersRheology.Weissenberg_a;
                                currentWeissenberg[1] = Control.PhysicalParametersRheology.Weissenberg_b;


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

                                    var agg = LsTrk.GetAgglomerator(LsTrk.SpeciesIdS.ToArray(), m_HMForder, this.Control.AgglomerationThreshold);

                                    DelComputeOperatorMatrix(SaddlePointMatrix, AffineDummy, this.CurrentSolution.Mapping,
                                    this.CurrentSolution.Mapping.Fields.ToArray(), agg.CellLengthScales, 0.0, 1);

                                    AggregationGridBasis[][] MgBasis = AggregationGridBasis.CreateSequence(this.MultigridSequence, this.CurrentSolution.Mapping.BasisS);
                                    //todo: AsyncCallback update
                                    MgBasis.UpdateXdgAggregationBasis(agg);
                                    MultigridOperator mgOp = new MultigridOperator(MgBasis, CurrentSolution.Mapping,
                                        SaddlePointMatrix, this.MassFact.GetMassMatrix(CurrentSolution.Mapping, false),
                                        this.MultigridOperatorConfig, null);

                                    MsrMatrix FullMatrix = mgOp.OperatorMatrix.ToMsrMatrix();

                                    MsrMatrix DiffMatrix;
                                    {
                                        int[] VelVarIdx = new int[] { 3, 4, 5 };

                                        long[] USubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                                        long[] USubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                                        int L = USubMatrixIdx_Row.Length;

                                        DiffMatrix = new MsrMatrix(L, L, 1, 1);
                                        FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(long[]), USubMatrixIdx_Col, default(long[]));
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


        protected override void SetInitial(double t) {
            base.SetInitial(t);

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

                this.LsTrk.UpdateTracker(Time, incremental: true);

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
                            DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                            mode = this.Control.VelocityBlockPrecondMode,
                            VarIndex = new int[] { d }
                        };
                    }
                    // configuration for pressure
                    configs[iLevel][D] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                        mode = this.Control.PressureBlockPrecondMode,
                        VarIndex = new int[] { D }
                    };
                    //configurations for stresses
                    for (int d = 3; d < 6; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(1, pStr - iLevel) },
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

                    GridRefinementController gridRefinementController = new GridRefinementController((GridData)this.GridData, BlockedCells);
                    bool AnyChange = gridRefinementController.ComputeGridChange(LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
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
                    long oldJ = this.GridData.CellPartitioning.TotalLength;

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
        public override void DataBackupBeforeBalancing(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
        }




        #endregion


        // =========
        // level-set
        // =========
        #region level-set

        ///// <summary>
        ///// Information of the current Fourier Level-Set
        ///// DFT_coeff
        ///// </summary>
        //FourierLevSetBase Fourier_LevSet;

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

        /// <summary>
        /// wrapper introduced due to API change
        /// </summary>
        class LevelSetTimeIntegratorWrapper : ISlaveTimeIntegrator {

            public LevelSetTimeIntegratorWrapper(XRheology_SolverMain __owner) {
                m_owner = __owner;
            }
            XRheology_SolverMain m_owner;

            public void Pop() {
                throw new NotImplementedException();
            }

            public void Push() {
                m_owner.PushLevelSetAndRelatedStuff();
            }

            public double Update(DGField[] CurrentState, double time, double dt, double UnderRelax, bool incremental) {
                return m_owner.DelUpdateLevelSet(CurrentState, time, dt, UnderRelax, incremental);
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
            var surfElemVol = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spcId, 0);
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
                EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"), 0);

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
        /// encapsulated handling of query values
        /// </summary>
        public void LogQueryValue(double phystime) {

            base.QueryResultTable.LogValue("time", phystime);

            

        }


        #endregion


    }
}




















