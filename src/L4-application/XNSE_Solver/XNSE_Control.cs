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
using System.Runtime.Serialization;
using System.Linq;
using System.Text;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.NSECommon;

using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.EllipticExtension;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.Timestepping;
using Newtonsoft.Json;
using BoSSS.Solution.EnergyCommon;

namespace BoSSS.Application.XNSE_Solver {


    /// <summary>
    /// 
    /// </summary>
    [DataContract]
    [Serializable]
    public class XNSE_Control : AppControlSolver {

        /// <summary>
        /// Ctor.
        /// </summary>
        public XNSE_Control() {
            base.LinearSolver.NoOfMultigridLevels = 1;
            base.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            //shift of Solver Information
            base.LinearSolver.MaxKrylovDim = 100; //Solver_MaxKrylovDim;
            base.LinearSolver.MaxSolverIterations = 2000; //Solver_MaxIterations
            base.LinearSolver.MinSolverIterations = 4; //Solver_MinIterations
            base.LinearSolver.ConvergenceCriterion = 1.0e-10; //Solver_ConvergenceCriterion
            base.LinearSolver.SolverCode = LinearSolverCode.classic_mumps; //LinearSolver
            base.NonLinearSolver.MaxSolverIterations = 2000; //Solver_MaxIterations
            base.NonLinearSolver.MinSolverIterations = 4; //Solver_MinIterations
            base.NonLinearSolver.ConvergenceCriterion = 1.0e-10; //Solver_ConvergenceCriterion
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard; //NonLinearSolver
        }

        /// <summary>
        /// Type of <see cref="XNSE_SolverMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XNSE_SolverMain);
        }

        /// <summary>
        /// 
        /// </summary>
        public override void SetDGdegree(int p) {
            SetFieldOptions(p, Math.Max(2, p));
        }


        /// <summary>
        /// 
        /// </summary>
        public void SetFieldOptions(int VelDegree, int LevSetDegree, FieldOpts.SaveToDBOpt SaveFilteredVelocity =  FieldOpts.SaveToDBOpt.TRUE, FieldOpts.SaveToDBOpt SaveCurvature = FieldOpts.SaveToDBOpt.TRUE) {
            FieldOptions.Add("Velocity*", new FieldOpts() {
                Degree = VelDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add("FilteredVelocity*", new FieldOpts() {
                SaveToDB = SaveFilteredVelocity
            });
            FieldOptions.Add("SurfaceForceDiagnostic*", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.FALSE
            });
            FieldOptions.Add(VariableNames.Pressure, new FieldOpts() {
                Degree = VelDegree - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add("Phi", new FieldOpts() {
                Degree = LevSetDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = LevSetDegree*2,
                SaveToDB = SaveCurvature
            });
            FieldOptions.Add(VariableNames.Temperature, new FieldOpts() {
                Degree = VelDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }

       

        /// <summary>
        /// Width of the narrow band.
        /// </summary>
        [DataMember]
        public int LS_TrackerWidth = 1;

        /// <summary>
        /// different implementations for the level indicator 
        /// </summary>
        public enum RefinementStrategy {

            /// <summary>
            /// same refinement level on near band
            /// </summary>
            constantInterface,

            /// <summary>
            /// additional refinement on cells with high curvature
            /// </summary>
            CurvatureRefined,

            /// <summary>
            /// additional refinement at contact line
            /// </summary>
            ContactLineRefined
        }

        /// <summary>
        /// See <see cref="LoggingValues"/>
        /// </summary>
        [DataMember]
        public RefinementStrategy RefineStrategy = RefinementStrategy.constantInterface;

        /// <summary>
        /// desired minimum refinement level at interface
        /// </summary>
        [DataMember]
        public int BaseRefinementLevel = 1;

        /// <summary>
        /// maximum refinement level including additional refinement (contact line, curvature, etc.)
        /// </summary>
        [DataMember]
        public int RefinementLevel = 1;


        /// <summary>
        /// additional refinement of the navier slip boundary 
        /// </summary>
        [DataMember]
        public bool RefineNavierSlipBoundary = false;


        [DataMember]
        public int ReInitPeriod = 0;

        /// <summary>
        /// Expert options regarding the spatial discretization.
        /// </summary>
        [DataMember]
        public DoNotTouchParameters AdvancedDiscretizationOptions = new DoNotTouchParameters();

        /// <summary>
        /// Viscosity, density and surface tension.
        /// </summary>
        [DataMember]
        public PhysicalParameters PhysicalParameters = new PhysicalParameters() {
            Material = true,
            IncludeConvection = false,
            mu_A = 1.0,
            mu_B = 1.0,
            rho_A = 1.0,
            rho_B = 1.0,
            Sigma = 0.0
        };       

        /// <summary>
        /// Only for debugging purpose:
        /// solver is turned of and residual of initial value/exact solution is evaluated, used to 
        /// test the consistency of the implementation.
        /// </summary>
        public bool SkipSolveAndEvaluateResidual = false;

        /// <summary>
        /// Data to be written in LogFile
        /// </summary>
        public enum LoggingValues {

            /// <summary>
            /// no data will be written
            /// </summary>
            None,

            /// <summary>
            /// for elemental test programm with line like interfaces
            /// </summary>
            LinelikeLS,

            /// <summary>
            /// for elemental test programm with circle like interfaces
            /// </summary>
            CirclelikeLS,

            /// <summary>
            /// for wavelike simulation as CapillaryWave, RT-Instability
            /// interface height (interface points)
            /// </summary>
            Wavelike,

            /// <summary>
            /// for the benchmark quantities of the Rising Bubble testcase
            /// </summary>
            RisingBubble,

            /// <summary>
            /// contact points and corresponding contact angle
            /// </summary>
            MovingContactLine,

            /// <summary>
            /// height of a rising capillary in a tube
            /// </summary>
            CapillaryHeight,

            /// <summary>
            /// Evaporative mass flux and speed of displacement (Line interface)
            /// </summary>
            EvaporationL,

            /// <summary>
            /// Evaporative mass flux and speed of displacement (circle interface)
            /// </summary>
            EvaporationC
        }

        /// <summary>
        /// See <see cref="LoggingValues"/>
        /// </summary>
        [DataMember]
        public LoggingValues LogValues = LoggingValues.None;

        [DataMember]
        public int LogPeriod = 1;

        public bool WriteInterfaceP = false;

        public bool TestMode = false;


        /// <summary>
        /// See <see cref="TimestepperInit"/>
        /// </summary>
        [DataMember]
        public TimeStepperInit Timestepper_BDFinit = TimeStepperInit.SingleInit;

        /// <summary>
        /// defines the number of incremental timesteps in one gloabl timestep (for incrementInit)
        /// </summary>
        public int incrementTimesteps = 1;

        /// <summary>
        /// See <see cref="LevelSetHandling"/>
        /// </summary>
        [DataMember]
        public LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

        /// <summary>
        /// underrelaxation of the level set movement in case of coupled iterative
        /// </summary>
        public double LSunderrelax = 1.0;


        /// <summary>
        /// See <see cref="LevelSetEvolution"/>.
        /// </summary>
        [DataMember]
        public LevelSetEvolution Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

        /// <summary>
        /// Options for the initialization of the Fourier Level-set
        /// </summary>
        [DataMember]
        public FourierLevSetControl FourierLevSetControl;

        /// <summary>
        /// array of additional parameter values for some testcases
        /// </summary>
        [DataMember]
        public double[] AdditionalParameters;


        /// <summary>
        /// The termination criterion for fully coupled/implicit level-set evolution.
        /// </summary>
        [DataMember]
        public double LevelSet_ConvergenceCriterion = 1.0e-6;

      

        /// <summary>
        /// See <see cref="ContinuityProjection"/>
        /// </summary>
        [DataMember]
        public ContinuityProjectionOption LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;

        /// <summary>
        /// Enforce the level-set to be globally conservative, by adding a constant to the level-set field
        /// </summary>
        public bool EnforceLevelSetConservation = false;


        /// <summary>
        /// if true, kinetic energy equation will be solved 
        /// </summary>
        [DataMember]
        public bool solveKineticEnergyEquation = false;

        /// <summary>
        /// discretization option for the visocus source terms of the kinetic energy equation
        /// </summary>
        [DataMember]
        public KineticEnergyViscousSourceTerms kinEViscousDiscretization;

        /// <summary>
        /// discretization option for the pressure source terms of the kinetic energy equation
        /// </summary>
        [DataMember]
        public KineticEnergyPressureSourceTerms kinEPressureDiscretization;

        /// <summary>
        /// switch for the pressure term in the Dissipation term
        /// </summary>
        [DataMember]
        public bool withDissipativePressure;

        /// <summary>
        /// Block-Precondition for the kinetic-Energy-block
        /// </summary>
        public MultigridOperator.Mode KineticEnergyeBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;

        /// <summary>
        /// If true, various energy properties will be evaluated in every cycle.
        /// </summary>
        [DataMember]
        public bool ComputeEnergyProperties = false;

        /// <summary>
        /// if true, the jump condition for mass, momentum and energy will be checked
        /// </summary>
        [DataMember]
        public bool CheckJumpConditions = false;

        /// <summary>
        /// if true, the mass conservation and the surface changerate is checked
        /// </summary>
        [DataMember]
        public bool CheckInterfaceProps = false;

        /// <summary>
        /// Registers all utility (also energy) fields to IOFields
        /// </summary>
        [DataMember]
        public bool RegisterUtilitiesToIOFields = false;
        
        /// <summary>
        /// average method for interface values
        /// </summary>
        public enum InterfaceAveraging {

            /// <summary>
            /// arithmetic mean
            /// </summary>
            mean,

            /// <summary>
            /// density weighted average
            /// </summary>
            density,

            /// <summary>
            /// viscosity weighted average
            /// </summary>
            viscosity

        }

        /// <summary>
        /// See <see cref="InterfaceAveraging"/>
        /// </summary>
        public InterfaceAveraging InterAverage = InterfaceAveraging.density;



        /// <summary>
        /// An explicit expression of the Level-set over time.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> Phi;

        /// <summary>
        /// Exact solution for velocity, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public IDictionary<string, Func<double[], double, double>[]> ExactSolutionVelocity;

        /// <summary>
        /// Exact solution, pressure, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public IDictionary<string, Func<double[], double, double>> ExactSolutionPressure;

        /// <summary>
        /// Control Options for ReInit
        /// </summary>
        public EllipticReInitAlgoControl ReInitControl = new EllipticReInitAlgoControl();

        /// <summary>
        /// Control Options for ExtVel
        /// </summary>
        public EllipticExtVelAlgoControl EllipticExtVelAlgoControl = new EllipticExtVelAlgoControl();



        /// <summary>
        /// switch for the computation of the coupled heat solver
        /// </summary>
        [DataMember]
        public bool solveCoupledHeatEquation = false;

        /// <summary>
        /// switch for advanced parameter Update for nonlinear solver
        /// </summary>
        [DataMember]
        public bool useSolutionParamUpdate = false;

        /// <summary>
        /// only available if no heat equation is solved
        /// </summary>
        public Func<double[], double, double> prescribedMassflux_Evaluator;

        [DataMember]
        public IBoundaryAndInitialData prescribedMassflux;

        /// <summary>
        /// implementations for the conductivity part (laplace operator) of the heat equation 
        /// </summary>
        [DataMember]
        public ConductivityInSpeciesBulk.ConductivityMode conductMode = ConductivityInSpeciesBulk.ConductivityMode.SIP;

        /// <summary>
        /// Block-Precondition for the Temperature-block
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode TemperatureBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;


        /// <summary>
        /// function for the disjoining pressure
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> DisjoiningPressureFunc;

        /// <summary>
        /// density, heat capacity and thermal conductivity
        /// </summary>
        [DataMember]
        public ThermalParameters ThermalParameters = new ThermalParameters() {
            rho_A = 1.0,
            rho_B = 1.0,
            c_A = 1.0,
            c_B = 1.0,
            k_A = 1.0,
            k_B = 1.0,
        };
    }
}
