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
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation;
using BoSSS.Application.XNSE_Solver.LoadBalancing;
using ilPSP;

namespace BoSSS.Application.XNSE_Solver {


    /// <summary>
    /// 
    /// </summary>
    [DataContract]
    [Serializable]
    public class XNSE_Control : SolverWithLevelSetUpdaterControl {

        /// <summary>
        /// Ctor.
        /// </summary>
        public XNSE_Control() {

            //base.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            //shift of Solver Information
            base.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig(); //LinearSolver
            base.NonLinearSolver.MaxSolverIterations = 2000; //Solver_MaxIterations
            base.NonLinearSolver.MinSolverIterations = 4; //Solver_MinIterations
            base.NonLinearSolver.ConvergenceCriterion = 0.0; //Solver_ConvergenceCriterion: solve as accurate as possible. Don't change this, Grüße von FK!
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton; //NonLinearSolver
            base.TimesteppingMode = AppControl._TimesteppingMode.Steady;


            // fk, 20mar23: Dynamic load balancing: for XNSE, turn on by default
            base.DynamicLoadBalancing_CellCostEstimators.Clear();
            base.DynamicLoadBalancing_CellCostEstimators.Add(new Loadbalancing.XNSECellCostEstimator());
            base.DynamicLoadBalancing_On = true;
            base.DynamicLoadBalancing_RedistributeAtStartup = true;

        }


        /// <summary>
        /// Activation of second level-set (fluid/solid boundary)
        /// </summary>
        [DataMember]
        virtual public bool UseImmersedBoundary {
            get;
            set;
        }

        /// <summary>
        /// Temporary. Suggestion: Move Rigid body benchmarks to FSI solver in future.
        /// Sets Parameter for Rigidbody.
        /// </summary>
        [DataMember]
        public XRigid Rigidbody = new XRigid();

        /// <summary>
        /// Sets Field options for residual fields,
        /// residual fields are now written to database.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <param name="k">Velocity Degree</param>
        public void SetOptionsResFields (int k) {
            this.FieldOptions.Add("Residual-MomentumX", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            this.FieldOptions.Add("Residual-MomentumY", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            this.FieldOptions.Add("Residual-MomentumZ", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            this.FieldOptions.Add("Residual-ContiEq", new FieldOpts() {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }

        public void SetMaximalRefinementLevel(int maxLvl) {
            this.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = maxLvl });
        }

        /// <summary>
        /// - default (false): preconditioning for velocity and pressure is determined by 
        ///   <see cref="VelocityBlockPrecondMode"/> and <see cref="PressureBlockPrecondMode"/>, respectively;
        /// - true: former options are ignored, Schur complement is used instead.
        /// </summary>
        [DataMember]
        public bool UseSchurBlockPrec = false;

        /// <summary>
        /// Type of <see cref="XNSE"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XNSE<XNSE_Control>);
        }

        /// <summary>
        /// 
        /// </summary>
        public override void SetDGdegree(int p) {
            SetFieldOptions(p, Math.Max(2, p));
        }

        /// <summary>
        /// Set Field Options, i.e. the DG degrees
        /// </summary>
        public void SetFieldOptions(int VelDegree, int LevSetDegree, FieldOpts.SaveToDBOpt SaveCurvature = FieldOpts.SaveToDBOpt.TRUE) {
            if(VelDegree < 1)
                throw new ArgumentOutOfRangeException("Velocity degree must be 1 at minimum.");
            
            FieldOptions.Add("Velocity*", new FieldOpts() {
                Degree = VelDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.Pressure, new FieldOpts() {
                Degree = VelDegree - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetDG, new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetCG, new FieldOpts() {
                Degree = LevSetDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetDGidx(1), new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetCGidx(1), new FieldOpts() {
                Degree = LevSetDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.Curvature, new FieldOpts() {
                Degree = LevSetDegree*2,
                SaveToDB = SaveCurvature
            });
        }

        /*
        public void SetDGdegree2(int p) {
            FieldOptions.Add(VariableNames.VelocityX, new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.VelocityY, new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.Pressure, new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }
        

        /// <summary>
        /// Allows to set DG degree of level set and flow solver 
        /// </summary>
        /// <param name="VelDegree"></param>
        /// <param name="LevSetDegree"></param>
        public void SetFieldOptions2(int VelDegree, int LevSetDegree) {
            FieldOptions.Add(VariableNames.VelocityX, new FieldOpts() {
                Degree = VelDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.VelocityY, new FieldOpts() {
                Degree = VelDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
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
                Degree = LevSetDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }
        */

        [DataMember]
        public string methodTagLS;

        /// <summary>
        /// 
        /// </summary>
        [Obsolete] // really?
        public void SetLevelSetMethod(int method, FourierLevSetControl _FourierControl = null) {

            LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            switch (method) {
                case 0: {
                        goto default;
                    }
                case 1: {
                        // fast marching with Curvature and default filtering 
                        methodTagLS = "FastMarchCurv";
                        Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                        FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;
                        AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
                        AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;
                        break;
                    }
                case 2: {
                        // Extension Velocity with Laplace Beltrami without filtering
                        methodTagLS = "ExtVelLB";
                        Option_LevelSetEvolution = LevelSetEvolution.ExtensionVelocity;
                        EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
                        //AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
                        AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
                        EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
                        //fullReInit = true;
                        break;
                    }
                case 3: {
                        // Extension Velocity with Curvature and default filtering 
                        methodTagLS = "ExtVelCurv";
                        Option_LevelSetEvolution = LevelSetEvolution.ExtensionVelocity;
                        EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
                        AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
                        AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;
                        EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
                        fullReInit = true;
                        break;
                    }
                case 4: {
                        methodTagLS = "Fourier";
                        FourierLevSetControl = _FourierControl;
                        Option_LevelSetEvolution = LevelSetEvolution.Fourier;
                        AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
                        FourierLevSetControl.Timestepper = FourierLevelSet_Timestepper.RungeKutta1901;
                        break;
                    }
                default: {
                        // (standard) fast marching with Laplace Beltrami without filtering
                        methodTagLS = "FastMarchLB";
                        Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                        FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;
                        AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
                        break;
                    }
            }
        }

        /* killed by fk:
         * the following adds a configuration redundancy, which is always a recipe for confusion
        /// <summary>
        /// switches off all plotCurrentState calls
        /// </summary>
        [DataMember]
        public bool switchOffPlotting = false;
        */


        /// <summary>
        /// different implementations for the level indicator 
        /// </summary>
        public enum RefinementStrategy {

            /// <summary>
            /// same refinement level on near band
            /// </summary>
            constantInterface,

            /// <summary>
            /// additional refinement on cells in phase A
            /// </summary>
            PhaseARefined,

            /// <summary>
            /// additional refinement on cells with high curvature
            /// </summary>
            CurvatureRefined,

            /// <summary>
            /// additional refinement at contact line
            /// </summary>
            ContactLineRefined,

            /// <summary>
            /// additional refinement at navier slip boundary
            /// </summary>
            NavierSlipRefined,

            /// <summary>
            /// additional refinement on near band cells for high velocity gradients
            /// </summary>
            VelocityGradient
        }

        /// <summary>
        /// For legacy solver <see cref="Legacy.XNSE_SolverMain"/>
        /// See <see cref="LoggingValues"/>
        /// </summary>
        [DataMember]
        public RefinementStrategy RefineStrategy = RefinementStrategy.constantInterface;

        /// <summary>
        /// For legacy solver <see cref="Legacy.XNSE_SolverMain"/>
        /// desired minimum refinement level at interface
        /// </summary>
        [DataMember]
        public int BaseRefinementLevel = 0;

        /// <summary>
        /// For legacy solver <see cref="Legacy.XNSE_SolverMain"/>
        /// maximum refinement level including additional refinement (contact line, curvature, etc.)
        /// </summary>
        [DataMember]
        public int RefinementLevel = 0;


        /// <summary>
        /// For legacy solver <see cref="Legacy.XNSE_SolverMain"/>
        /// additional refinement of the navier slip boundary 
        /// </summary>
        [DataMember]
        public bool RefineNavierSlipBoundary = false;


        /// <summary>
        /// option for clearing the velocities for restart
        /// </summary>
        [DataMember]
        public bool ClearVelocitiesOnRestart = false;

        [DataMember]
        public bool ReInitOnRestart = false;


        [DataMember]
        public bool adaptiveReInit = false;

        [DataMember]
        public bool InitSignedDistance = false;

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
            IncludeDiffusion = true,
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
        [DataMember]
        public bool SkipSolveAndEvaluateResidual = false;

        /// <summary>
        /// Terminates the simulation if the linear or nonlinear solver fails to converge
        /// </summary>
        [DataMember]
        public bool FailOnSolverFail = true;


        /// <summary>
        /// See <see cref="TimestepperInit"/>
        /// </summary>
        [DataMember]
        public TimeStepperInit Timestepper_BDFinit = TimeStepperInit.SingleInit;

        ///// <summary>
        ///// defines the number of incremental timesteps in one global timestep (for incrementInit)
        ///// </summary>
        //public int incrementTimesteps = 1;

       
        /// <summary>
        /// array of additional parameter values for some testcases
        /// </summary>
        [DataMember]
        public double[] AdditionalParameters;

        /// <summary>
        /// amplitude values for wave-like interfaces
        /// </summary>
        [DataMember]
        public double[] prescribedLSwaveData;

        /// <summary>
        /// Enforce the level-set to be globally conservative, by adding a constant to the level-set field
        /// </summary>
        [DataMember]
        public bool EnforceLevelSetConservation = false;


        /// <summary>
        /// if true, kinetic energy equation will be solved as postprocessing
        /// </summary>
        [DataMember]
        public bool solveKineticEnergyEquation = false;

        /// <summary>
        /// if false, the kinetic energy timestepping is one order higher than the flow solver
        /// </summary>
        [DataMember]
        public bool equalTimesteppingForKineticEnergy = true;

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
        /// average method for constructing the interface velocity
        /// </summary>
        public enum InterfaceVelocityAveraging {

            /// <summary>
            /// arithmetic mean
            /// </summary>
            mean = 1,

            /// <summary>
            /// density weighted average (recommended default value for most cases)
            /// </summary>
            density = 0,

            /// <summary>
            /// viscosity weighted average
            /// </summary>
            viscosity = 2,

            /// <summary>
            /// only take velocity from phase A
            /// </summary>
            phaseA = 3,

            /// <summary>
            /// only take velocity from phase B
            /// </summary>
            phaseB = 4

        }


        /// <summary>
        /// An explicit expression of the Level-set over time: \phi = f(x,y;t).
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> Phi;

        /// <summary>
        /// See <see cref="InterfaceAveraging"/>
        /// </summary>
        public InterfaceVelocityAveraging InterVelocAverage = InterfaceVelocityAveraging.density;

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
        /// Exact solution, temperature, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public IDictionary<string, Func<double[], double, double>> ExactSolutionTemperature;
        /// <summary>
        /// Exact solution, Mixture fraction, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public IDictionary<string, Func<double[], double, double>> ExactSolutionMixtureFraction;

        /// <summary>
        /// Time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        public ScalarFunctionTimeDep GetGravity(string species, int d) {
            bool bfound = this.InitialValues_EvaluatorsVec.TryGetValue(VariableNames.Gravity_d(d) + "#" + species, out var ret);
            if(!bfound)
                this.InitialValues_EvaluatorsVec.TryGetValue(VariableNames.Gravity_d(d), out ret);
            return ret;
        }

        /// <summary>
        /// Setting time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        public void SetGravity(string species, int d, IBoundaryAndInitialData g) {
            this.InitialValues[VariableNames.Gravity_d(d) + "#" + species] = g;
        }

        /// <summary>
        /// Setting time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        /// <remarks>
        /// Note: using the setter is not recommended when working with the job management system,
        /// since these values specified here cannot be serialized.
        /// Instead, <see cref="AppControl.InitialValues"/> or <see cref="SetGravity(string, int, IBoundaryAndInitialData)"/> should be used.
        /// </remarks>
        public void SetGravity(string species, int d, Func<double[], double, double> g) {
            this.InitialValues_Evaluators_TimeDep[VariableNames.Gravity_d(d) + "#" + species] = g;
        }

        /// <summary>
        /// Setting time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        /// <remarks>
        /// Note: using the setter is not recommended when working with the job management system,
        /// since these values specified here cannot be serialized.
        /// Instead, <see cref="AppControl.InitialValues"/> or <see cref="SetGravity(string, int, IBoundaryAndInitialData)"/> should be used.
        /// </remarks>
        public void SetGravity(string species, Func<double[], double, double>[] G) {
            for(int d = 0; d < G.Length; d++)
                this.InitialValues_Evaluators_TimeDep[VariableNames.Gravity_d(d) + "#" + species] = G[d];
        }


        /// <summary>
        /// Time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        public ScalarFunctionTimeDep GetVolumeForce(string species, int d) {
            bool bfound = this.InitialValues_EvaluatorsVec.TryGetValue(VariableNames.VolumeForce_d(d) + "#" + species, out var ret);
            if(!bfound)
                this.InitialValues_EvaluatorsVec.TryGetValue(VariableNames.VolumeForce_d(d), out ret);
            //Console.WriteLine("Using volume Force: " + (ret?.ToString() ?? "NIX"));
            return ret;
        }

        /// <summary>
        /// Setting time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        public void SetVolumeForce(string species, int d, IBoundaryAndInitialData g) {
            this.InitialValues[VariableNames.VolumeForce_d(d) + "#" + species] = g;
        }

        /// <summary>
        /// Setting time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        /// <remarks>
        /// Note: using the setter is not recommended when working with the job management system,
        /// since these values specified here cannot be serialized.
        /// Instead, <see cref="AppControl.InitialValues"/> or <see cref="SetVolumeForce(string, int, IBoundaryAndInitialData)"/> should be used.
        /// </remarks>
        public void SetVolumeForce(string species, int d, Func<double[], double, double> g) {
            this.InitialValues_Evaluators_TimeDep[VariableNames.VolumeForce_d(d) + "#" + species] = g;
        }

        /// <summary>
        /// Control Options for ReInit
        /// </summary>
        [DataMember]
        public EllipticReInitAlgoControl ReInitControl = new EllipticReInitAlgoControl();

        /// <summary>
        /// Control Options for ExtVel
        /// </summary>
        public EllipticExtVelAlgoControl EllipticExtVelAlgoControl = new EllipticExtVelAlgoControl();

        /// <summary>
        /// three-step reinitialization with preconditioning fast-marching
        /// </summary>
        [DataMember]
        public bool fullReInit = false;

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
        /// Contact lines and thin-films: 
        /// function for the disjoining pressure model
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
            alpha_A = 0.0,
            alpha_B = 0.0,
        };

        /// <summary>
        /// Used to active nonlinear solver even if convection is not included
        /// </summary>
        [DataMember]
        public bool NonlinearCouplingSolidFluid = false;

        ///// <summary>
        ///// 
        ///// </summary>
        //[DataMember]
        //public ClassifierType DynamicLoadbalancing_ClassifierType = ClassifierType.VoidCutNormal;


        /// <summary>
        /// Configuring <see cref="AppControl._TimesteppingMode.Steady"/> sets the <see cref="TimeSteppingScheme.ImplicitEuler"/>
        /// </summary>
        [JsonIgnore]
        public override _TimesteppingMode TimesteppingMode {
            get {
                return base.TimesteppingMode;
            }
            set {
                base.TimesteppingMode = value;
                if(value == _TimesteppingMode.Steady)
                    this.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        /// <summary>
        /// 
        /// </summary>
        public override bool Equals(object obj) {
            //System.Diagnostics. dbg_launch();
            if(!base.Equals(obj))
                return false;

            var other = obj as XNSE_Control;
            if(other == null)
                return false;

            if(!this.PhysicalParameters.Equals(other.PhysicalParameters))
                return false;

            if(!this.ThermalParameters.Equals(other.ThermalParameters))
                return false;


            return true;
        }

    }
}
