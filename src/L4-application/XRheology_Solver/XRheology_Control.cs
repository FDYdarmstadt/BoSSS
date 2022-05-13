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
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;

using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.EllipticExtension;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.XheatCommon;

namespace BoSSS.Application.XRheology_Solver {
    /// <summary>
    /// Control File For calculation with viscoelastic extra stress tensor
    /// </summary>
    [Serializable]
    [DataContract]
    public class XRheology_Control : XBase_Control {


        /// <summary>
        /// Ctor.
        /// </summary>
        public XRheology_Control() {
            base.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            //shift of Solver Information
            base.LinearSolver = LinearSolverCode.direct_mumps.GetConfig(); //LinearSolver
            base.NonLinearSolver.MaxSolverIterations = 50; //Solver_MaxIterations
            base.NonLinearSolver.MinSolverIterations = 50; //Solver_MinIterations
            base.NonLinearSolver.ConvergenceCriterion = 1.0e-10; //Solver_ConvergenceCriterion
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard; //NonLinearSolver
            base.NonLinearSolver.UnderRelax = 1.0; //UnderRelax
        }

        /// <summary>
        /// Type of <see cref="XRheology_SolverMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XRheology_SolverMain);
        }

        public void SetFieldOptions(int VelDegree, int LevSetDegree, FieldOpts.SaveToDBOpt SaveFilteredVelocity = FieldOpts.SaveToDBOpt.TRUE, FieldOpts.SaveToDBOpt SaveCurvature = FieldOpts.SaveToDBOpt.TRUE) {
            FieldOptions.Add(VariableNames.VelocityX, new FieldOpts() {Degree = VelDegree,SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.VelocityY, new FieldOpts() {Degree = VelDegree,SaveToDB = FieldOpts.SaveToDBOpt.TRUE});
            FieldOptions.Add("FilteredVelocityX", new FieldOpts() {SaveToDB = SaveFilteredVelocity});
            FieldOptions.Add("FilteredVelocityY", new FieldOpts() {SaveToDB = SaveFilteredVelocity});
            FieldOptions.Add(VariableNames.Pressure, new FieldOpts() {Degree = VelDegree - 1,SaveToDB = FieldOpts.SaveToDBOpt.TRUE});
            FieldOptions.Add(VariableNames.StressXX, new FieldOpts() {Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.StressXY, new FieldOpts() {Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.StressYY, new FieldOpts() {Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("PhiDG", new FieldOpts() {SaveToDB = FieldOpts.SaveToDBOpt.TRUE});
            FieldOptions.Add("Phi", new FieldOpts() {Degree = LevSetDegree,SaveToDB = FieldOpts.SaveToDBOpt.TRUE});
            FieldOptions.Add("Curvature", new FieldOpts() {Degree = LevSetDegree * 2,SaveToDB = SaveCurvature});
        }
 

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

        /// <summary>
        /// option for clearing the velocities for restart
        /// </summary>
        //[DataMember]
        //public bool ClearVelocitiesOnRestart = false;

        /// <summary>
        /// Only for debugging purpose:
        /// solver is turned of and residual of initial value/exact solution is evaluated, used to 
        /// test the consistency of the implementation.
        /// </summary>
        public bool SkipSolveAndEvaluateResidual = false;

        public bool FixedStreamwisePeriodicBC = false;

        ///// <summary>
        ///// Data to be written in LogFile
        ///// </summary>
        //public enum LoggingValues {

        //    /// <summary>
        //    /// no data will be written
        //    /// </summary>
        //    None,

        //    /// <summary>
        //    /// for elemental test programm with line like interfaces
        //    /// </summary>
        //    LinelikeLS,

        //    /// <summary>
        //    /// for elemental test programm with circle like interfaces
        //    /// </summary>
        //    CirclelikeLS,

        //    /// <summary>
        //    /// for wavelike simulation as CapillaryWave, RT-Instability
        //    /// interface height (interface points)
        //    /// </summary>
        //    Wavelike,

        //    /// <summary>
        //    /// for the benchmark quantities of the Rising Bubble testcase
        //    /// </summary>
        //    RisingBubble,

        //    /// <summary>
        //    /// contact points and corresponding contact angle
        //    /// </summary>
        //    MovingContactLine,

        //    /// <summary>
        //    /// height of a rising capillary in a tube
        //    /// </summary>
        //    CapillaryHeight
        //}



        public bool TestMode = false;

        public bool InterfaceTest = false;

        ///// <summary>
        ///// switch for the initialization of the <see cref="XdgBDFTimestepping"/> 
        ///// </summary>
        //public enum TimestepperInit {

        //    /// Initialization from a single timestep, i.e. if this time-stepper should use BDF4,
        //    /// it starts with BDF1, BDF2, BDF3 in the first, second and third time-step.
        //    SingleInit,

        //    /// same initialization for SingleInit, but the first timesteps 
        //    /// are computed with a smaller timestepsize
        //    IncrementInit,

        //    /// Initialization for a multi-step method, e.g. BDF4 requires 4 timesteps.
        //    /// can be used if an analytic solution is known or simulation is restarted form another session
        //    MultiInit
        //}

        /// <summary>
        /// See <see cref="TimestepperInit"/>
        /// </summary>
        [DataMember]
        public TimeStepperInit Timestepper_BDFinit = TimeStepperInit.SingleInit;

        /// <summary>
        /// defines the number of incremental timesteps in one gloabl timestep (for incrementInit)
        /// </summary>
        [DataMember]
        public int incrementTimesteps = 1;

        /// <summary>
        /// See <see cref="LevelSetHandling"/>
        /// </summary>
        [DataMember]
        public LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

        /// <summary>
        /// underrelaxation of the level set movement in case of coupled iterative
        /// </summary>
        [DataMember]
        public double LSunderrelax = 1.0;

        /// <summary>
        /// array of additional parameter values for some testcases
        /// </summary>
        [DataMember]
        public double[] AdditionalParameters;

        ///// <summary>
        ///// If iterative saddle-point solvers like GMRES or Orthonormalization are used, the maximum number of basis vectors
        ///// that are used to construct the accelerated solution.
        ///// </summary>
        //public int Solver_MaxKrylovDim = 100;

        ///// <summary>
        ///// If iterative saddle-point solvers are used, the termination criterion. 
        ///// </summary>
        //[DataMember]
        //public double Solver_ConvergenceCriterion = 1.0e-10;

        /// <summary>
        /// The termination criterion for fully coupled/implicit level-set evolution.
        /// </summary>
        [DataMember]
        public double LevelSet_ConvergenceCriterion = 1.0e-6;

        ///// <summary>
        ///// If iterative solvers are used, the maximum number of iterations.
        ///// </summary>
        //[DataMember]
        //public int Solver_MaxIterations = 2000;

        ///// <summary>
        ///// If iterative solvers are used, the minimum number of iterations.
        ///// </summary>
        //[DataMember]
        //public int Solver_MinIterations = 4;

        /// <summary>
        /// Block-Preconditiond for the velocity/momentum-block of the saddle-point system
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;

        /// <summary>
        /// Block-Preconditiond for the pressure/continuity-block of the saddle-point system
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode PressureBlockPrecondMode = MultigridOperator.Mode.IdMass;

        /// <summary>
        /// Block-Preconditiond for the stresses/constitutive-block of the system
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode StressBlockPrecondMode = MultigridOperator.Mode.Eye;

        /// <summary>
        /// Block-Preconditiond for the stresses/constitutive-block of the system
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode VelocityGradientBlockPrecondMode = MultigridOperator.Mode.Eye;


        /// <summary>
        /// If true, kinetic and surface energy will be evaluated in every cycle.
        /// </summary>
        [DataMember]
        public bool ComputeEnergy = false;

        /// <summary>
        /// If true, energy balance at the interface will be evaluated in every cycle.
        /// </summary>
        [DataMember]
        public bool ComputeInterfaceEnergy = false;

        /// <summary>
        /// Exact solution for GravityX source.
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> GravityX;
        /// <summary>
        /// Exact solution for GravityY source.
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> GravityY;
        /// <summary>
        /// Exact solution for GravityXX source.
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> GravityXX;
        /// <summary>
        /// Exact solution for GravityXY source.
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> GravityXY;
        /// <summary>
        /// Exact solution for GravityYY source.
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> GravityYY;

        /// <summary>
        /// Turn XDG for the velocity on/off; if off, only the pressure is approximated by XDG,
        /// the velocity is plain DG.
        /// </summary>
        //public bool UseXDG4Velocity = true;

        /// <summary>
        /// Exact solution for velocity, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>[]> ExactSolutionVelocity;

        /// <summary>
        /// Exact solution, pressure, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> ExactSolutionPressure;

        /// <summary>
        /// Exact solution for stresses, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> ExactSolutionStressXX;
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> ExactSolutionStressXY;
        [NonSerialized]
        public IDictionary<string, Func<double[], double, double>> ExactSolutionStressYY;

        /// <summary>
        /// Control Options for ReInit
        /// </summary>
        public EllipticReInitAlgoControl ReInitControl = new EllipticReInitAlgoControl();


        //RHEOLOGY PARAMETERS!!!
        //====================================================================================

        /// <summary>
        /// Use artificial Diffusion for smoothing singularities in stresses
        /// </summary>
        [DataMember]
        public bool UseArtificialDiffusion = false;

        /// <summary>
        /// Use finite differences Jacobian for Linearization
        /// </summary>
        [DataMember]
        public bool useJacobianForOperatorMatrix = false;

        /// <summary>
        /// adds a gravity source to the RHS
        /// </summary>
        [DataMember]
        public bool GravitySource = false;

        /// <summary>
        /// Raise Weissenberg Number?
        /// </summary>
        [DataMember]
        public bool RaiseWeissenberg = false;

        /// <summary>
        /// which increment?
        /// </summary>
        [DataMember]
        public double WeissenbergIncrement = 0.1;

        /// <summary>
        /// Use Persson Sensor to detect high energy modes of singularities
        /// </summary>
        [DataMember]
        public bool UsePerssonSensor = false;

        /// <summary>
        /// bound for perssonsensor should be around 1e-7 - 1e-8 that there is refinement or art. diffusion behind the cylinder!
        /// </summary>
        [DataMember]
        public double SensorLimit = 1e-7;

        /// <summary>
        /// Analysis of Operator Matrix (rank, cond...)?
        /// </summary>
        [DataMember]
        public bool OperatorMatrixAnalysis = false;

        /// <summary>
        /// Analytical solution for linearized problem A(u_ex) * u_new = b(u_ex)?
        /// </summary>
        [DataMember]
        public bool SetParamsAnalyticalSol = false;

        /// <summary>
        /// default velocity U function
        /// </summary>
        [NonSerialized]
        public Func<double[], double> VelFunctionU = X => 0;

        /// <summary>
        /// default velocity V function
        /// </summary>
        [NonSerialized]
        public Func<double[], double> VelFunctionV = X => 0;

        /// <summary>
        /// Exact solution for stresses.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double>[] ExSol_Stress;
    }
}









//// RHEOLOGY SOLVER OPTIONS
//// =============================================
///// <summary>
///// insert initial conditions
///// </summary>
//[DataMember]
//public bool SetInitialConditions = true;


///// <summary>
///// Convergence criterion stresses
///// </summary>
//[DataMember]
//public double ConvCritStress = 1E-10;


////=============================================================
//// NOT QUITE SURE WHAT THEY ARE FOR??????
///// <summary>
///// Convergence criterion VelocitySolver
///// </summary>
//[DataMember]
//public double VelocitySolver_ConvergenceCriterion = 1e-5;

///// <summary>
///// Convergence criterion StressSolver
///// </summary>
//[DataMember]
//public double StressSolver_ConvergenceCriterion = 1e-5;

///// <summary>
///// should velocity be solved
///// </summary>
//[DataMember]
//public bool solveVelocity = true;
////===========================================================






///// <summary>
///// periodic BC?
///// </summary>
//[DataMember]
//public bool FixedStreamwisePeriodicBC = false;

///// <summary>
///// Fixed SrcPressureGradient which should be used if periodic BC are applied
///// </summary>
//[DataMember]
//public double[] SrcPressureGrad = new double[] { -1, 0 };

///// <summary>
///// penalty factor in viscous part (SIP)
///// </summary>
//[DataMember]
//public double ViscousPenaltyScaling = 1;




////DEBUGGING PARAMETERS
////_____________________________________________________________________________________________



///// <summary>
///// Compute L2 Error of exact solution?
///// </summary>
//[DataMember]
//public bool ComputeL2Error = false;

///// <summary>
///// Analysis Level Operator Matrix
///// </summary>
//[DataMember]
//public int AnalysisLevel = 2;

///// <summary>
///// solver is turned off and residual of initial value/exact solution is evaluated, used to test the consistency of the implementation. We need an initial condition for pressure.
///// </summary>
//[DataMember]
//public bool SkipSolveAndEvaluateResidual = false;

///// <summary>
///// Initial pressure?
///// </summary>
//[DataMember]
//public bool SetInitialPressure = false;

///// <summary>
///// Analytical solution for linearized problem A(u_ex) * u_new = b(u_ex)?
///// </summary>
//public bool SetParamsAnalyticalSol = false;

///// <summary>
///// default velocity U function
///// </summary>
//[NonSerialized]
//public Func<double[], double> VelFunctionU = X => 0;

///// <summary>
///// default velocity V function
///// </summary>
//[NonSerialized]
//public Func<double[], double> VelFunctionV = X => 0;

///// <summary>
///// default pressure function
///// </summary>
//[NonSerialized]
//public Func<double[], double> PresFunction = X => 0;

///// <summary>
///// Exact solution for velocity.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double>[] ExSol_Velocity;


///// <summary>
///// Exact solution, pressure, for each species (either A or B).
///// </summary>
//[NonSerialized]
//public Func<double[], double, double> ExSol_Pressure;

///// <summary>
///// Exact solution for stresses.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double>[] ExSol_Stress;

///// <summary>
///// Exact solution for GravityX source.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double> GravityX;
///// <summary>
///// Exact solution for GravityY source.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double> GravityY;
///// <summary>
///// Exact solution for GravityXX source.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double> GravityXX;
///// <summary>
///// Exact solution for GravityXY source.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double> GravityXY;
///// <summary>
///// Exact solution for GravityYY source.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double> GravityYY;
///// <summary>
///// Exact solution for GravityDiv source.
///// </summary>
//[NonSerialized]
//public Func<double[], double, double> GravityDiv;

///// <summary>
///// polynomial degree for Unit Testing
///// </summary>
//[DataMember]
//public int deg;
///// <summary>
///// Grid resolution
///// </summary>
//[DataMember]
//public int grd;


//        //_____________________________________________________________________________________________
//    }
//}
