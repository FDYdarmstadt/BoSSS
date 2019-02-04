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

using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using System.Runtime.Serialization;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Control File For calculation with extra stress tensor
    /// </summary>
    [Serializable]
    [DataContract]
    public class RheologyControl : AppControl {

        /// <summary>
        /// Ctor.
        /// </summary>
        public RheologyControl() {
            base.LinearSolver.NoOfMultigridLevels = 1;
            //shift of Solver Information
            base.LinearSolver.MaxSolverIterations = 10; //MaxIter
            base.LinearSolver.MinSolverIterations = 1; //MinIter
            base.LinearSolver.ConvergenceCriterion = 1.0e-6; //ConvCritGMRES
            base.LinearSolver.SolverCode = LinearSolverConfig.Code.classic_mumps; //LinearSolver
            base.NonLinearSolver.MaxSolverIterations = 10; //MaxIter
            base.NonLinearSolver.MinSolverIterations = 1; //MinIter
            base.NonLinearSolver.ConvergenceCriterion = 1.0e-10; //ConvCrit
            base.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Picard; //NonLinearSolver
            base.NonLinearSolver.UnderRelax = 1.0; //UnderRelax
        }

        /// <summary>
        /// 
        /// </summary>
        public override Type GetSolverType() {
            return typeof(Rheology);
        }

        // CONSTITUTIVE EQUATIONS
        //_____________________________________________________________________________________________
        // 

        /// <summary>
        /// Reynoldsnumber of System (Re= rho * U * L / eta_0)
        /// </summary>
        [DataMember]
        public double Reynolds = 1;

        /// <summary>
        /// Weissenbergnumber of System (We= lambda_1 * U / L) 
        /// </summary>
        [DataMember]
        public double Weissenberg = 0.5;

        //Retardation vs Relaxation ratio (beta = lambda_2 / lambda_1 = eta_s / eta_0)
        [DataMember]
        public double beta = 0.11;

        // Relaxation factor for convective part in constitutive equation
        [DataMember]
        public double alpha = 0.5;

        //determines which implementation of objective Term should be used:
        // =1 only velocity gradient as param
        // =0 only stress tensor as param
        [DataMember]
        public double ObjectiveParam = 0;

        //_____________________________________________________________________________________________

        //SOLVING SYSTEM
        //_____________________________________________________________________________________________
        // Stokes (true) or Navier-Stokes (false) flow?
        [DataMember]
        public bool Stokes = true;

        //insert initial conditions
        [DataMember]
        public bool SetInitialConditions = true;

        //adds a gravity source to the RHS
        [DataMember]
        public bool GravitySource = false;

        //updating algorithm for u 
        [DataMember]
        public bool UpdateUAlg = false;
        [DataMember]
        public double ConvCritStress = 1E-10;

        // Raise Weissenberg Number? which increment?
        [DataMember]
        public bool RaiseWeissenberg = false;
        [DataMember]
        public double WeissenbergIncrement = 0.1;

        //Use Persson Sensor to detect high energy modes of singularities
        [DataMember]
        public bool UsePerssonSensor = false;
        //bound for perssonsensor should be around 1e-7 - 1e-8 that there is refinement or art. diffusion behind the cylinder!
        [DataMember]
        public double SensorLimit = 1e-7;

        //Use artificial Diffusion for smoothing singularities in stresses
        [DataMember]
        public bool UseArtificialDiffusion = false;

        // Fixed SrcPressureGradient which should be used if periodic BC are applied
        [DataMember]
        public bool FixedStreamwisePeriodicBC = true;
        [DataMember]
        public double[] SrcPressureGrad = new double[] { -1, 0 };

        // penalty factor in viscous part (SIP)
        [DataMember]
        public double ViscousPenaltyScaling = 10;

        // Penalty Values LDG (alpha, beta; Lit. e.g. Arnold et al. Unified Analysis of DG methods)
        [DataMember]
        public double[] Penalty1 = { 1.0, 1.0 }; //Penalty in Stress Divergence (beta)
        [DataMember]
        public double Penalty2 = 1.0; //Penalty in Constitutive Viscosity (alpha)
        [DataMember]
        public double[] PresPenalty1 = { 1.0, 1.0 }; //Penalty for pressure/conti (beta)
        [DataMember]
        public double PresPenalty2 = 1.0; //Penalty for pressure (alpha)
        [DataMember]
        public double StressPenalty = 1.0; //penalty for stress in objective term

        ////Iterations for nonlinear solver (NS)
        //public int MaxIter = 10;
        //public int MinIter = 1;
        //public double ConvCrit = 1E-10;
        //public double ConvCritGMRES = 1E-6;
        //public double UnderRelax = 1.0;

        ///// <summary>
        ///// Which linear solver should be used.
        ///// </summary>
        //public ISolverSmootherTemplate LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };

        ////Which nonliner Solver Method in iteration (Fixpunkt = Picard or Newton)
        //public NonlinearSolverMethod NonlinearMethod = NonlinearSolverMethod.Picard;

        // Block-Preconditiond for the velocity/momentum-block of the saddle-point system
        [DataMember]
        public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite; //.LeftInverse_DiagBlock; // SymPart_DiagBlockEquilib;

        //Block-Preconditiond for the pressure/continuity-block of the saddle-point system
        [DataMember]
        public MultigridOperator.Mode PressureBlockPrecondMode = MultigridOperator.Mode.Eye; // no SymPart_Diag-Präcon, because there may be no zero on the diagonal!!!

        // Block-Preconditiond for the stresses/constitutive-block of the system
        [DataMember]
        public MultigridOperator.Mode StressBlockPrecondMode = MultigridOperator.Mode.Eye;

        // Block-Preconditiond for the stresses/constitutive-block of the system
        [DataMember]
        public MultigridOperator.Mode VelocityGradientBlockPrecondMode = MultigridOperator.Mode.Eye;

        ////Aggregation levels for multigrid
        //public int MultigridNoOfLevels = 0; wird nicht verwendet und ist obsolet, siehe Konstruktor ...

        //Refinement level for adaptive mesh refinement
        [DataMember]
        public int RefinementLevel = 0;
        //_____________________________________________________________________________________________


        // TIMESTEPPING
        //_____________________________________________________________________________________________

        // Timestep (default is large for steady calculation)
        [DataMember]
        public double dt = 1E20;

        public enum TimesteppingScheme {


            ImplicitEuler = 1,

            CrankNicolson = 2,

            BDF2 = 3,

            BDF3 = 4,

            BDF4 = 5,

            BDF5 = 6,

            BDF6 = 7


        }
        [DataMember]
        public TimesteppingScheme Timestepper_Scheme;
        //_____________________________________________________________________________________________

        //DEBUGGING PARAMETERS
        //_____________________________________________________________________________________________

        // Analysis of Operator Matrix (rank, cond...)
        [DataMember]
        public bool OperatorMatrixAnalysis = false;

        //Compute L2 Error of exact solution
        [DataMember]
        public bool ComputeL2Error = false;

        //Compute body forces on wall
        [DataMember]
        public bool Bodyforces = false;

        //Analysis Level Operator Matrix
        [DataMember]
        public int AnalysisLevel = 2;

        // solver is turned off and residual of initial value/exact solution is evaluated, used to 
        // test the consistency of the implementation. We need an initial condition for pressure.
        [DataMember]
        public bool SkipSolveAndEvaluateResidual = false;
        [DataMember]
        public bool SetInitialPressure = false;

        //Analytical solution for linearized problem A(u_ex) * u_new = b(u_ex)
        public bool SetParamsAnalyticalSol = false;

        [NonSerialized]
        public Func<double[], double> VelFunctionU = X => 0;

        [NonSerialized]
        public Func<double[], double> VelFunctionV = X => 0;

        [NonSerialized]
        public Func<double[], double> PresFunction = X => 0;

        /// <summary>
        /// Exact solution for velocity.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double>[] ExSol_Velocity;


        /// <summary>
        /// Exact solution, pressure, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> ExSol_Pressure;

        /// <summary>
        /// Exact solution for stresses.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double>[] ExSol_Stress;

        /// <summary>
        /// Exact solution for Gravity source.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> GravityX;
        [NonSerialized]
        public Func<double[], double, double> GravityY;
        [NonSerialized]
        public Func<double[], double, double> GravityXX;
        [NonSerialized]
        public Func<double[], double, double> GravityXY;
        [NonSerialized]
        public Func<double[], double, double> GravityYY;
        [NonSerialized]
        public Func<double[], double, double> GravityDiv;

        /// <summary>
        /// Grid resolution and polynomial degree for Unit Testing
        /// </summary>
        [DataMember]
        public int deg;
        [DataMember]
        public int grd;


        //_____________________________________________________________________________________________



    }
}