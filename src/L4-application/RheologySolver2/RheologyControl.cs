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
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Control File For calculation with extra stress tensor
    /// </summary>
    public class RheologyControl : AppControl {

        /// <summary>
        /// Ctor.
        /// </summary>
        public RheologyControl() {
            base.NoOfMultigridLevels = 1;
        }

        // CONSTITUTIVE EQUATIONS
        //_____________________________________________________________________________________________
        // Reynoldsnumber of System (Re= rho * U * L / eta_0)
        public double Reynolds = 1;

        // Weissenbergnumber of System (We= lambda_1 * U / L)
        public double Weissenberg = 0.5;

        //Retardation vs Relaxation ratio (beta = lambda_2 / lambda_1 = eta_s / eta_0)
        public double beta = 0.11;

        // Relaxation factor for convective part in constitutive equation
        public double alpha = 0.5;

        //determines which implementation of objective Term should be used:
        // =1 only velocity gradient as param
        // =0 only stress tensor as param
        public double ObjectiveParam = 0;

        //_____________________________________________________________________________________________

        //SOLVING SYSTEM
        //_____________________________________________________________________________________________
        // Stokes (true) or Navier-Stokes (false) flow?
        public bool Stokes = true;

        //insert initial conditions
        public bool SetInitialConditions = true;

        //adds a gravity source to the RHS
        public bool GravitySource = false;

        //updating algorithm for u 
        public bool UpdateUAlg = false;
        public double ConvCritStress = 1E-10;

        // Raise Weissenberg Number? which increment?
        public bool RaiseWeissenberg = false;
        public double WeissenbergIncrement = 0.1;

        //Use Persson Sensor to detect high energy modes of singularities
        public bool UsePerssonSensor = false;
        //bound for perssonsensor should be around 1e-7 - 1e-8 that there is refinement or art. diffusion behind the cylinder!
        public double SensorLimit = 1e-7;

        //Use artificial Diffusion for smoothing singularities in stresses
        public bool UseArtificialDiffusion = false;

        // Fixed SrcPressureGradient which should be used if periodic BC are applied
        public bool FixedStreamwisePeriodicBC = true;
        public double[] SrcPressureGrad = new double[] { -1, 0 };

        // penalty factor in viscous part (SIP)
        public double ViscousPenaltyScaling = 10;

        // Penalty Values LDG (alpha, beta; Lit. e.g. Arnold et al. Unified Analysis of DG methods)
        public double[] Penalty1 = { 1.0, 1.0 }; //Penalty in Stress Divergence (beta)
        public double Penalty2 = 1.0; //Penalty in Constitutive Viscosity (alpha)
        public double[] PresPenalty1 = { 1.0, 1.0 }; //Penalty for pressure/conti (beta)
        public double PresPenalty2 = 1.0; //Penalty for pressure (alpha)
        public double StressPenalty = 1.0; //penalty for stress in objective term

        //Iterations for nonlinear solver (NS)
        public int MaxIter = 10;
        public int MinIter = 1;
        public double ConvCrit = 1E-10;
        public double ConvCritGMRES = 1E-6;
        public double UnderRelax = 1.0;

        /// <summary>
        /// Which linear solver should be used.
        /// </summary>
        public ISolverSmootherTemplate LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };

        //Which nonliner Solver Method in iteration (Fixpunkt = Picard or Newton)
        public NonlinearSolverMethod NonlinearMethod = NonlinearSolverMethod.Picard;

        // Block-Preconditiond for the velocity/momentum-block of the saddle-point system
        public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite; //.LeftInverse_DiagBlock; // SymPart_DiagBlockEquilib;

        //Block-Preconditiond for the pressure/continuity-block of the saddle-point system
        public MultigridOperator.Mode PressureBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;

        // Block-Preconditiond for the stresses/constitutive-block of the system
        public MultigridOperator.Mode StressBlockPrecondMode = MultigridOperator.Mode.Eye;

        // Block-Preconditiond for the stresses/constitutive-block of the system
        public MultigridOperator.Mode VelocityGradientBlockPrecondMode = MultigridOperator.Mode.Eye;

        //Aggregation levels for multigrid
        public int MultigridNoOfLevels = 0;

        //Refinement level for adaptive mesh refinement
        public int RefinementLevel = 0;
        //_____________________________________________________________________________________________


        // TIMESTEPPING
        //_____________________________________________________________________________________________

        // Timestep (default is large for steady calculation)
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

        public TimesteppingScheme Timestepper_Scheme;
        //_____________________________________________________________________________________________

        //DEBUGGING PARAMETERS
        //_____________________________________________________________________________________________

        // Analysis of Operator Matrix (rank, cond...)
        public bool OperatorMatrixAnalysis = false;

        //Compute L2 Error of exact solution
        public bool ComputeL2Error = false;

        //Compute body forces on wall
        public bool Bodyforces = false;

        //Analysis Level Operator Matrix
        public int AnalysisLevel = 2;

        // solver is turned off and residual of initial value/exact solution is evaluated, used to 
        // test the consistency of the implementation. We need an initial condition for pressure.
        public bool SkipSolveAndEvaluateResidual = false;
        public bool SetInitialPressure = false;

        //Analytical solution for linearized problem A(u_ex) * u_new = b(u_ex)
        public bool SetParamsAnalyticalSol = false;
        public Func<double[], double> VelFunctionU = X => 0;
        public Func<double[], double> VelFunctionV = X => 0;
        public Func<double[], double> PresFunction = X => 0;

        /// <summary>
        /// Exact solution for velocity.
        /// </summary>
        public Func<double[], double, double>[] ExSol_Velocity;


        /// <summary>
        /// Exact solution, pressure, for each species (either A or B).
        /// </summary>
        public Func<double[], double, double> ExSol_Pressure;

        /// <summary>
        /// Exact solution for stresses.
        /// </summary>
        public Func<double[], double, double>[] ExSol_Stress;

        /// <summary>
        /// Exact solution for Gravity source.
        /// </summary>
        public Func<double[], double, double> GravityX;
        public Func<double[], double, double> GravityY;
        public Func<double[], double, double> GravityXX;
        public Func<double[], double, double> GravityXY;
        public Func<double[], double, double> GravityYY;
        public Func<double[], double, double> GravityDiv;

        /// <summary>
        /// Grid resolution and polynomial degree for Unit Testing
        /// </summary>
        public int deg;
        public int grd;


        //_____________________________________________________________________________________________



    }
}