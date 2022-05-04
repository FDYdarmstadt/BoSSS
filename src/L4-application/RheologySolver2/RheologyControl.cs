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
using ilPSP.Utils;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Control File For calculation with viscoelastic extra stress tensor
    /// </summary>
    [Serializable]
    [DataContract]
    public class RheologyControl : AppControlSolver {

        /// <summary>
        /// Ctor.
        /// </summary>
        public RheologyControl() {
            base.NoOfMultigridLevels = 1;
            //shift of Solver Information
            base.LinearSolver = LinearSolverCode.direct_mumps.GetConfig(); //LinearSolver
            base.NonLinearSolver.MaxSolverIterations = 10; //MaxIter
            base.NonLinearSolver.MinSolverIterations = 1; //MinIter
            base.NonLinearSolver.ConvergenceCriterion = 1.0e-10; //ConvCrit
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton; //NonLinearSolver
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
        /// Reynolds-Number of System (Re= rho * U * L / eta_0)
        /// </summary>
        [DataMember]
        public double Reynolds = 1;

        /// <summary>
        /// Weissenberg-Number of System (We= lambda_1 * U / L) 
        /// </summary>
        [DataMember]
        public double Weissenberg = 0.5;

        /// <summary>
        /// Retardation vs Relaxation ratio (beta = lambda_2 / lambda_1 = eta_s / eta_0)
        /// - beta = 0: no Newtonian contribution in momentum equation
        /// - beta = 1: Newtonian fluid, all viscous effects computed in momentum equation
        /// </summary>
        [DataMember]
        public double beta = 0.11;

        /// <summary>
        /// Upwinding factor for convective part in constitutive equation
        /// </summary>
        [DataMember]
        public double alpha = 0.5;

        /// <summary>
        /// determines which implementation of objective Term should be used:
        /// =1 only velocity gradient as param
        /// =0 only stress tensor as param
        /// </summary>
        [DataMember]
        public double ObjectiveParam = 1.0;


        /// <summary>
        /// Giesekus model factor with which Giesekus term should be multiplied, if factor = 0 we have Oldroyd B
        /// </summary>
        [DataMember]
        public double giesekusfactor = 0.0;

        //_____________________________________________________________________________________________

        //SOLVING SYSTEM
        //_____________________________________________________________________________________________

        /// <summary>
        /// Stokes (true) or Navier-Stokes (false) flow?
        /// </summary>
        [DataMember]
        public bool Stokes = false;

        /// <summary>
        /// Stokes (true), but iterative for convection in constitutive
        /// </summary>
        [DataMember]
        public bool StokesConvection = false;

        /// <summary>
        /// insert initial conditions
        /// </summary>
        [DataMember]
        public bool SetInitialConditions = true;

        /// <summary>
        /// adds a gravity source to the RHS
        /// </summary>
        [DataMember]
        public bool GravitySource = false;

        // <summary>
        // updating algorithm for u
        // </summary>
        //[DataMember]
        //public bool UpdateUAlg = false;

        /// <summary>
        /// Convergence criterion stresses
        /// </summary>
        [DataMember]
        public double ConvCritStress = 1E-10;

        /// <summary>
        /// Switches Homotopy solver on/off:
        /// <see cref="Newton.UseHomotopy"/>, <see cref="BoSSS.Foundation.ISpatialOperator.HomotopyUpdate"/>
        /// </summary>
        [DataMember]
        public bool RaiseWeissenberg = false;

        
        /// <summary>
        /// Use finite differences Jacobian for Linearization
        /// </summary>
        [DataMember]
        public bool useFDJacobianForOperatorMatrix = false;

        /// <summary>
        /// periodic BC?
        /// </summary>
        [DataMember]
        public bool FixedStreamwisePeriodicBC = false;

        /// <summary>
        /// Fixed SrcPressureGradient which should be used if periodic BC are applied
        /// </summary>
        [DataMember]
        public double[] SrcPressureGrad = new double[] { -1, 0 };

        /// <summary>
        /// penalty factor in viscous part (SIP)
        /// </summary>
        [DataMember]
        public double ViscousPenaltyScaling = 1;

        /// <summary>
        /// Penalty Values LDG (alpha, beta; Lit. Cockburn (2002) Local DG Methods for the Stokes system)
        /// Penalty in Stress Divergence (beta)
        /// </summary>
        [DataMember]
        public double[] Penalty1 = { 0, 0 };

        /// <summary>
        /// Penalty in Constitutive Viscosity (alpha)
        /// </summary>
        [DataMember]
        public double Penalty2 = 1.0;

        /// <summary>
        /// Penalty for pressure/conti (beta)
        /// </summary>
        [DataMember]
        public double[] PresPenalty1 = { 0, 0 };

        /// <summary>
        /// Penalty for pressure (alpha)
        /// </summary>
        [DataMember]
        public double PresPenalty2 = 1.0;

        /// <summary>
        /// penalty for stress in objective term
        /// </summary>
        [DataMember]
        public double StressPenalty = 1.0;

        ///// <summary>
        ///// Block-Preconditiond for the velocity/momentum-block of the saddle-point system
        ///// </summary>
        //[DataMember]
        //public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.Eye;// LeftInverse_DiagBlock;// LeftInverse_DiagBlock;  // SymPart_DiagBlockEquilib_DropIndefinite;

        ///// <summary>
        ///// Block-Preconditiond for the velocity/momentum + pressure/conti-block of the saddle-point system
        ///// </summary>
        //[DataMember]
        //public MultigridOperator.Mode NSEBlockPrecondMode = MultigridOperator.Mode.Eye;// LeftInverse_DiagBlock; //.LeftInverse_DiagBlock; // SymPart_DiagBlockEquilib_DropIndefinite;

        ///// <summary>
        ///// Block-Preconditiond for the pressure/continuity-block of the saddle-point system
        ///// </summary>
        //[DataMember]
        //public MultigridOperator.Mode PressureBlockPrecondMode = MultigridOperator.Mode.Eye;// SymPart_DiagBlockEquilib;//LeftInverse_DiagBlock; // no SymPart_Diag-Präcon, because there may be no zero on the diagonal!!!

        ///// <summary>
        ///// Block-Preconditiond for the stresses/constitutive-block of the system
        ///// </summary>
        //[DataMember]
        //public MultigridOperator.Mode StressBlockPrecondMode = MultigridOperator.Mode.Eye;// SymPart_DiagBlockEquilib;//LeftInverse_DiagBlock;

        ///// <summary>
        ///// Block-Preconditiond for the stresses/constitutive-block of the system
        ///// </summary>
        //[DataMember]
        //public MultigridOperator.Mode VelocityGradientBlockPrecondMode = MultigridOperator.Mode.Eye;

        /// <summary>
        /// Refinement level for adaptive mesh refinement
        /// </summary>
        [DataMember]
        public int RefinementLevel = 0;
        //_____________________________________________________________________________________________


        // TIMESTEPPING
        //_____________________________________________________________________________________________

        /// <summary>
        /// Timestep (default is large for steady calculation)
        /// </summary>
        [DataMember]
        public double dt = 1E20;

        /// <summary>
        /// Timestepping scheme
        /// </summary>
        public enum TimesteppingScheme {

            /// <summary>
            /// ImplicitEuler
            /// </summary>
            ImplicitEuler = 1,
            /// <summary>
            /// CrankNicolson
            /// </summary>
            CrankNicolson = 2,
            /// <summary>
            /// BDF2
            /// </summary>
            BDF2 = 3,
            /// <summary>
            /// BDF3
            /// </summary>
            BDF3 = 4,
            /// <summary>
            /// BDF4
            /// </summary>
            BDF4 = 5,
            /// <summary>
            /// BDF5
            /// </summary>
            BDF5 = 6,
            /// <summary>
            /// BDF6
            /// </summary>
            BDF6 = 7


        }
        /// <summary>
        /// Timestepper scheme
        /// </summary>
        [DataMember]
        public TimesteppingScheme Timestepper_Scheme;
        //_____________________________________________________________________________________________

        //DEBUGGING PARAMETERS
        //_____________________________________________________________________________________________

        ///// <summary>
        ///// Analysis of Operator Matrix (rank, cond...)?
        ///// </summary>
        //[DataMember]
        //public bool OperatorMatrixAnalysis = false;

        /// <summary>
        /// Compute L2 Error of exact solution?
        /// </summary>
        [DataMember]
        public bool ComputeL2Error = false;

        /// <summary>
        /// Compute body forces on wall?
        /// </summary>
        [DataMember]
        public bool Bodyforces = false;

        /// <summary>
        /// Analysis Level Operator Matrix
        /// </summary>
        [DataMember]
        public int AnalysisLevel = 2;

        /// <summary>
        /// solver is turned off and residual of initial value/exact solution is evaluated, used to test the consistency of the implementation. We need an initial condition for pressure.
        /// </summary>
        [DataMember]
        public bool SkipSolveAndEvaluateResidual = false;

        /// <summary>
        /// Initial pressure?
        /// </summary>
        [DataMember]
        public bool SetInitialPressure = false;

        /// <summary>
        /// Analytical solution for linearized problem A(u_ex) * u_new = b(u_ex)?
        /// </summary>
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
        /// default pressure function
        /// </summary>
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
        /// Exact solution for GravityX source.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> GravityX;
        /// <summary>
        /// Exact solution for GravityY source.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> GravityY;
        /// <summary>
        /// Exact solution for GravityXX source.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> GravityXX;
        /// <summary>
        /// Exact solution for GravityXY source.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> GravityXY;
        /// <summary>
        /// Exact solution for GravityYY source.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> GravityYY;
        /// <summary>
        /// Exact solution for GravityDiv source.
        /// </summary>
        [NonSerialized]
        public Func<double[], double, double> GravityDiv;

        /// <summary>
        /// polynomial degree for Unit Testing
        /// </summary>
        [DataMember]
        public int deg;
        
        /// <summary>
        /// Grid resolution
        /// </summary>
        [DataMember]
        public int grd;


        //_____________________________________________________________________________________________

        /// <summary>
        /// %
        /// </summary>
        public override void SetDGdegree(int degree) {
            FieldOptions.Clear();

            FieldOptions.Add("VelocityX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("VelocityY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("Pressure", new FieldOpts() { Degree = degree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            FieldOptions.Add("StressXX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("StressXY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("StressYY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            FieldOptions.Add("ResidualMomentumX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("ResidualMomentumY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("ResidualConti", new FieldOpts() { Degree = degree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            FieldOptions.Add("ResidualStressXX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("ResidualStressXY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("ResidualStressYY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            FieldOptions.Add("PhiDG", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("Phi", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        }

        /// <summary>
        /// Dummy override
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        /// <summary>
        /// Mainly for use by the job manager in BoSSSpad, in order to check if a specific configuration is already computed 
        /// (i.e. an equal control object can be found in the database),
        /// or not.
        /// Therefore, we only check attributes that considered `essential' properties (e.g. <see cref="Weissenberg"/>, not <see cref="AppControl.Tags"/>),
        /// in order to avoid a discrimination which is `too sharp'.
        /// </summary>
        public override bool Equals(object obj) {
            if(!base.Equals(obj))
                return false; // checks initial values, etc.

            RheologyControl oCtrl = obj as RheologyControl;
            if(oCtrl == null)
                return false;

            if(this.Weissenberg != oCtrl.Weissenberg)
                return false;

            if(this.Reynolds != oCtrl.Reynolds)
                return false;

            if(this.Stokes != oCtrl.Stokes)
                return false;

            //if(this.UsePerssonSensor != oCtrl.UsePerssonSensor)
            //    return false;

            //if(this.UsePerssonSensor != oCtrl.UsePerssonSensor)
            //    return false;
 
            if(this.alpha != oCtrl.alpha)
                return false;

            if(this.beta != oCtrl.beta)
                return false;

            if(this.GravitySource != oCtrl.GravitySource)
                return false;

            if(!this.Penalty1.ListEquals(oCtrl.Penalty1))
                return false;

            if(this.Penalty2 != oCtrl.Penalty2)
                return false;

            if(!this.PresPenalty1.ListEquals(oCtrl.PresPenalty1))
                return false;

            if(this.PresPenalty2 != oCtrl.PresPenalty2)
                return false;

            if(this.PresPenalty2 != oCtrl.PresPenalty2)
                return false;

                        if(this.PresPenalty2 != oCtrl.PresPenalty2)
                return false;


            return true;
        }

    }
}