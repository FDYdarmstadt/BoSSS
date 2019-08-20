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

namespace BoSSS.Application.XRheology_Solver {
    /// <summary>
    /// Control File For calculation with viscoelastic extra stress tensor
    /// </summary>
    [Serializable]
    [DataContract]
    public class XRheology_Control : XNSE_Control {

        /// <summary>
        /// Ctor.
        /// </summary>
        public XRheology_Control() {
            base.LinearSolver.NoOfMultigridLevels = 1;
            base.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            //shift of Solver Information
            base.LinearSolver.MaxSolverIterations = 10; //MaxIter
            base.LinearSolver.MinSolverIterations = 1; //MinIter
            base.LinearSolver.ConvergenceCriterion = 1.0e-6; //ConvCritGMRES
            base.LinearSolver.SolverCode = LinearSolverConfig.Code.classic_mumps; //LinearSolver
            base.NonLinearSolver.MaxSolverIterations = 10; //MaxIter
            base.NonLinearSolver.MinSolverIterations = 1; //MinIter
            base.NonLinearSolver.ConvergenceCriterion = 1.0e-10; //ConvCrit
            base.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Newton; //NonLinearSolver
            base.NonLinearSolver.UnderRelax = 1.0; //UnderRelax
        }

        /// <summary>
        /// Type of <see cref="XRheology_SolverMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XRheology_SolverMain);
        }

        public void SetFieldOptions(int VelDegree, int LevSetDegree, FieldOpts.SaveToDBOpt SaveFilteredVelocity = FieldOpts.SaveToDBOpt.TRUE, FieldOpts.SaveToDBOpt SaveCurvature = FieldOpts.SaveToDBOpt.TRUE) {
            FieldOptions.Add(VariableNames.VelocityX, new FieldOpts() { Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.VelocityY, new FieldOpts() { Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.Pressure, new FieldOpts() { Degree = VelDegree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.StressXX, new FieldOpts() { Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.StressXY, new FieldOpts() { Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.StressYY, new FieldOpts() { Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("PhiDG", new FieldOpts() { SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("Phi", new FieldOpts() { Degree = LevSetDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("Curvature", new FieldOpts() { Degree = LevSetDegree * 2, SaveToDB = SaveCurvature });
            FieldOptions.Add(VariableNames.Temperature, new FieldOpts() { Degree = VelDegree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        }


        [DataMember]
        public PhysicalParameters PhysicalParameters = new PhysicalParameters() {
            Material = true,
            IncludeConvection = false,
            reynolds_A = 1.0,
            reynolds_B = 1.0,
            rho_A = 1.0,
            rho_B = 1.0,
            Sigma = 0.0,
            Weissenberg_a =0.0,
            Weissenberg_b = 0.0,
            beta_a = 1.0,
            beta_b = 1.0,
        };

        /// <summary>
        /// Expert options regarding the spatial discretization.
        /// </summary>
        [DataMember]
        public DoNotTouchParameters AdvancedDiscretizationOptions = new DoNotTouchParameters() {
            PenaltySafety = 1.0,
            alpha = 1.0,
            ObjectiveParam =1.0,
            //Penalty1 = { 0, 0 },
            Penalty2 = 1.0,
            //PresPenalty1 = {0,0},
            PresPenalty2 = 1.0,
            StressPenalty = 1.0,
            ViscosityMode = ViscosityMode.Viscoelastic

        };

        // RHEOLOGY SOLVER OPTIONS
        // =============================================
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

        /// <summary>
        /// Convergence criterion stresses
        /// </summary>
        [DataMember]
        public double ConvCritStress = 1E-10;

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

        //BLOCK PRECONDITIONING
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


        //DEBUGGING PARAMETERS
        //_____________________________________________________________________________________________

        /// <summary>
        /// Analysis of Operator Matrix (rank, cond...)?
        /// </summary>
        [DataMember]
        public bool OperatorMatrixAnalysis = false;

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
    }
}
