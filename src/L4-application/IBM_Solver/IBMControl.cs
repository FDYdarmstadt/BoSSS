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
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.IBM_Solver {

    /// <summary>
    /// 
    /// </summary>
    [DataContract]
    [Serializable]
    public class IBM_Control : AppControl {

        /// <summary>
        /// Ctor.
        /// </summary>
        public IBM_Control() {
            base.LinearSolver.NoOfMultigridLevels = 1;
            //shift of Solver Information
            base.LinearSolver.MaxKrylovDim = 30; //MaxKrylovDim
            base.LinearSolver.MaxSolverIterations = 2000; //MaxSolverIterations
            base.LinearSolver.MinSolverIterations = 2; //MinSolverIterations
            base.LinearSolver.ConvergenceCriterion = 1.0e-8; //Solver_ConvergenceCriterion
            base.LinearSolver.SolverCode = LinearSolverConfig.Code.classic_mumps; //public LinearSolverCodes LinearSolve = LinearSolverCodes.classic_mumps;
            base.NonLinearSolver.MaxSolverIterations = 2000; //MaxSolverIterations
            base.NonLinearSolver.MinSolverIterations = 2; //MinSolverIterations
            base.NonLinearSolver.ConvergenceCriterion = 1.0e-8; //Solver_ConvergenceCriterion
            base.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Picard; //public NonlinearSolverCodes NonlinearSolve = NonlinearSolverCodes.Picard;
        }

        /// <summary>
        /// Type of <see cref="IBM_SolverMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(IBM_SolverMain);
        }

        /// <summary>
        /// Expert options regarding the spatial discretization.
        /// </summary>
        [DataMember]
        public DoNotTouchParameters AdvancedDiscretizationOptions = new DoNotTouchParameters();

        /// <summary>
        /// desired minimum refinement level, 2 is minimum
        /// </summary>
        [DataMember]
        public int RefinementLevel = 2;

        /// <summary>
        /// reciprocal of the ratio between curvature and hmin
        /// </summary>
        [DataMember]
        public int maxCurvature = 2;

        /// <summary>
        /// Sets the DG polynomial degree 
        /// </summary>
        /// <param name="k">Degree for velocity; pressure  will be one order lower.</param>
        public override void SetDGdegree(int k) {
            if (k < 1)
                throw new ArgumentOutOfRangeException("DG polynomial degree must be at least 1.");

            base.FieldOptions.Clear();
            this.AddFieldOption("Velocity*", k);
            this.AddFieldOption("Pressure", k - 1);
            this.AddFieldOption("PhiDG", 2);
            this.AddFieldOption("Phi", 2);
        }

        /// <summary>
        /// Block-Preconditiond for the velocity-components of the saddel-point system
        /// </summary>
        public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;

        /// <summary>
        /// Enforce the level-set to be continuous; 
        /// </summary>
        public bool LevelSetSmoothing = false;

        [NonSerialized]
        Func<double[], double, double>[] m_ExSol_Velocity_Evaluator = null;

        /// <summary>
        /// Exact solution for velocity, in order to compute e.g. convergence rates.
        /// This property is non-serializeable, it cannot be used with control object files; if this is required, use 
        /// <see cref="ExSol_Velocity"/>.
        /// </summary>
        public Func<double[], double, double>[] ExSol_Velocity_Evaluator {
            get {
                if (m_ExSol_Velocity_Evaluator == null) {
                    if (ExSol_Velocity != null) {
                        m_ExSol_Velocity_Evaluator = ExSol_Velocity.Select<IBoundaryAndInitialData, Func<double[], double, double>>(a => a.Evaluate).ToArray();
                    }
                }
                return m_ExSol_Velocity_Evaluator;
            }
            set {
                ExSol_Velocity = null;
                m_ExSol_Velocity_Evaluator = value;
            }
        }

        /// <summary>
        /// Exact solution for velocity, in order to compute e.g. convergence rates.
        /// </summary>
        public IBoundaryAndInitialData[] ExSol_Velocity = null;


        [NonSerialized]
        Func<double[], double, double> m_ExSol_Pressure_Evaluator = null;


        /// <summary>
        /// Exact solution for pressure, in order to compute e.g. convergence rates.
        /// This property is non-serializeable, it cannot be used with control object files; if this is required, use 
        /// <see cref="ExSol_Pressure"/>.
        /// </summary>
        public Func<double[], double, double> ExSol_Pressure_Evaluator {
            get {
                if (m_ExSol_Pressure_Evaluator == null) {
                    if (ExSol_Pressure != null)
                        m_ExSol_Pressure_Evaluator = ExSol_Pressure.Evaluate;
                }
                return m_ExSol_Pressure_Evaluator;
            }
            set {
                ExSol_Pressure = null;
                m_ExSol_Pressure_Evaluator = value;
            }
        }

        /// <summary>
        /// Exact solution for pressure, in order to compute e.g. convergence rates.
        /// </summary>
        public IBoundaryAndInitialData ExSol_Pressure = null;

        /// <summary>
        /// Viscosity, density and surface tension. Note phase A is fluid and phase B particle.
        /// </summary>
        [DataMember]
        public PhysicalParameters PhysicalParameters = new PhysicalParameters() {
            IncludeConvection = true,
            rho_A = 1,
            mu_A = 1,
        };

        /// <summary>
        /// Radius of the circular particle immersed in the fluid
        /// </summary>
        [DataMember]
        public double particleRadius;


        /// <summary>
        /// 
        /// </summary>
        public enum TimesteppingScheme {


            ImplicitEuler = 1,

            CrankNicolson = 2,

            BDF2 = 3,

            BDF3 = 4,

            BDF4 = 5,

            BDF5 = 6,

            BDF6 = 7


        }

        /// <summary>
        /// See <see cref="TimesteppingScheme"/>
        /// </summary>
        [DataMember]
        public TimesteppingScheme Timestepper_Scheme = TimesteppingScheme.BDF2;

        /// <summary>
        /// Set true if periodic boundary conditions in streamwise direction are applied
        /// </summary>
        public bool FixedStreamwisePeriodicBC = false;

        /// <summary>
        /// Fixed SrcPressureGradient which should be used if periodic BC are applied
        /// </summary>
        [DataMember]
        public double[] SrcPressureGrad;

        public enum TimestepperInit {

            SingleInit,

            IncrementInit,

            MultiInit
        }

        [DataMember]
        public BoSSS.Solution.Timestepping.TimeStepperInit TimeStepper_Init = Solution.Timestepping.TimeStepperInit.SingleInit;
    }
}
