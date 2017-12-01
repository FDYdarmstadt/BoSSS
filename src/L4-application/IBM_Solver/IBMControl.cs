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
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.XdgTimestepping;
using System.Linq;

namespace BoSSS.Application.IBM_Solver {

    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class IBM_Control : AppControl {

        /// <summary>
        /// Ctor.
        /// </summary>
        public IBM_Control() {
            base.NoOfMultigridLevels = 1;
        }

        ///// <summary>
        ///// Expert options regarding the level set solver.
        ///// </summary>
        //public XVelocityProjection.Configuration LevelSetOptions = new XVelocityProjection.Configuration() {
        //};
        //public int LS_TrackerWidth = 1;

        /// <summary>
        /// Expert options regarding the spatial discretization.
        /// </summary>
        public DoNotTouchParameters AdvancedDiscretizationOptions = new DoNotTouchParameters();


        
        /// <summary>
        /// If iterative saddle-point solvers like GMRES or Orthonormalization are used, the maximum number of basis vectors
        /// that are used to construct the accelerated solution.
        /// </summary>
        public int MaxKrylovDim = 100;

        /// <summary>
        /// If iterative solvers are used, the maximum number of iterations.
        /// </summary>
        public int MaxSolverIterations = 2000;

        /// <summary>
        /// If iterative solvers are used, the maximum number of iterations.
        /// </summary>
        public int MinSolverIterations = 2;

        /// <summary>
        /// Convergence criterion for linear/nonlinear solver.
        /// </summary>
        public double Solver_ConvergenceCriterion = 1.0e-8;

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
                    if(ExSol_Velocity!= null) {
                        m_ExSol_Velocity_Evaluator = ExSol_Velocity.Select< IBoundaryAndInitialData, Func<double[], double, double>>( a => a.Evaluate).ToArray();
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
                if(m_ExSol_Pressure_Evaluator == null) {
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
        public PhysicalParameters PhysicalParameters = new PhysicalParameters()
        {
            IncludeConvection = true,
            rho_A = 1,
            mu_A = 1,
        };

        /// <summary>
        /// Radius of the circular particle immersed in the fluid
        /// </summary>
        public double particleRadius;

        public double MeshFactor;



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
        public TimesteppingScheme Timestepper_Scheme;

        /// <summary>
        /// Set true if periodic boundary conditions in streamwise direction are applied
        /// </summary>
        public bool FixedStreamwisePeriodicBC = false;

        /// <summary>
        /// Fixed SrcPressureGradient which should be used if periodic BC are applied
        /// </summary>
        public double[] SrcPressureGrad;


        /// <summary>
        /// Which direct solver should be used.
        /// </summary>
        public DirectSolver._whichSolver whichSolver = DirectSolver._whichSolver.PARDISO;


        public enum TimestepperInit {

            SingleInit,

            IncrementInit,

            MultiInit
        }

        public ISolverSmootherTemplate LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };

        public NonlinearSolverMethod NonlinearMethod = NonlinearSolverMethod.Picard;

        /// <summary>
        /// 
        /// </summary>
        public TimestepperInit Timestepper_Init = TimestepperInit.SingleInit;
    }
}
