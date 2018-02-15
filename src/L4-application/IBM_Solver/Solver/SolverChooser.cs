using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.Multigrid;

namespace BoSSS.Application.IBM_Solver {
    class SolverChooser {

        /// <summary>
        /// Choose solver depending on configurations made in the control file.
        /// </summary>
        /// <param name="nonlinSol"></param>
        /// <param name="linSol"></param>
        /// <param name="Timestepper"></param>
        public static void ChooseSolver(IBM_Control Control, ref XdgBDFTimestepping Timestepper) {

            // Set several solver options for Timestepper
            Timestepper.Config_SolverConvergenceCriterion = Control.Solver_ConvergenceCriterion;
            Timestepper.Config_MaxIterations = Control.MaxSolverIterations;
            Timestepper.Config_MinIterations = Control.MinSolverIterations;
            Timestepper.Config_MaxKrylovDim = Control.MaxKrylovDim;

            // Set nonlinear Solver
            switch (Control.NonlinearSolve) {
                case NonlinearSolverCodes.NewtonGMRES:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.Newton;
                    break;
                case NonlinearSolverCodes.Picard:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.Picard;
                    break;
                default:
                    throw new NotImplementedException("Nonlinear solver option not available");
            }


            switch (Control.LinearSolve) {
                case LinearSolverCodes.automatic:
                    AutomaticChoice();
                    break;

                case LinearSolverCodes.classic_mumps:
                    Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };
                    break;

                case LinearSolverCodes.classic_pardiso:
                    Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };
                    break;

                case LinearSolverCodes.exp_schwarz_directcoarse:
                    Timestepper.Config_linearSolver = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 1,
                        CoarseSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS }
                    };
                    break;

                case LinearSolverCodes.exp_schwarz_Kcycle_directcoarse:
                    Timestepper.Config_linearSolver = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = Control.NoOfMultigridLevels-1
                        },
                        Overlap = 1,
                        CoarseSolver = new ClassicMultigrid() {
                            CoarserLevelSolver= new ClassicMultigrid() {
                                CoarserLevelSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS}
                            }
                        }
                    };
                    break;

                case LinearSolverCodes.exp_softgmres:
                    Timestepper.Config_linearSolver = new SoftGMRES() {
                        MaxKrylovDim = Timestepper.Config_MaxKrylovDim,                       
                    };
                    break;

                default:
                    throw new NotImplementedException("Linear solver option not available");
            }
        }

        /// <summary>
        /// Automatic choice of linear solver depending on problem size, immersed boundary, polynomial degree, etc.
        /// </summary>
        static void AutomaticChoice() {
            throw new NotImplementedException("Option currently not available");
        }
    }
}
