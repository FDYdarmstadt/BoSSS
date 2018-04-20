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

            // Set to pseudo Picard if the Stokes equations should be solved
            if (Control.PhysicalParameters.IncludeConvection == false)
                Control.NonlinearSolve = NonlinearSolverCodes.Picard;

            // Set nonlinear Solver
            switch (Control.NonlinearSolve) {
                case NonlinearSolverCodes.NewtonGMRES:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.NewtonGMRES;
                    break;
                case NonlinearSolverCodes.Picard:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.Picard;
                    break;
                case NonlinearSolverCodes.Newton:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.Newton;
                    break;
                default:
                    throw new NotImplementedException("Nonlinear solver option not available");
            }

            switch (Control.LinearSolve) {
                case LinearSolverCodes.automatic:
                    AutomaticChoice(Control, Timestepper);
                    break;

                case LinearSolverCodes.classic_mumps:
                    Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };
                    break;

                case LinearSolverCodes.classic_pardiso:
                    Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };
                    break;

                case LinearSolverCodes.exp_schwarz_directcoarse_overlap:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    Timestepper.Config_linearSolver = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverCodes.exp_schwarz_directcoarse:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    Timestepper.Config_linearSolver = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverCodes.exp_schwarz_Kcycle_directcoarse:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    Timestepper.Config_linearSolver = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = Control.NoOfMultigridLevels - 1
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverCodes.exp_schwarz_Kcycle_directcoarse_overlap:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    Timestepper.Config_linearSolver = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = Control.NoOfMultigridLevels - 1
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverCodes.exp_softgmres:
                    Timestepper.Config_linearSolver = new SoftGMRES() {
                        MaxKrylovDim = Timestepper.Config_MaxKrylovDim,
                        m_Tolerance = Timestepper.Config_SolverConvergenceCriterion,
                    };
                    break;

                case LinearSolverCodes.exp_softgmres_schwarz_Kcycle_directcoarse_overlap:
                    Timestepper.Config_linearSolver = new SoftGMRES() {
                        MaxKrylovDim = Timestepper.Config_MaxKrylovDim,
                        m_Tolerance = Timestepper.Config_SolverConvergenceCriterion,
                        Precond = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                Depth = Control.NoOfMultigridLevels - 1
                            },
                            Overlap = 1,
                            CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                        },
                    };
                    break;

                case LinearSolverCodes.exp_softgmres_schwarz_directcoarse_overlap:
                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");
                    Timestepper.Config_linearSolver = new SoftGMRES() {
                        MaxKrylovDim = Timestepper.Config_MaxKrylovDim,
                        m_Tolerance = Timestepper.Config_SolverConvergenceCriterion,
                        Precond = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsPerProcess = 1,
                            },
                            Overlap = 1,
                            CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                        },
                    };
                    break;

                default:
                    throw new NotImplementedException("Linear solver option not available");
            }
        }

        /// <summary>
        /// Automatic choice of linear solver depending on problem size, immersed boundary, polynomial degree, etc.
        /// </summary>
        static void AutomaticChoice(IBM_Control Control, XdgBDFTimestepping Timestepper) {

            //int pV = Control.FieldOptions["VelocityX"].Degree;
            int pP = Control.FieldOptions["Pressure"].Degree;
            int pV = pP + 1;

            // Detecting variables for solver determination 
            var D = Timestepper.MultigridSequence[0].SpatialDimension;
            var cellsLoc = Timestepper.MultigridSequence[0].CellPartitioning.LocalLength;
            var cellsGlo = Timestepper.MultigridSequence[0].CellPartitioning.TotalLength;


            var size = Timestepper.MultigridSequence[0].CellPartitioning.MpiSize;

            // !!!!!!!!!!!UNTERSCHEIDUNG OB PICARD ODER NEWTON!!!!!!!!!!!!
            if (Timestepper.Config_NonlinearSolver == NonlinearSolverMethod.NewtonGMRES) {

                // Spatial Dimension
                switch (D) {
                    case 1:
                        break;
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                        break;

                    case 2:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                        break;

                    case 3:
                        var dofsPerCell3D = (3 * (pV * pV * pV + 6 * pV * pV + 11 * pV + 6) / 6 + 1 * (pP * pP * pP + 6 * pP * pP + 11 * pP + 6) / 6);
                        var dofsLoc = dofsPerCell3D * cellsLoc;
                        var dofsGlo = dofsPerCell3D * cellsGlo;

                        var PPP = (int)Math.Ceiling(dofsLoc / 6500.0);

                        Console.WriteLine("Analysing the problem yields " + PPP + " parts per process.");

                        if (dofsGlo > 10000) {

                            if (Control.NoOfMultigridLevels < 2)
                                throw new ApplicationException("At least 2 Multigridlevels are required");

                            Timestepper.Config_linearSolver = new Schwarz() {
                                m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                    NoOfPartsPerProcess = PPP,
                                },
                                Overlap = 1,
                                CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                            };
                        } else {
                            Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };
                        }
                        break;

                    default:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                }
            } else {
                // Spatial Dimension
                switch (D) {
                    case 1:
                        break;
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                        break;

                    case 2:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                        break;

                    case 3:
                        var dofsPerCell3D = (3 * (pV * pV * pV + 6 * pV * pV + 11 * pV + 6) / 6 + 1 * (pP * pP * pP + 6 * pP * pP + 11 * pP + 6) / 6);
                        var dofsLoc = dofsPerCell3D * cellsLoc;
                        var dofsGlo = dofsPerCell3D * cellsGlo;

                        if (dofsGlo > 10000) {

                            if (Control.NoOfMultigridLevels < 2)
                                throw new ApplicationException("At least 2 Multigridlevels are required");

                            Timestepper.Config_linearSolver = new SoftGMRES() {
                                MaxKrylovDim = Timestepper.Config_MaxKrylovDim,
                                m_Tolerance = Timestepper.Config_SolverConvergenceCriterion,
                                Precond = new Schwarz() {
                                    m_BlockingStrategy = new Schwarz.SimpleBlocking() {
                                        NoOfPartsPerProcess = (int)Math.Ceiling(dofsLoc / 6500.0),
                                    },
                                    Overlap = 1,
                                    CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                                },
                            };
                        } else {
                            Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };
                        }
                        break;

                    default:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                }
            }

            //Timestepper.

            // Wenn Gesamtproblem in 2D < 100000 DoFs -> Direct Solver
            // Wenn Gesamtproblem in 3D < 10000 DoFs -> Direct Solver 

            // Block Solve 3D ca. 6000 DoFs per Process -> Adjust Blocks per Process
            // Coarse Solve ca. 5000 bis 10000 DoFs. -> Adjust Multigrid Levels

        }

        /// <summary>
        /// Determines a solver sequence depending on MGlevels
        /// </summary>
        /// <param name="MGlevels"></param>
        /// <param name="CoarsestSolver"></param>
        /// <returns></returns>
        static ISolverSmootherTemplate DetermineMGSquence(int MGlevels) {
            ISolverSmootherTemplate solver;
            if (MGlevels > 0) {
                solver = new ClassicMultigrid() { CoarserLevelSolver = DetermineMGSquence(MGlevels - 1) };
            } else {
                solver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };
            }
            return solver;
        }
    }
}
