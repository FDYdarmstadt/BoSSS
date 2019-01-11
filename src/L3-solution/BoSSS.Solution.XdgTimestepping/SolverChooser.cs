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
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.LinSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Solution.AdvancedSolvers {
//namespace BoSSS.Solution {

    public class SolverChooser {

        public static NonlinearSolver GenerateNonLin(Control.NonLinearSolverConfig nc) {

        }

        public static ISolverSmootherTemplate GenerateLinear(Control.LinearSolverConfig lc) {

        }



        /// <summary>
        /// Choose solver depending on configurations made in the control file. Information about linear- and NonlinearSolver configuration is extracted from controlfile. This SolverChooser is optimized for Navier-Stokes Problems and thus needs additonal parameters like convection and viscosity. ToDo: linConfig is overtaken for linear and Nonlinear Solver -> set config of Nonlinear Solver
        /// </summary>
        /// <param name="Control"></param>
        /// <param name="Timestepper"></param>
        /// <param name="Convection"></param>
        /// <param name="muA"></param>
        public static void ChooseSolver(AppControl Control, XdgTimesteppingBase Timestepper, bool Convection, double muA) {

            if (Control == null) {
                throw new Exception("Control object empty");
            }
            //if (Timestepper == null) {
            //    throw new Exception("Timestepper object empty");
            //}

            LinearSolverConfig linConfig = Control.LinearSolver;
            NonLinearSolverConfig NonlinConfig = Control.NonLinearSolver;

            // Set several solver options for Timestepper
            Timestepper.Config_SolverConvergenceCriterion = linConfig.Solver_ConvergenceCriterion;
            Timestepper.Config_MaxIterations = linConfig.MaxSolverIterations;
            Timestepper.Config_MinIterations = linConfig.MinSolverIterations;
            Timestepper.Config_MaxKrylovDim = linConfig.MaxKrylovDim;

            // Set to pseudo Picard if the Stokes equations should be solved
            if (Convection == false)
                NonlinConfig.SolverCode = NonLinearSolverConfig.Code.Picard;
            
            ISolverSmootherTemplate templinearSolve = null;

            switch (linConfig.SolverCode) {
                case LinearSolverConfig.Code.automatic:
                    templinearSolve = AutomaticChoice(Control, Timestepper);
                    break;

                case LinearSolverConfig.Code.classic_mumps:
                    templinearSolve = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };
                    break;

                case LinearSolverConfig.Code.classic_pardiso:
                    templinearSolve = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_directcoarse_overlap:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_directcoarse:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_Kcycle_directcoarse:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = Control.NoOfMultigridLevels - 1
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_Kcycle_directcoarse_overlap:

                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = Control.NoOfMultigridLevels - 1
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(Control.NoOfMultigridLevels - 2)
                    };
                    break;

                case LinearSolverConfig.Code.exp_softgmres:
                    templinearSolve = new SoftGMRES() {
                        MaxKrylovDim = Timestepper.Config_MaxKrylovDim,
                        m_Tolerance = Timestepper.Config_SolverConvergenceCriterion,
                    };
                    break;

                case LinearSolverConfig.Code.exp_softgmres_schwarz_Kcycle_directcoarse_overlap:
                    templinearSolve = new SoftGMRES() {
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

                case LinearSolverConfig.Code.exp_softgmres_schwarz_directcoarse_overlap:
                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");
                    templinearSolve = new SoftGMRES() {
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

                case LinearSolverConfig.Code.exp_multigrid:
                    if (Control.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");
                    templinearSolve = new ILU() { };
                    break;

                case LinearSolverConfig.Code.exp_ILU:
                    templinearSolve = new ILU() { };
                    break;

                case LinearSolverConfig.Code.exp_Schur:
                    templinearSolve = new SchurPrecond()
                    {
                        SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
                    };
                    break;

                case LinearSolverConfig.Code.exp_Simple:
                    templinearSolve = new SchurPrecond()
                    {
                        SchurOpt = SchurPrecond.SchurOptions.SIMPLE
                    };
                    break;

                case LinearSolverConfig.Code.exp_AS_1000:
                    if (Timestepper.MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF 
                    {
                        templinearSolve = new Schwarz()
                        {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
                            {
                                //noofparts = 76,
                                NoOfPartsPerProcess = 213, // Warum 76
                            },
                            CoarseSolver = new DirectSolver()
                            {
                                WhichSolver = DirectSolver._whichSolver.MUMPS    //PARDISO
                            },
                            Overlap = 1
                        };
                    }
                    else  //2D --> 75088DoF
                    {
                        templinearSolve = new Schwarz()
                        {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
                            {
                                //noofparts = 213,
                                NoOfPartsPerProcess = 213,
                            },
                            CoarseSolver = new DirectSolver()
                            {
                                WhichSolver = DirectSolver._whichSolver.MUMPS    //PARDISO
                            },
                            Overlap = 1
                        };
                    }
                    break;

                case LinearSolverConfig.Code.exp_AS_5000:
                    if (Timestepper.MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF
                    {
                        templinearSolve = new Schwarz()
                        {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
                            {
                                //noofparts = 43,
                                NoOfPartsPerProcess = 43,
                            },
                            CoarseSolver = new DirectSolver()
                            {
                                WhichSolver = DirectSolver._whichSolver.MUMPS    //PARDISO
                            },
                            Overlap = 1
                        };
                    }
                    else  //2D --> 75088DoF
                    {
                        templinearSolve = new Schwarz()
                        {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
                            {
                                //noofparts = 16,
                                NoOfPartsPerProcess = 43,
                            },
                            CoarseSolver = new DirectSolver()
                            {
                                WhichSolver = DirectSolver._whichSolver.MUMPS    //PARDISO
                            },
                            Overlap = 1
                        };
                    } 
                        
                    break;

                case LinearSolverConfig.Code.exp_AS_10000:
                    if (Timestepper.MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF
                    {
                        templinearSolve = new Schwarz()
                        {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
                            {
                                //noofparts = 22,
                                NoOfPartsPerProcess = 22,

                            },
                            CoarseSolver = new DirectSolver()
                            {
                                WhichSolver = DirectSolver._whichSolver.MUMPS    //PARDISO
                            },
                            Overlap = 1
                        };
                    }
                    else  //2D --> 75088DoF
                    {
                        templinearSolve = new Schwarz()
                        {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
                            {
                                //noofparts = 8,
                                NoOfPartsPerProcess = 22, //

                            },
                            CoarseSolver = new DirectSolver()
                            {
                                WhichSolver = DirectSolver._whichSolver.MUMPS    //PARDISO
                            },
                            Overlap = 1
                        };
                    }
                        
                    break;

                case LinearSolverConfig.Code.exp_AS_MG:
                    templinearSolve = new Schwarz()
                    {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks()
                        {
                            //depth = asdepth,
                            Depth = 2,
                        },
                        CoarseSolver = new DirectSolver()
                        {
                            WhichSolver = DirectSolver._whichSolver.MUMPS    //PARDISO
                        },

                        Overlap = 1
                    };
                    break;


                case LinearSolverConfig.Code.exp_localPrec:
                    templinearSolve = new LocalizedOperatorPrec()
                    {
                        m_dt = Control.GetFixedTimestep(),
                        m_muA = muA,
                    };
                    break;

                case LinearSolverConfig.Code.classic_cg:
                    templinearSolve = new ilPSP.LinSolvers.monkey.CG() {
                        MaxIterations = 1000000,
                        Tolerance = 1.0e-10,
                        DevType = ilPSP.LinSolvers.monkey.DeviceType.Cuda
                    };
                    break;

                
                default:
                    throw new NotImplementedException("Linear solver option not available");
            }

            

            // Set nonlinear Solver
            switch (NonlinConfig.SolverCode)
            {
                case NonLinearSolverConfig.Code.NewtonGMRES:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.NewtonGMRES;
                    Timestepper.Config_linearSolver = templinearSolve;
                    break;
                case NonLinearSolverConfig.Code.Picard:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.Picard;
                    Timestepper.Config_linearSolver = templinearSolve;
                    break;
                case NonLinearSolverConfig.Code.Newton:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.Newton;
                    Timestepper.Config_linearSolver = templinearSolve;
                    break;
                case NonLinearSolverConfig.Code.PicardGMRES:
                    Timestepper.Config_NonlinearSolver = NonlinearSolverMethod.Picard;
                    Timestepper.Config_linearSolver = new SoftGMRES()
                    {
                        MaxKrylovDim = Timestepper.Config_MaxKrylovDim,
                        m_Tolerance = Timestepper.Config_SolverConvergenceCriterion,
                        Precond = templinearSolve,
                        m_MaxIterations = Timestepper.Config_MaxIterations,
                    };
                    break;
                default:
                    throw new NotImplementedException("Nonlinear solver option not available");
            }

        }

        /// <summary>
        /// Automatic choice of linear solver depending on problem size, immersed boundary, polynomial degree, etc.
        /// </summary>
        static ISolverSmootherTemplate AutomaticChoice(AppControl Control, XdgTimesteppingBase Timestepper) {

            //int pV = Control.FieldOptions["VelocityX"].Degree;
            int pP = Control.FieldOptions["Pressure"].Degree;
            int pV = pP + 1;

            // Detecting variables for solver determination 
            var D = Timestepper.MultigridSequence[0].SpatialDimension;
            var cellsLoc = Timestepper.MultigridSequence[0].CellPartitioning.LocalLength;
            var cellsGlo = Timestepper.MultigridSequence[0].CellPartitioning.TotalLength;

            ISolverSmootherTemplate tempsolve = null;


            var size = Timestepper.MultigridSequence[0].CellPartitioning.MpiSize;

            // !!!!!!!!!!!UNTERSCHEIDUNG OB PICARD ODER NEWTON!!!!!!!!!!!!
            if (Timestepper.Config_NonlinearSolver == NonlinearSolverMethod.NewtonGMRES) {

                // Spatial Dimension
                switch (D) {
                    case 1:
                        break;
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                        //break;

                    case 2:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                        //break;

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
                        //break;

                    case 2:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                        //break;

                    case 3:
                        var dofsPerCell3D = (3 * (pV * pV * pV + 6 * pV * pV + 11 * pV + 6) / 6 + 1 * (pP * pP * pP + 6 * pP * pP + 11 * pP + 6) / 6);
                        var dofsLoc = dofsPerCell3D * cellsLoc;
                        var dofsGlo = dofsPerCell3D * cellsGlo;

                        if (dofsGlo > 10000) {

                            if (Control.NoOfMultigridLevels < 2)
                                throw new ApplicationException("At least 2 Multigridlevels are required");

                            tempsolve = new SoftGMRES() {
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
                            tempsolve = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };
                        }
                        break;

                    default:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                }

                
            }

            return tempsolve;
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

        private class Converted (ilPSP.LinSolvers.ISparseSolver bla) : ISolverSmootherTemplate {
            ISparseSolver thisCG = new ilPSP.LinSolvers.monkey.CG();
            ISolverSmootherTemplate ISST;
            ISST.Converged
            ISST.Init(){
            }

           
        }

    }
}
