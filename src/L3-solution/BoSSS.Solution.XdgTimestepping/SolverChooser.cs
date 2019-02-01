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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using MPI.Wrappers;

namespace BoSSS.Solution.AdvancedSolvers {

    public class SolverFactory {

        /// <summary>
        /// This <see cref="SolverFactory"/> enables creation of linear and nonlinear solver objects. The configuration <see cref="LinearSolverConfig"/> and <see cref="NonLinearSolverConfig"/> can be set in Control-Files (defined in <see cref="AppControl"/>).
        /// </summary>
        /// <param name="nc"></param>
        /// <param name="lc"></param>
        public SolverFactory(NonLinearSolverConfig nc, LinearSolverConfig lc) {
            m_lc = lc;
            m_nc = nc;
        }

        /// <summary>
        /// This will return <code>linear</code> and <code>nonlinear</code> solver objects, which are configured according to <see cref="LinearSolverConfig"/> and <see cref="NonLinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>).
        /// </summary>
        /// <param name="nonlinSolver"></param>
        /// <param name="linsolver"></param>
        /// <param name="Timestepper"></param>
        /// <param name="ts_AssembleMatrixCallback"></param>
        /// <param name="ts_MultigridBasis"></param>
        /// <param name="LevelSetConvergenceReached"></param>
        /// <param name="PseudoNonlinear"></param>
        public void GenerateNonLin(out NonlinearSolver nonlinSolver, out ISolverSmootherTemplate linsolver, XdgTimesteppingBase Timestepper, OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, FixpointIterator.CoupledConvergenceReached LevelSetConvergenceReached, bool PseudoNonlinear) {

            if (m_nonlinsolver != null && (m_linsolver == null || m_precond == null))
                throw new NotImplementedException("an uncomplete nonlinear solver is overgiven.");

            linsolver = null;
            nonlinSolver = null;
            ISolverSmootherTemplate precondsolver = null;

            precondsolver = LinSolverDealer(Timestepper,true);
            linsolver = LinSolverDealer(Timestepper,false);
            Debug.Assert(linsolver != null);

            nonlinSolver = GenerateNonLin_body(Timestepper, ts_AssembleMatrixCallback, ts_MultigridBasis, LevelSetConvergenceReached, PseudoNonlinear, m_nc, m_lc, linsolver,precondsolver);

            Debug.Assert(nonlinSolver !=null);
            return;
        }

        /// <summary>
        /// Internal method to determine, it there is already a selfmade solver set. If there is none and <see cref="LinearSolverConfig.Code.automatic"/> is selected, then the nonlinear solver slection is considered in the choice of the linear solver, with <see cref="AutomaticChoice(Timestepper, m_nc, _lc)"/>.
        /// </summary>
        /// <param name="Timestepper"></param>
        /// <param name="IsPrecond"></param>
        /// <returns></returns>
        private ISolverSmootherTemplate LinSolverDealer(XdgTimesteppingBase Timestepper, bool IsPrecond) {
            LinearSolverConfig _lc;
            if (IsPrecond) {
                _lc = m_nc.PrecondSolver;
                if (m_precond != null) {
                    Check_precond();
#if DEBUG
                    Console.WriteLine("Preconditioner was handed over");
#endif
                    return m_precond;
                }
            } else {
                _lc = m_lc;
                if (m_linsolver != null) {
                    Check_linsolver();
#if DEBUG
                    Console.WriteLine("LinearSolver was handed over");
#endif
                    return m_linsolver;
                }
            }
            if (_lc != null) {
                return GenerateLinear_body(Timestepper, _lc , m_nc);
            } else {
                throw new ArgumentNullException("LinearSolver Config is null!");
            }
        }

        /// <summary>
        /// This one is the method-body of <see cref="GenerateNonLinear"/> and shall not be called from the outside. The parameters are mainly handed over to the NonLinearSolver object, which lives in <see cref="AdvancedSolvers.NonlinearSolver"/>.
        /// </summary>
        /// <param name="Timestepper"></param>
        /// <param name="ts_AssembleMatrixCallback"></param>
        /// <param name="ts_MultigridBasis"></param>
        /// <param name="LevelSetConvergenceReached"></param>
        /// <param name="PseudoNonlinear"></param>
        /// <param name="nc"></param>
        /// <param name="lc"></param>
        /// <param name="LinSolver"></param>
        /// <param name="PrecondSolver"></param>
        /// <returns></returns>
        private NonlinearSolver GenerateNonLin_body(XdgTimesteppingBase Timestepper, OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, FixpointIterator.CoupledConvergenceReached LevelSetConvergenceReached, bool PseudoNonlinear, NonLinearSolverConfig nc, LinearSolverConfig lc, ISolverSmootherTemplate LinSolver, ISolverSmootherTemplate PrecondSolver) {

            // +++++++++++++++++++++++++++++++++++++++++++++
            // the nonlinear solvers:
            // +++++++++++++++++++++++++++++++++++++++++++++

            NonlinearSolver nonlinSolver;

            // Set to pseudo Picard if the Stokes equations should be solved
            if (PseudoNonlinear == true)
                nc.SolverCode = NonLinearSolverConfig.Code.Picard;

            switch (nc.SolverCode) {
                    case NonLinearSolverConfig.Code.Picard:

                        nonlinSolver = new FixpointIterator(
                            ts_AssembleMatrixCallback,
                            ts_MultigridBasis,
                            Timestepper.Config_MultigridOperator) {
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            m_LinearSolver = LinSolver,
                            m_SessionPath = Timestepper.SessionPath, //is needed for Debug purposes, output of inter-timesteps
                            ConvCrit = nc.ConvergenceCriterion,
                            UnderRelax = nc.UnderRelax,
                        };

                        if (Timestepper.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                            ((FixpointIterator)nonlinSolver).CoupledIteration_Converged = LevelSetConvergenceReached;
                        }

                        break;

                //Besides NonLinearSolverConfig Newton needs also LinearSolverConfig
                //Newton uses MUMPS as linearsolver by default
                case NonLinearSolverConfig.Code.Newton:
                    
                        nonlinSolver = new Newton(
                            ts_AssembleMatrixCallback,
                            ts_MultigridBasis,
                            Timestepper.Config_MultigridOperator) {
                            maxKrylovDim = lc.MaxKrylovDim, 
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            ApproxJac = Newton.ApproxInvJacobianOptions.DirectSolver, //MUMPS is taken, todo: enable all linear solvers
                            Precond = PrecondSolver,
                            GMRESConvCrit = lc.ConvergenceCriterion,
                            ConvCrit = nc.ConvergenceCriterion,
                            m_SessionPath = Timestepper.SessionPath,
                        };
                        break;

                //in NewtonGMRES Newton is merged with GMRES, this is an optimized algorithm
                //NonLinearSolver and LinearSolver can not be separated in this case
                case NonLinearSolverConfig.Code.NewtonGMRES:

                        nonlinSolver = new Newton(
                            ts_AssembleMatrixCallback,
                            ts_MultigridBasis,
                            Timestepper.Config_MultigridOperator) {
                            maxKrylovDim = lc.MaxKrylovDim,
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            ApproxJac = Newton.ApproxInvJacobianOptions.GMRES,
                            Precond = PrecondSolver,
                            //Precond_solver = new RheologyJacobiPrecond() { m_We = 0.1},
                            GMRESConvCrit = lc.ConvergenceCriterion,
                            ConvCrit = nc.ConvergenceCriterion,
                            m_SessionPath = Timestepper.SessionPath,
                        };
                        break;

                    case NonLinearSolverConfig.Code.PicardGMRES:

                    nonlinSolver = new FixpointIterator(
                            ts_AssembleMatrixCallback,
                            ts_MultigridBasis,
                            Timestepper.Config_MultigridOperator) {
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            m_LinearSolver = LinSolver,
                            m_SessionPath = Timestepper.SessionPath, //is needed for Debug purposes, output of inter-timesteps
                            ConvCrit = nc.ConvergenceCriterion,
                            UnderRelax = nc.UnderRelax,
                            Precond = PrecondSolver,
                        };
                        break;

                default:
                        throw new NotImplementedException();
                }
#if DEBUG
            Console.WriteLine("nonlinear Solver: {0}",nc.SolverCode.ToString());
#endif
            return nonlinSolver;
        }

        /// <summary>
        /// This will return a <code>linear</code> solver object, which is configured according to <see cref="LinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>). 
        /// </summary>
        /// <param name="templinearSolve"></param>
        /// <param name="Timestepper"></param>
        public void GenerateLinear(out ISolverSmootherTemplate templinearSolve, XdgTimesteppingBase Timestepper) {
            if (m_linsolver != null) {
                templinearSolve = m_linsolver;
            } else {
                templinearSolve = GenerateLinear_body(Timestepper, m_lc, null);
            }
            Debug.Assert(templinearSolve!=null);
            return;
        }

        /// <summary>
        /// This one is the method-body of <see cref="GenerateLinear"/> and shall not be called from the outside. Some Solver aquire additional information, thus the timestepper is overgiven as well.
        /// </summary>
        /// <param name="Timestepper"></param>
        /// <param name="lc"></param>
        /// <returns></returns>
        private ISolverSmootherTemplate GenerateLinear_body(XdgTimesteppingBase Timestepper, LinearSolverConfig lc, NonLinearSolverConfig nc) {

            // +++++++++++++++++++++++++++++++++++++++++++++
            // the linear solvers:
            // +++++++++++++++++++++++++++++++++++++++++++++

            if (lc == null)
                throw new ArgumentNullException();

            ISolverSmootherTemplate templinearSolve=null;

            switch (lc.SolverCode) {
                case LinearSolverConfig.Code.automatic:
                    if (nc != null) {
                        templinearSolve = Automatic(Timestepper, nc, lc);
                    } else {
                        templinearSolve = AutomaticLinearOnly(Timestepper, lc);
                    }
                    break;

                case LinearSolverConfig.Code.classic_mumps:
                    templinearSolve = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.MUMPS,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverConfig.Code.classic_pardiso:
                    templinearSolve = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.PARDISO,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_directcoarse_overlap:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_directcoarse:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_Kcycle_directcoarse:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = lc.NoOfMultigridLevels - 1
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverConfig.Code.exp_schwarz_Kcycle_directcoarse_overlap:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = lc.NoOfMultigridLevels - 1
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverConfig.Code.exp_softgmres:

                    templinearSolve = new SoftGMRES() {
                        MaxKrylovDim = lc.MaxKrylovDim,
                        m_Tolerance = lc.ConvergenceCriterion,
                    };
                    break;

                case LinearSolverConfig.Code.exp_softgmres_schwarz_Kcycle_directcoarse_overlap:
                    templinearSolve = new SoftGMRES() {
                        MaxKrylovDim = lc.MaxKrylovDim,
                        m_Tolerance = lc.ConvergenceCriterion,
                        Precond = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                Depth = lc.NoOfMultigridLevels - 1
                            },
                            Overlap = 1,
                            CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                        },
                    };
                    break;

                case LinearSolverConfig.Code.exp_softgmres_schwarz_directcoarse_overlap:
                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");
                    templinearSolve = new SoftGMRES() {
                        MaxKrylovDim = lc.MaxKrylovDim,
                        m_Tolerance = lc.ConvergenceCriterion,
                        Precond = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsPerProcess = 1,
                            },
                            Overlap = 1,
                            CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2,lc)
                        },
                    };
                    break;

                case LinearSolverConfig.Code.exp_multigrid:
                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");
                    templinearSolve = new ILU() { };
                    break;

                case LinearSolverConfig.Code.exp_ILU:
                    templinearSolve = new ILU() { };
                    break;

                case LinearSolverConfig.Code.exp_Schur:
                    templinearSolve = new SchurPrecond() {
                        SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
                    };
                    break;

                case LinearSolverConfig.Code.exp_Simple:
                    templinearSolve = new SchurPrecond() {
                        SchurOpt = SchurPrecond.SchurOptions.SIMPLE
                    };
                    break;

                case LinearSolverConfig.Code.exp_AS_1000:
                    if (Timestepper.MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF 
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 76,
                                NoOfPartsPerProcess = 213, // Warum 76
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    } else  //2D --> 75088DoF
                      {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 213,
                                NoOfPartsPerProcess = 213,
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    }
                    break;

                case LinearSolverConfig.Code.exp_AS_5000:
                    if (Timestepper.MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 43,
                                NoOfPartsPerProcess = 43,
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    } else  //2D --> 75088DoF
                      {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 16,
                                NoOfPartsPerProcess = 43,
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    }

                    break;

                case LinearSolverConfig.Code.exp_AS_10000:
                    if (Timestepper.MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 22,
                                NoOfPartsPerProcess = 22,

                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    } else  //2D --> 75088DoF
                      {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 8,
                                NoOfPartsPerProcess = 22, //

                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    }

                    break;

                case LinearSolverConfig.Code.exp_AS_MG:
                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            //depth = asdepth,
                            Depth = 2,
                        },
                        CoarseSolver = new SparseSolver() {
                            WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                            LinConfig = lc
                        },

                        Overlap = 1
                    };
                    break;


                case LinearSolverConfig.Code.exp_localPrec:
                    templinearSolve = new LocalizedOperatorPrec() {
                        m_dt = lc.exp_localPrec_Min_dt,
                        m_muA = lc.exp_localPrec_muA,
                    };
                    break;

                case LinearSolverConfig.Code.classic_cg:
                    templinearSolve = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.CG,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverConfig.Code.exp_softpcg_mg:
                    //templinearSolve = new SpecialMultilevelSchwarz();
                    //Timestepper.Config_MultigridOperator
                    

                    break;


                default:
                    throw new NotImplementedException("Linear solver option not available");
            }
            Debug.Assert(templinearSolve != null);
#if DEBUG
            Console.WriteLine("linear Solver: {0}",lc.SolverCode.ToString());
#endif
            return templinearSolve;
        }

        /// <summary>
        /// Is get and set by <see cref="Selfmade_linsolver"/> and used by <see cref="GenerateLinear"/>.
        /// </summary>
        private ISolverSmootherTemplate m_linsolver;
        
        /// <summary>
        /// Is get and set by <see cref="Selfmade_nonlinsolver"/> and used by <see cref="GenerateNonLinear"/>.
        /// </summary>
        private NonlinearSolver m_nonlinsolver;

        /// <summary>
        /// Is get and set by <see cref="Selfmade_precond"/> and used by <see cref="GenerateLinear"/>.
        /// </summary>
        private ISolverSmootherTemplate m_precond;

        /// <summary>
        /// Internal linear solver configuration. Shall always be != null. 
        /// </summary>
        private LinearSolverConfig m_lc;

        /// <summary>
        /// Internal nonlinear solver configuration. Shall always be != null. 
        /// </summary>
        private NonLinearSolverConfig m_nc;

        /// <summary>
        /// For developers, who want full control over solvers: In <see cref="selfmade_linsolver"/> you can insert your own config of linear solver.
        /// Set the solver to <c>null</c> to enable solver generation from <see cref="LinearSolverConfig"/> again.
        /// </summary>
        public ISolverSmootherTemplate Selfmade_linsolver {
            set {
                m_linsolver = value;
            }
            get {
                return m_linsolver;
            }
        }

        /// <summary>
        /// For developers, who want full control over solvers: In <see cref="selfmade_nonlinsolver"/> you can insert your own config of nonlinear solver,
        /// which will overwrite the output of <see cref="GenerateNonLinear"/> with the overgiven solver.
        /// Set the solver to <c>null</c> to enable solver generation from <see cref="NonLinearSolverConfig"/> again.
        /// Note: The overgiven solver has to be completely defined (precond!=null and linsolve!=null) 
        /// </summary>
        public NonlinearSolver Selfmade_nonlinsolver {
            set {
                m_nonlinsolver = value;
            }
            get {
                return m_nonlinsolver;
            }
        }

        /// <summary>
        /// For developers, who want full control over solvers: In <see cref="Selfmade_precond"/> you can insert your own config of linear solver.
        /// Set the solver to <c>null</c> to enable solver generation from <see cref="NonLinearSolverConfig.Precond_solver"/> again.
        /// </summary>
        public ISolverSmootherTemplate Selfmade_precond {
            set {
                m_precond = value;
            }
            get {
                return m_precond;
            }
        }

        /// <summary>
        /// Automatic Chooser for LinearSolver standalone
        /// </summary>
        /// <param name="Timestepper"></param>
        /// <param name="lc"></param>
        private ISolverSmootherTemplate AutomaticLinearOnly(XdgTimesteppingBase Timestepper, LinearSolverConfig lc) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Automatic choice of linear solver depending on problem size, immersed boundary, polynomial degree, etc. In addition the nonlinearsolver config is considered as well.
        /// </summary>
        private ISolverSmootherTemplate Automatic(XdgTimesteppingBase Timestepper, NonLinearSolverConfig nc, LinearSolverConfig lc) {

            var D = Timestepper.MultigridSequence[0].SpatialDimension;

            //int pV = Control.FieldOptions["VelocityX"].Degree;
            int pP = Timestepper.Config_MultigridOperator[0][D].Degree; //Control.FieldOptions["Pressure"].Degree;
            int pV = pP + 1;

            // Detecting variables for solver determination 

            var cellsLoc = Timestepper.MultigridSequence[0].CellPartitioning.LocalLength;
            var cellsGlo = Timestepper.MultigridSequence[0].CellPartitioning.TotalLength;

            ISolverSmootherTemplate tempsolve = null;

            var size = Timestepper.MultigridSequence[0].CellPartitioning.MpiSize;

            // !!!!!!!!!!!UNTERSCHEIDUNG OB PICARD ODER NEWTON!!!!!!!!!!!!
            if (nc.SolverCode == NonLinearSolverConfig.Code.NewtonGMRES) {

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

                            if (lc.NoOfMultigridLevels < 2)
                                throw new ApplicationException("At least 2 Multigridlevels are required");

                            tempsolve = new Schwarz() {
                                m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                    NoOfPartsPerProcess = PPP,
                                },
                                Overlap = 1,
                                CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2,lc)
                            };
                        } else {
                            tempsolve = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,
                                LinConfig = lc
                            };
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

                            if (lc.NoOfMultigridLevels < 2)
                                throw new ApplicationException("At least 2 Multigridlevels are required");

                            tempsolve = new SoftGMRES() {
                                MaxKrylovDim = lc.MaxKrylovDim,
                                m_Tolerance = lc.ConvergenceCriterion,
                                Precond = new Schwarz() {
                                    m_BlockingStrategy = new Schwarz.SimpleBlocking() {
                                        NoOfPartsPerProcess = (int)Math.Ceiling(dofsLoc / 6500.0),
                                    },
                                    Overlap = 1,
                                    CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2,lc)
                                },
                            };
                        } else {
                            tempsolve = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,
                                LinConfig = lc
                            };
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
        private ISolverSmootherTemplate DetermineMGSquence(int MGlevels,LinearSolverConfig lc) {
            ISolverSmootherTemplate solver;
            if (MGlevels > 0) {
                solver = new ClassicMultigrid() { CoarserLevelSolver = DetermineMGSquence(MGlevels - 1, lc) };
            } else {
                solver = new SparseSolver() {
                    WhichSolver = SparseSolver._whichSolver.MUMPS,
                    LinConfig = lc
                };
            }
            return solver;
        }

        /// <summary>
        /// Updates solver configuration.
        /// </summary>
        /// <param name="nc"></param>
        /// <param name="lc"></param>
        public void Update(NonLinearSolverConfig nc, LinearSolverConfig lc) {
            this.m_lc = lc;
            this.m_nc = nc;
        }

        /// <summary>
        /// clears overgiven selfmade solvers
        /// </summary>
        public void Clear() {
            //this.m_lc = null;
            //this.m_nc = null;
            this.m_linsolver = null;
            this.m_nonlinsolver = null;
        }

        /// <summary>
        /// Checks overgiven selfmade linear solver
        /// </summary>
        /// <returns></returns>
        private bool Check_linsolver() {
            bool check = true;
            //test something ... m_linsolver
            return check;
        }

        /// <summary>
        /// Checks overgiven selfmade proconditioner solver
        /// </summary>
        /// <returns></returns>
        private bool Check_precond() {
            bool check = true;
            //test something ... m_precond
            return check;
        }

        /// <summary>
        /// Checks overgiven selfmade nonlinear solver
        /// </summary>
        /// <returns></returns>
        private bool Check_nonlinsolver() {
            bool check = true;
            //test something ... m_nonlinsolver
            return check;
        }



        private ISolverSmootherTemplate SpecialMultilevelSchwarz(MultigridOperator op,XdgTimesteppingBase TS, LinearSolverConfig _lc) {
            var solver = new SoftPCG() {
                m_MaxIterations = 500,
                m_Tolerance = 1.0e-12
            };
            //var solver = new OrthonormalizationScheme() {
            //    MaxIter = 500,
            //    Tolerance = 1.0e-10,
            //};
            //var solver = new SoftGMRES() {
            //    m_MaxIterations = 500,
            //    m_Tolerance = 1.0e-10,

            //};

            // my tests show that the ideal block size may be around 10'000
            int DirectKickIn = _lc.TargetBlockSize;

            //TS.MultigridSequence[0].CellPartitioning.TotalLength*TS.Config_MultigridOperator.
            MultigridOperator Current = op;
            ISolverSmootherTemplate[] MultigridChain = new ISolverSmootherTemplate[_lc.NoOfMultigridLevels];
            for (int iLevel = 0; iLevel < _lc.NoOfMultigridLevels; iLevel++) {
                
                int SysSize = Current.Mapping.TotalLength;
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == _lc.NoOfMultigridLevels - 1;
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (useDirect) {
                    MultigridChain[iLevel] = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.PARDISO,
                        TestSolution = false
                    };
                } else {

                    ClassicMultigrid MgLevel = new ClassicMultigrid() {
                        m_MaxIterations = 1,
                        m_Tolerance = 0.0 // termination controlled by top level PCG
                    };


                    MultigridChain[iLevel] = MgLevel;



                    ISolverSmootherTemplate pre, pst;
                    if (iLevel > 0) {

                        Schwarz swz1 = new Schwarz() {
                            m_MaxIterations = 1,
                            CoarseSolver = null,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsPerProcess = NoOfBlocks
                            },
                            Overlap = 0 // overlap does **NOT** seem to help
                        };

                        SoftPCG pcg1 = new SoftPCG() {
                            m_MinIterations = 5,
                            m_MaxIterations = 5
                        };

                        SoftPCG pcg2 = new SoftPCG() {
                            m_MinIterations = 5,
                            m_MaxIterations = 5
                        };

                        var preChain = new ISolverSmootherTemplate[] { swz1, pcg1 };
                        var pstChain = new ISolverSmootherTemplate[] { swz1, pcg2 };

                        pre = new SolverSquence() { SolverChain = preChain };
                        pst = new SolverSquence() { SolverChain = pstChain };
                    } else {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++
                        // top level - use only iterative (non-direct) solvers
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++

                        pre = new BlockJacobi() {
                            NoOfIterations = 3,
                            omega = 0.5
                        };

                        pst = new BlockJacobi() {
                            NoOfIterations = 3,
                            omega = 0.5
                        };

                        //preChain = new ISolverSmootherTemplate[] { pcg1 };
                        //pstChain = new ISolverSmootherTemplate[] { pcg2 };
                    }





                    //if (iLevel > 0) {
                    //    MgLevel.PreSmoother = pre;
                    //    MgLevel.PostSmoother = pst;
                    //} else {
                    //    //MgLevel.PreSmoother = pcg1;   // ganz schlechte Idee, konvergiert gegen FALSCHE lösung
                    //    //MgLevel.PostSmoother = pcg2;  // ganz schlechte Idee, konvergiert gegen FALSCHE lösung
                    //    MgLevel.PreSmoother = pre;
                    //    MgLevel.PostSmoother = pst;
                    //}
                    

                    MgLevel.PreSmoother = pre;
                    MgLevel.PostSmoother = pst;
                }

                if (iLevel > 0) {
                    ((ClassicMultigrid)(MultigridChain[iLevel - 1])).CoarserLevelSolver = MultigridChain[iLevel];
                }

                if (useDirect) {
                    Console.WriteLine("MG: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, SysSize, DirectKickIn);
                    break;
                }



                Current = Current.CoarserLevel;

            } // end of level loop


            solver.Precond = MultigridChain[0];
            //solver.PrecondS = new[] { MultigridChain[0] };

            return solver;
        }
    }
}