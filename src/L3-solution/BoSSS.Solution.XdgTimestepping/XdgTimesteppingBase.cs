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

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Multigrid;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ilPSP;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Aggregation;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// Callback-template for level-set updates.
    /// </summary>
    /// <param name="OpMtx">
    /// Output for the linear part.
    /// </param>
    /// <param name="OpAffine">
    /// Output for the affine part.
    /// </param>
    /// <param name="Mapping">
    /// Corresponds with row and columns of <paramref name="OpMtx"/>, resp. with <paramref name="OpAffine"/>.
    /// </param>
    /// <param name="CurrentState"></param>
    /// <param name="CurrentAgg"></param>
    /// <param name="time"></param>
    public delegate void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, MultiphaseCellAgglomerator CurrentAgg, double time);

    /// <summary>
    /// Callback-template for level-set updates.
    /// </summary>
    /// <param name="CurrentState"></param>
    /// <param name="time">
    /// Actual simulation time for the known value;
    /// </param>
    /// <param name="dt">
    /// Timestep size.
    /// </param>
    /// <param name="UnderRelax">
    /// </param>
    /// <param name="incremental">
    /// true for Splitting schmes with subdivided level-set evolution (e.g. Strang-Splitting)
    /// </param>
    /// <returns>
    /// Some kind of level-set-residual in order to check convergence in a fully coupled simulation
    /// (see <see cref="LevelSetHandling.Coupled_Iterative"/>)
    /// </returns>
    public delegate double DelUpdateLevelset(DGField[] CurrentState, double time, double dt, double UnderRelax, bool incremental);

    /// <summary>
    /// Callback-template for to obtain cut-cell metrics.
    /// </summary>
    public delegate CutCellMetrics DelUpdateCutCellMetrics();

    /// <summary>
    /// Callback-template for pushing the level-set in case of increment timestepping
    /// </summary>
    public delegate void DelPushLevelSetRelatedStuff();

    /// <summary>
    /// Controls the updating of the mass-matrix, resp. the temporal operator.
    /// </summary>
    public enum MassMatrixShapeandDependence {
        /// <summary>
        /// May be useful e.g. if one runs the FSI solver as a pure single-phase solver,
        /// i.e. if the Level-Set is outside the domain.
        /// </summary>
        IsIdentity = 0,

        /// <summary>
        /// For a mass-matrix which is not the identity, but constant over time.
        /// E.g., if the level-set is constant in time, the mass matrix does not change.
        /// </summary>
        IsNonIdentity = 1,

        /// <summary>
        /// For a mass matrix which changes over time, but does not depend on the solution.
        /// E.g., the interface motion is prescribed.
        /// </summary>
        IsTimeDependent = 2,

        /// <summary>
        /// The mass matrix depends on time and the solution, e.g. for multiphase flows,
        /// where the position of the fluid interface is influenced by the flow field.
        /// </summary>
        IsTimeAndSolutionDependent = 3
    }

    /// <summary>
    /// Which nonlinear solver should be used
    /// </summary>
    public enum NonlinearSolverMethod {
        /// <summary>
        /// Standard Fixpoint-Iteration. Convergence only linear.
        /// </summary>
        Picard = 0,

        /// <summary>
        /// Newtons method coupled with GMRES as default. Convergeces quadratically but needs a good approximation. Preconditioning of the Newton-GMRES is set to LinearSolver.
        /// </summary>
        Newton = 1,
    }

    /// <summary>
    /// Treatment of the level-set.
    /// </summary>
    public enum LevelSetHandling {
        /// <summary>
        /// Level-Set remains constant.
        /// </summary>
        None = 0,

        /// <summary>
        /// Level-Set is handled using Lie-Splitting.
        /// </summary>
        LieSplitting = 1,

        /// <summary>
        /// Level-Set is handled using Strang-Splitting.
        /// </summary>
        StrangSplitting = 2,

        /// <summary>
        /// Level-Set is updated once per time-step.
        /// </summary>
        Coupled_Once = 3,

        /// <summary>
        /// Level-Set is updated in very iteration, until convergence is reached.
        /// </summary>
        Coupled_Iterative = 4
    }

    public enum SpatialOperatorType {

        /// <summary>
        /// The spatial operator is linear, but time dependent (e.g. a Diffusion operator with time-dependent coefficient).
        /// </summary>
        LinearTimeDependent = 1,

        /// <summary>
        /// The spatial operator is non-linear, therefore it is also considered to be time-dependent;
        /// however, a linearization must be available.
        /// </summary>
        Nonlinear = 2
    }

    /// <summary>
    /// Common base-class for XDG timesteppers.
    /// </summary>
    abstract public class XdgTimesteppingBase {

        /// <summary>
        /// Configuration sanity checks, to be used by constructors of derived classes.
        /// </summary>
        protected void CommonConfigurationChecks() {
            if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                && this.Config_LevelSetHandling != LevelSetHandling.None)
                throw new ArgumentOutOfRangeException("illegal configuration");

            if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity && this.Config_MassScale != null) {
                // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                // i.e. if the Level-Set is outside the domain.

                foreach (var kv in this.Config_MassScale) {
                    SpeciesId spId = kv.Key;
                    double[] scaleVec = kv.Value.ToArray();
                    //for (int i = 0; i < scaleVec.Length; i++) {
                    //    if (scaleVec[i] != 1.0)
                    //        throw new ArithmeticException(string.Format("XDG time-stepping: illegal mass scale, mass matrix option {0} is set, but scaling factor for species {1}, variable no. {2} ({3}) is set to {4} (expecting 1.0).",
                    //            MassMatrixShapeandDependence.IsIdentity, this.m_LsTrk.GetSpeciesName(kv.Key), i, this.CurrentStateMapping.Fields[i].Identification, scaleVec[i]));
                    //}
                }

            }
        }


        public MassMatrixShapeandDependence Config_MassMatrixShapeandDependence {
            get;
            protected set;
        }


        public LevelSetHandling Config_LevelSetHandling {
            get;
            protected set;
        }

        public NonlinearSolverMethod Config_NonlinearSolver {
            get;
            set;
        }

        /// <summary>
        /// Convergence criterion if a nonlinear solver has to be used.
        /// </summary>
        public double Config_SolverConvergenceCriterion = 1.0e-8;

        public double Config_LevelSetConvergenceCriterion = 1.0e-6;

        /// <summary>
        /// Maximum number of iterations for an iterative solver (linear or nonlinear).
        /// </summary>
        public int Config_MaxIterations = 1000;

        /// <summary>
        /// Session path for writing in database
        /// </summary>
        public string SessionPath = "";

        /// <summary>
        /// Minimum number of iterations for an iterative solver (linear or nonlinear).
        /// </summary>
        public int Config_MinIterations = 4;

        /// <summary>
        /// Default Solver for the linear system
        /// </summary>
        public ISolverSmootherTemplate Config_linearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };

        /// <summary>
        /// Scaling of the mass matrix, for each species and each variable.
        /// </summary>
        public IDictionary<SpeciesId, IEnumerable<double>> Config_MassScale {
            get;
            protected set;
        }

        /// <summary>
        /// Whether the operator is linear, nonlinear.
        /// </summary>
        public SpatialOperatorType Config_SpatialOperatorType {
            get;
            protected set;
        }

        /// <summary>
        /// Callback routine to update the operator matrix.
        /// </summary>
        public DelComputeOperatorMatrix ComputeOperatorMatrix {
            get;
            protected set;
        }


        /// <summary>
        /// Callback routine to update the level set.
        /// </summary>
        public DelUpdateLevelset UpdateLevelset {
            get;
            protected set;
        }

        /// <summary>
        /// Callback routine to update the cut-cell metrics.
        /// </summary>
        public DelUpdateCutCellMetrics UpdateCutCellMetrics {
            get;
            protected set;
        }

        /// <summary>
        /// As usual the threshold for cell agglomeration.
        /// </summary>
        public double Config_AgglomerationThreshold {
            get;
            set;
        }

        /// <summary>
        /// True, if a nonlinear solver has to be used.
        /// </summary>
        protected virtual bool RequiresNonlinearSolver {
            get {
                // ask before changing


                if (Config_SpatialOperatorType == SpatialOperatorType.Nonlinear)
                    return true;

                if (Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative)
                    return true;

                //if (Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent)
                //    return true;

                return false;
            }

        }

        protected LevelSetTracker m_LsTrk;

        /// <summary>
        /// The sequence of aggregation grids, on which multi grid solvers work.
        /// </summary>
        protected AggregationGrid[] MultigridSequence;

        /// <summary>
        /// The aggregation grid basis.
        ///  - 1st index: grid level
        ///  - 2nd index: variable
        /// </summary>
        protected AggregationGridBasis[][] MultigridBasis;

        /*
        /// <summary>
        /// The XDG aggregation grid basis, for each multi grid level, which must be updated after each level set update.
        /// </summary>
        protected XdgAggregationBasis[] m_XdgAggregationBasis;
        */

        /// <summary>
        /// Last solver residuals.
        /// </summary>
        public CoordinateVector Residuals {
            get;
            protected set;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="ProblemMapping"></param>
        /// <param name="MultigridSequence"></param>
        /// <param name="useX">
        /// Create an XDG aggregation basis even if <paramref name="ProblemMapping"/> only contains DG basis.
        /// </param>
        /// <returns></returns>
        public static AggregationGridBasis[][] CreateAggregationGridBasis(UnsetteledCoordinateMapping ProblemMapping, LevelSetTracker _LsTrk, AggregationGrid[] MultigridSequence, bool useX) {

            Debug.Assert(object.ReferenceEquals(_LsTrk.GridDat, ProblemMapping.GridDat));
            Basis[] BasisS = ProblemMapping.BasisS.ToArray();
            Debug.Assert(BasisS.Where(b => !object.ReferenceEquals(b.GridDat, _LsTrk.GridDat)).Count() == 0);
            AggregationGridBasis[][] MultigridBasis = new AggregationGridBasis[MultigridSequence.Length][];



            XDGBasis maxXDGbasis = null;
            Basis maxDGbasis = null;
            foreach (var b in BasisS) {
                if (b is XDGBasis) {
                    XDGBasis xb = (XDGBasis)b;
                    if (maxXDGbasis == null || maxXDGbasis.Degree < xb.Degree)
                        maxXDGbasis = xb;
                }
                else {
                    if (maxDGbasis == null || maxDGbasis.Degree < b.Degree)
                        maxDGbasis = b;
                }
            }

            if (useX) {
                if (maxDGbasis != null) {
                    if (maxXDGbasis == null || maxXDGbasis.Degree < maxDGbasis.Degree) {
                        maxXDGbasis = new XDGBasis(_LsTrk, maxDGbasis.Degree);
                    }
                }
            }


            XdgAggregationBasis[] mgXdgB;
            if (maxXDGbasis != null)
                mgXdgB = MultigridSequence.Select(aggGrd => new XdgAggregationBasis(maxXDGbasis, aggGrd)).ToArray();
            else
                mgXdgB = null;


            AggregationGridBasis[] mgDgB;
            if (maxDGbasis != null)
                mgDgB = MultigridSequence.Select(aggGrd => new AggregationGridBasis(maxDGbasis, aggGrd)).ToArray();
            else
                mgDgB = null;

            for (int iLevel = 0; iLevel < MultigridSequence.Length; iLevel++) {
                MultigridBasis[iLevel] = new AggregationGridBasis[BasisS.Length];
                for (int r = 0; r < BasisS.Length; r++) {
                    if (BasisS[r] is XDGBasis || useX)
                        MultigridBasis[iLevel][r] = mgXdgB[iLevel];
                    else
                        MultigridBasis[iLevel][r] = mgDgB[iLevel];
                }
            }

            return MultigridBasis;
        }


        /// <summary>
        /// Initialization of the solver/preconditioner.
        /// </summary>
        protected void InitMultigrid(DGField[] Fields, bool useX) {
            this.MultigridBasis = CreateAggregationGridBasis(
                new UnsetteledCoordinateMapping(Fields.Select(f => f.Basis).ToArray()),
                this.m_LsTrk,
                this.MultigridSequence,
                useX);
        }

        /// <summary>
        /// Agglomerator for the currently set level-set position . 
        /// </summary>
        protected MultiphaseCellAgglomerator m_CurrentAgglomeration;

        /// <summary>
        /// Returns either a solver for the Navier-Stokes or the Stokes system.
        /// E.g. for testing purposes, one might also use a nonlinear solver on a Stokes system.
        /// </summary>
        /// <param name="nonlinSolver"></param>
        /// <param name="linearSolver"></param>
        protected string GetSolver(out NonlinearSolver nonlinSolver, out ISolverSmootherTemplate linearSolver) {
            nonlinSolver = null;
            linearSolver = null;

            string solverDescription = ""; // this.Control.option_solver;
            //solverDescription = string.Format("; No. of cells: {0}, p = {1}/{2}; MaxKrylovDim {3}, MaxIter {4};",
            //    this.GridData.Partitioning.TotalLength,
            //    this.CurrentVel[0].Basis.Degree,
            //    this.Pressure.Basis.Degree,
            //    MaxKrylovDim,
            //    MaxIterations);


            // define solver
            // -------------


            if (this.RequiresNonlinearSolver) {
                // +++++++++++++++++++++++++++++++++++++++++++++
                // the nonlinear solvers:
                // +++++++++++++++++++++++++++++++++++++++++++++

                if (Config_LevelSetHandling != LevelSetHandling.Coupled_Iterative) {
                    switch (Config_NonlinearSolver) {
                        case NonlinearSolverMethod.Picard:

                            nonlinSolver = new FixpointIterator(
                                this.AssembleMatrixCallback,
                                this.MultigridBasis,
                                this.Config_MultigridOperator) {
                                MaxIter = Config_MaxIterations,
                                MinIter = Config_MinIterations,
                                m_LinearSolver = Config_linearSolver,
                                m_SessionPath = SessionPath,
                                ConvCrit = Config_SolverConvergenceCriterion,
                            };
                            break;

                        case NonlinearSolverMethod.Newton:

                            nonlinSolver = new Newton(
                                this.AssembleMatrixCallback,
                                this.MultigridBasis,
                                this.Config_MultigridOperator) {
                                maxKrylovDim = 1000,
                                MaxIter = Config_MaxIterations,
                                MinIter = Config_MinIterations,
                                ApproxJac = Newton.ApproxInvJacobianOptions.GMRES,
                                Precond = Config_linearSolver,
                                GMRESConvCrit = Config_SolverConvergenceCriterion,
                                ConvCrit = Config_SolverConvergenceCriterion,
                                m_SessionPath = SessionPath,
                            };
                            break;

                        default:

                            throw new NotImplementedException();
                    }
                }
                else {
                    switch (Config_NonlinearSolver) {
                        case NonlinearSolverMethod.Picard:

                            nonlinSolver = new FixpointIterator(
                            this.AssembleMatrixCallback,
                            this.MultigridBasis,
                            this.Config_MultigridOperator) {
                                MaxIter = Config_MaxIterations,
                                MinIter = Config_MinIterations,
                                m_LinearSolver = Config_linearSolver,
                                ConvCrit = Config_SolverConvergenceCriterion,
                                LastIteration_Converged = LevelSetConvergenceReached,
                            };
                            break;

                        case NonlinearSolverMethod.Newton:

                            nonlinSolver = new Newton(
                                this.AssembleMatrixCallback,
                                this.MultigridBasis,
                                this.Config_MultigridOperator) {
                                maxKrylovDim = 1000,
                                MaxIter = Config_MaxIterations,
                                MinIter = Config_MinIterations,
                                ApproxJac = Newton.ApproxInvJacobianOptions.GMRES,
                                Precond = Config_linearSolver,
                                GMRESConvCrit = Config_SolverConvergenceCriterion,
                                ConvCrit = Config_SolverConvergenceCriterion,
                                m_SessionPath = SessionPath,
                            };
                            break;

                        default:

                            throw new NotImplementedException();
                    }
                }

            }
            else {

                // +++++++++++++++++++++++++++++++++++++++++++++
                // the linear solvers:
                // +++++++++++++++++++++++++++++++++++++++++++++


                linearSolver = Config_linearSolver;
            }



            // set callback for diagnostic output
            // ----------------------------------

            if (nonlinSolver != null) {
                nonlinSolver.IterationCallback += this.LogResis;
            }

            if (linearSolver != null && linearSolver is ISolverWithCallback) {
                ((ISolverWithCallback)linearSolver).IterationCallback = this.LogResis;
            }

            // return
            // ------


            return solverDescription;
        }


        /// <summary>
        /// Configuration for residual logging (provisional), see <see cref="LogResis(int, double[], double[], MultigridOperator)"/>.
        /// </summary>
        public ResidualLogger m_ResLogger;

        /// <summary>
        /// Configuration for residual logging (provisional), see <see cref="LogResis(int, double[], double[], MultigridOperator)"/>.
        /// 
        /// Names for the residual of each variable.
        /// </summary>
        public string[] m_ResidualNames;

        /// <summary>
        /// Configuration for residual logging (provisional), see <see cref="LogResis(int, double[], double[], MultigridOperator)"/>.
        /// 
        /// If true, the residual will we transformed back to the original XDG basis (before agglomeration and block preconditioning)
        /// before the L2-norm is computed.
        /// </summary>
        public bool m_TransformedResi = true;

        public double m_LastLevelSetResidual;

        protected bool LevelSetConvergenceReached(double Residual_Solver) {

            return (Residual_Solver < Config_SolverConvergenceCriterion && m_LastLevelSetResidual < Config_LevelSetConvergenceCriterion);
        }

        /// <summary>
        /// Logging of residuals (provisional).
        /// </summary>
        virtual protected void LogResis(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {

            if (m_ResLogger != null) {
                int NF = this.CurrentStateMapping.Fields.Count;
                m_ResLogger.IterationCounter = iterIndex;

                if (m_TransformedResi) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // transform current solution and residual back to the DG domain
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    //var X = new CoordinateVector(this.CurrentStateMapping.Fields.Select(f => f.CloneAs()).ToArray());
                    var R = this.Residuals;
                    R.Clear();

                    Mgop.TransformRhsFrom(R, currentRes);
                    //Mgop.TransformSolFrom(X, currentSol);
                    //this.m_Agglomerator.Extrapolate(X, X.Mapping);
                    this.m_CurrentAgglomeration.Extrapolate(R.Mapping);

                    for (int i = 0; i < NF; i++) {
                        double L2Res = R.Mapping.Fields[i].L2Norm();
                        m_ResLogger.CustomValue(L2Res, m_ResidualNames[i]);
                    }
                }
                else {

                    // +++++++++++++++++++++++
                    // un-transformed residual
                    // +++++++++++++++++++++++

                    var VarIdx = NF.ForLoop(i => Mgop.Mapping.GetSubvectorIndices(i));

                    for (int i = 0; i < VarIdx.Length; i++) {
                        double L2Res = 0.0;

                        foreach (int idx in VarIdx[i])
                            L2Res += currentRes[idx].Pow2();
                        L2Res = L2Res.MPISum().Sqrt(); // would be better to do the MPISum for all L2Res together,
                                                       //                                but this implementation is anyway inefficient....

                        m_ResLogger.CustomValue(L2Res, m_ResidualNames[i]);
                    }

                }

                if (Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                    m_ResLogger.CustomValue(m_LastLevelSetResidual, "LevelSet");
                }

                m_ResLogger.NextIteration(true);
            }
        }

        /// <summary>
        /// Sets <see cref="Config_MultigridOperator"/> to a default configuration.
        /// </summary>
        protected void SetConfig_MultigridOperator_Default(IEnumerable<DGField> Fields) {
            int NF = Fields.Count();

            // set the MultigridOperator configuration for each level:
            // it is not necessary to have exactly as many configurations as actual multigrid levels:
            // the last configuration entry will be used for all higher level
            MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
            for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[NF];

                // configurations for velocity
                for (int d = 0; d < NF; d++) {
                    int p = Fields.ElementAt(d).Basis.Degree;

                    configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, p - iLevel),
                        mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                        VarIndex = new int[] { d }
                    };
                }
            }

            this.Config_MultigridOperator = configs;
        }


        /// <summary>
        /// Configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        public virtual MultigridOperator.ChangeOfBasisConfig[][] Config_MultigridOperator {
            get;
            protected set;
        }

        /// <summary>
        /// Coordinate mapping of the current solution.
        /// </summary>
        abstract protected CoordinateMapping CurrentStateMapping {
            get;
        }

        /// <summary>
        /// Callback-routine  to update the linear resp. linearized system, 
        /// see <see cref="AssembleMatrixDel"/> resp. <see cref="NonlinearSolver.m_AssembleMatrix"/>.
        /// </summary>
        /// <param name="argCurSt">Input, current state of solution.</param>
        /// <param name="System">Output.</param>
        /// <param name="Affine">Output.</param>
        /// <param name="MassMatrix">
        /// Mass matrix including agglomeration, without any scaling,
        /// required for block-precond.
        /// </param>
        abstract protected void AssembleMatrixCallback(out BlockMsrMatrix System, out double[] Affine, out BlockMsrMatrix MassMatrix, DGField[] argCurSt);
    }
}
