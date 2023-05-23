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
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using MPI.Wrappers;
using System.Numerics;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using System.Runtime.Serialization;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Multigrid methods for linear systems, using pre- and post-smoother (<see cref="PreSmoother"/>, <see cref="PostSmoother"/>) as well as coarse-grid correction (<see cref="CoarserLevelSolver"/>).
    /// </summary>
    public class ClassicMultigrid : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination, ISubsystemSolver {


        /// <summary>
        /// Individual configuration of <see cref="OrthonormalizationMultigrid"/>
        /// </summary>
        [DataContract]
        [Serializable]
        public class Config : ISolverFactory {


            /// <summary>
            /// Ctor
            /// </summary>
            public Config() {
            }



            /// <summary>
            /// - True: the default value: <see cref="ClassicMultigrid.CoarserLevelSolver"/> is initialized and solved on coarser level
            /// - false: <see cref="ClassicMultigrid.CoarserLevelSolver"/> is initialized on the same level, but it may perform tis own restriction
            /// </summary>
            [DataMember]
            public bool CoarseOnLovwerLevel = true;

            /// <summary>
            /// - if set to 1, a this performs a V-cycle
            /// - if set to 2 or higher, a W-cycle, resp. WW-cycle, WWW-cycle, etc.
            /// </summary>
            [DataMember]
            public int m_omega = 1;

            /// <summary>
            /// 
            /// </summary>
            public string Name {
                get { return "ClassicMultigrid"; }
            }

            /// <summary>
            /// 
            /// </summary>
            public string Shortname {
                get { return "MG"; }
            }

            /// <summary>
            /// factory
            /// </summary>
            public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var instance = new ClassicMultigrid();
                instance.myConfig = this;
                instance.Init(level);
                return instance;
            }

            public bool Equals(ISolverFactory other) {
                var othr = other as ClassicMultigrid.Config;
                if (othr == null)
                    return false;

                if (othr.m_omega != this.m_omega)
                    return false;
                if (othr.CoarseOnLovwerLevel != this.CoarseOnLovwerLevel)
                    return false;

                return true;
            }
        }


        Config myConfig = new Config();

        /// <summary>
        /// Solver configuration
        /// </summary>
        public Config config {
            get {
                return myConfig;
            }
        }

        /// <summary>
        /// Creates a V- or W-cycle (depending on <see cref="ClassicMultigrid.Config.m_omega"/>) of multi-grid methods, 
        /// with each multi-grid solver being the coarse-level solver for the finer level.
        /// </summary>
        static public ISolverSmootherTemplate InitMultigridChain(MultigridOperator MgOp,
            Func<int, ISolverSmootherTemplate> PreSmootherFactory,
            Func<int, ISolverSmootherTemplate> PostSmootherFactory,
            Action<int, ClassicMultigrid> ParamsSeter,
            Func<ISolverSmootherTemplate> CoarsestSolverFactory) //
        {
            if(MgOp.CoarserLevel == null) {
                return CoarsestSolverFactory();
            } else {
                var MgTop = new ClassicMultigrid();
                ParamsSeter(MgOp.LevelIndex, MgTop);
                MgTop.PreSmoother = PreSmootherFactory(MgOp.LevelIndex);
                MgTop.PostSmoother = PostSmootherFactory(MgOp.LevelIndex);
                MgTop.CoarserLevelSolver = InitMultigridChain(MgOp.CoarserLevel, PreSmootherFactory, PostSmootherFactory, ParamsSeter, CoarsestSolverFactory);
                return MgTop;
            }
        }

        /// <summary>
        /// Creates a V- or W-cycle (depending on <see cref="ClassicMultigrid.Config.m_omega"/>) of multi-grid methods, 
        /// with each multi-grid solver being the coarse-level solver for the finer level.
        /// </summary>
        static public ISolverSmootherTemplate InitMultigridChain(IEnumerable<AggregationGridData> MultigridSequence,
            Func<int, ISolverSmootherTemplate> PreSmootherFactory,
            Func<int, ISolverSmootherTemplate> PostSmootherFactory,
            Action<int, ClassicMultigrid> ParamsSeter,
            Func<ISolverSmootherTemplate> CoarsestSolverFactory,
            int iLevel = 0) //
        {
            if(iLevel == MultigridSequence.Count() - 1) {
                return CoarsestSolverFactory();
            } else {
                var MgTop = new ClassicMultigrid();
                ParamsSeter(iLevel, MgTop);
                MgTop.PreSmoother = PreSmootherFactory(iLevel);
                MgTop.PostSmoother = PostSmootherFactory(iLevel);
                MgTop.CoarserLevelSolver = InitMultigridChain(MultigridSequence, PreSmootherFactory, PostSmootherFactory, ParamsSeter, CoarsestSolverFactory, iLevel + 1);
                return MgTop;
            }
        }
        
        
        /// <summary>
        /// The matrix at this level.
        /// </summary>
        public BlockMsrMatrix OpMatrix {
            get {
                return m_MgOperator.OperatorMatrix;
            }
        }


        IOperatorMappingPair m_MgOperator;

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        void InitImpl(IOperatorMappingPair op) {
            using(new FuncTrace()) {
                if(object.ReferenceEquals(op, this.m_MgOperator))
                    return; // already initialized
                else
                    this.Dispose();

                this.m_MgOperator = op;
                var Mtx = op.OperatorMatrix;
                var MgMap = op.DgMapping;

                if(!Mtx.RowPartitioning.EqualsPartition(MgMap))
                    throw new ArgumentException("Row partitioning mismatch.");
                if(!Mtx.ColPartition.EqualsPartition(MgMap))
                    throw new ArgumentException("Column partitioning mismatch.");


             
                // set operator
                // ============
                


               




                // initiate coarser level
                // ======================
                
                if (this.CoarserLevelSolver == null) {
                    //throw new NotSupportedException("Missing coarse level solver.");
                    Console.WriteLine("OrthonormalizationMultigrid: running without coarse solver.");
                } else {
                    if (op is MultigridOperator mgOp) {
                        if (myConfig.CoarseOnLovwerLevel && mgOp.CoarserLevel != null) {
                            this.CoarserLevelSolver.Init(mgOp.CoarserLevel);
                        } else {
                            Console.WriteLine("OrthonormalizationMultigrid: running coarse solver on same level.");
                            this.CoarserLevelSolver.Init(mgOp);
                        }
                    } else {
                        if (myConfig.CoarseOnLovwerLevel == false && this.CoarserLevelSolver is ISubsystemSolver ssCoarse) {
                            ssCoarse.Init(op);
                        } else {
                            throw new NotSupportedException($"Unable to initialize coarse-level-solver if operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }

                // init pre&post smoother
                // ======================
                if (PreSmoother != null) {
                    if (PreSmoother is ISubsystemSolver ssPreSmother) {
                        ssPreSmother.Init(op);
                    } else {
                        if (op is MultigridOperator mgOp) {
                            PreSmoother.Init(mgOp);
                        } else {
                            throw new NotSupportedException($"Unable to initialize pre-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }

                if (PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother)) {
                    //PostSmoother.Init(op);
                    if (PostSmoother is ISubsystemSolver ssPostSmother) {
                        ssPostSmother.Init(op);
                    } else {
                        if (op is MultigridOperator mgOp) {
                            PostSmoother.Init(mgOp);
                        } else {
                            throw new NotSupportedException($"Unable to initialize post-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                        }
                    }

                }
            }
        }
        
        public ISolverSmootherTemplate CoarserLevelSolver;
        public ISolverSmootherTemplate PreSmoother;
        public ISolverSmootherTemplate PostSmoother;

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }



        /// <summary>
        /// computes the residual on this level.
        /// </summary>
        public double Residual<V1,V2,V3>(V1 rl, V2 xl, V3 bl) 
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double>
        {
            OpMatrix.SpMV(-1.0, xl, 0.0, rl);
            rl.AccV(1.0, bl);
            return rl.MPI_L2Norm(OpMatrix.MPI_Comm);
        }

        /*

        int m_Gamma = 1;

        /// <summary>
        /// Number of coarse-level corrections:
        /// - 1 leads to a V-cycle
        /// - 2 or more leads to a W-cycle
        /// </summary>
        public int Gamma {
            get {
                return m_Gamma;
            }
            set {
                if(value < 1)
                    throw new ArgumentException();
                m_Gamma = value;
            }
        }
        */

        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
        }


        bool m_converged = false;
        int NoOfIterations = 0;

        /// <summary>
        /// ctor
        /// </summary>
        public ClassicMultigrid() {
            TerminationCriterion = (iIter, r0, ri) => (iIter <= 1, true);
        }


        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="bl">the right-hand-side of the problem</param>
        public void Solve<U, V>(U xl, V bl)
            where U : IList<double>
            where V : IList<double> {
            using (new FuncTrace()) {

                int N = xl.Count;
                
                double[] rl = new double[N];
                

                double iter0ResidualNorm = bl.MPI_L2Norm(m_MgOperator.DgMapping.MPI_Comm);
                double iterNorm = iter0ResidualNorm;

                for (int iIter = 0; true; iIter++) {
                    var term = TerminationCriterion(iIter, iter0ResidualNorm, iterNorm);
                    this.m_converged = term.bSuccess;
                    if (!term.bNotTerminate)
                        return;



                    if (PreSmoother != null)
                        PreSmoother.Solve(xl, bl); // Vorglättung
                  
                    iterNorm = Residual(rl, xl, bl); // Residual on this level

                    if (CoarserLevelSolver != null) {
                        if (config.CoarseOnLovwerLevel) {
                            var _MgOperator = this.m_MgOperator as MultigridOperator;

                            int NN = _MgOperator.CoarserLevel.Mapping.LocalLength;
                            double[] rlp1 = new double[NN];
                            _MgOperator.CoarserLevel.Restrict(rl, rlp1);

                            // Berechnung der Grobgitterkorrektur
                            double[] vlp1 = new double[NN];
                            for (int j = 0; j < config.m_omega; j++) {
                                this.CoarserLevelSolver.Solve(vlp1, rlp1);
                            }

                            // Prolongation der Grobgitterkorrektur
                            _MgOperator.CoarserLevel.Prolongate(1.0, xl, 1.0, vlp1);
                        } else {
                            var ssCoarse = CoarserLevelSolver as ISubsystemSolver;
                            for (int j = 0; j < config.m_omega; j++) {
                                ssCoarse.Solve(xl, bl);
                                if (j < config.m_omega - 1)
                                    iterNorm = Residual(rl, xl, bl);
                            }
                        }
                    }
                    // check residual after coarse grid correction
                    Residual(rl, xl, bl); // Residual on this level

                    // Nachglättung
                    if (PostSmoother != null)
                        PostSmoother.Solve(xl, bl);

                    // update residual
                    this.IterationCallback?.Invoke(iIter + 1, xl.ToArray(), rl, this.m_MgOperator as MultigridOperator);
                    iterNorm = Residual(rl, xl, bl); // Residual on this level
                    this.NoOfIterations++;
                }
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public int IterationsInNested {
            get {
                return
                    ((this.PreSmoother != null) ? (this.PreSmoother.IterationsInNested + this.PreSmoother.ThisLevelIterations) : 0)
                    + ((this.PostSmoother != null) ? (this.PostSmoother.IterationsInNested + this.PostSmoother.ThisLevelIterations) : 0)
                    + (this.CoarserLevelSolver.IterationsInNested + this.CoarserLevelSolver.ThisLevelIterations);
            }
        }

        public int ThisLevelIterations {
            get {
                return this.NoOfIterations;
            }
        }

        public bool Converged {
            get {
                return this.m_converged;
            }
        }

        public void ResetStat() {
            this.m_converged = false;
            this.NoOfIterations = 0;
            if(this.PreSmoother != null)
                this.PreSmoother.ResetStat();
            if(this.PostSmoother != null)
                this.PostSmoother.ResetStat();
            if(this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.ResetStat();
        }
        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

        public void Dispose() {
            if(PreSmoother != null)
                PreSmoother.Dispose();
            if (PostSmoother != null)
                PostSmoother.Dispose();
            if (CoarserLevelSolver != null)
                CoarserLevelSolver.Dispose();
            m_MgOperator = null;
        }

        public long UsedMemory() {
            long ret = 0;
            ret += PreSmoother?.UsedMemory() ?? 0;
            ret += PostSmoother?.UsedMemory() ?? 0;
            ret += CoarserLevelSolver?.UsedMemory() ?? 0;

            return ret;
        }

        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        public void Init(MultigridOperator op) {
            InitImpl(op);
        }
    }

}
