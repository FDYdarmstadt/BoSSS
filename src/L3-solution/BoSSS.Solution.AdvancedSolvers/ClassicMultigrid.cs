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
using MathNet.Numerics.Algorithms.LinearAlgebra;
using System.Numerics;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {


    public class ClassicMultigrid : ISolverSmootherTemplate, ISolverWithCallback {


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
        public BlockMsrMatrix OpMatrix;


        MultigridOperator m_MgOperator;

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(MultigridOperator op) {
            this.m_MgOperator = op;
            var Mtx = op.OperatorMatrix;
            var MgMap = op.Mapping;

            if(!Mtx.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if(!Mtx.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");


            double Dim = MgMap.ProblemMapping.GridDat.SpatialDimension;

            // set operator
            // ============
            if(op.CoarserLevel == null) {
                throw new NotSupportedException("Multigrid algorithm cannot be used as a solver on the finest level.");
            }
            this.OpMatrix = Mtx;

            
            /*
            // construct restriction and prolongation
            // ======================================

            var thisLevel = MgMap.AggBasis;
            var lowerLevel = MgMap.CoarserLevel.AggBasis;

            int[] this_p = MgMap.DgDegree;
            int[] low_p = MgMap.CoarserLevel.DgDegree;

            var __PrlgOperator = thisLevel.FromOtherLevelMatrix(MgMap.ProblemMapping, this_p, lowerLevel, low_p);
            Debug.Assert(__PrlgOperator.RowPartitioning.LocalLength == MgMap.LocalLength);
            Debug.Assert(__PrlgOperator.ColPartition.LocalLength == MgMap.CoarserLevel.LocalLength);
#if DEBUG
            var __RestOperator = lowerLevel.FromOtherLevelMatrix(MgMap.ProblemMapping, low_p, thisLevel, this_p);
            {
                var __RestOperator2 = __PrlgOperator.Transpose();
                __RestOperator2.Acc(-1.0, __RestOperator);
                double dingsNorm = __RestOperator2.InfNorm();
                Debug.Assert(dingsNorm <= 1.0e-10);
            }
#else 
            var __RestOperator = __PrlgOperator.Transpose();
#endif
            int iLevel = MgMap.LevelIndex;

            if(this.UseSmoothedInterpolation) {
                MsrMatrix DinvA;
                double specRadius = Utils.rhoDinvA(OpMatrix, out DinvA);
                Console.WriteLine("Level {0} spectral radius = {1}", iLevel, specRadius);
                double omega = ((Dim + 1) / Dim) * (1.0 / specRadius);

                DinvA.Scale(-omega);
                //DinvA.Clear();
                DinvA.AccEyeSp(1.0);

                PrologateOperator = MsrMatrix.Multiply(DinvA, __PrlgOperator);
                RestrictionOperator = this.PrologateOperator.Transpose();
            } else {
                PrologateOperator = __PrlgOperator;
                RestrictionOperator = __RestOperator;

                var RX = this.RestrictionOperator.CloneAs();
                RX.Acc(-1.0, m_MgOperator.RestrictionOperator);
                double tst = RX.InfNorm();
                Debug.Assert(tst <= 1.0e-10);

                var PX = this.PrologateOperator.CloneAs();
                PX.Acc(-1.0, m_MgOperator.PrologateOperator);
                double tst2 = PX.InfNorm();
                Debug.Assert(tst2 <= 1.0e-10);
            }
            //*/

            // initiate coarser level
            // ======================
            //var M1 = MsrMatrix.Multiply(this.OpMatrix, this.PrologateOperator);
            //var CoarseMtx = MsrMatrix.Multiply(this.RestrictionOperator, M1);
            if(this.CoarserLevelSolver == null)
                throw new NotSupportedException("Missing coarse level solver.");
            this.CoarserLevelSolver.Init(op.CoarserLevel);


            // init pre&post smoother
            // ======================
            if(this.PreSmoother != null)
                this.PreSmoother.Init(op);
            if(this.PostSmoother != null)
                this.PostSmoother.Init(op);
        }

        
        public ISolverSmootherTemplate CoarserLevelSolver;
        public ISolverSmootherTemplate PreSmoother;
        public ISolverSmootherTemplate PostSmoother;

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }


        public double m_Tolerance = 0.0;


        /// <summary>
        /// computes the residual on this level.
        /// </summary>
        public void Residual<V1,V2,V3>(V1 rl, V2 xl, V3 bl) 
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double>
        {
            OpMatrix.SpMV(-1.0, xl, 0.0, rl);
            rl.AccV(1.0, bl);
        }


        int m_Gamma = 1;


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

        public int m_MaxIterations = 1;

        bool m_converged = false;
        int NoOfIterations = 0;


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
                int NN = m_MgOperator.CoarserLevel.Mapping.LocalLength; // RestrictionOperator.RowPartitioning.LocalLength;
                double[] rl = new double[N];
                double[] rlp1 = new double[NN];

                for (int iIter = 0; iIter < this.m_MaxIterations; iIter++) {

                    var DGBasis = m_MgOperator.BaseGridProblemMapping.BasisS[0];
                    //XDGField ResB4Jacobi = new XDGField((XDGBasis)DGBasis, "Resi_B4Jacobi");
                    //m_MgOperator.TransformRhsFrom(ResB4Jacobi.CoordinatesAsVector, bl);

                    if (PreSmoother != null)
                            PreSmoother.Solve(xl, bl); // Vorglättung
                  
                    Residual(rl, xl, bl); // Residual on this level


                    if (this.IterationCallback != null)
                    {
                        //var _bl = bl.ToArray();
                        //this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                        this.IterationCallback(0, xl.ToArray(), rl, this.m_MgOperator);
                    }

                    //var ResAftJacobi = new XDGField((XDGBasis)DGBasis, "Resi_afJacobi");
                    //m_MgOperator.TransformRhsFrom(ResAftJacobi.CoordinatesAsVector, rl);

                    //Tecplot.Tecplot.PlotFields(new DGField[] { ResB4Jacobi, ResAftJacobi }, DGBasis.GridDat, "Resi", "Resi", 0, 4);

                    if (this.m_Tolerance > 0) {
                        double ResNorm = bl.L2NormPow2().MPISum().Sqrt();
                        if (ResNorm < this.m_Tolerance) {
                            m_converged = true;
                            return;
                        }
                    }
                    this.NoOfIterations++;

                    this.m_MgOperator.CoarserLevel.Restrict(rl, rlp1);

                    // Berechnung der Grobgitterkorrektur
                    double[] vlp1 = new double[NN];

                    for (int j = 0; j < m_Gamma; j++) {
                        this.CoarserLevelSolver.Solve(vlp1, rlp1);
                    }

                    // Prolongation der Grobgitterkorrektur
                    this.m_MgOperator.CoarserLevel.Prolongate(1.0, xl, 1.0, vlp1);

                    // check resudual after coarse grid correction
                    Residual(rl, xl, bl); // Residual on this level

                    
                    if (this.IterationCallback != null)
                    {
                        //var _bl = bl.ToArray();
                        //this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                        this.IterationCallback(iIter, xl.ToArray(), rl, this.m_MgOperator);
                    }

                    // Nachglättung
                    if (PostSmoother != null)
                        PostSmoother.Solve(xl, bl);

                }
            }
        }

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
        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }
    }

}
