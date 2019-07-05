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


    public class OrthonormalizationMultigrid : ISolverSmootherTemplate, ISolverWithCallback {

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

            if (!Mtx.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!Mtx.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");

            MxxHistory.Clear();
            SolHistory.Clear();

            double Dim = MgMap.ProblemMapping.GridDat.SpatialDimension;

            // set operator
            // ============
            if (op.CoarserLevel == null) {
                throw new NotSupportedException("Multigrid algorithm cannot be used as a solver on the finest level.");
            }
            this.OpMatrix = Mtx;


            // initiate coarser level
            // ======================
            if (this.CoarserLevelSolver == null)
                throw new NotSupportedException("Missing coarse level solver.");
            this.CoarserLevelSolver.Init(op.CoarserLevel);


            // init smoother
            // =============
            if (this.Smoother != null)
                this.Smoother.Init(op);
        }


        public ISolverSmootherTemplate CoarserLevelSolver;
        public ISolverSmootherTemplate Smoother;

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }


        public double m_Tolerance = 0.0;


        /// <summary>
        /// computes the residual on this level.
        /// </summary>
        public void Residual<V1, V2, V3>(V1 rl, V2 xl, V3 bl)
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double> {
            OpMatrix.SpMV(-1.0, xl, 0.0, rl);
            rl.AccV(1.0, bl);
        }


        public int m_MaxIterations = 1;

        bool m_converged = false;
        int NoOfIterations = 0;


        List<double[]> SolHistory = new List<double[]>();
        List<double[]> MxxHistory = new List<double[]>();


        void AddSol(ref double[] X) {
            Debug.Assert(SolHistory.Count == MxxHistory.Count);
            Debug.Assert(X.Length == OpMatrix._RowPartitioning.LocalLength);
            int L = X.Length;
            int KrylovDim = SolHistory.Count;

            double[] Mxx = new double[L];
            OpMatrix.SpMV(1.0, X, 0.0, Mxx);

            for (int i = 0; i < KrylovDim; i++) {
                Debug.Assert(!object.ReferenceEquals(Mxx, MxxHistory[i]));
                double beta = GenericBlas.InnerProd(Mxx, MxxHistory[i]).MPISum();
                Mxx.AccV(-beta, MxxHistory[i]);
                X.AccV(-beta, SolHistory[i]);
            }

            double gamma = 1.0 / GenericBlas.L2NormPow2(Mxx).MPISum().Sqrt();
            Mxx.ScaleV(gamma);
            X.ScaleV(gamma);

            SolHistory.Add(X);
            MxxHistory.Add(Mxx);
            X = null;
        }

        void MinimizeResidual<U,V>(U outX, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            Debug.Assert(SolHistory.Count == MxxHistory.Count);
            int KrylovDim = SolHistory.Count;

            double[] alpha = new double[KrylovDim];
            for (int i = 0; i < KrylovDim; i++) {
                alpha[i] = GenericBlas.InnerProd(MxxHistory[i], B).MPISum();
            }

            outX.ClearEntries();
            for (int i = 0; i < KrylovDim; i++)
                outX.AccV(alpha[i], SolHistory[i]);

        }

        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="_xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="bl">the right-hand-side of the problem</param>
        public void Solve<U, V>(U _xl, V bl)
            where U : IList<double>
            where V : IList<double> //
        {
            using (new FuncTrace()) {

                int L = _xl.Count;
                int Lc = m_MgOperator.CoarserLevel.Mapping.LocalLength; // RestrictionOperator.RowPartitioning.LocalLength;
                double[] rl = new double[L];
                double[] rlc = new double[Lc];

                /*
                if (this.IterationCallback != null) {
                    var _bl = bl.ToArray();
                    this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                    this.IterationCallback(0, xl.ToArray(), _bl, this.m_MgOperator);
                }
                */

                {
                    double[] Sol0 = _xl.ToArray();
                    AddSol(ref Sol0);
                }


                for (int iIter = 0; iIter < this.m_MaxIterations; iIter++) {

                    // pre-smoother
                    // ------------

                    {

                        Residual(rl, _xl, bl); // Residual on this level

                        // compute correction
                        double[] PreCorr = new double[L];
                        Smoother.Solve(PreCorr, rl); // Vorglättung

                        // orthonormalization and residual minimization
                        AddSol(ref PreCorr);
                        MinimizeResidual(_xl, bl);
                    }


                    // coarse grid correction
                    // ======================
                    {
                        Residual(rl, _xl, bl); // Residual on this level
                        this.m_MgOperator.CoarserLevel.Restrict(rl, rlc);

                        // Berechnung der Grobgitterkorrektur
                        double[] vlc = new double[Lc];
                        this.CoarserLevelSolver.Solve(vlc, rlc);

                        // Prolongation der Grobgitterkorrektur
                        double[] vl = new double[L];
                        this.m_MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);

                        // orthonormalization and residual minimization
                        AddSol(ref vl);
                        MinimizeResidual(_xl, bl);

                    }

                    // post-smoother
                    // -------------

                    {

                        Residual(rl, _xl, bl); // Residual on this level

                        // compute correction
                        double[] PreCorr = new double[L];
                        Smoother.Solve(PreCorr, rl); // Vorglättung

                        // orthonormalization and residual minimization
                        AddSol(ref PreCorr);
                        MinimizeResidual(_xl, bl);
                    }

                    /*
                    if (this.IterationCallback != null) {
                        var _bl = bl.ToArray();
                        this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                        this.IterationCallback(iIter + 1, xl.ToArray(), _bl, this.m_MgOperator);
                    }
                    */
                }
            }
        }

        public int IterationsInNested {
            get {
                return
                    ((this.Smoother != null) ? (this.Smoother.IterationsInNested + this.Smoother.ThisLevelIterations) : 0)
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
            if (this.Smoother != null)
                this.Smoother.ResetStat();
            if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.ResetStat();
        }
        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }
    }
}
