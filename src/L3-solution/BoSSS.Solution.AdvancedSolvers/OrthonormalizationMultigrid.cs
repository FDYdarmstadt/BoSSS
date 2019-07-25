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
            if (PreSmoother != null)
                PreSmoother.Init(op);
            if (PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                PostSmoother.Init(op);
        }


        public ISolverSmootherTemplate CoarserLevelSolver;
        public ISolverSmootherTemplate PreSmoother;
        public ISolverSmootherTemplate PostSmoother;

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }


        /// <summary>
        /// Threshold for convergence detection
        /// </summary>
        public double Tolerance = 1E-10;


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


        List<double[]> SolHistory = new List<double[]>();
        List<double[]> MxxHistory = new List<double[]>();

        
        void AddSol(ref double[] X) {

            var Xlo = X.CloneAs();
            var Xhi = X.CloneAs();

            int J = m_MgOperator.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int Np = m_MgOperator.Mapping.AggBasis[0].MaximalLength;
            if (J * Np != X.Length)
                throw new NotSupportedException("experimental stuff failed");

            int Ncut = 4;
            for(int j = 0; j < J; j++) {
                for(int n = 0; n < Np; n++) {
                    int l = m_MgOperator.Mapping.LocalUniqueIndex(0, j, n);
                    if (n < Ncut)
                        Xhi[l] = 0;
                    else
                        Xlo[l] = 0;
                }
            }

            __AddSol(ref Xlo);
            __AddSol(ref Xhi);

        }
        

        void __AddSol(ref double[] X) {
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

        double MinimizeResidual<U,V,W,Q>(U outX, V Sol0, W Res0, Q outRes)
            where U : IList<double> //
            where V : IList<double> //
            where W : IList<double> //
            where Q : IList<double> //
        {
            Debug.Assert(SolHistory.Count == MxxHistory.Count);
            int KrylovDim = SolHistory.Count;

            double[] alpha = new double[KrylovDim];
            for (int i = 0; i < KrylovDim; i++) {
                alpha[i] = GenericBlas.InnerProd(MxxHistory[i], Res0).MPISum();
            }

            outX.SetV(Sol0);
            outRes.SetV(Res0);
            for (int i = 0; i < KrylovDim; i++) {
                outX.AccV(alpha[i], SolHistory[i]);
                outRes.AccV(-alpha[i], MxxHistory[i]);
            }

            double ResNorm = outRes.L2NormPow2().MPISum().Sqrt();

            /*
            if(ResNorm < Tolerance) {
                alpha.SaveToStream(Console.Out, m_MgOperator.BaseGridProblemMapping.MPI_Comm);

                alpha.SaveToTextFile("ConvScheisse.txt", m_MgOperator.BaseGridProblemMapping.MPI_Comm);
            }
            */

            return ResNorm;

            /* we cannot do the following 
            // since the 'MxxHistory' vectors form an orthonormal system,
            // the L2-norm is the L2-Norm of the 'alpha'-coordinates (Parceval's equality)
            return alpha.L2Norm();            
            */
        }

        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="_xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="B">the right-hand-side of the problem</param>
        public void Solve<U, V>(U _xl, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (new FuncTrace()) {

                int L = _xl.Count;
                int Lc = m_MgOperator.CoarserLevel.Mapping.LocalLength; // RestrictionOperator.RowPartitioning.LocalLength;
                double[] rl = new double[L];
                double[] rlc = new double[Lc];

                
                if (this.IterationCallback != null) {
                    var _bl = B.ToArray();
                    this.OpMatrix.SpMV(-1.0, _xl, 1.0, _bl);
                    this.IterationCallback(0, _xl.ToArray(), _bl, this.m_MgOperator);
                }
                
                double[] Sol0 = _xl.ToArray();
                double[] Res0 = new double[L];
                Residual(Res0, Sol0, B);

                for (int iIter = 0; iIter < this.m_MaxIterations; iIter++) {

                    // pre-smoother
                    // ------------

                    {

                        Residual(rl, _xl, B); // Residual on this level

                        // compute correction
                        double[] PreCorr = new double[L];
                        PreSmoother.Solve(PreCorr, rl); // Vorglättung

                        // orthonormalization and residual minimization
                        AddSol(ref PreCorr);
                        double resNorm = MinimizeResidual(_xl, Sol0, Res0, rl);
                        if(resNorm < this.Tolerance) {
                            Converged = true;
                            return;
                        }
                    }


                    // coarse grid correction
                    // ----------------------
                    {
                        Residual(rl, _xl, B); // Residual on this level
                        this.m_MgOperator.CoarserLevel.Restrict(rl, rlc);

                        // Berechnung der Grobgitterkorrektur
                        double[] vlc = new double[Lc];
                        this.CoarserLevelSolver.Solve(vlc, rlc);

                        // Prolongation der Grobgitterkorrektur
                        double[] vl = new double[L];
                        this.m_MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);

                        // orthonormalization and residual minimization
                        AddSol(ref vl);
                        double resNorm = MinimizeResidual(_xl, Sol0, Res0, rl);
                        if (resNorm < this.Tolerance) {
                            Converged = true;
                            return;
                        }

                    }

                    // post-smoother
                    // -------------
                    
                    {
                        Residual(rl, _xl, B); // Residual on this level

                        // compute correction
                        double[] PreCorr = new double[L];
                        PostSmoother.Solve(PreCorr, rl); // Vorglättung

                        // orthonormalization and residual minimization
                        AddSol(ref PreCorr);
                        double resNorm = MinimizeResidual(_xl, Sol0, Res0, rl);
                        if (resNorm < this.Tolerance) {
                            Converged = true;
                            return;
                        }
                    }
                    
                    // iteration callback
                    // ------------------

                    this.ThisLevelIterations++;

                    if (this.IterationCallback != null) {
                        var _bl = B.ToArray();
                        this.OpMatrix.SpMV(-1.0, _xl, 1.0, _bl);
                        this.IterationCallback(iIter + 1, _xl.ToArray(), _bl, this.m_MgOperator);
                    }
                    
                }
            }
        }

        public int IterationsInNested {
            get {
                int iter = 0;

                if (PreSmoother != null)
                    iter += this.PreSmoother.IterationsInNested + this.PreSmoother.ThisLevelIterations;

                if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                    iter += this.PostSmoother.IterationsInNested + this.PostSmoother.ThisLevelIterations;

                iter += this.CoarserLevelSolver.IterationsInNested + this.CoarserLevelSolver.ThisLevelIterations;

                return iter;
            }
        }

        public int ThisLevelIterations {
            get;
            private set;
        }

        public bool Converged {
            get;
            private set;
        }

        public void ResetStat() {
            this.Converged = false;
            this.ThisLevelIterations = 0;
            if (this.PreSmoother != null)
                this.PreSmoother.ResetStat();
            if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                this.PostSmoother.ResetStat();
            if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.ResetStat();
        }
        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }
    }
}
