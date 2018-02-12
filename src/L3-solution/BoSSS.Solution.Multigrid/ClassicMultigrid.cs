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

namespace BoSSS.Solution.Multigrid {


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

        static public ISolverSmootherTemplate InitMultigridChain(IEnumerable<AggregationGrid> MultigridSequence,
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
            where V : IList<double> 
        {
            //

            int N = xl.Count;
            int NN = m_MgOperator.CoarserLevel.Mapping.LocalLength; // RestrictionOperator.RowPartitioning.LocalLength;
            double[] rl = new double[N];
            double[] rlp1 = new double[NN];

            if(this.IterationCallback != null) {
                var _bl = bl.ToArray();
                this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                this.IterationCallback(0, xl.ToArray(), _bl, this.m_MgOperator);
            }

            for(int iIter = 0; iIter < this.m_MaxIterations; iIter++) {

                var DGBasis = m_MgOperator.BaseGridProblemMapping.BasisS[0];
                //XDGField ResB4Jacobi = new XDGField((XDGBasis)DGBasis, "Resi_B4Jacobi");
                //m_MgOperator.TransformRhsFrom(ResB4Jacobi.CoordinatesAsVector, bl);

                if(PreSmoother != null)
                    PreSmoother.Solve(xl, bl); // Vorglättung
                Residual(rl, xl, bl); // Residual on this level

                //var ResAftJacobi = new XDGField((XDGBasis)DGBasis, "Resi_afJacobi");
                //m_MgOperator.TransformRhsFrom(ResAftJacobi.CoordinatesAsVector, rl);

                //Tecplot.Tecplot.PlotFields(new DGField[] { ResB4Jacobi, ResAftJacobi }, DGBasis.GridDat, "Resi", "Resi", 0, 4);

                if(this.m_Tolerance > 0) {
                    double ResNorm = bl.L2NormPow2().MPISum().Sqrt();
                    if(ResNorm < this.m_Tolerance) {
                        m_converged = true;
                        return;
                    }
                }
                this.NoOfIterations++;

                this.m_MgOperator.CoarserLevel.Restrict(rl, rlp1);

                // Berechnung der Grobgitterkorrektur
                double[] vlp1 = new double[NN];

                for(int j = 0; j < m_Gamma; j++) {
                    this.CoarserLevelSolver.Solve(vlp1, rlp1);
                }

                // Prolongation der Grobgitterkorrektur
                this.m_MgOperator.CoarserLevel.Prolongate(1.0, xl, 1.0, vlp1);

                
                // Nachglättung
                if(PostSmoother != null)
                    PostSmoother.Solve(xl, bl);

                if(this.IterationCallback != null) {
                    var _bl = bl.ToArray();
                    this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                    this.IterationCallback(iIter + 1, xl.ToArray(), _bl, this.m_MgOperator);
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
    }



    static class Utils {

        public static double rhoDinvA(MsrMatrix A, out MsrMatrix DinvA) {
            double rho;

            // extract diagonal matrix
            double[] Diag = A.GetDiagVector();
            int n = Diag.Length;

            // % form (D^-1) A
            //Diag = spdiags(1./Diag, 0, n, n);
            //DinvA = Diag*A;
            DinvA = new MsrMatrix(A.RowPartitioning, A.ColPartition);
            int i0 = A.RowPartitioning.i0;

            int Lr;
            int[] ColumnIdx = null;
            double[] Values = null;

            for (int i = 0; i < n; i++) {
                Lr = A.GetRow(i + i0, ref ColumnIdx, ref Values);
                for(int k = 0; k < Lr; k++) {
                    Values[k] *= 1.0 / Diag[i];
                }
                DinvA.SetRow(i + i0, ColumnIdx, Values, Lr);
            }
#if DEBUG
            for(int i = 0; i < n; i++) {
                Debug.Assert(Math.Abs(1.0 - DinvA[i + i0, i + i0]) < BLAS.MachineEps * 10);
            }
#endif

            // estimate the largest eigen value from Arnoldi iteration
            int kk = 20;
            var rand = new Random(0);
            double[] v0 = n.ForLoop(i => rand.NextDouble());
            double[][] V; double[,] H; int kact;
            arnoldi(out V, out H, out kact, DinvA, v0, kk, false);
            kk = Math.Min(H.GetLength(0), H.GetLength(1));
            H = H.GetSubMatrix(0, kk, 0, kk);

            rho = MaxAbsEigen(H);

            return rho;
        }

        /// <summary>
        /// Maximum of the absolute value of all Eigenvalues of <paramref name="H"/>.
        /// </summary>
        static public double MaxAbsEigen(double[,] H) {
            var linalg = new ManagedLinearAlgebraProvider();

            double Eigen;
            int N = H.GetLength(0);
            double[] Matrix = H.Resize(false);
            double[] EigenVect = new double[Matrix.Length];
            double[] diagM = new double[Matrix.Length];
            System.Numerics.Complex[] EigenValues = new System.Numerics.Complex[N];
            linalg.EigenDecomp(false, N, Matrix, EigenVect, EigenValues, diagM);
            Eigen = EigenValues.Select(ev => Complex.Abs(ev)).Max();
            return Eigen;
        }

        /// <summary>
        /// Arnoldi iteration 
        /// </summary>
        /// <param name="V">Output: Arnoldi vectors</param>
        /// <param name="H">Output: </param>
        /// <param name="kact">Output:</param>
        /// <param name="A">Input: (n-by-n) the matrix </param>
        /// <param name="v0">Input: n-vector</param>
        /// <param name="k">Input: number of Arnoldi steps requested</param>
        /// <param name="reorth">Input: (optional) set to 1 for reorthogonalization, (default), set to any other value to switch it off</param>
        /// <remarks>
        /// (c) Ren-Cang Li, rcli@uta.edu,  06/16/07
        /// </remarks>
        public static void arnoldi(out double[][] V, out double[,] H, out int kact, MsrMatrix A, double[] v0, int k, bool reorth = false) {
            //%
            //%             -----  kact=k -------
            //%      V      n-by-(k+1)  Arnoldi vectors
            //%      H      (k+1)-by-k
            //%             -----  kact=j<k -------
            //%      V      n-by-j  Arnoldi vectors
            //%      H      j-by-j

            double eps = BLAS.MachineEps;


            int n = A.RowPartitioning.LocalLength;
            if(A.ColPartition.LocalLength != A.RowPartitioning.LocalLength) {
                throw new ArgumentException("the sizes of input matrix incorrect");
            }

            V = (k + 1).ForLoop(i => new double[n]);
            H = new double[k + 1, k];

            double nrm2 = v0.L2NormPow2().MPISum().Sqrt();
            if(nrm2 == 0.0) {
                throw new ArgumentException("arnoldi: input v0 is a zero vector");
            }

            double tol = n * eps;

            V[0].SetV(v0, 1 / nrm2);   //v(:,1)=v0/nrm2;
            for(int j= 0; j < k; j++) {
                double[] vh = new double[n];
                A.SpMVpara(1.0, V[j], 0.0, vh);    //vh = A*V(:,j);   
                double nrmvh = vh.L2NormPow2().MPISum().Sqrt();

                //%   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                //%   by MGS
                for(int i= 0; i < j; i++) {
                    double hij = GenericBlas.InnerProd(V[i], vh).MPISum();
                    vh.AccV(-hij, V[i]); //vh = vh - hij*V(:,i);
                    H[i, j] = hij;
                }
                if(reorth) {
                    for(int i =0; i < j; i++) {
                        double tmp = GenericBlas.InnerProd(V[i], vh).MPISum();
                        vh.AccV(-tmp, V[i]);  //vh = vh - tmp*V(:,i);
                        H[i, j] = H[i, j] + tmp;
                    }
                }
                //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                H[j + 1, j] = vh.L2NormPow2().MPISum().Sqrt();
                V[j + 1].SetV(vh, 1.0 / H[j + 1, j]);

                if(H[j + 1, j] <= tol * nrmvh) {
                    //%             -----  kact<k -------
                    //%      V      n    -by- kact            Arnoldi vectors
                    //%      H      kact -by- kact
                    // Console.WriteLine("termination at step: " + j);
                    kact = j + 1;
                    V = V.GetSubVector(0, kact);
                    H = H.GetSubMatrix(0, kact, 0, kact);
                    return;
                }
            }
            kact = k;
            Debug.Assert(V.Length == kact + 1);
            Debug.Assert(V.Length == kact + 1);

            //%             -----  kact=k -------
            //%      V       n        -by-  (kact+1)  Arnoldi vectors
            //%      H      (kact+1)  -by-   kact
        }
    }

}
