using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MathNet.Numerics.Algorithms.LinearAlgebra;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {
    static class Utils {

        public static double rhoDinvA(MsrMatrix A, out MsrMatrix DinvA)
        {
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

            for (int i = 0; i < n; i++)
            {
                Lr = A.GetRow(i + i0, ref ColumnIdx, ref Values);
                for (int k = 0; k < Lr; k++)
                {
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
        static public double MaxAbsEigen(double[,] H)
        {
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
        public static void arnoldi(out double[][] V, out double[,] H, out int kact, MsrMatrix A, double[] v0, int k, bool reorth = false)
        {
            //%
            //%             -----  kact=k -------
            //%      V      n-by-(k+1)  Arnoldi vectors
            //%      H      (k+1)-by-k
            //%             -----  kact=j<k -------
            //%      V      n-by-j  Arnoldi vectors
            //%      H      j-by-j

            double eps = BLAS.MachineEps;


            int n = A.RowPartitioning.LocalLength;
            if (A.ColPartition.LocalLength != A.RowPartitioning.LocalLength)
            {
                throw new ArgumentException("the sizes of input matrix incorrect");
            }

            V = (k + 1).ForLoop(i => new double[n]);
            H = new double[k + 1, k];

            double nrm2 = v0.L2NormPow2().MPISum().Sqrt();
            if (nrm2 == 0.0)
            {
                throw new ArgumentException("arnoldi: input v0 is a zero vector");
            }

            double tol = n * eps;

            V[0].SetV(v0, 1 / nrm2);   //v(:,1)=v0/nrm2;
            for (int j = 0; j < k; j++)
            {
                double[] vh = new double[n];
                A.SpMVpara(1.0, V[j], 0.0, vh);    //vh = A*V(:,j);   
                double nrmvh = vh.L2NormPow2().MPISum().Sqrt();

                //%   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                //%   by MGS
                for (int i = 0; i < j; i++)
                {
                    double hij = GenericBlas.InnerProd(V[i], vh).MPISum();
                    vh.AccV(-hij, V[i]); //vh = vh - hij*V(:,i);
                    H[i, j] = hij;
                }
                if (reorth)
                {
                    for (int i = 0; i < j; i++)
                    {
                        double tmp = GenericBlas.InnerProd(V[i], vh).MPISum();
                        vh.AccV(-tmp, V[i]);  //vh = vh - tmp*V(:,i);
                        H[i, j] = H[i, j] + tmp;
                    }
                }
                //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                H[j + 1, j] = vh.L2NormPow2().MPISum().Sqrt();
                V[j + 1].SetV(vh, 1.0 / H[j + 1, j]);

                if (H[j + 1, j] <= tol * nrmvh)
                {
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
