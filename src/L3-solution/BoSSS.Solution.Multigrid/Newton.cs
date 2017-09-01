using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;
using ilPSP.Connectors.Matlab;

namespace BoSSS.Solution.Multigrid {
    public class Newton : NonlinearSolver {
        public Newton(AssembleMatrixDel __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq, MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig) : base(__AssembleMatrix, __AggBasisSeq, __MultigridOperatorConfig)
        {
        }

        public int MaxIter = 10;
        public int MinIter = 1;
        public double ConvCrit = 1e-9;
        public double GMRESConvCrit = 1e-9;

        public ISolverSmootherTemplate m_LinearSolver;

        public CoordinateVector m_SolutionVec;


        public override void SolverDriver<S>(CoordinateVector SolutionVec, S RHS)
        {
            m_SolutionVec = SolutionVec;

            int itc;
            itc = 0;
            double[] x, xt, xOld, f0, deltaX;
            double rat;
            
            // Eval_F0 
            base.Init(SolutionVec, RHS, out x, out f0);

            deltaX = new double[x.Length];
            xt = new double[x.Length];

            // fnorm
            double fnorm = f0.L2NormPow2().MPISum().Sqrt();
            double fNormo = 1;
            double errstep;
            double[] step;

            OnIterationCallback(itc, x.CloneAs(), f0.CloneAs(), this.CurrentLin);

            while (fnorm > ConvCrit && itc < MaxIter) {
                rat = fnorm / fNormo;
                fNormo = fnorm;
                itc++;

                // dKrylov
                step = Krylov(SolutionVec, x, f0, out errstep);

                // Line search starts here
                xOld = x;
                SolutionVec.AccV(1, step);
                double lambda = 1.0;
                base.Init(SolutionVec, RHS, out x, out f0);

                fnorm = f0.L2Norm();      

                OnIterationCallback(itc, x.CloneAs(), f0.CloneAs(), this.CurrentLin);

            }

            SolutionVec = m_SolutionVec;
        }

        public int krylovDim = 1;
        public int GMRESmaxIter = 5;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="SolutionVec">Current Point</param>
        /// <param name="f0">Function at current point</param>
        /// <param name="xinit">initial iterate</param>
        /// <param name="errstep">error of step</param>
        /// <returns></returns>
        public double[] NewtonGMRES(CoordinateVector SolutionVec, double[] f0, double[] xinit, out double errstep)
        {
            int n = f0.Length;
            int reorth = 1;

            // RHS of the linear equation system 
            double[] b = new double[n];
            b.AccV(-1, f0);

            double[] x = new double[n];
            double[] r = new double[n];

            int Nloc = base.CurrentLin.OperatorMatrix.RowPartitioning.LocalLength;
            int Ntot = base.CurrentLin.OperatorMatrix.RowPartitioning.TotalLength;

            r = b;

            //Initial solution
            if (xinit.L2Norm() != 0) {
                x = xinit;
                r = dirder(SolutionVec, x, f0);
                r.ScaleV(-1);
                r.AccV(-1, f0);
            }

            double[][] V = (GMRESmaxIter+1).ForLoop(i => new double[Nloc]); //   V(1:n,1:m+1) = zeros(n,m);
            MultidimensionalArray H = MultidimensionalArray.Create(GMRESmaxIter+1, GMRESmaxIter+1); //   H(1:m,1:m) = zeros(m,m);
            double[] c = new double[GMRESmaxIter + 1];
            double[] s = new double[GMRESmaxIter + 1];
            double[] y;
            double rho = r.L2Norm();
            errstep = rho;
            double[] g = new double[GMRESmaxIter+1];
            g[0] = rho;

           // Console.WriteLine("Error NewtonGMRES:   " + rho);

            // Termination of entry
            if (rho < GMRESConvCrit)
                return SolutionVec.ToArray();

            V[0].SetV(r, (1.0 / rho));
            double beta = rho;
            int k = 0;

            while ((rho > GMRESConvCrit) && k < GMRESmaxIter) {
                k++;

                // Call directional derivative
                V[k].SetV(dirder(SolutionVec, V[k - 1], f0));
                double normav = V[k].L2Norm();

                // Modified Gram-Schmidt
                for (int j = 1; j <= k; j++) {
                    H[j - 1, k - 1] = GenericBlas.InnerProd(V[k], V[j - 1]).MPISum();
                    V[k].AccV(-H[j - 1, k - 1], V[j - 1]);
                }
                H[k, k - 1] = V[k].L2Norm();
                double normav2 = H[k, k - 1];

                // Reorthogonalize ?
                if (reorth == 1 /*&& normav + 0.001 * normav2 == normav*/) {
                    for (int j = 1; j <= k; j++) {
                        double hr = GenericBlas.InnerProd(V[k], V[j - 1]).MPISum();
                        H[j - 1, k - 1] = H[j - 1, k - 1] + hr;
                        V[k].AccV(-hr, V[j - 1]);
                    }
                    H[k, k - 1] = V[k].MPI_L2Norm();
                }


                // Watch out for happy breakdown
                if (H[k, k - 1] != 0)
                    V[k].ScaleV(1 / H[k, k - 1]);



                // Form and store the information for the new Givens rotation
                if (k > 1) {
                   // for (int i = 1; i <= k; i++) {
                        H.SetColumn(k-1,givapp(c.GetSubVector(0, k - 1), s.GetSubVector(0, k - 1), H.GetColumn(k-1),k-1));
                    //}
                }


                // Don't divide by zero if solution has  been found
                var nu = Math.Sqrt(H[k - 1, k - 1].Pow2() + H[k, k - 1].Pow2());
                if (nu != 0) {
                    c[k - 1] = H[k - 1, k - 1] / nu;
                    s[k - 1] = H[k, k - 1] / nu;
                    H[k - 1, k - 1] = c[k - 1] * H[k - 1, k - 1] - s[k - 1] * H[k, k - 1];
                    H[k, k - 1] = 0;

                    // givapp for g
                    var w1 = c[k - 1] * g[k - 1] - s[k - 1] * g[k];
                    var w2 = s[k - 1] * g[k - 1] + c[k - 1] * g[k];
                    g[k - 1] = w1;
                    g[k] = w2;
                }

                rho = g[k].Abs();

               // Console.WriteLine("Error NewtonGMRES:   " + rho);

            }

            // update approximation and exit

            //y = H(1:i,1:i) \ s(1:i);    
            y = new double[k];
            H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { k-1 , k-1  })
                .Solve(y, g.GetSubVector(0, k));

            int totalIter = k;

            // x = x + V(:,1:i)*y;
            for (int ii = 0; ii <= k; ii++) {
                x.AccV(y[ii], V[ii]);
            }

            errstep = rho;

            return x;

        }

        public double[] Krylov(CoordinateVector SolutionVec, double[] x, double[] f0, out double errstep)
        {           

            double[] step = NewtonGMRES(SolutionVec, f0, new double[x.Length], out errstep);
            int kinn = 0;

            Console.WriteLine("Error Krylov:   " + errstep);

            while (kinn < krylovDim && errstep > ConvCrit) {
                kinn++;

                step = NewtonGMRES(SolutionVec , f0, step, out errstep);

                Console.WriteLine("Error Krylov:   " + errstep);
            }

            return step;

        }

        /// <summary>
        /// Finite difference directional derivative Approximate f'(x) w
        /// C.T.Kelley, April 1, 2003
        /// This code comes with no guarantee or warranty of any kind.
        /// </summary>
        /// <param name="SolutionVec">Solution point</param>
        /// <param name="w">Direction</param>
        /// <param name="f0">f0, usally has been calculated earlier</param>
        /// <returns></returns>
        public double[] dirder(CoordinateVector SolutionVec, double[] w, double[] f0)
        {
            double epsnew = 1E-7;
            int n = SolutionVec.Length;

            // Scale the step
            if (w.L2Norm() == 0) {
                SolutionVec.Clear();
                return SolutionVec.ToArray();
            }

            // Scale the difference increment
            double[] xs = new double[n];
            for (int i = 0; i < n; i++) {
                xs[i] = SolutionVec[i] * w[i];
            }
            xs.ScaleV(w.L2Norm());
            if (xs.L2Norm() == 0) {
                //epsnew = epsnew*Math.Max(xs.,)
            }
            epsnew = epsnew / w.L2Norm();

            double[] temp = new double[SolutionVec.Length];
            temp.CopyEntries(SolutionVec);

            SolutionVec.AccV(epsnew, w);
            double[] fx = new double[f0.Length];

            // Evaluate at SolutionVec
            Init(SolutionVec, RHSRaw, out xs, out fx);

            SolutionVec.CopyEntries(temp);

            // (f1 - f0) / epsnew
            fx.AccV(-1, f0);
            fx.ScaleV(1 / epsnew);

            return fx;

        }

        /// <summary>
        /// Apply a sequence of k Givens rotations, used within gmres codes.  
        /// C.T.Kelley, April 1, 2003
        /// This code comes with no guarantee or warranty of any kind.
        /// </summary>
        /// <param name="c"></param>
        /// <param name="s"></param>
        /// <param name="vin"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public double[] givapp(double[] c, double[] s, double[] vin, int k)
        {
            double[] vrot = vin;
            double w1, w2;

            for (int i = 1; i < k; i++) {
                w1 = c[i - 1] * vrot[i - 1] - s[i - 1] * vrot[i];
                w2 = s[i - 1] * vrot[i - 1] + c[i - 1] * vrot[i];
                vrot[i - 1] = w1;
                vrot[i] = w2;
            }
            return vrot;
        }


    }
}
