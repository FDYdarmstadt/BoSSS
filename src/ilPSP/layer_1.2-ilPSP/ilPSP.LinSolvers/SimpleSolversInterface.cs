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

using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ilPSP.LinSolvers {

    /// <summary>
    /// Simplified interfaces to solvers.
    /// </summary>
    static public class SimpleSolversInterface {

        /// <summary>
        /// Solves a symmetric, definite system using the conjugate gradient algorithm.
        /// </summary>
        /// <param name="Matrix">the matrix of the linear problem.</param>
        /// <param name="X">Output, (hopefully) the solution.</param>
        /// <param name="B">Input, the right-hand-side.</param>
        /// <param name="MaxIterations"></param>
        /// <param name="Tolerance"></param>
        /// <returns>Actual number of iterations.</returns>
        static public int Solve_CG<V, W>(this IMutableMatrixEx Matrix, V X, W B, int MaxIterations = 100000, double Tolerance = 1.0e-10)
            where V : IList<double>
            where W : IList<double> //
        {
            using (var slv = new ilPSP.LinSolvers.monkey.CG()) {
                slv.MaxIterations = MaxIterations;
                slv.Tolerance = Tolerance;
                slv.DevType = monkey.DeviceType.CPU;

                slv.DefineMatrix(Matrix);

                var SolRes = slv.Solve(X, B.ToArray());
                return SolRes.NoOfIterations;
            }
        }

        static string CheckResidual<V, W>(this IMutableMatrixEx Matrix, V X, W B, Type t)
            where V : IList<double>
            where W : IList<double> //
        {
            double[] Residual = B.ToArray();
            if (object.ReferenceEquals(B, Residual))
                throw new ApplicationException("ToArray does not work as expected.");
            double RhsNorm = B.L2NormPow2().MPISum(Matrix.MPI_Comm).Sqrt();
            double MatrixInfNorm = Matrix.InfNorm();
            Matrix.SpMV(-1.0, X, 1.0, Residual);

            double ResidualNorm = Residual.L2NormPow2().MPISum(Matrix.MPI_Comm).Sqrt();
            double SolutionNorm = X.L2NormPow2().MPISum(Matrix.MPI_Comm).Sqrt();
            double Denom = Math.Max(MatrixInfNorm, Math.Max(RhsNorm, Math.Max(SolutionNorm, Math.Sqrt(BLAS.MachineEps))));
            double RelResidualNorm = ResidualNorm / Denom;


            if (RelResidualNorm > 1.0e-10) {
                string ErrMsg;
                using (var stw = new System.IO.StringWriter()) {
                    stw.WriteLine("Linear Solver: High residual from direct solver " + t + ".");
                    stw.WriteLine("    L2 Norm of RHS:         " + RhsNorm);
                    stw.WriteLine("    L2 Norm of Solution:    " + SolutionNorm);
                    stw.WriteLine("    L2 Norm of Residual:    " + ResidualNorm);
                    stw.WriteLine("    Relative Residual norm: " + RelResidualNorm);
                    stw.WriteLine("    Matrix Inf norm:        " + MatrixInfNorm);

                    ErrMsg = stw.ToString();
                }
                return ErrMsg;
            } else {
                return null;
            }
        }

        /// <summary>
        /// Solves a linear system using a direct solver.
        /// </summary>
        /// <param name="Matrix">the matrix of the linear problem.</param>
        /// <param name="B">Input, the right-hand-side.</param>
        /// <remarks>
        /// The solver checks the residual norm and tries a different solver library if the norm is deemed to high.
        /// (E.g. on some processor architectures, <see cref="PARDISO.PARDISOSolver"/> is known to fail from time to time,
        /// especially in OpenMP-Mode).
        /// </remarks>
        static public double[] Solve_Direct<W>(this IMutableMatrixEx Matrix, W B)
            where W : IList<double> //
        {
            var X = new double[B.Count];
            Solve_Direct(Matrix, X, B);
            return X;
        }

        /// <summary>
        /// Solves a linear system using a direct solver.
        /// </summary>
        /// <param name="Matrix">the matrix of the linear problem.</param>
        /// <param name="X">Output, (hopefully) the solution.</param>
        /// <param name="B">Input, the right-hand-side.</param>
        /// <remarks>
        /// The solver checks the residual norm and tries a different solver library if the norm is deemed to high.
        /// (E.g. on some processor architectures, <see cref="PARDISO.PARDISOSolver"/> is known to fail from time to time,
        /// especially in OpenMP-Mode).
        /// </remarks>
        static public void Solve_Direct<V, W>(this IMutableMatrixEx Matrix, V X, W B)
        where V : IList<double>
        where W : IList<double> //
    {
            var SolverFallbackSeq = new Func<ISparseSolver>[] {
                () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver() { Parallelism = Parallelism.OMP },
                () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver() { Parallelism = Parallelism.SEQ },
                () => new ilPSP.LinSolvers.MUMPS.MUMPSSolver() { Parallelism = Parallelism.SEQ }
            };

            if (Matrix.MPI_Comm.Equals(csMPI.Raw._COMM.SELF)) // solve matrices on SELF always sequentially!
                SolverFallbackSeq = SolverFallbackSeq.Skip(1).ToArray();


            //*

            string LastError = null;
            foreach (var f in SolverFallbackSeq) {
                using (var slv = f()) {
                    //Console.WriteLine("Using solver: " + slv.GetType() + " ...");
                    slv.DefineMatrix(Matrix);
                    var SolRes = slv.Solve(X, B.ToArray());

                    LastError = CheckResidual(Matrix, X, B, slv.GetType());
                    if (LastError == null) {
                        //Console.WriteLine("residual is fine.");
                        return;
                    } else {
                        //Console.WriteLine("some error.");
                    }
                }
            }

            if (LastError != null)
                throw new ArithmeticException(LastError);

        }

        /// <summary>
        /// Condition number estimate by MUMPS; Rem.: MUMPS manual does not tell in which norm; it does not seem to be as reliable as MATLAB condest
        /// </summary>
        public static double Condest_MUMPS(this IMutableMatrixEx Mtx) {
            using (new FuncTrace()) {
                using (var slv = new ilPSP.LinSolvers.MUMPS.MUMPSSolver() { Parallelism = Parallelism.MPI }) {
                    slv.Statistics = MUMPS.MUMPSStatistics.AllStatistics;

                    slv.DefineMatrix(Mtx);
                    double[] dummyRHS = new double[Mtx.RowPartitioning.LocalLength];
                    double[] dummySol = new double[dummyRHS.Length];
                    dummyRHS.SetAll(1.11);

                    slv.Solve(dummySol, dummyRHS);

                    return slv.LastCondNo;
                }
            }
        }

        /// <summary>
        /// Computes the minimal Eigenvalue and related Eigenvector 
        /// using an Arnoldi iteration 
        /// on the inverse matrix.
        /// Note that the inverse is not computed; instead the required products of inverse matrix times a vector are 
        /// computed using the PARDISO solver.
        /// </summary>
        static public (double lambdaMin, double[] V) MinimalEigen(this IMutableMatrixEx Mtx, double tol = 1.0e-6) {

            int L = Mtx.RowPartitioning.LocalLength;
            using (var slv = new ilPSP.LinSolvers.PARDISO.PARDISOSolver()) {
                slv.CacheFactorization = true;
                slv.DefineMatrix(Mtx);

                double invMinLambda = 0;
                double[] Evect = new double[L];
                double[] tmp = new double[L];

                Evect.FillRandom();
                for (int i = 0; i < 100; i++) {
                    double norm = Evect.MPI_L2Norm(Mtx.MPI_Comm);
                    Evect.ScaleV(1.0 / norm);
                    slv.Solve(tmp, Evect);


                    string Correcto = Mtx.CheckResidual(tmp, Evect, slv.GetType());
                    if (Correcto != null) {
                        // PARDISO seems unable to solve the system:
                        // use the fallback levels in other routine:

                        Mtx.Solve_Direct(tmp, Evect);
                    }

                    // Monitor change in Eigenvalue:
                    double prev_invMinLambda = invMinLambda;
                    invMinLambda = tmp.MPI_InnerProd(Evect, Mtx.MPI_Comm);
                    double ChangeNorm = Math.Abs(invMinLambda - prev_invMinLambda) / Math.Max(invMinLambda.Abs(), prev_invMinLambda.Abs());
                    //Console.WriteLine(i + " -- Change norm is: " + ChangeNorm);

                    // prepare for next loop:
                    var a = Evect;
                    Evect = tmp;
                    tmp = a;

                    if (ChangeNorm < tol)
                        break;
                }

                double normFin = Evect.MPI_L2Norm(Mtx.MPI_Comm);
                Evect.ScaleV(1.0 / normFin);
                return (1.0 / invMinLambda, Evect);
            }
        }

        /// <summary>
        /// Computes the maximal Eigenvalue and Eigenvector using an Arnoldi iteration
        /// </summary>
        /// <param name="Mtx"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        static public (double lambdaMax, double[] V) MaximalEigen(this IMutableMatrixEx Mtx, double tol = 1.0e-6) {
            int L = Mtx.RowPartitioning.LocalLength;
            double[] Evect = new double[L];
            double[] tmp = new double[L];
            double MaxLambda = 0;
            Evect.FillRandom();

            var comm = Mtx.MPI_Comm;

            for (int i = 0; i < 100; i++) {
                double norm = Evect.MPI_L2Norm(comm);
                Evect.ScaleV(1.0 / norm);
                Mtx.SpMV(1.0, Evect, 0.0, tmp);
                double prev_MaxLambda = MaxLambda;
                MaxLambda = tmp.MPI_InnerProd(Evect, comm);
                double ChangeNorm = Math.Abs(MaxLambda - prev_MaxLambda) / Math.Max(MaxLambda.Abs(), prev_MaxLambda.Abs());
                var a = Evect;
                Evect = tmp;
                tmp = a;

                if (ChangeNorm < tol)
                    break;
            }
            double normFin = Evect.MPI_L2Norm(comm);
            Evect.ScaleV(1.0 / normFin);
            return (MaxLambda, Evect);
        }

        /// <summary>
        /// condition number estimate employing the Arnoldi iterations <see cref="MinimalEigen"/> and <see cref="MaximalEigen"/>
        /// </summary>
        static public double condestArnoldi(this IMutableMatrixEx Mtx, double tol = 1.0e-6) {
            (double lambdaMax, _) = MaximalEigen(Mtx, tol);
            (double lambdaMin, _) = MinimalEigen(Mtx, tol);
            return Math.Abs(lambdaMax / lambdaMin);
        }
    }
}
