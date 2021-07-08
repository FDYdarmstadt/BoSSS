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
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;

namespace ilPSP.Connectors.Matlab {
    
    /// <summary>
    /// shortcuts to some MATLAB routines
    /// </summary>
    static public class Extensions {

        /// <summary>
        /// Evaluation of the condition number of a full matrix
        /// <paramref name="M"/>.
        /// </summary>
        /// <param name="M">A full square matrix</param>
        /// <param name="workingPath"></param>
        /// <returns>
        /// The condition number of <paramref name="M"/>
        /// </returns>
        public static double cond(this IMatrix M, string workingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector( WorkingPath: workingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutMatrix(M, "Matrix");
                connector.Cmd("cond = cond(Matrix)");
                connector.GetMatrix(output, "cond");

                connector.Execute(false);

                return output[0, 0];
            }
        }

        /*
        /// <summary>
        /// Evaluation of the eigenvalues of a full matrix
        /// <paramref name="M"/>.
        /// </summary>
        /// <param name="M">A full square matrix</param>
        /// <param name="workingPath"></param>
        /// <returns>
        /// The eigenvalues of <paramref name="M"/> in ascending order
        /// </returns>
        public static double[] eig(this IMatrix M, string workingPath = null) {
            if (M.NoOfCols != M.NoOfRows) {
                throw new ArgumentException("Matrix must be square");
            }
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: workingPath)) {
                MultidimensionalArray output = MultidimensionalArray.Create(M.NoOfCols, 1);
                connector.PutMatrix(M, "Matrix");
                connector.Cmd("eig = eig(Matrix)");
                connector.GetMatrix(output, "eig");

                connector.Execute(false);

                return output.Storage;
            }
        }
        */

        /// <summary>
        /// MATLAB function 'condest' (condition number estimation for sparse
        /// matrices)
        /// </summary>
        public static double condest(this IMutableMatrixEx M, string __WorkingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd("cond = condest(Matrix)");
                connector.GetMatrix(output, "cond");

                connector.Execute(false);

                return output[0, 0];
            }
        }

        /// <summary>
        /// Assumes that there are Matrices on Each process with Self_Communicator (alias local matrices).
        /// </summary>
        /// <param name="M"></param>
        /// <param name="__WorkingPath"></param>
        /// <param name="printout">switch on printing of results</param>
        /// <returns>condition number for each local matrix</returns>
        public static double[] Multicondest(this IMutableMatrixEx M, string __WorkingPath = null, bool printout = false) {
            int MPISize = -1;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Size(comm, out MPISize);

            if (M.RowPartitioning.MpiSize != 1)
                throw new NotSupportedException("input is no local matrix!");
            if (M == null)
                throw new ArgumentNullException();
            double[] retcondest = new double[MPISize];

            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(MPISize, 1);
                connector.PutSparseMatrixRankExclusive(M, "Matrix");
                connector.CreateOutputBuffer(MPISize,"cond");
                for (int iRank=0; iRank < MPISize; iRank++) {
                    string cmdline = String.Format("cond({0}) = condest(Matrix_{1});",iRank+1,iRank);
                    connector.Cmd(cmdline);
                }
                connector.GetMatrix(output, "cond");
                connector.Execute(false);
                retcondest=MPISize.ForLoop(iRank => output[iRank, 0]);
            }
            if (printout) {
                for(int iRank=0; iRank<MPISize;iRank++)
                    Console.WriteLine("condest on rank {0}: {1}", iRank, retcondest[iRank]);
            }
            return retcondest;
        }

        private static void CreateOutputBuffer(this BatchmodeConnector bmc, int Length, string name) {
            string cmdline = name+"=[0";
            for(int i=0; i < Length-1; i++) {
                cmdline += ";0";
            }
            cmdline += "];";
            bmc.Cmd(cmdline);
        }

        /// <summary>
        /// Tests, via a Cholesky factorization, if a symmetric matrix is positive definite.
        /// </summary>
        public static bool IsPosDef(this IMutableMatrixEx M, string __WorkingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd("[V,r]=chol(0.5*(Matrix+Matrix'));");
                connector.GetMatrix(output, "r");

                connector.Execute(false);

                return output[0, 0] == 0;
            }
        }

        /// <summary>
        /// Tests, via a Cholesky factorization, if a symmetric matrix is negative definite.
        /// </summary>
        public static bool IsNegDef(this IMutableMatrixEx M, string __WorkingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd("[V,r]=chol(-0.5*(Matrix+Matrix'));");
                connector.GetMatrix(output, "r");

                connector.Execute(false);

                return output[0, 0] == 0;
            }
        }

        /// <summary>
        /// Tests, via a Cholesky factorization, if a symmetric matrix is positive or negative definite.
        /// </summary>
        public static bool IsDefinite(this IMutableMatrixEx M, string __WorkingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 2);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd("[V,pr]=chol( 0.5*(Matrix+Matrix'));");
                connector.Cmd("[V,nr]=chol(-0.5*(Matrix+Matrix'));");
                connector.Cmd("ret=[pr,nr]");
                connector.GetMatrix(output, "ret");

                connector.Execute(false);

                return output[0, 0] == 0 || output[0, 1] == 0;
            }
        }

        /// <summary>
        /// MATLAB function 'cond' (condition number for sparse matrices)
        /// </summary>
        public static double cond(this IMutableMatrixEx M, string __WorkingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd("cond = cond(Matrix)");
                connector.GetMatrix(output, "cond");

                connector.Execute(false);

                return output[0, 0];
            }
        }



        /// <summary>
        /// MATLAB function 'det' (determinante of a sparse matrix);
        /// </summary>
        public static double det(this IMutableMatrixEx M, string __WorkingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd("detM = det(Matrix)");
                connector.GetMatrix(output, "detM");

                connector.Execute(false);

                return output[0, 0];
            }

        }

        /// <summary>
        /// MATLAB function 'rank' (rank of a matrix);
        /// </summary>
        public static double rank(this IMutableMatrixEx M, string __WorkingPath = null) {
            if (M == null)
                throw new ArgumentNullException();
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd("rank = rank(full(Matrix))");
                connector.GetMatrix(output, "rank");

                connector.Execute(false);

                return output[0, 0];
            }

        }


        /// <summary>
        /// MATLAB 'backslash' solver for the system
        /// <paramref name="M"/>*<paramref name="X"/> = <paramref name="RHS"/>.
        /// </summary>
        /// <param name="M">
        /// Matrix of the linear system.
        /// </param>
        /// <param name="RHS">
        /// Input, the right hand side of the linear system.
        /// </param>
        /// <param name="X">
        /// Output, the solution.
        /// </param>
        /// <param name="__WorkingPath"></param>
        public static void SolveMATLAB<T1,T2>(this IMutableMatrixEx M, T1 X, T2 RHS, string __WorkingPath = null)
            where T1 : IList<double> 
            where T2 : IList<double> //
        {
            if (M == null)
                throw new ArgumentNullException();
            if (X == null)
                throw new ArgumentNullException();
            if (RHS == null)
                throw new ArgumentNullException();
            if (M.RowPartitioning.LocalLength != RHS.Count)
                throw new ArgumentException("Mismatch between number of rows and length of right-hand-side.");
            if (M.ColPartition.LocalLength != X.Count)
                throw new ArgumentException("Mismatch between number of columns and length of unknown vector.");

            MultidimensionalArray Xwrapper = MultidimensionalArray.Create(X.Count, 1);

            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {

                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.PutVector(RHS, "RHS");

                connector.Cmd("X = Matrix \\ RHS ;");
                connector.GetMatrix(Xwrapper, "X");

                connector.Execute(false);

                Xwrapper.GetColumn(0, X);
            }
        }


        /// <summary>
        /// MATLAB function 'eigs' (eigenvalues of a sparse matrix);
        /// </summary>
        /// <param name="C">
        /// options, ('lm': largest magnitude, etc.) see MATLAB documentation.
        /// </param>
        /// <param name="K">
        /// Number of eigenvalues.
        /// </param>
        /// <param name="M">
        /// A sparse Matrix.
        /// </param>
        /// <param name="__WorkingPath"></param>
        public static double[] eigs(this IMutableMatrixEx M, int K, string C, string __WorkingPath = null) {
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {
                if (M == null)
                    throw new ArgumentNullException();
                if(K < 1)
                    throw new ArgumentOutOfRangeException();

                MultidimensionalArray output = MultidimensionalArray.Create(K, 1);
                connector.PutSparseMatrix(M, "Matrix");
                connector.Cmd(string.Format("EV = eigs(Matrix,{0},'{1}');", K, C));
                connector.GetMatrix(output, "EV");

                connector.Execute(false);

                return output.GetColumn(0);
            }
        }


        /// <summary>
        /// MATLAB function 'eigs' (eigenvalues of a dense matrix);
        /// </summary>
        /// <param name="M">
        /// A sparse Matrix.
        /// </param>
        /// <param name="__WorkingPath"></param>
        public static double[] eigs(this IMatrix M, string __WorkingPath = null) {
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {
                if (M == null)
                    throw new ArgumentNullException();
                if(M.NoOfCols != M.NoOfRows) {
                    throw new ArgumentException("Not supported for non-symmetrical matrices.");
                }
                int N = M.NoOfCols;

                MultidimensionalArray output = MultidimensionalArray.Create(N, 1);
                connector.PutMatrix(M, "Matrix");
                connector.Cmd(string.Format("EV = eigs(Matrix,{0});", M.NoOfCols));
                connector.GetMatrix(output, "EV");

                connector.Execute(false);

                return output.GetColumn(0);
            }
        }
        
        /// <summary>
        /// MATLAB function 'eigs' (eigenvalues and eigenvectors);
        /// </summary>
        /// <param name="M">
        /// A sparse Matrix.
        /// </param>
        /// <param name="__WorkingPath"></param>
        public static (double[] EigenVals, MultidimensionalArray EigenVect) eigsV(this IMatrix M, string __WorkingPath = null) {
            using (var connector = new BatchmodeConnector(WorkingPath: __WorkingPath)) {
                if (M == null)
                    throw new ArgumentNullException();
                if(M.NoOfCols != M.NoOfRows) {
                    throw new ArgumentException("Not supported for non-symmetrical matrices.");
                }
                int N = M.NoOfCols;

                MultidimensionalArray eVect = MultidimensionalArray.Create(N, N);
                MultidimensionalArray eVal = MultidimensionalArray.Create(N, 1);
                connector.PutMatrix(M, "Matrix");
                connector.Cmd(string.Format("[V,D] = eigs(Matrix,{0});", M.NoOfCols));
                connector.Cmd(string.Format("ev = diag(D);"));
                connector.GetMatrix(eVect, "V");
                connector.GetMatrix(eVal, "ev");

                connector.Execute(false);

                return (eVal.GetColumn(0), eVect);
            }
        }
    }
}
