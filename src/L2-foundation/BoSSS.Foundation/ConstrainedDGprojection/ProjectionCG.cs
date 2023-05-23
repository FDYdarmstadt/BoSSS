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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.LinSolvers;
using MPI.Wrappers;

namespace BoSSS.Foundation.ConstrainedDGprojection {
    public static class SolverUtils {
        /// <summary>
        /// convert global BlockMsrMatrix to local Msr matrix
        /// </summary>
        /// <param name="Matrix"></param>
        /// <returns></returns>
        public static MsrMatrix GetLocalMatrix(BlockMsrMatrix Matrix) {
            int rank;
            MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out rank);

            // extraction of local Matrix block, ILU will be executed process local
            var LocBlocki0s = new List<long>();
            var LocBlockLen = new List<int>();
            long IdxOffset = Matrix._RowPartitioning.i0;
            for (int i = 0; i < Matrix._RowPartitioning.LocalNoOfBlocks; i++) {
                long iBlock = i + Matrix._RowPartitioning.FirstBlock;
                long i0 = Matrix._RowPartitioning.GetBlockI0(iBlock) - IdxOffset;
                int Len = Matrix._RowPartitioning.GetBlockLen(iBlock);
                LocBlocki0s.Add(i0);
                LocBlockLen.Add(Len);
            }
            long[] RowISrc = Matrix._RowPartitioning.LocalLength.ForLoop(i => i + IdxOffset);

            var part = new BlockPartitioning(Matrix._RowPartitioning.LocalLength, LocBlocki0s, LocBlockLen, csMPI.Raw._COMM.SELF);
            BlockMsrMatrix localMatrix = new BlockMsrMatrix(part);

            Matrix.WriteSubMatrixTo(localMatrix, RowISrc, default(long[]), RowISrc, default(long[]));
            return localMatrix.ToMsrMatrix();
        }

        public static void CheckAndRepairMatrix(BlockMsrMatrix M) {
            int rows = M.RowPartitioning.LocalLength;
            long i0 = M.RowPartitioning.i0;
            long iE = M.RowPartitioning.iE;
            for(long iRow = i0; iRow < iE; iRow++) {
                if (M.GetDiagonalElement(iRow) == 0)
                    M.SetDiagonalElement(iRow,1.0);
            }
        }

        /// <summary>
        /// constructs solver depending on size of matrix (DOF). Acts on MPI_Comm.WORLD
        /// </summary>
        /// <param name="DOF"></param>
        /// <returns></returns>
        public static ISparseSolver GlobalSolverFactory(long DOF) {
            bool UseDirect = DOF < 1E6;
            if (UseDirect)
                //return new ilPSP.LinSolvers.MUMPS.MUMPSSolver() {
                //    Parallelism = Parallelism.MPI
                //};
                return new ilPSP.LinSolvers.PARDISO.PARDISOSolver() {
                    Parallelism = Parallelism.OMP,
                    SymmIndefPivot = true,
                    CacheFactorization = false
                };
            else
                return new myCG();
        }

        /// <summary>
        /// constructs solver which acts on MPI_Comm.SELF
        /// </summary>
        /// <returns></returns>
        public static ISparseSolver PatchSolverFactory() {
            return new ilPSP.LinSolvers.PARDISO.PARDISOSolver() {
                Parallelism = Parallelism.SEQ,
                SymmIndefPivot = true,
                CacheFactorization = false
            };
            //return new ilPSP.LinSolvers.MUMPS.MUMPSSolver() {
            //    Parallelism = Parallelism.SEQ,
            //};
        }
    }

    public class myCG : ISparseSolver {

        public void DefineMatrix(IMutableMatrixEx matrix) {
            BlockMsrMatrix bmsrM;
            try {
                bmsrM = (BlockMsrMatrix)matrix;
            } catch (InvalidCastException) {
                throw new NotSupportedException("this solver only supports block msr matrix as input!");
            }
            m_matrix = bmsrM;

            //SolverUtils.CheckAndRepairMatrix(m_matrix);

            m_ILU = new ilPSP.LinSolvers.HYPRE.Euclid() {
                Level = 0,
            };
            //if (m_localMsrMatrix != null && m_localMsrMatrix.RowPartitioning.LocalLength > 0) ((ilPSP.LinSolvers.HYPRE.Solver)m_ILU).DefineMatrix(m_matrix);

            m_PCG = new ilPSP.LinSolvers.HYPRE.PCG() {
                NestedPrecond = m_ILU,
                PrintLevel = DiagnosticLevel,
                ConvergenceType = ConvergenceTypes.Relative,
                MaxIterations = 100,
                Tolerance = 1E-6,
            };
            m_PCG.DefineMatrix(m_matrix);
        }

        public int DiagnosticLevel = 0;
        BlockMsrMatrix m_matrix;
        ISparseSolver m_PCG;
        IImplicitPrecond m_ILU;
        //MsrMatrix m_localMsrMatrix;

        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            var StartTime = DateTime.Now;
            m_PCG.Solve(x, rhs);
            return new SolverResult() {
                Converged = true,
                NoOfIterations = 1,
                RunTime = DateTime.Now - StartTime
            };
        }
        public void Dispose() {
            m_PCG.Dispose();
            m_ILU.Dispose();
            m_matrix = null;
        }
    }
}
