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
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using MPI.Wrappers;
using BoSSS.Platform.Utils;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Block-Jacobi smoother, maybe only useful in combination with the multi-grid solver (<see cref="ClassicMultigrid"/>).
    /// </summary>
    public class BlockJacobi : ISolverSmootherTemplate, ISolverWithCallback {

        MultigridOperator m_MultigridOp;


        public void Init(MultigridOperator op) {
            BlockMsrMatrix M = op.OperatorMatrix;
            var MgMap = op.Mapping;
            this.m_MultigridOp = op;
            

            if(!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if(!M.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");

            Mtx = M;
            int L = M.RowPartitioning.LocalLength;

            /*
            diag = new double[L];
            int i0 = Mtx.RowPartitioning.i0;

            for(int i = 0; i < L; i++) {
                diag[i] = Mtx[i0 + i, i0 + i];
            }
             */

            //if (op.Mapping.MaximalLength != op.Mapping.MinimalLength)
            //    // 'BlockDiagonalMatrix' should be completely replaced by 'BlockMsrMatrix'
            //    throw new NotImplementedException("todo - Block Jacobi for variable block Sizes");

            Diag = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
            invDiag = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
            int Jloc = MgMap.LocalNoOfBlocks;
            int j0 = MgMap.FirstBlock;
            MultidimensionalArray temp = null;
            for (int j = 0; j < Jloc; j++) {
                int jBlock = j + j0;
                int Nblk = MgMap.GetBlockLen(jBlock);
                int i0 = MgMap.GetBlockI0(jBlock);

                if (temp == null || temp.NoOfCols != Nblk)
                    temp = MultidimensionalArray.Create(Nblk, Nblk);

                M.ReadBlock(i0, i0, temp);
                Diag.AccBlock(i0, i0, 1.0, temp, 0.0);

                temp.Invert();

                invDiag.AccBlock(i0, i0, 1.0, temp, 0.0);
            }
#if DEBUG
            invDiag.CheckForNanOrInfM();
#endif
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        BlockMsrMatrix Mtx;
        BlockMsrMatrix Diag;
        BlockMsrMatrix invDiag;
        //double[] diag;

        /// <summary>
        /// Jacobi-Damping
        /// </summary>
        public double omega = 1.0; // jacobi - under-relax

        public int NoOfIterations = 4;

        public double m_Tolerance = 0.0;


        /// <summary>
        /// Jacobi iteration
        /// </summary>
        public void Solve<U, V>(U xl, V bl)
            where U : IList<double>
            where V : IList<double> {
            int L = xl.Count;
            double[] ql = new double[L];

            for(int iIter = 0; iIter < NoOfIterations; iIter++) {
                ql.SetV(bl);
                Mtx.SpMV(-1.0, xl, 1.0, ql);
                if(this.m_Tolerance > 0) {
                    double ResNorm = ql.L2NormPow2().MPISum().Sqrt();
                    if(ResNorm < this.m_Tolerance) {
                        m_Converged = true;
                        return;
                    }
                }
                Diag.SpMV(1.0, xl, 1.0, ql);


                //xl.ScaleV(1.0 - omega);
                invDiag.SpMV(omega, ql, 1.0 - omega, xl);

                //for(int i = 0; i < L; i++) {
                //    xl[i] = omega * ((ql[i] + diag[i] * xl[i]) / diag[i]) 
                //        + (1.0 - omega) * xl[i];
                //}
                


                if(this.IterationCallback != null) {
                    double[] _xl = xl.ToArray();
                    double[] _bl = bl.ToArray();
                    Mtx.SpMV(-1.0, _xl, 1.0, _bl);
                    this.IterationCallback(iIter, _xl, _bl, this.m_MultigridOp);
                }
            }
        }

        bool m_Converged = false;
        int m_ThisLevelIterations = 0;

        public int IterationsInNested {
            get { return 0; }
        }

        public int ThisLevelIterations {
            get { return this.m_ThisLevelIterations; }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public void ResetStat() {
            m_Converged = false;
            m_ThisLevelIterations = 0;

        }
    }
}
