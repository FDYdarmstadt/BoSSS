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
using MPI.Wrappers;
using BoSSS.Platform;
using BoSSS.Platform.Utils;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Parallel ILU from HYPRE library
    /// </summary>
    public class HypreAMR : ISubsystemSolver, ISolverWithCallback, IDisposable {

        public int IterationsInNested {
            get { return 0; }
        }

        public int ThisLevelIterations {
            get { return this.m_ThisLevelIterations; }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get { return m_IterationCallback; }
            set { m_IterationCallback = value; }
        }

        public void ResetStat() {
            m_Converged = false;
            m_ThisLevelIterations = 0;
        }

        Action<int, double[], double[], MultigridOperator> m_IterationCallback;
        bool m_Converged = false;
        int m_ThisLevelIterations = 0;

        private int FillInLevel {
            get {
                int D = m_mgop.DgMapping.SpatialDimension;
                switch(D) {
                    case 1:
                    case 2:
                    return 2;
                    case 3:
                    return 1;
                    default:
                    throw new NotSupportedException("spatial Dimension (" + D + ") not supported.");
                }
            }
        }

        void InitImpl(IOperatorMappingPair op) {
            var Mtx = op.OperatorMatrix;
            var MgMap = op.DgMapping;
            this.m_mgop = op;

            if (!Mtx.RowPartitioning.EqualsPartition(MgMap))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!Mtx.ColPartition.EqualsPartition(MgMap))
                throw new ArgumentException("Column partitioning mismatch.");

            m_matrix = Mtx.CloneAs();

            AMG = new ilPSP.LinSolvers.HYPRE.BoomerAMG {
                RelaxType = ilPSP.LinSolvers.HYPRE.RelaxType.GaussSeidel,
                PrintLevel = 2,
                Tolerance = 1E-8,
                MaxIterations = 1000,
            };

            // Zeros on diagonal elements because of saddle point structure
            //long i0 = m_ILUmatrix._RowPartitioning.i0;
            //long iE = m_ILUmatrix._RowPartitioning.iE;
            //for (long iRow = i0; iRow < iE; iRow++) {
            //    if (m_ILUmatrix.GetDiagonalElement(iRow) == 0) {
            //        m_ILUmatrix.SetDiagonalElement(iRow, 1);
            //    }
            //}
            if (m_matrix != null && m_matrix.RowPartitioning.LocalLength > 0) 
                AMG.DefineMatrix(m_matrix);
        }

        ilPSP.LinSolvers.HYPRE.Solver AMG;
        BlockMsrMatrix m_matrix;
        IOperatorMappingPair m_mgop;

        public void Solve<P, Q>(P X, Q B)
            where P : IList<double>
            where Q : IList<double> //
        {
            var SolverResult = AMG.Solve(X, B);
            m_ThisLevelIterations += SolverResult.NoOfIterations;
        }


        public void Dispose() {
            AMG.Dispose();
            m_matrix = null;
            m_mgop = null;
        }

        public object Clone() {
            var clone = new HypreAMR();
            if(this.IterationCallback !=null) clone.IterationCallback = this.IterationCallback.CloneAs();
            return clone;
        }

        public long UsedMemory() {
            return 0;
        }

        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        public void Init(MultigridOperator op) {
            InitImpl(op);
        }
    }
}

