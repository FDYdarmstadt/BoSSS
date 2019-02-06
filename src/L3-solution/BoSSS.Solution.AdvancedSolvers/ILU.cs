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
    /// Algorithm based on "Numerik linearer Gleichungssysteme, 2. Auflage, Vieweg, Andreas Meister"
    /// </summary>
    public class ILU : ISolverSmootherTemplate, ISolverWithCallback {

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

        bool m_Converged = false;
        int m_ThisLevelIterations = 0;

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get {
                throw new NotImplementedException();
            }
            set {
                throw new NotImplementedException();
            }
        }

        MultigridOperator m_mgop;

        public void Init(MultigridOperator op) {
            var M = op.OperatorMatrix;
            var MgMap = op.Mapping;
            this.m_mgop = op;

            if (!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!M.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");

            Mtx = M;
            int L = M.RowPartitioning.LocalLength;

            diag = new double[L];
            int i0 = Mtx.RowPartitioning.i0;
        }

        BlockMsrMatrix Mtx;
        double[] diag;

        public void Solve<P, Q>(P X, Q B)
            where P : IList<double>
            where Q : IList<double> {

            int n = Mtx.NoOfCols;
            double sum = 0;

            var tempMtx = Mtx;

            // Zeros on diagonal elements because of saddle point structure
            for (int bla = 0; bla < n; bla++) {
                if (Mtx.GetDiagonalElement(bla) == 0)
                    throw new Exception("One or more diagonal elements are zero, ILU cannot work");
                Mtx.SetDiagonalElement(bla, 1);
            }

            // ILU decomposition of matrix
            for (int k = 0; k < n - 1; k++) {
                for (int i = k; i < n; i++) {
                    if (tempMtx[i, k] == 0) { i = n; }
                    else {
                        Mtx[i, k] = Mtx[i, k] / Mtx[k, k];
                        for (int j = k + 1; j < n; j++) {
                            if (tempMtx[i, j] == 0) {
                                j = n;
                            }
                            else {
                                Mtx[i, j] = Mtx[i, j] - Mtx[i, k] * Mtx[k, j];
                            }
                        }
                    }
                }
            }


            // LU decomposition of matrix        
            //for (int i = 0; i < n; i++) {
            //    for (int j = i; j < n; j++) {
            //        sum = 0;
            //        for (int k = 0; k < i; k++)
            //            sum += Mtx[i, k] * Mtx[k, j];
            //        Mtx[i, j] = tempMtx[i, j] - sum;
            //    }
            //    for (int j = i + 1; j < n; j++) {
            //        sum = 0;
            //        for (int k = 0; k < i; k++)
            //            sum += Mtx[j, k] * Mtx[k, i];
            //        Mtx[j, i] = (1 / Mtx[i, i]) * (tempMtx[j, i] - sum);
            //    }
            //}

            // find solution of Ly = b
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                sum = 0;
                for (int k = 0; k < i; k++)
                    sum += Mtx[i, k] * y[k];
                y[i] = B[i] - sum;
            }
            // find solution of Ux = y
            for (int i = n - 1; i >= 0; i--) {
                sum = 0;
                for (int k = i + 1; k < n; k++)
                    sum += Mtx[i, k] * X[k];
                X[i] = (1 / Mtx[i, i]) * (y[i] - sum);
            }

        }


    }
}

