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
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {
    public class Jacobi : ISubsystemSolver {

        IOperatorMappingPair m_mgop;


        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        public void Init(MultigridOperator op) {
            InitImpl(op);
        }
        void InitImpl(IOperatorMappingPair op) {
            using(new FuncTrace()) {
                if(object.ReferenceEquals(op, m_mgop))
                    return; // already initialized
                else
                    this.Dispose();

                var M = op.OperatorMatrix;
                var MgMap = op.DgMapping;
                this.m_mgop = op;

                if(!M.RowPartitioning.EqualsPartition(MgMap))
                    throw new ArgumentException("Row partitioning mismatch.");
                if(!M.ColPartition.EqualsPartition(MgMap))
                    throw new ArgumentException("Column partitioning mismatch.");

                Mtx = M;
                int L = M.RowPartitioning.LocalLength;

                diag = new double[L];
                long i0 = Mtx.RowPartitioning.i0;

                for(int i = 0; i < L; i++) {
                    diag[i] = Mtx[i0 + i, i0 + i];
                }
            }
        }

        //public Action<int, double[], double[], MultigridOperator> IterationCallback {
        //    get;
        //    set;
        //}

        BlockMsrMatrix Mtx;
        double[] diag;

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
            where V : IList<double> 
        {
            int L = xl.Count;
            double[] ql = new double[L];

            for(int iIter = 0; iIter < NoOfIterations; iIter++) {
                this.m_ThisLevelIterations++;
                ql.SetV(bl);
                Mtx.SpMV(-1.0, xl, 1.0, ql);
                if(this.m_Tolerance > 0) {
                    double ResNorm = ql.L2NormPow2().MPISum().Sqrt();
                    if(ResNorm < this.m_Tolerance) {
                        m_Converged = true;
                        return;
                    }
                }

                for(int i = 0; i < L; i++) {
                    xl[i] = omega * ((ql[i] + diag[i] * xl[i]) / diag[i]) + (1.0 - omega) * xl[i];
                }

                //if(this.IterationCallback != null) {
                //    double[] _xl = xl.ToArray();
                //    double[] _bl = bl.ToArray();
                //    Mtx.SpMV(-1.0, _xl, 1.0, _bl);
                //    this.IterationCallback(iIter, _xl, _bl, this.m_mgop);
                //}
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

        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

        public void Dispose() {
            throw new NotImplementedException();
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }
}
