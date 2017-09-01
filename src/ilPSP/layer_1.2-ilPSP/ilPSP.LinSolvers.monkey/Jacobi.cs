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
using System.Text;

namespace ilPSP.LinSolvers.monkey {

    /// <summary>
    /// A Jacobi solver - maybe useful as preconditioner, or as smoother
    /// </summary>
    public class Jacobi : Solver {

        /// <summary>
        ///  "one over diagonal elements" of the matrix
        /// </summary>
        VectorBase CreateInvDiag() {
            VectorBase InvDiag = Device.CreateVector(m_Matrix.RowPartitioning);
            int i0 = (int)m_Matrix.RowPartitioning.i0;
            int L = (int)m_Matrix.RowPartitioning.LocalLength;
            for (int i = 0; i < L; i++) {
                InvDiag[i] = 1.0 / m_Matrix.GetDiagonalElement(i + i0);
            }
            return InvDiag;
        }

        double m_UnderRelaxationFactor = 1.0;

        /// <summary>
        /// gets/sets the realxation factor for the Jacobi iteration; Default is 1.0;<br/>
        /// Recommended range is between 0 (excluding) and 1.0 (including);
        /// </summary>
        public double UnderRelaxationFactor {
            get { return m_UnderRelaxationFactor; }
            set {
                if (m_UnderRelaxationFactor <= 0)
                    throw new ArgumentOutOfRangeException("must be positive");
                m_UnderRelaxationFactor = value;
            }
        }

        /// <summary>
        /// executes the Jacobi iteration
        /// </summary>
        protected override void CallSolver(VectorBase x, VectorBase rhs, ref SolverResult stats) {

            // create objects 
            // ==============

            VectorBase.CommVector _xComm = x.CreateCommVector(m_Matrix);
            VectorBase tmp = Device.CreateVector(x.Part);
            VectorBase InvDiag = CreateInvDiag();

            // lock objects
            // ============
            x.Lock();
            m_Matrix.Lock();
            rhs.Lock();
            tmp.Lock();
            InvDiag.Lock();


            // iterate
            // =======
            stats.Converged = false;
            stats.NoOfIterations = 0;
            double residualNorm = double.MaxValue;
            double r_0 = double.NaN;
            while (true) {

                // loop termination
                // ================
                if (stats.NoOfIterations >= m_MinIterations) { // do at least the minimum number of iterations

                    if (residualNorm <= m_Tolerance) {
                        // success
                        stats.Converged = true;
                        break;
                    }

                    if (stats.NoOfIterations >= m_MaxIterations)
                        // terminate
                        break;
                }

                // Jacobi iteration
                // ================
                m_Matrix.SpMV_Expert(-1.0, _xComm, 0.0, tmp); // tmp = -M*x
                tmp.Acc(1.0, rhs);                            // tmp = -M*x + rhs
                if (m_UnderRelaxationFactor != 1.0)
                    tmp.Scale(m_UnderRelaxationFactor);

                double r = Math.Sqrt(tmp.TwoNormSquare());
                if (stats.NoOfIterations == 0)
                    r_0 = r;
                if (m_ConvergenceType == ConvergenceTypes.Absolute)
                    residualNorm = r;
                else
                    residualNorm = r / r_0;
                //Console.WriteLine("JACOBI: " + residualNorm);

                if (m_UnderRelaxationFactor != 1.0)
                    tmp.Scale(m_UnderRelaxationFactor);
                tmp.MultiplyElementWise(InvDiag);

                x.Acc(1.0, tmp);

                stats.NoOfIterations++;
            }


            // unlock
            // ======
            x.Unlock();
            m_Matrix.Unlock();
            rhs.Unlock();
            InvDiag.Unlock();
        }
    }

}