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
using System.Runtime.InteropServices;
using MPI.Wrappers;
using System.Diagnostics;

namespace ilPSP.LinSolvers.monkey {

    /// <summary>
    /// solver diagnostic
    /// </summary>
    /// <param name="IterNo">
    /// iteration number
    /// </param>
    /// <param name="x_approx">
    /// approximate solution so far
    /// </param>
    /// <param name="res_curr">
    /// residual in current iteration step
    /// </param>
    /// <remarks>
    /// usually, this feature causes massive overhead and should be used only for diagnostic purpose
    /// </remarks>
    public delegate void IterationCallback(int IterNo, double[] x_approx, double[] res_curr);
        
    /// <summary>
    /// a preconditioned conjugated gradient solver
    /// </summary>
    public class PCG : Solver, IImplicitPrecondSupport {
        
        IterationCallback m_IterationCallback = null;
        
        /// <summary>
        /// if unequal to null, the residual and current approximative
        /// solution is passed to the callback function;
        /// This feature causes massive overhead and should be used only for diagnostic purpose;
        /// </summary>
        public IterationCallback IterationCallback {
            get { return m_IterationCallback; }
            set { m_IterationCallback = value; }
        }


        /// <summary>
        /// implementation of the CG algorithm
        /// </summary>
        /// <param name="x"></param>
        /// <param name="rhs"></param>
        /// <param name="stats"></param>
        protected override void CallSolver(VectorBase x, VectorBase rhs, ref SolverResult stats) {

            VectorBase P = Device.CreateVector(x.Part);
            VectorBase.CommVector commP = P.CreateCommVector(m_Matrix);

            VectorBase R = rhs; // rhs is only needed once, so we can use it to store residuals
            VectorBase V = Device.CreateVector(x.Part);
            VectorBase Z = Device.CreateVector(x.Part);

            // lock objects
            // ============
            x.Lock();
            P.Lock();
            R.Lock();
            V.Lock();
            Z.Lock();

            m_Matrix.Lock();


            // configure Precond
            // =================

            if (m_NestedPrecond != null) {
                m_NestedPrecond.CreateTempObjects(Z, R, m_Matrix,this.Device);
            }
            

            // compute P0, R0
            // ==============
            P.Swap(x);
            m_Matrix.SpMV_Expert(-1.0, commP, 1.0, R);
            P.Swap(x);
            if (m_NestedPrecond != null) {
                m_NestedPrecond.DoPrecond();
                P.CopyFrom(Z);
            } else {
                P.CopyFrom(R);
            }

            double alpha = R.InnerProd(P);
            double alpha_0 = alpha;
            double ResNorm;

            if (m_ConvergenceType == ConvergenceTypes.Absolute)
                ResNorm = Math.Sqrt(alpha);
            else if (m_ConvergenceType == ConvergenceTypes.Relative)
                ResNorm = 1.0;
            else {
                throw new NotImplementedException("unknown convergence type: " + m_ConvergenceType.ToString());
            }

            //long total = 0;
            //long gemv = 0;
            //long rest = 0;
            //long st, en;

            // iterate
            // =======
            stats.Converged = false;
            stats.NoOfIterations = 1; // one iteration has allready been performed (P0, R0)
            for (int n = m_MaxIterations - 2; n >= 0; n--) {

                if (ResNorm <= m_Tolerance && stats.NoOfIterations >= base.m_MinIterations) {
                    stats.Converged = true;
                    break;
                }

                if (Math.Abs(alpha) <= double.Epsilon)
                    // numerical breakdown
                    break;

                m_Matrix.SpMV_Expert(1.0, commP, 0, V);
                double lambda = alpha / V.InnerProd(P);

                x.Acc(lambda, P);
                
                R.Acc(-lambda, V);

                if (m_IterationCallback != null) {
                    // pass approx. sol and residual to callback function

                    x.Unlock();
                    R.Unlock();

                    double[] x_approx = new double[x.Part.LocalLength];
                    x.CopyTo(x_approx, 0);

                    double[] R_curr = new double[R.Part.LocalLength];
                    R.CopyTo(R_curr, 0);

                    m_IterationCallback(stats.NoOfIterations, x_approx, R_curr);


                    x.Lock();
                    R.Lock();

                }

                if (m_NestedPrecond != null) {
                    Z.Clear();
                    m_NestedPrecond.DoPrecond();
                } else {
                    Z.CopyFrom(R);
                }

                double alpha_neu = R.InnerProd(Z);
                
                // compute residual norm
                if (m_ConvergenceType == ConvergenceTypes.Absolute)
                    ResNorm = Math.Sqrt(alpha);
                else
                    ResNorm = Math.Sqrt(alpha/alpha_0);
                ResNorm = Math.Sqrt(R.TwoNormSquare());
                
                P.Scale(alpha_neu / alpha);
                P.Acc(1.0, Z);
                                
                alpha = alpha_neu;
                stats.NoOfIterations++;
            }

            // unlock objects
            // ==============

            if (m_NestedPrecond != null)
                m_NestedPrecond.ReleaseTempObjects();

            x.Unlock();
            P.Unlock();
            R.Unlock();
            V.Unlock();

            m_Matrix.Unlock();
            
            commP.Dispose();
            P.Dispose();
        }

        #region IImplicitPrecondSupport Members

        /// <summary>
        /// see <see cref="IImplicitPrecondSupport.SupportedPrecond"/>
        /// </summary>
        public ICollection<Type> SupportedPrecond {
            get { return new Type[] { typeof(IMonkeyImplicitPrecond) }; }
        }

        IMonkeyImplicitPrecond m_NestedPrecond;

        /// <summary>
        /// see <see cref="IImplicitPrecondSupport.NestedPrecond"/>
        /// </summary>
        public IImplicitPrecond NestedPrecond {
            get {
                return m_NestedPrecond;
            }
            set {
                m_NestedPrecond = (IMonkeyImplicitPrecond)value;
            }
        }

        #endregion

        /// <summary>
        /// see <see cref="Solver.DefineMatrix"/>
        /// </summary>
        public override void DefineMatrix(IMutableMatrixEx M) {
            base.DefineMatrix(M);
            if (m_NestedPrecond != null)
                m_NestedPrecond.Initialize(M, base.Device, base.MatrixType);
        }



    }
}
