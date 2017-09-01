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

namespace ilPSP.LinSolvers.monkey {
    
    /// <summary>
    /// Jacobi Preconditioner
    /// </summary>
    public class JacobiPrecond : Precond {
        
        /// <summary>
        ///  "one over diagonal elements" of the matrix
        /// </summary>
        VectorBase CreateInvDiag(IMutableMatrixEx Matrix) {
            VectorBase InvDiag = m_Device.CreateVector(Matrix.RowPartitioning);
            int i0 = (int)Matrix.RowPartitioning.i0;
            int L = (int)Matrix.RowPartitioning.LocalLength;
            for (int i = 0; i < L; i++) {
                InvDiag[i] = 1.0 / Matrix[i + i0,i + i0];
            }
            return InvDiag;
        }

        double m_UnderRelaxationFactor = 1.0;

        /// <summary>
        /// gets/sets the relaxation factor for the Jacobi iteration; Default is 1.0;<br/>
        /// Recommended range is between 0 (excluding) and 1.0 (including, default);
        /// </summary>
        public double UnderRelaxationFactor {
            get { return m_UnderRelaxationFactor; }
            set {
                if (m_UnderRelaxationFactor <= 0)
                    throw new ArgumentOutOfRangeException("must be positive");
                m_UnderRelaxationFactor = value;
            }
        }
        
        #region IMonkeyImplicitPrecond Members

        VectorBase.CommVector _xComm;
        VectorBase tmp;
        
        /// <summary>
        /// containes the inverese of the diagonal elements of <see cref="m_Matrix"/>
        /// </summary>
        VectorBase m_InvDiag;
                
        /// <summary>
        /// input data for the preconditioner
        /// </summary>
        VectorBase m_PcInput;
        
        /// <summary>
        /// output object of the preconditioning 
        /// </summary>
        VectorBase m_PcOutput;
        
        /// <summary>
        /// 
        /// </summary>
        MatrixBase m_Matrix;

        /// <summary>
        /// see <see cref="IMonkeyImplicitPrecond.CreateTempObjects"/>
        /// </summary>
        public override void CreateTempObjects(VectorBase x, VectorBase b, MatrixBase mtx, Device dev) {
            if (!x.IsLocked)
                throw new ArgumentException("x must be locked.", "x");
            if (!b.IsLocked)
                throw new ArgumentException("b must be locked.", "b");
            if (!mtx.IsLocked)
                throw new ArgumentException("mtx must be locked.", "mtx");

            m_Matrix = mtx;
            m_PcInput = b;
            m_PcOutput = x;

            // create objects 
            // ==============

            _xComm = x.CreateCommVector(mtx);
            tmp = m_Device.CreateVector(m_PcOutput.Part);

            // lock objects
            // ============
            tmp.Lock();
            m_InvDiag.Lock();
        }

        /// <summary>
        /// see <see cref="IMonkeyImplicitPrecond.ReleaseTempObjects"/>
        /// </summary>
        public override void ReleaseTempObjects() {
            // unlock
            // ======
            tmp.Unlock();
            tmp.Dispose();
            tmp = null;

            m_InvDiag.Unlock();

            _xComm.Dispose();
            _xComm = null;
        }

        /// <summary>
        /// see <see cref="IMonkeyImplicitPrecond.DoPrecond"/>
        /// </summary>
        override public void DoPrecond() {
            
            m_PcOutput.CopyFrom(m_PcInput);

            // iterate
            // =======
            for( int i = 1; i < m_NoOfIter; i++) {
                                
                // Jacobi iteration
                // ================
                //((CPU.RefMatrix)m_Matrix).useSingle = true;
                m_Matrix.SpMV_Expert(-1.0, _xComm, 0.0, tmp); // tmp = -M*x
                //((CPU.RefMatrix)m_Matrix).useSingle = true;
                tmp.Acc(1.0, m_PcInput);                      // tmp = -M*x + m_PcInput
                if (m_UnderRelaxationFactor != 1.0)
                    tmp.Scale(m_UnderRelaxationFactor);

                
                if (m_UnderRelaxationFactor != 1.0)
                    tmp.Scale(m_UnderRelaxationFactor);
                tmp.MultiplyElementWise(m_InvDiag);

                m_PcOutput.Acc(1.0, tmp);               //m_PcOutput = m_PcOutput + tmp
            }
        }
        
        int m_NoOfIter = 3;

        /// <summary>
        /// fixed number of Jacobi itearations in all preconditioner steps
        /// </summary>
        /// <remarks>
        /// must be greater or equal to 1; With one iteration, 
        /// this preconditioner is just the identity matrix.
        /// </remarks>
        public int NoOfIter {
            get { return m_NoOfIter; }
            set {
                if (value <= 1)
                    throw new ArgumentOutOfRangeException("number of iterations must be greater or equal to 1.");
                m_NoOfIter = value; 
            }
        }

        /// <summary>
        /// see <see cref="IMonkeyImplicitPrecond.Initialize"/>
        /// </summary>
        public override void Initialize(IMutableMatrixEx OrigMatrix, Device dev, MatrixType mt) {
            base.Initialize(OrigMatrix, dev, mt);
            m_InvDiag = CreateInvDiag(OrigMatrix);
        }

        #endregion
    }
}
