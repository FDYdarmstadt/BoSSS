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
    /// This class represents Neumann Preconditioner implementing the abstract class Precond
    /// </summary>
    class NeumannSeries : Precond {

        //the number of iterations for calculating the neumann series
        int m_NoOfIter = 3;
        VectorBase.CommVector _yComm;

        VectorBase mPcOutput;
        VectorBase mPcInput;
        MatrixBase m_matrix;

        VectorBase yOld;
        VectorBase yNew;

        
        /// <summary>
        /// the number of iterations used for calculating the neumann series
        /// </summary>
        public int NoOfIter {
            get { return m_NoOfIter; }
            set { m_NoOfIter = value; }
        }
                
        /// <summary>
        /// The used temporary objects are unlocked and set to null 
        /// </summary>
        public override void ReleaseTempObjects() {

            //unlock
            yOld.Unlock();
            yOld.Dispose();
            yOld = null;

            _yComm.Dispose();
            _yComm = null;

            yNew.Unlock();
            yNew.Dispose();
            yNew = null;
        }

        /// <summary>
        /// this function implements the preconditioning according to the neumann preconditioner
        /// </summary>
        public override void DoPrecond() {

            mPcOutput.CopyFrom(mPcInput);
           
            for (int i = 1; i < m_NoOfIter; i++) {
                //_yNew = _yOld - m_matrix*_yOld;
                m_matrix.SpMV_Expert(-1.0, _yComm, 0.0, yNew);  //yNew = -m_matrix * yComm = -m_Matrix*yOld
                yNew.Acc(1.0, yOld);                            //yNew = yOld - m_Matrix*yOld;
                
                mPcOutput.Acc(1.0, yNew);                       //mPcOutput = mPcOutput + yNew
                yOld.Swap(yNew);
               
            }
        }

        /// <summary>
        /// This function creates the needed temporary objects according to the given parameters
        /// </summary>
        /// <param name="pc_output">The output vector</param>
        /// <param name="pc_input">The input vector</param>
        /// <param name="mtx">The matrix that is to be multiplied</param>
        /// <param name="dev"></param>
        public override void CreateTempObjects(VectorBase pc_output, VectorBase pc_input, MatrixBase mtx, Device dev) {

            if (!pc_output.IsLocked)
                throw new ArgumentException("pc_output must be locked.", "pc_output");
            if (!pc_input.IsLocked)
                throw new ArgumentException("pc_input must be locked.", "pc_input");
            if (!mtx.IsLocked)
                throw new ArgumentException("mtx must be locked.", "mtx");

            mPcInput = pc_input;
            mPcOutput = pc_output;
            m_matrix = mtx;

            //the temporary objects
            yOld = m_Device.CreateVector(mPcInput.Part);
            yNew = m_Device.CreateVector(mPcInput.Part);
            _yComm = yOld.CreateCommVector(m_matrix);
            //lock the temporary objects:
            yOld.Lock();
            yNew.Lock();
        }
    }
}
