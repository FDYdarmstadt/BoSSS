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
    /// baseclass for all preconditioners in monkey 
    /// </summary>
    abstract public class Precond : IMonkeyImplicitPrecond {
        #region IMonkeyImplicitPrecond Members

        /// <summary>
        /// see <see cref="IMonkeyImplicitPrecond.Initialize"/>
        /// </summary>
        virtual public void Initialize(IMutableMatrixEx OrigMatrix, Device dev, MatrixType mt) {
            m_Device = dev;
            m_MatrixType = mt;
        }

        /// <summary>
        /// factory for creating matrix and vector objects;
        /// </summary>
        protected Device m_Device;

        /// <summary>
        /// 'recommendation' for the storage format of matrices that are created by this
        /// preconditioner
        /// </summary>
        protected MatrixType m_MatrixType;

        /// <summary>
        /// called after by the hosting solver before it starts, 
        /// to provide configuration and the possibility th create temporary objects
        /// </summary>
        /// <param name="pc_output">
        /// Here, the hosting solver expects the result of the precondioning;
        /// see <see cref="DoPrecond"/>;
        /// </param>
        /// <param name="pc_input">
        /// Here, the hosting solver will place the input vector for the precondioning;
        /// see <see cref="DoPrecond"/>;
        /// </param>
        /// <param name="mtx"></param>
        /// <param name="dev"></param>
        abstract public void CreateTempObjects(VectorBase pc_output, VectorBase pc_input, MatrixBase mtx, Device dev);

        /// <summary>
        /// called after by the hosting solver after it is finished, to free temporary objects 
        /// </summary>
        abstract public void ReleaseTempObjects();

        /// <summary>
        /// Computes the preconditioning, i.e. \f$ y = M \cdot x\f$ ,
        /// with the input (<em>x</em>) and output vectors (<em>y</em>)
        /// specified
        /// in <see cref="CreateTempObjects"/>;
        /// </summary>
        abstract public void DoPrecond();

        #endregion
    }
}
