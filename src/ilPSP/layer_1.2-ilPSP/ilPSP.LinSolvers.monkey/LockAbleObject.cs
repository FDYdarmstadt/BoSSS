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

namespace ilPSP.LinSolvers.monkey {

    /// <summary>
    /// Baseclass for matrices, vectors and other data objects which can be uploaded into a GPU;
    /// </summary>
    /// <remarks>
    /// Normally, this object (i.e. its content) is stored in the main memory, but a call
    /// to <see cref="Lock"/> copies/clones its content in some mammer and may store it into
    /// the dedicated memory of some accellerator, e.g. a GPU;<br/>
    /// A call to <see cref="Unlock"/> undoes this thing, i.e. downloads the content to main memory
    /// and frees the memory on the accelerator device.
    /// </remarks>
    public abstract class LockAbleObject : IDisposable {

        /// <summary>
        /// <see cref="IsLocked"/>
        /// </summary>
        protected bool m_IsLocked = false;

        /// <summary>
        /// tue is this object is locked, see <see cref="Lock"/> and <see cref="Unlock"/>;
        /// </summary>
        virtual public bool IsLocked {
            get { return m_IsLocked; }
        }

        /// <summary>
        /// puts this object into 'locked' - state; This can be used to implement e.g.
        /// uploading of onto a GPU.
        /// </summary>
        virtual public void Lock() {
            if (m_IsLocked == true)
                throw new ApplicationException("matrix is allready locked.");
            m_IsLocked = true;
        }


        /// <summary>
        /// unlocks this object; When unlocked, the objects properties can 
        /// be altered, which is not possible in locked state
        /// </summary>
        virtual public void Unlock() {
            if (m_IsLocked == false)
                throw new ApplicationException("matrix must be locked.");
            m_IsLocked = false;
        }

        #region IDisposable Members

        /// <summary>
        /// override to support disposal of unmanaged resources.
        /// </summary>
        virtual public void Dispose() {
        }

        #endregion
    }
}
