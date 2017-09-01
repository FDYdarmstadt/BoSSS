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
using ilPSP.Utils;
using System.Linq;

namespace ilPSP.LinSolvers.monkey.mtCPU {
    
    /// <summary>
    /// factory for creating multithreaded implementations: <see cref="MtVector"/> and <see cref="MtMatrix"/>;
    /// </summary>
    public class MtDevice : Device {

        /// <summary>
        /// see <see cref="Device.CreateVector(IPartitioning)"/>;
        /// </summary>
        public override VectorBase CreateVector(IPartitioning p) {
            if (p.IsMutable)
                throw new NotSupportedException();
            return new MtVector(p);
        }

        /// <summary>
        /// see <see cref="Device.CreateMatrix(MsrMatrix,MatrixType)"/>;
        /// </summary>
        public override MatrixBase CreateMatrix(MsrMatrix M, MatrixType matType) {
            return new MtMatrix(M);
        }

        /// <summary>
        /// see <see cref="Device.CreateVector{T}(IPartitioning,T,out bool)"/>;
        /// </summary>
        public override VectorBase CreateVector<T>(IPartitioning p, T content, out bool shallowInit) {
            double[] vals = null;
            if (typeof(T).Equals(typeof(double[]))) {
                shallowInit = true;
                vals = content as double[];
            } else {
                shallowInit = false;
                vals = content.ToArray();
            }
            return new MtVector(p, vals);
        }
    }
}
