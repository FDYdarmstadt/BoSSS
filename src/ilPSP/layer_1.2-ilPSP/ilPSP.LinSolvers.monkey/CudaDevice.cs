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
using System.IO;
using System.Reflection;
using ilPSP.Utils;
using ilPSP.Tracing;
using System.Linq;

namespace ilPSP.LinSolvers.monkey.CUDA {

    /// <summary>
    /// a monkey device that uses CUDA
    /// </summary>
    /// <remarks>
    /// All objects
    /// (matrices: <see cref="CudaMatrix"/> and derivatives, vectors: <see cref="CudaVector"/>) are created within 
    /// this CUDA context, and therefore, (methods of these objects) must be called from the same thread. <br/>
    /// Using multiple <see cref="CudaDevice"/>-objects in one thread is neither tested nor supported.
    /// </remarks>
    public class  CudaDevice : Device {


        private CudaEnviroment m_Env;

        /// <summary>
        /// CUDA Environment (use which GPU on which MPI process, ...);
        /// </summary>
        public CudaEnviroment Env {
            get { return m_Env; }
        }


        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="cuEnv">Distribution of processes and CUDA devices</param>
        public CudaDevice(CudaEnviroment cuEnv) {
            m_Env = cuEnv;
        }


        /// <summary>
        /// 
        /// </summary>
        public override VectorBase CreateVector(IPartitioning p) {
            if (p.IsMutable)
                throw new NotSupportedException();
            return new CudaVector(p,m_Env);
        }

        /// <summary>
        /// see <see cref="Device.CreateMatrix(MsrMatrix,MatrixType)"/>;
        /// </summary>
        public override MatrixBase CreateMatrix(MsrMatrix M, MatrixType matType) {
            using (var f = new FuncTrace()) {
                f.Info("desired matrix type: " + matType);
                switch (matType) {
                    case MatrixType.CCBCSR:
                        return new CudaCCBCSRMatrix(M, this.Env);
                    case MatrixType.CSR:
                        return new CudaCSRMatrix(M, this.Env);
                    case MatrixType.ELLPACK:
                        return new CudaELLPACKmodMatrix(M, this.Env);
                    case MatrixType.Auto:
                    case MatrixType.ELLPACKcache:
                        return new CudaELLPACKcacheMatrix(M, this.Env);
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// see <see cref="Device.CreateVector{T}(IPartitioning,T,out bool)"/>
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
            return new CudaVector(p, vals, m_Env);
        }
    }
}
