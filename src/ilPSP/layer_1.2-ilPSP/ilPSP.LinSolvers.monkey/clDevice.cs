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

namespace ilPSP.LinSolvers.monkey.CL
{
    /// <summary>
    /// OpenCL device
    /// </summary>
    public class clDevice : Device, IDisposable
    {
        private bool disposed = false;

        internal clEnvironment env;
        internal cl_command_queue cq;
        internal cl_program matrixProgram, vectorProgram;


        /// <summary>
        /// Create OpenCL device
        /// </summary>
        /// <param name="clEnv">OpenCL environment</param>
        public clDevice(clEnvironment clEnv)
        {
            env = clEnv;

            cq = cl.CreateCommandQueue(env.context, env.device);
            vectorProgram = loadProgram(clVectorSource.source);
            matrixProgram = loadProgram(clMatrixSource.source);
        }

        /// <summary>
        /// Destructor
        /// </summary>
        ~clDevice()
        {
            Dispose();
        }

        private cl_program loadProgram(string[] src)
        {
            
            cl_program program = cl.CreateProgramWithSource(env.context, src);
            try
            {
                cl.BuildProgram(program);
            }
            catch (ApplicationException e)
            {
                cl_program_build_info_return info = cl.GetProgramBuildInfo(program, env.device);
                Console.WriteLine(info.status);
                Console.WriteLine(info.log);
                throw e;
            }
            return program;
        }

        /// <summary>
        /// Free resources
        /// </summary>
        public void Dispose()
        {
            if (!disposed)
            {
                cl.ReleaseProgram(matrixProgram);
                cl.ReleaseProgram(vectorProgram);
                cl.ReleaseCommandQueue(cq);
                disposed = true;
            }
        }

        /// <summary>
        /// Create vector for this device
        /// </summary>
        /// <param name="p">Parition</param>
        /// <returns>The vector</returns>
        public override VectorBase CreateVector(IPartitioning p) {
            if (p.IsMutable)
                throw new NotSupportedException();
            if (disposed)
                throw new ApplicationException("object is disposed.");

            return new clVector(p, this);
        }

        /// <summary>
        /// Create vector for this device
        /// </summary>
        /// <typeparam name="T">Type for data storage</typeparam>
        /// <param name="p">Partition</param>
        /// <param name="content">Data storage</param>
        /// <param name="shallowInit">Data is copied or referenced</param>
        /// <returns>The vector</returns>
        public override VectorBase CreateVector<T>(IPartitioning p, T content, out bool shallowInit) {
            if (p.IsMutable)
                throw new NotSupportedException();
            if (disposed)
                throw new ApplicationException("object is disposed.");

            double[] vals = null;
            if (typeof(T).Equals(typeof(double[]))) {
                shallowInit = true;
                vals = content as double[];
            } else {
                shallowInit = false;
                vals = content.ToArray();
            }
            return new clVector(p, vals, this);
        }

        /// <summary>
        /// 
        /// </summary>
        public override MatrixBase CreateMatrix(MsrMatrix M, MatrixType mt)
        {
            if (disposed)
                throw new ApplicationException("object is disposed.");

            switch (mt)
            {
                case MatrixType.CCBCSR:
                    return new clCCBCSRMatrix(M, this);
                case MatrixType.CSR:
                    return new clCSRMatrix(M, this);
                case MatrixType.ELLPACK:
                    return new clELLPACKmodMatrix(M, this);
                case MatrixType.Auto:
                case MatrixType.ELLPACKcache:
                    return new clELLPACKcacheMatrix(M, this);
                default:
                    throw new NotImplementedException();
            }
        }
    }
}
