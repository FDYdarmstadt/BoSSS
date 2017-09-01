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
    /// encoding of driver (subclass of <see cref="Device"/>)
    /// </summary>
    public enum DeviceType {

        /// <summary>
        /// automatic selection of driver
        /// </summary>
        Auto,

        /// <summary>
        /// enforce usage of CUDA - device will be (<see cref="CUDA.CudaDevice"/>);
        /// </summary>
        Cuda,

        /// <summary>
        /// enforce usage of OpenCL - device will be (<see cref="CL.clDevice"/>);
        /// </summary>
        OpenCL,

        /// <summary>
        /// enforce usage of CPU - device will be (<see cref="ilPSP.LinSolvers.monkey.CPU.ReferenceDevice"/>);
        /// </summary>
        CPU,

        /// <summary>
        /// enforce usage of CPU - device will be (<see cref="mtCPU.MtDevice"/>);
        /// </summary>
        MultiThreadCPU
    }

    /// <summary>
    /// Data structure of matrix, only used with CUDA
    /// </summary>
    public enum MatrixType {
        /// <summary>
        /// automatic selection of matrix format
        /// </summary>
        Auto,

        /// <summary>
        /// classical CSR, not very perfomant, only recommended if
        /// the number of nonzeros per row varies a lot, see <see cref="MatrixBase.CSR"/>
        /// </summary>
        CSR,

        /// <summary>
        /// block matrix format (Block - CSR), recommended if the matrix contains lots of dense block matrices,
        /// <see cref="MatrixBase.CCBCSR"/>
        /// </summary>
        CCBCSR,

        /// <summary>
        /// classical ELLPACK (see <see cref="MatrixBase.ELLPACKmod"/>),
        /// modified for the GPU, usually <see cref="ELLPACKcache"/>
        /// works better
        /// </summary>
        ELLPACK,

        /// <summary>
        /// ELLPACK with manual caching - ususally best (see <see cref="MatrixBase.ManualCacheELLPACK"/>)
        /// </summary>
        ELLPACKcache
    }

    
    /// <summary>
    /// an abstract factory for creating monkey-vectors (<see cref="VectorBase"/>) and -matrices(<see cref="MatrixBase"/>);
    /// </summary>
    /// <remarks>
    /// solvers use factories of this type to create data objects which may 'live' either on 
    /// main processor or on accelerator (e.g. GPU) devices.
    /// </remarks>
    public abstract class Device {

        /// <summary>
        /// creates an empty vector with MPI-partitioning <paramref name="p"/>
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        abstract public VectorBase CreateVector(IPartitioning p);

        /// <summary>
        /// creates a matrix and loads the contents of <paramref name="M"/>
        /// </summary>
        /// <param name="M"></param>
        /// <param name="matType">
        /// matrix storage format, may be ignored on some implementations (like e.g. <see cref="CPU.ReferenceDevice"/>);
        /// </param>
        /// <returns></returns>
        abstract public MatrixBase CreateMatrix(MsrMatrix M, MatrixType matType);

        /// <summary>
        /// creates a vector with MPI-partitioning <paramref name="p"/>
        /// and loades the content of the list <paramref name="content"/>.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="p"></param>
        /// <param name="content"></param>
        /// <param name="shallowInit">
        /// if true, a shallow copy has been performed, i.e. the returned vector did not  
        /// allocate its own memory, but uses the memory of <paramref name="content"/>.
        /// </param>
        /// <returns></returns>
        abstract public VectorBase CreateVector<T>(IPartitioning p, T content, out bool shallowInit) where T : IList<double>;
    }
}
