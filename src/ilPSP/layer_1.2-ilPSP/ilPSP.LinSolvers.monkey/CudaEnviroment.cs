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
using log4net;

namespace ilPSP.LinSolvers.monkey.CUDA {
    
    
    /// <summary>
    /// CUDA Environment - which MPI process uses which GPU
    /// </summary>
    public class CudaEnviroment : IDisposable {

        static ILog Logger = LogManager.GetLogger(typeof(CudaEnviroment));

        /// <summary>
        /// inits the cuda driver, if possible
        /// </summary>
        static CudaEnviroment() {
            try {
                cu.Init(0);
            } catch (Exception) {
                m_CudaAvailable = false;
            }
            m_CudaAvailable = true;


            // check that all Data Types have the correct size (sometimes the Runtime makes a mess with this): https://bugzilla.novell.com/show_bug.cgi?id=697898
            int d_PtrSz = 4;

            if (Marshal.SizeOf(typeof(CUdevice)) != sizeof(int))
                throw new ApplicationException();

            if (Marshal.SizeOf(typeof(CUcontext)) != IntPtr.Size)
                throw new ApplicationException();

            if (Marshal.SizeOf(typeof(CUmodule)) != IntPtr.Size)
                throw new ApplicationException();

            if (Marshal.SizeOf(typeof(CUfunction)) != IntPtr.Size)
                throw new ApplicationException();

            if (Marshal.SizeOf(typeof(CUstream)) != IntPtr.Size)
                throw new ApplicationException();

            if (Marshal.SizeOf(typeof(CUdeviceptr)) != d_PtrSz)
                throw new ApplicationException();

        }

        static bool m_CudaAvailable = false;

        /// <summary>
        /// true, if this computer has some CUDA driver installed
        /// </summary>
        public static bool CudaAvailable {
            get { return m_CudaAvailable; } 
        }
		
        /// <summary>
        /// ctor
        /// </summary>
        public CudaEnviroment(MPIEnviroment _MpiEnv) {
            m_MpiEnv = _MpiEnv;
            if (!CudaAvailable)
                throw new ApplicationException("No CUDA driver on this machine.");

            int NoOfDev;
            cu.DeviceGetCount(out NoOfDev);

            if (NoOfDev < _MpiEnv.ProcessesOnMySMP)
                throw new ApplicationException("not enough CUDA devices; There must be at least one CUDA device for each MPI process;");
            
            //cu.DeviceGet(out m_CUDAdev, 0); // provisorisch
            //Console.WriteLine("Cuda: prowisoischä dev-allok.");
            cu.DeviceGet(out m_CUDAdev, _MpiEnv.ProcessRankOnSMP); // korrekt
            
            int major, minor;
            cu.DeviceComputeCapability(out major, out minor, m_CUDAdev);
            if((major < 1) || (major <= 1 && minor <= 2)) 
                throw new ApplicationException("requiring at least compute capability 1.3 - got " + major + "." + minor + " for device #" + _MpiEnv.ProcessRankOnSMP + " on host '" + MpiEnv.Hostname + "';");

            cu.CtxCreate(out m_context, CUDAdev, CUctx_flags_enum.CU_CTX_MAP_HOST);
            disposed = false;
            String name;
            cu.DeviceGetName(out name, CUDAdev);

            //uint membytes;
            //cu.DeviceTotalMem(out membytes, m_Env.CUDAdev);
            //membytes /= 1024 * 1024;

            Logger.Info("Process " + MpiEnv.MPI_Rank + " on " + MpiEnv.Hostname + " using " + name + ".");
        }

        MPIEnviroment m_MpiEnv;

        /// <summary>
        /// the used MPI environment
        /// </summary>
        public MPIEnviroment MpiEnv {
            get { return m_MpiEnv; }
        }
        
        CUdevice m_CUDAdev;

        /// <summary>
        /// CUDA device handle for the current process
        /// </summary>
        public CUdevice CUDAdev {
            get { return m_CUDAdev; }
        }


        CUmodule CudaMatrixKernelDP = default(CUmodule);

        SortedDictionary<string, CUfunction> CudaMatrixKernelDP_Functions = null;

        /// <summary>
        /// 
        /// </summary>
        internal CUfunction Get_CudaMatrixKernelDP_Function(string funcname) {
            if (CudaMatrixKernelDP_Functions == null) {
                IntPtr CudaMatrixKernelDP_ptx = Marshal.StringToHGlobalAnsi(Properties.Resources.CudaMatrixKernelDP);
                cu.ModuleLoadData(out CudaMatrixKernelDP, CudaMatrixKernelDP_ptx);
                Marshal.FreeHGlobal(CudaMatrixKernelDP_ptx);

                //cu.ModuleGetFunction(out sparseMultiply, module, funcName);
                //cu.ModuleGetFunction(out cuaccext, module, "accumulateExternal");

                CudaMatrixKernelDP_Functions = new SortedDictionary<string, CUfunction>();
            }

            CUfunction ret;
            if (!CudaMatrixKernelDP_Functions.TryGetValue(funcname, out ret)) {
                cu.ModuleGetFunction(out ret, CudaMatrixKernelDP, funcname);
            }
            return ret;
        }        

        CUmodule CudaVectorKernelDP = default(CUmodule);

        SortedDictionary<string, CUfunction> CudaVectorKernelDP_Functions = null;

        /// <summary>
        /// 
        /// </summary>
        internal CUfunction Get_CudaVectorKernelDP_Function(string funcname) {
            if (CudaVectorKernelDP_Functions == null) {
                IntPtr CudaVectorKernelDP_ptx = Marshal.StringToHGlobalAnsi(Properties.Resources.CudaVectorKernelDP);
                cu.ModuleLoadData(out CudaVectorKernelDP, CudaVectorKernelDP_ptx);
                Marshal.FreeHGlobal(CudaVectorKernelDP_ptx);

                CudaVectorKernelDP_Functions = new SortedDictionary<string, CUfunction>();
            }

            CUfunction ret;
            if (!CudaMatrixKernelDP_Functions.TryGetValue(funcname, out ret)) {
                cu.ModuleGetFunction(out ret, CudaVectorKernelDP, funcname);
            }
            return ret;
        }


        private CUcontext m_context;

        /// <summary>
        /// gets access to the currently used CUDA context handle
        /// </summary>
        public CUcontext Context {
            get { return m_context; }
        }
        

        /// <summary>
        /// dtor
        /// </summary>
        ~CudaEnviroment() {
            //Dispose();
        }

        bool disposed = true;

        /// <summary>
        /// frees the CUDA context
        /// </summary>
        public void Dispose() {
            if (disposed)
                return;
            
            disposed = true;

            cu.ModuleUnload(CudaMatrixKernelDP);
            CudaMatrixKernelDP = default(CUmodule);
            CudaMatrixKernelDP_Functions = null;

            cu.ModuleUnload(CudaVectorKernelDP);
            CudaVectorKernelDP = default(CUmodule);
            CudaVectorKernelDP_Functions = null;

            cu.CtxDestroy(m_context);
            m_context = default(CUcontext);
        }
    }
}
