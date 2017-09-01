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
using MPI.Wrappers.Utils;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Threading;

namespace ilPSP.LinSolvers.monkey.CUDA
{

    //public class Kaos {
    //    [DllImport("dl", CharSet = CharSet.Ansi)]
    //    static extern IntPtr dlopen(string filename, int flag);
    //    [DllImport("dl", CharSet = CharSet.Ansi)]
    //    static extern IntPtr dlsym(IntPtr handle, string symbol);
    //    private delegate int _cuInit(uint Flags);
    //    public static void La() {
    //        IntPtr p = dlopen("libcuda.so", 1);
    //        IntPtr f = dlsym(p, "cuInit");
    //        _cuInit init = (_cuInit)Marshal.GetDelegateForFunctionPointer(f, typeof(_cuInit));
    //        Console.Write("cuinit res: ");
    //        Console.Write(init(0));
    //    }
    //}


    #region CUDA structures and data types
    /// <summary> see CUDA doc </summary>
    //[StructLayout(LayoutKind.Explicit)] // Explicit layout seems to produce datatypes of wrong size in Mono 2.10.1: https://bugzilla.novell.com/show_bug.cgi?id=697898
    [StructLayout(LayoutKind.Sequential)]
    public struct CUdevice {
        //[FieldOffset(0)]
        /// <summary> % </summary>
        public int d;
    }

    /// <summary> see CUDA doc </summary>
    //[StructLayout(LayoutKind.Explicit)] // Explicit layout seems to produce datatypes of wrong size in Mono 2.10.1: https://bugzilla.novell.com/show_bug.cgi?id=697898
    [StructLayout(LayoutKind.Sequential)]
    public struct CUcontext {
        //[FieldOffset(0)]
        /// <summary> % </summary>
        public IntPtr p;
    }

    /// <summary> see CUDA doc </summary>
    //[StructLayout(LayoutKind.Explicit)] // Explicit layout seems to produce datatypes of wrong size in Mono 2.10.1: https://bugzilla.novell.com/show_bug.cgi?id=697898
    [StructLayout(LayoutKind.Sequential)]
    public struct CUmodule {
        //[FieldOffset(0)]
        /// <summary> % </summary>
        public IntPtr p;
    }

    /// <summary> see CUDA doc </summary>
    //[StructLayout(LayoutKind.Explicit)] // Explicit layout seems to produce datatypes of wrong size in Mono 2.10.1: https://bugzilla.novell.com/show_bug.cgi?id=697898
    [StructLayout(LayoutKind.Sequential)]
    public struct CUfunction {
        //[FieldOffset(0)]
        /// <summary> % </summary>
        public IntPtr p;
    }

    /// <summary> see CUDA doc </summary>
    //[StructLayout(LayoutKind.Explicit)] // Explicit layout seems to produce datatypes of wrong size in Mono 2.10.1: https://bugzilla.novell.com/show_bug.cgi?id=697898
    [StructLayout(LayoutKind.Sequential)]
    public struct CUstream {
        //[FieldOffset(0)]
        /// <summary> % </summary>
        public IntPtr p;
    }

    /// <summary> see CUDA doc </summary>
    [StructLayout(LayoutKind.Sequential)]
    //[StructLayout(LayoutKind.Explicit)] // Explicit layout seems to produce datatypes of wrong size in Mono 2.10.1: https://bugzilla.novell.com/show_bug.cgi?id=697898
    public struct CUdeviceptr
    {
        //[FieldOffset(0)]
        /// <summary> % </summary>
        public uint p;
    }

    /// <summary>
    /// CUDA device attributes
    /// </summary>
    public enum CUdevice_attribute {
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 1,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X = 2,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y = 3,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z = 4,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X = 5,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y = 6,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z = 7,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK = 8,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_SHARED_MEMORY_PER_BLOCK = 8,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY = 9,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_WARP_SIZE = 10,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_PITCH = 11,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK = 12,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_REGISTERS_PER_BLOCK = 12,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_CLOCK_RATE = 13,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT = 14,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_GPU_OVERLAP = 15,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT = 16,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT = 17,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_INTEGRATED = 18,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY = 19,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_COMPUTE_MODE = 20,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_WIDTH = 21,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_WIDTH = 22,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_HEIGHT = 23,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_WIDTH = 24,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_HEIGHT = 25,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_DEPTH = 26,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_WIDTH = 27,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_HEIGHT = 28,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_NUMSLICES = 29,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_SURFACE_ALIGNMENT = 30,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS = 31,
        /// <summary> na was wohl? </summary>
        CU_DEVICE_ATTRIBUTE_ECC_ENABLED = 32 
    }

    /// <summary>
    /// CUDA driver success/error codes
    /// </summary>
    public enum CUresult {
        /// <summary>No errors</summary>
        CUDA_SUCCESS = 0,
        /// <summary>Invalid value </summary>
        CUDA_ERROR_INVALID_VALUE = 1,
        /// <summary> Out of memory</summary>
        CUDA_ERROR_OUT_OF_MEMORY = 2,
        /// <summary>Driver not initialized </summary>
        CUDA_ERROR_NOT_INITIALIZED = 3,
        /// <summary>Driver de-initialized </summary>
        CUDA_ERROR_DEINITIALIZED = 4,
        /// <summary>No CUDA-capable device available </summary>
        CUDA_ERROR_NO_DEVICE = 100,
        /// <summary>Invalid device </summary>
        CUDA_ERROR_INVALID_DEVICE = 101,
        /// <summary>Invalid kernel image </summary>
        CUDA_ERROR_INVALID_IMAGE = 200,
        /// <summary>Invalid context </summary>
        CUDA_ERROR_INVALID_CONTEXT = 201,
        /// <summary>Context already current </summary>
        CUDA_ERROR_CONTEXT_ALREADY_CURRENT = 202,
        /// <summary>Map failed </summary>
        CUDA_ERROR_MAP_FAILED = 205,
        /// <summary>Unmap failed </summary>
        CUDA_ERROR_UNMAP_FAILED = 206,
        /// <summary>Array is mapped </summary>
        CUDA_ERROR_ARRAY_IS_MAPPED = 207,
        /// <summary>Already mapped </summary>
        CUDA_ERROR_ALREADY_MAPPED = 208,
        /// <summary>No binary for GPU </summary>
        CUDA_ERROR_NO_BINARY_FOR_GPU = 209,
        /// <summary>Already acquired </summary>
        CUDA_ERROR_ALREADY_ACQUIRED = 210,
        /// <summary>Not mapped </summary>
        CUDA_ERROR_NOT_MAPPED = 211,
        /// <summary>Mapped resource not available for access as an array</summary>
        CUDA_ERROR_NOT_MAPPED_AS_ARRAY = 212,
        /// <summary>Mapped resource not available for access as a pointer </summary>
        CUDA_ERROR_NOT_MAPPED_AS_POINTER = 213,
        /// <summary>Uncorrectable ECC error detected </summary>
        CUDA_ERROR_ECC_UNCORRECTABLE = 214,
        /// <summary>CUlimit not supported by device </summary>
        CUDA_ERROR_UNSUPPORTED_LIMIT = 215,
        /// <summary>Invalid source </summary>
        CUDA_ERROR_INVALID_SOURCE = 300,
        /// <summary>File not found </summary>
        CUDA_ERROR_FILE_NOT_FOUND = 301,
        /// <summary>Link to a shared object failed to resolve </summary>
        CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND = 302,
        /// <summary>Shared object initialization failed </summary>
        CUDA_ERROR_SHARED_OBJECT_INIT_FAILED = 303,
        /// <summary>Invalid handle </summary>
        CUDA_ERROR_INVALID_HANDLE = 400,
        /// <summary>Not found </summary>
        CUDA_ERROR_NOT_FOUND = 500,
        /// <summary>CUDA not ready </summary>
        CUDA_ERROR_NOT_READY = 600,
        /// <summary>Launch failed </summary>
        CUDA_ERROR_LAUNCH_FAILED = 700,
        /// <summary>Launch exceeded resources </summary>
        CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES = 701,
        /// <summary> </summary>
        CUDA_ERROR_LAUNCH_TIMEOUT = 702,
        /// <summary>Launch with incompatible texturing </summary>
        CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING = 703,
        /// <summary> Attempted to retrieve 64-bit pointer via 32-bit API function</summary>
        CUDA_ERROR_POINTER_IS_64BIT = 800,
        /// <summary> Attempted to retrieve 64-bit size via 32-bit API function </summary>
        CUDA_ERROR_SIZE_IS_64BIT = 801,
        /// <summary> Unknown error</summary>
        CUDA_ERROR_UNKNOWN = 999
    }

    /// <summary>
    /// Function properties 
    /// </summary>
    public enum CUfunction_attribute {
        /**
         * The number of threads beyond which a launch of the function would fail.
         * This number depends on both the function and the device on which the
         * function is currently loaded.
         */
        CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 0,

        /**
         * The size in bytes of statically-allocated shared memory required by
         * this function. This does not include dynamically-allocated shared
         * memory requested by the user at runtime.
         */
        CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES = 1,

        /**
         * The size in bytes of user-allocated constant memory required by this
         * function.
         */
        CU_FUNC_ATTRIBUTE_CONST_SIZE_BYTES = 2,

        /**
         * The size in bytes of thread local memory used by this function.
         */
        CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES = 3,

        /**
         * The number of registers used by each thread of this function.
         */
        CU_FUNC_ATTRIBUTE_NUM_REGS = 4,

        /**
         * The PTX virtual architecture version for which the function was compiled.
         */
        CU_FUNC_ATTRIBUTE_PTX_VERSION = 5,

        /**
         * The binary version for which the function was compiled.
         */
        CU_FUNC_ATTRIBUTE_BINARY_VERSION = 6,

        /// <summary>
        /// dummy
        /// </summary>
        CU_FUNC_ATTRIBUTE_MAX
    }

    /// <summary> see CUDA doc; </summary>
    public enum CUfunc_cache {
        /// <summary> see CUDA doc; </summary>
        CU_FUNC_CACHE_PREFER_NONE = 0x00,
        /// <summary> see CUDA doc; </summary>
        CU_FUNC_CACHE_PREFER_SHARED = 0x01,
        /// <summary> see CUDA doc; </summary>
        CU_FUNC_CACHE_PREFER_L1 = 0x02
    }

    /// <summary> see CUDA doc; </summary>
    public enum CUmem_host_alloc {
        /// <summary> see CUDA doc; </summary>
        CU_MEMHOSTALLOC_PORTABLE = 0x01,
        /// <summary> see CUDA doc; </summary>
        CU_MEMHOSTALLOC_DEVICEMAP = 0x02,
        /// <summary> see CUDA doc; </summary>
        CU_MEMHOSTALLOC_WRITECOMBINED = 0x04
    }

    /// <summary> see CUDA doc; </summary>
    public enum CUctx_flags_enum {
        /// <summary> Automatic scheduling; </summary>
        CU_CTX_SCHED_AUTO  = 0,
        /// <summary> Set spin as default scheduling </summary>
        CU_CTX_SCHED_SPIN  = 1,
        /// <summary> Set yield as default scheduling </summary>
        CU_CTX_SCHED_YIELD = 2,
        /// <summary> Scheduling Mask </summary>
        CU_CTX_SCHED_MASK  = 0x3,
        /// <summary> Use blocking synchronization </summary>
        CU_CTX_BLOCKING_SYNC = 4,
        /// <summary> Support mapped pinned allocations </summary>
        CU_CTX_MAP_HOST = 8,
        /// <summary> Keep local memory allocation after launch </summary>
        CU_CTX_LMEM_RESIZE_TO_MAX = 16,
        /// <summary> Flags Mask </summary>
        CU_CTX_FLAGS_MASK = 0x1f
    }
    #endregion

    /// <summary>
    /// wrappers to cuda driver interface
    /// </summary>
    public class cu : DynLibLoader
    {
	
		//[MethodImpl(MethodImplOptions.NoInlining)]
		//static void PoorManDebug() {
        //    string _name ;
        //    {
        //        StackFrame fr = new StackFrame(2, true);
        //        StackTrace st = new StackTrace(fr);
		//
        //        _MethodBase m = fr.GetMethod();
        //        _name = m.DeclaringType.FullName + "." + m.Name;
        //    }
		//	
		//	Console.WriteLine(_name);
		//}
		
		

        public cu() : base(
              new string[] { "nvcuda.dll", "libcuda.so" },
              new string[2][][],
              new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity },
              new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
              new int[] { -1, -1 })
        {

        }

        private static cu my;// = new cu();

        static bool m_CudaSupported = true;

        /// <summary>
        /// true if this system supports CUDA
        /// </summary>
        public static bool CudaSupported {
            get { return cu.m_CudaSupported; }
        }
        
        /// <summary>
        /// tries to load CUDa library
        /// </summary>
        static cu() {
            try {
                my = new cu();
                m_CudaSupported = true;
            } catch (Exception ) {
                m_CudaSupported = false;
            }
        }

		[MethodImpl(MethodImplOptions.NoInlining)]
        private static void testResult(CUresult r)
        {
			//Console.WriteLine("thr = '" + Thread.CurrentThread.Priority + "'");
            //{
            //    StackFrame fr = new StackFrame(1, true);
			
            //    _MethodBase m = fr.GetMethod();
            //    string _name = m.DeclaringType.FullName + "." + m.Name;
            //    Console.WriteLine(_name);
            //}
			if (r != CUresult.CUDA_SUCCESS)
                throw new ApplicationException("CUDA error: " + r.ToString());
        }

#pragma warning disable        649
        // ------------------------------------------
        #region Initialization
        private delegate CUresult _cuInit(uint Flags);

        _cuInit cuInit;

        /// <summary> see CUDA doc; </summary>
        public static void Init(uint Flags)
        {
            testResult(my.cuInit(Flags));
        }
        #endregion
        // ------------------------------------------
        #region Device Management
        private delegate CUresult _cuDeviceGetCount(out int count);
        private delegate CUresult _cuDeviceGet(out CUdevice device, int ordinal);
        private delegate CUresult _cuDeviceGetName(IntPtr name, int len, CUdevice device);
        private delegate CUresult _cuDeviceComputeCapability(out int major, out int minor, CUdevice dev);
        private delegate CUresult _cuDeviceTotalMem(out uint bytes, CUdevice dev);
        private delegate CUresult _cuDeviceGetAttribute( out int pi, CUdevice_attribute attrib, CUdevice dev);


        _cuDeviceGetCount cuDeviceGetCount;
        _cuDeviceGet cuDeviceGet;
        _cuDeviceGetName cuDeviceGetName;
        _cuDeviceComputeCapability cuDeviceComputeCapability;
        _cuDeviceTotalMem cuDeviceTotalMem;
        _cuDeviceGetAttribute cuDeviceGetAttribute;

        /// <summary> see CUDA doc; </summary>
        public static void DeviceGetAttribute(out int pi, CUdevice_attribute attrib, CUdevice dev) {
            testResult(my.cuDeviceGetAttribute(out pi,attrib, dev));
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceGetCount(out int count) {
            testResult(my.cuDeviceGetCount(out count));
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceGet(out CUdevice device, int ordinal)
        {
            testResult(my.cuDeviceGet(out device, ordinal));
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceGetName(out string name, CUdevice dev)
        {
            int maxlen = 1024;
            IntPtr mem = Marshal.AllocHGlobal(maxlen);
            CUresult res;
            res = my.cuDeviceGetName(mem, maxlen, dev);
            name = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            testResult(res);
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceComputeCapability(out int major, out int minor, CUdevice device)
        {
            testResult(my.cuDeviceComputeCapability(out major, out minor, device));
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceTotalMem(out uint bytes, CUdevice dev)
        {
            testResult(my.cuDeviceTotalMem(out bytes, dev));
        }
        #endregion
        // ------------------------------------------
        #region Context Management
        private delegate CUresult _cuCtxCreate(out CUcontext pctx, uint flags, CUdevice dev);
        private delegate CUresult _cuCtxDestroy(CUcontext pctx);
        private delegate CUresult _cuCtxSynchronize();

        _cuCtxCreate cuCtxCreate;
        _cuCtxDestroy cuCtxDestroy;
        _cuCtxSynchronize cuCtxSynchronize;

        /// <summary> see CUDA doc; </summary>
        public static void CtxCreate(out CUcontext pctx, CUdevice dev, params CUctx_flags_enum[] flags)
        {
            uint _flags = 0;
            foreach (CUctx_flags_enum x in flags)
            {
                _flags |= (uint)x;
            }
            testResult(my.cuCtxCreate(out pctx, _flags, dev));
        }

        /// <summary> see CUDA doc; </summary>
        public static void CtxDestroy(CUcontext ctx)
        {
            testResult(my.cuCtxDestroy(ctx));
        }

        /// <summary> see CUDA doc; </summary>
        public static void CtxSynchronize()
        {
            testResult(my.cuCtxSynchronize());
        }
        #endregion
        // ------------------------------------------
        #region Module Management
        private delegate CUresult _cuModuleLoad(out CUmodule module, IntPtr fname);
        private delegate CUresult _cuModuleLoadData(out CUmodule module, IntPtr img);
        private delegate CUresult _cuModuleGetFunction(out CUfunction hfunc, CUmodule hmodule, IntPtr name);
        private delegate CUresult _cuModuleUnload(CUmodule module);

        _cuModuleLoad cuModuleLoad;
        _cuModuleLoadData cuModuleLoadData;
        _cuModuleGetFunction cuModuleGetFunction;
        _cuModuleUnload cuModuleUnload;

        /// <summary> see CUDA doc; </summary>
        public static void ModuleUnload(CUmodule module) {
            testResult(my.cuModuleUnload(module));
        }

        /// <summary> see CUDA doc; </summary>
        public static void ModuleLoad(out CUmodule module, string fname)
        {
            IntPtr _fname = Marshal.StringToHGlobalAnsi(fname);
            CUresult res = my.cuModuleLoad(out module, _fname);
            Marshal.FreeHGlobal(_fname);
            testResult(res);
        }

        /// <summary> see CUDA doc; </summary>
        public static void ModuleLoadData(out CUmodule module, IntPtr img)
        {
            CUresult res = my.cuModuleLoadData(out module, img);
            testResult(res);
        }

        /// <summary> see CUDA doc; </summary>
        public static void ModuleGetFunction(out CUfunction func, CUmodule module, string name)
        {
            //IntPtr _name = Marshal.StringToHGlobalAnsi(name);
			IntPtr _name = Marshal.StringToHGlobalAnsi(name);
            CUresult res = my.cuModuleGetFunction(out func, module, _name);
            //Marshal.FreeHGlobal(_name);
			Marshal.FreeHGlobal(_name);
            testResult(res);
        }
        #endregion
        // ------------------------------------------
        #region Memory Management
        private delegate CUresult _cuMemAlloc(out CUdeviceptr dptr, uint bytesize);
        private delegate CUresult _cuMemHostAlloc(out IntPtr hptr, uint bytesize, uint flags);
        private delegate CUresult _cuMemHostGetDevicePointer(out CUdeviceptr dptr, IntPtr p, uint flags);
        private delegate CUresult _cuMemFreeHost(IntPtr p);
        private delegate CUresult _cuMemFree(CUdeviceptr dptr);
        private delegate CUresult _cuMemsetD8(CUdeviceptr dstDevice, byte uc, uint N);

        _cuMemAlloc cuMemAlloc;
        _cuMemHostAlloc cuMemHostAlloc;
        _cuMemHostGetDevicePointer cuMemHostGetDevicePointer;
        _cuMemFreeHost cuMemFreeHost;
        _cuMemFree cuMemFree;
        _cuMemsetD8 cuMemsetD8;

        /// <summary> see CUDA doc; </summary>
        public static void MemAlloc(out CUdeviceptr dptr, uint bytesize) {
            testResult(my.cuMemAlloc(out dptr, bytesize));
        }
        
        /// <summary> see CUDA doc; </summary>
        public static void MemsetD8(CUdeviceptr dstDevice, byte uc, uint N) {
            testResult(my.cuMemsetD8(dstDevice,uc,N));
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemHostAlloc(out IntPtr hptr, uint bytesize, params CUmem_host_alloc[] flags)
        {
            uint _flags = 0;
            foreach (CUmem_host_alloc x in flags)
            {
                _flags |= (uint)x;
            }
            testResult(my.cuMemHostAlloc(out hptr, bytesize, _flags));
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemHostGetDevicePointer(out CUdeviceptr dptr, IntPtr p, uint flags)
        {
            testResult(my.cuMemHostGetDevicePointer(out dptr, p, flags));
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemFreeHost(IntPtr p)
        {
            testResult(my.cuMemFreeHost(p));
        }

        /// <summary> see CUDA doc; </summary>
        static public void MemFree(CUdeviceptr dptr)
        {
            testResult(my.cuMemFree(dptr));
        }
        #endregion
        // ------------------------------------------
        #region Memory Transfer
        private delegate CUresult _cuMemcpyHtoD(CUdeviceptr dstDevice, IntPtr srcHost, uint byteCount);
        private delegate CUresult _cuMemcpyDtoH(IntPtr dstHost, CUdeviceptr srcDevice, uint byteCount);
        private delegate CUresult _cuMemcpyDtoD(CUdeviceptr dstDevice, CUdeviceptr srcDevice, uint byteCount);

        _cuMemcpyHtoD cuMemcpyHtoD;
        _cuMemcpyDtoH cuMemcpyDtoH;
        _cuMemcpyDtoD cuMemcpyDtoD;

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyHtoD(CUdeviceptr dstDevice, double[] srcHost, uint ByteCount) {
            unsafe {
                fixed (double* p = srcHost) {
                    testResult(my.cuMemcpyHtoD(dstDevice, (IntPtr)p, ByteCount));
                }
            }
        }
        
        /// <summary> see CUDA doc; </summary>
        public static void MemcpyHtoD(CUdeviceptr dstDevice, IntPtr srcHost, uint ByteCount) {
            testResult(my.cuMemcpyHtoD(dstDevice, srcHost, ByteCount));
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyHtoD(CUdeviceptr dstDevice, int[] srcHost, uint ByteCount)
        {
            unsafe
            {
                fixed (int* p = srcHost)
                {
                    testResult(my.cuMemcpyHtoD(dstDevice, (IntPtr)p, ByteCount));
                }
            }
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyHtoD(CUdeviceptr dstDevice, ushort[] srcHost, uint ByteCount)
        {
            unsafe
            {
                fixed (ushort* p = srcHost)
                {
                    testResult(my.cuMemcpyHtoD(dstDevice, (IntPtr)p, ByteCount));
                }
            }
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyDtoH(double[] dstHost, CUdeviceptr srcDevice, uint ByteCount)
        {
            unsafe
            {
                fixed (double* p = dstHost)
                {
                    testResult(my.cuMemcpyDtoH((IntPtr)p, srcDevice, ByteCount));
                }
            }
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyDtoH(int[] dstHost, CUdeviceptr srcDevice, uint ByteCount)
        {
            unsafe
            {
                fixed (int* p = dstHost)
                {
                    testResult(my.cuMemcpyDtoH((IntPtr)p, srcDevice, ByteCount));
                }
            }
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyDtoH(ushort[] dstHost, CUdeviceptr srcDevice, uint ByteCount)
        {
            unsafe
            {
                fixed (ushort* p = dstHost)
                {
                    testResult(my.cuMemcpyDtoH((IntPtr)p, srcDevice, ByteCount));
                }
            }
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyDtoD(CUdeviceptr dstDevice, CUdeviceptr srcDevice, uint ByteCount)
        {
            testResult(my.cuMemcpyDtoD(dstDevice, srcDevice, ByteCount));
        }
        #endregion
        // ------------------------------------------
        #region Kernel Parameter Management
        private delegate CUresult _cuParamSeti(CUfunction hfunc, int offset, uint value);
        private delegate CUresult _cuParamSetf(CUfunction hfunc, int offset, float value);
        private delegate CUresult _cuParamSetv(CUfunction hfunc, int offset, IntPtr ptr, uint numbytes);
        private delegate CUresult _cuParamSetSize(CUfunction hfunc, uint numbytes);

        _cuParamSeti cuParamSeti;
        _cuParamSetf cuParamSetf;
        _cuParamSetv cuParamSetv;
        _cuParamSetSize cuParamSetSize;

        /// <summary> see CUDA doc; </summary>
        static public void ParamSeti(CUfunction hfunc, int offset, int value)
        {
            testResult(my.cuParamSeti(hfunc, offset, (uint)value));
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetf(CUfunction hfunc, int offset, float value)
        {
            testResult(my.cuParamSetf(hfunc, offset, value));
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetd(CUfunction hfunc, int offset, double value)
        {
            unsafe
            {
                double* p = &value;
                testResult(my.cuParamSetv(hfunc, offset, (IntPtr)p, sizeof(double)));
            }
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetl(CUfunction hfunc, int offset, long value)
        {
            unsafe
            {
                long* p = &value;
                testResult(my.cuParamSetv(hfunc, offset, (IntPtr)p, sizeof(long)));
            }
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetp(CUfunction hfunc, int offset, CUdeviceptr ptr)
        {
            ParamSetl(hfunc, offset, (long)ptr.p);
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetSize(CUfunction hfunc, uint numbytes)
        {
            testResult(my.cuParamSetSize(hfunc, numbytes));
        }

        #endregion
        // ------------------------------------------
        #region Kernel Execution Management
        private delegate CUresult _cuFuncSetBlockShape(CUfunction hfunc, int x, int y, int z);
        private delegate CUresult _cuFuncSetSharedSize(CUfunction hfunc, uint bytes);
        private delegate CUresult _cuFuncGetAttribute(out int value, CUfunction_attribute attrib, CUfunction hfunc);
        private delegate CUresult _cuFuncSetCacheConfig(CUfunction hfunc, CUfunc_cache config);
        private delegate CUresult _cuLaunchGrid(CUfunction hfunc, int grid_width, int grid_height);
        private delegate CUresult _cuLaunchGridAsync(CUfunction hfunc, int grid_width, int grid_height, CUstream hStream);

        _cuFuncSetBlockShape cuFuncSetBlockShape;
        _cuFuncSetSharedSize cuFuncSetSharedSize;
        _cuFuncGetAttribute cuFuncGetAttribute;
        _cuFuncSetCacheConfig cuFuncSetCacheConfig;
        _cuLaunchGrid cuLaunchGrid;
        _cuLaunchGridAsync cuLaunchGridAsync;

        /// <summary> see CUDA doc; </summary>
        public static void FuncSetBlockShape(CUfunction hfunc, int x, int y, int z)
        {
            testResult(my.cuFuncSetBlockShape(hfunc, x, y, z));
        }

        /// <summary> see CUDA doc; </summary>
        static public void FuncSetSharedSize(CUfunction hfunc, uint bytes)
        {
            testResult(my.cuFuncSetSharedSize(hfunc, bytes));
        }

        /// <summary> see CUDA doc; </summary>
        static public void FuncGetAttribute(out int value, CUfunction_attribute attrib, CUfunction hfunc)
        {
            testResult(my.cuFuncGetAttribute(out value, attrib, hfunc));
        }

        /// <summary> 
        /// see CUDA doc; 
        /// </summary>
        public static void FuncSetCacheConfig(CUfunction hfunc, CUfunc_cache config)        {
            testResult(my.cuFuncSetCacheConfig(hfunc, config));
        }

        /// <summary> see CUDA doc; </summary>
        public static void LaunchGrid(CUfunction f, int grid_width, int grid_height)
        {
            testResult(my.cuLaunchGrid(f, grid_width, grid_height));
        }

        /// <summary> see CUDA doc; </summary>
        static public void LaunchGridAsync(CUfunction f, int grid_width, int grid_height, CUstream hStream)
        {
            testResult(my.cuLaunchGridAsync(f, grid_width, grid_height, hStream));
        }
        #endregion
        // ------------------------------------------
        #region Stream Management
        private delegate CUresult _cuStreamCreate(out CUstream hStream, uint Flags);
        private delegate CUresult _cuStreamQuery(CUstream hStream);
        private delegate CUresult _cuStreamSynchronize(CUstream hStream);
        private delegate CUresult _cuStreamDestroy(CUstream hStream);

        _cuStreamCreate cuStreamCreate;
        _cuStreamQuery cuStreamQuery;
        _cuStreamSynchronize cuStreamSynchronize;
        _cuStreamDestroy cuStreamDestroy;

        /// <summary> see CUDA doc; </summary>
        public static void StreamCreate(out CUstream hStream, uint Flags)
        {
            testResult(my.cuStreamCreate(out hStream, Flags));
        }

        /// <summary> see CUDA doc; </summary>
        public static CUresult StreamQuery(CUstream hStream)
        {
            return my.cuStreamQuery(hStream);
        }

        /// <summary> see CUDA doc; </summary>
        public static void StreamSynchronize(CUstream hStream)
        {
            testResult(my.cuStreamSynchronize(hStream));
        }

        /// <summary> see CUDA doc; </summary>
        public static void StreamDestroy(CUstream hStream)
        {
            testResult(my.cuStreamDestroy(hStream));
        }
        #endregion
#pragma warning restore        649
    }
}
