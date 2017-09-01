using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.IO;
using ilPSP.LinSolvers.monkey.CUDA;


namespace vectorAddDrvSharp {

    /*
    [StructLayout(LayoutKind.Explicit)]
    struct CUdevice {
        [FieldOffset(0)]
        int d;
    }

    [StructLayout(LayoutKind.Explicit)]
    struct CUcontext {
        [FieldOffset(0)]
        IntPtr p;
    }

    [StructLayout(LayoutKind.Explicit)]
    struct CUmodule {
        [FieldOffset(0)]
        IntPtr p;
    }

    [StructLayout(LayoutKind.Explicit)]
    struct CUfunction {
        [FieldOffset(0)]
        IntPtr p;
    }


    enum CUresult {

        CUDA_SUCCESS = 0,   ///< No errors
        CUDA_ERROR_INVALID_VALUE = 1,   ///< Invalid value
        CUDA_ERROR_OUT_OF_MEMORY = 2,   ///< Out of memory
        CUDA_ERROR_NOT_INITIALIZED = 3,   ///< Driver not initialized
        CUDA_ERROR_DEINITIALIZED = 4,   ///< Driver deinitialized

        CUDA_ERROR_NO_DEVICE = 100, ///< No CUDA-capable device available
        CUDA_ERROR_INVALID_DEVICE = 101, ///< Invalid device

        CUDA_ERROR_INVALID_IMAGE = 200, ///< Invalid kernel image
        CUDA_ERROR_INVALID_CONTEXT = 201, ///< Invalid context
        CUDA_ERROR_CONTEXT_ALREADY_CURRENT = 202, ///< Context already current
        CUDA_ERROR_MAP_FAILED = 205, ///< Map failed
        CUDA_ERROR_UNMAP_FAILED = 206, ///< Unmap failed
        CUDA_ERROR_ARRAY_IS_MAPPED = 207, ///< Array is mapped
        CUDA_ERROR_ALREADY_MAPPED = 208, ///< Already mapped
        CUDA_ERROR_NO_BINARY_FOR_GPU = 209, ///< No binary for GPU
        CUDA_ERROR_ALREADY_ACQUIRED = 210, ///< Already acquired
        CUDA_ERROR_NOT_MAPPED = 211, ///< Not mapped
        CUDA_ERROR_NOT_MAPPED_AS_ARRAY = 212, ///< Mapped resource not available for access as an array
        CUDA_ERROR_NOT_MAPPED_AS_POINTER = 213, ///< Mapped resource not available for access as a pointer
        CUDA_ERROR_ECC_UNCORRECTABLE = 214, ///< Uncorrectable ECC error detected
        CUDA_ERROR_UNSUPPORTED_LIMIT = 215, ///< CUlimit not supported by device

        CUDA_ERROR_INVALID_SOURCE = 300, ///< Invalid source
        CUDA_ERROR_FILE_NOT_FOUND = 301, ///< File not found
        CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND = 302, ///< Link to a shared object failed to resolve
        CUDA_ERROR_SHARED_OBJECT_INIT_FAILED = 303, ///< Shared object initialization failed

        CUDA_ERROR_INVALID_HANDLE = 400, ///< Invalid handle

        CUDA_ERROR_NOT_FOUND = 500, ///< Not found

        CUDA_ERROR_NOT_READY = 600, ///< CUDA not ready

        CUDA_ERROR_LAUNCH_FAILED = 700, ///< Launch failed
        CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES = 701, ///< Launch exceeded resources
        CUDA_ERROR_LAUNCH_TIMEOUT = 702, ///< Launch exceeded timeout
        CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING = 703, ///< Launch with incompatible texturing

        CUDA_ERROR_POINTER_IS_64BIT = 800, ///< Attempted to retrieve 64-bit pointer via 32-bit API function
        CUDA_ERROR_SIZE_IS_64BIT = 801, ///< Attempted to retrieve 64-bit size via 32-bit API function

        CUDA_ERROR_UNKNOWN = 999  ///< Unknown error
    }

    [StructLayout(LayoutKind.Explicit)]
    struct CUdeviceptr {
        [FieldOffset(0)]
        public uint p;
    }

    class cu {
        [DllImport("nvcuda")]
        public static extern CUresult cuInit(uint Flags);

        [DllImport("nvcuda")]
        public static extern CUresult cuDeviceGet(out CUdevice device, int ordinal);

        [DllImport("nvcuda")]
        public static extern CUresult cuCtxCreate(out CUcontext pctx, uint flags, CUdevice dev);

        [DllImport("nvcuda")]
        unsafe private static extern CUresult cuDeviceGetName(char* name, int len, CUdevice dev);


        public static CUresult cuDeviceGetName(out string name, CUdevice dev) {
            int maxlen = 1024;
            IntPtr mem = Marshal.AllocHGlobal(maxlen);
            CUresult res;
            unsafe {
                res = cuDeviceGetName((char*)mem, maxlen, dev);
            }
            name = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            return res;
        }

        [DllImport("nvcuda")]
        private static extern CUresult cuModuleLoad(out CUmodule module, IntPtr fname);

        public static CUresult cuModuleLoad(out CUmodule module, string fname) {
            IntPtr _fname = Marshal.StringToHGlobalAnsi(fname);
            CUresult res = cuModuleLoad(out module, _fname);
            Marshal.FreeHGlobal(_fname);
            return res;
        }

        [DllImport("nvcuda")]
        private static extern CUresult cuModuleGetFunction(out CUfunction hfunc, CUmodule hmod, IntPtr name);

        public static CUresult cuModuleGetFunction(out CUfunction func, CUmodule module, string name) {
            IntPtr _name = Marshal.StringToHGlobalAnsi(name);
            CUresult res = cuModuleGetFunction(out func, module, _name);
            Marshal.FreeHGlobal(_name);
            return res;
        }

        [DllImport("nvcuda")]
        public static extern CUresult cuMemAlloc(out CUdeviceptr dptr, uint bytesize);
       
        //[DllImport("nvcuda")]
        //unsafe public static extern CUresult cuMemAlloc(ulong*  dptr, uint bytesize);

        [DllImport("nvcuda")]
        unsafe public static extern CUresult cuMemcpyHtoD(CUdeviceptr dstDevice, void* srcHost, uint ByteCount);

        [DllImport("nvcuda")]
        unsafe public static extern CUresult cuMemcpyDtoH(void* dstHost, CUdeviceptr srcDevice, uint ByteCount);

        [DllImport("nvcuda")]
        public static extern CUresult cuMemcpyHtoD(CUdeviceptr dstDevice, double[] srcHost, uint ByteCount);

        [DllImport("nvcuda")]
        public static extern CUresult cuMemcpyDtoH(double[] dstHost, CUdeviceptr srcDevice, uint ByteCount);

        [DllImport("nvcuda")]
        public static extern CUresult cuMemcpyHtoD(CUdeviceptr dstDevice, int[] srcHost, uint ByteCount);

        [DllImport("nvcuda")]
        public static extern CUresult cuMemcpyDtoH(int[] dstHost, CUdeviceptr srcDevice, uint ByteCount);

        [DllImport("nvcuda")]
        public static extern CUresult cuParamSeti(CUfunction hfunc, int offset, uint value);

        [DllImport("nvcuda")]
        public static extern CUresult cuParamSetf(CUfunction hfunc, int offset, float value);

        [DllImport("nvcuda")]
        public static extern CUresult cuParamSetv(CUfunction hfunc, int offset, ref long ptr, uint numbytes);

        [DllImport("nvcuda")]
        public static extern CUresult cuParamSetv(CUfunction hfunc, int offset, ref double ptr, uint numbytes);

        [DllImport("nvcuda")]
        public static extern CUresult cuParamSetSize(CUfunction hfunc, uint numbytes);

        [DllImport("nvcuda")]
        public static extern CUresult cuFuncSetBlockShape(CUfunction hfunc, int x, int y, int z);

        [DllImport("nvcuda")]
        public static extern CUresult cuFuncSetSharedSize(CUfunction hfunc, uint bytes);

        [DllImport("nvcuda")]
        public static extern CUresult cuLaunchGrid(CUfunction f, int grid_width, int grid_height);

        [DllImport("nvcuda")]
        public static extern CUresult cuCtxSynchronize();

        [DllImport("nvcuda")]
        public static extern CUresult cuMemFree(CUdeviceptr dptr);
    }*/


    class cu2 : MPI.Wrappers.Utils.DynLibLoader {

        cu2() : base(
                new string[] { "nvcuda.dll", "libcuda.so" },
                new GetNameMangling[] { MPI.Wrappers.Utils.DynLibLoader.Identity, MPI.Wrappers.Utils.DynLibLoader.Identity },
                new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                new int[] { -1, -1 })
        {
        
        }

        private static cu2 my = new cu2();

        private static void TestResult(CUresult r)
        {
            if (r != CUresult.CUDA_SUCCESS)
                throw new ApplicationException("CUDA error: " + r.ToString());
        }

        // ------------------------------------------
        #region Initialization
        private delegate CUresult _cuInit(uint Flags);

#pragma warning disable        649
        _cuInit cuInit;
#pragma warning restore        649

        /// <summary> see CUDA doc; </summary>
        public static void Init(uint Flags)
        {
            TestResult(my.cuInit(Flags));
        }
        #endregion
        // ------------------------------------------
        #region Device Management
        private delegate CUresult _cuDeviceGetCount(out int count);
        private delegate CUresult _cuDeviceGet(out CUdevice device, int ordinal);
        private delegate CUresult _cuDeviceGetName(IntPtr name, int len, CUdevice device);
        private delegate CUresult _cuDeviceComputeCapability(out int major, out int minor, CUdevice dev);
        private delegate CUresult _cuDeviceTotalMem(out uint bytes, CUdevice dev);

#pragma warning disable        649
        _cuDeviceGetCount cuDeviceGetCount;
        _cuDeviceGet cuDeviceGet;
        _cuDeviceGetName cuDeviceGetName;
        _cuDeviceComputeCapability cuDeviceComputeCapability;
        _cuDeviceTotalMem cuDeviceTotalMem;
#pragma warning restore        649

        /// <summary> see CUDA doc; </summary>
        public static void DeviceGetCount(out int count)
        {
            TestResult(my.cuDeviceGetCount(out count));
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceGet(out CUdevice device, int ordinal)
        {
            TestResult(my.cuDeviceGet(out device, ordinal));
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceGetName(out string name, CUdevice dev)
        {
            int maxlen = 1024;
            IntPtr mem = Marshal.AllocHGlobal(maxlen);
            CUresult res;
            unsafe
            {
                res = my.cuDeviceGetName(mem, maxlen, dev);
            }
            name = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            TestResult(res);
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceComputeCapability(out int major, out int minor, CUdevice device)
        {
            TestResult(my.cuDeviceComputeCapability(out major, out minor, device));
        }

        /// <summary> see CUDA doc; </summary>
        public static void DeviceTotalMem(out uint bytes, CUdevice dev)
        {
            TestResult(my.cuDeviceTotalMem(out bytes, dev));
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
            TestResult(my.cuCtxCreate(out pctx, _flags, dev));
        }

        /// <summary> see CUDA doc; </summary>
        public static void CtxDestroy(CUcontext ctx)
        {
            TestResult(my.cuCtxDestroy(ctx));
        }

        /// <summary> see CUDA doc; </summary>
        public static void CtxSynchronize()
        {
            TestResult(my.cuCtxSynchronize());
        }
        #endregion
        // ------------------------------------------
        #region Module Management
        private delegate CUresult _cuModuleLoad(out CUmodule module, IntPtr fname);
        private delegate CUresult _cuModuleLoadData(out CUmodule module, IntPtr img);
        private delegate CUresult _cuModuleGetFunction(out CUfunction hfunc, CUmodule hmodule, IntPtr name);

        _cuModuleLoad cuModuleLoad;
        _cuModuleLoadData cuModuleLoadData;
        _cuModuleGetFunction cuModuleGetFunction;

        /// <summary> see CUDA doc; </summary>
        public static void ModuleLoad(out CUmodule module, string fname)
        {
            IntPtr _fname = Marshal.StringToHGlobalAnsi(fname);
            CUresult res = my.cuModuleLoad(out module, _fname);
            Marshal.FreeHGlobal(_fname);
            TestResult(res);
        }

        /// <summary> see CUDA doc; </summary>
        public static void ModuleLoadData(out CUmodule module, IntPtr img)
        {
            CUresult res = my.cuModuleLoadData(out module, img);
            TestResult(res);
        }

        /// <summary> see CUDA doc; </summary>
        public static void ModuleGetFunction(out CUfunction func, CUmodule module, string name)
        {
            IntPtr _name = Marshal.StringToHGlobalAnsi(name);
            CUresult res = my.cuModuleGetFunction(out func, module, _name);
            Marshal.FreeHGlobal(_name);
            TestResult(res);
        }
        #endregion
        // ------------------------------------------
        #region Memory Management
        private delegate CUresult _cuMemAlloc(out CUdeviceptr dptr, uint bytesize);
        private delegate CUresult _cuMemHostAlloc(out IntPtr hptr, uint bytesize, uint flags);
        private delegate CUresult _cuMemHostGetDevicePointer(out CUdeviceptr dptr, IntPtr p, uint flags);
        private delegate CUresult _cuMemFreeHost(IntPtr p);
        private delegate CUresult _cuMemFree(CUdeviceptr dptr);

        _cuMemAlloc cuMemAlloc;
        _cuMemHostAlloc cuMemHostAlloc;
        _cuMemHostGetDevicePointer cuMemHostGetDevicePointer;
        _cuMemFreeHost cuMemFreeHost;
        _cuMemFree cuMemFree;

        /// <summary> see CUDA doc; </summary>
        public static void MemAlloc(out CUdeviceptr dptr, uint bytesize)
        {
            TestResult(my.cuMemAlloc(out dptr, bytesize));
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemHostAlloc(out IntPtr hptr, uint bytesize, params CUmem_host_alloc[] flags)
        {
            uint _flags = 0;
            foreach (CUmem_host_alloc x in flags)
            {
                _flags |= (uint)x;
            }
            TestResult(my.cuMemHostAlloc(out hptr, bytesize, _flags));
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemHostGetDevicePointer(out CUdeviceptr dptr, IntPtr p, uint flags)
        {
            TestResult(my.cuMemHostGetDevicePointer(out dptr, p, flags));
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemFreeHost(IntPtr p)
        {
            TestResult(my.cuMemFreeHost(p));
        }

        /// <summary> see CUDA doc; </summary>
        static public void MemFree(CUdeviceptr dptr)
        {
            TestResult(my.cuMemFree(dptr));
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
        public static void MemcpyHtoD(CUdeviceptr dstDevice, double[] srcHost, uint ByteCount)
        {
            unsafe
            {
                fixed (double* p = srcHost)
                {
                    TestResult(my.cuMemcpyHtoD(dstDevice, (IntPtr)p, ByteCount));
                }
            }
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyHtoD(CUdeviceptr dstDevice, int[] srcHost, uint ByteCount)
        {
            unsafe
            {
                fixed (int* p = srcHost)
                {
                    TestResult(my.cuMemcpyHtoD(dstDevice, (IntPtr)p, ByteCount));
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
                    TestResult(my.cuMemcpyHtoD(dstDevice, (IntPtr)p, ByteCount));
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
                    TestResult(my.cuMemcpyDtoH((IntPtr)p, srcDevice, ByteCount));
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
                    TestResult(my.cuMemcpyDtoH((IntPtr)p, srcDevice, ByteCount));
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
                    TestResult(my.cuMemcpyDtoH((IntPtr)p, srcDevice, ByteCount));
                }
            }
        }

        /// <summary> see CUDA doc; </summary>
        public static void MemcpyDtoD(CUdeviceptr dstDevice, CUdeviceptr srcDevice, uint ByteCount)
        {
            TestResult(my.cuMemcpyDtoD(dstDevice, srcDevice, ByteCount));
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
        static public void ParamSeti(CUfunction hfunc, int offset, uint value)
        {
            TestResult(my.cuParamSeti(hfunc, offset, value));
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetf(CUfunction hfunc, int offset, float value)
        {
            TestResult(my.cuParamSetf(hfunc, offset, value));
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetd(CUfunction hfunc, int offset, double value)
        {
            unsafe
            {
                double* p = &value;
                TestResult(my.cuParamSetv(hfunc, offset, (IntPtr)p, sizeof(double)));
            }
        }

        /// <summary> see CUDA doc; </summary>
        static public void ParamSetl(CUfunction hfunc, int offset, long value)
        {
            unsafe
            {
                long* p = &value;
                TestResult(my.cuParamSetv(hfunc, offset, (IntPtr)p, sizeof(double)));
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
            TestResult(my.cuParamSetSize(hfunc, numbytes));
        }

        #endregion
        // ------------------------------------------
        #region Kernel Execution Management
        private delegate CUresult _cuFuncSetBlockShape(CUfunction hfunc, int x, int y, int z);
        private delegate CUresult _cuFuncSetSharedSize(CUfunction hfunc, uint bytes);
        //private delegate CUresult _cuFuncSetCacheConfig(CUfunction hfunc, CUfunc_cache config);
        private delegate CUresult _cuLaunchGrid(CUfunction hfunc, int grid_width, int grid_height);
        private delegate CUresult _cuLaunchGridAsync(CUfunction hfunc, int grid_width, int grid_height, CUstream hStream);

        _cuFuncSetBlockShape cuFuncSetBlockShape;
        _cuFuncSetSharedSize cuFuncSetSharedSize;
        //_cuFuncSetCacheConfig cuFuncSetCacheConfig;
        _cuLaunchGrid cuLaunchGrid;
        _cuLaunchGridAsync cuLaunchGridAsync;

        /// <summary> see CUDA doc; </summary>
        public static void FuncSetBlockShape(CUfunction hfunc, int x, int y, int z)
        {
            TestResult(my.cuFuncSetBlockShape(hfunc, x, y, z));
        }

        /// <summary> see CUDA doc; </summary>
        static public void FuncSetSharedSize(CUfunction hfunc, uint bytes)
        {
            TestResult(my.cuFuncSetSharedSize(hfunc, bytes));
        }

        /// <summary> see CUDA doc; </summary>
        //public static void FuncSetCacheConfig(CUfunction hfunc, CUfunc_cache config)
        //{
        //    TestResult(my.cuFuncSetCacheConfig(hfunc, config));
        //}

        /// <summary> see CUDA doc; </summary>
        public static void LaunchGrid(CUfunction f, int grid_width, int grid_height)
        {
            TestResult(my.cuLaunchGrid(f, grid_width, grid_height));
        }

        /// <summary> see CUDA doc; </summary>
        static public void LaunchGridAsync(CUfunction f, int grid_width, int grid_height, CUstream hStream)
        {
            TestResult(my.cuLaunchGridAsync(f, grid_width, grid_height, hStream));
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
            TestResult(my.cuStreamCreate(out hStream, Flags));
        }

        /// <summary> see CUDA doc; </summary>
        public static CUresult StreamQuery(CUstream hStream)
        {
            return my.cuStreamQuery(hStream);
        }

        /// <summary> see CUDA doc; </summary>
        public static void StreamSynchronize(CUstream hStream)
        {
            TestResult(my.cuStreamSynchronize(hStream));
        }

        /// <summary> see CUDA doc; </summary>
        public static void StreamDestroy(CUstream hStream)
        {
            TestResult(my.cuStreamDestroy(hStream));
        }
        #endregion
    }

    class vectorAdd {

        static int Main(string[] args) {

            //Console.WriteLine("CPU pointer size: " + IntPtr.Size);

            cu2.Init(0);
            int NoOfDev;
            cu2.DeviceGetCount(out NoOfDev);
            Console.WriteLine("No of devices: " + NoOfDev);
            for (int i = 0; i < NoOfDev; i++)
            {
                CUdevice sdevice;
                String sname;
                cu2.DeviceGet(out sdevice, i);
                cu2.DeviceGetName(out sname, sdevice);
                Console.WriteLine("Device " + i + ": " + sname);
            }
            
            // Select device by number and create context
            // The 480 GTX should be 0 and the Quadro NVS then 1
            CUdevice device;
            cu2.DeviceGet(out device, 0);

            CUcontext context;
            cu2.CtxCreate(out context, device);

            // Load kernel function from pre-compiled PTX file
            CUmodule module;
            string path;
            {
                System.Reflection.Assembly a = System.Reflection.Assembly.GetEntryAssembly();
                DirectoryInfo di = new DirectoryInfo(Path.GetDirectoryName(a.Location));
                path = Path.Combine(di.ToString(), "CudaMatrixKernelDP.ptx");
            }
            cu2.ModuleLoad(out module, path);

            CUfunction func;
            cu2.ModuleGetFunction(out func, module, "sparseMultiply");

            // Prepare parameters, create data, allocate and copy memory
            // Number of rows
            uint rows = 256000;
            // Number of matrix entries
            uint values = rows*4;
            // Size of one CUDA block
            int blocksize = 128;
            // Number of blocks
            int blockcount = (int)Math.Ceiling((float)rows / blocksize);
            // Offset for parameters
            int offset = 0;
            double alpha = 2.0;
            double beta = 0.0;

            // Device memory allocation
            CUdeviceptr d_val, d_colIdx, d_rowStart, d_result, d_x;
            cu2.MemAlloc(out d_val, values * sizeof(double));
            cu2.MemAlloc(out d_colIdx, values * sizeof(int));
            cu2.MemAlloc(out d_rowStart, (rows + 1) * sizeof(int));
            cu2.MemAlloc(out d_result, rows * sizeof(double));
            cu2.MemAlloc(out d_x, rows * sizeof(double));

            //ulong test = 0xFFFFFFFFFFFFFFFF;
            //unsafe {
            //    cu2.MemAlloc(&test, 122));
            //}


            // Host memory allocation
            double[] h_val = new double[values];
            int[] h_colIdx = new int[values];
            int[] h_rowStart = new int[rows + 1];
            double[] h_result = new double[rows];
            double[] h_x = new double[rows];

            // Create sparse matrix with random values in random columns
            Random rnd = new Random(DateTime.Now.Millisecond);
            for (uint i = 0; i < values; i++) {
                h_val[i] = rnd.NextDouble();
                h_colIdx[i] = (int)Math.Abs((int)rnd.Next()) % (int)rows;
            }
            // For simplicity put 4 values in one row
            // Initialize x with random values and result with 0
            for (uint i = 0; i < rows; i++) {
                h_rowStart[i] = (int)(i * 4);
                h_x[i] = rnd.NextDouble();
                h_result[i] = 0.0;
            }
            // Last value / end of last row
            h_rowStart[rows] = (int)values;

            // Copy the 4 necessary arrays to device memory
            cu2.MemcpyHtoD(d_val, h_val, values * sizeof(double));
            cu2.MemcpyHtoD(d_colIdx, h_colIdx, values * sizeof(int));
            cu2.MemcpyHtoD(d_rowStart, h_rowStart, (rows + 1) * sizeof(int));
            cu2.MemcpyHtoD(d_x, h_x, rows * sizeof(double));

            // Set kernel parameters, launch kernel and copy results back
            // Device pointers are casted to void* before given to cuParamSetv, which takes a pointer to arbitrary data
            // In each step, offset is the accumulated size of all previous parameters
            
            cu2.ParamSetp(func, offset, d_val);
            offset += sizeof(long);
            
            cu2.ParamSetp(func, offset, d_colIdx);
            offset += sizeof(long);
            
            cu2.ParamSetp(func, offset, d_rowStart);
            offset += sizeof(long);
            
            cu2.ParamSetp(func, offset, d_result);
            offset += sizeof(long);
            
            cu2.ParamSetp(func, offset, d_x);
            offset += sizeof(long);
            
            cu2.ParamSetd(func, offset, alpha);
            offset += sizeof(double);
            
            cu2.ParamSetd(func, offset, beta);
            offset += sizeof(double);
            
            cu2.ParamSeti(func, offset, rows);
            offset += sizeof(uint);
            
            // Set overall parameter size
            cu2.ParamSetSize(func, (uint)offset);
            // Set block size
            cu2.FuncSetBlockShape(func, blocksize, 1, 1);
            // Set size of dynamically allocated shared memory ("extern __shared__" in kernel code)
            cu2.FuncSetSharedSize(func, (uint)((blocksize + 1) * sizeof(int)));
            // Launch kernel with number of blocks
            cu2.LaunchGrid(func, blockcount, 1);
            // Syncronize (not needed, because launch isn't async)
            cu2.CtxSynchronize();

            // Copy results back
            cu2.MemcpyDtoH(h_result, d_result, rows * sizeof(double));
            
            ////// Free memory and validate results
            //// Free device memory
            cu2.MemFree(d_val);
            cu2.MemFree(d_colIdx);
            cu2.MemFree(d_rowStart);
            cu2.MemFree(d_result);
            cu2.MemFree(d_x);

            // Calculate reference data on host
            int pos = 0;
            double rowacc;
            double[] x_result = new double[rows];
            for (uint i = 0; i < rows; i++) {
                rowacc = 0.0;
                for (int j = h_rowStart[i]; j < h_rowStart[i + 1]; j++) {
                    rowacc += h_val[pos] * h_x[h_colIdx[pos]];
                    pos++;
                }
                x_result[i] = rowacc * alpha;
            }

            // Compare host and device data
            int flag = 0;
            for (uint i = 0; i < rows; i++) {
                if (Math.Abs(h_result[i] - x_result[i]) > 0.000001) {
                    flag = 1;
                    //Console.WriteLine("Error in row " + i + ">> CPU: " + x_result[i] + " , GPU: " + h_result[i]);
                }
            }
            if (flag == 0) {
                Console.WriteLine("No Errors.");
            }

            return 0;
        }

    }
}