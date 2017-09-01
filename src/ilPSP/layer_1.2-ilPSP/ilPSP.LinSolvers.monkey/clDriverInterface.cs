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

namespace ilPSP.LinSolvers.monkey.CL
{
    // ------------------------------------------
    #region OpenCL data types
    /// <summary>
    /// OpenCL context identifier
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    public struct cl_context
    {
        [FieldOffset(0)]
        internal IntPtr p;
    }

    /// <summary>
    /// OpenCL device identifier
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    public struct cl_device_id
    {
        [FieldOffset(0)]
        internal IntPtr p;
    }

    /// <summary>
    /// OpenCL platform identifier
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    public struct cl_platform_id
    {
        [FieldOffset(0)]
        internal IntPtr p;

        ///// <summary>
        ///// Tests if two platforms are euqal.
        ///// </summary>
        ///// <param name="a">Platform a</param>
        ///// <param name="b">Platform b</param>
        ///// <returns>If they are euqal</returns>
        //public static bool operator ==(cl_platform_id a, cl_platform_id b)
        //{
        //    if (a == null && b == null)
        //    {
        //        return true;
        //    }
        //    if (a == null || b == null)
        //    {
        //        return false;
        //    }
        //    return a.p == b.p;
        //}

        ///// <summary>
        ///// Tests if two platforms are not euqal.
        ///// </summary>
        ///// <param name="a">Platform a</param>
        ///// <param name="b">Platform b</param>
        ///// <returns>If they are not euqal</returns>
        //public static bool operator !=(cl_platform_id a, cl_platform_id b)
        //{
        //    return !(a == b);
        //}
    }

    /// <summary>
    /// OpenCL command queue
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    public struct cl_command_queue
    {
        [FieldOffset(0)]
        internal IntPtr p;
    }

    /// <summary>
    /// OpenCL device memory buffer
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    public struct cl_mem
    {
        [FieldOffset(0)]
        internal IntPtr p;
    }

    /// <summary>
    /// OpenCL program
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    public struct cl_program
    {
        [FieldOffset(0)]
        internal IntPtr p;
    }

    /// <summary>
    /// OpenCL kernel
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    public struct cl_kernel
    {
        [FieldOffset(0)]
        internal IntPtr p;
    }

    /// <summary>
    /// OpenCL event identifier
    /// </summary>
    [StructLayout(LayoutKind.Explicit)]
    struct cl_event_id
    {
        [FieldOffset(0)]
        internal IntPtr p;
    }

    /// <summary>
    /// OpenCL event
    /// </summary>
    /// Capsulates an event identifier
    public class cl_event
    {
        internal cl_event_id id;

        /// <summary>
        /// Constructor creates event from event identifier
        /// </summary>
        /// <param name="id">Event</param>
        internal cl_event(cl_event_id id)
        {
            this.id = id;
        }

        /// <summary>
        /// Destructor
        /// </summary>
        ~cl_event()
        {
            cl.ReleaseEvent(this);
        }

        internal static cl_event_id[] convertList(cl_event[] list)
        {
            if (list == null)
            {
                return null;
            }

            cl_event_id[] ret = new cl_event_id[list.Length];
            for (int i = 0; i < list.Length; i++)
            {
                ret[i] = list[i].id;
            }
            return ret;
        }

        /// <summary>
        /// Return the command of this event.
        /// </summary>
        /// <returns>The command type</returns>
        public cl_command_type getCommandType()
        {
            return (cl_command_type)cl.GetEventInfoUInt(this, cl_event_info.CL_EVENT_COMMAND_TYPE);
        }

        /// <summary>
        /// Return the status of this event.
        /// </summary>
        /// <returns>The command status</returns>
        public cl_status getCommandStatus()
        {
            return (cl_status)cl.GetEventInfoUInt(this, cl_event_info.CL_EVENT_COMMAND_EXECUTION_STATUS);
        }

        /// <summary>
        /// Return the command queue of this event.
        /// </summary>
        /// <returns>The command-queue</returns>
        public cl_command_queue getCommandQueue()
        {
            cl_command_queue cq;
            cq.p = cl.GetEventInfoIntPtr(this, cl_event_info.CL_EVENT_COMMAND_QUEUE);
            return cq;
        }

        /// <summary>
        /// Wait for this event to complete.
        /// </summary>
        public void waitFor()
        {
            cl.WaitForEvent(this);
        }

#if OPENCL_1_1
        /// <summary>
        /// Set Callback function.
        /// </summary>
        /// <param name="command_exec_callback_type">The command execution status for which the callback is registered.</param>
        /// <param name="pfn_event_notify">Callback, see <see cref="cl_set_event_callback"/></param>
        /// <param name="user_data">A pointer to user supplied data.</param>
        public void setCallback(cl_command_execution_status command_exec_callback_type, cl_set_event_callback pfn_event_notify, IntPtr user_data = default(IntPtr))
        {
            cl.SetEventCallback(this, command_exec_callback_type, pfn_event_notify, user_data);
        }
#endif
    }
    #endregion
    // ------------------------------------------
    #region OpenCL enumerations
    /// <summary>
    /// OpenCL device type
    /// </summary>
    public enum cl_device_type : ulong
    {
        /// <summary>Default device</summary>
        CL_DEVICE_TYPE_DEFAULT                       = (1 << 0),
        /// <summary>CPU</summary>
        CL_DEVICE_TYPE_CPU                           = (1 << 1),
        /// <summary>GPU</summary>
        CL_DEVICE_TYPE_GPU                           = (1 << 2),
        /// <summary>Accelerator (e.g. Cell BE)</summary>
        CL_DEVICE_TYPE_ACCELERATOR                   = (1 << 3),
        /// <summary>All devices</summary>
        CL_DEVICE_TYPE_ALL                           = 0xFFFFFFFF
    }

    /// <summary>
    /// OpenCL device properties
    /// </summary>
    public enum cl_device_info : uint
    {
        /// <summary>The OpenCL device type.</summary>
        CL_DEVICE_TYPE                               = 0x1000,
        /// <summary>A unique device vendor identifier. An example of a unique device identifier could be the PCIe ID.</summary>
        CL_DEVICE_VENDOR_ID                          = 0x1001,
        /// <summary>The number of parallel compute cores on the OpenCL device. The minimum value is 1.</summary>
        CL_DEVICE_MAX_COMPUTE_UNITS                  = 0x1002,
        /// <summary>Maximum dimensions that specify the global and local work-item IDs used by the data parallel execution model. The minimum value is 3.</summary>
        CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS           = 0x1003,
        /// <summary>Maximum number of work-items in a work-group executing a kernel using the data parallel execution model. The minimum value is 1.</summary>
        CL_DEVICE_MAX_WORK_GROUP_SIZE                = 0x1004,
        /// <summary>Maximum number of work-items that can be specified in each dimension of the work-group to clEnqueueNDRangeKernel. Returns n size_t entries, where n is the value returned by the query for CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS. The minimum value is (1, 1, 1).</summary>
        CL_DEVICE_MAX_WORK_ITEM_SIZES                = 0x1005,
        /// <summary>Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR        = 0x1006,
        /// <summary>Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT       = 0x1007,
        /// <summary>Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT         = 0x1008,
        /// <summary>Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG        = 0x1009,
        /// <summary>Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT       = 0x100A,
        /// <summary>Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE      = 0x100B,
        /// <summary>Maximum configured clock frequency of the device in MHz.</summary>
        CL_DEVICE_MAX_CLOCK_FREQUENCY                = 0x100C,
        /// <summary>The default compute device address space size specified as an unsigned integer value in bits. Currently supported values are 32 or 64 bits.</summary>
        CL_DEVICE_ADDRESS_BITS                       = 0x100D,
        /// <summary>Max number of simultaneous image objects that can be read by a kernel. The minimum value is 128 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_MAX_READ_IMAGE_ARGS                = 0x100E,
        /// <summary>Max number of simultaneous image objects that can be written to by a kernel. The minimum value is 8 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_MAX_WRITE_IMAGE_ARGS               = 0x100F,
        /// <summary>Max size of memory object allocation in bytes. The minimum value is max (1/4th of CL_DEVICE_GLOBAL_MEM_SIZE, 128*1024*1024).</summary>
        CL_DEVICE_MAX_MEM_ALLOC_SIZE                 = 0x1010,
        /// <summary>Max width of 2D image in pixels. The minimum value is 8192 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_IMAGE2D_MAX_WIDTH                  = 0x1011,
        /// <summary>Max height of 2D image in pixels. The minimum value is 8192 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_IMAGE2D_MAX_HEIGHT                 = 0x1012,
        /// <summary>Max width of 3D image in pixels. The minimum value is 2048 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_IMAGE3D_MAX_WIDTH                  = 0x1013,
        /// <summary>Max height of 3D image in pixels. The minimum value is 2048 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_IMAGE3D_MAX_HEIGHT                 = 0x1014,
        /// <summary>Max depth of 3D image in pixels. The minimum value is 2048 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_IMAGE3D_MAX_DEPTH                  = 0x1015,
        /// <summary>Is CL_TRUE if images are supported by the OpenCL device and CL_FALSE otherwise.</summary>
        CL_DEVICE_IMAGE_SUPPORT                      = 0x1016,
        /// <summary>Max size in bytes of the arguments that can be passed to a kernel. The minimum value is 1024. For this minimum value, only a maximum of 128 arguments can be passed to a kernel.</summary>
        CL_DEVICE_MAX_PARAMETER_SIZE                 = 0x1017,
        /// <summary>Maximum number of samplers that can be used in a kernel. The minimum value is 16 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        CL_DEVICE_MAX_SAMPLERS                       = 0x1018,
        /// <summary>Describes the alignment in bits of the base address of any allocated memory object.</summary>
        CL_DEVICE_MEM_BASE_ADDR_ALIGN                = 0x1019,
        /// <summary>The smallest alignment in bytes which can be used for any data type.</summary>
        CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE           = 0x101A,
        /// <summary>Describes single precision floating-point capability of the device.</summary>
        CL_DEVICE_SINGLE_FP_CONFIG                   = 0x101B,
        /// <summary>Type of global memory cache supported. Valid values are: CL_NONE, CL_READ_ONLY_CACHE, and CL_READ_WRITE_CACHE.</summary>
        CL_DEVICE_GLOBAL_MEM_CACHE_TYPE              = 0x101C,
        /// <summary>Size of global memory cache line in bytes.</summary>
        CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE          = 0x101D,
        /// <summary>Size of global memory cache in bytes.</summary>
        CL_DEVICE_GLOBAL_MEM_CACHE_SIZE              = 0x101E,
        /// <summary>Size of global device memory in bytes.</summary>
        CL_DEVICE_GLOBAL_MEM_SIZE                    = 0x101F,
        /// <summary>Max size in bytes of a constant buffer allocation. The minimum value is 64 KB.</summary>
        CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE           = 0x1020,
        /// <summary>Max number of arguments declared with the __constant qualifier in a kernel. The minimum value is 8.</summary>
        CL_DEVICE_MAX_CONSTANT_ARGS                  = 0x1021,
        /// <summary>Type of local memory supported. This can be set to CL_LOCAL implying dedicated local memory storage such as SRAM, or CL_GLOBAL.</summary>
        CL_DEVICE_LOCAL_MEM_TYPE                     = 0x1022,
        /// <summary>Size of local memory arena in bytes. The minimum value is 32 KB.</summary>
        CL_DEVICE_LOCAL_MEM_SIZE                     = 0x1023,
        /// <summary>Is CL_TRUE if the device implements error correction for all accesses to compute device memory (global and constant). Is CL_FALSE if the device does not implement such error correction.</summary>
        CL_DEVICE_ERROR_CORRECTION_SUPPORT           = 0x1024,
        /// <summary>Describes the resolution of device timer. This is measured in nanoseconds.</summary>
        CL_DEVICE_PROFILING_TIMER_RESOLUTION         = 0x1025,
        /// <summary>Is CL_TRUE if the OpenCL device is a little endian device and CL_FALSE otherwise.</summary>
        CL_DEVICE_ENDIAN_LITTLE                      = 0x1026,
        /// <summary>Is CL_TRUE if the device is available and CL_FALSE if the device is not available.</summary>
        CL_DEVICE_AVAILABLE                          = 0x1027,
        /// <summary>Is CL_FALSE if the implementation does not have a compiler available to compile the program source. Is CL_TRUE if the compiler is available. This can be CL_FALSE for the embedded platform profile only.</summary>
        CL_DEVICE_COMPILER_AVAILABLE                 = 0x1028,
        /// <summary>Describes the execution capabilities of the device.</summary>
        CL_DEVICE_EXECUTION_CAPABILITIES             = 0x1029,
        /// <summary>Describes the command-queue properties supported by the device.</summary>
        CL_DEVICE_QUEUE_PROPERTIES                   = 0x102A,
        /// <summary>Device name string.</summary>
        CL_DEVICE_NAME                               = 0x102B,
        /// <summary>Vendor name string.</summary>
        CL_DEVICE_VENDOR                             = 0x102C,
        /// <summary>OpenCL software driver version string in the form major_number.minor_number.</summary>
        CL_DRIVER_VERSION                            = 0x102D,
        /// <summary>OpenCL profile string. Returns the profile name supported by the device.</summary>
        CL_DEVICE_PROFILE                            = 0x102E,
        /// <summary>OpenCL version string. Returns the OpenCL version supported by the device.</summary>
        CL_DEVICE_VERSION                            = 0x102F,
        /// <summary>Returns a space-separated list of extension names (the extension names themselves do not contain any spaces).</summary>
        CL_DEVICE_EXTENSIONS                         = 0x1030,
        /// <summary>The platform associated with this device.</summary>
        CL_DEVICE_PLATFORM                           = 0x1031,
        /// <summary>Describes the OPTIONAL double precision floating-point capability of the OpenCL device.</summary>
        CL_DEVICE_DOUBLE_FP_CONFIG                   = 0x1032,
        /// <summary>Describes the OPTIONAL half precision floating-point capability of the OpenCL device.</summary>
        CL_DEVICE_HALF_FP_CONFIG                     = 0x1033,
        /// <summary>Preferred native vector width size for built-in scalar types that can be put into vectors. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF        = 0x1034,
        /// <summary>Returns the native ISA vector width. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR           = 0x1036,
        /// <summary>Returns the native ISA vector width. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT          = 0x1037,
        /// <summary>Returns the native ISA vector width. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_NATIVE_VECTOR_WIDTH_INT            = 0x1038,
        /// <summary>Returns the native ISA vector width. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG           = 0x1039,
        /// <summary>Returns the native ISA vector width. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT          = 0x103A,
        /// <summary>Returns the native ISA vector width. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE         = 0x103B,
        /// <summary>Returns the native ISA vector width. The vector width is defined as the number of scalar elements that can be stored in the vector.</summary>
        CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF           = 0x103C,
#if OPENCL_1_1
        /// <summary>OpenCL C version string. Returns the highest OpenCL C version supported by the compiler for this device.</summary>
        CL_DEVICE_OPENCL_C_VERSION                   = 0x103D,
        /// <summary>Is CL_TRUE if the device and the host have a unified memory subsystem and is CL_FALSE otherwise.</summary>
        CL_DEVICE_HOST_UNIFIED_MEMORY                = 0x1035,
#endif
    }

    /// <summary>
    /// OpenCL device properties structure
    /// </summary>
    public struct cl_device_info_return
    {
        /// <summary>The default compute device address space size specified as an unsigned integer value in bits. Currently supported values are 32 or 64 bits.</summary>
        public uint address_bits;
        /// <summary>Is CL_TRUE if the device is available and CL_FALSE if the device is not available.</summary>
        public bool available;
        /// <summary>Is CL_FALSE if the implementation does not have a compiler available to compile the program source. Is CL_TRUE if the compiler is available. This can be CL_FALSE for the embedded platform profile only.</summary>
        public bool compiler_available;
        /// <summary>Describes the OPTIONAL double precision floating-point capability of the OpenCL device.</summary>
        public cl_device_fp_config double_fp_config;
        /// <summary>Is CL_TRUE if the OpenCL device is a little endian device and CL_FALSE otherwise.</summary>
        public bool endian_little;
        /// <summary>Is CL_TRUE if the device implements error correction for all accesses to compute device memory (global and constant). Is CL_FALSE if the device does not implement such error correction.</summary>
        public bool error_correction_support;
        /// <summary>Describes the execution capabilities of the device.</summary>
        public cl_device_exec_capabilities execution_capabilities;
        /// <summary>Returns a space-separated list of extension names (the extension names themselves do not contain any spaces).</summary>
        public string extensions;
        /// <summary>Size of global memory cache in bytes.</summary>
        public ulong global_mem_cache_size;
        /// <summary>Type of global memory cache supported. Valid values are: CL_NONE, CL_READ_ONLY_CACHE, and CL_READ_WRITE_CACHE.</summary>
        public cl_device_mem_cache_type global_mem_cache_type;
        /// <summary>Size of global memory cache line in bytes.</summary>
        public uint global_mem_cacheline_size;
        /// <summary>Size of global device memory in bytes.</summary>
        public ulong global_mem_size;
        // <summary>Describes the OPTIONAL half precision floating-point capability of the OpenCL device.</summary>
        //public cl_device_fp_config half_fp_config;
        /// <summary>Is CL_TRUE if images are supported by the OpenCL device and CL_FALSE otherwise.</summary>
        public bool image_support;
        /// <summary>Max height of 2D image in pixels. The minimum value is 8192 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint image2d_max_height;
        /// <summary>Max width of 2D image in pixels. The minimum value is 8192 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint image2d_max_width;
        /// <summary>Max depth of 3D image in pixels. The minimum value is 2048 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint image3d_max_depth;
        /// <summary>Max height of 3D image in pixels. The minimum value is 2048 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint image3d_max_height;
        /// <summary>Max width of 3D image in pixels. The minimum value is 2048 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint image3d_max_width;
        /// <summary>Size of local memory arena in bytes. The minimum value is 32 KB.</summary>
        public ulong local_mem_size;
        /// <summary>Type of local memory supported. This can be set to CL_LOCAL implying dedicated local memory storage such as SRAM, or CL_GLOBAL.</summary>
        public cl_device_local_mem_type local_mem_type;
        /// <summary>Maximum configured clock frequency of the device in MHz.</summary>
        public uint max_clock_frequency;
        /// <summary>The number of parallel compute cores on the OpenCL device. The minimum value is 1.</summary>
        public uint max_compute_units;
        /// <summary>Max number of arguments declared with the __constant qualifier in a kernel. The minimum value is 8.</summary>
        public uint max_constant_args;
        /// <summary>Max size in bytes of a constant buffer allocation. The minimum value is 64 KB.</summary>
        public ulong max_constant_buffer_size;
        /// <summary>Max size of memory object allocation in bytes. The minimum value is max (1/4th of CL_DEVICE_GLOBAL_MEM_SIZE, 128*1024*1024).</summary>
        public ulong max_mem_alloc_size;
        /// <summary>Max size in bytes of the arguments that can be passed to a kernel. The minimum value is 1024. For this minimum value, only a maximum of 128 arguments can be passed to a kernel.</summary>
        public uint max_parameter_size;
        /// <summary>Max number of simultaneous image objects that can be read by a kernel. The minimum value is 128 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint max_read_image_args;
        /// <summary>Maximum number of samplers that can be used in a kernel. The minimum value is 16 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint max_samplers;
        /// <summary>Maximum number of work-items in a work-group executing a kernel using the data parallel execution model. The minimum value is 1.</summary>
        public uint max_work_group_size;
        /// <summary>Maximum dimensions that specify the global and local work-item IDs used by the data parallel execution model. The minimum value is 3.</summary>
        public uint max_work_item_dimensions;
        /// <summary>Maximum number of work-items that can be specified in each dimension of the work-group to clEnqueueNDRangeKernel. Returns n size_t entries, where n is the value returned by the query for CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS. The minimum value is (1, 1, 1).</summary>
        public uint[] max_work_item_sizes;
        /// <summary>Max number of simultaneous image objects that can be written to by a kernel. The minimum value is 8 if CL_DEVICE_IMAGE_SUPPORT is CL_TRUE.</summary>
        public uint max_write_image_args;
        /// <summary>Describes the alignment in bits of the base address of any allocated memory object.</summary>
        public uint mem_base_addr_align;
        /// <summary>The smallest alignment in bytes which can be used for any data type.</summary>
        public uint min_data_type_align_size;
        /// <summary>Device name string.</summary>
        public string name;
        /// <summary>The platform associated with this device.</summary>
        public cl_platform_id platform;
        /// <summary>OpenCL profile string. Returns the profile name supported by the device.</summary>
        public string profile;
        /// <summary>Describes the resolution of device timer. This is measured in nanoseconds.</summary>
        public uint profiling_timer_resultion;
        /// <summary>Describes the command-queue properties supported by the device.</summary>
        public cl_command_queue_properties queue_properties;
        /// <summary>Describes single precision floating-point capability of the device.</summary>
        public cl_device_fp_config single_fp_config;
        /// <summary>The OpenCL device type.</summary>
        public cl_device_type type;
        /// <summary>Vendor name string.</summary>
        public string vendor;
        /// <summary>A unique device vendor identifier. An example of a unique device identifier could be the PCIe ID.</summary>
        public uint vendor_id;
        /// <summary>OpenCL version string. Returns the OpenCL version supported by the device.</summary>
        public string version;
#if OPENCL_1_1
        /// <summary>Is CL_TRUE if the device and the host have a unified memory subsystem and is CL_FALSE otherwise.</summary>
        public bool host_unified_memory;
        /// <summary>OpenCL C version string. Returns the highest OpenCL C version supported by the compiler for this device.</summary>
        public string opencl_c_version;
#endif
    }

    /// <summary>
    /// Floating-point capability
    /// </summary>
    public enum cl_device_fp_config : ulong
    {
        /// <summary>Denorms are supported.</summary>
        CL_FP_DENORM                                 = (1 << 0),
        /// <summary>CL_FP_INF_NAN - INF and NaNs are supported.</summary>
        CL_FP_INF_NAN                                = (1 << 1),
        /// <summary>Round to nearest even rounding mode supported.</summary>
        CL_FP_ROUND_TO_NEAREST                       = (1 << 2),
        /// <summary>Round to zero rounding mode supported.</summary>
        CL_FP_ROUND_TO_ZERO                          = (1 << 3),
        /// <summary>Round to +ve and -ve infinity rounding modes supported.</summary>
        CL_FP_ROUND_TO_INF                           = (1 << 4),
        /// <summary>IEEE754-2008 fused multiply-add is supported. </summary>
        CL_FP_FMA                                    = (1 << 5),
        /// <summary>Basic floating-point operations (such as addition, subtraction, multiplication) are implemented in software.</summary>
        CL_FP_SOFT_FLOAT                             = (1 << 6)
    }

    /// <summary>
    /// Local memory type
    /// </summary>
    public enum cl_device_local_mem_type : uint
    {
        /// <summary>Dedicated local memory storage such as SRAM</summary>
        CL_LOCAL                                     = 0x1,
        /// <summary>Inside global memory</summary>
        CL_GLOBAL                                    = 0x2
    }

    /// <summary>
    /// Global memory cache type
    /// </summary>
    public enum cl_device_mem_cache_type : uint
    {
        /// <summary>No Cache</summary>
        CL_NONE                                      = 0x0,
        /// <summary>Read-Only-Cache</summary>
        CL_READ_ONLY_CACHE                           = 0x1,
        /// <summary>Write-Only-Cache</summary>
        CL_READ_WRITE_CACHE                          = 0x2
    }

    /// <summary>
    /// Execution capability
    /// </summary>
    public enum cl_device_exec_capabilities : ulong
    {
        /// <summary>The OpenCL device can execute OpenCL kernels.</summary>
        CL_EXEC_KERNEL                               = (1 << 0),
        /// <summary>The OpenCL device can execute native kernels.</summary>
        CL_EXEC_NATIVE_KERNEL                        = (1 << 1)
    }

    /// <summary>
    /// OpenCL error codes
    /// </summary>
    /// Also used as OpenCL command execution status
    public enum cl_status : int
    {
        /// <summary>Command has been enqueued in the command-queue.</summary>
        CL_QUEUED                                    = 0x3,
        /// <summary>Enqueued command has been submitted by the host to the device associated with the command-queue.</summary>
        CL_SUBMITTED                                 = 0x2,
        /// <summary>Device is currently executing this command.</summary>
        CL_RUNNING                                   = 0x1,
        /// <summary>Command has completed with no errors.</summary>
        /// Identical to CL_COMPLETE
        CL_SUCCESS                                   = 0,
        /// <summary>Device not found</summary>
        CL_DEVICE_NOT_FOUND                          = -1,
        /// <summary>Device not available</summary>
        CL_DEVICE_NOT_AVAILABLE                      = -2,
        /// <summary>Compiler not available</summary>
        CL_COMPILER_NOT_AVAILABLE                    = -3,
        /// <summary>Memory allocation failure</summary>
        CL_MEM_OBJECT_ALLOCATION_FAILURE             = -4,
        /// <summary>Out of resources</summary>
        CL_OUT_OF_RESOURCES                          = -5,
        /// <summary>Out of host memory</summary>
        CL_OUT_OF_HOST_MEMORY                        = -6,
        /// <summary>Profiling not available</summary>
        CL_PROFILING_INFO_NOT_AVAILABLE              = -7,
        /// <summary>Memory copy overlap</summary>
        CL_MEM_COPY_OVERLAP                          = -8,
        /// <summary>Image format mismatch</summary>
        CL_IMAGE_FORMAT_MISMATCH                     = -9,
        /// <summary>Image format not supported</summary>
        CL_IMAGE_FORMAT_NOT_SUPPORTED                = -10,
        /// <summary>Build program failure</summary>
        CL_BUILD_PROGRAM_FAILURE                     = -11,
        /// <summary>Map failure</summary>
        CL_MAP_FAILURE                               = -12,
        /// <summary>Misaligned sub-buffer offset</summary>
        CL_MISALIGNED_SUB_BUFFER_OFFSET              = -13,
        /// <summary>Exec status error for events in wait-list</summary>
        CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST = -14,
        /// <summary>Invalid value</summary>
        CL_INVALID_VALUE                             = -30,
        /// <summary>Invalid device type</summary>
        CL_INVALID_DEVICE_TYPE                       = -31,
        /// <summary>Invalid platform</summary>
        CL_INVALID_PLATFORM                          = -32,
        /// <summary>Invalid device</summary>
        CL_INVALID_DEVICE                            = -33,
        /// <summary>Invalid context</summary>
        CL_INVALID_CONTEXT                           = -34,
        /// <summary>Invalid command queue properties</summary>
        CL_INVALID_QUEUE_PROPERTIES                  = -35,
        /// <summary>Invalid command queue</summary>
        CL_INVALID_COMMAND_QUEUE                     = -36,
        /// <summary>Invalid host pointer</summary>
        CL_INVALID_HOST_PTR                          = -37,
        /// <summary>Invalid memory object</summary>
        CL_INVALID_MEM_OBJECT                        = -38,
        /// <summary>Invalid image format descriptor</summary>
        CL_INVALID_IMAGE_FORMAT_DESCRIPTOR           = -39,
        /// <summary>Invalid image size</summary>
        CL_INVALID_IMAGE_SIZE                        = -40,
        /// <summary>Invalid sampler</summary>
        CL_INVALID_SAMPLER                           = -41,
        /// <summary>Invalid binary</summary>
        CL_INVALID_BINARY                            = -42,
        /// <summary>Invalid build options</summary>
        CL_INVALID_BUILD_OPTIONS                     = -43,
        /// <summary>Invalid programm</summary>
        CL_INVALID_PROGRAM                           = -44,
        /// <summary>Invalid program executable</summary>
        CL_INVALID_PROGRAM_EXECUTABLE                = -45,
        /// <summary>Invalid kernel name</summary>
        CL_INVALID_KERNEL_NAME                       = -46,
        /// <summary>Invalid kernel definition</summary>
        CL_INVALID_KERNEL_DEFINITION                 = -47,
        /// <summary>Invalid kernel</summary>
        CL_INVALID_KERNEL                            = -48,
        /// <summary>Invalid argument index</summary>
        CL_INVALID_ARG_INDEX                         = -49,
        /// <summary>Invalid argument value</summary>
        CL_INVALID_ARG_VALUE                         = -50,
        /// <summary>Invalid argument size</summary>
        CL_INVALID_ARG_SIZE                          = -51,
        /// <summary>Invalid kernel arguments</summary>
        CL_INVALID_KERNEL_ARGS                       = -52,
        /// <summary>Invalid work dimension</summary>
        CL_INVALID_WORK_DIMENSION                    = -53,
        /// <summary>Invalid work group size</summary>
        CL_INVALID_WORK_GROUP_SIZE                   = -54,
        /// <summary>Invalid work item size</summary>
        CL_INVALID_WORK_ITEM_SIZE                    = -55,
        /// <summary>Invalid global offset</summary>
        CL_INVALID_GLOBAL_OFFSET                     = -56,
        /// <summary>Invalid event wait-list</summary>
        CL_INVALID_EVENT_WAIT_LIST                   = -57,
        /// <summary>Invalid event</summary>
        CL_INVALID_EVENT                             = -58,
        /// <summary>Invalid operation</summary>
        CL_INVALID_OPERATION                         = -59,
        /// <summary>Invalid OpenGL buffer</summary>
        CL_INVALID_GL_OBJECT                         = -60,
        /// <summary>Invalid buffer size</summary>
        CL_INVALID_BUFFER_SIZE                       = -61,
        /// <summary>Invalid MIP level</summary>
        CL_INVALID_MIP_LEVEL                         = -62,
        /// <summary>Invalid global work size</summary>
        CL_INVALID_GLOBAL_WORK_SIZE                  = -63
    }

    /// <summary>
    /// OpenCL context info
    /// </summary>
    public enum cl_context_info : uint
    {
        /// <summary>Return the context reference count.</summary>
        /// The reference count returned should be considered immediately stale. It is unsuitable for general use in applications. This feature is provided for identifying memory leaks.
        CL_CONTEXT_REFERENCE_COUNT                   = 0x1080,
        /// <summary>Return the list of devices in context.</summary>
        CL_CONTEXT_DEVICES                           = 0x1081,
        /// <summary>Return the properties argument specified in <see cref="cl.CreateContext"/>.</summary>
        CL_CONTEXT_PROPERTIES                        = 0x1082,
        /// <summary>Return the number of devices in context.</summary>
        CL_CONTEXT_NUM_DEVICES                       = 0x1083,
    }

    /// <summary>
    /// OpenCL context info structure
    /// </summary>
    public struct cl_context_info_return
    {
        /// <summary>The context reference count.</summary>
        public int reference_count;
        /// <summary>The list of devices in context.</summary>
        public cl_device_id[] devices;
        /// <summary>The properties argument specified in <see cref="cl.CreateContext"/>.</summary>
        public cl_context_properties[] properties;
    }

    /// <summary>
    /// OpenCL platform info
    /// </summary>
    public enum cl_platform_info : uint
    {
        /// <summary>OpenCL profile string. Returns the profile name supported by the implementation.</summary>
        /// The profile name returned can be one of the following strings:
        /// <list type="bullet">
        ///     <item>
        ///         <term>FULL_PROFILE</term>
        ///         <description>if the implementation supports the OpenCL specification (functionality defined as part of the core specification and does not require any extensions to be supported).</description>
        ///     </item>
        ///     <item>
        ///         <term>EMBEDDED_PROFILE</term>
        ///         <description>if the implementation supports the OpenCL embedded profile. The embedded profile is defined to be a subset for each version of OpenCL.</description>
        ///     </item>
        /// </list>
        CL_PLATFORM_PROFILE                          = 0x0900,

        /// <summary>OpenCL version string. Returns the OpenCL version supported by the implementation.</summary>
        /// This version string has the following format: OpenCL &lt;space&gt;&lt;major_version.minor_version&gt;&lt;space&gt;&lt;platform-specific information&gt;
        CL_PLATFORM_VERSION                          = 0x0901,

        /// <summary>Platform name string.</summary>
        CL_PLATFORM_NAME                             = 0x0902,

        /// <summary>Platform vendor string.</summary>
        CL_PLATFORM_VENDOR                           = 0x0903,

        /// <summary>Returns a space-separated list of extension names (the extension names themselves do not contain any spaces) supported by the platform.</summary>
        /// Extensions defined here must be supported by all devices associated with this platform.
        CL_PLATFORM_EXTENSIONS                       = 0x0904
    }

    /// <summary>
    /// OpenCL platform info structure
    /// </summary>
    public struct cl_platform_info_return
    {
        /// <summary>OpenCL profile string. Returns the profile name supported by the implementation.</summary>
        /// The profile name returned can be one of the following strings:
        /// <list type="bullet">
        ///     <item>
        ///         <term>FULL_PROFILE</term>
        ///         <description>if the implementation supports the OpenCL specification (functionality defined as part of the core specification and does not require any extensions to be supported).</description>
        ///     </item>
        ///     <item>
        ///         <term>EMBEDDED_PROFILE</term>
        ///         <description>if the implementation supports the OpenCL embedded profile. The embedded profile is defined to be a subset for each version of OpenCL.</description>
        ///     </item>
        /// </list>
        public string profile;

        /// <summary>OpenCL version string. Returns the OpenCL version supported by the implementation.</summary>
        /// This version string has the following format: OpenCL &lt;space&gt;&lt;major_version.minor_version&gt;&lt;space&gt;&lt;platform-specific information&gt;
        public string version;

        /// <summary>Platform name string.</summary>
        public string name;

        /// <summary>Platform vendor string.</summary>
        public string vendor;

        /// <summary>Returns a space-separated list of extension names (the extension names themselves do not contain any spaces) supported by the platform.</summary>
        /// Extensions defined here must be supported by all devices associated with this platform.
        public string extensions;
    }

    /// <summary>
    /// OpenCL context properties
    /// </summary>
    public enum cl_context_properties
    {
        /// <summary>No doc</summary>
        CL_CONTEXT_PLATFORM                          = 0x1084
    }

    /// <summary>
    /// OpenCL command queue properties
    /// </summary>
    [Flags]
    public enum cl_command_queue_properties : ulong
    {
        /// <summary>Determines whether the commands queued in the command-queue are executed in-order or out-of-order.</summary>
        /// If set, the commands in the command-queue are executed out-of-order. Otherwise, commands are executed in-order.
        CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE       = (1 << 0),

        /// <summary>Enable or disable profiling of commands in the command-queue.</summary>
        /// If set, the profiling of commands is enabled. Otherwise profiling of commands is disabled. See clGetEventProfilingInfo for more information.
        CL_QUEUE_PROFILING_ENABLE                    = (1 << 1)
    }

    /// <summary>
    /// OpenCL device memory flags
    /// </summary>
    [Flags]
    public enum cl_mem_flags : ulong
    {
        /// <summary>This flag specifies that the memory object will be read and written by a kernel. This is the default.</summary>
        CL_MEM_READ_WRITE                            = (1 << 0),

        /// <summary>This flags specifies that the memory object will be written but not read by a kernel.</summary>
        /// Reading from a buffer or image object created with CL_MEM_WRITE_ONLY inside a kernel is undefined.
        CL_MEM_WRITE_ONLY                            = (1 << 1),

        /// <summary>This flag specifies that the memory object is a read-only memory object when used inside a kernel.</summary>
        /// Writing to a buffer or image object created with CL_MEM_READ_ONLY inside a kernel is undefined.
        CL_MEM_READ_ONLY                             = (1 << 2),

        /// <summary>This flag is valid only if host_ptr is not NULL. If specified, it indicates that the application wants the OpenCL implementation to use memory referenced by host_ptr as the storage bits for the memory object.</summary>
        /// OpenCL implementations are allowed to cache the buffer contents pointed to by host_ptr in device memory. This cached copy can be used when kernels are executed on a device.
        /// The result of OpenCL commands that operate on multiple buffer objects created with the same host_ptr or overlapping host regions is considered to be undefined.
        CL_MEM_USE_HOST_PTR                          = (1 << 3),

        /// <summary>
        /// This flag specifies that the application wants the OpenCL implementation to allocate memory from host accessible memory.
        /// </summary>
        /// CL_MEM_ALLOC_HOST_PTR and CL_MEM_USE_HOST_PTR are mutually exclusive.
        CL_MEM_ALLOC_HOST_PTR                        = (1 << 4),

        /// <summary>This flag is valid only if host_ptr is not NULL. If specified, it indicates that the application wants the OpenCL implementation to allocate memory for the memory object and copy the data from memory referenced by host_ptr.</summary>
        /// CL_MEM_COPY_HOST_PTR and CL_MEM_USE_HOST_PTR are mutually exclusive.
        /// CL_MEM_COPY_HOST_PTR can be used with CL_MEM_ALLOC_HOST_PTR to initialize the contents of the cl_mem object allocated using host-accessible (e.g. PCIe) memory. 
        CL_MEM_COPY_HOST_PTR                         = (1 << 5)
    }

    /// <summary>
    /// OpenCL memory mapping flags
    /// </summary>
    [Flags]
    public enum cl_map_flags : ulong
    {
        /// <summary>Map for reading</summary>
        CL_MAP_READ                                  = (1 << 0),
        /// <summary>Map for writing</summary>
        CL_MAP_WRITE                                 = (1 << 1)
    }

    /// <summary>
    /// OpenCL program build info
    /// </summary>
    public enum cl_program_build_info : uint
    {
        /// <summary>Returns the build status of program for a specific device as given by device.</summary>
        CL_PROGRAM_BUILD_STATUS                      = 0x1181,

        /// <summary>Return the build options specified by the options argument in <see cref="cl.BuildProgram"/> for device.</summary>
        /// If build status of program for device is CL_BUILD_NONE, an empty string is returned.
        CL_PROGRAM_BUILD_OPTIONS                     = 0x1182,

        /// <summary>Return the build log when <see cref="cl.BuildProgram"/> was called for device.</summary>
        /// If build status of program for device is CL_BUILD_NONE, an empty string is returned.
        CL_PROGRAM_BUILD_LOG                         = 0x1183
    }

    /// <summary>
    /// OpenCL program build info structure
    /// </summary>
    public struct cl_program_build_info_return
    {
        /// <summary>Returns the build status of program for a specific device as given by device.</summary>
        public cl_build_status status;

        /// <summary>Return the build options specified by the options argument in <see cref="cl.BuildProgram"/> for device.</summary>
        /// If build status of program for device is CL_BUILD_NONE, an empty string is returned.
        public string options;

        /// <summary>Return the build log when <see cref="cl.BuildProgram"/> was called for device.</summary>
        /// If build status of program for device is CL_BUILD_NONE, an empty string is returned.
        public string log;
    }

    /// <summary>
    /// OpenCL program build status
    /// </summary>
    public enum cl_build_status : int
    {
        /// <summary>The build status returned if no build has been performed on the specified program object for device.</summary>
        CL_BUILD_SUCCESS                             = 0,
        /// <summary>The build status returned if no build has been performed on the specified program object for device.</summary>
        CL_BUILD_NONE                                = -1,
        /// <summary>The build status returned if the last call to clBuildProgram on the specified program object for device generated an error.</summary>
        CL_BUILD_ERROR                               = -2,
        /// <summary>The build status returned if the last call to clBuildProgram on the specified program object for device has not finished.</summary>
        CL_BUILD_IN_PROGRESS                         = -3
    }

    /// <summary>
    /// OpenCL kernel info
    /// </summary>
    public enum cl_kernel_work_group_info : uint
    {
        /// <summary>This provides a mechanism for the application to query the maximum work-group size that can be used to execute a kernel on a specific device given by device.</summary>
        /// The OpenCL implementation uses the resource requirements of the kernel (register usage etc.) to determine what this work-group size should be.
        CL_KERNEL_WORK_GROUP_SIZE                    = 0x11B0,
        
        /// <summary>Returns the work-group size specified by the __attribute__((reqd_work_gr oup_size(X, Y, Z))) qualifier.</summary>
        /// See Function Qualifiers. If the work-group size is not specified using the above attribute qualifier (0, 0, 0) is returned.
        CL_KERNEL_COMPILE_WORK_GROUP_SIZE            = 0x11B1,
        
        /// <summary>Returns the amount of local memory in bytes being used by a kernel.</summary>
        /// This includes local memory that may be needed by an implementation to execute the kernel, variables declared inside the kernel with the __local address qualifier and local memory to be allocated for arguments to the kernel declared as pointers with the __local address qualifier and whose size is specified with clSetKernelArg.
        /// If the local memory size, for any pointer argument to the kernel declared with the __local address qualifier, is not specified, its size is assumed to be 0.
        CL_KERNEL_LOCAL_MEM_SIZE                     = 0x11B2,

#if OPENCL_1_1
        /// <summary>Returns the preferred multiple of workgroup size for launch.</summary>
        /// This is a performance hint. Specifying a workgroup size that is not a multiple of the value returned by this query as the value of the local work size argument to
        /// clEnqueueNDRangeKernelwill not fail to enqueue the kernel for execution unless the work-group size specified is larger than the device maximum.
        CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE = 0x11B3,
        
        /// <summary>Returns the minimum amount of private memory, in bytes, used by each workitem in the kernel.</summary>
        /// This value may include any private memory needed by an implementation to execute the kernel, including that used by the language built-ins and variable declared inside the kernel with the __private qualifier.
        CL_KERNEL_PRIVATE_MEM_SIZE                   = 0x11B4,
#endif
    }

    /// <summary>
    /// OpenCL kernel info structure
    /// </summary>
    public struct cl_kernel_work_group_info_return
    {
        /// <summary>This provides a mechanism for the application to query the maximum work-group size that can be used to execute a kernel on a specific device given by device.</summary>
        /// The OpenCL implementation uses the resource requirements of the kernel (register usage etc.) to determine what this work-group size should be.
        public uint work_group_size;

        /// <summary>Returns the work-group size specified by the __attribute__((reqd_work_gr oup_size(X, Y, Z))) qualifier.</summary>
        /// See Function Qualifiers. If the work-group size is not specified using the above attribute qualifier (0, 0, 0) is returned.
        public uint[] compile_work_group_size;

        /// <summary>Returns the amount of local memory in bytes being used by a kernel.</summary>
        /// This includes local memory that may be needed by an implementation to execute the kernel, variables declared inside the kernel with the __local address qualifier and local memory to be allocated for arguments to the kernel declared as pointers with the __local address qualifier and whose size is specified with clSetKernelArg.
        /// If the local memory size, for any pointer argument to the kernel declared with the __local address qualifier, is not specified, its size is assumed to be 0.
        public ulong local_mem_size;

#if OPENCL_1_1
        /// <summary>Returns the preferred multiple of workgroup size for launch.</summary>
        /// This is a performance hint. Specifying a workgroup size that is not a multiple of the value returned by this query as the value of the local work size argument to
        /// clEnqueueNDRangeKernelwill not fail to enqueue the kernel for execution unless the work-group size specified is larger than the device maximum.
        public uint preferred_work_group_size_multiple;

        /// <summary>Returns the minimum amount of private memory, in bytes, used by each workitem in the kernel.</summary>
        /// This value may include any private memory needed by an implementation to execute the kernel, including that used by the language built-ins and variable declared inside the kernel with the __private qualifier. 
        public ulong private_mem_size;
#endif
    }

    /// <summary>
    /// OpenCL event info
    /// </summary>
    public enum cl_event_info : uint
    {
        /// <summary>Return the command-queue associated with event. For user event objects, a NULL value is returned.</summary>
        CL_EVENT_COMMAND_QUEUE                       = 0x11D0,
        /// <summary>Return the command associated with event.</summary>
        CL_EVENT_COMMAND_TYPE                        = 0x11D1,
        /// <summary>Return the event reference count.</summary>
        ///  The reference count returned should be considered immediately stale. It is unsuitable for general use in applications. This feature is provided for identifying memory leaks. 
        CL_EVENT_REFERENCE_COUNT                     = 0x11D2,
        /// <summary>Return the execution status of the command identified by event.</summary>
        CL_EVENT_COMMAND_EXECUTION_STATUS            = 0x11D3,
#if OPENCL_1_1
        /// <summary>Return the context associated with event.</summary>
        CL_EVENT_CONTEXT                             = 0x11D4,
#endif
    }

    /// <summary>
    /// OpenCL event info structure
    /// </summary>
    public struct cl_event_info_return
    {
        /// <summary>Return the command-queue associated with event. For user event objects, a NULL value is returned.</summary>
        public cl_command_queue queue;
        /// <summary>Return the command associated with event.</summary>
        public cl_command_type type;
        /// <summary>Return the event reference count.</summary>
        /// The reference count returned should be considered immediately stale. It is unsuitable for general use in applications. This feature is provided for identifying memory leaks. 
        public uint reference_count;
        /// <summary>Return the execution status of the command identified by event.</summary>
        public cl_status status;
#if OPENCL_1_1
        /// <summary>Return the context associated with event.</summary>
        public cl_context context;
#endif
    }

    /// <summary>
    /// OpenCL command type
    /// </summary>
    public enum cl_command_type : uint
    {
        /// <summary>OpenCL kernel.</summary>
        CL_COMMAND_NDRANGE_KERNEL                    = 0x11F0,
        /// <summary>Task</summary>
        CL_COMMAND_TASK                              = 0x11F1,
        /// <summary>Native kernel</summary>
        CL_COMMAND_NATIVE_KERNEL                     = 0x11F2,
        /// <summary>Read buffer operation</summary>
        CL_COMMAND_READ_BUFFER                       = 0x11F3,
        /// <summary>Write buffer operation</summary>
        CL_COMMAND_WRITE_BUFFER                      = 0x11F4,
        /// <summary>Copy buffer operation</summary>
        CL_COMMAND_COPY_BUFFER                       = 0x11F5,
        /// <summary>Read image operation</summary>
        CL_COMMAND_READ_IMAGE                        = 0x11F6,
        /// <summary>Write image operation</summary>
        CL_COMMAND_WRITE_IMAGE                       = 0x11F7,
        /// <summary>Copy image operation</summary>
        CL_COMMAND_COPY_IMAGE                        = 0x11F8,
        /// <summary>Copy image to buffer</summary>
        CL_COMMAND_COPY_IMAGE_TO_BUFFER              = 0x11F9,
        /// <summary>Copy buffer to image</summary>
        CL_COMMAND_COPY_BUFFER_TO_IMAGE              = 0x11FA,
        /// <summary>Map buffer operation</summary>
        CL_COMMAND_MAP_BUFFER                        = 0x11FB,
        /// <summary>map image oepration</summary>
        CL_COMMAND_MAP_IMAGE                         = 0x11FC,
        /// <summary>Unmap memory object</summary>
        CL_COMMAND_UNMAP_MEM_OBJECT                  = 0x11FD,
        /// <summary>Marker</summary>
        CL_COMMAND_MARKER                            = 0x11FE,
        /// <summary>Get OpenGL objects</summary>
        CL_COMMAND_ACQUIRE_GL_OBJECTS                = 0x11FF,
        /// <summary>Release OpenGL objects</summary>
        CL_COMMAND_RELEASE_GL_OBJECTS                = 0x1200,
        /// <summary>Read buffer region</summary>
        CL_COMMAND_READ_BUFFER_RECT                  = 0x1201,
        /// <summary>Write buffer region</summary>
        CL_COMMAND_WRITE_BUFFER_RECT                 = 0x1202,
        /// <summary>Copy buffer region</summary>
        CL_COMMAND_COPY_BUFFER_RECT                  = 0x1203,
        /// <summary>User event</summary>
        CL_COMMAND_USER                              = 0x1204
    }

    /// <summary>
    /// OpenCL program info
    /// </summary>
    public enum cl_program_info : uint
    {
        /// <summary>Return the program reference count.</summary>
        CL_PROGRAM_REFERENCE_COUNT                   = 0x1160,
        /// <summary>Return the context specified when the program object is created.</summary>
        CL_PROGRAM_CONTEXT                           = 0x1161,
        /// <summary>Return the number of devices associated with program.</summary>
        CL_PROGRAM_NUM_DEVICES = 0x1162,
        /// <summary>Return the list of devices associated with the program object.</summary>
        /// This can be the devices associated with context on which the program object has been created or can be a subset of devices that are specified when a progam object is created using clCreateProgramWithBinary.
        CL_PROGRAM_DEVICES = 0x1163,
        /// <summary>Return the program source code specified by clCreateProgramWithSource.</summary>
        CL_PROGRAM_SOURCE                            = 0x1164,
        /// <summary>Returns an array that contains the size in bytes of the program binary for each device associated with program.</summary>
        CL_PROGRAM_BINARY_SIZES                      = 0x1165,
        /// <summary>Return the program binaries for all devices associated with program.</summary>
        CL_PROGRAM_BINARIES                          = 0x1166
    }

    /// <summary>
    /// OpenCL program info structure
    /// </summary>
    public struct cl_program_info_return
    {
        /// <summary>Return the program reference count.</summary>
        public uint reference_count;
        /// <summary>Return the number of devices associated with program.</summary>
        public cl_context context;
        /// <summary>Return the list of devices associated with the program object.</summary>
        /// This can be the devices associated with context on which the program object has been created or can be a subset of devices that are specified when a progam object is created using clCreateProgramWithBinary.
        public cl_device_id[] devices;
        /// <summary>Return the program source code specified by clCreateProgramWithSource.</summary>
        public string source;
        /// <summary>Return the program binaries for all devices associated with program.</summary>
        public byte[][] binaries;
    }

    /// <summary>
    /// OpenCL kernel info
    /// </summary>
    public enum cl_kernel_info : uint
    {
        /// <summary>Return the kernel function name.</summary>
        CL_KERNEL_FUNCTION_NAME                      = 0x1190,
        /// <summary>Return the number of arguments to kernel.</summary>
        CL_KERNEL_NUM_ARGS                           = 0x1191,
        /// <summary>Return the kernel reference count.</summary>
        CL_KERNEL_REFERENCE_COUNT                    = 0x1192,
        /// <summary>Return the context associated with kernel.</summary>
        CL_KERNEL_CONTEXT                            = 0x1193,
        /// <summary>Return the program object associated with kernel.</summary>
        CL_KERNEL_PROGRAM                            = 0x1194
    }

    /// <summary>
    /// OpenCL kernel info structure
    /// </summary>
    public struct cl_kernel_info_return
    {
        /// <summary>Return the kernel function name.</summary>
        public string kernel_func_name;
        /// <summary>Return the number of arguments to kernel.</summary>
        public uint kernel_num_args;
        /// <summary>Return the kernel reference count.</summary>
        public uint reference_count;
        /// <summary>Return the context associated with kernel.</summary>
        public cl_context kernel_context;
        /// <summary>Return the program object associated with kernel.</summary>
        public cl_program kernel_program;
    }
    #endregion
    // ------------------------------------------
    #region Callback Delegates
    /// <summary>
    /// A callback function that can be registered by the application.
    /// </summary>
    /// This callback function will be used by the OpenCL implementation to report information on errors that occur in this context.
    /// This callback function may be called asynchronously by the OpenCL implementation. It is the application's responsibility to ensure that the callback function is thread-safe.
    /// If pfn_notify is NULL, no callback function is registered. The parameters to this callback function are:
    /// <param name="errinfo">A pointer to an error string.</param>
    /// <param name="private_info">Represent a pointer to binary data that is returned by the OpenCL implementation that can be used to log additional information helpful in debugging the error.</param>
    /// <param name="cb">Size of private_info.</param>
    /// <param name="user_data">A pointer to user supplied data.</param>
    /// NOTE: There are a number of cases where error notifications need to be delivered due to an error that occurs outside a context.
    /// Such notifications may not be delivered through the pfn_notify callback. Where these notifications go is implementation-defined.
    public delegate void cl_create_context_callback(IntPtr errinfo, IntPtr private_info, IntPtr cb, IntPtr user_data);

    /// <summary>
    /// A function pointer to a notification routine.
    /// </summary>
    /// The notification routine is a callback function that an application can register and which will be called when the program executable has been built (successfully or unsuccessfully).
    /// If pfn_notify is not NULL, clBuildProgram does not need to wait for the build to complete and can return immediately.
    /// If pfn_notify is NULL, clBuildProgram does not return until the build has completed.
    /// This callback function may be called asynchronously by the OpenCL implementation. It is the application's responsibility to ensure that the callback function is thread-safe. 
    /// <param name="program">The program.</param>
    /// <param name="user_data">A pointer to user supplied data.</param>
    public delegate void cl_build_program_callback(cl_program program, IntPtr user_data);

    /// <summary>
    /// The event callback function that can be registered by the application.
    /// </summary>
    /// This callback function may be called asynchronously by the OpenCL implementation. It is the application's responsibility to ensure that the callback function is thread-safe.
    /// The parameters to this callback function are: 
    /// <param name="e">The event object for which the callback function is invoked.</param>
    /// <param name="event_command_exec_status">The execution status of command for which this callback function is invoked.</param>
    /// <param name="user_data">A pointer to user supplied data.</param>
    public delegate void cl_set_event_callback(cl_event e, cl_status event_command_exec_status, IntPtr user_data);
    #endregion
    // ------------------------------------------

    /// <summary>
    /// OpenCL API interface
    /// </summary>
    public class cl : DynLibLoader
    {
        public cl()
            : base(
                new string[] { "opencl.dll", "libopencl.so", "libOpenCL.so" },
                new string[3][][],
                new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity, DynLibLoader.Identity },
                new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix },
                new int[] { -1, -1, -1 })
        {

        }

        private static cl my;

        static cl()
        {
            my = new cl();
        }

        private static void testResult(cl_status r)
        {
            if (r != cl_status.CL_SUCCESS)
                throw new ApplicationException("OpenCL error: " + r.ToString());
        }
        
#pragma warning disable        649
        // ------------------------------------------
        #region Platform Management
        private delegate cl_status _clGetPlatformIDs(uint numEntries, IntPtr platforms, out uint num_platforms);
        private delegate cl_status _clGetPlatformInfo(cl_platform_id platform, cl_platform_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);

        _clGetPlatformIDs clGetPlatformIDs;
        _clGetPlatformInfo clGetPlatformInfo;

        /// <summary>
        /// Obtain the list of platforms available.
        /// </summary>
        /// <returns>Returns a list of OpenCL platforms found. The cl_platform_id values returned in platforms can be used to identify a specific OpenCL platform.</returns>
        public static cl_platform_id[] GetPlatformIDs()
        {
            cl_status r;
            uint num_platforms, ignored;
            testResult(my.clGetPlatformIDs(0, IntPtr.Zero, out num_platforms));
            if (num_platforms == 0)
            {
                return null;
            }
            
            cl_platform_id[] platforms = new cl_platform_id[num_platforms];
            unsafe
            {
                fixed (cl_platform_id* ptr = platforms)
                {
                    r = my.clGetPlatformIDs(num_platforms, (IntPtr)ptr, out ignored);
                }
            }
            
            testResult(r);
            return platforms;
        }

        /// <summary>
        /// Get specific information about the OpenCL platform.
        /// </summary>
        /// <param name="platform">The platform ID returned by <see cref="GetPlatformIDs"/>.</param>
        /// <returns>Returns the information.</returns>
        public static cl_platform_info_return GetPlatformInfo(cl_platform_id platform)
        {
            cl_platform_info_return ret = new cl_platform_info_return();

            ret.extensions = GetPlatformInfoString(platform, cl_platform_info.CL_PLATFORM_EXTENSIONS);
            ret.name = GetPlatformInfoString(platform, cl_platform_info.CL_PLATFORM_NAME);
            ret.profile = GetPlatformInfoString(platform, cl_platform_info.CL_PLATFORM_PROFILE);
            ret.vendor = GetPlatformInfoString(platform, cl_platform_info.CL_PLATFORM_VENDOR);
            ret.version = GetPlatformInfoString(platform, cl_platform_info.CL_PLATFORM_VERSION);

            return ret;
        }

        internal static string GetPlatformInfoString(cl_platform_id platform, cl_platform_info param_name)
        {
            IntPtr len, ignored;

            testResult(my.clGetPlatformInfo(platform, param_name, IntPtr.Zero, IntPtr.Zero, out len));
            IntPtr mem = Marshal.AllocHGlobal((int)len);
            cl_status r = my.clGetPlatformInfo(platform, param_name, len, mem, out ignored);
            string param_value = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            testResult(r);

            return param_value;
        }
        #endregion
        // ------------------------------------------
        #region Context Management
        private delegate cl_context _clCreateContext(IntPtr properties, uint num_devices, IntPtr devices, cl_create_context_callback pfn_notify, IntPtr user_data, out cl_status errcode_ret);
        private delegate cl_context _clCreateContextFromType(IntPtr properties, cl_device_type device_type, cl_create_context_callback pfn_notify, IntPtr user_data, out cl_status errcode_ret);
        private delegate cl_status _clGetContextInfo(cl_context context, cl_context_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);
        private delegate cl_status _clRetainContext(cl_context context);
        private delegate cl_status _clReleaseContext(cl_context context);

        _clCreateContext clCreateContext;
        _clCreateContextFromType clCreateContextFromType;
        _clGetContextInfo clGetContextInfo;
        _clRetainContext clRetainContext;
        _clReleaseContext clReleaseContext;
		
        /// <summary>
        /// Creates an OpenCL context.
        /// </summary>
        /// An OpenCL context is created with one or more devices. Contexts are used by the OpenCL runtime for managing objects such as command-queues, memory, program and kernel objects and for executing kernels on one or more devices specified in the context.
        /// <param name="platform">Specifies the platform to use.</param>
        /// <param name="devices">A list of unique devices returned by <see cref="GetDeviceIDs"/> for a platform. </param>
        /// <param name="pfn_notify">Callback, see <see cref="cl_create_context_callback"/>.</param>
        /// <param name="user_data">Passed as the user_data argument when pfn_notify is called. user_data can be NULL.</param>
        /// <returns>Returns the context.</returns>
        public static cl_context CreateContext(cl_platform_id platform, cl_device_id[] devices, cl_create_context_callback pfn_notify = null, IntPtr user_data = default(IntPtr))
        {
            cl_context ctx;
            cl_status r;
            IntPtr[] prop = new IntPtr[3];
            prop[0] = (IntPtr)cl_context_properties.CL_CONTEXT_PLATFORM;
            prop[1] = platform.p;
            prop[2] = IntPtr.Zero;

            unsafe
            {
                fixed (IntPtr* platformAddress = prop)
                {
                    fixed (cl_device_id* deviceAddress = devices)
                    {
                        ctx = my.clCreateContext((IntPtr)platformAddress, (uint)devices.Length, (IntPtr)deviceAddress, pfn_notify, user_data, out r);
                    }
                }
            }
            testResult(r);
            return ctx;
        }

        /// <summary>
        /// Create an OpenCL context from a device type that identifies the specific device(s) to use. 
        /// </summary>
        /// <param name="platform">Specifies the platform to use.</param>
        /// <param name="device_type">Identifies the type of device.</param>
        /// <param name="pfn_notify">Callback, see <see cref="cl_create_context_callback"/>.</param>
        /// <param name="user_data">Passed as the user_data argument when pfn_notify is called. user_data can be NULL.</param>
        /// <returns>Returns the context.</returns>
        public static cl_context CreateContextFromType(cl_platform_id platform, cl_device_type device_type, cl_create_context_callback pfn_notify = null, IntPtr user_data = default(IntPtr))
        {
            cl_context ctx;
            cl_status r;
            IntPtr[] prop = new IntPtr[3];
            prop[0] = (IntPtr)cl_context_properties.CL_CONTEXT_PLATFORM;
            prop[1] = platform.p;
            prop[2] = IntPtr.Zero;

            unsafe
            {
                fixed (IntPtr* platformAddress = prop)
                {
                    ctx = my.clCreateContextFromType((IntPtr)platformAddress, device_type, pfn_notify, user_data, out r);
                }
            }
            testResult(r);
            return ctx;
        }

        /// <summary>
        /// Query information about a context.
        /// </summary>
        /// <param name="context">Specifies the OpenCL context being queried.</param>
        /// <returns>Returns a struct with all information.</returns>
        public static cl_context_info_return GetContextInfo(cl_context context)
        {
            cl_context_info_return ret = new cl_context_info_return();
            IntPtr size, ignored;
            int count;

            unsafe
            {
                size = (IntPtr)sizeof(uint);
                int* pref = &ret.reference_count;
                testResult(my.clGetContextInfo(context, cl_context_info.CL_CONTEXT_REFERENCE_COUNT, size, (IntPtr)pref, out ignored));

                testResult(my.clGetContextInfo(context, cl_context_info.CL_CONTEXT_DEVICES, IntPtr.Zero, IntPtr.Zero, out size));
                count = (int)size / IntPtr.Size;
                ret.devices = new cl_device_id[count];
                fixed (cl_device_id* p = ret.devices)
                {
                    testResult(my.clGetContextInfo(context, cl_context_info.CL_CONTEXT_DEVICES, size, (IntPtr)p, out ignored));
                }

                testResult(my.clGetContextInfo(context, cl_context_info.CL_CONTEXT_PROPERTIES, IntPtr.Zero, IntPtr.Zero, out size));
                count = (int)size / IntPtr.Size;
                ret.properties = new cl_context_properties[count];
                fixed (cl_context_properties* p = ret.properties)
                {
                    testResult(my.clGetContextInfo(context, cl_context_info.CL_CONTEXT_DEVICES, size, (IntPtr)p, out ignored));
                }
            }

            return ret;
        }

        /// <summary>
        /// Increment the context reference count. 
        /// </summary>
        /// clCreateContext and clCreateContextFromType perform an implicit retain.
        /// <param name="context">The context to retain.</param>
        public static void RetainContext(cl_context context)
        {
            testResult(my.clRetainContext(context));
        }

        /// <summary>
        /// Decrement the context reference count.
        /// </summary>
        /// <param name="context">The context to release.</param>
        public static void ReleaseContext(cl_context context)
        {
            testResult(my.clReleaseContext(context));
        }
        #endregion
        // ------------------------------------------
        #region Device Management
        private delegate cl_status _clGetDeviceIDs(cl_platform_id platform, cl_device_type device_type, uint num_entries, IntPtr devices, out uint num_devices);
        private delegate cl_status _clGetDeviceInfo(cl_device_id device, cl_device_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);
        
        _clGetDeviceIDs clGetDeviceIDs;
        _clGetDeviceInfo clGetDeviceInfo;

        /// <summary>
        /// Obtain the list of devices available on a platform.
        /// </summary>
        /// <param name="platform">Refers to the platform ID returned by <see cref="GetPlatformIDs"/>.</param>
        /// <param name="device_type">The device_type can be used to query specific OpenCL devices or all OpenCL devices available.</param>
        /// <returns>A list of OpenCL devices found.</returns>
        public static cl_device_id[] GetDeviceIDs(cl_platform_id platform, cl_device_type device_type = cl_device_type.CL_DEVICE_TYPE_ALL)
        {
            uint num_devices, ignored;
            testResult(my.clGetDeviceIDs(platform, device_type, 0, IntPtr.Zero, out num_devices));
            if (num_devices == 0)
            {
                return null;
            }

            cl_device_id[] devices = new cl_device_id[num_devices];
            unsafe
            {
                fixed (cl_device_id* ptr = devices)
                {
                    testResult(my.clGetDeviceIDs(platform, device_type, num_devices, (IntPtr)ptr, out ignored));
                }
            }
            return devices;
        }

        /// <summary>
        /// Get information about an OpenCL device.
        /// </summary>
        /// <param name="device">Refers to the device returned by clGetDeviceIDs.</param>
        /// <returns>A structure with device information.</returns>
        public static cl_device_info_return GetDeviceInfo(cl_device_id device)
        {
            cl_device_info_return ret = new cl_device_info_return();
            IntPtr size;

            ret.address_bits = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_ADDRESS_BITS);
            ret.available = GetDeviceInfoBool(device, cl_device_info.CL_DEVICE_AVAILABLE);
            ret.compiler_available = GetDeviceInfoBool(device, cl_device_info.CL_DEVICE_COMPILER_AVAILABLE);
            ret.double_fp_config = (cl_device_fp_config)GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_DOUBLE_FP_CONFIG);
            //ret.half_fp_config = (cl_device_fp_config)GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_HALF_FP_CONFIG);
            ret.single_fp_config = (cl_device_fp_config)GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_SINGLE_FP_CONFIG);
            ret.endian_little = GetDeviceInfoBool(device, cl_device_info.CL_DEVICE_ENDIAN_LITTLE);
            ret.error_correction_support = GetDeviceInfoBool(device, cl_device_info.CL_DEVICE_ERROR_CORRECTION_SUPPORT);
            ret.execution_capabilities = (cl_device_exec_capabilities)GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_EXECUTION_CAPABILITIES);
            ret.extensions = GetDeviceInfoString(device, cl_device_info.CL_DEVICE_EXTENSIONS);
            ret.global_mem_cache_size = GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_GLOBAL_MEM_CACHE_SIZE);
            ret.global_mem_cache_type = (cl_device_mem_cache_type)GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_GLOBAL_MEM_CACHE_TYPE);
            ret.global_mem_cacheline_size = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE);
            ret.global_mem_size = GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_GLOBAL_MEM_SIZE);
            ret.image_support = GetDeviceInfoBool(device, cl_device_info.CL_DEVICE_IMAGE_SUPPORT);
            ret.image2d_max_height = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_IMAGE2D_MAX_HEIGHT).ToUInt32();
            ret.image2d_max_width = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_IMAGE2D_MAX_WIDTH).ToUInt32();
            ret.image3d_max_depth = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_IMAGE3D_MAX_DEPTH).ToUInt32();
            ret.image3d_max_height = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_IMAGE3D_MAX_HEIGHT).ToUInt32();
            ret.image3d_max_width = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_IMAGE3D_MAX_WIDTH).ToUInt32();
            ret.local_mem_size = GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_LOCAL_MEM_SIZE);
            ret.local_mem_type = (cl_device_local_mem_type)GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_LOCAL_MEM_TYPE);
            ret.max_clock_frequency = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MAX_CLOCK_FREQUENCY);
            ret.max_compute_units = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MAX_COMPUTE_UNITS);
            ret.max_constant_args = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MAX_CONSTANT_ARGS);
            ret.max_constant_buffer_size = GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE);
            ret.max_mem_alloc_size = GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_MAX_MEM_ALLOC_SIZE);
            ret.max_parameter_size = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_MAX_PARAMETER_SIZE).ToUInt32();
            ret.max_read_image_args = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MAX_READ_IMAGE_ARGS);
            ret.max_samplers = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MAX_SAMPLERS);
            ret.max_work_group_size = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_MAX_WORK_GROUP_SIZE).ToUInt32();
            ret.max_work_item_dimensions = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);

            UIntPtr[] val = new UIntPtr[ret.max_work_item_dimensions];
            size = (IntPtr)(ret.max_work_item_dimensions * IntPtr.Size);
            unsafe
            {
                fixed (UIntPtr* p = val)
                {
                    testResult(my.clGetDeviceInfo(device, cl_device_info.CL_DEVICE_MAX_WORK_ITEM_SIZES, size, (IntPtr)p, out size));
                }
            }
            ret.max_work_item_sizes = new uint[ret.max_work_item_dimensions];
            for (uint i = 0; i < ret.max_work_item_dimensions; i++)
            {
                ret.max_work_item_sizes[i] = val[i].ToUInt32();
            }

            ret.max_write_image_args = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MAX_WRITE_IMAGE_ARGS);
            ret.mem_base_addr_align = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MEM_BASE_ADDR_ALIGN);
            ret.min_data_type_align_size = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE);
            ret.name = GetDeviceInfoString(device, cl_device_info.CL_DEVICE_NAME);
            ret.platform.p = GetDeviceInfoIntPtr(device, cl_device_info.CL_DEVICE_PLATFORM);
            ret.profile = GetDeviceInfoString(device, cl_device_info.CL_DEVICE_PROFILE);
            ret.profiling_timer_resultion = GetDeviceInfoUIntPtr(device, cl_device_info.CL_DEVICE_PROFILING_TIMER_RESOLUTION).ToUInt32();
            ret.queue_properties = (cl_command_queue_properties)GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_QUEUE_PROPERTIES);
            ret.type = (cl_device_type)GetDeviceInfoULong(device, cl_device_info.CL_DEVICE_TYPE);
            ret.vendor = GetDeviceInfoString(device, cl_device_info.CL_DEVICE_VENDOR);
            ret.vendor_id = GetDeviceInfoUInt(device, cl_device_info.CL_DEVICE_VENDOR_ID);
            ret.version = GetDeviceInfoString(device, cl_device_info.CL_DEVICE_VERSION);
#if OPENCL_1_1
            ret.host_unified_memory = GetDeviceInfoBool(device, cl_device_info.CL_DEVICE_HOST_UNIFIED_MEMORY);
            ret.opencl_c_version = GetDeviceInfoString(device, cl_device_info.CL_DEVICE_OPENCL_C_VERSION);
#endif

            return ret;
        }

        internal static string GetDeviceInfoString(cl_device_id device, cl_device_info param_name)
        {
            IntPtr len, ignored;
            cl_status r;

            testResult(my.clGetDeviceInfo(device, param_name, IntPtr.Zero, IntPtr.Zero, out len));
            IntPtr mem = Marshal.AllocHGlobal((int)len);
            r = my.clGetDeviceInfo(device, param_name, len, mem, out ignored);
            string ret = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            testResult(r);

            return ret;
        }

        internal static uint GetDeviceInfoUInt(cl_device_id device, cl_device_info param_name)
        {
            uint ret;
            IntPtr ignored, size = (IntPtr)sizeof(uint);
            unsafe
            {
                uint* p = &ret;
                testResult(my.clGetDeviceInfo(device, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        internal static ulong GetDeviceInfoULong(cl_device_id device, cl_device_info param_name)
        {
            ulong ret;
            IntPtr ignored, size = (IntPtr)sizeof(ulong);
            unsafe
            {
                ulong* p = &ret;
                testResult(my.clGetDeviceInfo(device, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        internal static IntPtr GetDeviceInfoIntPtr(cl_device_id device, cl_device_info param_name)
        {
            IntPtr ret, ignored, size = (IntPtr)IntPtr.Size;
            unsafe
            {
                IntPtr* p = &ret;
                testResult(my.clGetDeviceInfo(device, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        internal static UIntPtr GetDeviceInfoUIntPtr(cl_device_id device, cl_device_info param_name)
        {
            UIntPtr ret;
            IntPtr ignored, size = (IntPtr)UIntPtr.Size;
            unsafe
            {
                UIntPtr* p = &ret;
                testResult(my.clGetDeviceInfo(device, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        internal static bool GetDeviceInfoBool(cl_device_id device, cl_device_info param_name)
        {
            return GetDeviceInfoUInt(device, param_name) == 0 ? false : true;
        }
        #endregion
        // ------------------------------------------
        #region Command-Queue Management
        private delegate cl_command_queue _clCreateCommandQueue(cl_context context, cl_device_id device, cl_command_queue_properties properties, out cl_status errcode_ret);
        private delegate cl_status _clReleaseCommandQueue(cl_command_queue queue);
        private delegate cl_status _clEnqueueMarker(cl_command_queue queue, out cl_event_id event_ret);
        private delegate cl_status _clEnqueueWaitForEvents(cl_command_queue queue, uint num_events, IntPtr event_list);
        private delegate cl_status _clEnqueueBarrier(cl_command_queue queue);
        private delegate cl_status _clFlush(cl_command_queue queue);
        private delegate cl_status _clFinish(cl_command_queue queue);

        _clCreateCommandQueue clCreateCommandQueue;
        _clReleaseCommandQueue clReleaseCommandQueue;
        _clEnqueueMarker clEnqueueMarker;
        _clEnqueueWaitForEvents clEnqueueWaitForEvents;
        _clEnqueueBarrier clEnqueueBarrier;
        _clFlush clFlush;
        _clFinish clFinish;

        /// <summary>
        /// Create a command-queue on a specific device.
        /// </summary>
        /// <param name="context">Must be a valid OpenCL context. </param>
        /// <param name="device">Must be a device associated with context.</param>
        /// <param name="properties">Specifies a list of properties for the command-queue.</param>
        /// <returns>Returns the command queue.</returns>
        public static cl_command_queue CreateCommandQueue(cl_context context, cl_device_id device, cl_command_queue_properties properties = 0)
        {
            cl_status r;
            cl_command_queue queue = my.clCreateCommandQueue(context, device, properties, out r);
            testResult(r);
            return queue;
        }

        /// <summary>
        /// Decrements the command_queue reference count.
        /// </summary>
        /// <param name="queue">Specifies the command-queue to release. </param>
        public static void ReleaseCommandQueue(cl_command_queue queue)
        {
            testResult(my.clReleaseCommandQueue(queue));
        }

        /// <summary>
        /// Enqueues a marker command.
        /// </summary>
        /// Enqueues a marker command to command_queue. The marker command is not completed until all commands enqueued before it have completed.
        /// The marker command returns an event which can be waited on, i.e. this event can be waited on to ensure that all commands which have been queued before the market command have been completed.
        /// <param name="queue">A valid command-queue.</param>
        /// <returns>Returns an event for this marker.</returns>
        public static cl_event EnqueueMarker(cl_command_queue queue)
        {
            cl_event_id event_ret;
            testResult(my.clEnqueueMarker(queue, out event_ret));
            return new cl_event(event_ret);
        }

        /// <summary>
        /// Enqueues a wait for a specific event or a list of events to complete before any future commands queued in the command-queue are executed.
        /// </summary>
        /// Events specified in event_list act as synchronization points. The context associated with events in event_list and command_queue must be the same.
        /// <param name="queue">A valid command-queue.</param>
        /// <param name="event_list">A list of events.</param>
        public static void EnqueueWaitForEvents(cl_command_queue queue, cl_event[] event_list)
        {
            if (event_list == null)
            {
                return;
            }

            uint num_events_in_event_list = (uint)event_list.Length;
            cl_event_id[] list = cl_event.convertList(event_list);

            unsafe
            {
                fixed (cl_event_id* eventAddress = list)
                {
                    testResult(my.clEnqueueWaitForEvents(queue, num_events_in_event_list, (IntPtr)eventAddress));
                }
            }
        }

        /// <summary>
        /// A synchronization point that enqueues a barrier operation.
        /// </summary>
        /// <param name="queue">A valid command-queue.</param>
        public static void EnqueueBarrier(cl_command_queue queue)
        {
            testResult(my.clEnqueueBarrier(queue));
        }

        /// <summary>
        /// Issues all previously queued OpenCL commands in a command-queue to the device associated with the command-queue.
        /// </summary>
        /// <param name="queue">Specifies the command-queue.</param>
        public static void Flush(cl_command_queue queue)
        {
            testResult(my.clFlush(queue));
        }

        /// <summary>
        /// Blocks until all previously queued OpenCL commands in a command-queue are issued to the associated device and have completed.
        /// </summary>
        /// <param name="queue">Specifies the command-queue.</param>
        public static void Finish(cl_command_queue queue)
        {
            testResult(my.clFinish(queue));
        }
        #endregion
        // ------------------------------------------
        #region Memory Management
        private delegate cl_mem _clCreateBuffer(cl_context context, cl_mem_flags flags, IntPtr size, IntPtr host_ptr, out cl_status errcode_ret);
        private delegate cl_status _clEnqueueReadBuffer(cl_command_queue command_queue, cl_mem buffer, uint blocking_read, IntPtr offset, IntPtr cb, IntPtr ptr, uint num_events_in_wait_list, IntPtr wait_list, out cl_event_id event_ret);
        private delegate cl_status _clEnqueueWriteBuffer(cl_command_queue command_queue, cl_mem buffer, uint blocking_write, IntPtr offset, IntPtr cb, IntPtr ptr, uint num_events_in_wait_list, IntPtr wait_list, out cl_event_id event_ret);
        private delegate cl_status _clEnqueueCopyBuffer(cl_command_queue command_queue, cl_mem src_buffer, cl_mem dest_buffer, IntPtr src_offset, IntPtr dest_offset, IntPtr cb, uint num_events_in_wait_list, IntPtr wait_list, out cl_event_id event_ret);
        private delegate IntPtr _clEnqueueMapBuffer(cl_command_queue command_queue, cl_mem buffer, uint blocking_map, cl_map_flags map_flags, IntPtr offset, IntPtr cb, uint num_events_in_wait_list, IntPtr wait_list, out cl_event_id event_ret, out cl_status errcode_ret);
        private delegate cl_status _clEnqueueUnmapMemObject(cl_command_queue command_queue, cl_mem memobj, IntPtr mapped_ptr, uint num_events_in_wait_list, IntPtr wait_list, out cl_event_id event_ret);
        private delegate cl_status _clReleaseMemObject(cl_mem memobj);

        _clCreateBuffer clCreateBuffer;
        _clEnqueueReadBuffer clEnqueueReadBuffer;
        _clEnqueueWriteBuffer clEnqueueWriteBuffer;
        _clEnqueueCopyBuffer clEnqueueCopyBuffer;
        _clEnqueueMapBuffer clEnqueueMapBuffer;
        _clEnqueueUnmapMemObject clEnqueueUnmapMemObject;
        _clReleaseMemObject clReleaseMemObject;

        /// <summary>
        /// Creates a buffer object.
        /// </summary>
        /// <param name="context">A valid OpenCL context used to create the buffer object.</param>
        /// <param name="flags">A bit-field that is used to specify allocation and usage information such as the memory arena that should be used to allocate the buffer object and how it will be used.</param>
        /// <param name="size">The size in bytes of the buffer memory object to be allocated.</param>
        /// <returns>Returns a valid non-zero buffer object if successfull.</returns>
        public static cl_mem CreateBuffer(cl_context context, cl_mem_flags flags, uint size)
        {
            cl_status r;
            cl_mem mem;
            mem = my.clCreateBuffer(context, flags, (IntPtr)size, IntPtr.Zero, out r);
            testResult(r);
            return mem;
        }

        /// <summary>
        /// Creates a buffer object.
        /// </summary>
        /// <param name="context">A valid OpenCL context used to create the buffer object.</param>
        /// <param name="flags">A bit-field that is used to specify allocation and usage information such as the memory arena that should be used to allocate the buffer object and how it will be used.</param>
        /// <param name="size">The size in bytes of the buffer memory object to be allocated.</param>
        /// <param name="host_ptr">A pointer to the buffer data that may already be allocated by the application. The size of the buffer that host_ptr points to must be greater than or equal to the size bytes.</param>
        /// <returns>Returns a valid non-zero buffer object if successfull.</returns>
        public static cl_mem CreateBuffer(cl_context context, cl_mem_flags flags, uint size, IntPtr host_ptr)
        {
            cl_status r;
            cl_mem mem = my.clCreateBuffer(context, flags, (IntPtr)size, host_ptr, out r);
            testResult(r);
            return mem;
        }

        /// <summary>
        /// Creates a buffer object.
        /// </summary>
        /// <param name="context">A valid OpenCL context used to create the buffer object.</param>
        /// <param name="flags">A bit-field that is used to specify allocation and usage information such as the memory arena that should be used to allocate the buffer object and how it will be used.</param>
        /// <param name="size">The size in bytes of the buffer memory object to be allocated.</param>
        /// <param name="host_ptr">A pointer to the buffer data that may already be allocated by the application. The size of the buffer that host_ptr points to must be greater than or equal to the size bytes.</param>
        /// <returns>Returns a valid non-zero buffer object if successfull.</returns>
        public static cl_mem CreateBuffer(cl_context context, cl_mem_flags flags, uint size, Array host_ptr)
        {
            cl_status r;
            GCHandle handle = GCHandle.Alloc(host_ptr, GCHandleType.Pinned);
            IntPtr address = Marshal.UnsafeAddrOfPinnedArrayElement(host_ptr, 0);
            cl_mem mem = my.clCreateBuffer(context, flags, (IntPtr)size, address, out r);
            handle.Free();
            testResult(r);
            return mem;
        }

        /// <summary>
        /// Enqueue commands to read from a buffer object to host memory.
        /// </summary>
        /// <param name="command_queue">Refers to the command-queue in which the read command will be queued. command_queue and buffer must be created with the same OpenCL context.</param>
        /// <param name="buffer">Refers to a valid buffer object.</param>
        /// <param name="blocking_read">
        ///     Indicates if the read operations are blocking or non-blocking.
        ///     If blocking_read is true i.e. the read command is blocking, clEnqueueReadBuffer does not return until the buffer data has been read and copied into memory pointed to by ptr.
        ///     If blocking_read is false i.e. the read command is non-blocking, clEnqueueReadBuffer queues a non-blocking read command and returns.
        ///     The contents of the buffer that ptr points to cannot be used until the read command has completed. The event argument returns an event object which can be used to query the execution status of the read command.
        ///     When the read command has completed, the contents of the buffer that ptr points to can be used by the application.
        /// </param>
        /// <param name="offset">The offset in bytes in the buffer object to read from.</param>
        /// <param name="cb">The size in bytes of data being read.</param>
        /// <param name="ptr">The pointer to buffer in host memory where data is to be read into.</param>
        /// <param name="wait_list">
        ///     event_wait_list specifies events that need to complete before this particular command can be executed.
        ///     If event_wait_list is null, then this particular command does not wait on any event to complete. The events specified in event_wait_list act as synchronization points.
        ///     The context associated with events in event_wait_list and command_queue must be the same. 
        /// </param>
        /// <returns>Returns an event object that identifies this particular read command and can be used to query or queue a wait for this particular command to complete.</returns>
        public static cl_event EnqueueReadBuffer(cl_command_queue command_queue, cl_mem buffer, bool blocking_read, uint offset, uint cb, IntPtr ptr, cl_event[] wait_list = null)
        {
            cl_event_id event_ret;
            uint _blocking_read = (blocking_read ? 1u : 0u);
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);

            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    testResult(my.clEnqueueReadBuffer(command_queue, buffer, _blocking_read, (IntPtr)offset, (IntPtr)cb, ptr, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret));
                }
            }

            return new cl_event(event_ret);
        }

        /// <summary>
        /// Enqueue commands to read from a buffer object to host memory.
        /// </summary>
        /// <param name="command_queue">Refers to the command-queue in which the read command will be queued. command_queue and buffer must be created with the same OpenCL context.</param>
        /// <param name="buffer">Refers to a valid buffer object.</param>
        /// <param name="blocking_read">
        ///     Indicates if the read operations are blocking or non-blocking.
        ///     If blocking_read is true i.e. the read command is blocking, clEnqueueReadBuffer does not return until the buffer data has been read and copied into memory pointed to by ptr.
        ///     If blocking_read is false i.e. the read command is non-blocking, clEnqueueReadBuffer queues a non-blocking read command and returns.
        ///     The contents of the buffer that ptr points to cannot be used until the read command has completed. The event argument returns an event object which can be used to query the execution status of the read command.
        ///     When the read command has completed, the contents of the buffer that ptr points to can be used by the application.
        /// </param>
        /// <param name="offset">The offset in bytes in the buffer object to read from.</param>
        /// <param name="cb">The size in bytes of data being read.</param>
        /// <param name="ptr">The pointer to buffer in host memory where data is to be read into.</param>
        /// <param name="wait_list">
        ///     event_wait_list specifies events that need to complete before this particular command can be executed.
        ///     If event_wait_list is null, then this particular command does not wait on any event to complete. The events specified in event_wait_list act as synchronization points.
        ///     The context associated with events in event_wait_list and command_queue must be the same. 
        /// </param>
        /// <returns>Returns an event object that identifies this particular read command and can be used to query or queue a wait for this particular command to complete.</returns>
        public static cl_event EnqueueReadBuffer(cl_command_queue command_queue, cl_mem buffer, bool blocking_read, uint offset, uint cb, Array ptr, cl_event[] wait_list = null)
        {
            cl_event_id event_ret;
            uint _blocking_read = (blocking_read ? 1u : 0u);
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);
            
            GCHandle arrayHandle = GCHandle.Alloc(ptr, GCHandleType.Pinned);
            IntPtr arrayAddress = Marshal.UnsafeAddrOfPinnedArrayElement(ptr, 0);
            cl_status r;
            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    r = my.clEnqueueReadBuffer(command_queue, buffer, _blocking_read, (IntPtr)offset, (IntPtr)cb, arrayAddress, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret);
                }
            }
            arrayHandle.Free();
            testResult(r);

            return new cl_event(event_ret);
        }

        /// <summary>
        /// Enqueue commands to write to a buffer object from host memory.
        /// </summary>
        /// <param name="command_queue">Refers to the command-queue in which the write command will be queued. command_queue and buffer must be created with the same OpenCL context.</param>
        /// <param name="buffer">Refers to a valid buffer object.</param>
        /// <param name="blocking_write">
        ///     Indicates if the write operations are blocking or non-blocking.
        ///     If blocking_write is true, the OpenCL implementation copies the data referred to by ptr and enqueues the write operation in the command-queue. The memory pointed to by ptr can be reused by the application after the clEnqueueWriteBuffer call returns.
        ///     If blocking_write is false, the OpenCL implementation will use ptr to perform a nonblocking write. As the write is non-blocking the implementation can return immediately. The memory pointed to by ptr cannot be reused by the application after the call returns.
        ///     The event argument returns an event object which can be used to query the execution status of the write command. When the write command has completed, the memory pointed to by ptr can then be reused by the application. 
        /// </param>
        /// <param name="offset">The offset in bytes in the buffer object to write to.</param>
        /// <param name="cb">The size in bytes of data being written.</param>
        /// <param name="ptr">The pointer to buffer in host memory where data is to be written from.</param>
        /// <param name="wait_list">
        ///     event_wait_list specifies events that need to complete before this particular command can be executed.
        ///     If event_wait_list is null, then this particular command does not wait on any event to complete. The events specified in event_wait_list act as synchronization points.
        ///     The context associated with events in event_wait_list and command_queue must be the same. 
        /// </param>
        /// <returns>Returns an event object that identifies this particular write command and can be used to query or queue a wait for this particular command to complete.</returns>
        public static cl_event EnqueueWriteBuffer(cl_command_queue command_queue, cl_mem buffer, bool blocking_write, uint offset, uint cb, IntPtr ptr, cl_event[] wait_list = null)
        {
            cl_event_id event_ret;
            uint _blocking_write = (blocking_write ? 1u : 0u);
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);
            
            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    testResult(my.clEnqueueWriteBuffer(command_queue, buffer, _blocking_write, (IntPtr)offset, (IntPtr)cb, ptr, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret));
                }
            }

            return new cl_event(event_ret);
        }

        /// <summary>
        /// Enqueue commands to write to a buffer object from host memory.
        /// </summary>
        /// <param name="command_queue">Refers to the command-queue in which the write command will be queued. command_queue and buffer must be created with the same OpenCL context.</param>
        /// <param name="buffer">Refers to a valid buffer object.</param>
        /// <param name="blocking_write">
        ///     Indicates if the write operations are blocking or non-blocking.
        ///     If blocking_write is true, the OpenCL implementation copies the data referred to by ptr and enqueues the write operation in the command-queue. The memory pointed to by ptr can be reused by the application after the clEnqueueWriteBuffer call returns.
        ///     If blocking_write is false, the OpenCL implementation will use ptr to perform a nonblocking write. As the write is non-blocking the implementation can return immediately. The memory pointed to by ptr cannot be reused by the application after the call returns.
        ///     The event argument returns an event object which can be used to query the execution status of the write command. When the write command has completed, the memory pointed to by ptr can then be reused by the application. 
        /// </param>
        /// <param name="offset">The offset in bytes in the buffer object to write to.</param>
        /// <param name="cb">The size in bytes of data being written.</param>
        /// <param name="ptr">The pointer to buffer in host memory where data is to be written from.</param>
        /// <param name="wait_list">
        ///     event_wait_list specifies events that need to complete before this particular command can be executed.
        ///     If event_wait_list is null, then this particular command does not wait on any event to complete. The events specified in event_wait_list act as synchronization points.
        ///     The context associated with events in event_wait_list and command_queue must be the same. 
        /// </param>
        /// <returns>Returns an event object that identifies this particular write command and can be used to query or queue a wait for this particular command to complete.</returns>
        public static cl_event EnqueueWriteBuffer(cl_command_queue command_queue, cl_mem buffer, bool blocking_write, uint offset, uint cb, Array ptr, cl_event[] wait_list = null)
        {
            cl_event_id event_ret;
            uint _blocking_write = (blocking_write ? 1u : 0u);
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);
            
            GCHandle arrayHandle = GCHandle.Alloc(ptr, GCHandleType.Pinned);
            IntPtr arrayAddress = Marshal.UnsafeAddrOfPinnedArrayElement(ptr, 0);
            cl_status r;

            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    r = my.clEnqueueWriteBuffer(command_queue, buffer, _blocking_write, (IntPtr)offset, (IntPtr)cb, arrayAddress, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret);
                }
            }

            arrayHandle.Free();
            testResult(r);

            return new cl_event(event_ret);
        }

        /// <summary>
        /// Enqueues a command to copy a buffer object to another buffer object.
        /// </summary>
        /// <param name="command_queue">The command-queue in which the copy command will be queued. The OpenCL context associated with command_queue, src_buffer, and dst_buffer must be the same.</param>
        /// <param name="src_buffer">Source buffer</param>
        /// <param name="dest_buffer">Destination buffer</param>
        /// <param name="src_offset">The offset where to begin copying data from src_buffer.</param>
        /// <param name="dest_offset">The offset where to begin copying data into dst_buffer.</param>
        /// <param name="cb">Refers to the size in bytes to copy.</param>
        /// <param name="wait_list">
        ///     event_wait_list specifies events that need to complete before this particular command can be executed.
        ///     If event_wait_list is null, then this particular command does not wait on any event to complete. The events specified in event_wait_list act as synchronization points.
        ///     The context associated with events in event_wait_list and command_queue must be the same. 
        /// </param>
        /// <returns>Returns an event object that identifies this particular copy command and can be used toquery or queue a wait for this particular command to complete.</returns>
        public static cl_event EnqueueCopyBuffer(cl_command_queue command_queue, cl_mem src_buffer, cl_mem dest_buffer, uint src_offset, uint dest_offset, uint cb, cl_event[] wait_list = null)
        {
            cl_event_id event_ret;
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);

            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    testResult(my.clEnqueueCopyBuffer(command_queue, src_buffer, dest_buffer, (IntPtr)src_offset, (IntPtr)dest_offset, (IntPtr)cb, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret));
                }
            }

            return new cl_event(event_ret);
        }

        /// <summary>
        /// Enqueues a command to map a region of the buffer object given by buffer into the host address space and returns a pointer to this mapped region.
        /// </summary>
        /// <param name="command_queue">Must be a valid command-queue.</param>
        /// <param name="buffer">A valid buffer object. The OpenCL context associated with command_queue and buffer must be the same.</param>
        /// <param name="host_ptr">Host address where buffer is mapped into</param>
        /// <param name="blocking_map">
        ///     Indicates if the map operation is blocking or non-blocking.
        ///     If blocking_map is true, clEnqueueMapBuffer does not return until the specified region in buffer can be mapped.
        ///     If blocking_map is false i.e. map operation is non-blocking, the pointer to the mapped region returned by clEnqueueMapBuffer cannot be used until the map command has completed.
        ///     The event argument returns an event object which can be used to query the execution status of the map command. When the map command is completed, the application can access the contents of the mapped region using the pointer returned by clEnqueueMapBuffer.
        /// </param>
        /// <param name="map_flags">Is a bit-field and can be set to CL_MAP_READ to indicate that the region specified by (offset, cb) in the buffer object is being mapped for reading, and/or CL_MAP_WRITE to indicate that the region specified by (offset, cb) in the buffer object is being mapped for writing.</param>
        /// <param name="offset">The offset in bytes of the region in the buffer object that is being mapped.</param>
        /// <param name="cb">The size of the region in the buffer object that is being mapped.</param>
        /// <param name="wait_list">
        ///     event_wait_list specifies events that need to complete before this particular command can be executed.
        ///     If event_wait_list is null, then this particular command does not wait on any event to complete. The events specified in event_wait_list act as synchronization points.
        ///     The context associated with events in event_wait_list and command_queue must be the same. 
        /// </param>
        /// <returns>Returns an event object that identifies this particular copy command and can be used toquery or queue a wait for this particular command to complete.</returns>
        public static cl_event EnqueueMapBuffer(cl_command_queue command_queue, cl_mem buffer, out IntPtr host_ptr, bool blocking_map, cl_map_flags map_flags, uint offset, uint cb, cl_event[] wait_list = null)
        {
            cl_event_id event_ret;
            cl_status r;
            uint _blocking_map = (blocking_map ? 1u : 0u);
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);

            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    host_ptr = my.clEnqueueMapBuffer(command_queue, buffer, _blocking_map, map_flags, (IntPtr)offset, (IntPtr)cb, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret, out r);
                }
            }

            testResult(r);

            return new cl_event(event_ret);
        }

        /// <summary>
        /// Enqueues a command to unmap a previously mapped region of a memory object.
        /// </summary>
        /// <param name="command_queue">Must be a valid command-queue.</param>
        /// <param name="memobj">A valid memory object. The OpenCL context associated with command_queue and memobj must be the same.</param>
        /// <param name="mapped_ptr">The host address returned by a previous call to clEnqueueMapBuffer or clEnqueueMapImage for memobj.</param>
        /// <param name="wait_list">
        ///     Specify events that need to complete before clEnqueueUnmapMemObject can be executed.
        ///     If event_wait_list is null, then clEnqueueUnmapMemObject does not wait on any event to complete. The events specified in event_wait_list act as synchronization points.
        ///     The context associated with events in event_wait_list and command_queue must be the same.
        /// </param>
        /// <returns></returns>
        public static cl_event EnqueueUnmapMemObject(cl_command_queue command_queue, cl_mem memobj, IntPtr mapped_ptr, cl_event[] wait_list = null)
        {
            cl_event_id event_ret;
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);
            
            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    testResult(my.clEnqueueUnmapMemObject(command_queue, memobj, mapped_ptr, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret));
                }
            }

            return new cl_event(event_ret);
        }

        /// <summary>
        /// Decrements the memory object reference count.
        /// </summary>
        /// <param name="memobj">A valid memory object.</param>
        /// After the memobj reference count becomes zero and commands queued for execution on a command-queue(s) that use memobj have finished, the memory object is deleted.
        public static void ReleaseMemObject(cl_mem memobj)
        {
            testResult(my.clReleaseMemObject(memobj));
        }
        #endregion
        // ------------------------------------------
        #region Kernel Management
        private delegate cl_program _clCreateProgramWithSource(cl_context context, uint count, IntPtr source, IntPtr lengths, out cl_status errcode_ret);
        private delegate cl_program _clCreateProgramWithBinary(cl_context context, IntPtr device_list, IntPtr lengths, IntPtr binaries, IntPtr binary_status, out cl_status errcode_ret);
        private delegate cl_status _clBuildProgram(cl_program program, uint num_devices, IntPtr device_list, IntPtr options, cl_build_program_callback pfn_notify, IntPtr user_data);
        private delegate cl_status _clGetProgramInfo(cl_program program, cl_program_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);
        private delegate cl_status _clGetProgramBuildInfo(cl_program program, cl_device_id device, cl_program_build_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);
        private delegate cl_kernel _clCreateKernel(cl_program program, IntPtr kernel_name, out cl_status errcode_ret);
        private delegate cl_status _clCreateKernelsInProgram(cl_program program, uint num_kernels, IntPtr kernels, out uint num_kernels_ret);
        private delegate cl_status _clGetKernelInfo(cl_kernel kernel, cl_kernel_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);
        private delegate cl_status _clGetKernelWorkGroupInfo(cl_kernel kernel, cl_device_id device, cl_kernel_work_group_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);
        private delegate cl_status _clSetKernelArg(cl_kernel kernel, uint arg_index, IntPtr arg_size, IntPtr arg_value);
        private delegate cl_status _clEnqueueNDRangeKernel(cl_command_queue command_queue, cl_kernel kernel, uint work_dim, IntPtr global_work_offset, IntPtr global_work_size, IntPtr local_work_size, uint num_events_in_wait_list, IntPtr wait_list, out cl_event_id event_ret);
        private delegate cl_status _clReleaseKernel(cl_kernel kernel);
        private delegate cl_status _clReleaseProgram(cl_program program);

        _clCreateProgramWithSource clCreateProgramWithSource;
        _clCreateProgramWithBinary clCreateProgramWithBinary;
        _clBuildProgram clBuildProgram;
        _clGetProgramInfo clGetProgramInfo;
        _clGetProgramBuildInfo clGetProgramBuildInfo;
        _clCreateKernel clCreateKernel;
        _clCreateKernelsInProgram clCreateKernelsInProgram;
        _clGetKernelInfo clGetKernelInfo;
        _clGetKernelWorkGroupInfo clGetKernelWorkGroupInfo;
        _clSetKernelArg clSetKernelArg;
        _clEnqueueNDRangeKernel clEnqueueNDRangeKernel;
        _clReleaseKernel clReleaseKernel;
        _clReleaseProgram clReleaseProgram;

        /// <summary>
        /// Creates a program object for a context, and loads the source code specified by the text strings in the strings array into the program object.
        /// </summary>
        /// <param name="context">Must be a valid OpenCL context.</param>
        /// <param name="source">An array of strings that make up the source code.</param>
        /// <returns>Returns the program.</returns>
        public static cl_program CreateProgramWithSource(cl_context context, string[] source)
        {
            int size = source.Length;
            IntPtr[] _source = new IntPtr[size];
            for (int i = 0; i < size; i++)
            {
                _source[i] = Marshal.StringToHGlobalAnsi(source[i]);
            }

            IntPtr _sourcePtr = Marshal.AllocHGlobal(IntPtr.Size * size);
            Marshal.Copy(_source, 0, _sourcePtr, size);
            cl_status r;
            cl_program program = my.clCreateProgramWithSource(context, (uint)size, _sourcePtr, IntPtr.Zero, out r);

            for (int i = 0; i < size; i++)
            {
                Marshal.FreeHGlobal(_source[i]);
            }
            Marshal.FreeHGlobal(_sourcePtr);
            testResult(r);
            return program;
        }

        /// <summary>
        /// Creates a program object for a context, and loads specified binary data into the program object.
        /// </summary>
        /// <param name="context">Must be a valid OpenCL context.</param>
        /// <param name="device_list">A pointer to a list of devices that are in context. device_list must be a non-NULL value. The binaries are loaded for devices specified in this list.</param>
        /// <param name="binaries">An array to program binaries to be loaded for devices specified by device_list. For each device given by device_list[i], the program binary for that device is given as byte[] in binaries[i].</param>
        /// <param name="binary_status">Returns whether the program binary for each device specified in device_list was loaded successfully or not.</param>
        /// <returns>Returns a valid non-zero program object.</returns>
        public static cl_program CreateProgramWithBinary(cl_context context, cl_device_id[] device_list, byte[][] binaries, out cl_status[] binary_status)
        {
            int size = device_list.Length;

            binary_status = new cl_status[size];
            IntPtr[] _binaries = new IntPtr[size];
            IntPtr[] lengths = new IntPtr[size];
            for (int i = 0; i < size; i++)
            {
                _binaries[i] = Marshal.AllocHGlobal(binaries[i].Length);
                Marshal.Copy(binaries[i], 0, _binaries[i], binaries[i].Length);
                lengths[i] = (IntPtr)binaries[i].Length;
            }

            IntPtr _binariesPtr = Marshal.AllocHGlobal(IntPtr.Size * size);
            Marshal.Copy(_binaries, 0, _binariesPtr, size);
            cl_status r;
            cl_program program;
            unsafe
            {
                fixed (cl_device_id* pdevice = device_list)
                {
                    fixed (IntPtr* plengths = lengths)
                    {
                        fixed (cl_status* pstatus = binary_status)
                        {
                            program = my.clCreateProgramWithBinary(context, (IntPtr)pdevice, (IntPtr)plengths, _binariesPtr, (IntPtr)pstatus, out r);
                        }
                    }
                }
            }

            for (int i = 0; i < size; i++)
            {
                Marshal.FreeHGlobal(_binaries[i]);
            }
            Marshal.FreeHGlobal(_binariesPtr);
            testResult(r);
            return program;
        }

        /// <summary>
        /// Builds (compiles and links) a program executable from the program source or binary.
        /// </summary>
        /// <param name="program">The program object.</param>
        /// <param name="devices">
        ///     A list of devices associated with program.
        ///     If device_list is NULL value, the program executable is built for all devices associated with program for which a source or binary has been loaded.
        ///     If device_list is a non-NULL value, the program executable is built for devices specified in this list for which a source or binary has been loaded.
        /// </param>
        /// <param name="options">A string that describes the build options to be used for building the program executable.</param>
        /// <param name="pfn_notify">Callback, see <see cref="cl_build_program_callback"/></param>
        /// <param name="user_data">Passed as an argument when pfn_notify is called. user_data can be NULL.</param>
        public static void BuildProgram(cl_program program, cl_device_id[] devices = null, string options = null, cl_build_program_callback pfn_notify = null, IntPtr user_data = default(IntPtr))
        {
            cl_status r;
            IntPtr _options = options == null ? IntPtr.Zero : Marshal.StringToHGlobalAnsi(options);
            uint num_devices = devices == null ? 0 : (uint)devices.Length;

            unsafe
            {
                fixed (cl_device_id* p = devices)
                {
                    r = my.clBuildProgram(program, num_devices, (IntPtr)p, _options, pfn_notify, user_data);
                }
            }

            if (options != null)
            {
                Marshal.FreeHGlobal(_options);
            }

            testResult(r);
        }

        /// <summary>
        /// Returns information about the program object.
        /// </summary>
        /// <param name="program">Specifies the program object being queried.</param>
        /// <returns>A structure with the information.</returns>
        public static cl_program_info_return GetProgramInfo(cl_program program)
        {
            cl_program_info_return ret = new cl_program_info_return();

            ret.reference_count = GetProgramInfoUInt(program, cl_program_info.CL_PROGRAM_REFERENCE_COUNT);
            ret.context.p = GetProgramInfoIntPtr(program, cl_program_info.CL_PROGRAM_CONTEXT);
            ret.devices = GetProgramInfoDevices(program);
            ret.source = GetProgramInfoString(program, cl_program_info.CL_PROGRAM_SOURCE);
            ret.binaries = GetProgramInfoBinaries(program);

            return ret;
        }

        internal static byte[][] GetProgramInfoBinaries(cl_program program)
        {
            int count = (int)GetProgramInfoUInt(program, cl_program_info.CL_PROGRAM_NUM_DEVICES);

            IntPtr[] lengths = new IntPtr[count];
            IntPtr ignored, size = (IntPtr)(IntPtr.Size * count);
            unsafe
            {
                fixed (IntPtr* plengths = lengths)
                {
                    testResult(my.clGetProgramInfo(program, cl_program_info.CL_PROGRAM_BINARY_SIZES, size, (IntPtr)plengths, out ignored));
                }
            }

            IntPtr[] _binaries = new IntPtr[count];
            for (int i = 0; i < count; i++)
            {
                _binaries[i] = Marshal.AllocHGlobal(lengths[i].ToInt32());
            }
            IntPtr _binariesPtr = Marshal.AllocHGlobal(IntPtr.Size * count);
            Marshal.Copy(_binaries, 0, _binariesPtr, count);
            
            cl_status r = my.clGetProgramInfo(program, cl_program_info.CL_PROGRAM_BINARIES, size, _binariesPtr, out ignored);

            byte[][] binaries = new byte[count][];
            for (int i = 0; i < count; i++)
            {
                int length = lengths[i].ToInt32();
                if (length > 0)
                {
                    binaries[i] = new byte[length];
                    Marshal.Copy(_binaries[i], binaries[i], 0, length);
                }
                Marshal.FreeHGlobal(_binaries[i]);
            }
            Marshal.FreeHGlobal(_binariesPtr);
            testResult(r);
            return binaries;
        }

        internal static cl_device_id[] GetProgramInfoDevices(cl_program program)
        {
            cl_device_id[] devices;

            uint count = GetProgramInfoUInt(program, cl_program_info.CL_PROGRAM_NUM_DEVICES);
            if (count == 0)
            {
                return null;
            }

            devices = new cl_device_id[count];
            IntPtr ignored, size = (IntPtr)(IntPtr.Size * count);
            unsafe
            {
                fixed (cl_device_id* pdevices = devices)
                {
                    testResult(my.clGetProgramInfo(program, cl_program_info.CL_PROGRAM_DEVICES, size, (IntPtr)pdevices, out ignored));
                }
            }

            return devices;
        }

        internal static string GetProgramInfoString(cl_program program, cl_program_info param_name)
        {
            IntPtr size, ignored;
            testResult(my.clGetProgramInfo(program, param_name, IntPtr.Zero, IntPtr.Zero, out size));
            IntPtr mem = Marshal.AllocHGlobal((int)size);
            cl_status r = my.clGetProgramInfo(program, param_name, size, mem, out ignored);
            string ret = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            testResult(r);
            return ret;
        }

        internal static IntPtr GetProgramInfoIntPtr(cl_program program, cl_program_info param_name)
        {
            IntPtr ret, ignored, size = (IntPtr)IntPtr.Size;
            unsafe
            {
                IntPtr* p = &ret;
                testResult(my.clGetProgramInfo(program, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        internal static uint GetProgramInfoUInt(cl_program program, cl_program_info param_name)
        {
            uint ret;
            IntPtr ignored, size = (IntPtr)sizeof(uint);
            unsafe
            {
                uint* p = &ret;
                testResult(my.clGetProgramInfo(program, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        /// <summary>
        /// Returns build information for each device in the program object.
        /// </summary>
        /// <param name="program">Specifies the program object being queried.</param>
        /// <param name="device">Specifies the device for which build information is being queried. device must be a valid device associated with program.</param>
        /// <returns>Returns the information.</returns>
        public static cl_program_build_info_return GetProgramBuildInfo(cl_program program, cl_device_id device)
        {
            cl_program_build_info_return ret = new cl_program_build_info_return();

            ret.log = GetProgramBuildInfoString(program, device, cl_program_build_info.CL_PROGRAM_BUILD_LOG);
            ret.options = GetProgramBuildInfoString(program, device, cl_program_build_info.CL_PROGRAM_BUILD_OPTIONS);
            ret.status = GetProgramBuildStatus(program, device);

            return ret;
        }

        internal static string GetProgramBuildInfoString(cl_program program, cl_device_id device, cl_program_build_info param_name)
        {
            IntPtr size, ignored;
            testResult(my.clGetProgramBuildInfo(program, device, param_name, IntPtr.Zero, IntPtr.Zero, out size));
            IntPtr mem = Marshal.AllocHGlobal((int)size);
            cl_status r = my.clGetProgramBuildInfo(program, device, param_name, size, mem, out ignored);
            string ret = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            testResult(r);
            return ret;
        }

        /// <summary>
        /// Returns the build status of program for a specific device as given by device.
        /// </summary>
        /// <param name="program">Specifies the program object being queried.</param>
        /// <param name="device">Specifies the device for which build information is being queried. device must be a valid device associated with program.</param>
        /// <returns>Returns the build status.</returns>
        public static cl_build_status GetProgramBuildStatus(cl_program program, cl_device_id device)
        {
            int ret;
            IntPtr ignored, size = (IntPtr)sizeof(int);
            unsafe
            {
                int* p = &ret;
                testResult(my.clGetProgramBuildInfo(program, device, cl_program_build_info.CL_PROGRAM_BUILD_STATUS, size, (IntPtr)p, out ignored));
            }
            return (cl_build_status)ret;
        }

        /// <summary>
        /// Creates a kernel object.
        /// </summary>
        /// <param name="program">A program object with a successfully built executable.</param>
        /// <param name="kernel_name">A function name in the program declared with the __kernel qualifier.</param>
        /// <returns>Returns the kernel.</returns>
        public static cl_kernel CreateKernel(cl_program program, string kernel_name)
        {
            IntPtr _kernel_name = Marshal.StringToHGlobalAnsi(kernel_name);
            cl_status r;
            cl_kernel kernel = my.clCreateKernel(program, _kernel_name, out r);
            Marshal.FreeHGlobal(_kernel_name);
            testResult(r);
            return kernel;
        }

        /// <summary>
        /// Creates kernel objects for all kernel functions in a program object.
        /// </summary>
        /// <param name="program">A program object with a successfully built executable.</param>
        /// <returns>Returns an array with kernels or null, if there are no kernels in this program.</returns>
        public static cl_kernel[] CreateKernelsInProgram(cl_program program)
        {
            uint numKernels, ignored;
            testResult(my.clCreateKernelsInProgram(program, 0, IntPtr.Zero, out numKernels));
            if (numKernels > 0)
            {
                cl_kernel[] kernels = new cl_kernel[numKernels];
                unsafe
                {
                    fixed (cl_kernel* p = kernels)
                    {
                        testResult(my.clCreateKernelsInProgram(program, numKernels, (IntPtr)p, out ignored));
                    }
                }
                return kernels;
            }
            else
            {
                return null;
            }
        }

        /// <summary>
        /// Returns information about the kernel object.
        /// </summary>
        /// <param name="kernel">Specifies the kernel object being queried.</param>
        /// <returns>A structure with the information.</returns>
        public static cl_kernel_info_return GetKernelInfo(cl_kernel kernel)
        {
            cl_kernel_info_return ret = new cl_kernel_info_return();

            ret.kernel_func_name = GetKernelInfoString(kernel, cl_kernel_info.CL_KERNEL_FUNCTION_NAME);
            ret.kernel_num_args = GetKernelInfoUInt(kernel, cl_kernel_info.CL_KERNEL_NUM_ARGS);
            ret.reference_count = GetKernelInfoUInt(kernel, cl_kernel_info.CL_KERNEL_REFERENCE_COUNT);
            ret.kernel_context.p = GetKernelInfoIntPtr(kernel, cl_kernel_info.CL_KERNEL_CONTEXT);
            ret.kernel_program.p = GetKernelInfoIntPtr(kernel, cl_kernel_info.CL_KERNEL_PROGRAM);

            return ret;
        }

        internal static string GetKernelInfoString(cl_kernel kernel, cl_kernel_info param_name)
        {
            IntPtr size, ignored;
            testResult(my.clGetKernelInfo(kernel, param_name, IntPtr.Zero, IntPtr.Zero, out size));
            IntPtr mem = Marshal.AllocHGlobal((int)size);
            cl_status r = my.clGetKernelInfo(kernel, param_name, size, mem, out ignored);
            string ret = Marshal.PtrToStringAnsi(mem);
            Marshal.FreeHGlobal(mem);
            testResult(r);
            return ret;
        }

        internal static uint GetKernelInfoUInt(cl_kernel kernel, cl_kernel_info param_name)
        {
            uint ret;
            IntPtr ignored, size = (IntPtr)sizeof(uint);
            unsafe
            {
                uint* p = &ret;
                testResult(my.clGetKernelInfo(kernel, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        internal static IntPtr GetKernelInfoIntPtr(cl_kernel kernel, cl_kernel_info param_name)
        {
            IntPtr ret, ignored, size = (IntPtr)IntPtr.Size;
            unsafe
            {
                IntPtr* p = &ret;
                testResult(my.clGetKernelInfo(kernel, param_name, size, (IntPtr)p, out ignored));
            }
            return ret;
        }

        /// <summary>
        /// Returns information about the kernel object that may be specific to a device. 
        /// </summary>
        /// <param name="kernel">Specifies the kernel object being queried.</param>
        /// <param name="device">Identifies a specific device in the list of devices associated with kernel.</param>
        /// <returns>A structure with kernel information.</returns>
        public static cl_kernel_work_group_info_return GetKernelWorkGroupInfo(cl_kernel kernel, cl_device_id device)
        {
            cl_kernel_work_group_info_return ret = new cl_kernel_work_group_info_return();
            IntPtr ignored;

            unsafe
            {
                UIntPtr val;
                IntPtr p = (IntPtr)(&val), size = (IntPtr)IntPtr.Size;
                testResult(my.clGetKernelWorkGroupInfo(kernel, device, cl_kernel_work_group_info.CL_KERNEL_WORK_GROUP_SIZE, size, p, out ignored));
                ret.work_group_size = val.ToUInt32();
            }

            unsafe
            {
                IntPtr size = (IntPtr)(IntPtr.Size * 3);
                UIntPtr[] array = new UIntPtr[3];
                fixed (UIntPtr* p = array)
                {
                    testResult(my.clGetKernelWorkGroupInfo(kernel, device, cl_kernel_work_group_info.CL_KERNEL_COMPILE_WORK_GROUP_SIZE, size, (IntPtr)p, out ignored));
                }
                ret.compile_work_group_size = new uint[3];
                ret.compile_work_group_size[0] = array[0].ToUInt32();
                ret.compile_work_group_size[1] = array[1].ToUInt32();
                ret.compile_work_group_size[2] = array[2].ToUInt32();
            }

            unsafe
            {
                ulong val;
                IntPtr size = (IntPtr)sizeof(ulong);
                testResult(my.clGetKernelWorkGroupInfo(kernel, device, cl_kernel_work_group_info.CL_KERNEL_LOCAL_MEM_SIZE, size, (IntPtr)(&val), out ignored));
                ret.local_mem_size = val;
            }

#if OPENCL_1_1
            unsafe
            {
                UIntPtr val;
                IntPtr p = (IntPtr)(&val), size = (IntPtr)IntPtr.Size;
                testResult(my.clGetKernelWorkGroupInfo(kernel, device, cl_kernel_work_group_info.CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, size, p, out ignored));
                ret.preferred_work_group_size_multiple = val.ToUInt32();
            }

            unsafe
            {
                ulong val;
                IntPtr size = (IntPtr)sizeof(ulong);
                testResult(my.clGetKernelWorkGroupInfo(kernel, device, cl_kernel_work_group_info.CL_KERNEL_PRIVATE_MEM_SIZE, size, (IntPtr)(&val), out ignored));
                ret.private_mem_size = val;
            }
#endif

            return ret;
        }

        /// <summary>
        /// Used to set the argument value for a specific argument of a kernel. 
        /// </summary>
        /// <param name="kernel">A valid kernel object.</param>
        /// <param name="arg_index">The argument index. Arguments to the kernel are referred by indices that go from 0 for the leftmost argument to n - 1, where n is the total number of arguments declared by a kernel.</param>
        /// <param name="value">A buffer as parameter.</param>
        public static void SetKernelArg(cl_kernel kernel, uint arg_index, cl_mem value)
        {
            unsafe
            {
                cl_mem* p = &value;
                testResult(my.clSetKernelArg(kernel, arg_index, (IntPtr)IntPtr.Size, (IntPtr)p));
            }
        }

        /// <summary>
        /// Used to set the argument value for a specific argument of a kernel. 
        /// </summary>
        /// <param name="kernel">A valid kernel object.</param>
        /// <param name="arg_index">The argument index. Arguments to the kernel are referred by indices that go from 0 for the leftmost argument to n - 1, where n is the total number of arguments declared by a kernel.</param>
        /// <param name="value">An integer as parameter.</param>
        public static void SetKernelArg(cl_kernel kernel, uint arg_index, int value)
        {
            unsafe
            {
                int* p = &value;
                testResult(my.clSetKernelArg(kernel, arg_index, (IntPtr)sizeof(int), (IntPtr)p));
            }
        }

        /// <summary>
        /// Used to set the argument value for a specific argument of a kernel. 
        /// </summary>
        /// <param name="kernel">A valid kernel object.</param>
        /// <param name="arg_index">The argument index. Arguments to the kernel are referred by indices that go from 0 for the leftmost argument to n - 1, where n is the total number of arguments declared by a kernel.</param>
        /// <param name="value">A double as parameter.</param>
        public static void SetKernelArg(cl_kernel kernel, uint arg_index, double value)
        {
            unsafe
            {
                double* p = &value;
                testResult(my.clSetKernelArg(kernel, arg_index, (IntPtr)sizeof(double), (IntPtr)p));
            }
        }

        /// <summary>
        /// Uses to set the size of a local memory argument of a kernel.
        /// </summary>
        /// <param name="kernel">A valid kernel object.</param>
        /// <param name="arg_index">The argument index. Arguments to the kernel are referred by indices that go from 0 for the leftmost argument to n - 1, where n is the total number of arguments declared by a kernel.</param>
        /// <param name="size">Size of the local memory.</param>
        public static void SetKernelArgLocalSize(cl_kernel kernel, uint arg_index, uint size)
        {
            testResult(my.clSetKernelArg(kernel, arg_index, (IntPtr)size, IntPtr.Zero));
        }

        /// <summary>
        /// Enqueues a command to execute a kernel on a device.
        /// </summary>
        /// <param name="command_queue">A valid command-queue. The kernel will be queued for execution on the device associated with command_queue.</param>
        /// <param name="kernel">A valid kernel object. The OpenCL context associated with kernel and command_queue must be the same.</param>
        /// <param name="work_dim">The number of dimensions used to specify the global work-items and work-items in the work-group.</param>
        /// <param name="global_work_size">An array of work_dim unsigned values that describe the number of global work-items in work_dim dimensions that will execute the kernel function. The total number of global work-items is computed as global_work_size[0] *...* global_work_size[work_dim - 1].</param>
        /// <param name="local_work_size">An array of work_dim unsigned values that describe the number of work-items that make up a work-group (also referred to as the size of the work-group) that will execute the kernel specified by kernel. The total number of work-items in a work-group is computed as local_work_size[0] *... * local_work_size[work_dim - 1].</param>
        /// <param name="wait_list">Specify events that need to complete before this particular command can be executed.</param>
        /// <returns>Returns an event object that identifies this particular kernel execution instance.</returns>
        public static cl_event EnqueueNDRangeKernel(cl_command_queue command_queue, cl_kernel kernel, uint work_dim, int[] global_work_size, int[] local_work_size, cl_event[] wait_list = null)
        {
            if (work_dim != global_work_size.Length || work_dim != local_work_size.Length)
            {
                throw new ApplicationException("OpenCL error: Inconsistent work size dimensions.");
            }

            cl_event_id event_ret;
            uint num_events_in_wait_list = wait_list == null ? 0 : (uint)wait_list.Length;
            cl_event_id[] list = cl_event.convertList(wait_list);
            IntPtr[] global = new IntPtr[work_dim];
            IntPtr[] local = new IntPtr[work_dim];
            for (int i = 0; i < work_dim; i++)
            {
                global[i] = (IntPtr)global_work_size[i];
                local[i] = (IntPtr)local_work_size[i];
            }
            
            unsafe
            {
                fixed (cl_event_id* waitAddress = list)
                {
                    fixed (IntPtr* globalAddress = global, localAddress = local)
                    {
                        testResult(my.clEnqueueNDRangeKernel(command_queue, kernel, 1, IntPtr.Zero, (IntPtr)globalAddress, (IntPtr)localAddress, num_events_in_wait_list, (IntPtr)waitAddress, out event_ret));
                    }
                }
            }
            
            return new cl_event(event_ret);
        }

        /// <summary>
        /// Decrements the kernel reference count.
        /// </summary>
        /// <param name="kernel">The kernel to release.</param>
        public static void ReleaseKernel(cl_kernel kernel)
        {
            testResult(my.clReleaseKernel(kernel));
        }

        /// <summary>
        /// Decrements the program reference count.
        /// </summary>
        /// <param name="program">The program to release.</param>
        public static void ReleaseProgram(cl_program program)
        {
            testResult(my.clReleaseProgram(program));
        }
        #endregion
        // ------------------------------------------
        #region Event Management
        private delegate cl_status _clWaitForEvents(uint num_events, IntPtr event_list);
        private delegate cl_status _clRetainEvent(cl_event_id e);
        private delegate cl_status _clReleaseEvent(cl_event_id e);
        private delegate cl_status _clGetEventInfo(cl_event_id e, cl_event_info param_name, IntPtr param_value_size, IntPtr param_value, out IntPtr param_value_size_ret);
#if OPENCL_1_1
        private delegate cl_error _clSetEventCallback(cl_event_id e, cl_command_execution_status command_exec_callback_type, cl_set_event_callback pfn_event_notify, IntPtr user_data);
        private delegate cl_event_id _clCreateUserEvent(cl_context context, out cl_error errcode_ret);
        private delegate cl_error _clSetUserEventStatus(cl_event_id e, cl_command_execution_status execution_status);
#endif

        _clWaitForEvents clWaitForEvents;
        _clRetainEvent clRetainEvent;
        _clReleaseEvent clReleaseEvent;
        _clGetEventInfo clGetEventInfo;
#if OPENCL_1_1
        _clSetEventCallback clSetEventCallback;
        _clCreateUserEvent clCreateUserEvent;
        _clSetUserEventStatus clSetUserEventStatus;
#endif

        /// <summary>
        /// Waits on the host thread for commands identified by event objects to complete.
        /// </summary>
        /// <param name="event_list">The events specified in event_list act as synchronization points.</param>
        public static void WaitForEvents(cl_event[] event_list)
        {
            if (event_list == null)
            {
                return;
            }

            uint num_events_in_event_list = (uint)event_list.Length;
            cl_event_id[] list = cl_event.convertList(event_list);

            unsafe
            {
                fixed (cl_event_id* eventAddress = list)
                {
                    testResult(my.clWaitForEvents(num_events_in_event_list, (IntPtr)eventAddress));
                }
            }
        }

        /// <summary>
        /// Waits on the host thread for command identified by event object to complete.
        /// </summary>
        /// <param name="e">The event specified in e acts as synchronization point.</param>
        public static void WaitForEvent(cl_event e)
        {
            if (e == null)
            {
                return;
            }

            unsafe
            {
                fixed (cl_event_id* eventAddress = &e.id)
                {
                    testResult(my.clWaitForEvents(1, (IntPtr)eventAddress));
                }
            }
        }

        /// <summary>
        /// Increments the event reference count.
        /// </summary>
        /// The OpenCL commands that return an event perform an implicit retain.
        /// <param name="e">Event object being retained.</param>
        public static void RetainEvent(cl_event e)
        {
            testResult(my.clRetainEvent(e.id));
        }

        /// <summary>
        /// Decrements the event reference count.
        /// </summary>
        /// Decrements the event reference count.
        /// The event object is deleted once the reference count becomes zero, the specific command identified by this event has completed (or terminated)
        /// and there are no commands in the command-queues of a context that require a wait for this event to complete.
        /// <param name="e">Event object being released.</param>
        public static void ReleaseEvent(cl_event e)
        {
            testResult(my.clReleaseEvent(e.id));
        }

        /// <summary>
        /// Returns information about the event object.
        /// </summary>
        /// Using clGetEventInfo to determine if a command identified by event has finished execution (i.e. CL_EVENT_COMMAND_EXECUTION_STATUS returns CL_COMPLETE) is not a synchronization point.
        /// There are no guarantees that the memory objects being modified by command associated with event will be visible to other enqueued commands.
        /// <param name="e">Specifies the event object being queried.</param>
        /// <returns>A structure with the information</returns>
        public static cl_event_info_return GetEventInfo(cl_event e)
        {
            cl_event_info_return ret = new cl_event_info_return();

            ret.queue.p = GetEventInfoIntPtr(e, cl_event_info.CL_EVENT_COMMAND_QUEUE);
            ret.type = (cl_command_type)GetEventInfoUInt(e, cl_event_info.CL_EVENT_COMMAND_TYPE);
            ret.status = (cl_status)GetEventInfoUInt(e, cl_event_info.CL_EVENT_COMMAND_EXECUTION_STATUS);
            ret.reference_count = GetEventInfoUInt(e, cl_event_info.CL_EVENT_REFERENCE_COUNT);
#if OPENCL_1_1
            ret.context.p = GetEventInfoIntPtr(e, cl_event_info.CL_EVENT_CONTEXT);
#endif
            
            return ret;
        }

        internal static IntPtr GetEventInfoIntPtr(cl_event e, cl_event_info param_name)
        {
            IntPtr ignored, val, size = (IntPtr)IntPtr.Size;

            unsafe
            {
                IntPtr* pval = &val;
                testResult(my.clGetEventInfo(e.id, param_name, size, (IntPtr)pval, out ignored));
            }

            return val;
        }

        internal static uint GetEventInfoUInt(cl_event e, cl_event_info param_name)
        {
            IntPtr ignored, size = (IntPtr)sizeof(uint);
            uint val;

            unsafe
            {
                uint* pval = &val;
                testResult(my.clGetEventInfo(e.id, param_name, size, (IntPtr)pval, out ignored));
            }

            return val;
        }

#if OPENCL_1_1
        /// <summary>
        /// Registers a user callback function for a specific command execution status.
        /// </summary>
        /// <param name="e">A valid event object.</param>
        /// <param name="command_exec_callback_type">
        ///         The command execution status for which the callback is registered.
        ///         The command execution callback value for which a callback can be registered is CL_COMPLETE.
        ///         There is no guarantee that the callback functions registered for various execution status values for an event will be called in the exact order that the execution status of a command changes.
        ///         The callback function registered for a command_exec_callback_type value of CL_COMPLETE will be called when the command has completed successfully or is abnormally terminated.
        /// </param>
        /// <param name="pfn_event_notify">Callback, see <see cref="cl_set_event_callback"/>.</param>
        /// <param name="user_data">Will be passed as the user_data argument when pfn_notify is called. user_data can be NULL.</param>
        public static void SetEventCallback(cl_event e, cl_command_execution_status command_exec_callback_type, cl_set_event_callback pfn_event_notify, IntPtr user_data = default(IntPtr))
        {
            testResult(my.clSetEventCallback(e.id, command_exec_callback_type, pfn_event_notify, user_data));
        }

        /// <summary>
        /// Creates a user event object.
        /// </summary>
        /// <param name="context">A valid OpenCL context used to create the user event object.</param>
        /// User events allow applications to enqueue commands that wait on a user event to finish before the command is executed by the device.
        /// The execution status of the user event object created is set to CL_SUBMITTED.
        /// <returns></returns>
        public static cl_event CreateUserEvent(cl_context context)
        {
            cl_error r;
            cl_event_id e = my.clCreateUserEvent(context, out r);
            testResult(r);
            return new cl_event(e);
        }

        /// <summary>
        /// Sets the execution status of a user event object.
        /// </summary>
        /// <param name="e">A user event object created using <see cref="clCreateUserEvent"/>.</param>
        /// <param name="execution_status">
        ///     Specifies the new execution status to be set and can be CL_COMPLETE or a negative integer value to indicate an error.
        ///     A negative integer value causes all enqueued commands that wait on this user event to be terminated.
        ///     clSetUserEventStatus can only be called once to change the execution status of event. 
        /// </param>
        public static void SetUserEventStatus(cl_event e, cl_command_execution_status execution_status)
        {
            testResult(my.clSetUserEventStatus(e.id, execution_status));
        }
#endif
        #endregion
        // ------------------------------------------
#pragma warning restore        649
    }
}
