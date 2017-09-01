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

namespace ilPSP.LinSolvers.monkey.CL
{
    /// <summary>
    /// OpenCL enrivonment
    /// </summary>
    public class clEnvironment
    {
        internal MPIEnviroment MpiEnv;
        internal cl_device_id device;
        internal cl_context context;

        /// <summary>
        /// Manufacturer of OpenCL Device
        /// </summary>
        public enum Vendor {
            /// <summary>
            /// Advanced Micro Devices
            /// </summary>
            AMD,

            /// <summary>
            /// NVIDIA
            /// </summary>
            NVIDIA,

            /// <summary>
            /// Other or unkown OpenCL Vendor
            /// </summary>
            Other
        }

        /// <summary>
        /// see <see cref="Vendor"/>;
        /// </summary>
        public Vendor _vendor {
            get;
            private set;
        }

        /// <summary>
        /// Creates OpenCL environment on this process
        /// </summary>
        /// <param name="_MpiEnv">MPI environment of this process</param>
        public clEnvironment(MPIEnviroment _MpiEnv)
        {
            MpiEnv = _MpiEnv;

            cl_platform_id[] platforms = cl.GetPlatformIDs();
            cl_platform_info_return pinfo = cl.GetPlatformInfo(platforms[0]);
            string vnd = pinfo.vendor;
            switch (vnd)
            {
                case "NVIDIA Corporation":
                    _vendor = Vendor.NVIDIA;
                    break;
                case "Advanced Micro Devices, Inc.":
                    _vendor = Vendor.AMD;
                    break;
                default:
                    _vendor = Vendor.Other;
                    break;
            }
            
            Console.WriteLine(vnd + " : " + _vendor.ToString());
            Console.WriteLine(pinfo.version);

            if (_vendor == Vendor.AMD)
            {
                clVectorSource.source[0] = clVectorSource.source[0].Replace("cl_khr_fp64", "cl_amd_fp64");
                clMatrixSource.source[0] = clMatrixSource.source[0].Replace("cl_khr_fp64", "cl_amd_fp64");
                Console.WriteLine("Adjusted AMD fp64 extenstion");
            }
            
            cl_device_id[] devices = cl.GetDeviceIDs(platforms[0], cl_device_type.CL_DEVICE_TYPE_GPU);
            int numDevices = devices.Length;

            if (numDevices < _MpiEnv.ProcessesOnMySMP)
                throw new ApplicationException("not enougth OpenCL devices; There must be at least one OpenCL device for each MPI process;");

            device = devices[0];
            context = cl.CreateContext(platforms[0], devices);

            //cl_device_info_return dinfo = cl.GetDeviceInfo(device);
            //Console.WriteLine("Max work group size: " + dinfo.max_work_group_size);
            //Console.WriteLine("Process " + _MpiEnv.ProcessRankOnSMP + " running on device " + dinfo.name + ", " + dinfo.version);
            //Console.WriteLine(dinfo.extensions);
        }

        /// <summary>
        /// Destructor
        /// </summary>
        ~clEnvironment()
        {
            cl.ReleaseContext(context);
        }
    }
}
