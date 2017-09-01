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

namespace MPI.Wrappers {

    /// <summary>
    /// contains the *_c2f ( C - to - FORTRAN) and *_f2c -- routines for OpenMPI
    /// </summary>
    class OPENMPI_Converter : Utils.DynLibLoader {

        /// <summary>
        /// ctor
        /// </summary>
        public OPENMPI_Converter()
            : base(
                new string[] { "libmpi.so", "/usr/local/opt/open-mpi/lib/libmpi.dylib" },
                new string[2][][],
                new GetNameMangling[] { Utils.DynLibLoader.Identity, FortranMPIdriver.MacOsMangling },
                new PlatformID[] { PlatformID.Unix, PlatformID.MacOSX },
                new int[] { -1, -1 }) {
        }

        /// <summary>
        /// copies the content of an byte-array to an IntPtr - structure
        /// </summary>
        /// <param name="arr">
        /// number of elements must be equal to the size, in bytes of a pointer (4 on 32-bit mode, 8 on 64-bit mode)
        /// </param>
        /// <returns></returns>
        static internal IntPtr ByteArray2IntPtr(byte[] arr) {
            if (arr.Length != IntPtr.Size)
                throw new ArgumentException("wrong length of input array - must match 'IntPtr.Size';");
            IntPtr ret;
            unsafe {
                fixed (byte* pcComm = &arr[0]) {
                    byte* pComm = (byte*)(&ret);
                    for (int i = 0; i < IntPtr.Size; i++)
                        pComm[i] = pcComm[i];
                }
            }
            return ret;
        }

        /// <summary>
        /// copies the content of an  IntPtr - structure to a byte-array
        /// </summary>
        /// <param name="ptr"></param>
        /// <returns></returns>
        static internal byte[] IntPtr2ByteArray(IntPtr ptr) {
            byte[] ret = new byte[IntPtr.Size];
            unsafe {
                fixed (byte* pcComm = &ret[0]) {
                    byte* pComm = (byte*)(&ptr);
                    for (int i = 0; i < IntPtr.Size; i++)
                        pcComm[i] = pComm[i];
                }
            }
            return ret;
        }


#pragma warning disable 649
        public delegate int _MPI_Comm_c2f(IntPtr C_Comm);
        /// <summary> </summary>
        public _MPI_Comm_c2f MPI_Comm_c2f;
#pragma warning restore 649

#pragma warning disable 649
        public delegate IntPtr _MPI_Comm_f2c(MPI_Comm F_Comm);
        /// <summary> </summary>
        public _MPI_Comm_f2c MPI_Comm_f2c;
#pragma warning restore 649

#pragma warning disable 649
        public delegate int _MPI_Op_c2f(IntPtr C_Op);
        /// <summary> </summary>
        public _MPI_Op_c2f MPI_Op_c2f;
#pragma warning restore 649

#pragma warning disable 649
        public delegate IntPtr _MPI_Op_f2c(MPI_Op F_Op);
        /// <summary> </summary>
        public _MPI_Op_f2c MPI_Op_f2c;
#pragma warning restore 649

#pragma warning disable 649
        public delegate int _MPI_Type_c2f(IntPtr C_Datatype);
        /// <summary> </summary>
        public _MPI_Type_c2f MPI_Type_c2f;
#pragma warning restore 649

#pragma warning disable 649
        public delegate IntPtr _MPI_Type_f2c(MPI_Datatype F_Datatype);
        /// <summary> </summary>
        public _MPI_Type_f2c MPI_Type_f2c;
#pragma warning restore 649

#pragma warning disable 649
        public delegate int _MPI_Request_c2f(IntPtr C_Datatype);
        /// <summary> </summary>
        public _MPI_Request_c2f MPI_Request_c2f;
#pragma warning restore 649

#pragma warning disable 649
        public delegate IntPtr _MPI_Request_f2c(MPI_Request F_Datatype);
        /// <summary> </summary>
        public _MPI_Request_f2c MPI_Request_f2c;
#pragma warning restore 649

    }

}
