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
using System.Runtime.InteropServices;

namespace MPI.Wrappers {

    /// <summary>
    /// predefined communicators and conversion between FORTRAN and C
    /// Communicator handles
    /// </summary>
    public interface IMPI_CommConstants {

        /// <summary>
        /// MPI_COMM_WORLD
        /// </summary>
        MPI_Comm WORLD {
            get;
        }

        /// <summary>
        /// MPI_COMM_SELF
        /// </summary>
        MPI_Comm SELF {
            get;
        }

        /// <summary>
        /// MPI_COMM_NULL
        /// </summary>
        MPI_Comm NULL {
            get;
        }

        /// <summary>
        /// gets the size, in bytes of the C MPI Communicator handle
        /// </summary>
        int GetSizeof_C_MPI_comm();

        /// <summary>
        /// gets the size, in bytes of the FORTRAN MPI Communicator handle
        /// </summary>
        int Get_Sizeof_F77_MPI_COMM();

        /// <summary>
        /// conversion of an C-MPI communicator handle (which can be 4 or 8
        /// byte, depending on  the MPI implementation (MPICH or OpenMPI) and
        /// the 'bitness' (32 Bit or 64 Bit) of the platform) to a standard
        /// 4-byte <see cref="MPI_Comm"/>
        /// </summary>
        MPI_Comm Comm_c2f(byte[] C_MPI_Comm);

        /// <summary>
        /// inverse function of <see cref="Comm_c2f"/>
        /// </summary>
        byte[] Comm_f2c(MPI_Comm F_Comm);
    }

    class MPICH_MPI_Comm : IMPI_CommConstants {

        /// <summary>
        /// Predefined communicator containing all of the MPI processes. See
        /// <see cref="IMPI_CommConstants.WORLD"/>.
        /// </summary>
        public MPI_Comm WORLD {
            get {
                MPI_Comm ret;
                ret.m1 = 0x44000000;
                return ret;
            }
        }

        /// <summary>
        /// Predefined communicator containing only the calling process. See 
        /// <see cref="IMPI_CommConstants.SELF"/>.
        /// </summary>
        public MPI_Comm SELF {
            get {
                MPI_Comm ret;
                ret.m1 = 0x44000001;
                return ret;
            }
        }

        /// <summary>
        /// Predefined communicator representing "no communicator". In the
        /// higher-level interface, this is represented by a <c>null</c>
        /// <see cref="MPI_Comm"/> object.
        /// </summary>
        public MPI_Comm NULL {
            get {
                MPI_Comm ret;
                ret.m1 = 0x04000000;
                return ret;
            }
        }

        /// <summary>
        /// see <see cref="IMPI_CommConstants.GetSizeof_C_MPI_comm"/>
        /// </summary>
        public int GetSizeof_C_MPI_comm() {
            return 4;
        }

        /// <summary>
        /// see <see cref="IMPI_CommConstants.Get_Sizeof_F77_MPI_COMM"/>
        /// </summary>
        public int Get_Sizeof_F77_MPI_COMM() {
            return 4;
        }

        /// <summary>
        /// 
        /// </summary>
        public MPI_Comm Comm_c2f(byte[] C_MPI_Comm) {
            if (C_MPI_Comm.Length != GetSizeof_C_MPI_comm())
                throw new ArgumentException("wrong number of bytes; length of argument must match the value that is returend by 'GetSizeof_C_MPI_comm()';");
            MPI_Comm F_Comm = default(MPI_Comm);
            unsafe {
                fixed (byte* pcComm = &C_MPI_Comm[0]) {
                    byte* pComm = (byte*)&(F_Comm.m1);
                    for (int i = 0; i < C_MPI_Comm.Length; i++) {
                        *pComm = pcComm[i];
                        pComm++;
                    }
                }
            }
            return F_Comm;
        }

        /// <summary>
        /// 
        /// </summary>
        public byte[] Comm_f2c(MPI_Comm F_Comm) {
            byte[] ret = new byte[this.GetSizeof_C_MPI_comm()];
            unsafe {
                byte* pComm = (byte*)&(F_Comm.m1);
                for (int i = 0; i < ret.Length; i++) {
                    ret[i] = *pComm;
                    pComm++;
                }
            }
            return ret;
        }
    }


    class OpenMPI_MPI_Comm : IMPI_CommConstants {

        OPENMPI_Converter m_conv;

        public OpenMPI_MPI_Comm(OPENMPI_Converter conv) {
            m_conv = conv;
        }

        bool deleyedInitdone = false;

        private void DeleaedInit() {
            if (deleyedInitdone)
                return;
            deleyedInitdone = true;
            LoadSymb("ompi_mpi_comm_world", ref m_WORLD);
            LoadSymb("ompi_mpi_comm_self", ref m_SELF);
            LoadSymb("ompi_mpi_comm_null", ref m_NULL);
        }

        private void LoadSymb(string name, ref MPI_Comm sym) {
            string errstr;
            IntPtr addr = Utils.DynamicLibraries.LoadSymbol(m_conv.LibHandle, name, out errstr);
            if (addr == IntPtr.Zero)
                throw new ApplicationException("OpenMPI error: unable to load symbol '" + name + "' from library '" + m_conv.CurrentLibraryName + "', Error string: >" + errstr + "<;");

            sym.m1 = m_conv.MPI_Comm_c2f(addr);
            //Console.WriteLine("val of '" + name + "' is: " + addr + ", fortran value is " + sym.m1);
            IntPtr test = m_conv.MPI_Comm_f2c(sym);
            if (test != addr)
                throw new ApplicationException("fuck");
            //else
            //    Console.WriteLine("back passed.");
        }


        #region IMPI_Comm Members

        MPI_Comm m_WORLD;
        public MPI_Comm WORLD {
            get {
                DeleaedInit();
                return m_WORLD;
            }
        }

        MPI_Comm m_SELF;
        public MPI_Comm SELF {
            get {
                DeleaedInit();
                return m_SELF;
            }
        }

        MPI_Comm m_NULL;
        public MPI_Comm NULL {
            get {
                DeleaedInit();
                return m_NULL;
            }
        }

        /// <summary>
        /// the size of a pointer
        /// </summary>
        public int GetSizeof_C_MPI_comm() {
            return IntPtr.Size;
        }

        /// <summary>
        /// the ((size of the FORTRAN type INTEGER) == (sizoof(int) in C)) - hopefully 4 bytes;
        /// </summary>
        /// <returns></returns>
        public int Get_Sizeof_F77_MPI_COMM() {
            return 4;
        }

        public MPI_Comm Comm_c2f(byte[] C_MPI_Comm) {
            IntPtr _c_MPI_Comm = OPENMPI_Converter.ByteArray2IntPtr(C_MPI_Comm);
            MPI_Comm ret;
            ret.m1 = m_conv.MPI_Comm_c2f(_c_MPI_Comm);
            return ret;
        }

        public byte[] Comm_f2c(MPI_Comm F_Comm) {
            IntPtr ret = m_conv.MPI_Comm_f2c(F_Comm);
            return OPENMPI_Converter.IntPtr2ByteArray(ret);
        }

        #endregion
    }


    class Platform_Native_Comm : IMPI_CommConstants {


        [DllImport("Platform_Native")]
        static extern int BoSSS_Get_MPI_COMM_WORLD();

        public MPI_Comm WORLD {
            get {
                MPI_Comm ret = default(MPI_Comm);
                ret.m1 = BoSSS_Get_MPI_COMM_WORLD();
                return ret;
            }
        }

        [DllImport("Platform_Native")]
        static extern int BoSSS_Get_MPI_COMM_SELF();


        public MPI_Comm SELF {
            get {
                MPI_Comm ret = default(MPI_Comm);
                ret.m1 = BoSSS_Get_MPI_COMM_SELF();
                return ret;
            }
        }

        [DllImport("Platform_Native")]
        static extern int BoSSS_Get_MPI_COMM_NULL();


        public MPI_Comm NULL {
            get {
                MPI_Comm ret = default(MPI_Comm);
                ret.m1 = BoSSS_Get_MPI_COMM_NULL();
                return ret;
            }
        }

        public int GetSizeof_C_MPI_comm() {
            throw new NotImplementedException();
        }

        public int Get_Sizeof_F77_MPI_COMM() {
            throw new NotImplementedException();
        }

        public MPI_Comm Comm_c2f(byte[] C_MPI_Comm) {
            throw new NotImplementedException();
        }

        public byte[] Comm_f2c(MPI_Comm F_Comm) {
            throw new NotImplementedException();
        }
    }

}
