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
    /// predefined MPI operations
    /// </summary>
    public interface IMPI_OpConstants {

        /// <summary> </summary>
        MPI_Op MAX {
            get;
        }
        /// <summary> </summary>
        MPI_Op MIN {
            get;
        }
        /// <summary> </summary>
        MPI_Op SUM {
            get;
        }
        /// <summary> </summary>
        MPI_Op PROD {
            get;
        }
        /// <summary> </summary>
        MPI_Op LAND {
            get;
        }
        /// <summary> </summary>
        MPI_Op BAND {
            get;
        }
        /// <summary> </summary>
        MPI_Op LOR {
            get;
        }
        /// <summary> </summary>
        MPI_Op BOR {
            get;
        }
        /// <summary> </summary>
        MPI_Op LXOR {
            get;
        }
        /// <summary> </summary>
        MPI_Op BXOR {
            get;
        }
        /// <summary> </summary>
        MPI_Op MINLOC {
            get;
        }
        /// <summary> </summary>
        MPI_Op MAXLOC {
            get;
        }
        /// <summary> </summary>
        MPI_Op REPLACE {
            get;
        }

        /// <summary>
        /// conversion of an C-MPI operation handle (which can be 4 or 8 byte, depending on 
        /// the MPI implementation (MPICH or OpenMPI) and the Bitness (32 Bit or 64 Bit) of the platform)
        /// to a standard 4-byte <see cref="MPI_Comm"/>;<br/>
        /// </summary>
        MPI_Op Op_c2f(byte[] C_MPI_Comm);

        /// <summary>
        /// inverse function of <see cref="Op_c2f"/>
        /// </summary>
        byte[] Op_f2c(MPI_Op F_Comm);

        /// <summary>
        /// sizeof C-MPI_Op handele on current platform, in bytes
        /// </summary>
        int GetSizeof_C_MPI_Op();

        /// <summary>
        /// sizeof FORTRAN-MPI_Op handele on current platform, in bytes
        /// </summary>
        int GetSizeof_F_MPI_Op();
    }

    /// <summary>
    /// MPICH operation codes;
    /// equal for MPICH and MS-MPI;
    /// </summary>
    class MPICH_MPI_op : IMPI_OpConstants {


        #region IMPI_Op Members

        public MPI_Op MAX {
            get {
                MPI_Op r;
                r.m1 = 0x58000001;
                return r;
            }
        }

        public MPI_Op MIN {
            get {
                MPI_Op r;
                r.m1 = 0x58000002;
                return r;
            }
        }

        public MPI_Op SUM {
            get {
                MPI_Op r;
                r.m1 = 0x58000003;
                return r;
            }
        }

        public MPI_Op PROD {
            get {
                MPI_Op r;
                r.m1 = 0x58000004;
                return r;
            }
        }

        public MPI_Op LAND {
            get {
                MPI_Op r;
                r.m1 = 0x58000005;
                return r;
            }
        }

        public MPI_Op BAND {
            get {
                MPI_Op r;
                r.m1 = 0x58000006;
                return r;
            }
        }

        public MPI_Op LOR {
            get {
                MPI_Op r;
                r.m1 = 0x58000007;
                return r;
            }
        }

        public MPI_Op BOR {
            get {
                MPI_Op r;
                r.m1 = 0x58000008;
                return r;
            }
        }

        public MPI_Op LXOR {
            get {
                MPI_Op r;
                r.m1 = 0x58000009;
                return r;
            }
        }

        public MPI_Op BXOR {
            get {
                MPI_Op r;
                r.m1 = 0x5800000a;
                return r;
            }
        }

        public MPI_Op MINLOC {
            get {
                MPI_Op r;
                r.m1 = 0x5800000b;
                return r;
            }
        }

        public MPI_Op MAXLOC {
            get {
                MPI_Op r;
                r.m1 = 0x5800000c;
                return r;
            }
        }

        public MPI_Op REPLACE {
            get {
                MPI_Op r;
                r.m1 = 0x5800000d;
                return r;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public MPI_Op Op_c2f(byte[] C_MPI_Comm) {
            if (C_MPI_Comm.Length != GetSizeof_C_MPI_Op())
                throw new ArgumentException("wrong number of bytes; length of argument must match the value that is returend by 'GetSizeof_C_MPI_comm()';");
            MPI_Op F_Comm = default(MPI_Op);
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
        public byte[] Op_f2c(MPI_Op F_Comm) {
            byte[] ret = new byte[this.GetSizeof_C_MPI_Op()];
            unsafe {
                byte* pComm = (byte*)&(F_Comm.m1);
                for (int i = 0; i < ret.Length; i++) {
                    ret[i] = *pComm;
                    pComm++;
                }
            }
            return ret;
        }

        public int GetSizeof_C_MPI_Op() {
            return sizeof(int);
        }

        public int GetSizeof_F_MPI_Op() {
            return sizeof(int);
        }

        #endregion
    }



    class OpenMPI_MPI_Op : IMPI_OpConstants {

        OPENMPI_Converter m_conv;

        public OpenMPI_MPI_Op(OPENMPI_Converter conv) {
            m_conv = conv;
        }

        bool deleyedInitdone = false;

        private void DeleaedInit() {
            if (deleyedInitdone)
                return;
            deleyedInitdone = true;
            LoadSymb("ompi_mpi_op_max", ref m_MPI_MAX);
            LoadSymb("ompi_mpi_op_min", ref m_MPI_MIN);
            LoadSymb("ompi_mpi_op_sum", ref m_MPI_SUM);
            LoadSymb("ompi_mpi_op_prod", ref m_MPI_PROD);
            LoadSymb("ompi_mpi_op_land", ref m_MPI_LAND);
            LoadSymb("ompi_mpi_op_band", ref m_MPI_BAND);
            LoadSymb("ompi_mpi_op_lor", ref m_MPI_LOR);
            LoadSymb("ompi_mpi_op_bor", ref m_MPI_BOR);
            LoadSymb("ompi_mpi_op_lxor", ref m_MPI_LXOR);
            LoadSymb("ompi_mpi_op_bxor", ref m_MPI_BXOR);
            LoadSymb("ompi_mpi_op_minloc", ref m_MPI_MINLOC);
            LoadSymb("ompi_mpi_op_maxloc", ref m_MPI_MAXLOC);
            LoadSymb("ompi_mpi_op_replace", ref m_MPI_REPLACE);
        }

        private void LoadSymb(string name, ref MPI_Op sym) {

            string errstr;
            IntPtr addr = Utils.DynamicLibraries.LoadSymbol(m_conv.LibHandle, name, out errstr);
            if (addr == IntPtr.Zero)
                throw new ApplicationException("OpenMPI error: unable to load symbol '" + name + "' from library '" + m_conv.CurrentLibraryName + "', Error string: >" + errstr + "<;");
            //Console.WriteLine("val of '" + name + "' is: " + addr);
            sym.m1 = m_conv.MPI_Op_c2f(addr);
        }

        MPI_Op m_MPI_MAX;
        MPI_Op m_MPI_MIN;
        MPI_Op m_MPI_SUM;
        MPI_Op m_MPI_PROD;
        MPI_Op m_MPI_LAND;
        MPI_Op m_MPI_BAND;
        MPI_Op m_MPI_LOR;
        MPI_Op m_MPI_BOR;
        MPI_Op m_MPI_LXOR;
        MPI_Op m_MPI_BXOR;
        MPI_Op m_MPI_MINLOC;
        MPI_Op m_MPI_MAXLOC;
        MPI_Op m_MPI_REPLACE;

        #region IMPI_Op Members


        public MPI_Op MAX {
            get {
                DeleaedInit();
                return m_MPI_MAX;
            }
        }

        public MPI_Op MIN {
            get {
                DeleaedInit();
                return m_MPI_MIN;
            }
        }

        public MPI_Op SUM {
            get {
                DeleaedInit();
                return m_MPI_SUM;
            }
        }

        public MPI_Op PROD {
            get {
                DeleaedInit();
                return m_MPI_PROD;
            }
        }

        public MPI_Op LAND {
            get {
                DeleaedInit();
                return m_MPI_LAND;
            }
        }

        public MPI_Op BAND {
            get {
                DeleaedInit();
                return m_MPI_BAND;
            }
        }

        public MPI_Op LOR {
            get {
                DeleaedInit();
                return m_MPI_LOR;
            }
        }

        public MPI_Op BOR {
            get {
                DeleaedInit();
                return m_MPI_BOR;
            }
        }

        public MPI_Op LXOR {
            get {
                DeleaedInit();
                return m_MPI_LXOR;
            }
        }

        public MPI_Op BXOR {
            get {
                DeleaedInit();
                return m_MPI_BXOR;
            }
        }

        public MPI_Op MINLOC {
            get {
                DeleaedInit();
                return m_MPI_MINLOC;
            }
        }

        public MPI_Op MAXLOC {
            get {
                DeleaedInit();
                return m_MPI_MAXLOC;
            }
        }

        public MPI_Op REPLACE {
            get {
                DeleaedInit();
                return m_MPI_REPLACE;
            }
        }

        public int GetSizeof_C_MPI_Op() {
            return IntPtr.Size;
        }

        public int GetSizeof_F_MPI_Op() {
            return sizeof(int);
        }

        public MPI_Op Op_c2f(byte[] C_MPI_Op) {
            IntPtr _c_op = OPENMPI_Converter.ByteArray2IntPtr(C_MPI_Op);
            MPI_Op ret;
            ret.m1 = m_conv.MPI_Op_c2f(_c_op);
            return ret;
        }

        public byte[] Op_f2c(MPI_Op F_Op) {
            IntPtr c_op = m_conv.MPI_Op_f2c(F_Op);
            return OPENMPI_Converter.IntPtr2ByteArray(c_op);
        }

        #endregion
    }

    class Platform_Native_MPI_Op : IMPI_OpConstants {

        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_MAX();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_MIN();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_SUM();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_PROD();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_LAND();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_BAND();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_LOR();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_BOR();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_LXOR();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_BXOR();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_MINLOC();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_MAXLOC();


        public MPI_Op MAX {
            get {
                return new MPI_Op(BoSSS_Get_MPI_MAX());
            }
        }

        public MPI_Op MIN {
            get {
                return new MPI_Op(BoSSS_Get_MPI_MIN());
            }
        }

        public MPI_Op SUM {
            get {
                return new MPI_Op(BoSSS_Get_MPI_SUM());
            }
        }

        public MPI_Op PROD {
            get {
                return new MPI_Op(BoSSS_Get_MPI_MAX());
            }
        }

        public MPI_Op LAND {
            get {
                return new MPI_Op(BoSSS_Get_MPI_LAND());
            }
        }

        public MPI_Op BAND {
            get {
                return new MPI_Op(BoSSS_Get_MPI_BAND());
            }
        }

        public MPI_Op LOR {
            get {
                return new MPI_Op(BoSSS_Get_MPI_LOR());
            }
        }

        public MPI_Op BOR {
            get {
                return new MPI_Op(BoSSS_Get_MPI_BOR());
            }
        }

        public MPI_Op LXOR {
            get {
                return new MPI_Op(BoSSS_Get_MPI_LXOR());
            }
        }

        public MPI_Op BXOR {
            get {
                return new MPI_Op(BoSSS_Get_MPI_BXOR());
            }
        }

        public MPI_Op MINLOC {
            get {
                return new MPI_Op(BoSSS_Get_MPI_MINLOC());
            }
        }

        public MPI_Op MAXLOC {
            get {
                return new MPI_Op(BoSSS_Get_MPI_MAXLOC());
            }
        }

        public MPI_Op REPLACE {
            get {
                throw new NotImplementedException();
                //return new MPI_Op(BoSSS_Get_MPI_RE());
            }
        }

        public MPI_Op Op_c2f(byte[] C_MPI_Comm) {
            throw new NotImplementedException();
        }

        public byte[] Op_f2c(MPI_Op F_Comm) {
            throw new NotImplementedException();
        }

        public int GetSizeof_C_MPI_Op() {
            throw new NotImplementedException();
        }

        public int GetSizeof_F_MPI_Op() {
            throw new NotImplementedException();
        }
    }

}
