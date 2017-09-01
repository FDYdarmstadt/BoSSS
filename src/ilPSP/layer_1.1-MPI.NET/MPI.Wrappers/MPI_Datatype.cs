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
    /// predefined MPI data-type handles
    /// </summary>
    public interface IMPI_DatatypeConstants {
        /// <summary> </summary>
        MPI_Datatype CHAR {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype SIGNED_CHAR {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype UNSIGNED_CHAR {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype BYTE {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype WCHAR {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype SHORT {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype UNSIGNED_SHORT {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype INT {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype UNSIGNED {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype LONG {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype UNSIGNED_LONG {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype FLOAT {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype DOUBLE {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype LONG_DOUBLE {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype LONG_LONG_INT {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype UNSIGNED_LONG_LONG {
            get;
        }
        /// <summary> </summary>
        MPI_Datatype LONG_LONG {
            get;
        }
    }

    class MPICH_MPI_Datatype : IMPI_DatatypeConstants {
        public MPI_Datatype CHAR {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000101;
                return r;
            }
        }
        public MPI_Datatype SIGNED_CHAR {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000118;
                return r;
            }
        }
        public MPI_Datatype UNSIGNED_CHAR {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000102;
                return r;
            }
        }
        public MPI_Datatype BYTE {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c00010d;
                return r;
            }
        }
        public MPI_Datatype WCHAR {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c00020e;
                return r;
            }
        }
        public MPI_Datatype SHORT {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000203;
                return r;
            }
        }
        public MPI_Datatype UNSIGNED_SHORT {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000204;
                return r;
            }
        }
        public MPI_Datatype INT {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000405;
                return r;
            }
        }
        public MPI_Datatype UNSIGNED {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000406;
                return r;
            }
        }
        public MPI_Datatype LONG {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000407;
                return r;
            }
        }
        public MPI_Datatype UNSIGNED_LONG {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000408;
                return r;
            }
        }
        public MPI_Datatype FLOAT {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c00040a;
                return r;
            }
        }
        public MPI_Datatype DOUBLE {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c00080b;
                return r;
            }
        }
        public MPI_Datatype LONG_DOUBLE {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c00080c;
                return r;
            }
        }
        public MPI_Datatype LONG_LONG_INT {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000809;
                return r;
            }
        }
        public MPI_Datatype UNSIGNED_LONG_LONG {
            get {
                MPI_Datatype r;
                r.m1 = 0x4c000819;
                return r;
            }
        }
        public MPI_Datatype LONG_LONG {
            get {
                return LONG_LONG_INT;
            }
        }
    }


    class OpenMPI_MPI_Datatype : IMPI_DatatypeConstants {
        OPENMPI_Converter m_conv;

        public OpenMPI_MPI_Datatype(OPENMPI_Converter conv) {
            m_conv = conv;
        }

        bool deleyedInitdone = false;

        private void DeleaedInit() {
            if (deleyedInitdone)
                return;
            deleyedInitdone = true;
            LoadSymb("ompi_mpi_char", ref m_MPI_CHAR);
            LoadSymb("ompi_mpi_signed_char", ref m_MPI_SIGNED_CHAR);
            LoadSymb("ompi_mpi_unsigned_char", ref m_MPI_UNSIGNED_CHAR);
            LoadSymb("ompi_mpi_byte", ref m_MPI_BYTE);
            LoadSymb("ompi_mpi_wchar", ref m_MPI_WCHAR);
            LoadSymb("ompi_mpi_short", ref m_MPI_SHORT);
            LoadSymb("ompi_mpi_unsigned_short", ref m_MPI_UNSIGNED_SHORT);
            LoadSymb("ompi_mpi_int", ref m_MPI_INT);
            LoadSymb("ompi_mpi_unsigned", ref m_MPI_UNSIGNED);
            LoadSymb("ompi_mpi_long", ref m_MPI_LONG);
            LoadSymb("ompi_mpi_unsigned_long", ref m_MPI_UNSIGNED_LONG);
            LoadSymb("ompi_mpi_float", ref m_MPI_FLOAT);
            LoadSymb("ompi_mpi_double", ref m_MPI_DOUBLE);
            LoadSymb("ompi_mpi_long_double", ref m_MPI_LONG_DOUBLE);
            LoadSymb("ompi_mpi_long_long_int", ref m_MPI_LONG_LONG_INT);
            LoadSymb("ompi_mpi_unsigned_long_long", ref m_MPI_UNSIGNED_LONG_LONG);
        }

        private void LoadSymb(string name, ref MPI_Datatype sym) {
            string errstr;
            IntPtr addr = Utils.DynamicLibraries.LoadSymbol(m_conv.LibHandle, name, out errstr);
            if (addr == IntPtr.Zero)
                throw new ApplicationException("OpenMPI error: unable to load symbol '" + name + "' from library '" + m_conv.CurrentLibraryName + "', Error string: >" + errstr + "<;");
            //Console.WriteLine("val of '" + name + "' is: " + addr);
            sym.m1 = m_conv.MPI_Type_c2f(addr);
        }

        MPI_Datatype m_MPI_CHAR;
        MPI_Datatype m_MPI_SIGNED_CHAR;
        MPI_Datatype m_MPI_UNSIGNED_CHAR;
        MPI_Datatype m_MPI_BYTE;
        MPI_Datatype m_MPI_WCHAR;
        MPI_Datatype m_MPI_SHORT;
        MPI_Datatype m_MPI_UNSIGNED_SHORT;
        MPI_Datatype m_MPI_INT;
        MPI_Datatype m_MPI_UNSIGNED;
        MPI_Datatype m_MPI_LONG;
        MPI_Datatype m_MPI_UNSIGNED_LONG;
        MPI_Datatype m_MPI_FLOAT;
        MPI_Datatype m_MPI_DOUBLE;
        MPI_Datatype m_MPI_LONG_DOUBLE;
        MPI_Datatype m_MPI_LONG_LONG_INT;
        MPI_Datatype m_MPI_UNSIGNED_LONG_LONG;


        public MPI_Datatype CHAR {
            get {
                DeleaedInit();
                return m_MPI_CHAR;
            }
        }
        public MPI_Datatype SIGNED_CHAR {
            get {
                DeleaedInit();
                return m_MPI_SIGNED_CHAR;
            }
        }
        public MPI_Datatype UNSIGNED_CHAR {
            get {
                DeleaedInit();
                return m_MPI_UNSIGNED_CHAR;
            }
        }
        public MPI_Datatype BYTE {
            get {
                DeleaedInit();
                return m_MPI_BYTE;
            }
        }
        public MPI_Datatype WCHAR {
            get {
                DeleaedInit();
                return m_MPI_WCHAR;
            }
        }
        public MPI_Datatype SHORT {
            get {
                DeleaedInit();
                return m_MPI_SHORT;
            }
        }
        public MPI_Datatype UNSIGNED_SHORT {
            get {
                DeleaedInit();
                return m_MPI_UNSIGNED_SHORT;
            }
        }
        public MPI_Datatype INT {
            get {
                DeleaedInit();
                return m_MPI_INT;
            }
        }
        public MPI_Datatype UNSIGNED {
            get {
                DeleaedInit();
                return m_MPI_UNSIGNED;
            }
        }
        public MPI_Datatype LONG {
            get {
                DeleaedInit();
                return m_MPI_LONG;
            }
        }
        public MPI_Datatype UNSIGNED_LONG {
            get {
                DeleaedInit();
                return m_MPI_UNSIGNED_LONG;
            }
        }
        public MPI_Datatype FLOAT {
            get {
                DeleaedInit();
                return m_MPI_FLOAT;
            }
        }
        public MPI_Datatype DOUBLE {
            get {
                DeleaedInit();
                return m_MPI_DOUBLE;
            }
        }
        public MPI_Datatype LONG_DOUBLE {
            get {
                DeleaedInit();
                return m_MPI_LONG_DOUBLE;
            }
        }
        public MPI_Datatype LONG_LONG_INT {
            get {
                DeleaedInit();
                return m_MPI_LONG_LONG_INT;
            }
        }
        public MPI_Datatype UNSIGNED_LONG_LONG {
            get {
                DeleaedInit();
                return m_MPI_UNSIGNED_LONG_LONG;
            }
        }
        public MPI_Datatype LONG_LONG {
            get {
                return LONG_LONG_INT;
            }
        }
    }

    class Platform_Native_MPI_Datatype : IMPI_DatatypeConstants {

        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_CHAR();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_UNSIGNED_CHAR();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_BYTE();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_SHORT();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_UNSIGNED_SHORT();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_INT();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_UNSIGNED();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_LONG();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_UNSIGNED_LONG();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_FLOAT();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_DOUBLE();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_LONG_DOUBLE();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_Datatype_LONG_LONG_INT();

        public MPI_Datatype CHAR {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_CHAR());
            }
        }

        public MPI_Datatype SIGNED_CHAR {
            get {
                throw new NotImplementedException();
                //return new MPI_Datatype(BoSSS_Get_MPI_Datatype_Sis());
            }
        }

        public MPI_Datatype UNSIGNED_CHAR {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_UNSIGNED_CHAR());
            }
        }

        public MPI_Datatype BYTE {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_BYTE());
            }
        }

        public MPI_Datatype WCHAR {
            get {
                throw new NotImplementedException();
                //return new MPI_Datatype(BoSSS_Get_MPI_Datatype_WC());
            }
        }

        public MPI_Datatype SHORT {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_SHORT());
            }
        }

        public MPI_Datatype UNSIGNED_SHORT {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_UNSIGNED_SHORT());
            }
        }

        public MPI_Datatype INT {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_INT());
            }
        }

        public MPI_Datatype UNSIGNED {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_UNSIGNED());
            }
        }

        public MPI_Datatype LONG {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_LONG());
            }
        }

        public MPI_Datatype UNSIGNED_LONG {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_UNSIGNED_LONG());
            }
        }

        public MPI_Datatype FLOAT {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_FLOAT());
            }
        }

        public MPI_Datatype DOUBLE {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_DOUBLE());
            }
        }

        public MPI_Datatype LONG_DOUBLE {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_LONG_DOUBLE());
            }
        }

        public MPI_Datatype LONG_LONG_INT {
            get {
                return new MPI_Datatype(BoSSS_Get_MPI_Datatype_LONG_LONG_INT());
            }
        }

        public MPI_Datatype UNSIGNED_LONG_LONG {
            get {
                throw new NotImplementedException();
                //return new MPI_Datatype(BoSSS_Get_MPI_Datatype_CHAR());
            }
        }

        public MPI_Datatype LONG_LONG {
            get {
                throw new NotImplementedException();
                //return new MPI_Datatype(BoSSS_Get_MPI_Datatype_LONG_LONG_INT());
            }
        }
    }
}
