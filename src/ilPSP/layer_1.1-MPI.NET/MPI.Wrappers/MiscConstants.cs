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
    /// misc MPI constants
    /// </summary>
    public interface IMiscConstants {
        /// <summary> match any source rank </summary>
        int ANY_SOURCE {
            get;
        }      /*  */
        /// <summary> rank of null process </summary>
        int PROC_NULL {
            get;
        }      /*  */
        /// <summary>  </summary>
        int ROOT {
            get;
        }
        /// <summary>  match any message tag </summary>
        int ANY_TAG {
            get;
        }      /* */
        /// <summary>  max proc. name length </summary>
        int MAX_PROCESSOR_NAME {
            get;
        }     /* */
        /// <summary>  max error message length </summary>
        int MAX_ERROR_STRING {
            get;
        }     /* */
        /// <summary> max object name length </summary>
        int MAX_OBJECT_NAME {
            get;
        }      /*  */
        /// <summary> undefined stuff </summary>
        int UNDEFINED {
            get;
        }  /*  */
        /// <summary> cartesian topology </summary>
        int CART {
            get;
        }       /*  */
        /// <summary> graph topology </summary>
        int GRAPH {
            get;
        }       /*  */
        /// <summary> invalid key value </summary>
        int KEYVAL_INVALID {
            get;
        }      /*  */

        /// <summary> NULL Request</summary>
        MPI_Request MPI_REQUEST_NULL {
            get;
        }
    }

    class OPENMPI_MiscConstants : IMiscConstants {

        OPENMPI_Converter m_conv;

        public OPENMPI_MiscConstants(OPENMPI_Converter conv) {
            m_conv = conv;
        }

        bool deleyedInitdone = false;

        private void DeleaedInit() {
            if (deleyedInitdone)
                return;
            deleyedInitdone = true;
            LoadSymb("ompi_request_null", ref m_MPI_REQUEST_NULL);
        }

        private void LoadSymb(string name, ref MPI_Request sym) {
            string errstr;
            IntPtr addr = Utils.DynamicLibraries.LoadSymbol(m_conv.LibHandle, name, out errstr);
            if (addr == IntPtr.Zero)
                throw new ApplicationException("OpenMPI error: unable to load symbol '" + name + "' from library '" + m_conv.CurrentLibraryName + "', Error string: >" + errstr + "<;");
            //Console.WriteLine("val of '" + name + "' is: " + addr);
            sym.m1 = m_conv.MPI_Request_c2f(addr);
        }
        /// <summary> match any source rank </summary>
        public int ANY_SOURCE {
            get {
                return -1;
            }
        }      /*  */
        /// <summary> rank of null process </summary>
        public int PROC_NULL {
            get {
                return -2;
            }
        }      /*  */
        /// <summary>  </summary>
        public int ROOT {
            get {
                return -4;
            }
        }
        /// <summary>  match any message tag </summary>
        public int ANY_TAG {
            get {
                return -1;
            }
        }      /* */
        /// <summary>  max proc. name length </summary>
        public int MAX_PROCESSOR_NAME {
            get {
                return 255;
            }
        }     /* */
        /// <summary>  max error message length </summary>
        public int MAX_ERROR_STRING {
            get {
                return 255;
            }
        }     /* */
        /// <summary> max object name length </summary>
        public int MAX_OBJECT_NAME {
            get {
                return 63;
            }
        }      /*  */
        /// <summary> undefined stuff </summary>
        public int UNDEFINED {
            get {
                return -32766;
            }
        }  /*  */
        /// <summary> cartesian topology </summary>
        public int CART {
            get {
                return 1;
            }
        }       /*  */
        /// <summary> graph topology </summary>
        public int GRAPH {
            get {
                return 2;
            }
        }       /*  */
        /// <summary> invalid key value </summary>
        public int KEYVAL_INVALID {
            get {
                return -1;
            }
        }      /*  */

        MPI_Request m_MPI_REQUEST_NULL;
        public MPI_Request MPI_REQUEST_NULL {
            get {
                return m_MPI_REQUEST_NULL;
            }
        }
    }

    class MPICH_MiscConstants : IMiscConstants {

        /// <summary> match any source rank </summary>
        public int ANY_SOURCE {
            get {
                return -2;
            }
        }      /*  */
        /// <summary> rank of null process </summary>
        public int PROC_NULL {
            get {
                return -1;
            }
        }      /*  */
        /// <summary>  </summary>
        public int ROOT {
            get {
                return -3;
            }
        }
        /// <summary>  match any message tag </summary>
        public int ANY_TAG {
            get {
                return -1;
            }
        }      /* */
        /// <summary>  max proc. name length </summary>
        public int MAX_PROCESSOR_NAME {
            get {
                return 128 - 1;
            }
        }     /* */
        /// <summary>  max error message length </summary>
        public int MAX_ERROR_STRING {
            get {
                return 511;
            }
        }     /* */
        /// <summary> max object name length </summary>
        public int MAX_OBJECT_NAME {
            get {
                return 127;
            }
        }      /*  */
        /// <summary> undefined stuff </summary>
        public int UNDEFINED {
            get {
                return -32766;
            }
        }  /*  */
        /// <summary> cartesian topology </summary>
        public int CART {
            get {
                return 1;
            }
        }       /*  */
        /// <summary> graph topology </summary>
        public int GRAPH {
            get {
                return 2;
            }
        }       /*  */
        /// <summary> invalid key value </summary>
        public int KEYVAL_INVALID {
            get {
                return 603979776;
            }
        }      /*  */

        public MPI_Request MPI_REQUEST_NULL {
            get {
                MPI_Request r;
                r.m1 = 0x2c000000;
                return r;
            }
        }
    }


    class Platform_Native_MiscConstants : IMiscConstants {

        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_ANY_SOURCE();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_PROC_NULL();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_ROOT();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_ANY_TAG();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_MAX_PROCESSOR_NAME();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_MAX_ERROR_STRING();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_MAX_OBJECT_NAME();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_UNDEFINED();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_CART();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_GRAPH();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_KEYVAL_INVALID();
        [DllImport("Platform_Native")]
        extern static int BoSSS_Get_MPI_REQUEST_NULL();


        /// <summary> match any source rank </summary>
        public int ANY_SOURCE {
            get {
                return BoSSS_Get_MPI_ANY_SOURCE();
            }
        }
        /// <summary> rank of null process </summary>
        public int PROC_NULL {
            get {
                return BoSSS_Get_MPI_PROC_NULL();
            }
        }
        /// <summary>  </summary>
        public int ROOT {
            get {
                return BoSSS_Get_MPI_ROOT();
            }
        }
        /// <summary>  match any message tag </summary>
        public int ANY_TAG {
            get {
                return BoSSS_Get_MPI_ANY_TAG();
            }
        }
        /// <summary>  max proc. name length </summary>
        public int MAX_PROCESSOR_NAME {
            get {
                return BoSSS_Get_MPI_MAX_PROCESSOR_NAME();
            }
        }
        /// <summary>  max error message length </summary>
        public int MAX_ERROR_STRING {
            get {
                return BoSSS_Get_MPI_MAX_ERROR_STRING();
            }
        }
        /// <summary> max object name length </summary>
        public int MAX_OBJECT_NAME {
            get {
                return BoSSS_Get_MPI_MAX_OBJECT_NAME();
            }
        }
        /// <summary> undefined stuff </summary>
        public int UNDEFINED {
            get {
                return BoSSS_Get_MPI_UNDEFINED();
            }
        }
        /// <summary> cartesian topology </summary>
        public int CART {
            get {
                return BoSSS_Get_MPI_CART();
            }
        }
        /// <summary> graph topology </summary>
        public int GRAPH {
            get {
                return BoSSS_Get_MPI_GRAPH();
            }
        }
        /// <summary> invalid key value </summary>
        public int KEYVAL_INVALID {
            get {
                return BoSSS_Get_MPI_KEYVAL_INVALID();
            }
        }

        public MPI_Request MPI_REQUEST_NULL {
            get {
                MPI_Request ret = default(MPI_Request);
                ret.m1 = BoSSS_Get_MPI_REQUEST_NULL();
                return ret;
            }
        }
    }
}
