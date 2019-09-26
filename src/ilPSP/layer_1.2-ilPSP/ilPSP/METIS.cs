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
using MPI.Wrappers.Utils;

namespace ilPSP.Kraypis {

    /// <summary>
    /// wrappers to METIS functions (API version 5.1.0)
    /// </summary>
    /// <remarks>
    /// METIS can downloaded from http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
    /// </remarks>
    static public class METIS {

        static UnsafeMETIS m_METIS; 

        static METIS()
        {
            m_METIS = new UnsafeMETIS();
        }

        /// <summary>
        /// Length of the METIS options array
        /// </summary>
        public const int METIS_NOPTIONS = 40;

        /// <summary>
        /// Indices into METIS options array
        /// </summary>
        public enum OptionCodes {
            METIS_OPTION_PTYPE = 0,
            METIS_OPTION_OBJTYPE = 1,
            METIS_OPTION_CTYPE = 2,
            METIS_OPTION_IPTYPE = 3,
            METIS_OPTION_RTYPE = 4,
            METIS_OPTION_DBGLVL = 5,
            METIS_OPTION_NITER = 6,
            METIS_OPTION_NCUTS = 7,
            METIS_OPTION_SEED = 8,
            METIS_OPTION_NO2HOP = 9,
            METIS_OPTION_MINCONN = 10,
            METIS_OPTION_CONTIG = 11,
            METIS_OPTION_COMPRESS = 12,
            METIS_OPTION_CCORDER = 13,
            METIS_OPTION_PFACTOR = 14,
            METIS_OPTION_NSEPS = 15,
            METIS_OPTION_UFACTOR = 16,
            METIS_OPTION_NUMBERING = 17
        }

        /// <summary>
        /// Return codes returned by METIS functions
        /// </summary>
        public enum ReturnCodes {
            METIS_OK = 1,    /*!< Returned normally */
            METIS_ERROR_INPUT = -2,   /*!< Returned due to erroneous inputs and/or options */
            METIS_ERROR_MEMORY = -3,   /*!< Returned due to insufficient memory */
            METIS_ERROR = -4    /*!< Some other errors */
        }

        static public int PARTGRAPHKWAY(ref int nvtxs, ref int ncon, int[] xadj,
                                                int[] adjncy, int[] vwgt, int[] vsize,
                                                int[] adjwgt, ref int nparts, double[] tpwgts,
                                                double[] ubvec, int[] options, ref int objval, int[] part)
        {
            return m_METIS.PartGraphKway(ref nvtxs, ref ncon, xadj,
                                                adjncy, vwgt, vsize,
                                                adjwgt, ref nparts, tpwgts,
                                                ubvec, options, ref objval, part);
        }

        static public int PARTGRAPHRECURSIVE(ref int nvtxs, ref int ncon, int[] xadj,
                                                int[] adjncy, int[] vwgt, int[] vsize,
                                                int[] adjwgt, ref int nparts, double[] tpwgts,
                                                double[] ubvec, int[] options, ref int objval, int[] part)
        {
            return m_METIS.PartGraphRecursive(ref nvtxs, ref ncon, xadj,
                                                adjncy, vwgt, vsize,
                                                adjwgt, ref nparts, tpwgts,
                                                ubvec, options, ref objval, part);
        }


        static public int SETDEFAULTOPTIONS(int[] options) {
            return m_METIS.SetDefaultOptions(options);
        }


    }

    public sealed class UnsafeMETIS : DynLibLoader
    {
        // workaround for .NET bug:
        // https://connect.microsoft.com/VisualStudio/feedback/details/635365/runtimehelpers-initializearray-fails-on-64b-framework
        static PlatformID[] Helper()
        {
            PlatformID[] p = new PlatformID[2];
            p[0] = PlatformID.Win32NT;
            p[1] = PlatformID.Unix;
            return p;
        }

        /// <summary>
        /// ctor
        /// </summary>
        public UnsafeMETIS() :
            base(new string[] { "metis.dll", "libBoSSSnative_seq.so" },
                  new string[2][][],
                  new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.BoSSS_Prefix },
                  Helper(), //new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix, PlatformID.Unix, PlatformID.Unix },
                  new int[] { -1, -1 })
        { }

#pragma warning disable 649
        _PartGraphKway METIS_PartGraphKway;
        _PartGraphRecursive METIS_PartGraphRecursive;
        _SetDefaultOptions METIS_SetDefaultOptions;
#pragma warning restore 649

        /// <summary>
        /// see METIS manual;
        /// </summary>
        public unsafe delegate int _PartGraphKway(ref int nvtxs, ref int ncon, int[] xadj,
                                                int[] adjncy, int[] vwgt, int[] vsize,
                                                int[] adjwgt, ref int nparts, double[] tpwgts,
                                                double[] ubvec, int[] options, ref int objval, int[] part);

        public unsafe _PartGraphKway PartGraphKway
        {
            get { return METIS_PartGraphKway; }
        }

        

        /// <summary>
        /// see METIS manual;
        /// </summary>
        public unsafe delegate int _PartGraphRecursive(ref int nvtxs, ref int ncon, int[] xadj,
                                                int[] adjncy, int[] vwgt, int[] vsize,
                                                int[] adjwgt, ref int nparts, double[] tpwgts,
                                                double[] ubvec, int[] options, ref int objval, int[] part);

        public unsafe _PartGraphRecursive PartGraphRecursive
        {
            get { return METIS_PartGraphRecursive; }
        }


        /// <summary>
        /// see METIS manual;
        /// </summary>
        public unsafe delegate int _SetDefaultOptions(int[] options);

        public unsafe _SetDefaultOptions SetDefaultOptions
        {
            get { return METIS_SetDefaultOptions; }
        }
    }
    
}