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

using System.Runtime.InteropServices;

namespace ilPSP.Kraypis {

    /// <summary>
    /// wrappers to METIS functions (API version 5.1.0)
    /// </summary>
    /// <remarks>
    /// METIS can downloaded from http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
    /// </remarks>
    public class METIS {

        /// <summary>
        /// Length of the METIS options array
        /// </summary>
        public const int METIS_NOPTIONS = 40;

        /// <summary>
        /// Indices into METIS options array
        /// </summary>
        public enum OPTIONS_CODES {
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
        /// see METIS manual;
        /// </summary>
        [DllImport("metis", EntryPoint = "METIS_PartGraphKway")]
        public static extern void PartGraphKway(ref int nvtxs, ref int ncon,
                                                int[] xadj, int[] adjncy,
                                                int[] vwgt, int[] vsize, int[] adjwgt,
                                                ref int nparts, double[] tpwgts, double[] ubvec,
                                                int[] options, ref int objval, int[] part);

        /// <summary>
        /// see METIS manual;
        /// </summary>
        [DllImport("metis", EntryPoint = "METIS_PartGraphRecursive")]
        public static extern void PartGraphRecursive(ref int nvtxs, ref int ncon,
                                                     int[] xadj, int[] adjncy,
                                                     int[] vwgt, int[] vsize, int[] adjwgt,
                                                     ref int nparts, double[] tpwgts, double[] ubvec,
                                                     int[] options, ref int objval, int[] part);
    }
}
