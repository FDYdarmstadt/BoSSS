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
using MPI.Wrappers;

namespace ilPSP.Kraypis {

    /// <summary>
    /// wrappers to METIS functions
    /// </summary>
    /// <remarks>
    /// <b>IMPORTANT: Licencing issues:</b><br/>
    /// METIS and ParMETIS <see cref="ParMETIS_V3"/> do not ship with free
    /// license, neither their source nor  binaries compiled from it can be
    /// shipped with this software;<br/>
    /// On some Linux distributions, like Debian (and derived distributions
    /// like Ubuntu) and Gentoo, this software can be installed from the local
    /// package repository (because they own a special permission from the
    /// METIS/ParMETIS developers).<br/>
    /// For Windows systems, the sources of the METIS software can be
    /// downloaded from http://glaros.dtc.umn.edu/gkhome/metis/metis/overview,
    /// they must be compiled individually. Visual Studio Project files are provided.
    /// </remarks>
    public class METIS {

        /// <summary>
        /// see METIS manual;
        /// </summary>
        [DllImport("metis", EntryPoint = "METIS_PartGraphKway")]
        public static extern void PartGraphKway(ref int nvtxs,
                                                int[] xadj, int[] adjncy,
                                                int[] vwgt, int[] adjwgt,
                                                ref int wgtflag, ref int numflag, ref int nparts,
                                                int[] options, ref int edgecut, int[] part);
        /// <summary>
        /// see METIS manual;
        /// </summary>
        [DllImport("metis", EntryPoint = "METIS_PartGraphRecursive")]
        public static extern void PartGraphRecursive(ref int nvtxs,
                                                     int[] xadj, int[] adjncy,
                                                     int[] vwgt, int[] adjwgt,
                                                     ref int wgtflag, ref int numflag, ref int nparts,
                                                     int[] options, ref int edgecut, int[] part);
    }
}
