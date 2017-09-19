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

using MPI.Wrappers;
using System.Runtime.InteropServices;

namespace ilPSP.Kraypis {

    /// <summary>
    /// Wrappers to ParMETIS functions
    /// </summary>
    /// <remarks>
    /// <b>IMPORTANT: Licensing issues:</b><br/>
    /// ParMETIS does not ship with a free license,
    /// neither their source nor binaries compiled from it can be shipped with
    /// this software;<br/>
    /// On some Linux distributions, like Debian (and derived distributions
    /// like Ubuntu) and Gentoo, this software can be installed from the local
    /// package repository (because they own a special permission from the
    /// METIS/ParMETIS developers).<br/>
    /// For Windows systems, the sources of the METIS software can be
    /// downloaded from
    /// http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview, they must be
    /// compiled individually. Instructions are given here:
    /// https://github.com/FDYdarmstadt/BoSSS/wiki/Compiling-ParMETIS-on-Windows
    /// </remarks>
    public class ParMETIS {

        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        [DllImport("parmetis", EntryPoint = "ParMETIS_V3_PartKway")]
        static extern int V3_PartKway(int[] vtxdist,
                                      int[] xadj,
                                      int[] adjncy,
                                      int[] vwgt,
                                      int[] adjwgt,
                                      ref int wgtflag,
                                      ref int numflag,
                                      ref int ncon,
                                      ref int nparts,
                                      float[] tpwgts,
                                      float[] ubvec,
                                      int[] options,
                                      ref int edgecut,
                                      int[] part,
                                      byte[] C_MPI_Comm);


        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        [DllImport("parmetis", EntryPoint = "ParMETIS_V3_PartGeomKway")]
        static extern int V3_PartGeomKway(int[] vtxdist,
                                          int[] xadj,
                                          int[] adjncy,
                                          int[] vwgt,
                                          int[] adjwgt,
                                          int[] wgtflag,
                                          ref int numflag,
                                          ref int ndims,
                                          float[] xyz,
                                          ref int ncon,
                                          ref int nparts,
                                          float[] tpwgts,
                                          float[] ubvec,
                                          int[] options,
                                          int[] edgecut,
                                          int[] part,
                                          byte[] C_MPI_Comm);


        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        [DllImport("parmetis", EntryPoint = "ParMETIS_V3_RefineKway")]
        static extern int ParMETIS_V3_RefineKway(int[] vtxdist,
                                                 int[] xadj,
                                                 int[] adjncy,
                                                 int[] vwgt,
                                                 int[] adjwgt,
                                                 ref int wgtflag,
                                                 ref int numflag,
                                                 ref int ncon,
                                                 ref int nparts,
                                                 float[] tpwgts,
                                                 float[] ubvec,
                                                 int[] options,
                                                 ref int edgecut,
                                                 int[] part,
                                                 byte[] MPI_Comm);

        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        public static METIS.ReturnCodes V3_PartKway(
            int[] vtxdist,
            int[] xadj,
            int[] adjncy,
            int[] vwgt,
            int[] adjwgt,
            ref int wgtflag,
            ref int numflag,
            ref int ncon,
            ref int nparts,
            float[] tpwgts,
            float[] ubvec,
            int[] options,
            ref int edgecut,
            int[] part,
            MPI_Comm C_MPI_Comm) {

            return (METIS.ReturnCodes)V3_PartKway(
                vtxdist, xadj, adjncy, vwgt, adjwgt, ref wgtflag, ref numflag, ref ncon, ref nparts, tpwgts, ubvec, options, ref edgecut, part, csMPI.Raw._COMM.Comm_f2c(C_MPI_Comm));
        }

        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        public static METIS.ReturnCodes V3_PartGeomKway(
            int[] vtxdist,
            int[] xadj,
            int[] adjncy,
            int[] vwgt,
            int[] adjwgt,
            int[] wgtflag,
            ref int numflag,
            ref int ndims,
            float[] xyz,
            ref int ncon,
            ref int nparts,
            float[] tpwgts,
            float[] ubvec,
            int[] options,
            int[] edgecut,
            int[] part,
            ref MPI_Comm MPI_Comm) {

            return (METIS.ReturnCodes)V3_PartGeomKway(
                vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, ref numflag, ref ndims, xyz, ref ncon, ref nparts, tpwgts, ubvec, options, edgecut, part, csMPI.Raw._COMM.Comm_f2c(MPI_Comm));
        }

        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        public static METIS.ReturnCodes V3_RefineKway(
            int[] vtxdist,
            int[] xadj,
            int[] adjncy,
            int[] vwgt,
            int[] adjwgt,
            ref int wgtflag,
            ref int numflag,
            ref int ncon,
            ref int nparts,
            float[] tpwgts,
            float[] ubvec,
            int[] options,
            ref int edgecut,
            int[] part,
            MPI_Comm MPI_Comm) {

            return (METIS.ReturnCodes)ParMETIS_V3_RefineKway(
                vtxdist, xadj, adjncy, vwgt, adjwgt, ref wgtflag, ref numflag, ref ncon, ref nparts, tpwgts, ubvec, options, ref edgecut, part, csMPI.Raw._COMM.Comm_f2c(MPI_Comm));
        }
    }
}
