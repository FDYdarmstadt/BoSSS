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

        /// <summary>
        /// see METIS manual;
        /// </summary>
        [DllImport("metis", EntryPoint = "METIS_WPartGraphRecursive")]
        public static extern void WPartGraphRecursive(ref int n,
                                                      int[] xadj, int[] adjncy,
                                                      int[] vwgt, int[] adjwgt,
                                                      ref int wgtflag, ref int numflag, ref int nparts,
                                                      float[] tpwgts,
                                                      int[] options, ref int edgecut, int[] part);
        /// <summary>
        /// see METIS manual;
        /// </summary>
        [DllImport("metis", EntryPoint = "METIS_WPartGraphKway")]
        public static extern void WPartGraphKway(ref int n,
                                                 int[] xadj, int[] adjncy,
                                                 int[] vwgt, int[] adjwgt,
                                                 ref int wgtflag, ref int numflag, ref int nparts,
                                                 float[] tpwgts,
                                                 int[] options, ref int edgecut, int[] part);
        /// <summary>
        /// see METIS manual;
        /// </summary>
        [DllImport("metis", EntryPoint = "METIS_WPartGraphVKway")]
        public static extern void WPartGraphVKway(ref int n,
                                                  int[] xadj, int[] adjncy,
                                                  int[] vwgt, int[] adjwgt,
                                                  ref int wgtflag, ref int numflag, ref int nparts,
                                                  float[] tpwgts,
                                                  int[] options, ref int edgecut, int[] part);
    }

    /// <summary>
    /// Wrappers to ParMETIS functions
    /// </summary>
    /// <remarks>
    /// <b>IMPORTANT: Licensing issues:</b><br/>
    /// METIS <see cref="METIS"/> and ParMETIS do not ship with free license,
    /// neither their source nor  binaries compiled from it can be shipped with
    /// this software;<br/>
    /// On some Linux distributions, like Debian (and derived distributions
    /// like Ubuntu) and Gentoo, this software can be installed from the local
    /// package repository (because they own a special permission from the
    /// METIS/ParMETIS developers).<br/>
    /// For Windows systems, the sources of the METIS software can be
    /// downloaded from
    /// http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview, they must be
    /// compiled individually. Visual Studio Project files are provided.
    /// </remarks>
    public class ParMETIS_V3 {

        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        [DllImport("parmetis", EntryPoint = "ParMETIS_V3_PartKway")]
        static extern void V3_PartKway(int[] vtxdist,
                                       int[] xadj,
                                       int[] adjncy,
                                       int[] vwgt,
                                       int[] adjwgt,
                                       int[] wgtflag,
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
        /// This routine is used to compute a k-way partitioning of a graph on
        /// p processors using the multilevel k-way multi-constraint
        /// partitioning algorithm.
        /// </summary>
        /// <param name="vtxdist">
        /// This array describes how the vertices of the graph are distributed
        /// among the processors. Its contents are identical for every
        /// processor.
        /// </param>
        /// <param name="xadj">
        /// Stores indices into <paramref name="adjncy"/>
        /// </param>
        /// <param name="adjncy">
        /// Store the (local) adjacency structure of the graph at each
        /// processor
        /// </param>
        /// <param name="vwgt">
        /// Stores the weights of the vertices
        /// </param>
        /// <param name="adjwgt">
        /// Stores the weights of the edges
        /// </param>
        /// <param name="wgtflag">
        /// This is used to indicate if the graph is weighted. Can take one of
        /// four values:
        /// <list type="bullet">
        ///     <item>
        ///     0: No weights (<paramref name="vwgt"/> and
        ///     <paramref name="adjwgt"/> are both null
        ///     </item>
        ///     <item>
        ///     1: Weights on the edges only (<paramref name="vwgt"/> is null)
        ///     </item>
        ///     <item>
        ///     2: Weights on the vertices only (<paramref name="adjwgt"/> is
        ///     null)
        ///     </item>
        ///     <item>
        ///     3: Weights on both the vertices and edges
        ///     </item>
        /// </list>
        /// </param>
        /// <param name="numflag">
        /// This is used to indicate the numbering scheme that is used for the
        /// <paramref name="vtxdist"/>, <paramref name="xadj"/>,
        /// <paramref name="adjncy"/>, and <paramref name="part"/> arrays. Can
        /// take one of two values:
        /// <list type="bullet">
        ///     <item>0: C-style numbering that starts from 0</item>
        ///     <item>1: Fortran-style numbering that starts from 1.</item>
        /// </list>
        /// </param>
        /// <param name="ncon">
        /// This is used to specify the number of weights that each vertex has.
        /// It is also the number of balance constraints that must be satisfied.
        /// </param>
        /// <param name="nparts">
        /// This is used to specify the number of sub-domains that are desired.
        /// Note that the number of sub-domains is independent of the number of
        /// processors that call this routine.
        /// </param>
        /// <param name="tpwgts">
        /// An array of size <paramref name="ncon"/>x<paramref name="nparts"/>
        /// that is used to specify the fraction of vertex weight that should
        /// be distributed to each sub-domain for each balance constraint. If
        /// all of the sub-domains are to be of the same size for every vertex
        /// weight, then each of the
        /// <paramref name="nparts"/>*<paramref name="nparts"/> elements should
        /// be set to a value of 1 over <paramref name="nparts"/>. If
        /// <paramref name="ncon"/> is greater than 1, the target sub-domain
        /// weights for each sub-domain are stored contiguously (similar to the
        /// <paramref name="vwgt"/> array). Note that the sum of all of the
        /// <paramref name="tpwgts"/> for a given vertex weight should be one.
        /// </param>
        /// <param name="ubvec">
        /// An array of size <paramref name="ncon"/> that is used to specify
        /// the imbalance tolerance for each vertex weight, with 1 being
        /// perfect balance and <paramref name="nparts"/> being perfect
        /// imbalance. A value of 1.05 for each of the <paramref name="ncon"/>
        /// weights is recommended
        /// </param>
        /// <param name="options">
        /// This is an array of integers that is used to pass additional
        /// parameters for the routine. The first element (i.e., options[0])
        /// can take either the value of 0 or 1. If it is 0, then the default
        /// values are used, otherwise the remaining two elements of options
        /// are interpreted as follows:
        /// <list type="bullet">
        ///     <item>
        ///     options[1]: This specifies the level of information to be
        ///     returned during the execution of the algorithm. Timing
        ///     information can be obtained by setting this to 1. Additional
        ///     options for this parameter can be obtained by looking at
        ///     parmetis.h. The numerical values there should be added to
        ///     obtain the correct value. The default value is 0
        ///     </item>
        ///     <item>
        ///     options[2]: This is the random number seed for the routine
        ///     </item>
        /// </list>
        /// </param>
        /// <param name="edgecut">
        /// Upon successful completion, the number of edges that are cut by the
        /// partitioning is written to this parameter.
        /// </param>
        /// <param name="part">
        /// This is an array of size equal to the number of locally-stored
        /// vertices. Upon successful completion the partition vector of the
        /// locally-stored vertices is written to this array
        /// </param>
        /// <param name="C_MPI_Comm">
        /// This is a pointer to the MPI communicator of the processes that
        /// call PARMETIS. For most programs this will point to MPI_COMM_WORLD
        /// </param>
        public static void V3_PartKway(int[] vtxdist,
                                       int[] xadj,
                                       int[] adjncy,
                                       int[] vwgt,
                                       int[] adjwgt,
                                       int[] wgtflag,
                                       ref int numflag,
                                       ref int ncon,
                                       ref int nparts,
                                       float[] tpwgts,
                                       float[] ubvec,
                                       int[] options,
                                       ref int edgecut,
                                       int[] part,
                                       MPI_Comm C_MPI_Comm) {
            V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, ref numflag, ref ncon, ref nparts, tpwgts, ubvec, options, ref edgecut, part, csMPI.Raw._COMM.Comm_f2c(C_MPI_Comm));
        }

        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        [DllImport("parmetis", EntryPoint = "ParMETIS_V3_PartGeomKway")]
        static extern void V3_PartGeomKway(int[] vtxdist,
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
        public static void V3_PartGeomKway(
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

            V3_PartGeomKway(
                vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, ref numflag, ref ndims, xyz, ref ncon, ref nparts, tpwgts, ubvec, options, edgecut, part, csMPI.Raw._COMM.Comm_f2c(MPI_Comm));
        }

        /// <summary>
        /// see ParMETIS manual;
        /// </summary>
        [DllImport("parmetis", EntryPoint = "ParMETIS_V3_RefineKway")]
        static extern void ParMETIS_V3_RefineKway(
            int[] vtxdist,
            int[] xadj,
            int[] adjncy,
            int[] vwgt,
            int[] adjwgt,
            int[] wgtflag,
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
        public static void V3_RefineKway(
            int[] vtxdist,
            int[] xadj,
            int[] adjncy,
            int[] vwgt,
            int[] adjwgt,
            int[] wgtflag,
            ref int numflag,
            ref int ncon,
            ref int nparts,
            float[] tpwgts,
            float[] ubvec,
            int[] options,
            ref int edgecut,
            int[] part,
            MPI_Comm MPI_Comm) {

            ParMETIS_V3_RefineKway(
                vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, ref numflag, ref ncon, ref nparts, tpwgts, ubvec, options, ref edgecut, part, csMPI.Raw._COMM.Comm_f2c(MPI_Comm));
        }
    }
}
