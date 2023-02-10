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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Kraypis;

namespace BoSSS.Foundation.ConstrainedDGprojection {

    /// <summary>
    /// Abstract base class, template for different strategies 
    /// </summary>
    public abstract class BlockingStrategy {

        /// <summary>
        /// Returns lists which form the blocks of the additive-Schwarz domain decomposition.
        /// </summary>
        /// <returns>
        /// - outer enumeration: corresponds to domain-decomposition blocks
        /// - inner index: indices within the sub-blocks
        /// - content: local cell indices which form the respective additive-Schwarz block 
        /// </returns>
        abstract internal IEnumerable<List<int>> GetBlocking(GridData grdDat, CellMask mask);

        /// <summary>
        /// Number of blocs returned by <see cref="GetBlocking"/>
        /// </summary>
        internal abstract int GetNoOfBlocks();
    }


    /// <summary>
    /// creates a fixed number of blocks by using METIS
    /// </summary>
    public class METISBlockingStrategy : BlockingStrategy {

        /// <summary>
        /// Number of parts/additive Schwarz blocks on current MPI process (can be different on other processors)
        /// </summary>
        public int NoOfPartsOnCurrentProcess = 1;

        internal override IEnumerable<List<int>> GetBlocking(GridData grdDat, CellMask mask) {

            //if (cache != null) {
            //    return cache.Select(orgList => new List<int>(orgList)).ToArray();
            //}

            SubGrid sbgrd = new SubGrid(mask);
            int JComp = mask.NoOfItemsLocally; // number of local cells

            int[] xadj = new int[JComp + 1];
            List<int> adjncy = new List<int>();
            for (int j = 0; j < JComp; j++) {
                Debug.Assert(xadj[j] == adjncy.Count);

                int j_loc = sbgrd.SubgridIndex2LocalCellIndex[j];
                int[] neigh_j = grdDat.iLogicalCells.CellNeighbours[j_loc];
                int nCnt = 0;
                foreach (int jNeigh in neigh_j) {
                    //adjncy.AddRange(neigh_j);
                    int jNeigh_sbgrd = sbgrd.LocalCellIndex2SubgridIndex[jNeigh];
                    if (jNeigh_sbgrd >= 0 && jNeigh_sbgrd < JComp) {
                        adjncy.Add(jNeigh_sbgrd);
                        nCnt++;
                    } else {
                        //Console.WriteLine("Skipping external cell");
                    }
                }
                xadj[j + 1] = xadj[j] + nCnt;
            }
            Debug.Assert(xadj[JComp] == adjncy.Count);

            int MPIrank, MPIsize;
            MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out MPIrank);
            MPI.Wrappers.csMPI.Raw.Comm_Size(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out MPIsize);
            //if (MPIrank == 1)
            //    NoOfParts = 1;
            // dbg_launch();

            int[] part = new int[JComp];
            {
                if (NoOfPartsOnCurrentProcess > 1) {
                    int ncon = 1;
                    int edgecut = 0;
                    int[] options = new int[METIS.METIS_NOPTIONS];
                    METIS.SETDEFAULTOPTIONS(options);

                    options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = 1; // 
                    options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
                    options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)
                    options[(int)METIS.OptionCodes.METIS_OPTION_NUMBERING] = 0;

                    //Console.WriteLine("xadj.Count = {0}", xadj.Count());
                    //for (int i = 0; i < JComp + 1; i++) {
                    //    Console.WriteLine("xadj[{0}] = {1}", i, xadj[i]);
                    //}
                    Debug.Assert(xadj.Where(idx => idx > adjncy.Count).Count() == 0);
                    //Console.WriteLine("JComp = {0}", JComp);
                    //for (int j = 0; j < adjncy.Count; j++) {
                    //    if (adjncy.ElementAt(j) >= JComp)
                    //        Console.WriteLine("adjncy.ElemntAt({0}) = {1}", j, adjncy.ElementAt(j));
                    //}
                    Debug.Assert(adjncy.Where(j => j >= JComp).Count() == 0);

                    METIS.PARTGRAPHKWAY(
                            ref JComp, ref ncon,
                            xadj,
                            adjncy.ToArray(),
                            null,
                            null,
                            null,
                            ref NoOfPartsOnCurrentProcess,
                            null,
                            null,
                            options,
                            ref edgecut,
                            part);
                } else {
                    part.SetAll(0);
                }
            }

            {
                List<List<int>> _Blocks = NoOfPartsOnCurrentProcess.ForLoop(i => new List<int>((int)Math.Ceiling(1.1 * JComp / NoOfPartsOnCurrentProcess))).ToList();
                for (int j = 0; j < JComp; j++) { // loop over cells...
                    _Blocks[part[j]].Add(sbgrd.SubgridIndex2LocalCellIndex[j]); // cell `j` belongs to block `part[j]`
                }

                // remove empty blocks:
                for (int iB = 0; iB < _Blocks.Count; iB++) {
                    if (_Blocks[iB].Count <= 0) {
                        _Blocks.RemoveAt(iB);
                        iB--;
                    }
                }

                if (_Blocks.Count < NoOfPartsOnCurrentProcess)
                    Console.WriteLine("METIS WARNING: requested " + NoOfPartsOnCurrentProcess + " blocks, but got " + _Blocks.Count);

                return _Blocks.ToArray().Select(orgList => new List<int>(orgList)).ToArray();
            }
        }

        //List<int>[] cache;


        /// <summary>
        /// %
        /// </summary>
        internal override int GetNoOfBlocks() {
            return NoOfPartsOnCurrentProcess;
        }
    }


}
