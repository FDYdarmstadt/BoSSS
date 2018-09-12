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

using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Extensions for <see cref="GridCommons"/> and <see cref="GridData"/>.
    /// </summary>
    static public class GridCommons_extensions {

        /// <summary>
        /// Creates fields that mark the boundary conditions.
        /// </summary>
        public static DGField[] BoundaryMark(this IGridData m_gridData) {
            List<DGField> demo_fields = new List<DGField>();

            // mark boundary cells
            {
                SinglePhaseField boundary = new SinglePhaseField(new Basis(m_gridData, 0), "AllBndCells");
                int[,] Edges = m_gridData.iLogicalEdges.CellIndices;
                int E = Edges.GetLength(0);
                for (int ee = 0; ee < E; ee++) {
                    int j = Edges[ee, 0];
                    if (Edges[ee, 1] < 0)
                        boundary.SetMeanValue(j, boundary.GetMeanValue(j) + 1);
                }
                demo_fields.Add(boundary);
            }

            {
                DGField[] boundaries = new DGField[m_gridData.EdgeTagNames.Keys.Max() + 1];
                foreach (byte edgeTag in m_gridData.EdgeTagNames.Keys) {
                    if (edgeTag != 0)
                        boundaries[edgeTag] = new SinglePhaseField(new Basis(m_gridData, 0), "Boundary_" + m_gridData.EdgeTagNames[edgeTag]);
                }
                boundaries[0] = new SinglePhaseField(new Basis(m_gridData, 0), "UnspecifiedBoundary");

                for (int ii = 0; ii < boundaries.Length; ii++)
                    if (boundaries[ii] == null) {
                        boundaries[ii] = new SinglePhaseField(new Basis(m_gridData, 0), "notused-" + ii);
                    }


                int[,] Edges = m_gridData.iGeomEdges.CellIndices;
                byte[] EdgeTags = m_gridData.iGeomEdges.EdgeTags;

                int E = Edges.GetLength(0);
                for (int ee = 0; ee < E; ee++) {
                    int j = Edges[ee, 0];
                    if (Edges[ee, 1] < 0) {
                        int tag = EdgeTags[ee];

                        boundaries[tag].SetMeanValue(j, boundaries[tag].GetMeanValue(j) + 1);
                    }
                }
                demo_fields.AddRange(boundaries);
            }


            return demo_fields.ToArray();

        }
    }
}
