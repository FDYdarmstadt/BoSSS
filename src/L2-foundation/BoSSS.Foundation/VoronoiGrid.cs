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
using System.Runtime.Serialization;
using BoSSS.Platform.LinAlg;
using ilPSP;


namespace BoSSS.Foundation.Grid.Voronoi {

    [Serializable]
    public class VoronoiGrid : Aggregation.AggregationGrid
    {
        VoronoiInfo info;

        VoronoiNodes nodes;

        public VoronoiNodes Nodes {
            get { return nodes; }
        }

        public VoronoiBoundary Boundary {
            get { return info.Boundary; }
        }

        public VoronoiInfo Info {
            get { return info; }
        }

        public VoronoiGrid(IGrid pGrid,
            int[][] logialToGeometricalCellMap,
            MultidimensionalArray nodePositions,
            VoronoiInfo info
            )
            : base(pGrid, logialToGeometricalCellMap)
        {
            SetVoronoiNodes(nodePositions);
            this.info = info;
        }

        void SetVoronoiNodes(MultidimensionalArray nodePositions)
        {
            long[] voronoiNodeIds = GetGlobalCellIds(this.AggregationCells);
            nodes = new VoronoiNodes(nodePositions, voronoiNodeIds);
        }

        long[] GetGlobalCellIds(AggCell[] cells)
        {
            long[] ids = new long[cells.Length];
            for (int i = 0; i < cells.Length; ++i)
            {
                ids[i] = cells[i].GlobalID;
            }
            return ids;
        }

        public VoronoiGrid(IGrid pGrid,
            int[][] logialToGeometricalCellMap,
            VoronoiNodes nodes,
            VoronoiInfo voronoiInfo)
            : base(pGrid, logialToGeometricalCellMap)
        {
            this.nodes = nodes;
            info = voronoiInfo;
        }

        VoronoiGrid() { }
    }

    public class VoronoiInfo
    {
        public VoronoiBoundary Boundary;
        public int NumberOfLloydIterations = 10;
        public int FirstCellNode_indice = 0;
    }

    public class VoronoiBoundary
    {
        public Vector[] Polygon;
    }
}
