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

        public VoronoiInfo Info {
            get { return info; }
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

        public double EdgeVelocity(int jEdge, double[] x, Vector normal)
        {
            int jCellIn = this.iGridData.iGeomEdges.CellIndices[jEdge, 1];
            int jCellOut = this.iGridData.iGeomEdges.CellIndices[jEdge, 0];
            int jCell_in = this.iGridData.iGeomCells.GeomCell2LogicalCell[jCellIn];
            int jCell_ot = this.iGridData.iGeomCells.GeomCell2LogicalCell[jCellOut];
            MultidimensionalArray positions = Nodes.Positions;
            MultidimensionalArray velocities = Nodes.Velocity;

            double[] posOt = positions.GetRow(jCell_ot);
            double[] posIn = positions.GetRow(jCell_in);
            double[] velOt = velocities.GetRow(jCell_ot);
            double[] velIn = velocities.GetRow(jCell_in);
            double result = VoronoiEdge.NormalVelocity(posOt, velOt, posIn, velIn, x, normal);
            return result;
        }
    }

    public class VoronoiInfo
    {
        public Vector[] BoundingBox;
        public Vector[] Boundary;
    }
}
