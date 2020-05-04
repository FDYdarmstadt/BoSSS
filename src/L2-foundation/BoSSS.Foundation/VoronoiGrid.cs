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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using ilPSP;


namespace BoSSS.Foundation.Grid.Voronoi {

    [Serializable]
    public class VoronoiGrid : Aggregation.AggregationGrid
    {
        VoronoiBoundary boundary;

        VoronoiNodes nodes;

        public VoronoiNodes Nodes {
            get { return nodes; }
        }

        public VoronoiBoundary Boundary {
            get { return boundary; }
        }

        public int FirstCornerNodeIndice;

        VoronoiGrid() { }

        public VoronoiGrid(GridCommons pGrid,
            int[][] logialToGeometricalCellMap,
            VoronoiNodes nodes,
            VoronoiBoundary boundary)
            : base(pGrid, logialToGeometricalCellMap)
        {
            this.nodes = nodes;
            this.boundary = boundary;
        }

        public double NormalEdgeVelocity(int jEdge, double[] x, Vector normal)
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

            //transform if Edge is periodic
            if (this.iGridData.iGeomEdges.EdgeTags[jEdge] >= GridCommons.FIRST_PERIODIC_BC_TAG) 
            {
                int periodicEdgeTag = this.iGridData.iGeomEdges.EdgeTags[jEdge] - GridCommons.FIRST_PERIODIC_BC_TAG;
                AffineTrafo PerT = ((GridCommons)ParentGrid).PeriodicTrafo[periodicEdgeTag];
                posIn = PerT.Transform(posIn);
            };

            double result = VoronoiEdge.NormalVelocity(posOt, velOt, posIn, velIn, x, normal);
            return result;
        }

        public Vector CellVelocity(int jCell)
        {
            int jCellIn = this.iGridData.iGeomCells.GeomCell2LogicalCell[jCell];
            var velocity = new Vector(Nodes.Velocity[jCellIn, 0], Nodes.Velocity[jCellIn, 1]);
            return velocity;
        }
    }

    public class VoronoiBoundary
    {
        public Vector[] BoundingBox;

        public Vector[] Polygon;

        public byte[] EdgeTags;

        public IDictionary<byte, string> EdgeTagNames;

        public byte GetEdgeTagOfPolygonEdge(int index)
        {
            return EdgeTags[index];
        }
    }
}
