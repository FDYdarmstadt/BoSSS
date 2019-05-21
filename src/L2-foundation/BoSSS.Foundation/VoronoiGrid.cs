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
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using ilPSP;


namespace BoSSS.Foundation.Grid.Voronoi {

    [Serializable]
    public class VoronoiGrid : Aggregation.AggregationGrid
    {
        //Generating Voronoi node of each cell jCell
        MultidimensionalArray voronoiNodes;

        public MultidimensionalArray VoronoiNodes {
            get { return voronoiNodes;}
        }
        
        //Velocity of each cell jCell
        public MultidimensionalArray NodeVelocity;

        public VoronoiGrid(IGrid pGrid,
            int[][] AggregationCells,
            MultidimensionalArray voronoiNodes)
            : base(pGrid, AggregationCells)
        {
            this.voronoiNodes = voronoiNodes;
            this.NodeVelocity = MultidimensionalArray.Create(NumberOfCells, 2);
        }

        VoronoiGrid() { }
    }
}
