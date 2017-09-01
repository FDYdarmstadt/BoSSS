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
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using ilPSP;

namespace BoSSS.Foundation {
    
    /// <summary>
    /// Extension methods 
    /// </summary>
    public static class DGFieldExtensions {


        /// <summary>
        /// Utility function to evaluate DG fields together with global nodes.
        /// </summary>
        /// <param name="f">DG field to evaluate.</param>
        /// <param name="NS">node set.</param>
        /// <param name="cm">optional cell mask.</param>
        /// <returns>
        /// Global nodes and function values at these nodes;
        /// 1st index: cell index;
        /// 2nd index: node index;
        /// 3rd index: spatial direction
        /// </returns>
        public static MultidimensionalArray Evaluate(this DGField f, NodeSet NS, CellMask cm = null) {
            IGridData g = f.GridDat;
            if(cm == null)
                cm = CellMask.GetFullMask(g);
            int D = g.SpatialDimension;

            MultidimensionalArray ret = MultidimensionalArray.Create(new int[] { cm.NoOfItemsLocally, NS.NoOfNodes, D + 1 });

            int jsub = 0;
            foreach(Chunk ck in cm) {
                var globNodes = ret.ExtractSubArrayShallow(new int[] { jsub, 0, 0 }, new int[] { jsub + ck.Len - 1, NS.NoOfNodes - 1, D - 1 });
                var fldValues = ret.ExtractSubArrayShallow(new int[] { jsub, 0, D }, new int[] { jsub + ck.Len - 1, NS.NoOfNodes - 1, D - 1 });

                g.TransformLocal2Global(NS, ck.i0, ck.Len, globNodes);
                f.Evaluate(ck.i0, ck.Len, NS, fldValues);
            }
            
            return ret;
        }


    }
}
