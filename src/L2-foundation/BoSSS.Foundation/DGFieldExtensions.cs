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
using System.Runtime.CompilerServices;
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
            if (cm == null)
                cm = CellMask.GetFullMask(g);
            int D = g.SpatialDimension;

            MultidimensionalArray ret = MultidimensionalArray.Create(new int[] { cm.NoOfItemsLocally, NS.NoOfNodes, D + 1 });

            int jsub = 0;
            foreach (Chunk ck in cm) {
                var globNodes = ret.ExtractSubArrayShallow(new int[] { jsub, 0, 0 }, new int[] { jsub + ck.Len - 1, NS.NoOfNodes - 1, D - 1 });
                var fldValues = ret.ExtractSubArrayShallow(new int[] { jsub, 0, D }, new int[] { jsub + ck.Len - 1, NS.NoOfNodes - 1, D - 1 });

                g.TransformLocal2Global(NS, ck.i0, ck.Len, globNodes);
                f.Evaluate(ck.i0, ck.Len, NS, fldValues);
            }

            return ret;
        }

        /// <summary>
        /// Central Difference Divergence, <see cref="DGField.DivergenceByFlux{T}"/>
        /// </summary>
        public static T DivergenceByFlux<T>(this VectorField<T> Vec, SubGrid optionalSubGrid = null, SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) where T : DGField {
            T myDiv = Vec[0].CloneAs() as T;
            myDiv.Identification = "div";
            myDiv.Clear();
            myDiv.DivergenceByFlux(1.0, Vec, optionalSubGrid, bndMode);

            return myDiv;
        }

        /// <summary>
        /// Broken Divergence, <see cref="DGField.Divergence{T}"/>
        /// </summary>
        public static T Divergence<T>(this VectorField<T> Vec, CellMask em = null) where T : DGField {
            T myDiv = Vec[0].CloneAs() as T;
            myDiv.Identification = "div";
            myDiv.Clear();
            myDiv.Divergence(1.0, Vec, em);

            return myDiv;
        }

        /// <summary>
        /// Central Difference Divergence, <see cref="DGField.DivergenceByFlux{T}"/>
        /// </summary>
        public static T DivergenceByFlux<T>(this IEnumerable<T> Vec, SubGrid optionalSubGrid = null, SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) where T : DGField {
            return (new VectorField<T>(Vec.ToArray())).DivergenceByFlux();

        }

        /// <summary>
        /// Broken Divergence, <see cref="DGField.Divergence{T}"/>
        /// </summary>
        public static T Divergence<T>(this IEnumerable<T> Vec, CellMask em = null) where T : DGField {
            return (new VectorField<T>(Vec.ToArray())).Divergence();
        }

        /// <summary>
        /// Central Difference Curl, <see cref="VectorField{T}.Curl3DByFlux"/>
        /// </summary>
        public static VectorField<T> Curl3DByFlux<T>(this VectorField<T> Vec, SubGrid optionalSubGrid = null, SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) where T : DGField {
            var Curl = Vec.CloneAs();
            Curl.Clear();
            Curl.Curl3DByFlux(1.0, Vec, optionalSubGrid, bndMode);

            return Curl;
        }

        /// <summary>
        /// Broken Curl, <see cref="VectorField{T}.Curl3D"/>
        /// </summary>
        public static VectorField<T> Curl3D<T>(this VectorField<T> Vec, CellMask em = null) where T : DGField {
            var Curl = Vec.CloneAs();
            Curl.Clear();
            Curl.Curl3D(1.0, Vec, em);

            return Curl;
        }

        /// <summary>
        /// Central Difference Curl, <see cref="VectorField{T}.Curl3DByFlux"/>
        /// </summary>
        public static VectorField<T> Curl3DByFlux<T>(this IEnumerable<T> Vec, SubGrid optionalSubGrid = null, SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) where T : DGField {
            return (new VectorField<T>(Vec.ToArray())).Curl3DByFlux(optionalSubGrid, bndMode);
        }

        /// <summary>
        /// Broken Curl, <see cref="VectorField{T}.Curl3D"/>
        /// </summary>
        public static VectorField<T> Curl3D<T>(this IEnumerable<T> Vec, CellMask em = null) where T : DGField {
            return (new VectorField<T>(Vec.ToArray())).Curl3D(em);
        }


        /// <summary>
        /// Central Difference 2D Curl, <see cref="DGField.Curl2DByFlux{T}"/>
        /// </summary>
        public static T Curl2DByFlux<T>(this VectorField<T> Vec, SubGrid optionalSubGrid = null, SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) where T : DGField {
            T Curl = Vec[0].CloneAs() as T;
            Curl.Clear();
            Curl.Clear();
            Curl.Curl2DByFlux(1.0, Vec, optionalSubGrid, bndMode);

            return Curl;
        }

        /// <summary>
        /// Broken 2D Curl, <see cref="DGField.Curl2D{T}"/>
        /// </summary>
        public static T Curl2D<T>(this VectorField<T> Vec, CellMask em = null) where T : DGField {
            T Curl = Vec[0].CloneAs() as T;
            Curl.Clear();
            Curl.Clear();
            Curl.Curl2D(1.0, Vec, em);

            return Curl;
        }

        /// <summary>
        /// Central Difference 2D Curl, <see cref="DGField.Curl2DByFlux{T}"/>
        /// </summary>
        public static T Curl2DByFlux<T>(this IEnumerable<T> Vec, SubGrid optionalSubGrid = null, SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) where T : DGField {
            return (new VectorField<T>(Vec.ToArray())).Curl2DByFlux(optionalSubGrid, bndMode);
        }

        /// <summary>
        /// Broken 2D Curl, <see cref="DGField.Curl2D{T}"/>
        /// </summary>
        public static T Curl2D<T>(this IEnumerable<T> Vec, CellMask em = null) where T : DGField {
            return (new VectorField<T>(Vec.ToArray())).Curl2D(em);
        }

        /// <summary>
        /// Central Difference Derivative, <see cref="DGField.Derivative"/>
        /// </summary>
        public static T Derivative<T>(this T f, int d, CellMask em = null) where T : DGField {
            var R = f.CloneAs() as T;
            R.Clear();
            R.Identification = "∂" + (f.Identification??"_") + "/∂x" + d;
            R.Derivative(1.0, f, d, em);
            return R;
        }

        /// <summary>
        /// Central Difference Derivative, <see cref="DGField.DerivativeByFlux"/>
        /// </summary>
        public static T DerivativeByFlux<T>(this T f, int d, SubGrid optionalSubGrid = null, SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) where T : DGField {
            var R = f.CloneAs() as T;
            R.Clear();
            R.Identification = "∂" + (f.Identification??"_") + "/∂x" + d;
            R.DerivativeByFlux(1.0, f, d, optionalSubGrid, bndMode);
            return R;
        }

      
    }
}
