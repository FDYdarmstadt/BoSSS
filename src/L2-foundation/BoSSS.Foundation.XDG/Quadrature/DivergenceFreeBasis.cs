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
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// Represents a basis of vector-valued, divergence-free polynomials.
    /// </summary>
    public partial class DivergenceFreeBasis : PolynomialList {

        ///// <summary>
        ///// Cache for
        ///// <see cref="GetPolynomials(GridData, RefElement, RefElement[], int)"/>.
        ///// </summary>
        //private static Dictionary<RefElement, Polynomial[]> cache =
        //    new Dictionary<RefElement, Polynomial[]>();

        /// <summary>
        /// ref elm.
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }


        /// <summary>
        /// Constructors
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="simplex"></param>
        /// <param name="order"></param>
        public DivergenceFreeBasis(GridData gridData, RefElement simplex, int order)
            : base(GetPolynomials(gridData, simplex, order)) //
        {
            this.RefElement = simplex;
        }

        

        /// <summary>
        /// Retrieves all polynomial for the given <paramref name="element"/>
        /// that must be an element of <paramref name="refElements"/>
        /// </summary>
        /// <param name="g"></param>
        /// <param name="element"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        private static IEnumerable<Polynomial> GetPolynomials(
            GridData g, RefElement element, int p) {

            //var ret = new IEnumerable<Polynomial>[refElements.Length];

            int D = element.SpatialDimension;

            // hotfix to get the correct spatial dimension on all polynomials:
            var P = GetPolynomials(element, p).ToArray();
            for (int iP = 0; iP < P.Length; iP++) {
                if (P[iP].Coeff.Length == 0) {
                    Debug.Assert(P[iP].Exponents.GetLength(0) == 0);
                    P[iP].Exponents = new int[0, D];
                }
            }

            //for (int i = 0; i < ret.Length; i++) {
            //    if (object.ReferenceEquals(refElements[i], element)) {
            //        ret[i] = P;
            //    } else {
            //        ret[i] = new Polynomial[P.Count()];
            //    }
            //}

            return P;
        }

        /// <summary>
        /// Uses <see cref="GetPolynomials2D"/> or
        /// <see cref="GetPolynomials3D"/> to retrieve the polynomials with
        /// order less than or equal <paramref name="order"/>
        /// </summary>
        /// <param name="simplex"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        private static IEnumerable<Polynomial> GetPolynomials(RefElement simplex, int order) {
            int numberOfPolynomials;
            if (simplex.SpatialDimension == 2) {
                numberOfPolynomials = (order + 2) * (order + 3) / 2 - 1;
            } else if (simplex.SpatialDimension == 3) {
                numberOfPolynomials = (order + 1) * (order + 2) * (2 * order + 9) / 6;
            } else {
                throw new Exception();
            }
            numberOfPolynomials *= simplex.SpatialDimension;

            Polynomial[] allPolys;
            if (simplex is Square) {
                allPolys = PolynomialsSquare;
            } else if (simplex is Triangle) {
                allPolys = PolynomialsTriangle;
            } else if (simplex is Cube) {
                allPolys = PolynomialsCube;
            } else if (simplex is Tetra) {
                allPolys = PolynomialsTetra;
            } else {
                throw new Exception();
            }

            if (numberOfPolynomials > allPolys.Length)
                throw new NotSupportedException("Divergence-free basis for reference element '" + simplex.GetType().ToString() + "' is not specified for degree " + order + ", max. supported degree is " + (allPolys.Max(poly => poly.AbsoluteDegree)) + ".");

#if DEBUG
            for (int i = 0; i < Math.Min(allPolys.Length, numberOfPolynomials); i++) {
                Debug.Assert(allPolys[i].AbsoluteDegree <= order);
            }
            int D = simplex.SpatialDimension;
            Debug.Assert(allPolys.Length % D == 0);
            for (int i = numberOfPolynomials/D; i < allPolys.Length/D; i++) {
                int DegMax = 0;
                for (int d = 0; d < D; d++) {
                    DegMax = Math.Max(DegMax, allPolys[i * D + d].AbsoluteDegree);
                }
                Debug.Assert(DegMax > order);
            }
#endif

            Polynomial[] ret = new Polynomial[numberOfPolynomials];
            Array.Copy(allPolys, 0, ret, 0, numberOfPolynomials);
            return ret;
        }
    }
}
