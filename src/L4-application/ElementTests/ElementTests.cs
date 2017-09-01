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
using System.Runtime.CompilerServices;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using NUnit.Framework;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Application.ElementTests {

    [TestFixture]
    class ElementTests {

        /// <summary>
        /// Ensures bootstrapping has been conducted correctly.
        /// </summary>
        [TestFixtureSetUp]
        public void SetUp() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
        }

        [Test]
        [MethodImpl(MethodImplOptions.NoOptimization)]
        public void Bla([Values(0, 1, 2, 3)] int refElementIndex, [Values(0, 1, 2)] int foreignTypeIndex) {
            RefElement Element = ListElementsMain.Elements[refElementIndex];
            RefElement.ExchangeFormats ft = ListElementsMain.ForeignTypes[foreignTypeIndex];

            foreach (var type in Element.SupportedCellTypes) {
                try {
                    int fnum;
                    string fname;
                    Element.GetForeignElementType(type, ft, out fname, out fnum);
                } catch (NotSupportedException) {
                    // Combination of element type and foreign convection
                    // simply does not exist, swallow exception.
                }

                // Check any of these methods throws an exception
                MultidimensionalArray InterpolationNodes = Element.GetInterpolationNodes(type);
                int[] NodeType = Element.GetInterpolationNodes_NodeType(type);
                int[] EntityIndex = Element.GetInterpolationNodes_EntityIndices(type);
                int NoOfNodes = InterpolationNodes.GetLength(0);
                var R = Element.GetForeignElementMapping(type, ft);
            }
        }
    }
}
