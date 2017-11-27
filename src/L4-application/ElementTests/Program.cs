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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution;

namespace BoSSS.Application.ElementTests {

    class ListElementsMain : BoSSS.Solution.Application {

        public static RefElement[] Elements = new RefElement[] {
            Line.Instance,
            Square.Instance,
            Triangle.Instance,
            Cube.Instance,
            Tetra.Instance
        };

        public static RefElement.ExchangeFormats[] ForeignTypes =
            new RefElement.ExchangeFormats[] {
                RefElement.ExchangeFormats.Gmsh,
                RefElement.ExchangeFormats.CGNS,
                RefElement.ExchangeFormats.GambitNeutral };

        public static void Main(string[] args) {
            BoSSS.Solution.Application._Main(args, true, null, delegate() {
                return new ListElementsMain();
            });
        }

        public void List() {
            foreach(var Element in Elements) {
                var Types = Element.SupportedCellTypes;

                foreach (var ElmType in Types) {
                    // write cell info
                    Console.WriteLine("{0}, minor = {1}", Element.GetType().Name, ElmType);
                    foreach (var ft in ForeignTypes) {
                        Console.Write("  in {0}: ", ft.ToString());
                        try {
                            int fnum;
                            string fname;
                            Element.GetForeignElementType(ElmType, ft, out fname, out fnum);
                            //Console.Write("number {0}, name '{1}'", fnum, fname);
                        } catch (Exception ee) {
                            Console.Write("'{0}', '{1}'", ee.GetType().Name, ee.Message);
                            //Console.Write("not impl.");
                        }

                        Console.WriteLine();
                    }


                    // write header
                    Console.Write("#\t");
                    string[] xyz = new string[] { "x", "y", "z" };
                    int D = Element.SpatialDimension;
                    for (int d = 0; d < D; d++) {
                        Console.Write("{0}\t", xyz[d]);
                    }
                    Console.Write("N.type\t");
                    Console.Write("Ent.I.\t");
                    foreach (var ft in ForeignTypes) {
                        Console.Write("{0}\t", ft.ToString());
                    }
                    Console.WriteLine();

                    // write node list
                    MultidimensionalArray InterpolationNodes = Element.GetInterpolationNodes(ElmType);
                    int[] NodeType = Element.GetInterpolationNodes_NodeType(ElmType);
                    int[] EntityIndex = Element.GetInterpolationNodes_EntityIndices(ElmType);
                    int NoOfNodes = InterpolationNodes.GetLength(0);
                    var NodeIdxTrans = new Dictionary<RefElement.ExchangeFormats, int[]>();
                    foreach (var ft in ForeignTypes) {
                        try {
                            var R = Element.GetForeignElementMapping(ElmType, ft);
                            NodeIdxTrans.Add(ft, R);
                        } catch (Exception) {
                            NodeIdxTrans.Add(ft, null);
                        }
                    }
                    for (int k = 0; k < NoOfNodes; k++) {
                        Console.Write(k);
                        Console.Write("\t");

                        for (int d = 0; d < D; d++) {
                            Console.Write(InterpolationNodes[k, d].ToString("0.####"));
                            Console.Write("\t");
                        }
                        Console.Write(NodeType[k]);
                        Console.Write("\t");
                        Console.Write(EntityIndex[k]);
                        Console.Write("\t");

                        foreach (var ft in ForeignTypes) {
                            var R = NodeIdxTrans[ft];
                            //if (R == null)
                            //    Console.Write("n.a.\t");
                            //else
                            //    Console.Write("{0}\t", R[k]);
                        }

                        Console.WriteLine();
                    }

                }

                Console.WriteLine("--------------------------");
            }
        }

        protected override GridCommons CreateOrLoadGrid() {
            List();
            return Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, 1, 5), GenericBlas.Linspace(-1, 1, 5));
        }

        protected override void CreateEquationsAndSolvers(GridUpdateData L) {
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            base.TerminationKey = true;
            return 0.0;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
        }
    }
}
