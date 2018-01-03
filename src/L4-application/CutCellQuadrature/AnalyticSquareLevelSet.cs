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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CutCellQuadrature {

    public class AnalyticSquareLevelSet : ILevelSet {

        private GridData gridData;

        private double offsetX;

        private double offsetY;

        public Object Clone() {
            return new AnalyticSquareLevelSet(this.gridData) {
                offsetX = this.offsetX,
                offsetY = this.offsetY
            };
        }

        public AnalyticSquareLevelSet(GridData gridData) {
            this.gridData = gridData;
        }

        public void SetOffset(double offsetX, double offsetY) {
            this.offsetX = offsetX;
            this.offsetY = offsetY;
        }

        private MultidimensionalArray GetGlobalNodes(int j0, int Len, NodeSet nodeSet) {
            int noOfNodes = nodeSet.GetLength(0);
            int D = nodeSet.GetLength(1);

            MultidimensionalArray nodes = MultidimensionalArray.Create(Len, noOfNodes, D);
            gridData.TransformLocal2Global(nodeSet, j0, Len, nodes, 0);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    nodes[i, j, 0] -= offsetX;
                    nodes[i, j, 1] -= offsetY;
                }
            }

            return nodes;
        }

        private static double L2Distance(double x1, double y1, double x2, double y2) {
            return Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        }

        private static double LInfinityNorm(double x, double y) {
            return Math.Max(Math.Abs(x), Math.Abs(y));
        }

        private static double L1Norm(double x, double y) {
            return Math.Abs(x) + Math.Abs(y);
        }

        private static double c = Math.Sqrt(2.0) / 2.0;

        #region ILevelSet Members

        public void Evaluate(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            MultidimensionalArray nodes = GetGlobalNodes(j0, Len, NS);
            int noOfNodes = nodes.GetLength(1);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    double x = nodes[i, j, 0];
                    double y = nodes[i, j, 1];

                    if (y > x + c && y > -x + c) {
                        result[i, j] = -L2Distance(x, y, 0, c);
                    }
                    else if (y > x + c && y < -x - c) {
                        result[i, j] = -L2Distance(x, y, -c, 0);
                    }
                    else if (y < -x - c && y < x - c) {
                        result[i, j] = -L2Distance(x, y, 0, -c);
                    }
                    else if (y < x - c && y > -x + c) {
                        result[i, j] = -L2Distance(x, y, c, 0);
                    }
                    else {
                        result[i, j] = (c - L1Norm(x, y)) * c;
                    }
                }
            }
        }

        public void EvaluateGradient(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            MultidimensionalArray nodes = GetGlobalNodes(j0, Len, NS);
            int noOfNodes = nodes.GetLength(1);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    double x = nodes[i, j, 0];
                    double y = nodes[i, j, 1];

                    if (y > x + c && y > -x + c) {
                        result[i, j, 0] = -(x + c) / L2Distance(x, y, -c, -c);
                        result[i, j, 1] = -(y + c) / L2Distance(x, y, -c, -c);
                    }
                    else if (y > x + c && y < -x - c) {
                        result[i, j, 0] = -(x + c) / L2Distance(x, y, -c, c);
                        result[i, j, 1] = -(y - c) / L2Distance(x, y, -c, c);
                    }
                    else if (y < -x - c && y < x - c) {
                        result[i, j, 0] = -(x - c) / L2Distance(x, y, c, -c);
                        result[i, j, 1] = -(y + c) / L2Distance(x, y, c, -c);
                    }
                    else if (y < x - c && y > -x + c) {
                        result[i, j, 0] = -(x - c) / L2Distance(x, y, c, c);
                        result[i, j, 1] = -(y - c) / L2Distance(x, y, c, c);
                    }
                    else {
                        result[i, j, 0] = Math.Sign(x) * c;
                        result[i, j, 1] = Math.Sign(y) * c;
                    }
                }
            }
        }

        public void EvaluateHessian(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
            MultidimensionalArray nodes = GetGlobalNodes(j0, Len, Ns);
            int noOfNodes = nodes.GetLength(1);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    double x = nodes[i, j, 0];
                    double y = nodes[i, j, 1];

                    if (y > x + c && y > -x + c) {
                        double d = Math.Pow(L2Distance(x, y, -c, -c), 3);
                        result[i, j, 0, 0] = -(y + c) * (y + c) / d;
                        result[i, j, 1, 0] = (x + c) * (y + c) / d;
                        result[i, j, 0, 1] = result[i, j, 1, 0];
                        result[i, j, 1, 1] = -(x + c) * (x + c) / d;
                    }
                    else if (y > x + c && y < -x - c) {
                        double d = Math.Pow(L2Distance(x, y, -c, c), 3);
                        result[i, j, 0, 0] = -(y - c) * (y - c) / d;
                        result[i, j, 1, 0] = (x + c) * (y - c) / d;
                        result[i, j, 0, 1] = result[i, j, 1, 0];
                        result[i, j, 1, 1] = -(x + c) * (x + c) / d;
                    }
                    else if (y < -x - c && y < x - c) {
                        double d = Math.Pow(L2Distance(x, y, c, -c), 3);
                        result[i, j, 0, 0] = -(y + c) * (y + c) / d;
                        result[i, j, 1, 0] = (x - c) * (y + c) / d;
                        result[i, j, 0, 1] = result[i, j, 1, 0];
                        result[i, j, 1, 1] = -(x - c) * (x - c) / d;
                    }
                    else if (y < x - c && y > -x + c) {
                        double d = Math.Pow(L2Distance(x, y, c, c), 3);
                        result[i, j, 0, 0] = -(y - c) * (y - c) / d;
                        result[i, j, 1, 0] = (x - c) * (y - c) / d;
                        result[i, j, 0, 1] = result[i, j, 1, 0];
                        result[i, j, 1, 1] = -(x - c) * (x - c) / d;
                    }
                    else {
                        result[i, j, 0, 0] = 0.0;
                        result[i, j, 1, 0] = 0.0;
                        result[i, j, 0, 1] = 0.0;
                        result[i, j, 1, 1] = 0.0;
                    }
                }
            }
        }

        public void EvaluateTotalCurvature(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
            MultidimensionalArray nodes = GetGlobalNodes(j0, Len, Ns);
            int noOfNodes = nodes.GetLength(1);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    double x = nodes[i, j, 0];
                    double y = nodes[i, j, 1];

                    if (x == y) {
                        result[i, j] = double.PositiveInfinity;
                    }
                    else {
                        result[i, j] = 0.0;
                    }
                }
            }
        }
        #endregion
    }
}
