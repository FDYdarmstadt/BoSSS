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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CutCellQuadrature {

    class AnalyticStarLevelSet : IAnalyticLevelSet {

        private GridData gridData;

        private double offsetX;

        private double offsetY;

        public Object Clone() {
            return new AnalyticStarLevelSet(this.gridData) {
                offsetX = this.offsetX,
                offsetY = this.offsetY
            };
        }


        public AnalyticStarLevelSet(GridData gridData) {
            this.gridData = gridData;
        }

        public void SetParameters(double offsetX, double offsetY) {
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

        #region IAnalyticLevelSet Members

        public IEnumerable<double> GetRoots(LineSegment lineSegment, int element) {
            double left = -1.0;
            double right = 1.0;
            MultidimensionalArray levelSetValue = MultidimensionalArray.Create(1, 1);

            double[] point = lineSegment.GetPointOnSegment(left);
            NodeSet pointAsArray = new NodeSet(gridData.Cells.GetRefElement(element), point);
            Evaluate(element, 1, pointAsArray, levelSetValue);
            double valueLeft = levelSetValue[0, 0];

            point = lineSegment.GetPointOnSegment(right);
            pointAsArray = new NodeSet(gridData.Cells.GetRefElement(element), point);
            Evaluate(element, 1, pointAsArray, levelSetValue);
            double valueRight = levelSetValue[0, 0];

            if (valueLeft.Sign() == valueRight.Sign()) {
                yield break;
            }

            if (valueLeft > valueRight) {
                double temp = valueRight;
                valueRight = valueLeft;
                valueLeft = temp;

                left = 1.0;
                right = -1.0;
            }

            while (true) {
                double center = 0.5 * (left + right);
                point = lineSegment.GetPointOnSegment(center);
                pointAsArray = new NodeSet(gridData.Cells.GetRefElement(element), point);
            
                Evaluate(element, 1, pointAsArray, levelSetValue);

                double valueCenter = levelSetValue[0, 0];

                if (valueCenter.Abs() < 1e-8) {
                    yield return center;
                    break;
                } else if (valueCenter > 0.0) {
                    right = center;
                    valueRight = valueCenter;
                } else if (valueCenter < 0.0) {
                    left = center;
                    valueLeft = center;
                }
            }
        }

        #endregion

        #region ILevelSet Members

        public void Evaluate(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
            MultidimensionalArray nodes = GetGlobalNodes(j0, Len, Ns);
            int noOfNodes = nodes.GetLength(1);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    double x = nodes[i, j, 0];
                    double y = nodes[i, j, 1];

                    double r = Math.Sqrt(x * x + y * y);
                    double theta = Math.Atan2(y, x);
                    double value = 1.0 + 0.5 * Math.Cos(5.0 * theta) - r;
                    if (value > 0.9 || r < 0.05) {
                        value = 0.9;
                    }
                    result[i, j] = value;
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

                    double theta = Math.Atan2(y, x);
                    double r = Math.Sqrt(x * x + y * y);
                    //if (r < 0.001) {
                    //    result[i, j, 0] = 0.0;
                    //    result[i, j, 1] = 0.0;
                    //} else {
                    result[i, j, 0] = -2.5 * Math.Sin(5.0 * theta) * (theta - x * y / r / r) - x / r;
                    result[i, j, 1] = 2.5 * Math.Sin(5.0 * theta) * y * y / r / r - y / r;
                    //}
                }
            }
        }

        public void EvaluateHessian(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
            throw new NotImplementedException();
        }

        public void EvaluateTotalCurvature(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result) {
            throw new NotImplementedException();
        }

        #endregion
    }
}
