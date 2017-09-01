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

    public class AnalyticSphereLevelSet : IAnalyticLevelSet {

        private GridData gridData;

        private double offsetX;

        private double offsetY;

        private double offsetZ;

        public AnalyticSphereLevelSet(GridData gridData) {
            this.gridData = gridData;
        }

        public void SetOffset(double offsetX, double offsetY, double offsetZ) {
            this.offsetX = offsetX;
            this.offsetY = offsetY;
            this.offsetZ = offsetZ;
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
                    nodes[i, j, 2] -= offsetZ;
                }
            }

            return nodes;
        }

        #region IAnalyticLevelSet Members

        public IEnumerable<double> GetRoots(LineSegment lineSegment, int element) {
            int D = lineSegment.SpatialDimension;

            MultidimensionalArray localSegmentCoordinates = MultidimensionalArray.Create(2, D);
            localSegmentCoordinates[0, 0] = lineSegment.Start[0];
            localSegmentCoordinates[0, 1] = lineSegment.Start[1];
            localSegmentCoordinates[0, 2] = lineSegment.Start[2];
            localSegmentCoordinates[1, 0] = lineSegment.End[0];
            localSegmentCoordinates[1, 1] = lineSegment.End[1];
            localSegmentCoordinates[1, 2] = lineSegment.End[2];

            MultidimensionalArray globalSegmentCoordinates = MultidimensionalArray.Create(1, 2, D);
            gridData.TransformLocal2Global(localSegmentCoordinates, element, 1, globalSegmentCoordinates, 0);

            double[] start = new double[D];
            start[0] = globalSegmentCoordinates[0, 0, 0];
            start[1] = globalSegmentCoordinates[0, 0, 1];
            start[2] = globalSegmentCoordinates[0, 0, 2];

            double[] end = new double[D];
            end[0] = globalSegmentCoordinates[0, 1, 0];
            end[1] = globalSegmentCoordinates[0, 1, 1];
            end[2] = globalSegmentCoordinates[0, 1, 2];

            LineSegment globalSegment = new LineSegment(
                D,
                lineSegment.RefElement,
                start,
                end);



            double p = 0.0;
            double q = 0.0;
            double divisor = 0.0;
            {
                // x
                double a = 0.5 * (globalSegment.End[0] - globalSegment.Start[0]);
                double b = 0.5 * (globalSegment.Start[0] + globalSegment.End[0]) - offsetX;

                p += 2.0 * a * b;
                q += b * b;
                divisor += a * a;

                // y
                a = 0.5 * (globalSegment.End[1] - globalSegment.Start[1]);
                b = 0.5 * (globalSegment.Start[1] + globalSegment.End[1]) - offsetY;

                p += 2.0 * a * b;
                q += b * b;
                divisor += a * a;

                // 7
                a = 0.5 * (globalSegment.End[2] - globalSegment.Start[2]);
                b = 0.5 * (globalSegment.Start[2] + globalSegment.End[2]) - offsetZ;

                p += 2.0 * a * b;
                q += b * b;
                divisor += a * a;
            }

            p /= divisor;
            q = (q - 1.0) / divisor;

            double radicand = 0.25 * p * p - q;
            //double EPS = 1e-10;
            //if (radicand.Abs() < EPS) {
            //    radicand = 0.0;
            //}

            if (radicand < 0.0) {
                yield break;
            } else if (radicand == 0.0) {
                double root = -0.5 * p;
                if (Math.Abs(root) <= 1.0) {
                    yield return root;
                }
            } else {
                double root = -0.5 * p - Math.Sqrt(radicand);
                if (Math.Abs(root) <= 1.0) {
                    yield return root;
                }

                root = -0.5 * p + Math.Sqrt(radicand);
                if (Math.Abs(root) <= 1.0) {
                    yield return root;
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
                    double z = nodes[i, j, 2];

                    result[i, j] = 1.0 - Math.Sqrt(x * x + y * y + z * z);
                }
            }
        }

        public void EvaluateGradient(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
            MultidimensionalArray nodes = GetGlobalNodes(j0, Len, Ns);
            int noOfNodes = nodes.GetLength(1);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    double x = nodes[i, j, 0];
                    double y = nodes[i, j, 1];
                    double z = nodes[i, j, 2];

                    double d = Math.Sqrt(x * x + y * y + z * z);
                    result[i, j, 0] = -x / d;
                    result[i, j, 1] = -y / d;
                    result[i, j, 2] = -z / d;
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
                    double z = nodes[i, j, 2];

                    double d = Math.Pow(x * x + y * y + z * z, 1.5);
                    result[i, j, 0, 0] = (-y * y - z * z) / d;
                    result[i, j, 1, 0] = x * y / d;
                    result[i, j, 2, 0] = x * z / d;

                    result[i, j, 0, 1] = result[i, j, 1, 0];
                    result[i, j, 1, 1] = (-x * x - z * z) / d;
                    result[i, j, 2, 1] = y * z / d;

                    result[i, j, 0, 2] = result[i, j, 2, 0];
                    result[i, j, 1, 2] = result[i, j, 2, 1];
                    result[i, j, 2, 2] = (-x * x - y * y) / d;
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
                    double z = nodes[i, j, 2];

                    result[i, j] = -1.0 / Math.Sqrt(x * x + y * y + z * z);
                }
            }
        }

        #endregion
    }
}
