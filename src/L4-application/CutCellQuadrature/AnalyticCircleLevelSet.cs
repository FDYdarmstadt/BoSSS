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

    public class AnalyticCircleLevelSet : IAnalyticLevelSet {

        private GridData gridData;

        private double radius;

        private double offsetX;

        private double offsetY;

        public AnalyticCircleLevelSet(GridData gridData) {
            this.gridData = gridData;
        }

        public void SetParameters(double radius, double offsetX, double offsetY) {
            this.radius = radius;
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

        #region ILevelSet Members

        public void Evaluate(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {

            MultidimensionalArray nodes = GetGlobalNodes(j0, Len, Ns);
            int noOfNodes = nodes.GetLength(1);

            for (int i = 0; i < Len; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    double x = nodes[i, j, 0];
                    double y = nodes[i, j, 1];

                    result[i, j] = radius - Math.Sqrt(x * x + y * y);
                    //result[i, j] = radius - (x * x + y * y);
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

                    double norm = Math.Sqrt(x * x + y * y);

                    // Gradient not defined, just perturb it
                    if (norm == 0) {
                        x += 1e-15;
                        y += 1e-15;
                        norm = Math.Sqrt(x * x + y * y);
                    }

                    result[i, j, 0] = -x / norm;
                    result[i, j, 1] = -y / norm;

                    //result[i, j, 0] = -2.0 * x;
                    //result[i, j, 1] = -2.0 * y;
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

                    double c = Math.Pow(x * x + y * y, -1.5);
                    result[i, j, 0, 0] = -y * y * c;
                    result[i, j, 1, 0] = x * y * c;
                    result[i, j, 0, 1] = result[i, j, 1, 0];
                    result[i, j, 1, 1] = -x * x * c;
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

                    double norm = Math.Sqrt(x * x + y * y);

                    // Small perturbation at the singularity
                    if (norm == 0) {
                        norm = 1e-15;
                    }
                    
                    result[i, j] = -1/norm;
                }
            }
        }
        

        public IEnumerable<double> GetRoots(LineSegment lineSegment, int element) {
            int D = lineSegment.SpatialDimension;

            MultidimensionalArray localSegmentCoordinates = MultidimensionalArray.Create(2, D);
            localSegmentCoordinates[0, 0] = lineSegment.Start[0];
            localSegmentCoordinates[0, 1] = lineSegment.Start[1];
            localSegmentCoordinates[1, 0] = lineSegment.End[0];
            localSegmentCoordinates[1, 1] = lineSegment.End[1];

            MultidimensionalArray globalSegmentCoordinates = MultidimensionalArray.Create(1, 2, D);
            gridData.TransformLocal2Global(localSegmentCoordinates, element, 1, globalSegmentCoordinates, 0);

            double[] start = new double[D];
            double[] end = new double[D];
            start[0] = globalSegmentCoordinates[0, 0, 0];
            start[1] = globalSegmentCoordinates[0, 0, 1];
            end[0] = globalSegmentCoordinates[0, 1, 0];
            end[1] = globalSegmentCoordinates[0, 1, 1];
            LineSegment globalSegment = new LineSegment(D, gridData.Cells.RefElements[0], start, end);

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
            }

            p /= divisor;
            q = (q - radius * radius) / divisor;

            double radicand = 0.25 * p * p - q;
            if (radicand < 0.0) {
                yield break;
            } else if (radicand == 0.0) {
                double root = -0.5 * p;
                if (Math.Abs(root) < 1.0) {
                    yield return root;
                }
            } else {
                double root = -0.5 * p - Math.Sqrt(radicand);
                if (Math.Abs(root) < 1.0) {
                    yield return root;
                }

                root = -0.5 * p + Math.Sqrt(radicand);
                if (Math.Abs(root) < 1.0) {
                    yield return root;
                }
            }
        }

        #endregion
    }
}
