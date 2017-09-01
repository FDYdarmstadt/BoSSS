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
using System.Collections;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Computes the integral over a specific edge of the computational domain
    /// </summary>
    public class EdgeIntegral {

        /// <summary>
        /// The underlying quadrature object used for the evaluation of the
        /// integral
        /// </summary>
        private MyEdgeQuadrature edgeQuadrature;

        /// <summary>
        /// An execution mask which only contains boundary edges with a
        /// specified edge tag
        /// </summary>
        private EdgeMask relevantEdgesMask;

        /// <summary>
        /// Initializes the integral
        /// </summary>
        /// <param name="grdDat">Information about the grid</param>
        /// <param name="edgeTag">
        /// The edge tag of the boundary edges to be considered in the integral
        /// </param>
        /// <param name="flux">
        /// The integrand. To be more specific,
        /// <see cref="INonlinearFlux.BorderEdgeFlux"/> will be evaluated on
        /// every relevant edge to compute the value of the integrand
        /// </param>
        /// <param name="mapping">
        /// The coordinate mapping that maps field values to the arguments of
        /// <paramref name="flux"/> (see
        /// <see cref="IEquationComponent.ArgumentOrdering"/>)
        /// </param>
        /// <param name="integrationOrder">
        /// The desired order of the applied quadrature rule
        /// </param>
        public EdgeIntegral(GridData grdDat, byte edgeTag, INonlinearFlux flux, CoordinateMapping mapping, int integrationOrder) {
            byte[] edgeTags = grdDat.Edges.EdgeTags;
            BitArray mask = new BitArray(edgeTags.Length);

            int numberOfRelevantEdges = 0;
            for (int i = 0; i < edgeTags.Length; i++) {
                if (edgeTag == edgeTags[i]) {
                    mask[i] = true;
                    numberOfRelevantEdges++;
                }
            }

            relevantEdgesMask = new EdgeMask(grdDat, mask);
            // Just use first ref element for the moment
            edgeQuadrature = new MyEdgeQuadrature(grdDat, integrationOrder, flux, mapping, 0, relevantEdgesMask);
        }

        /// <summary>
        /// Variant of
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// where the integration order is set to the maximum order of any
        /// <see cref="DGField"/> in <paramref name="mapping"/>.
        /// </summary>
        /// <param name="grdDat">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="edgeTag">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="flux">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="mapping">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        public EdgeIntegral(GridData grdDat, byte edgeTag, INonlinearFlux flux, CoordinateMapping mapping)
            : this(grdDat, edgeTag, flux, mapping, mapping.BasisS.Max(basis => basis.Degree)) {
        }

        /// <summary>
        /// Variant of 
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping)"/>
        /// where the edge tag is derived from <paramref name="edgeTagName"/>
        /// </summary>
        /// <param name="grdDat">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="edgeTagName">
        /// The name of the considered edge
        /// </param>
        /// <param name="flux">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="mapping">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        public EdgeIntegral(GridData grdDat, string edgeTagName, INonlinearFlux flux, CoordinateMapping mapping)
            : this(grdDat, GetEdgeTag(grdDat, edgeTagName), flux, mapping, mapping.BasisS.Max(basis => basis.Degree)) {
        }

        /// <summary>
        /// Variant of
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// where the edge tag is derived from <paramref name="edgeTagName"/>
        /// </summary>
        /// <param name="grdDat">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="edgeTagName">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="flux">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="mapping">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        /// <param name="integrationOrder">
        /// <see cref="EdgeIntegral(GridData, byte, INonlinearFlux, CoordinateMapping, int)"/>
        /// </param>
        public EdgeIntegral(GridData grdDat, string edgeTagName, INonlinearFlux flux, CoordinateMapping mapping, int integrationOrder)
            : this(grdDat, GetEdgeTag(grdDat, edgeTagName), flux, mapping, integrationOrder) {
        }

        /// <summary>
        /// Evaluates the integral
        /// </summary>
        /// <returns>
        /// The integral over the considered edges using the current field
        /// values. Involves MPI Communication
        /// </returns>
        public double Evaluate() {
            edgeQuadrature.Result = 0.0;
            edgeQuadrature.Execute();
            double localResult = edgeQuadrature.Result;
            double globalResult;
            unsafe
            {
                csMPI.Raw.Allreduce(
                    (IntPtr)(&localResult),
                    (IntPtr)(&globalResult),
                    1,
                    csMPI.Raw._DATATYPE.DOUBLE,
                    csMPI.Raw._OP.SUM,
                    csMPI.Raw._COMM.WORLD);
            }

            return globalResult;
        }

        /// <summary>
        /// Determines the edge tag associated with <paramref name="edgeTagName"/>
        /// </summary>
        /// <param name="grdDat">The omnipresent context</param>
        /// <param name="edgeTagName">
        /// The name of the edge to be considered
        /// </param>
        /// <returns></returns>
        private static byte GetEdgeTag(GridData grdDat, string edgeTagName) {
            try {
                return grdDat.Grid.EdgeTagNames.First(item => item.Value.Equals(edgeTagName)).Key;
            } catch (InvalidOperationException e) {
                throw new ArgumentException("An edge tag with name \"" + edgeTagName + "\" does not exist", "edgeTagName", e);
            }
        }

        /// <summary>
        /// Quadrature over a given edge of a given flux function
        /// </summary>
        private class MyEdgeQuadrature : EdgeQuadrature {

            
            /// <summary>
            /// The fields that need to be evaluated for the evaluation of
            /// <see cref="flux"/>
            /// </summary>
            private DGField[] evaluators;

            /// <summary>
            /// A flux function representing the integrand
            /// </summary>
            private INonlinearFlux flux;

            /// <summary>
            /// The result of the integration
            /// </summary>
            public double Result;

            /// <summary>
            /// Initializes the quadrature
            /// </summary>
            /// <param name="grdDat">
            /// The omnipresent Grid
            /// </param>
            /// <param name="quadratureOrder">
            /// The order of the applied quadrature rule
            /// </param>
            /// <param name="flux">
            /// A flux function representing the integrand
            /// </param>
            /// <param name="mapping">
            /// The coordinate mapping for the arguments of
            /// <paramref name="flux"/>.
            /// </param>
            /// <param name="iKref">
            /// Reference element index
            /// </param>
            /// <param name="edgeMask">
            /// An optional restriction of the domain
            /// </param>
            public MyEdgeQuadrature(GridData grdDat, int quadratureOrder, INonlinearFlux flux, CoordinateMapping mapping, int iKref, EdgeMask edgeMask)
                : base(new int[] { 1 }, grdDat, (new EdgeQuadratureScheme(true, edgeMask)).Compile(grdDat, quadratureOrder)) {
                evaluators = mapping.Fields.ToArray();
                this.flux = flux;
            }


            /// <summary>
            /// Evaluates the integrand by applying <see cref="flux"/>.
            /// </summary>
            protected override void Evaluate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = QR.Nodes;
                int D = gridData.SpatialDimension;
                int noOfFields = evaluators.Length;
                int NoOfNodes = QuadNodes.NoOfNodes;

                MultidimensionalArray fluxValues = MultidimensionalArray.Create(Length, NoOfNodes);

                MultidimensionalArray[] fieldValues = new MultidimensionalArray[evaluators.Length];
                for (int i = 0; i < evaluators.Length; i++) {
                    fieldValues[i] = MultidimensionalArray.Create(Length, NoOfNodes);
                }

                
                MultidimensionalArray globalCoordinates = MultidimensionalArray.Create(Length, NoOfNodes, D);

                int noOfNodes = QR.NoOfNodes;
                
                //MultidimensionalArray normals = MultidimensionalArray.Create(Length, noOfNodes, D);
                //gridData.Edges.GetNormals(i0, Length, m_NodeSetFamily[0].NodeSet, normals);
                var Normals = gridData.iGeomEdges.NormalsCache.GetNormals_Edge(QR.Nodes, i0, Length);


                // Loop over edges
                for (int i = 0; i < Length; i++) {
                    int cell = gridData.iGeomEdges.CellIndices[i + i0, 0];
                    int edge = gridData.iGeomEdges.FaceIndices[i + i0, 0];
                    int iTrf = gridData.iGeomEdges.Edge2CellTrafoIndex[i + i0, 0];

                    NodeSet CellNodes = QR.Nodes.GetVolumeNodeSet(base.GridDat, iTrf);

                    // Evaluate all fields
                    for (int j = 0; j < noOfFields; j++) {
                        evaluators[j].Evaluate(cell, 1, CellNodes, fieldValues[j], i, 0.0);
                    }

                    // Compute global coordinates
                    gridData.TransformLocal2Global(CellNodes, cell, 1, globalCoordinates, i);

                    // Evaluate integrand
                    flux.BorderEdgeFlux(
                        0.0,
                        i0 + i,
                        globalCoordinates,
                        Normals,
                        false,
                        gridData.iGeomEdges.EdgeTags,
                        i0 + i,
                        fieldValues,
                        i,
                        1,
                        fluxValues);

                    // Save results
                    for (int j = 0; j < NoOfNodes; j++) {
                        EvalResult[i, j, 0] = fluxValues[i, j];
                    }
                }
                
                
            }

            /// <summary>
            /// Sums up the integrals over the respective edges
            /// </summary>
            /// <param name="i0">
            /// <see cref="Quadrature{U,V}.SaveIntegrationResults"/>
            /// </param>
            /// <param name="Length">
            /// <see cref="Quadrature{U,V}.SaveIntegrationResults"/>
            /// </param>
            /// <param name="ResultsOfIntegration">
            /// <see cref="Quadrature{U,V}.SaveIntegrationResults"/>
            /// </param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++) {
                    Result += ResultsOfIntegration[i, 0];
                }
            }
        }

        /// <summary>
        /// Utility class that can be used instead of implementing
        /// <see cref="INonlinearFlux"/> by hand. Here, only
        /// <see cref="NonlinearFlux.BorderEdgeFlux(double, double[], double[], byte, double[], int)"/>
        /// and
        /// <see cref="NonlinearFlux.ArgumentOrdering"/> have to be overridden.
        /// </summary>
        public abstract class EdgeFlux : NonlinearFlux {

            /// <summary>
            /// Not implemented, not intended for overriding;
            /// </summary>
            protected sealed override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
                throw new NotImplementedException();
            }

            /// <summary>
            /// Not implemented, not intended for overriding;
            /// </summary>
            protected sealed override void Flux(double time, double[] x, double[] U, double[] output) {
                throw new NotImplementedException();
            }
        }

    }
}
