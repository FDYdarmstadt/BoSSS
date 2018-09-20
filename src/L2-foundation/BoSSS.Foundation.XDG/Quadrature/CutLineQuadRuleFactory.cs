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
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;
using System.Linq;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// Different type of jumps assumed.
    /// </summary>
    public enum JumpTypes {

        /// <summary>
        /// Function is taken as is; jump is already contained in the function
        /// </summary>
        Implicit,

        /// <summary>
        /// Function is assumed zero for negative level set
        /// </summary>
        Heaviside,

        /// <summary>
        /// Function is assumed zero for positive level set
        /// </summary>
        OneMinusHeaviside,

        /// <summary>
        /// Function is multiplied by -1 for negative level and 1 for positive
        /// level set.
        /// </summary>
        Sign
    }

    /// <summary>
    /// Quadrature rule factory for Gaussian integrals over edges of
    /// two-dimensional domains that are intersected by the level set (i.e.,
    /// integrals over cut segments)
    /// </summary>
    public class CutLineQuadRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule>, IObserver<LevelSetTracker.LevelSetRegions> {

        /// <summary>
        /// The line reference element.
        /// </summary>
        private static readonly Line lineSimplex = Line.Instance;

        ///// <summary>
        ///// Tracks the level set location
        ///// </summary>
        //private LevelSetTracker tracker;

        /// <summary>
        /// Index of the considered level set.
        /// </summary>
        private int levelSetIndex;

        /// <summary>
        /// A collection of line segments bounding the volume reference element
        /// </summary>
        private LineSegment[] referenceLineSegments;

        /// <summary>
        /// Determines the domain of integration, i.e. whether the rule should
        /// be valid for positive level set values, negative level set values,
        /// or even for both
        /// </summary>
        private JumpTypes jumpType;

        /// <summary>
        /// Stores the integration of the set of quadrature rule request last
        /// (cf. <see cref="GetQuadRuleSet"/>)
        /// </summary>
        private int lastOrder = -1;

        /// <summary>
        /// Cache for quadrature rules
        /// <list type="bullet">
        ///     <item>key: local cell index</item>
        ///     <item>value: cell boundary rule</item>
        /// </list>
        /// </summary>
        private Dictionary<int, CellBoundaryQuadRule> cache =
            new Dictionary<int, CellBoundaryQuadRule>();

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return cache.Keys.ToArray();
        }

        /// <summary>
        /// Tolerance for the sign of the level set function. Using the example
        /// of <see cref="JumpTypes.Heaviside"/>, values phi with
        /// phi - Tolerance > 0 are considered 'positive'.
        /// </summary>
        public double Tolerance {
            get;
            private set;
        }

        /// <summary>
        /// Constructs a new factory
        /// </summary>
        /// <param name="tracker">
        /// Tracker for the considered level set
        /// </param>
        /// <param name="Kref">
        /// Reference element index
        /// </param>
        /// <param name="lsData"></param>
        /// <param name="rootFindingAlgorithm">
        /// Selected root-finding algorithm for the line segments. Default is
        /// <see cref="LineSegment.DefaultRootFindingAlgorithm"/>
        /// </param>
        /// <param name="jumpType">
        /// Determines the domain of integration, i.e. whether the rule should
        /// be valid for positive level set values, negative level set values,
        /// or even for both
        /// </param>
        /// <param name="tolerace">
        /// Tolerance for the sign of the level set function. Using the example
        /// of <see cref="JumpTypes.Heaviside"/>, values phi with
        /// phi - Tolerance > 0 are considered 'positive'.
        /// </param>
        public CutLineQuadRuleFactory(
            LevelSetTracker.LevelSetData lsData,
            RefElement Kref,
            LineSegment.IRootFindingAlgorithm rootFindingAlgorithm = null,
            JumpTypes jumpType = JumpTypes.Heaviside,
            double tolerace = 1.0e-13) {

            this.RefElement = Kref;
            this.RootFindingAlgorithm = rootFindingAlgorithm ?? LineSegment.DefaultRootFindingAlgorithm;
            this.referenceLineSegments = GetReferenceLineSegments();
            this.jumpType = jumpType;
            this.levelSetData = lsData;

            this.iKref = lsData.GridDat.Grid.RefElements.IndexOf(this.RefElement, (A, B) => object.ReferenceEquals(A, B));
            if (this.iKref < 0)
                throw new ArgumentException("Reference element cannot be found in the provided grid.");

            
            this.levelSetIndex = lsData.LevelSetIndex;
            if (tolerace < 0.0)
                throw new ArgumentOutOfRangeException();
            this.Tolerance = tolerace;

            emptyrule = CellBoundaryQuadRule.CreateEmpty(this.RefElement, 1,  levelSetData.GridDat.SpatialDimension, referenceLineSegments.Length);
                        // create a rule with just one node and weight zero;
                        // this should avoid some special-case handling for empty rules
            emptyrule.NumbersOfNodesPerFace[0] = 1;
            emptyrule.Nodes.LockForever();
            
            //tracker.Subscribe(this);
        }

        CellBoundaryQuadRule emptyrule;

        LevelSetTracker.LevelSetData levelSetData;

        /// <summary>
        /// The selected root-finding algorithm
        /// </summary>
        public LineSegment.IRootFindingAlgorithm RootFindingAlgorithm {
            get;
            private set;
        }

        #region IQuadRuleFactory<CellBoundaryQuadRule> Members

        /// <summary>
        /// Returns a <see cref="CellBoundaryQuadRule"/> for each cell in the
        /// given <paramref name="mask"/>. This rule consists of Gaussian
        /// quadrature rules on each continuous line segment on the edges of a
        /// given cell (where continuous means 'not intersected by the zero
        /// level set'
        /// </summary>
        /// <param name="mask"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (!(mask is CellMask))
                throw new ArgumentException("This works on cell basis, so a volume mask is required.");
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");
            //Console.WriteLine("boundary order: " + order);

            if (mask == null) {
                mask = CellMask.GetFullMask(levelSetData.GridDat);
            }

            if (order != lastOrder) {
                cache.Clear();
            }

            double[] EdgeToVolumeTransformationDeterminants = this.RefElement.FaceTrafoGramianSqrt;

            QuadRule baseRule = lineSimplex.GetQuadratureRule(order);
            int D = levelSetData.GridDat.SpatialDimension;
            var _Cells = levelSetData.GridDat.Cells;

            var result = new List<ChunkRulePair<CellBoundaryQuadRule>>(mask.NoOfItemsLocally);
            foreach (Chunk chunk in mask) {
                for (int i = 0; i < chunk.Len; i++) {
                    int cell = i + chunk.i0;

                    if (cache.ContainsKey(cell)) {
                        Debug.Assert(cache[cell].Nodes.IsLocked, "Quadrule with non-locked nodes in cache.");
                        result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                            Chunk.GetSingleElementChunk(cell), cache[cell]));
                        continue;
                    }

                    List<double[]> nodes = new List<double[]>();
                    List<double> weights = new List<double>();
                    int[] noOfNodesPerEdge = new int[referenceLineSegments.Length];

                    for (int e = 0; e < referenceLineSegments.Length; e++) {
                        LineSegment referenceSegment = referenceLineSegments[e];
                        int iKref = _Cells.GetRefElementIndex(cell);
                        double[] roots = referenceSegment.GetRoots(levelSetData.LevelSet, cell, iKref);
                        double edgeDet = EdgeToVolumeTransformationDeterminants[e];

                        LineSegment[] subSegments = referenceSegment.Split(roots);

                        for (int k = 0; k < subSegments.Length; k++) {
                            // Evaluate sub segment at center to determine sign
                            NodeSet _point = new NodeSet(this.RefElement, subSegments[k].GetPointOnSegment(0.0));
                            
                            double weightFactor = subSegments[k].Length / referenceSegment.Length;

                            if (weightFactor < 1.0e-14)
                                // segment has a length of approximately 0.0 => no need to care about it.
                                continue;

                            weightFactor *= edgeDet;

                            if(jumpType != JumpTypes.Implicit) {
                                //using (tracker.GridDat.NSC.CreateLock(MultidimensionalArray.CreateWrapper(point, 1, D), this.iKref, -1.0)) {
                                MultidimensionalArray levelSetValue = this.levelSetData.GetLevSetValues(_point, cell, 1);

                                switch(jumpType) {
                                    case JumpTypes.Heaviside:
                                    if(levelSetValue[0, 0] <= -Tolerance) {
                                        continue;
                                    }
                                    break;

                                    case JumpTypes.OneMinusHeaviside:
                                    if(levelSetValue[0, 0] >= Tolerance) {
                                        continue;
                                    }
                                    break;

                                    case JumpTypes.Sign:
                                    weightFactor *= levelSetValue[0, 0].Sign();
                                    break;

                                    default:
                                    throw new NotImplementedException();
                                }
                                //}
                            }

                            for (int m = 0; m < baseRule.NoOfNodes; m++) {
                                // Base rule _always_ is a line rule, thus Nodes[*, _0_]
                                double[] point = subSegments[k].GetPointOnSegment(baseRule.Nodes[m, 0]);

                                weights.Add(weightFactor * baseRule.Weights[m]);
                                nodes.Add(point);

                                noOfNodesPerEdge[e]++;
                            }
                        }
                    }

                    if (weights.Count == 0) {
                        

                        result.Add(new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(cell), emptyrule));
                        continue;
                    }

                    NodeSet localNodes = new NodeSet(this.RefElement, nodes.Count, D);
                    for (int j = 0; j < nodes.Count; j++) {
                        for (int d = 0; d < D; d++) {
                            localNodes[j, d] = nodes[j][d];
                        }
                    }
                    localNodes.LockForever();

                    CellBoundaryQuadRule subdividedRule = new CellBoundaryQuadRule() {
                        OrderOfPrecision = order,
                        Weights = MultidimensionalArray.Create(weights.Count),
                        Nodes = localNodes,
                        NumbersOfNodesPerFace = noOfNodesPerEdge
                    };
                    subdividedRule.Weights.SetSubVector(weights, -1);

                    cache.Add(cell, subdividedRule);

                    result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                        Chunk.GetSingleElementChunk(cell), subdividedRule));
                }
            }

            return result;
        }

        /// <summary>
        /// Reference element index
        /// </summary>
        private int iKref;

        /// <summary>
        /// Selected reference element
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        #endregion

        /// <summary>
        /// Determines the line segments bounding the selected reference
        /// element (cf. <see cref="RefElement"/>) in reference coordinates
        /// </summary>
        /// <returns></returns>
        private LineSegment[] GetReferenceLineSegments() {
            Stack<RefElement> simplexHierarchy = new Stack<RefElement>();
            RefElement currentSimplex = RefElement;
            int spatialDimension = RefElement.SpatialDimension;

            int n = 0;
            while (currentSimplex.GetType() != lineSimplex.GetType()) {
                // If n > 2, the edge of the edge of a simplex is not a line.
                // This is hardly likely in 2d/3d.
                if (n > 2) {
                    throw new ApplicationException("Something went terribly wrong.");
                }

                simplexHierarchy.Push(currentSimplex);

                currentSimplex = currentSimplex.FaceRefElement;
                n++;
            }

            MultidimensionalArray vertexCoordinates = MultidimensionalArray.Create(2, 1);
            vertexCoordinates[0, 0] = -1.0;
            vertexCoordinates[1, 0] = 1.0;

            while (simplexHierarchy.Count > 0) {
                currentSimplex = simplexHierarchy.Pop();

                int noOfVertices = vertexCoordinates.GetLength(0);
                int D = currentSimplex.SpatialDimension;
                MultidimensionalArray volumeCoordinates = MultidimensionalArray.Create(
                    noOfVertices * currentSimplex.NoOfFaces, currentSimplex.SpatialDimension);

                for (int e = 0; e < currentSimplex.NoOfFaces; e++) {
                    MultidimensionalArray coordinates = MultidimensionalArray.Create(noOfVertices, D);
                    currentSimplex.TransformFaceCoordinates(e, vertexCoordinates, coordinates);

                    for (int i = 0; i < noOfVertices; i++) {
                        for (int d = 0; d < D; d++) {
                            volumeCoordinates[e * noOfVertices + i, d] = coordinates[i, d];
                        }
                    }
                }

                vertexCoordinates = volumeCoordinates;
            }

            Debug.Assert(
                vertexCoordinates.GetLength(0) % 2 == 0,
                "Even number of vertices expected");
            int initialNumberOfLineSegments = vertexCoordinates.GetLength(0) / 2;

            List<LineSegment> lineSegments = new List<LineSegment>(initialNumberOfLineSegments);
            for (int i = 0; i < initialNumberOfLineSegments; i++) {

                var p0 = vertexCoordinates.GetRow(2 * i + 0);
                int iP0 = this.RefElement.Vertices.FindRow(p0, 1.0e-8);
                var p1 = vertexCoordinates.GetRow(2 * i + 1);
                int iP1 = this.RefElement.Vertices.FindRow(p1, 1.0e-8);

                LineSegment newSegment = new LineSegment(spatialDimension, this.RefElement,
                    p0, p1, iP0, iP1,
                    RootFindingAlgorithm);

                if (!lineSegments.Contains(newSegment)) {
                    lineSegments.Add(newSegment);
                    //tracker.Subscribe(newSegment);
                }
            }

            foreach (LineSegment segment in lineSegments) {
                LevelSet levelSetField = levelSetData.LevelSet as LevelSet;
                if (levelSetField != null) {
                    segment.ProjectBasisPolynomials(levelSetField.Basis);
                }
            }

            return lineSegments.ToArray();
        }

        #region IObserver<LevelSetData> Members

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void OnCompleted() {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="error"></param>
        public void OnError(Exception error) {
        }

        /// <summary>
        /// Clears the cache.
        /// </summary>
        /// <param name="value"></param>
        public void OnNext(LevelSetTracker.LevelSetRegions value) {
            cache.Clear();
        }

        #endregion
    }
}
