using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG.Quadrature {
    class LevelSetOnEdgeRule {
        GridData grddat;

        LevelSetTracker.LevelSetData levelSetData;

        int D;

        public LevelSetOnEdgeRule(LevelSetTracker.LevelSetData levelSetData) {
            this.levelSetData = levelSetData;
            grddat = levelSetData.GridDat;
            D = grddat.SpatialDimension;
        }

        public bool IsSpecialCell(int j) {
            (int iLevSet, int iFace)[][] CoIncFaces = levelSetData.Region.m_LevSetCoincidingFaces;
            if (CoIncFaces == null)
                return false;
            if (CoIncFaces[j] == null)
                return false;

            foreach (var t in CoIncFaces[j]) {
                int levelSetIndex = levelSetData.LevelSetIndex;
                if (t.iLevSet == levelSetIndex)
                    return true; // jetzt geht der Spass los!
            }
            return false;
        }

        int GetSpecialFaceIndex(int j) {
            (int iLevSet, int iFace)[][] CoIncFaces = levelSetData.Region.m_LevSetCoincidingFaces;
            foreach (var t in CoIncFaces[j]) {
                int levelSetIndex = levelSetData.LevelSetIndex;
                if (t.iLevSet == levelSetIndex)
                    return t.iFace; // jetzt geht der Spass los!
            }
            throw new Exception("Face does not have registered special Face");
        }

        /// <summary>
        /// Call this function when the IsSpecialCell returns true.
        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// special case: level-set coincides with a cell-face; 
        /// cell is either full or empty 
        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// </summary>
        public (ChunkRulePair<QuadRule> surfaceRule, ChunkRulePair<QuadRule> volumeRule) ComboQuadRule(int intOrder, int jCell, int speciesSign = 1) {
            ChunkRulePair<QuadRule> surfaceRule = SurfaceQuadRule(intOrder, jCell);
            ChunkRulePair<QuadRule> volumeRule = VolumeQuadRule(intOrder, jCell, speciesSign);
            return (surfaceRule, volumeRule);
        }

        /// <summary>
        /// Call this function when the IsSpecialCell returns true.
        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// special case: level-set coincides with a cell-face; 
        /// cell is either full or empty 
        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// Note:
        /// Only one of the two cells on the special edge should have a surface rule, the other cell should be empty
        /// Otherwise, the level-set components are integrates twice (or never)!
        ///
        /// The convention which we use is:
        ///   - the cell with the **lower global index** has the full rule
        ///   - the other cell is empty
        /// We use the global index here, so that the result is "stable" even if we are at an MPI boundary.
        /// </summary>
        /// <param name="intOrder"></param>
        /// <param name="jCell"></param>
        /// <returns></returns>
        public ChunkRulePair<QuadRule> SurfaceQuadRule(int intOrder, int jCell) {
            int SpecialFace = GetSpecialFaceIndex(jCell);
            RefElement refElement = grddat.Cells.GetRefElement(jCell);

            // determine whether this cell is 
            // inside or outside 
            // w.r.t. the edge that corresponds with face 
            int iEdge = grddat.Cells.GetEdgesForFace(jCell, SpecialFace, out int InOrOut, out int[] FurtherEdges);
            if (FurtherEdges != null && FurtherEdges.Length > 0) {
                // throw new NotSupportedException("Hanging node on a edge which coincides with the level set - this should be avoided.");
                // Console.WriteLine("Hanging node on a edge which coincides with the level set - this should be avoided.");
            }
            int J = grddat.CellPartitioning.LocalLength;

            // surface rule -- classification, using conformality of edges
            // always use conformal cell, if both conformal choose lower global index
            bool EmptyOrFool; // true: full rule; false: empty rule
            {
                int OtherCell = grddat.Edges.CellIndices[iEdge, InOrOut == 0 ? 1 : 0];
                bool jConform = InOrOut == 0 ? grddat.Edges.IsEdgeConformalWithCell1(iEdge) : grddat.Edges.IsEdgeConformalWithCell2(iEdge);
                bool OtherConform = InOrOut == 0 ? grddat.Edges.IsEdgeConformalWithCell2(iEdge) : grddat.Edges.IsEdgeConformalWithCell1(iEdge);

                // if both cell faces have no hanging nodes use globally lower index
                if (jConform & OtherConform) {
                    long jCellGlob = jCell + grddat.CellPartitioning.i0;
                    long OtherCellGlob = OtherCell < J ? OtherCell + grddat.CellPartitioning.i0 : grddat.Parallel.GlobalIndicesExternalCells[OtherCell - J];
                    EmptyOrFool = jCellGlob < OtherCellGlob;
                } else if (jConform & !OtherConform) {
                    EmptyOrFool = true; // this cells face consists of one edge only, empty rule
                } else if (!jConform & OtherConform) {
                    EmptyOrFool = false; // this cells face consists of more than one edge only, use full rule, so that each partial edge contains the rule of full degree!
                } else {
                    throw new NotSupportedException(String.Format("Error in cell {0}, {1}: Only one cell should have the hanging node.", jCell, OtherCell));
                }
            }

            ChunkRulePair<QuadRule> surfaceRule;
            if (EmptyOrFool) {

                var FaceRule = refElement.FaceRefElement.GetQuadratureRule(intOrder);
                int K = FaceRule.NoOfNodes;
                NodeSet VolumeNodes = new NodeSet(refElement, K, D, true);
                refElement.TransformFaceCoordinates(SpecialFace, FaceRule.Nodes, VolumeNodes);
                VolumeNodes.LockForever();

                double gTrF = refElement.FaceTrafoGramianSqrt[SpecialFace];
                var metrics = this.levelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(VolumeNodes, jCell, 1);

                QuadRule qr_l = new QuadRule() {
                    OrderOfPrecision = intOrder,
                    Weights = MultidimensionalArray.Create(K),
                    Nodes = VolumeNodes
                };

                for (int k = 0; k < K; k++) {
                    qr_l.Weights[k] = FaceRule.Weights[k] * gTrF / metrics[0, k];
                }

                surfaceRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);

            } else {
                // use empty rule
                QuadRule empty = new QuadRule() {
                    OrderOfPrecision = intOrder,
                    Weights = MultidimensionalArray.Create(1),
                    Nodes = refElement.Center
                };

                surfaceRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), empty);
            }
            return surfaceRule;
        }

        /// <summary>
        /// Call this function when the IsSpecialCell returns true.
        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// special case: level-set coincides with a cell-face; 
        /// cell is either full or empty 
        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// Note:
        /// Only one of the two cells on the special edge should have a surface rule, the other cell should be empty
        /// Otherwise, the level-set components are integrates twice (or never)!
        ///
        /// The convention which we use is:
        ///   - the cell with the **lower global index** has the full rule
        ///   - the other cell is empty
        /// We use the global index here, so that the result is "stable" even if we are at an MPI boundary.
        /// </summary>
        public ChunkRulePair<QuadRule> VolumeQuadRule(int intOrder, int jCell, int speciesSign) {
            RefElement refElement = grddat.Cells.GetRefElement(jCell);

            //
            // Volume rule
            // -----------

            ChunkRulePair<QuadRule> volumeRule;
            // we assume that the cell is either completely empty or completely occupied
            bool EmptyOrFoolVol; // true: full rule; false: empty rule
            {
                var LsVals = levelSetData.GetLevSetValues(refElement.Center, jCell, 1);
                EmptyOrFoolVol = LsVals[0, 0] * speciesSign > 0;
            }

            if (EmptyOrFoolVol) {

                var VolRule = refElement.GetQuadratureRule(intOrder);

                volumeRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), VolRule);

            } else {
                // use empty rule
                QuadRule empty = new QuadRule() {
                    OrderOfPrecision = intOrder,
                    Weights = MultidimensionalArray.Create(1),
                    Nodes = refElement.Center
                };

                volumeRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), empty);
            }
            return volumeRule;
        }
    }
}
