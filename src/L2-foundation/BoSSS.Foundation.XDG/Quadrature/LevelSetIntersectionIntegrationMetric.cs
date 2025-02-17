using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG {


    /// <summary>
    /// integration metric for $D-2$-dimensional integral on the intersection of two level-sets. This is a point in 2D and a line in 3D,
    /// i.e., for a cell $K_j$ this is an integral of the type:
    /// \[
    ///    \int_{K_j  \cap \mathfrak{I}_2 \cap \mathfrak{I}_2}  f \mathrm{dl} .
    /// \]
    /// </summary>
    internal class LevelSetIntersectionIntegrationMetric : IIntegrationMetric {


        public LevelSetIntersectionIntegrationMetric(LevelSetTracker.LevelSetData levelSetData1, LevelSetTracker.LevelSetData levelSetData2) {
            m_levelSetData1 = levelSetData1;
            m_levelSetData2 = levelSetData2;

            if(!object.ReferenceEquals(levelSetData1.GridDat, levelSetData2.GridDat))
                throw new ArgumentException("Grid Data mismatch");
            if(object.ReferenceEquals(levelSetData1, levelSetData2))
                throw new ArgumentException("must be two different level-sets");

        }

        readonly LevelSetTracker.LevelSetData m_levelSetData1;
        readonly LevelSetTracker.LevelSetData m_levelSetData2;

        /// <summary>
        /// For the integral over the level-set intersection, 
        /// the metric depends on the node (if the element transformation applies some shearing transformation),
        /// but only in 3D
        /// </summary>
        public bool AlwaysUsePerNodeScaling {
            get {
                if(m_levelSetData1.GridDat.SpatialDimension == 2)
                    return false;
                if(m_levelSetData1.GridDat.SpatialDimension == 3)
                    return true;
                throw new NotSupportedException();
            }
        }


        /// <summary>
        /// in 2D, the level-set intersection are just points -- one needs a count measure, where the integration metric is always 1
        /// </summary>
        public MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int jCell0, int L) {
            if(m_levelSetData1.GridDat.SpatialDimension == 2) {
                var ret = MultidimensionalArray.Create(L);
                ret.SetAll(1.0);
                return ret;
            } else {
                throw new NotSupportedException();
            }

        }

        /// <summary>
        /// In 3D, the integral transformation metric is the local stretching of a tangent.
        /// </summary>
        public MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int jCell0, int L) {
            if(!object.ReferenceEquals(gridData, m_levelSetData1.GridDat))
                throw new ArgumentException();
            if(!object.ReferenceEquals(gridData, m_levelSetData2.GridDat))
                throw new ArgumentException();

            //throw new NotImplementedException("todo");

            int[,] E2C = gridData.iGeomEdges.CellIndices;
            int[,] E2Ctrafo = gridData.iGeomEdges.Edge2CellTrafoIndex;
            byte[,] E2F = gridData.iGeomEdges.FaceIndices;
            int D = gridData.SpatialDimension;
            int K = qr.NoOfNodes;

            var metric = MultidimensionalArray.Create(L, K);

            var LevSetNormals1 = this.m_levelSetData1.GetLevelSetReferenceNormals(qr.Nodes, jCell0, 1).ExtractSubArrayShallow(0, -1, -1);
            var LevSetNormals2 = this.m_levelSetData2.GetLevelSetReferenceNormals(qr.Nodes, jCell0, 1).ExtractSubArrayShallow(0, -1, -1);
            var Tangents = TempBuffer.GetTempMultidimensionalarray(out int iTangengBuffer, L, 2 * K, D);

            for(int j = 0; j < L; j++) {

                var LevSetNormals1_j = LevSetNormals1.ExtractSubArrayShallow(j, -1, -1);
                var LevSetNormals2_j = LevSetNormals1.ExtractSubArrayShallow(j, -1, -1);
                var Tangents_j = Tangents.ExtractSubArrayShallow(j, -1, -1);

                for(int k = 0; k < K; k++) {
                    var SurfNormal1 = LevSetNormals1_j.GetRowPt(k);
                    var SurfNormal2 = LevSetNormals2_j.GetRowPt(k);
                    var Node = qr.Nodes.GetRowPt(k);
                    var Tangent = SurfNormal1.CrossProduct(SurfNormal2);
                    Tangent.NormalizeInPlace();

                    Tangents.SetRowPt(k, Tangent + Node);
                    Tangents.SetRowPt(k + K, Node);
                }
            }

            var TangentsTransformed = TempBuffer.GetTempMultidimensionalarray(out int iTangentsTransformedBuffer, L, 2 * K, D);
            gridData.TransformLocal2Global(Tangents, jCell0, L, TangentsTransformed);

            for(int j = 0; j < L; j++) {
                var TangentsTransformed_j = TangentsTransformed.ExtractSubArrayShallow(j, -1, -1);

                // the metric is the local stretching of the tangent on the boundary line
                for(int k = 0; k < K; k++) {
                    var TangentTransformed = TangentsTransformed_j.GetRowPt(k) - TangentsTransformed_j.GetRow(k + K);

                    metric[j, k] = TangentTransformed.L2Norm();
                }
            }

            TempBuffer.FreeTempBuffer(iTangengBuffer);
            TempBuffer.FreeTempBuffer(iTangentsTransformedBuffer);

            return metric;

        }
    }
}
