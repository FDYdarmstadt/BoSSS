using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature {

    /// <summary>
    /// Integration metric for $D - 2$-dimensional integral the boundary of surface elements for each edge.
    /// This is a point in 2D and a line in 3D,
    /// i.e., for an edge $e_i$ which is subset of some cells $ K_j $ boundary, i.e., $ e_i \subset K_j $  an integral of the type
    /// \[
    ///    \int_{e_j \cap \mathfrak{I}}  f \mathrm{dl} .
    /// \]
    /// </summary>
    public class SurfaceElementEdgeIntegrationMetric : IIntegrationMetric {


        public SurfaceElementEdgeIntegrationMetric(LevelSetTracker.LevelSetData levelSetData) {
            m_levelSetData = levelSetData;
        }

        readonly LevelSetTracker.LevelSetData m_levelSetData;

        /// <summary>
        /// For the integral over the surface element boundary, 
        /// the metric depends on the node (if the element transformation applies some shearing transformation),
        /// but only in 3D
        /// </summary>
        public bool AlwaysUsePerNodeScaling {
            get {
                if(m_levelSetData.GridDat.SpatialDimension == 2)
                    return false;
                if(m_levelSetData.GridDat.SpatialDimension == 3)
                    return true;
                throw new NotSupportedException();
            }
        }


        /// <summary>
        /// in 2D, the boundary of surface elements are just points -- one needs a count measure, where the integration metric is always 1
        /// </summary>
        public MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int iEdge, int L) {
            if(m_levelSetData.GridDat.SpatialDimension == 2) {
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
        public MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int iEdge, int L) {
            if(!object.ReferenceEquals(gridData, m_levelSetData.GridDat))
                throw new ArgumentException();

            //throw new NotImplementedException("todo");

            int[,] E2C = gridData.iGeomEdges.CellIndices;
            int[,] E2Ctrafo = gridData.iGeomEdges.Edge2CellTrafoIndex;
            byte[,] E2F = gridData.iGeomEdges.FaceIndices;
            int D = gridData.SpatialDimension;
            int K = qr.NoOfNodes;

            var metric = MultidimensionalArray.Create(L, K);
            var Tangents = TempBuffer.GetTempMultidimensionalarray(out int iTangentsBuf, 2 * K, D);

            for(int i = 0; i < L; i++) {
                // The metric must be computed in a cell which is conformal with the edge:
                // the weights of the rule are line measures in the reference coordinates of the edge, resp. of a cell conformal with it
                // (EdgeRuleFromCellBoundaryFactory takes the rule from such a cell). On a hanging face, cell 0 may be the coarse cell, which is
                // not conformal; its tangent stretching is too large by the refinement ratio (factor 2 on the line integrals for 2:1).
                int iSide = gridData.iGeomEdges.IsEdgeConformalWithCell1(iEdge + i) ? 0 : 1;
                int jCell = E2C[iEdge + i, iSide];
                int iFace = E2F[iEdge + i, iSide];
                var Kref = gridData.iGeomCells.GetRefElement(jCell);

                var NodesCell = qr.Nodes.GetVolumeNodeSet(gridData, E2Ctrafo[iEdge + i, iSide], false);
                var LevSetNormals = this.m_levelSetData.GetLevelSetReferenceNormals(NodesCell, jCell, 1).ExtractSubArrayShallow(0, -1, -1);
                var CellNormal = Kref.FaceNormals.GetRowPt(iFace);



                for(int k = 0; k < K; k++) {
                    if(qr.Weights[k] != 0) { // avoid empty rules
                        var SurfNormal = LevSetNormals.GetRowPt(k);
                        var Node = NodesCell.GetRowPt(k);

                        var Tangent = CellNormal.CrossProduct(SurfNormal);
                        Tangent.NormalizeInPlace();

                        Tangents.SetRowPt(k, Tangent + Node);
                        Tangents.SetRowPt(k + K, Node);
                    }
                }

                var TangentsTransformed = gridData.TransformLocal2Global(Tangents, jCell);

                // the metric is the local stretching of the tangent on the boundary line
                for(int k = 0; k < K; k++) {
                    if(qr.Weights[k] != 0) { // avoid empty rules
                        var TangentTransformed = TangentsTransformed.GetRowPt(k) - TangentsTransformed.GetRow(k + K);
                        metric[i, k] = TangentTransformed.L2Norm();
                    }
                }
            }

            TempBuffer.FreeTempBuffer(iTangentsBuf);

            return metric;
        }
    }




}


