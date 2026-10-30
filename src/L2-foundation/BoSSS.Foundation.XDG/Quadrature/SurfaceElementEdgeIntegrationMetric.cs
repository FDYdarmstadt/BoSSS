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
            var TangentLengthInCell1 = new double[K];

            for(int i = 0; i < L; i++) {
                int jCell0 = E2C[iEdge + i, 0];
                int iFace = E2F[iEdge + i, 0];
                var Kref = gridData.iGeomCells.GetRefElement(jCell0);

                // On a hanging face, cell 0 may be the coarse cell, which is not conformal with the edge; the rule is then taken from cell 1
                // (see EdgeRuleFromCellBoundaryFactory), i.e. its weights are line measures in the reference coordinates of cell 1.
                // The tangent below is a unit vector in the reference coordinates of cell 0, which is longer in those of cell 1
                // (by the refinement ratio, i.e. twice as long for 2:1); the metric is divided by this length.
                // Only the geometry of cell 1 is used, since it may be an external cell, where the level set is not available.
                bool scaleToCell1 = !gridData.iGeomEdges.IsEdgeConformalWithCell1(iEdge + i);
                MultidimensionalArray M0 = null, M1 = null;
                if(scaleToCell1) {
                    M0 = gridData.iGeomEdges.Edge2CellTrafos[E2Ctrafo[iEdge + i, 0]].Matrix;
                    M1 = gridData.iGeomEdges.Edge2CellTrafos[E2Ctrafo[iEdge + i, 1]].Matrix;
                }

                var NodesCell = qr.Nodes.GetVolumeNodeSet(gridData, E2Ctrafo[iEdge + i, 0], false);
                var LevSetNormals = this.m_levelSetData.GetLevelSetReferenceNormals(NodesCell, jCell0, 1).ExtractSubArrayShallow(0, -1, -1);
                var CellNormal = Kref.FaceNormals.GetRowPt(iFace);



                for(int k = 0; k < K; k++) {
                    TangentLengthInCell1[k] = 1.0;
                    if(qr.Weights[k] != 0) { // avoid empty rules
                        var SurfNormal = LevSetNormals.GetRowPt(k);
                        var Node = NodesCell.GetRowPt(k);

                        var Tangent = CellNormal.CrossProduct(SurfNormal);
                        Tangent.NormalizeInPlace();

                        Tangents.SetRowPt(k, Tangent + Node);
                        Tangents.SetRowPt(k + K, Node);

                        if(scaleToCell1)
                            TangentLengthInCell1[k] = LengthInOtherCell(M0, M1, Tangent);
                    }
                }

                var TangentsTransformed = gridData.TransformLocal2Global(Tangents, jCell0);

                // the metric is the local stretching of the tangent on the boundary line
                for(int k = 0; k < K; k++) {
                    if(qr.Weights[k] != 0) { // avoid empty rules
                        var TangentTransformed = TangentsTransformed.GetRowPt(k) - TangentsTransformed.GetRow(k + K);
                        metric[i, k] = TangentTransformed.L2Norm() / TangentLengthInCell1[k];
                    }
                }
            }

            TempBuffer.FreeTempBuffer(iTangentsBuf);

            return metric;
        }

        /// <summary>
        /// Length of a vector <paramref name="t"/> on an edge, given in the reference coordinates of cell 0, in the reference coordinates of cell 1:
        /// with the matrices <paramref name="M0"/>, <paramref name="M1"/> of the edge-to-cell transformations, the vector in edge coordinates is
        /// $ \tau = (M_0^T M_0)^{-1} M_0^T t $ and the result is $ | M_1 \tau | $.
        /// </summary>
        static double LengthInOtherCell(MultidimensionalArray M0, MultidimensionalArray M1, Vector t) {
            int D = M0.GetLength(0);
            if(M0.GetLength(1) != 2)
                throw new NotSupportedException("only implemented for edges of 3D cells");

            double g00 = 0, g01 = 0, g11 = 0, r0 = 0, r1 = 0;
            for(int d = 0; d < D; d++) {
                g00 += M0[d, 0] * M0[d, 0];
                g01 += M0[d, 0] * M0[d, 1];
                g11 += M0[d, 1] * M0[d, 1];
                r0 += M0[d, 0] * t[d];
                r1 += M0[d, 1] * t[d];
            }
            double det = g00 * g11 - g01 * g01;
            double tau0 = (g11 * r0 - g01 * r1) / det;
            double tau1 = (g00 * r1 - g01 * r0) / det;

            double acc = 0;
            for(int d = 0; d < D; d++) {
                double v = M1[d, 0] * tau0 + M1[d, 1] * tau1;
                acc += v * v;
            }
            return Math.Sqrt(acc);
        }
    }




}


