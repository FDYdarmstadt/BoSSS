using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP.Tracing;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;

namespace BoSSS.Solution.Statistic {

    /// <summary>
    /// extension methods for <see cref="BoSSS.Foundation.Field"/>
    /// </summary>
    public static class Field_Extension {

        /// <summary>
        /// the broken H1-seminorm of DG field <paramref name="f"/>:
        /// <latex mode="display">
        ///  \left| f \right|_{H^1} 
        ///  :=
        ///  \left\|  \nabla_h f \right\|_{(L^2(\Omega))^D}
        /// </latex>
        /// </summary>
        /// <param name="f"></param>
        /// <param name="CM">optional restriction of domain.</param>
        /// <returns></returns>
        public static double H1SemiNorm(this Field f, CellMask CM = null) {
            using (new FuncTrace()) {

                Field tmp = f.CloneAs();
                int D = f.Basis.Context.GridDat.SpatialDimension;

                double Norm = 0;
                for (int d = 0; d < D; d++) {
                    tmp.Clear();
                    tmp.Derivative(1.0, f, d, CM);

                    double nd = tmp.L2Norm(CM);
                    Norm += nd * nd;
                }

                return Math.Sqrt(Norm);
            }
        }

        /// <summary>
        /// the jump-seminorm of DG-field <paramref name="f"/>:
        /// <latex mode="display">
        ///   
        /// </latex>
        /// </summary>
        /// <param name="f"></param>
        /// <param name="CM"></param>
        /// <returns></returns>
        public static double JumpSemiNorm(this Field f, EdgeMask CM = null) {
            using (new FuncTrace()) {

                var qi = new EdgeQuadratureScheme(CM);

                var quad = new JumpSemiNormIntegrator(f, qi);
                quad.Execute();
                return Math.Sqrt(quad.overallResult);
            }
        }

        class JumpSemiNormIntegrator : EdgeQuadrature {

            public JumpSemiNormIntegrator(Field _f, EdgeQuadratureScheme eq) :
                base(new int[] { 1 }, _f.Basis.Context, eq, _f.Basis.Degree * 2) {
                f = _f;
            }

            Field f;

            protected override NodeSetController.NodeSetContainer[] CreateNodeSetFamily(Platform.MultidimensionalArray NodesUntransformed) {
                Simplex volSimplx = m_Context.Grid.GridSimplex;
                int NoOfEdgesPerSimplex = volSimplx.NoOfEdges;
                int NoOfTrafos = m_Context.GridDat.InterCellTransformations.Length;
                int NoOfNodes = NodesUntransformed.GetLength(0);
                int D = m_Context.Grid.SpatialDimension;

                //// ======
                //// create
                //// ======

                //for (int i = 0; i < m_TestFunctions_1stCell.GetLength(0); i++) {
                //    m_NodesUntransfomedVolCoord_1stEdge[i] = MultidimensionalArray.CreateEmpty(2);
                //    for (int _k = 0; _k < m_NodesUntransfomedVolCoord_2ndEdge.GetLength(1); _k++)
                //        m_NodesUntransfomedVolCoord_2ndEdge[i, _k] = MultidimensionalArray.CreateEmpty(2);
                //}


                //// ==========================
                //// allocate storage structure
                //// ==========================

                //for (int e = 0; e < NoOfEdgesPerSimplex; e++) {
                //    m_NodesUntransfomedVolCoord_1stEdge[e].Allocate(NoOfNodes, D);

                //    for (int k = 0; k < NoOfTrafos; k++) {
                //        m_NodesUntransfomedVolCoord_2ndEdge[e, k].Allocate(NoOfNodes, D);
                //    }
                //}

                // ================================================================================
                // transform quadrature nodes from edge simplex to volume simplex coordinate system
                // ================================================================================

                NodeSetController.NodeSetContainer[] nodeSetFamily = new NodeSetController.NodeSetContainer[NoOfEdgesPerSimplex * (NoOfTrafos + 1)];

                int NoOftrafos = m_Context.GridDat.InterCellTransformations.Length;
                for (int e = 0; e < NoOfEdgesPerSimplex; e++) {
                    var nodes = MultidimensionalArray.Create(NoOfNodes, D);
                    volSimplx.EdgeToVolumeCoordinates(e, NodesUntransformed, nodes);
                    nodeSetFamily[e] = m_Context.NSC.CreateContainer(nodes, 0.0);

                    for (int _k = 0; _k < NoOftrafos; _k++) {
                        var nodes2 = MultidimensionalArray.Create(NoOfNodes, D);
                        m_Context.GridDat.InterCellTransformations[_k].Transform(nodes, nodes2);
                        nodeSetFamily[NodeSetIndex(e, _k)] = m_Context.NSC.CreateContainer(nodes2, 0.0);
                    }
                }
                return nodeSetFamily;
            }

            private int NodeSetIndex(int SecEdge, int ict) {
                GridData grd = m_Context.GridDat;
                Simplex volSimplx = m_Context.Grid.GridSimplex;
                int NoOfIctS = grd.InterCellTransformations.Length;
                int NoOfEdgesPerSimplex = volSimplx.NoOfEdges;
                return (NoOfEdgesPerSimplex + NoOfIctS * SecEdge + ict);
            }
            int GetNodeSetIndex2ndCell(int EdgNo1stCell, int Trafo) {
                int r = m_Context.Grid.GridSimplex.NoOfEdges;

                int NoOfTrafos = m_Context.GridDat.InterCellTransformations.Length;
                r += EdgNo1stCell * NoOfTrafos;
                r += Trafo;
                return r;
            }

            protected override void Evaluate(int i0, int Length, int NoOfNodes, MultidimensionalArray EvalResult) {
                var Edges = m_Context.GridDat.Edges;
                var EdgeIndices = m_Context.GridDat.EdgeIndices;
                byte[] Trafo22ndCell = m_Context.GridDat.TrafoTo2ndCell;

                var res = EvalResult.ExtractSubArrayShallow(-1, -1, 0);
                int I1 = res.GetLength(0), I2 = res.GetLength(1);

                for (int i = 0; i < Length; i++) {
                    int iEdge = i + i0;

                    int jCell1 = Edges[iEdge, 0];
                    int jCell2 = Edges[iEdge, 1];
                    var edge1 = EdgeIndices[iEdge, 0];

                    //m_Context.GridDat.EdgeIndices

                    if (jCell2 >= 0) {
                        f.Evaluate(jCell1, 1, edge1, res, i, 0.0);
                        f.Evaluate(jCell2, 1, GetNodeSetIndex2ndCell(edge1, Trafo22ndCell[iEdge]), res, i, -1.0);
                    } else {
                        for (int i2 = 0; i2 < I2; i2++) {
                            res[i, i2] = 0;
                        }
                    }
                }

                for (int i1 = 0; i1 < I1; i1++) {
                    for (int i2 = 0; i2 < I2; i2++) {
                        double x = res[i1, i2];
                        res[i1, i2] = x * x;
                    }
                }

            }

            public double overallResult = 0;

            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                var hF = m_Context.GridDat.h_max_Edge;

                for (int i = 0; i < Length; i++) {
                    overallResult += (1.0 / hF[i + i0]) * ResultsOfIntegration[i, 0];
                }
            }
        }


    }
}
