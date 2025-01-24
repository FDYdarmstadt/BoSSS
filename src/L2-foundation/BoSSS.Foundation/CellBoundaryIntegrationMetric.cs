using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace BoSSS.Foundation.Quadrature {


    /// <summary>
    /// integration metric for area integrals on cell boundaries
    /// </summary>
    public class CellBoundaryIntegrationMetric : IIntegrationMetric {

        /// <summary>
        /// For cell boundaries, the scaling is constant for each face of the cell for linear elements, therefore always false;
        /// c.f. <see cref="IIntegrationMetric.AlwaysUsePerNodeScaling"/>;
        /// </summary>
        public bool AlwaysUsePerNodeScaling => false;


        /// <summary>
        /// integration metrics for linear mapped cells
        /// </summary>
        /// <returns>
        ///  - 1st index: cell 
        ///  - 2nd index: face 
        /// </returns>
        public MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int i0, int L) {

            CellBoundaryQuadRule CurrentRule = qr as CellBoundaryQuadRule;

            int NoOFFaces = CurrentRule.RefElement.NoOfFaces;
            var ret = MultidimensionalArray.Create(L, NoOFFaces);
            var Cell2Edge = gridData.iLogicalCells.Cells2Edges;
            var FaceIndices = gridData.iGeomEdges.FaceIndices;
            var Scalings = gridData.iGeomEdges.SqrtGramian;
            var faceGrm = gridData.iGeomCells.GetRefElement(i0).FaceTrafoGramianSqrt;

            for(int i = 0; i < L; i++) {
                int jCell = i + i0;

                int[] Cell2edge_j = Cell2Edge[jCell];
                Debug.Assert(gridData.iGeomCells.IsCellAffineLinear(jCell) == true);

                for(int iF = 0; iF < NoOFFaces; iF++) {
                    ret[i, iF] = 0;
                }

                for(int e = Cell2edge_j.Length - 1; e >= 0; e--) {
                    int iEdge = Math.Abs(Cell2edge_j[e]) - 1;
                    int iFace = FaceIndices[iEdge, Cell2edge_j[e] > 0 ? 0 : 1];

                    ret[i, iFace] += Scalings[iEdge] / faceGrm[iFace];
                }
            }

#if DEBUG
            MultidimensionalArray NonLinScaling = GetScalingsForNonlinElements(gridData, qr, i0, L);
            int[] NodesPerFace = CurrentRule.NumbersOfNodesPerFace;
            for(int l = 0; l < L; l++) {
                int k = 0;
                for(int iFace = 0; iFace < NoOFFaces; iFace++) {
                    for(int n = 0; n < NodesPerFace[iFace]; n++) {
                        Debug.Assert(Math.Abs(NonLinScaling[l, k] - ret[l, iFace]) <= 1.0e-8);
                        k++;
                    }
                }
            }
#endif

            return ret;
        }

        /// <summary>
        /// integration metrics for non-linear mapped cells
        /// </summary>
        /// <returns>
        ///  - 1st index: cell
        ///  - 2nd index: quadrature node
        /// </returns>
        public MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int i0, int L) {

            CellBoundaryQuadRule CurrentRule = qr as CellBoundaryQuadRule;


            // Nodes in the edge coordinate system
            int[] NodesPerFace = CurrentRule.NumbersOfNodesPerFace;

            // Cell-Reference element:
            // we can assume that cells {i0, ... , i0 + Len - 1}
            // map to the same reference element
            var RefElmCell = gridData.iGeomCells.GetRefElement(i0);

            // bla
            int N = qr.Nodes.NoOfNodes;
            int D = gridData.SpatialDimension;
            var Trafos = gridData.iGeomEdges.Edge2CellTrafos;
            int NoOfFaces = RefElmCell.NoOfFaces;
            Debug.Assert(NoOfFaces == NodesPerFace.Length);

            // normals in reference space
            MultidimensionalArray RefNormals = RefElmCell.FaceNormals;

            // return value
            MultidimensionalArray _R = MultidimensionalArray.Create(L, N);



            NodeSet[] NodesForFace = new NodeSet[NoOfFaces];
            {
                int n = 0;
                for(int iFace = 0; iFace < NoOfFaces; iFace++) {
                    int Ne = NodesPerFace[iFace];
                    if(Ne > 0)
                        NodesForFace[iFace] = new NodeSet(RefElmCell, qr.Nodes.ExtractSubArrayShallow(new int[] { n, 0 }, new int[] { n + Ne - 1, D - 1 }), true);
                    n += Ne;
                }
                Debug.Assert(n == N);
            }

            MultidimensionalArray Jac = gridData.Jacobian.GetValue_Cell(qr.Nodes, i0, L);
            MultidimensionalArray DetJ = gridData.JacobianDeterminat.GetValue_Cell(qr.Nodes, i0, L);

            // loop over cells:
            for(int j = 0; j < L; j++) {
                int jCell = j + i0;  // cell index


                // return value for Edge # jEdge
                var R = _R.ExtractSubArrayShallow(new int[] { j, 0 }, new int[] { j - 1, N - 1 });

                //// transform nodes to cell coordinates
                //int TrafoIdx = grd.Edges.TrafoIndex[jEdge, 0];
                //var Trafo = grd.Edges.EdgeToCellTrafos[TrafoIdx];
                //Trafo.Transform(NodesEdge, Nodes);
                ////RefElmCell.


                // the Jacobian and it's determinat
                //RefElmCell.JacobianDetTransformation(
                //    Nodes, DetJ, 0, grd.Grid.Cells[jCell].Type, grd.Grid.Cells[jCell].TransformationParams);
                //RefElmCell.JacobianOfTransformation(
                //    Nodes, Jac, 0, grd.Grid.Cells[jCell].Type, grd.Grid.Cells[jCell].TransformationParams);


                int n0 = 0;
                for(int RefFace = 0; RefFace < NoOfFaces; RefFace++) {
                    int Ne = NodesPerFace[RefFace];
                    if(Ne <= 0)
                        continue;

                    // compute normals in physical space
                    MultidimensionalArray CurvedNormals = MultidimensionalArray.Create(Ne, D);
                    gridData.iGeomEdges.GetNormalsForCell(NodesForFace[RefFace], jCell, RefFace, CurvedNormals);


                    for(int _n = 0; _n < Ne; _n++) { // loop over quadrature nodes
                        int n = n0 + _n;
                        double r = 0.0;
                        double sc = 0.0;
                        for(int d = 0; d < D; d++) {//calculate the scalar product between the normals in ref. and phyic. domain (pro quad. node)
                            double acc = 0.0;
                            for(int s = 0; s < D; s++) {
                                acc += RefNormals[RefFace, s] * Jac[j, n, d, s];
                            }
                            sc += acc * CurvedNormals[_n, d];
                        }

                        r += DetJ[j, n] * (1 / (sc));

                        R[n] = r;
                    }

                    n0 += Ne;
                }

            }

            return _R;
        }


    }
}
