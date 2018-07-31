using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding {
    public class GridImportTest {

        public static void Main() {

            // Checks
            // ======

            Common_.BoSSSInitialize();
            Debug.Assert(nFaces == faces.Length);
            Debug.Assert(nInternalFaces == neighbour.Length);
            Debug.Assert(nFaces == owner.Length);
            Debug.Assert(nPoints == points.Length);


            // Build Cells-to-Faces correlation
            // ================================

            List<int>[] Cells2Faces = new List<int>[nCells];
            for(int j = 0; j < nCells; j++) {
                Cells2Faces[j] = new List<int>();
            }

            for(int iFace = 0; iFace < nFaces; iFace++) {
                Debug.Assert(Cells2Faces[owner[iFace]].Contains(iFace) == false);
                Cells2Faces[owner[iFace]].Add(iFace);

            }

            for(int iInternalFace = 0; iInternalFace < nInternalFaces; iInternalFace++) {
                int NeighCell = neighbour[iInternalFace];

                Debug.Assert(Cells2Faces[NeighCell].Contains(iInternalFace) == false);
                Cells2Faces[NeighCell].Add(iInternalFace);
            }

            // convert to BoSSS grid
            // =====================

            Cell[] bosss_cells = new Cell[nCells];

            for(int iCell = 0; iCell < nCells; iCell++) {
                var Cell2Faces = Cells2Faces[iCell];

                if(    Cell2Faces.Count == 6
                   && !Cell2Faces.Select(iFace => faces[iFace].Length == 4).Contains(false)) {
                    // +++++++++++++++
                    // found some Cube
                    // +++++++++++++++

                    int[][] cube_faces = Cell2Faces.Select(iFace => faces[iFace]).ToArray();
                    Debug.Assert(cube_faces.Length == 6);

                    HashSet<int> NodesUnsorted = new HashSet<int>();
                    foreach(int iFace in Cell2Faces) {
                        foreach(int iPoint in faces[iFace]) {
                            NodesUnsorted.Add(iPoint);
                        }
                    }
                    if(NodesUnsorted.Count != 8) {
                        throw new NotImplementedException("Degenerate cubes are not supported.");
                    }

                    // create cell
                    Cell cl = new Cell();
                    bosss_cells[iCell] = cl;
                    cl.GlobalID = iCell;
                    cl.Type = CellType.Cube_8;
                    cl.NodeIndices = new int[8];
                    cl.TransformationParams = MultidimensionalArray.Create(8, 3);

                    // Find point indices for bottom, i.e. 'y == -1' (in ref. coordinates):
                    cl.NodeIndices[0] = cube_faces[0][0];
                    cl.NodeIndices[1] = cube_faces[0][1];
                    cl.NodeIndices[6] = cube_faces[0][2];
                    cl.NodeIndices[3] = cube_faces[0][3];


                    // Find point indices for 'x == -1'-wall (in ref. coordinates):
                    bool found = false;
                    for(int iCF = 1; iCF < 6; iCF++) {
                        found = ContainsEdge(cube_faces[iCF], cl.NodeIndices[0], cl.NodeIndices[3], out cl.NodeIndices[2], out cl.NodeIndices[7]);
                        if (found)
                            break;
                    }
                    if (found == false)
                        throw new NotSupportedException("weired cell topology in cell " + iCell + "; assuming a cube/hex;");

                    // Find point indices for 'x == +1'-wall (in ref. coordinates):
                    found = false;
                    for(int iCF = 1; iCF < 6; iCF++) {
                        found = ContainsEdge(cube_faces[iCF], cl.NodeIndices[1], cl.NodeIndices[6], out cl.NodeIndices[4], out cl.NodeIndices[5]);
                        if (found)
                            break;
                    }
                    if (found == false)
                        throw new NotSupportedException("weired cell topology in cell " + iCell + "; assuming a cube/hex;");


                    // set point coordinates
                    for(int i = 0; i < 8; i++) {
                        cl.TransformationParams.SetRow(i, points[cl.NodeIndices[i]]);
                    }

                    // Jacobian Test
                    JacobianTest(cl, out bool PositiveJacobianFlag, out bool NegativeJacobianFlag, out bool Linear);
                    if(NegativeJacobianFlag) {
                        // flip top/bottom
                        int _0 = cl.NodeIndices[0];
                        int _1 = cl.NodeIndices[1];
                        int _6 = cl.NodeIndices[6];
                        int _3 = cl.NodeIndices[3];

                        cl.NodeIndices[0] = cl.NodeIndices[2];
                        cl.NodeIndices[1] = cl.NodeIndices[4];
                        cl.NodeIndices[6] = cl.NodeIndices[5];
                        cl.NodeIndices[3] = cl.NodeIndices[7];

                        cl.NodeIndices[2] = _0;
                        cl.NodeIndices[4] = _1;
                        cl.NodeIndices[5] = _6;
                        cl.NodeIndices[7] = _3;
                        
                        for (int i = 0; i < 8; i++) {
                            cl.TransformationParams.SetRow(i, points[cl.NodeIndices[i]]);
                        }

                        JacobianTest(cl, out PositiveJacobianFlag, out NegativeJacobianFlag, out Linear);
                    }

                    if(NegativeJacobianFlag || !PositiveJacobianFlag) {
                        throw new NotSupportedException("Found degenerate Jacobian in cell " + iCell);
                    }
                    Debug.Assert(PositiveJacobianFlag == true);
                    Debug.Assert(NegativeJacobianFlag == false);

                    if (Linear) {
                        cl.Type = CellType.Cube_Linear;
                        cl.TransformationParams = cl.TransformationParams.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { 3, 2 }).CloneAs();
                    }
                    

                } else {
                    throw new NotImplementedException("At the moment, only support for hex cells.");
                }
            }

            // final test
            // ==========

            GridCommons grd = new GridCommons(new RefElement[] { Cube.Instance }, new RefElement[] { Square.Instance }) {
                Cells = bosss_cells,
                Description = "imported form OpenFOAM"
            };

            var gDat = new GridData(grd);



            // Finalize
            // ========
            Common_.BoSSSFinalize();
        }


        static void JacobianTest(Cell cl, out bool PositiveJacobianFlag, out bool NegativeJacobianFlag, out bool Linear) {
            RefElement Kref = cl.Type.GetRefElement();

            int D = 3; // for OpenFOAM, we assume 3D;

            // get nodes within ref element, to test the Jacobian
            int deg = Kref.GetInterpolationDegree(cl.Type);
            if (deg > 1)
                deg--;
            deg *= 3;

            NodeSet TestNodes = Kref.GetQuadratureRule(2 * deg).Nodes;

            // evaluate derivatives of nodal polynomials for transformation
            PolynomialList[] Deriv = Kref.GetInterpolationPolynomials1stDeriv(cl.Type);
            MultidimensionalArray[] DerivEval = new MultidimensionalArray[D];
            for (int d = 0; d < D; d++) {
                DerivEval[d] = Deriv[d].Values.GetValues(TestNodes);
            }

            // evaluate Jacobian matrix
            MultidimensionalArray Jacobi = MultidimensionalArray.Create(TestNodes.NoOfNodes, D, D); // temporary storage for Jacobian matrix

            Debug.Assert(cl.TransformationParams.Dimension == 2);
            Debug.Assert(cl.TransformationParams.GetLength(1) == 3);

            for (int d1 = 0; d1 < D; d1++) {
                Debug.Assert(cl.TransformationParams.GetLength(0) == Deriv[d1].Count);
                MultidimensionalArray JacobiCol = Jacobi.ExtractSubArrayShallow(-1, -1, d1);
                JacobiCol.Multiply(1.0, DerivEval[d1], cl.TransformationParams, 0.0, "kd", "kn", "nd");
            }

            // do the tests
            NegativeJacobianFlag = false;
            PositiveJacobianFlag = false;
            double MinAbsJacobiDet = double.MaxValue;
            double MaxAbsJacobiDet = 0.0;
            for (int n = 0; n < TestNodes.NoOfNodes; n++) {
                double detJac = Jacobi.ExtractSubArrayShallow(n, -1, -1).Determinant();
                if (detJac <= 0) {
                    NegativeJacobianFlag = true;
                }
                if (detJac > 0) {
                    PositiveJacobianFlag = true;
                }

                MinAbsJacobiDet = Math.Min(MinAbsJacobiDet, detJac.Abs());
                MaxAbsJacobiDet = Math.Max(MaxAbsJacobiDet, detJac.Abs());
            }

            if ((MaxAbsJacobiDet - MinAbsJacobiDet) / MinAbsJacobiDet <= 1.0e-8)
                Linear = true;
            else
                Linear = false;

        }


        static bool ContainsEdge(int[] face, int iPt1, int iPt2, out int iPrev, out int iNext) {
            int L = face.Length;
            for (int i = 0; i < L; i++) {
                int v1 = face[i];
                int v2 = face[(i + 1) % L];
                int _iPrev = face[i > 0 ? (i - 1) : (L - 1)];
                int _iNext = face[(i + 2) % L];

                if (v1 == iPt1 && v2 == iPt2) {
                    iPrev = _iPrev;
                    iNext = _iNext;
                    return true;
                }
                if (v2 == iPt1 && v1 == iPt2) {
                    iNext = _iPrev;
                    iPrev= _iNext;
                    return true;
                }
            }
            iNext = int.MinValue;
            iPrev = int.MinValue;
            return false;
        }


        // ---------------------------------
        // test data from OpenFOAM tutorials
        //

        static int nPoints = 32;
        static int nCells = 9;
        static int nFaces = 42;
        static int nInternalFaces = 12;


        static int[][] faces = new int[][] {
            new int[] {1, 5, 21, 17},
            new int[] {4, 20, 21, 5},
            new int[] {2, 6, 22, 18},
            new int[] {5, 21, 22, 6},
            new int[] {6, 22, 23, 7},
            new int[] {5, 9, 25, 21},
            new int[] {8, 24, 25, 9},
            new int[] {6, 10, 26, 22},
            new int[] {9, 25, 26, 10},
            new int[] {10, 26, 27, 11},
            new int[] {9, 13, 29, 25},
            new int[] {10, 14, 30, 26},
            new int[] {12, 28, 29, 13},
            new int[] {13, 29, 30, 14},
            new int[] {14, 30, 31, 15},
            new int[] {0, 16, 20, 4},
            new int[] {4, 20, 24, 8},
            new int[] {8, 24, 28, 12},
            new int[] {3, 7, 23, 19},
            new int[] {7, 11, 27, 23},
            new int[] {11, 15, 31, 27},
            new int[] {0, 1, 17, 16},
            new int[] {1, 2, 18, 17},
            new int[] {2, 3, 19, 18},
            new int[] {0, 4, 5, 1},
            new int[] {4, 8, 9, 5},
            new int[] {8, 12, 13, 9},
            new int[] {1, 5, 6, 2},
            new int[] {5, 9, 10, 6},
            new int[] {9, 13, 14, 10},
            new int[] {2, 6, 7, 3},
            new int[] {6, 10, 11, 7},
            new int[] {10, 14, 15, 11},
            new int[] {16, 17, 21, 20},
            new int[] {20, 21, 25, 24},
            new int[] {24, 25, 29, 28},
            new int[] {17, 18, 22, 21},
            new int[] {21, 22, 26, 25},
            new int[] {25, 26, 30, 29},
            new int[] {18, 19, 23, 22},
            new int[] {22, 23, 27, 26},
            new int[] {26, 27, 31, 30}
        };

        static int[] neighbour = new int[] {
            1,
            3,
            2,
            4,
            5,
            4,
            6,
            5,
            7,
            8,
            7,
            8
        };

        static int[] owner = new int[] {
            0,
            0,
            1,
            1,
            2,
            3,
            3,
            4,
            4,
            5,
            6,
            7,
            6,
            7,
            8,
            0,
            3,
            6,
            2,
            5,
            8,
            0,
            1,
            2,
            0,
            3,
            6,
            1,
            4,
            7,
            2,
            5,
            8,
            0,
            3,
            6,
            1,
            4,
            7,
            2,
            5,
            8
        };

        static double[][] points = new double[][] {
            new double[] {0,0,0},
            new double[] {0.03333333333,0,0},
            new double[] {0.06666666667,0,0},
            new double[] {0.1,0,0},
            new double[] {0,0.03333333333,0},
            new double[] {0.03333333333,0.03333333333,0},
            new double[] {0.06666666667,0.03333333333,0},
            new double[] {0.1,0.03333333333,0},
            new double[] {0,0.06666666667,0},
            new double[] {0.03333333333,0.06666666667,0},
            new double[] {0.06666666667,0.06666666667,0},
            new double[] {0.1,0.06666666667,0},
            new double[] {0,0.1,0},
            new double[] {0.03333333333,0.1,0},
            new double[] {0.06666666667,0.1,0},
            new double[] {0.1,0.1,0},
            new double[] {0,0,0.01},
            new double[] {0.03333333333,0,0.01},
            new double[] {0.06666666667,0,0.01},
            new double[] {0.1,0,0.01},
            new double[] {0,0.03333333333,0.01},
            new double[] {0.03333333333,0.03333333333,0.01},
            new double[] {0.06666666667,0.03333333333,0.01},
            new double[] {0.1,0.03333333333,0.01},
            new double[] {0,0.06666666667,0.01},
            new double[] {0.03333333333,0.06666666667,0.01},
            new double[] {0.06666666667,0.06666666667,0.01},
            new double[] {0.1,0.06666666667,0.01},
            new double[] {0,0.1,0.01},
            new double[] {0.03333333333,0.1,0.01},
            new double[] {0.06666666667,0.1,0.01},
            new double[] {0.1,0.1,0.01}
        };
    }
}
