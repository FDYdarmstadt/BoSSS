using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding {
    public class GridImportTest {

        public static void Main() {
            Debug.Assert(nFaces == faces.Length);
            Debug.Assert(nInternalFaces == neighbour.Length);
            Debug.Assert(nFaces == owner.Length);
            Debug.Assert(nPoints == points.Length);

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

            //
            //
           
            for(int iCell = 0; iCell < nCells; iCell++) {
                var Cell2Faces = Cells2Faces[iCell];

                if(    Cell2Faces.Count == 6
                   && !Cell2Faces.Select(iFace => faces[iFace].Length == 4).Contains(false)) {
                    // +++++++++++++++
                    // found some Cube
                    // +++++++++++++++

                    int[] face1 = faces[Cell2Faces[0]];
                    int[] face2 = faces[Cell2Faces[1]];
                    int[] face3 = faces[Cell2Faces[2]];
                    int[] face4 = faces[Cell2Faces[3]];
                    int[] face5 = faces[Cell2Faces[4]];
                    int[] face6 = faces[Cell2Faces[5]];

                    HashSet<int> NodesUnsorted = new HashSet<int>();
                    foreach(int iFace in Cell2Faces) {
                        foreach(int iPoint in faces[iFace]) {
                            NodesUnsorted.Add(iPoint);
                        }
                    }
                    if(NodesUnsorted.Count != 8) {
                        throw new NotImplementedException("Degenerate cubes are not supported.");
                    }

                    Cell cl = new Cell();
                    cl.GlobalID = iCell;
                    cl.Type = CellType.Cube_8;
                    cl.NodeIndices = new int[8];

                    // Bottom
                    cl.NodeIndices[0] = face1[0];
                    cl.NodeIndices[1] = face1[1];
                    cl.NodeIndices[6] = face1[2];
                    cl.NodeIndices[3] = face1[3];




                    
                } else {
                    throw new NotImplementedException("At the moment, only support for hex cells.");
                }

            }
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
