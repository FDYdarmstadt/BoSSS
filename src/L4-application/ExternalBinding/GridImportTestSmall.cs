using BoSSS.Foundation.Grid;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.ExternalBinding {


    /// <summary>
    /// Basic Testing for external language binding.
    /// </summary>
    [TestFixture]
    static public class GridImportTestSmall {

        
        internal static int[][] faces = new int[][] {
            // new int[] {1, 5, 13, 9},
            // new int[] {2, 6, 14, 10},
            // new int[] {0, 8, 12, 4},
            // new int[] {3, 7, 15, 11},
            // new int[] {0, 1, 9, 8},
            // new int[] {1, 2, 10, 9},
            // new int[] {2, 3, 11, 10},
            // new int[] {4, 12, 13, 5},
            // new int[] {5, 13, 14, 6},
            // new int[] {6, 14, 15, 7},
            // new int[] {0, 4, 5, 1},
            // new int[] {1, 5, 6, 2},
            // new int[] {2, 6, 7, 3},
            // new int[] {8, 9, 13, 12},
            // new int[] {9, 10, 14, 13},
            // new int[] {10, 11, 15, 14}
            new int[] {1, 7, 19, 13},
            new int[] {2, 8, 20, 14},
            new int[] {3, 9, 21, 15},
            new int[] {4, 10, 22, 16},
            new int[] {0, 12, 18, 6},
            new int[] {5, 11, 23, 17},
            new int[] {0, 1, 13, 12},
            new int[] {1, 2, 14, 13},
            new int[] {2, 3, 15, 14},
            new int[] {3, 4, 16, 15},
            new int[] {4, 5, 17, 16},
            new int[] {6, 18, 19, 7},
            new int[] {7, 19, 20, 8},
            new int[] {8, 20, 21, 9},
            new int[] {9, 21, 22, 10},
            new int[] {10, 22, 23, 11},
            new int[] {0, 6, 7, 1},
            new int[] {1, 7, 8, 2},
            new int[] {2, 8, 9, 3},
            new int[] {3, 9, 10, 4},
            new int[] {4, 10, 11, 5},
            new int[] {12, 13, 19, 18},
            new int[] {13, 14, 20, 19},
            new int[] {14, 15, 21, 20},
            new int[] {15, 16, 22, 21},
            new int[] {16, 17, 23, 22}
        };

        internal static int[] neighbour = new int[] {
            // 1, 2
            1,2,3,4
        };

        internal static int[] owner = new int[] {
            // 0,
            // 1,
            // 0,
            // 2,
            // 0,
            // 1,
            // 2,
            // 0,
            // 1,
            // 2,
            // 0,
            // 1,
            // 2,
            // 0,
            // 1,
            // 20
            0,
1,
2,
3,
0,
4,
0,
1,
2,
3,
4,
0,
1,
2,
3,
4,
0,
1,
2,
3,
4,
0,
1,
2,
3,
4
        };

        internal static double[,] points = new double[,] {
            // {0, 0, 0},
            // {1.666666667, 0, 0},
            // {3.333333333, 0, 0},
            // {5, 0, 0},
            // {0, 5, 0},
            // {1.666666667, 5, 0},
            // {3.333333333, 5, 0},
            // {5, 5, 0},
            // {0, 0, 1},
            // {1.666666667, 0, 1},
            // {3.333333333, 0, 1},
            // {5, 0, 1},
            // {0, 5, 1},
            // {1.666666667, 5, 1},
            // {3.333333333, 5, 1},
            // {5, 5, 1}
{0, 0, 0},
{1, 0, 0},
{2, 0, 0},
{3, 0, 0},
{4, 0, 0},
{5, 0, 0},
{0, 5, 0},
{1, 5, 0},
{2, 5, 0},
{3, 5, 0},
{4, 5, 0},
{5, 5, 0},
{0, 0, 1},
{1, 0, 1},
{2, 0, 1},
{3, 0, 1},
{4, 0, 1},
{5, 0, 1},
{0, 5, 1},
{1, 5, 1},
{2, 5, 1},
{3, 5, 1},
{4, 5, 1},
{5, 5, 1}
            };

        internal static string[] names = new string[] {"left", "right", "empty"};
        internal static int[] patchIDs = new int[] {0,
                                                       1,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2
        };
        

        /*
        internal static int[][] faces = new int[][] {
            new int[] {1, 6 ,16 ,11 },
            new int[] {2, 7, 17, 12 },
            new int[] {3, 8, 18, 13 },
            new int[] {0, 10, 15, 5 },
            new int[] {4, 9, 19, 14 },
            new int[] {0, 1, 11, 10 },
            new int[] {1, 2, 12, 11 },
            new int[] {2, 3, 13, 12 },
            new int[] {3, 4, 14, 13 },
            new int[] {5, 15, 16, 6 },
            new int[] {6, 16, 17, 7 },
            new int[] {7, 17, 18, 8 },
            new int[] {8, 18, 19, 9 },
            new int[] {0, 5, 6, 1 },
            new int[] {1, 6, 7, 2 },
            new int[] {2, 7, 8, 3 },
            new int[] {3, 8, 9, 4 },
            new int[] {10, 11, 16, 15 },
            new int[] {11, 12, 17, 16 },
            new int[] {12, 13, 18, 17 },
            new int[] {13, 14, 19, 18 },
        };

        internal static double[,] points = new double[,] {
            { 0, 0, 0 },
            { 1.25, 0, 0  },
            { 2.5, 0, 0  },
            { 3.75, 0, 0  },
            { 5, 0, 0  },
            { 0, 1, 0  },
            { 1.25, 1, 0  },
            { 2.5, 1, 0  },
            { 3.75, 1, 0  },
            { 5, 1, 0  },
            { 0, 0, 1  },
            { 1.25, 0, 1  },
            { 2.5, 0, 1  },
            { 3.75, 0, 1  },
            { 5, 0, 1  },
            { 0, 1, 1  },
            { 1.25, 1, 1  },
            { 2.5, 1, 1  },
            { 3.75, 1 ,1  },
            { 5, 1, 1  }
        };

        internal static int[] owner = new int[] {
            0,
            1,
            2,
            0,
            3,
            0,
            1,
            2,
            3,
            0,
            1,
            2,
            3,
            0,
            1,
            2,
            3,
            0,
            1,
            2,
            3 
        };


        internal static int[] neighbour = new int[] {
            1,
            2,
            3
        };

        internal static string[] names = new string[] { "left", "right", "empty" };
        internal static int[] patchIDs = new int[] { 0,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2 
        };
        */

        public static void GridImportMain() {
            
            //int nPoints = points.GetLength(0);
            //int nFaces = owner.Length;
            //int nInternalFaces = neighbour.Length;


            Init();
            ConvertFOAMGrid();
            Cleanup();

        }

        /// <summary>
        /// test for <see cref="OpenFOAMGrid"/>
        /// </summary>
        [Test]
        public static void ConvertFOAMGrid() {

            int nCells = Math.Max(owner.Max(), neighbour.Max()) + 1;
            var g = GenerateFOAMGrid();

            Assert.AreEqual(g.GridData.iLogicalCells.Count, nCells, "Mismatch in expected number of cells.");
        }

        public static OpenFOAMGrid GenerateFOAMGrid() {

            int nCells = Math.Max(owner.Max(), neighbour.Max()) + 1;
            var g = new OpenFOAMGrid(nCells, faces, neighbour, owner, points, names, patchIDs, -1);

            return g;
        }

        
        static Initializer MyInit;

        /// <summary>
        /// MPI Init
        /// </summary>
        [OneTimeSetUp]
        public static void Init() {
            MyInit = new Initializer();
            Console.WriteLine("Test");
            MyInit.BoSSSInitialize();
        }

        /// <summary>
        /// MPI shutdown
        /// </summary>
        [OneTimeTearDown]
        public static void Cleanup() {
            MyInit.BoSSSFinalize();
        }

    }
}
