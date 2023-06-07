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
    static public class GridImportTest {

        
        internal static int[][] faces = new int[][] {
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

        internal static int[] neighbour = new int[] {
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

        internal static int[] owner = new int[] {
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

        internal static double[,] points = new double[,] {
             {0,0,0},
             {0.03333333333,0,0},
             {0.06666666667,0,0},
             {0.1,0,0},
             {0,0.03333333333,0},
             {0.03333333333,0.03333333333,0},
             {0.06666666667,0.03333333333,0},
             {0.1,0.03333333333,0},
             {0,0.06666666667,0},
             {0.03333333333,0.06666666667,0},
             {0.06666666667,0.06666666667,0},
             {0.1,0.06666666667,0},
             {0,0.1,0},
             {0.03333333333,0.1,0},
             {0.06666666667,0.1,0},
             {0.1,0.1,0},
             {0,0,0.01},
             {0.03333333333,0,0.01},
             {0.06666666667,0,0.01},
             {0.1,0,0.01},
             {0,0.03333333333,0.01},
             {0.03333333333,0.03333333333,0.01},
             {0.06666666667,0.03333333333,0.01},
             {0.1,0.03333333333,0.01},
             {0,0.06666666667,0.01},
             {0.03333333333,0.06666666667,0.01},
             {0.06666666667,0.06666666667,0.01},
             {0.1,0.06666666667,0.01},
             {0,0.1,0.01},
             {0.03333333333,0.1,0.01},
             {0.06666666667,0.1,0.01},
             {0.1,0.1,0.01}
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

        static int Max(int a, int b) {
            if (a > b)
            {
                return a;
            } else {
                return b;
            }
        }
        /// <summary>
        /// test for <see cref="OpenFOAMGrid"/>
        /// </summary>
        [Test]
        public static void ConvertFOAMGrid() {

            int nCells = Max(owner.Max(), neighbour.Max()) + 1;
            var g = GenerateFOAMGrid();

            Assert.AreEqual(g.GridData.iLogicalCells.Count, nCells, "Mismatch in expected number of cells.");
        }

        public static OpenFOAMGrid GenerateFOAMGrid() {

            int nCells = Max(owner.Max(), neighbour.Max()) + 1;
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
