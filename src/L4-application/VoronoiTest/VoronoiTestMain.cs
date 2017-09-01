using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.Grid.Voronoi;
using ilPSP;
using BoSSS.Solution.Gnuplot;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Aggregation;

namespace VoronoiTest {
    class VoronoitestMain {
        public static void SetUp() {
            bool MpiInit;
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out MpiInit);
        }

        public static void Cleanup() {
            Console.Out.Dispose();
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }



        static void Main(string[] args) {
            // =======
            // startup
            // =======
            ilPSP.Connectors.Matlab.BatchmodeConnector.Flav = ilPSP.Connectors.Matlab.BatchmodeConnector.Flavor.Octave;
            ilPSP.Connectors.Matlab.BatchmodeConnector.MatlabExecuteable = @"C:\cygwin64\bin\bash.exe";
            SetUp();

            /*
            int J = 100; // number of cells
            int D = 2; //   spatial dimension

            var vG = new VoronoiGrid();
            vG.DelaunayVertices = MultidimensionalArray.Create(J, D);

            // random Delaunay vertices
            Random rnd = new Random(2345);
            for(int j = 0; j < J; j++) {
                for(int d = 0; d < D; d++) {
                    vG.DelaunayVertices[j, d] = rnd.NextDouble();
                }
            }

            // run matlab
            vG.CreateWithMatlab();
            */

            var fineGrdIO = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, 1, 20), GenericBlas.Linspace(-1, 1, 20));
            var fineGrd = new GridData(fineGrdIO);

            var coarseGrd = CoarseningAlgorithms.Coarsen(fineGrd);


            // ========
            // teardown
            // ========
            Cleanup();

        }
    }
}
