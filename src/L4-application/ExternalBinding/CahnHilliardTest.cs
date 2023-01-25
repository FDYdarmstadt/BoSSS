using System;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using NUnit.Framework;

namespace BoSSS.Application.ExternalBinding {

    [TestFixture]
    static public class CahnHilliardTest {

        [Test]
        public static void Main() {

            Init();
            GridImportTest.ConvertFOAMGrid();
            Console.WriteLine("Running Cahn-Hilliard Test");
            var chOp = new FixedOperators();
            // OpenFOAMGrid grd = GridImportTestSmall.GenerateFOAMGrid();
            OpenFOAMGrid smallGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/small/polyMesh/");
            OpenFOAMGrid mediumGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/medium/polyMesh/");
            OpenFOAMGrid largeGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/large/polyMesh/");
            var norms = new List<double>();
            var preNorms = new List<double>();
            var normRelChanges = new List<double>();
            var jumpNorms = new List<double>();
            // OpenFOAMGrid grd = GridImportTest.GenerateFOAMGrid();
            // var EdgeValues = new List<List<double>>();
            // foreach (var val in new double[]{1, -1, 0}){
            //     EdgeValues.Add(new List<double>{val});
            // }
            // OpenFoamPatchField cPtch = new(grd, 1, new int[]{1,2,3}, new string[]{"dirichlet","dirichlet","neumann"}, new double[]{1,-1,0});
            int i = 0;
            // foreach (var grd in new List<OpenFOAMGrid>{smallGrd, mediumGrd, largeGrd}){
            foreach (var grd in new List<OpenFOAMGrid>{smallGrd, mediumGrd}){
                OpenFoamDGField f = new OpenFoamDGField(grd, 2, 2);
                OpenFoamMatrix mtx = new OpenFoamMatrix(grd, f);
                OpenFoamPatchField cPtch;
                int[] safeEts = new int[] { 1, 2, 3 };
                // int* ets = (int*)safeEts[0];
                string[] safeEtyps = new string[] { "neumann", "neumann", "neumann" };
                // int* eTyps = (int*)safeEtyps[0];
                // double[] safeVals = new double[]{ 1.0, -1.0, 0 };
                double[] safeVals = new double[] { 0.0, -0.0, 0 };
                cPtch = new OpenFoamPatchField(grd, 1, safeEts, safeEtyps, safeVals);

                double[] safeValsU = new double[] { 1.0, -1.0, 0, 0, 0, 0, 0, 0, 0 };
                OpenFoamPatchField uPtch = new OpenFoamPatchField(grd, 3, safeEts, safeEtyps, safeValsU);
                // unsafe {
                //     int[] safeEts = new[]{1, 2, 3};
                //     // int* ets = (int*)safeEts[0];
                //     int[] safeEtyps = new[]{1, 1, 0};
                //     // int* eTyps = (int*)safeEtyps[0];
                //     double[] safeVals = new[]{1.0, -1.0, 0};
                //     fixed (double* eVals = safeVals){
                //         fixed (int* ets = safeEts)
                //         {
                //             fixed (int* eTyps = safeEtyps)
                //             {
                //                 cPtch = new(grd, 1, 3, ets, eTyps, eVals);
                //             }
                //         }
                //     }
                // }
                OpenFoamDGField U = new OpenFoamDGField(grd, 2, 3);
                // Console.ReadLine();

                int noOfTotalCells = grd.GridData.Grid.NumberOfCells;
                ScalarFunction func()
                {
                    // double rMin = 2.0e-3 / Math.Sqrt(noOfTotalCells) * 3.0 / Math.Sqrt(2);
                    double radius = 0.5e-3;
                    // double radius = rMin * 1.3;
                    return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x - 1.0e-3, 2) + Math.Pow(z - 0.0e-3, 2)) + Math.Pow(radius, 1)) * 50000))).Vectorize();
                    // return ((_3D)((x, y, z) => tanh(((x - 0.0011) + 0.01 * z)*3750))).Vectorize();
                    // return ((_3D)((x, y, z) => Math.Tanh(((x - 0.0011) + 0.01 * z) * 5500))).Vectorize();
                    // return ((_3D)((x, y, z) => tanh(((x - 2.5) + 0.1 * (y - 2.5))/1))).Vectorize();
                }

                double preNorm = chOp.Norm(mtx, func());
                chOp.CahnHilliard(mtx, U, cPtch, uPtch, func());
                double postNorm = chOp.Norm();
                double normRelChange = Math.Abs((postNorm-preNorm)/preNorm);
                double jumpNorm = chOp.JumpNorm();
                preNorms.Add(preNorm);
                norms.Add(postNorm);
                jumpNorms.Add(jumpNorm);
                normRelChanges.Add(normRelChange);

                // move all plt files into their own directory before starting the next calculation
                var source = new System.IO.DirectoryInfo("./");
                System.IO.FileInfo[] files = source.GetFiles("*.plt");
                foreach (var file in files) {
                    file.MoveTo(new List<string>{"./small/", "./medium/", "./large/"}[i] + file.Name, true);
                }
                i++;
            }
            Console.WriteLine("preNorms:");
            preNorms.ForEach(i => Console.WriteLine("{0}", i));
            Console.WriteLine("postNorms:");
            norms.ForEach(i => Console.WriteLine("{0}", i));
            Console.WriteLine("normrelchanges:");
            normRelChanges.ForEach(i => Console.WriteLine("{0}", i));
            Console.WriteLine("jumpNorms:");
            jumpNorms.ForEach(i => Console.WriteLine("{0}", i));

            // make sure it works better with a finer grid
            Assert.IsTrue(normRelChanges[0] > normRelChanges[1]);
            Assert.IsTrue(jumpNorms[0] > jumpNorms[1]);

            // also have some absolute constraints in place
            Assert.IsTrue(normRelChanges[1] < 1e-2);
            Assert.IsTrue(jumpNorms[1] < 1e-3);

            Cleanup();

        }

        static Initializer MyInit;

        /// <summary>
        /// MPI Init
        /// </summary>
        public static void Init() {
            MyInit = new Initializer();
            MyInit.BoSSSInitialize();
        }

        /// <summary>
        /// MPI shutdown
        /// </summary>
        public static void Cleanup() {
            MyInit.BoSSSFinalize();
        }
    }
}
