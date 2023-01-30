using System;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Statistic;
using ilPSP.Utils;
using NUnit.Framework;

namespace BoSSS.Application.ExternalBinding {

    [TestFixture]
    static public class CahnHilliardTest {

        // [Test]
        public static void DropletTest() {

            Init();

            GridImportTest.ConvertFOAMGrid();
            Console.WriteLine("Running Cahn-Hilliard Test");
            var chOp = new FixedOperators();
            // OpenFOAMGrid grd = GridImportTestSmall.GenerateFOAMGrid();
            OpenFOAMGrid smallGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/small/polyMesh/");
            OpenFOAMGrid mediumGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/medium/polyMesh/");
            OpenFOAMGrid largeGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/large/polyMesh/");
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
                string[] safeEtyps = new string[] { "neumann", "neumann", "neumann" };
                double[] safeVals = new double[] { 0.0, 0.0, 0 };
                cPtch = new OpenFoamPatchField(grd, 1, safeEts, safeEtyps, safeVals);

                double[] safeValsU = new double[] { 1.0, -1.0, 0, 0, 0, 0, 0, 0, 0 };
                OpenFoamPatchField uPtch = new OpenFoamPatchField(grd, 3, safeEts, safeEtyps, safeValsU);
                OpenFoamDGField U = new OpenFoamDGField(grd, 2, 3);

                int noOfTotalCells = grd.GridData.Grid.NumberOfCells;
                ScalarFunction func()
                {
                    double radius = 7;
                    // return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x - 1.0e-3, 2) + Math.Pow(z - 0.0e-3, 2)) + Math.Pow(radius, 1)) * 50000))).Vectorize();
                    return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x, 2) + Math.Pow(z, 2)) + Math.Pow(radius, 1)) * Math.Sqrt(2)))).Vectorize();
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
                var targetPath = new List<string>{"./small/", "./medium/", "./large/"}[i];
                System.IO.Directory.CreateDirectory(targetPath);
                foreach (var file in files) {
                    file.MoveTo(targetPath + file.Name, true);
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
            Assert.IsTrue(normRelChanges[1] > normRelChanges[2]);
            Assert.IsTrue(jumpNorms[1] > jumpNorms[2]);

            // also have some absolute constraints in place
            Assert.IsTrue(normRelChanges[2] < 1e-2);
            Assert.IsTrue(jumpNorms[2] < 1e-3);

            Cleanup();

        }

        [Test]
        public static void ConvergenceTest() {

            Init();

            GridImportTest.ConvertFOAMGrid();
            Console.WriteLine("Running Cahn-Hilliard Test");
            var chOp = new FixedOperators();
            // OpenFOAMGrid grd = GridImportTestSmall.GenerateFOAMGrid();
            OpenFOAMGrid smallGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/small/polyMesh/");
            OpenFOAMGrid mediumGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/medium/polyMesh/");
            OpenFOAMGrid largeGrd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/large/polyMesh/");
            // OpenFOAMGrid smallGrd = GridImportFromDirectory.GenerateFOAMGrid("/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dummyConvAnalysis/small/constant/polyMesh/");
            // OpenFOAMGrid mediumGrd = GridImportFromDirectory.GenerateFOAMGrid("/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dummyConvAnalysis/medium/constant/polyMesh/");
            // OpenFOAMGrid largeGrd = GridImportFromDirectory.GenerateFOAMGrid("/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dummyConvAnalysis/large/constant/polyMesh/");
            int i = 0;
            List<IEnumerable<DGField>> solutionOnDifferentResolutions = new List<IEnumerable<DGField>>();
            foreach (var grd in new List<OpenFOAMGrid>{smallGrd, mediumGrd, largeGrd}){
            // foreach (var grd in new List<OpenFOAMGrid>{smallGrd, mediumGrd}){
                OpenFoamDGField f = new OpenFoamDGField(grd, 2, 2);
                OpenFoamMatrix mtx = new OpenFoamMatrix(grd, f);
                OpenFoamPatchField cPtch;
                int[] safeEts = new int[] { 1, 2, 3 };
                string[] safeEtyps = new string[] { "neumann", "neumann", "neumann" };
                double[] safeVals = new double[] { 0.0, 0.0, 0 };
                cPtch = new OpenFoamPatchField(grd, 1, safeEts, safeEtyps, safeVals);

                double[] safeValsU = new double[] { 1.0, -1.0, 0, 0, 0, 0, 0, 0, 0 };
                OpenFoamPatchField uPtch = new OpenFoamPatchField(grd, 3, safeEts, safeEtyps, safeValsU);
                OpenFoamDGField U = new OpenFoamDGField(grd, 2, 3);

                int noOfTotalCells = grd.GridData.Grid.NumberOfCells;
                ScalarFunction func()
                {
                    double radius = 7;
                    // return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x - 1.0e-3, 2) + Math.Pow(z - 0.0e-3, 2)) + Math.Pow(radius, 1)) * 50000))).Vectorize();
                    return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x, 2) + Math.Pow(z, 2)) + Math.Pow(radius, 1)) * Math.Sqrt(2)))).Vectorize();
                }

                chOp.CahnHilliard(mtx, U, cPtch, uPtch, func());

                var field = new SinglePhaseField(mtx.ColMap.BasisS[0], "c");
                // DGField field = mtx.Fields[0].Fields[0] as SinglePhaseField;
                field.Acc(1.0, mtx.Fields[0].Fields[0] as SinglePhaseField);
                solutionOnDifferentResolutions.Add(new DGField[]{field});
                Console.WriteLine(field == null);

                // move all plt files into their own directory before starting the next calculation
                var source = new System.IO.DirectoryInfo("./");
                System.IO.FileInfo[] files = source.GetFiles("*.plt");
                var targetPath = new List<string>{"./small/", "./medium/", "./large/"}[i];
                System.IO.Directory.CreateDirectory(targetPath);
                foreach (var file in files) {
                    file.MoveTo(targetPath + file.Name, true);
                }
                i++;
            }
            DGFieldComparison.ComputeErrors(
                // solutionOnDifferentResolutions, out var hS, out var DOFs, out var errorS, NormType.L2_embedded);
                solutionOnDifferentResolutions, out var hS, out var DOFs, out var errorS, NormType.L2_approximate);
            var slope = (errorS["c"][1] - errorS["c"][0])/(hS[1] - hS[0]);
            Console.WriteLine("Slope: " + slope);
            Assert.IsTrue(slope >= 3.0);

            Cleanup();

        }

        public static void Main() {

            ConvergenceTest();
            // Init();

            // GridImportTest.ConvertFOAMGrid();
            // Console.WriteLine("Running Cahn-Hilliard Test");
            // var chOp = new FixedOperators();
            // // OpenFOAMGrid grd = GridImportTestSmall.GenerateFOAMGrid();
            // // OpenFOAMGrid grd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/large/polyMesh/");
            // OpenFOAMGrid grd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/medium/polyMesh/");
            // // OpenFOAMGrid grd = GridImportFromDirectory.GenerateFOAMGrid("./meshes/big/small/polyMesh/");
            // var norms = new List<double>();
            // var preNorms = new List<double>();
            // var normRelChanges = new List<double>();
            // var jumpNorms = new List<double>();
            // OpenFoamDGField f = new OpenFoamDGField(grd, 2, 2);
            // OpenFoamMatrix mtx = new OpenFoamMatrix(grd, f);
            // OpenFoamPatchField cPtch;
            // int[] safeEts = new int[] { 1, 2, 3 };
            // // int* ets = (int*)safeEts[0];
            // string[] safeEtyps = new string[] { "neumann", "neumann", "neumann" };
            // // string[] safeEtyps = new string[] { "dirichlet", "dirichlet", "neumann" };
            // // int* eTyps = (int*)safeEtyps[0];
            // double[] safeVals = new double[] { 0.0, 0.0, 0 };
            // // double[] safeVals = new double[] { -1.0, -1.0, 0 };
            // cPtch = new OpenFoamPatchField(grd, 1, safeEts, safeEtyps, safeVals);

            // double[] safeValsU = new double[] { 1.0, -1.0, 0, 0, 0, 0, 0, 0, 0 };
            // OpenFoamPatchField uPtch = new OpenFoamPatchField(grd, 3, safeEts, safeEtyps, safeValsU);
            // OpenFoamDGField U = new OpenFoamDGField(grd, 2, 3);

            // int noOfTotalCells = grd.GridData.Grid.NumberOfCells;
            // ScalarFunction func()
            // {
            //     double radius = 7;
            //     // double radius = rMin * 1.3;
            //     // return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x - 1.0e-3, 2) + Math.Pow(z - 0.0e-3, 2)) + Math.Pow(radius, 1)) * 50000))).Vectorize();
            //     return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x, 2) + Math.Pow(z - 0.0e-3, 2)) + Math.Pow(radius, 1)) * Math.Sqrt(2)))).Vectorize();
            // }

            // double preNorm = chOp.Norm(mtx, func());
            // chOp.CahnHilliard(mtx, U, cPtch, uPtch, func());
            // double postNorm = chOp.Norm();
            // double normRelChange = Math.Abs((postNorm - preNorm) / preNorm);
            // double jumpNorm = chOp.JumpNorm();
            // preNorms.Add(preNorm);
            // norms.Add(postNorm);
            // jumpNorms.Add(jumpNorm);
            // normRelChanges.Add(normRelChange);

            // // move all plt files into their own directory before starting the next calculation
            // // var source = new System.IO.DirectoryInfo("./");
            // // System.IO.FileInfo[] files = source.GetFiles("*.plt");
            // // var targetPath = "./big/";
            // // System.IO.Directory.CreateDirectory(targetPath);
            // // foreach (var file in files)
            // // {
            // //     file.MoveTo(targetPath + file.Name, true);
            // // }
            // // i++;

            // Console.WriteLine("preNorms:");
            // preNorms.ForEach(i => Console.WriteLine("{0}", i));
            // Console.WriteLine("postNorms:");
            // norms.ForEach(i => Console.WriteLine("{0}", i));
            // Console.WriteLine("normrelchanges:");
            // normRelChanges.ForEach(i => Console.WriteLine("{0}", i));
            // Console.WriteLine("jumpNorms:");
            // jumpNorms.ForEach(i => Console.WriteLine("{0}", i));

            // Cleanup();

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
