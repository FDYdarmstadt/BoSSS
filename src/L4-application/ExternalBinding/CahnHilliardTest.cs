using ilPSP;
using System;
using System.Linq;
using System.IO;
using System.Reflection;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Statistic;
using BoSSS.Solution.ASCIIExport;
using ilPSP.Utils;
using NUnit.Framework;

namespace BoSSS.Application.ExternalBinding {

    [TestFixture]
    static public class CahnHilliardTest {

        static double LogLogRegression(IEnumerable<int> _xValues, IEnumerable<double> _yValues)
        {
            double[] xValues = _xValues.Select(x => Math.Log10(x)).ToArray();
            double[] yValues = _yValues.Select(y => Math.Log10(y)).ToArray();

            double xAvg = xValues.Average();
            double yAvg = yValues.Average();

            double v1 = 0.0;
            double v2 = 0.0;

            for (int i = 0; i < yValues.Length; i++)
            {
                v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
                v2 += Math.Pow(xValues[i] - xAvg, 2);
            }

            double a = v1 / v2;
            double b = yAvg - a * xAvg;

            return a;
        }

        public static SinglePhaseField RunDropletTest(string GridPath, string PlotTargetDir = "./plots/", FixedOperators chOp = null) {
            Init();
            //GridImportTest.ConvertFOAMGrid();
            Console.WriteLine("Running Cahn-Hilliard Droplet Test");
            if (chOp == null) {
                chOp = new FixedOperators();
            }
            OpenFOAMGrid grd = GridImportFromDirectory.GenerateFOAMGrid(GridPath);
            OpenFoamDGField f = new OpenFoamDGField(grd, 2, 2);
            OpenFoamMatrix mtx = new OpenFoamMatrix(grd, f);
            OpenFoamPatchField cPtch;
            int[] safeEts = new int[] { 1, 2, 3, 4, 5 };
            string[] safeEtyps = new string[] { "neumann", "neumann", "neumann", "neumann", "neumann" };
            double[] safeVals = new double[] { 0.0, 0.0, 0.0 , 0.0 , 0.0 };
            cPtch = new OpenFoamPatchField(grd, 1, safeEts, safeEtyps, safeVals);

            double[] safeValsU = new double[] { 0.01*15, 0.01*(-15), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            OpenFoamPatchField uPtch = new OpenFoamPatchField(grd, 3, safeEts, safeEtyps, safeValsU);
            OpenFoamDGField U = new OpenFoamDGField(grd, 2, 3);

            int noOfTotalCells = grd.GridData.Grid.NumberOfCells;
            // ScalarFunction func()
            // {
            //     double radius = 7;
            //     // return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x - 1.0e-3, 2) + Math.Pow(z - 0.0e-3, 2)) + Math.Pow(radius, 1)) * 50000))).Vectorize();
            //     return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x, 2) + Math.Pow(z, 2)) + Math.Pow(radius, 1)) * Math.Sqrt(2)))).Vectorize();
            // }

            var _chParams = new CahnHilliardParameters(_cahn: 0.1, _diffusion: 0.1, _stationary: false, _dt: 0.2, _endT: 0.2*1.1);
            chOp.CahnHilliardInternal(mtx, null, U, cPtch, uPtch, chParams: _chParams);

            var field = new SinglePhaseField(mtx.ColMap.BasisS[0], "c");
            field.Acc(1.0, mtx.Fields[0].Fields[0] as SinglePhaseField);

            // move all plt files into their own directory before starting the next calculation
            var source = new System.IO.DirectoryInfo("./");
            System.IO.FileInfo[] files = source.GetFiles("*.plt");
            var targetPath = PlotTargetDir;
            System.IO.Directory.CreateDirectory(targetPath);
            foreach (var file in files)
            {
                try {
                    file.MoveTo(targetPath + "/" + file.Name);
                } catch (Exception e){}
            }
            return field;
        }

        [NUnitFileToCopyHack(
        "OpenFOAM/meshes/1D/small/polyMesh/*",

        "OpenFOAM/meshes/1D/medium/polyMesh/*",

        "OpenFOAM/meshes/1D/large/polyMesh/*"
            )]
        [Test]
        public static void Test1D() {
            Init();

            string smallGrd = "./meshes/1D/small/polyMesh/";
            string mediumGrd = "./meshes/1D/medium/polyMesh/";
            string largeGrd = "./meshes/1D/large/polyMesh/";
            var chOp = new FixedOperators();

            // OpenFOAMGrid grd = GridImportFromDirectory.GenerateFOAMGrid(smallGrd);
            List<int> DOFs = new List<int>();
            List<double> errors = new List<double>();
            int i = 0;
            foreach (var grdStr in new List<string>{smallGrd, mediumGrd, largeGrd}){
            // foreach (var grdStr in new List<string>{mediumGrd, largeGrd}){
                OpenFOAMGrid grd = GridImportFromDirectory.GenerateFOAMGrid(grdStr);
                OpenFoamDGField f = new OpenFoamDGField(grd, 2, 2);
                OpenFoamMatrix mtx = new OpenFoamMatrix(grd, f);
                OpenFoamPatchField cPtch;
                int[] safeEts = new int[] { 1, 2, 3 };
                string[] safeEtyps = new string[] { "neumann", "neumann", "neumann" };
                double[] safeVals = new double[] { 0, 0, 0 };
                cPtch = new OpenFoamPatchField(grd, 1, safeEts, safeEtyps, safeVals);

                double[] safeValsU = new double[] { 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0 };
                OpenFoamPatchField uPtch = new OpenFoamPatchField(grd, 3, safeEts, safeEtyps, safeValsU);
                OpenFoamDGField U = new OpenFoamDGField(grd, 2, 3);

                int noOfTotalCells = grd.GridData.Grid.NumberOfCells;

                double cahn = 1.0e-1;

                ScalarFunction func()
                {
                    return ((_3D)((x, y, z) => Math.Sign(x))).Vectorize();
                }
                ScalarFunction Ufunc()
                {
                    return ((_3D)((x, y, z) => 0)).Vectorize();
                }

                var _chParams = new CahnHilliardParameters(_cahn: cahn, _diffusion: 1.0, _stationary: true);
                chOp.CahnHilliardInternal(mtx, null, U, cPtch, uPtch, func: func(), Ufunc: Ufunc(), chParams: _chParams);

                var field = new SinglePhaseField(mtx.ColMap.BasisS[0], "c");
                field.Acc(1.0, mtx.Fields[0].Fields[0] as SinglePhaseField);
                var other = new SinglePhaseField(mtx.ColMap.BasisS[0], "c");
                other.ProjectField(((x, y, z) => Math.Tanh((x) / (Math.Sqrt(2) * cahn))));
                double error = field.L2Error(other);
                errors.Add(error);
                DOFs.Add(grd.NumberOfCells);

                // move all plt files into their own directory before starting the next calculation
                var source = new System.IO.DirectoryInfo("./");
                System.IO.FileInfo[] files = source.GetFiles("*.plt");
                var targetPath = new List<string>{"./small1D/", "./medium1D/", "./large1D/"}[i];
                System.IO.Directory.CreateDirectory(targetPath);
                foreach (var file in files) {
                    string targetFile = targetPath + file.Name;
                    System.IO.File.Delete(targetFile);
                    file.MoveTo(targetFile);
                }

                // export 1D data to csv
                var csvExp = new CSVExportDriver(grd.GetGridData(), false, false, 3);
                csvExp.PlotFields(targetPath + "out", 1e5, field);
                i++;
            }
            Console.WriteLine("L2Errors: ");
            errors.ForEach(j => Console.WriteLine("{0}", j));

            var slope = -LogLogRegression(DOFs, errors);
            Console.WriteLine("Slope: " + slope);
            Console.WriteLine("DOFs[0]: " + DOFs[0]);
            Console.WriteLine("DOFs[1]: " + DOFs[1]);
            Console.WriteLine("DOFs[2]: " + DOFs[2]);
            Console.WriteLine("errorsS[0]: " + errors[0]);
            Console.WriteLine("errorsS[1]: " + errors[1]);
            Console.WriteLine("errorsS[2]: " + errors[2]);
            Assert.IsTrue(slope >= 2.99);
        }


        [NUnitFileToCopyHack(
        "OpenFOAM/meshes/big/small/polyMesh/*",

        "OpenFOAM/meshes/big/medium/polyMesh/*",

        "OpenFOAM/meshes/big/large/polyMesh/*"
            )]
        [Test]
        public static void DropletTest() {

            Init();

            // GridImportTest.ConvertFOAMGrid();

            Console.WriteLine("Running Cahn-Hilliard Droplet Test");
            // OpenFOAMGrid grd = GridImportTestSmall.GenerateFOAMGrid();
            string smallGrd = "./meshes/big/small/polyMesh/";
            string mediumGrd = "./meshes/big/medium/polyMesh/";
            string largeGrd = "./meshes/big/large/polyMesh/";
            var chOp = new FixedOperators();
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
            foreach (var grd in new List<string>{smallGrd, mediumGrd}){
            // i++;
            // foreach (var grd in new List<string>{mediumGrd}){

                RunDropletTest(grd, new List<string>{"./small/", "./medium/", "./large/"}[i], chOp);
                double normRelChange = chOp.NormRelChange();
                double jumpNorm = chOp.JumpNorm();
                normRelChanges.Add(normRelChange);
                jumpNorms.Add(jumpNorm);

                // move all plt files into their own directory before starting the next calculation
                var source = new System.IO.DirectoryInfo("./");
                System.IO.FileInfo[] files = source.GetFiles("*.plt");
                var targetPath = new List<string>{"./small/", "./medium/", "./large/"}[i];
                System.IO.Directory.CreateDirectory(targetPath);
                // foreach (var file in files) {
                //     file.MoveTo(targetPath + file.Name, true);
                // }
                i++;
            }

            Console.WriteLine("normrelchanges:");
            normRelChanges.ForEach(j => Console.WriteLine("{0}", j));
            Console.WriteLine("jumpNorms:");
            jumpNorms.ForEach(j => Console.WriteLine("{0}", j));

            // make sure it works better with a finer grid
            // Assert.IsTrue(normRelChanges[0] > normRelChanges[1]);
            Assert.IsTrue(jumpNorms[0] > jumpNorms[1]);
            // Assert.IsTrue(normRelChanges[1] > normRelChanges[2]);
            // Assert.IsTrue(jumpNorms[1] > jumpNorms[2]);

            // also have some absolute constraints in place
            Assert.IsTrue(normRelChanges[1] < 1e-2);
            // Assert.IsTrue(jumpNorms[2] < 1e-3);

            Cleanup();

        }

#if !DEBUG // this test takes too much time in debug mode
        [NUnitFileToCopyHack(
        "OpenFOAM/meshes/big/small/polyMesh/*",

        "OpenFOAM/meshes/big/medium/polyMesh/*",

        "OpenFOAM/meshes/big/large/polyMesh/*"
            )]
        [Test]
#endif
        public static void ConvergenceTest() {

            Console.WriteLine("Running Cahn-Hilliard Convergence Test");
            var chOp = new FixedOperators();
            // OpenFOAMGrid grd = GridImportTestSmall.GenerateFOAMGrid();
            string currentDirectory = "";
            string smallGrd = currentDirectory + "./meshes/big/small/polyMesh/";
            string mediumGrd = currentDirectory + "./meshes/big/medium/polyMesh/";
            string largeGrd = currentDirectory + "./meshes/big/large/polyMesh/";
            // string smallGrd = "/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dummyConvAnalysis/small/constant/polyMesh/";
            // string mediumGrd = "/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dummyConvAnalysis/medium/constant/polyMesh/";
            // string largeGrd = "/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dummyConvAnalysis/large/constant/polyMesh/";
            int i = 0;
            List<IEnumerable<DGField>> solutionOnDifferentResolutions = new List<IEnumerable<DGField>>();
            List<DGField> solutionOnDifferentResolutions2 = new List<DGField>();
            foreach (var grd in new List<string>{smallGrd, mediumGrd, largeGrd}){
            // foreach (var grd in new List<OpenFOAMGrid>{smallGrd, mediumGrd}){
                var field = RunDropletTest(grd, new List<string>{"./small/", "./medium/", "./large/"}[i], chOp);
                solutionOnDifferentResolutions.Add(new DGField[]{field});
                solutionOnDifferentResolutions2.Add(field);
                i++;
            }
            DGFieldComparison.ComputeErrors(
                solutionOnDifferentResolutions, out var hS, out var DOFs, out var errorS, NormType.L2_embedded);
            var slope = (Math.Log(errorS["c"][1]) - Math.Log(errorS["c"][0]))/(Math.Log(Math.Sqrt(1.0/DOFs["c"][1])) - Math.Log(Math.Sqrt(1.0/DOFs["c"][0])));
            Console.WriteLine("Slope: " + slope);
            Console.WriteLine("DOFs[0]: " + DOFs["c"][0]);
            Console.WriteLine("DOFs[1]: " + DOFs["c"][1]);
            Console.WriteLine("errorsS[0]: " + errorS["c"][0]);
            Console.WriteLine("errorsS[1]: " + errorS["c"][1]);
            Assert.IsTrue(slope >= 2.4); // TODO should this be larger than 3?

            Cleanup();
        }

        public static void Main() {

            // GridImportTest.ConvertFOAMGrid();
            Init();

            string examplesDirectory   = "../../../../../../examples/OpenFOAM/meshes/big/";
            string examplesDirectory1D = "../../../../../../examples/OpenFOAM/meshes/1D/";
            bool overwrite = true;
            try {
                Directory.CreateDirectory("./meshes/");
                Directory.CreateDirectory("./meshes/big/");
                Directory.CreateDirectory("./meshes/big/small/");
                Directory.CreateDirectory("./meshes/big/small/polyMesh/");
                Directory.CreateDirectory("./meshes/big/medium/");
                Directory.CreateDirectory("./meshes/big/medium/polyMesh/");
                Directory.CreateDirectory("./meshes/big/large/");
                Directory.CreateDirectory("./meshes/big/large/polyMesh/");

                File.Copy(Path.Combine(examplesDirectory, "small/polyMesh/boundarySmall"), "./meshes/big/small/polyMesh/boundary", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "small/polyMesh/facesSmall"), "./meshes/big/small/polyMesh/faces", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "small/polyMesh/neighbourSmall"), "./meshes/big/small/polyMesh/neighbour", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "small/polyMesh/ownerSmall"), "./meshes/big/small/polyMesh/owner", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "small/polyMesh/pointsSmall"), "./meshes/big/small/polyMesh/points", overwrite);

                File.Copy(Path.Combine(examplesDirectory, "medium/polyMesh/boundaryMedium"), "./meshes/big/medium/polyMesh/boundary", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "medium/polyMesh/facesMedium"), "./meshes/big/medium/polyMesh/faces", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "medium/polyMesh/neighbourMedium"), "./meshes/big/medium/polyMesh/neighbour", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "medium/polyMesh/ownerMedium"), "./meshes/big/medium/polyMesh/owner", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "medium/polyMesh/pointsMedium"), "./meshes/big/medium/polyMesh/points", overwrite);

                File.Copy(Path.Combine(examplesDirectory, "large/polyMesh/boundaryLarge"), "./meshes/big/large/polyMesh/boundary", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "large/polyMesh/facesLarge"), "./meshes/big/large/polyMesh/faces", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "large/polyMesh/neighbourLarge"), "./meshes/big/large/polyMesh/neighbour", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "large/polyMesh/ownerLarge"), "./meshes/big/large/polyMesh/owner", overwrite);
                File.Copy(Path.Combine(examplesDirectory, "large/polyMesh/pointsLarge"), "./meshes/big/large/polyMesh/points", overwrite);
            } catch (Exception e) { Console.WriteLine("Continuing despite " + e); }

            try {
                Directory.CreateDirectory("./meshes/");
                Directory.CreateDirectory("./meshes/1D/");
                Directory.CreateDirectory("./meshes/1D/small/");
                Directory.CreateDirectory("./meshes/1D/small/polyMesh/");
                Directory.CreateDirectory("./meshes/1D/medium/");
                Directory.CreateDirectory("./meshes/1D/medium/polyMesh/");
                Directory.CreateDirectory("./meshes/1D/large/");
                Directory.CreateDirectory("./meshes/1D/large/polyMesh/");

                File.Copy(Path.Combine(examplesDirectory1D, "small/polyMesh/boundarySmall1D"), "./meshes/1D/small/polyMesh/boundary", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "small/polyMesh/facesSmall1D"), "./meshes/1D/small/polyMesh/faces", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "small/polyMesh/neighbourSmall1D"), "./meshes/1D/small/polyMesh/neighbour", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "small/polyMesh/ownerSmall1D"), "./meshes/1D/small/polyMesh/owner", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "small/polyMesh/pointsSmall1D"), "./meshes/1D/small/polyMesh/points", overwrite);

                File.Copy(Path.Combine(examplesDirectory1D, "medium/polyMesh/boundaryMedium1D"), "./meshes/1D/medium/polyMesh/boundary", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "medium/polyMesh/facesMedium1D"), "./meshes/1D/medium/polyMesh/faces", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "medium/polyMesh/neighbourMedium1D"), "./meshes/1D/medium/polyMesh/neighbour", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "medium/polyMesh/ownerMedium1D"), "./meshes/1D/medium/polyMesh/owner", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "medium/polyMesh/pointsMedium1D"), "./meshes/1D/medium/polyMesh/points", overwrite);

                File.Copy(Path.Combine(examplesDirectory1D, "large/polyMesh/boundaryLarge1D"), "./meshes/1D/large/polyMesh/boundary", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "large/polyMesh/facesLarge1D"), "./meshes/1D/large/polyMesh/faces", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "large/polyMesh/neighbourLarge1D"), "./meshes/1D/large/polyMesh/neighbour", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "large/polyMesh/ownerLarge1D"), "./meshes/1D/large/polyMesh/owner", overwrite);
                File.Copy(Path.Combine(examplesDirectory1D, "large/polyMesh/pointsLarge1D"), "./meshes/1D/large/polyMesh/points", overwrite);
            } catch (Exception e) { Console.WriteLine("Continuing despite " + e); }

            string currentDirectory = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);
            string grd = Path.Combine(currentDirectory, "meshes", "big", "medium", "polyMesh");

            // RunDropletTest(grd);

            // DropletTest();
            // ConvergenceTest();
            Test1D();

        }

        static Initializer MyInit;

        /// <summary>
        /// MPI Init
        /// </summary>
        public static void Init() {
            MyInit = new Initializer();
            MyInit.BoSSSInitialize();

            // recover directory structure
            try {
                Directory.CreateDirectory("./meshes/");
                Directory.CreateDirectory("./meshes/big/");
                Directory.CreateDirectory("./meshes/big/small/");
                Directory.CreateDirectory("./meshes/big/small/polyMesh/");
                Directory.CreateDirectory("./meshes/big/medium/");
                Directory.CreateDirectory("./meshes/big/medium/polyMesh/");
                Directory.CreateDirectory("./meshes/big/large/");
                Directory.CreateDirectory("./meshes/big/large/polyMesh/");
                File.Copy("boundarySmall", "./meshes/big/small/polyMesh/boundary");
                File.Copy("facesSmall", "./meshes/big/small/polyMesh/faces");
                File.Copy("neighbourSmall", "./meshes/big/small/polyMesh/neighbour");
                File.Copy("ownerSmall", "./meshes/big/small/polyMesh/owner");
                File.Copy("pointsSmall", "./meshes/big/small/polyMesh/points");

                File.Copy("boundaryMedium", "./meshes/big/medium/polyMesh/boundary");
                File.Copy("facesMedium", "./meshes/big/medium/polyMesh/faces");
                File.Copy("neighbourMedium", "./meshes/big/medium/polyMesh/neighbour");
                File.Copy("ownerMedium", "./meshes/big/medium/polyMesh/owner");
                File.Copy("pointsMedium", "./meshes/big/medium/polyMesh/points");

                File.Copy("boundaryLarge", "./meshes/big/large/polyMesh/boundary");
                File.Copy("facesLarge", "./meshes/big/large/polyMesh/faces");
                File.Copy("neighbourLarge", "./meshes/big/large/polyMesh/neighbour");
                File.Copy("ownerLarge", "./meshes/big/large/polyMesh/owner");
                File.Copy("pointsLarge", "./meshes/big/large/polyMesh/points");

            } catch (Exception e) { Console.WriteLine("Continuing despite " + e); }

            try {
                Directory.CreateDirectory("./meshes/");
                Directory.CreateDirectory("./meshes/1D/");
                Directory.CreateDirectory("./meshes/1D/small/");
                Directory.CreateDirectory("./meshes/1D/small/polyMesh/");
                Directory.CreateDirectory("./meshes/1D/medium/");
                Directory.CreateDirectory("./meshes/1D/medium/polyMesh/");
                Directory.CreateDirectory("./meshes/1D/large/");
                Directory.CreateDirectory("./meshes/1D/large/polyMesh/");

                File.Copy("boundarySmall1D", "./meshes/1D/small/polyMesh/boundary");
                File.Copy("facesSmall1D", "./meshes/1D/small/polyMesh/faces");
                File.Copy("neighbourSmall1D", "./meshes/1D/small/polyMesh/neighbour");
                File.Copy("ownerSmall1D", "./meshes/1D/small/polyMesh/owner");
                File.Copy("pointsSmall1D", "./meshes/1D/small/polyMesh/points");

                File.Copy("boundaryMedium1D", "./meshes/1D/medium/polyMesh/boundary");
                File.Copy("facesMedium1D", "./meshes/1D/medium/polyMesh/faces");
                File.Copy("neighbourMedium1D", "./meshes/1D/medium/polyMesh/neighbour");
                File.Copy("ownerMedium1D", "./meshes/1D/medium/polyMesh/owner");
                File.Copy("pointsMedium1D", "./meshes/1D/medium/polyMesh/points");

                File.Copy("boundaryLarge1D", "./meshes/1D/large/polyMesh/boundary");
                File.Copy("facesLarge1D", "./meshes/1D/large/polyMesh/faces");
                File.Copy("neighbourLarge1D", "./meshes/1D/large/polyMesh/neighbour");
                File.Copy("ownerLarge1D", "./meshes/1D/large/polyMesh/owner");
                File.Copy("pointsLarge1D", "./meshes/1D/large/polyMesh/points");

            } catch (Exception e) { Console.WriteLine("Continuing despite " + e); }
        }

        /// <summary>
        /// MPI shutdown
        /// </summary>
        public static void Cleanup() {
            MyInit.BoSSSFinalize();
        }
    }
}
