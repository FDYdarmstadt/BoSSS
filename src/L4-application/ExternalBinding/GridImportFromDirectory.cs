using ilPSP;
using BoSSS.Foundation.Grid;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Linq;
using System.Text.RegularExpressions;

namespace BoSSS.Application.ExternalBinding {


    /// <summary>
    /// Basic Testing for external language binding.
    /// </summary>
    [TestFixture]
    static public class GridImportFromDirectory {

        // static string polyMeshDir = "/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dropletInShearFlowMedium/constant/polyMesh/";
        // static string polyMeshDir = "/home/klingenberg/Documents-work/programming/foam-dg/foam-dg/run/dropletInShearFlowSmall/constant/polyMesh/";

        static void RemoveComments(ref string input){
            // see https://stackoverflow.com/questions/3524317/regex-to-strip-line-comments-from-c-sharp/3524689#3524689
            var blockComments = @"/\*(.*?)\*/";
            var lineComments = @"//(.*?)\r?\n";
            var strings = @"""((\\[^\n]|[^""\n])*)""";
            var verbatimStrings = @"@(""[^""]*"")+";

            string noComments = Regex.Replace(input,
                                          blockComments + "|" + lineComments + "|" + strings + "|" + verbatimStrings,
                                          me =>
                                          {
                                              if (me.Value.StartsWith("/*") || me.Value.StartsWith("//"))
                                                return me.Value.StartsWith("//") ? System.Environment.NewLine : "";
                                            // Keep the literal strings
                                            return me.Value;
                                        },
                                        RegexOptions.Singleline);

        input = noComments;
        }

        static void RemoveHeader(ref string input){
            var foamHeader = @"FoamFile.*?{(.*?)}";
            // var foamHeader = @"FoamFile\n{";
            string noHeader = Regex.Replace(input,
                                            foamHeader,
                                            "",
                                            RegexOptions.Singleline);

        input = noHeader;
        }

        static void RemoveEmptyLines(ref string input){
            var pat = @"^\n";
            string noEL = Regex.Replace(input,
                                            pat,
                                            "");

        input = noEL;
        }

        static void Clean(ref string input){
            RemoveComments(ref input);
            RemoveHeader(ref input);
            RemoveEmptyLines(ref input);
        }

        internal static int[][] getNestedIntArray(string filename, string polyMeshDir) {
            var ret = new List<int[]>();
            string text = File.ReadAllText(polyMeshDir + filename);
            Clean(ref text);
            string[] lines = text.Split('\n');
            var pat = @"^[0-9]*\((.*?)\)";
            foreach (var line in lines){
                foreach (Match match in Regex.Matches(line, pat)){
                    string[] strArr = match.Groups[1].Value.Split();
                    int len = strArr.Length;
                    int[] arr = new int[len];
                    for (int i = 0; i < len; i++){
                        arr[i] = int.Parse(strArr[i]);
                    }
                    ret.Add(arr);
                }
            }
            return ret.ToArray();
        }

        internal static int[] getIntArray(string filename, string polyMeshDir) {
            var ret = new List<int>();
            string text = File.ReadAllText(polyMeshDir + filename);
            Clean(ref text);
            text = Regex.Replace(text, "\n", " ");
            var pat = @"[0-9]*\s*\((.*?)\)";
            foreach (Match match in Regex.Matches(text, pat, RegexOptions.Singleline)){
                string[] strArr = match.Groups[1].Value.Split(' ');
                int len = strArr.Length;
                int[] arr = new int[len];
                foreach (var elem in strArr){
                    try{
                        ret.Add(int.Parse(elem));
                    } catch (System.FormatException e){}
                }
            }
            return ret.ToArray();
        }

        internal static double[,] getNestedDoubleArray(string filename, string polyMeshDir) {
            var ret = new List<double[]>();
            string text = File.ReadAllText(polyMeshDir + filename);
            Clean(ref text);
            string[] lines = text.Split('\n');
            var pat = @"^[0-9]*\((.*?)\)";
            int dim = 3;
            foreach (var line in lines){
                foreach (Match match in Regex.Matches(line, pat)){
                    string[] strArr = match.Groups[1].Value.Split();
                    dim = strArr.Length;
                    double[] arr = new double[dim];
                    for (int i = 0; i < dim; i++){
                        arr[i] = double.Parse(strArr[i]);
                    }
                    ret.Add(arr);
                }
            }
            int len = ret.Count;
            double[,] retArr = new double[len, dim];
            for (int i = 0; i < len; i++){
                for (int j = 0; j < dim; j++){
                    retArr[i,j] = ret[i][j];
                }
            }
            return retArr;
        }


        internal static int[][] getFaces(string polyMeshDir) {
            return getNestedIntArray("faces", polyMeshDir);
        }

        internal static int[] getNeighbour(string polyMeshDir) {
            return getIntArray("neighbour", polyMeshDir);
        }

        internal static int[] getOwner(string polyMeshDir){
            return getIntArray("owner", polyMeshDir);
        }

        internal static double[,] getPoints(string polyMeshDir){
            return getNestedDoubleArray("points", polyMeshDir);
        }

        internal static string[] names = new string[] {"left", "right", "empty"};

        internal static int[] getPatchIDs(int len){
            int[] ret = new int[len];
            ret[0] = 0;
            ret[1] = 1;
            for (int i = 2; i < len; i++){
                ret[i] = 2;
            }
            return ret;
        }
        

        public static void GridImportMain() {
            
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
        [NUnitFileToCopyHack(
        "OpenFOAM/meshes/big/small/polyMesh/boundarySmall",
        "OpenFOAM/meshes/big/small/polyMesh/facesSmall",
        "OpenFOAM/meshes/big/small/polyMesh/neighbourSmall",
        "OpenFOAM/meshes/big/small/polyMesh/ownerSmall",
        "OpenFOAM/meshes/big/small/polyMesh/pointsSmall",

        "OpenFOAM/meshes/big/medium/polyMesh/boundaryMedium",
        "OpenFOAM/meshes/big/medium/polyMesh/facesMedium",
        "OpenFOAM/meshes/big/medium/polyMesh/neighbourMedium",
        "OpenFOAM/meshes/big/medium/polyMesh/ownerMedium",
        "OpenFOAM/meshes/big/medium/polyMesh/pointsMedium",

        "OpenFOAM/meshes/big/large/polyMesh/boundaryLarge",
        "OpenFOAM/meshes/big/large/polyMesh/facesLarge",
        "OpenFOAM/meshes/big/large/polyMesh/neighbourLarge",
        "OpenFOAM/meshes/big/large/polyMesh/ownerLarge",
        "OpenFOAM/meshes/big/large/polyMesh/pointsLarge"
            )]
        [Test]
        public static void ConvertFOAMGrid() {

            // string currentDirectory = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);
            string polyMeshDir = "./meshes/big/small/polyMesh/";
            var owner = getOwner(polyMeshDir);
            var neighbour = getNeighbour(polyMeshDir);
            int nCells = Max(owner.Max(), neighbour.Max()) + 1;
            var g = GenerateFOAMGrid(polyMeshDir);

            Assert.AreEqual(g.GridData.iLogicalCells.Count, nCells, "Mismatch in expected number of cells.");
        }

        public static OpenFOAMGrid GenerateFOAMGrid(string polyMeshDir) {

            var faces = getFaces(polyMeshDir);
            var owner = getOwner(polyMeshDir);
            var neighbour = getNeighbour(polyMeshDir);
            var points = getPoints(polyMeshDir);
            int nCells = Max(owner.Max(), neighbour.Max()) + 1;
            var patchIDs = getPatchIDs(faces.Length);
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

            // recover directory structure
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
