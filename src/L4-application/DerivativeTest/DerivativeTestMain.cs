/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.IO;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using NUnit.Framework;
using BoSSS.Solution.GridImport;
using MPI.Wrappers;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution;

namespace BoSSS.Application.DerivativeTest {

    [TestFixture]
    public class Tests {

        static int CHUNK_DATA_LIMIT_bkup;

        [TestFixtureSetUp]
        public static void SetUp() {
            bool MpiInit;
            CHUNK_DATA_LIMIT_bkup = Quadrature_Bulksize.CHUNK_DATA_LIMIT;
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out MpiInit);
        }

        [TestFixtureTearDown]
        public static void Cleanup() {
            //Console.Out.Dispose();
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        [Test]
#if DEBUG
        public static void DerivativeTest_BuildInGrid([Range(1, 10)] int gridCase, [Values(2,10000000)] int bulksize_limit, [Values(1024,1024*1024*128)] int cache_size) {
#else
        public static void DerivativeTest_BuildInGrid([Range(1, 14)] int gridCase, [Values(1, 500, 10000000)] int bulksize_limit, [Values(1024, 1024 * 1024 * 128)] int cache_size) {
#endif
            DerivativeTestMain.GRID_CASE = gridCase;
            DerivativeTestMain p = null;
            Quadrature_Bulksize.CHUNK_DATA_LIMIT = bulksize_limit;
            BoSSS.Foundation.Caching.Cache.MaxMem = cache_size;

            BoSSS.Solution.Application._Main(new string[0], true, null, delegate() {
                p = new DerivativeTestMain();
                return p;
            });

            Assert.IsTrue(p.m_passed);
        }

        /// <summary>
        /// Filenames of test grids.
        /// </summary>
        string[] m_testFiles {
            get {
                List<string> R = new List<string>();
                R.AddRange(Directory.GetFiles("../../TestGrids/", "*.msh").Select(f => Path.GetFileName(f)));
                R.AddRange(Directory.GetFiles("../../TestGrids/", "*.cgns").Select(f => Path.GetFileName(f)));

                // blacklist that isn't working
                R.Remove("Ringleb6th.msh");
                R.Remove("Ringleb5th.msh");
                R.Remove("Ringleb4th.msh");

                // Too much for DEBUG-configration.
#if DEBUG
                R.Remove("kubus.cgns");
                R.Remove("kubus.cgns");
                R.Remove("WallMountedCube.cgns");
                R.Remove("Tria_coarse3rd.msh");
                R.Remove("Tria_coarse2nd.msh");
                R.Remove("Ringleb3rd.msh");
                R.Remove("wedding2D_v16.cgns");
#endif

                return R.ToArray();
            }
        }


        /// <summary>
        /// Test using grids imported from gmsh/cgns
        /// </summary>
        [Test]
        public static void DerivativeTest_GridImport([ValueSource("m_testFiles")]string File) {
            DerivativeTestMain.GRID_CASE = 50;
            DerivativeTestMain.GRID_FILE = Path.Combine("../../TestGrids/", File);
            DerivativeTestMain p = null;
            Quadrature_Bulksize.CHUNK_DATA_LIMIT = CHUNK_DATA_LIMIT_bkup; // might have been changed by other test, needs re-set
            if(CHUNK_DATA_LIMIT_bkup < 1)
                throw new ApplicationException();

            BoSSS.Solution.Application._Main(new string[0], true, null, delegate() {
                p = new DerivativeTestMain();
                return p;
            });

            Assert.IsTrue(p.m_passed);
        }


    }




    class DerivativeTestMain : BoSSS.Solution.Application {

        public static int GRID_CASE = 13;
        public static string GRID_FILE = "..\\..\\TestGrids\\wedding2D_v16.cgns";
        
        static void Main(string[] args) {
            //Quadrature_Bulksize.BULKSIZE_LIMIT_OVERRIDE = 1;
            BoSSS.Solution.Application.InitMPI(args);

            // Build-In Grids
            // ==============

            
            for(int i = 9; i <= 9; i++) {
                BoSSS.Solution.Application._Main(args, true, null, delegate() {
                    var R = new DerivativeTestMain();
                    GRID_CASE = i;
                    return R;
                });
            }
            //*/

            // gmsh Grids
            // ==========
            
            /*
            string[] gmshMeshFiles = Directory.GetFiles(@"../../TestGrids/", "WallMountedCube.cgns");
            //string[] gmshMeshFiles = Directory.GetFiles(@"../../TestGrids/", "ring.cgns");
            foreach(string gmf in gmshMeshFiles) {
                Console.Write(Path.GetFileName(gmf) + " ... ");

                BoSSS.Solution.Application._Main(args, true, null, delegate() {
                    var R = new DerivativeTestMain();
                    GRID_CASE = 50;
                    GRID_FILE = gmf;
                    return R;
                });
            }
            //*/

            Console.WriteLine("Number of cache hits:   " + BoSSS.Foundation.Caching.Cache.Hits);
            Console.WriteLine("Number of cache misses: " + BoSSS.Foundation.Caching.Cache.Misses);

            BoSSS.Solution.Application.FinalizeMPI();
        }

        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
        }

        SinglePhaseField f1;
        SinglePhaseField f2;
        SinglePhaseField[] f1Gradient_Analytical;
        SinglePhaseField[] f1Gradient_Numerical;
        SinglePhaseField[] f2Gradient_Analytical;
        SinglePhaseField[] f2Gradient_Numerical;
        SinglePhaseField Laplace_f1_Numerical;
        SinglePhaseField Laplace_f2_Numerical;
        SinglePhaseField Laplace_f1_Analytical;
        SinglePhaseField Laplace_f2_Analytical;

        protected override void CreateFields() {
            int GridDeg = this.Grid.Cells.Select(cl => this.Grid.GetRefElement(cl.Type).GetInterpolationDegree(cl.Type)).Max();
            int D = this.GridData.SpatialDimension;
            if(this.Grid.GetRefElement(this.Grid.Cells[0].Type) == Square.Instance
                || this.Grid.GetRefElement(this.Grid.Cells[0].Type) == Cube.Instance) {
                // hack: otherwise DG deg gets to high
                GridDeg = (int)Math.Round(Math.Pow(GridDeg, 1.0 / D));
            }

            Basis b = new Basis(this.GridData, 2 + GridDeg);
           

            Console.WriteLine("Grid degree is "+ GridDeg + " => Using DG order: " + b.Degree);
                       

            f1 = new SinglePhaseField(b, "f1");
            f2 = new SinglePhaseField(b, "f2");
            Laplace_f1_Numerical = new SinglePhaseField(b, "Laplace_f1_Numerical");
            Laplace_f2_Numerical = new SinglePhaseField(b, "Laplace_f2_Numerical");
            Laplace_f1_Analytical = new SinglePhaseField(b, "Laplace_f1_Analytical");
            Laplace_f2_Analytical = new SinglePhaseField(b, "Laplace_f2_Analytical");

            f1Gradient_Analytical = new SinglePhaseField[D];
            f2Gradient_Analytical = new SinglePhaseField[D];
            f1Gradient_Numerical = new SinglePhaseField[D];
            f2Gradient_Numerical = new SinglePhaseField[D];

            for (int d = 0; d < D; d++) {
                f1Gradient_Analytical[d] = new SinglePhaseField(b, string.Format("df1_dx{0}_Analytical", d));
                f2Gradient_Analytical[d] = new SinglePhaseField(b, string.Format("df2_dx{0}_Analytical", d));
                f1Gradient_Numerical[d] = new SinglePhaseField(b, string.Format("df1_dx{0}_Numerical", d));
                f2Gradient_Numerical[d] = new SinglePhaseField(b, string.Format("df2_dx{0}_Numerical", d));
            }
        }

        double EdgeArea = -1;
        double CellVolume = -1;

        protected override GridCommons CreateOrLoadGrid() {

            GridCommons grd;
            switch (GRID_CASE) {
                case 1:
                grd = Grid1D.LineGrid(GenericBlas.Linspace(-4, 4, 5));
                break;


                case 2: {
                    grd = Grid1D.LineGrid(GenericBlas.Linspace(-4, 4, 20));
                    break;
                }

                case 3: {
                    double[] xnodes = new double[] { -2, 0, 2 };
                    double[] ynodes = new double[] { -2, 0, 2 };
                    double dx = xnodes[1] - xnodes[0];
                    double dy = ynodes[1] - ynodes[0];
                    //this.CellVolume = dx * dy;
                    //if(Math.Abs(dx - dy) <= 1.0e-12)
                    //    EdgeArea = dx;
                    grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, periodicX: false, periodicY: false, type: CellType.Square_4);
                    break;
                }

                case 4: {
                    double[] xnodes = GenericBlas.Linspace(-1, 5, 9);
                    double[] ynodes = GenericBlas.Linspace(-1, 5, 13);
                    double dx = xnodes[1] - xnodes[0];
                    double dy = ynodes[1] - ynodes[0];
                    this.CellVolume = dx * dy;
                    if (Math.Abs(dx - dy) <= 1.0e-12)
                        EdgeArea = dx;
                    grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, periodicX: false, periodicY: false, type: CellType.Square_4);
                    break;
                }

                case 5: {
                    double[] xnodes = GenericBlas.Linspace(-1, 1, 8);
                    double[] ynodes = GenericBlas.Linspace(-1, 1, 13);
                    grd = Grid2D.UnstructuredTriangleGrid(xnodes, ynodes, JitterScale: 0.5);
                    break;
                }

                case 6: {
                    grd = Circle();
                    break;
                }

                case 7: {
                    // test periodicity

                    grd = Grid2D.CurvedSquareGrid(GenericBlas.Linspace(1, 2, 4), GenericBlas.Linspace(0, 0.25, 10), CellType.Square_9, PeriodicS: true);
                    AltRefSol = true;
                    break;
                }

                case 8: {
                    double[] rNodes = GenericBlas.Linspace(1, 4, 8);
                    double[] sNodes = GenericBlas.Linspace(0, 0.5, 15);
                    grd = Grid2D.CurvedSquareGrid(rNodes, sNodes, CellType.Square_4, PeriodicS: false);
                    break;
                }

                case 9: {
                    double[] xNodes1 = GenericBlas.Linspace(-1, 0.3, 7);
                    double[] yNodes1 = GenericBlas.Linspace(-1, 1, 13);
                    double[] xNodes2 = GenericBlas.Linspace(0.3, 1, 5);
                    double[] yNodes2 = GenericBlas.Linspace(-1, 1, 25);
                    double[] xNodes3 = GenericBlas.Linspace(-1, 1, 8);
                    double[] yNodes3 = GenericBlas.Linspace(-2, -1, 5);

                    var grd1 = Grid2D.Cartesian2DGrid(xNodes1, yNodes1, type: CellType.Square_Linear);
                    var grd2 = Grid2D.Cartesian2DGrid(xNodes2, yNodes2, type: CellType.Square_Linear);
                    var grd3 = Grid2D.Cartesian2DGrid(xNodes3, yNodes3, type: CellType.Square_Linear);
                    var grdJ = GridCommons.MergeLogically(grd1, GridCommons.MergeLogically(grd2, grd3));
                    grd = GridCommons.Seal(grdJ, 4);

                    break;
                }

                case 10: {
                    
                    double[] xNodes1 = GenericBlas.Linspace(-1, 0.3, 4);
                    double[] xNodes2 = GenericBlas.Linspace(0.3, 1, 5);

                    double[] yNodes1 = GenericBlas.Linspace(-1, 1, 9);
                    double[] yNodes2 = GenericBlas.Linspace(-1, 1, 5);

                    double[] zNodes1 = GenericBlas.Linspace(-1, 1, 5);
                    double[] zNodes2 = GenericBlas.Linspace(-1, 1, 3);

                    var grd1 = Grid3D.Cartesian3DGrid(xNodes1, yNodes1, zNodes1);
                    var grd2 = Grid3D.Cartesian3DGrid(xNodes2, yNodes2, zNodes2);
                    var grdJ = GridCommons.MergeLogically(grd1, grd2);
                    grd = GridCommons.Seal(grdJ, 4);
                                                  

                    break;

                }


                case 11: {
                    double[] rNodes = GenericBlas.Linspace(1, 4, 8);
                    double[] sNodes = GenericBlas.Linspace(0, 0.5, 15);
                    grd = Grid2D.CurvedSquareGrid(rNodes, sNodes, CellType.Square_9, PeriodicS: false);
                    break;
                }

                case 12: {
                    double[] rNodes = GenericBlas.Linspace(1, 4, 13);
                    double[] sNodes = GenericBlas.Linspace(0, 0.5, 25);
                    grd = Grid2D.CurvedSquareGrid(rNodes, sNodes, CellType.Square_16, PeriodicS: false);
                    break;
                }

                case 13: {
                    double[] xnodes = GenericBlas.Linspace(-1, 1, 7);
                    double[] ynodes = GenericBlas.Linspace(-1, 1, 9);
                    double[] znodes = GenericBlas.Linspace(-1, 1, 8);
                    grd = Grid3D.Cartesian3DGrid(xnodes, ynodes, znodes, periodicX: false, periodicY: false, periodicZ: false);
                    break;
                }

                case 14: {
                    double[] rNodes = GenericBlas.Linspace(1, 2, 4);
                    double[] sNodes = GenericBlas.Linspace(0, 0.5, 4);
                    double[] zNodes = GenericBlas.Linspace(-1, 1, 5);
                    grd = Grid3D.CylinderGrid(rNodes, sNodes, zNodes, CellType.Cube_27, PeriodicS: false, PeriodicZ: false);
                    break;
                }

                case 15: {
                    grd = Grid2D.Trapezoidal2dGrid(4, 2, 2, GenericBlas.Linspace(0, 1, 2));
                    break;
                }


                case 50: {
                    // gmsh grid import test

                    Console.WriteLine("Loading file: '" + GRID_FILE + "'...");
                    grd = GridImporter.Import(GRID_FILE);
                    //Console.WriteLine("done. " + grd.NoOfUpdateCells.MPISum() + " cells loaded.");

                    //Plot2dGridGnuplot(grd);

                    HashSet<CellType> cellTypes = new HashSet<CellType>();
                    foreach (var cell in grd.Cells) {
                        if (!cellTypes.Contains(cell.Type))
                            cellTypes.Add(cell.Type);
                    }
                    Console.Write("Cell types: ");
                    foreach (var ct in cellTypes) {
                        Console.Write(ct);
                        Console.Write(" ");
                    }
                    Console.WriteLine();


                    break;
                }

                default:
                throw new NotSupportedException();

            }
            return grd;
        }

        static GridCommons Circle() {
            int R = 7;
            bool aber = true;

            CellType ct = CellType.Triangle_6;
            NodeSet triNodes = Triangle.Instance.GetInterpolationNodes(ct);
            //MultidimensionalArray triCornerNodes = triNodes.ExtractSubArrayShallow(new int[] {3,0}, new int[] { 5, 1});

            GridCommons grd = new Grid2D(Triangle.Instance);
            grd.Cells = new Cell[R];

            int off = ct == CellType.Triangle_6 ? 3 : 0;

            for(int r = 0; r < R; r++) {
                double dAlpha = 2*Math.PI / ((double)(R));
                double alpha = dAlpha * r;
                double beta = alpha + dAlpha;
                double gamma = alpha + dAlpha * 0.5;

                MultidimensionalArray nodesGlobal = MultidimensionalArray.Create(3 + off, 2);

                nodesGlobal[off + 0, 0] = 0;
                nodesGlobal[off + 0, 1] = 0;
                nodesGlobal[off + 1, 0] = Math.Cos(alpha);
                nodesGlobal[off + 1, 1] = Math.Sin(alpha);
                nodesGlobal[off + 2, 0] = Math.Cos(beta);
                nodesGlobal[off + 2, 1] = Math.Sin(beta);

                if(ct == CellType.Triangle_6) {
                    nodesGlobal[0, 0] = nodesGlobal[4, 0] * (aber ? 0.3 : 0.5);
                    nodesGlobal[0, 1] = nodesGlobal[4, 1] * (aber ? 0.3 : 0.5);
                    nodesGlobal[1, 0] = (aber ? Math.Cos(gamma) : (nodesGlobal[4, 0] + nodesGlobal[5, 0]) * 0.5);
                    nodesGlobal[1, 1] = (aber ? Math.Sin(gamma) : (nodesGlobal[4, 1] + nodesGlobal[5, 1]) * 0.5);
                    nodesGlobal[2, 0] = nodesGlobal[5, 0] * (aber ? 0.3 : 0.5);
                    nodesGlobal[2, 1] = nodesGlobal[5, 1] * (aber ? 0.3 : 0.5);
                }

                Cell cl = new Cell();
                cl.GlobalID = r;
                cl.Type = ct;
                cl.TransformationParams = nodesGlobal.CloneAs();
                cl.TransformationParams.LockForever();
                cl.NodeIndices = new int[] { R, r, (r + 1) % R };

                grd.Cells[r] = cl;
            }

            return grd;
        }

        bool AltRefSol = false;

        protected override void SetInitial() {

            if (this.GridData.SpatialDimension == 3) {

                f1.ProjectField((x, y, z) => (3 * x + z));
                f1Gradient_Analytical[0].ProjectField((x, y, z) => 3.0);
                f1Gradient_Analytical[1].ProjectField((x, y, z) => 0.0);
                f1Gradient_Analytical[2].ProjectField((x, y, z) => 1.0);

                f2.ProjectField((x, y, z) => z + 2 * y);
                f2Gradient_Analytical[0].ProjectField((x, y, z) => 0.0);
                f2Gradient_Analytical[1].ProjectField((x, y, z) => 2.0);
                f2Gradient_Analytical[2].ProjectField((x, y, z) => 1.0);

                Laplace_f1_Analytical.ProjectField((x, y, z) => 0.0);
                Laplace_f2_Analytical.ProjectField((x, y, z) => 0.0);

            } else if (this.GridData.SpatialDimension == 2) {
                if (AltRefSol == false) {
                    f1.ProjectField((x, y) => (3 * x));
                    f1Gradient_Analytical[0].ProjectField((x, y) => 3);
                    f1Gradient_Analytical[1].ProjectField((x, y) => 0);

                    f2.ProjectField((x, y) => x + 2 * y);
                    f2Gradient_Analytical[0].ProjectField((x, y) => 1.0);
                    f2Gradient_Analytical[1].ProjectField((x, y) => 2.0);

                    Laplace_f1_Analytical.ProjectField((x, y) => 0.0);
                    Laplace_f2_Analytical.ProjectField((x, y) => 0.0);
                } else {
                    f1.ProjectField((x, y) => Math.Sin(Math.Atan(y / x) * 4.0));
                    f1Gradient_Analytical[0].ProjectField((x, y) => (-4 * Math.Cos(4 * Math.Atan(y / x)) * y / (x * x) / (1 + (y * y) / (x * x))));
                    f1Gradient_Analytical[1].ProjectField((x, y) => (4 * Math.Cos(4 * Math.Atan(y / x)) / x / (1 + (y * y) / (x * x))));
                }
            } else if (this.GridData.SpatialDimension == 1) {
                f1.ProjectField((x) => (3 * x));
                f1Gradient_Analytical[0].ProjectField((_1D)((x) => 3));

                f2.ProjectField((x) => x * x);
                f2Gradient_Analytical[0].ProjectField((_1D)((x) => 2 * x));

                Laplace_f1_Analytical.ProjectField((_1D)((x) => 0.0));
                Laplace_f2_Analytical.ProjectField((_1D)((x) => 2.0));

            } else
                throw new NotImplementedException();
        }
        
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(
                ArrayTools.Cat<DGField>(f1Gradient_Analytical, f1Gradient_Numerical, f1, GridData.BoundaryMark(), Laplace_f1_Numerical, Laplace_f2_Numerical),
                "derivatives", 0.0, superSampling);
        }

        internal bool m_passed = true;
        
        void DerivativeByFluxLinear(SinglePhaseField fin, SinglePhaseField fres, int d, SinglePhaseField fBnd) {
            var Op = (new LinearDerivFlx(d)).Operator();

            MsrMatrix OpMtx = new MsrMatrix(fres.Mapping, fin.Mapping);
            double[] OpAff = new double[fres.Mapping.LocalLength];
            Op.ComputeMatrixEx(fin.Mapping, new DGField[] { fBnd }, fres.Mapping,
                OpMtx, OpAff, OnlyAffine: false);
           
            fres.Clear();
            fres.CoordinateVector.Acc(1.0, OpAff);
            OpMtx.SpMVpara(1.0, fin.CoordinateVector, 1.0, fres.CoordinateVector);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            base.EndTime = 0.0;
            base.NoOfTimesteps = 0;

            int D = this.GridData.SpatialDimension;
            int J = this.GridData.Cells.NoOfLocalUpdatedCells;

            Console.WriteLine("DerivativeTest.exe, test case #" + GRID_CASE + " ******************************");

            //var Fix = this.GridData.iGeomEdges.FaceIndices;
            //for(int iEdge = 0; iEdge < Fix.GetLength(0); iEdge++) {
            //    Debug.Assert(Fix[iEdge, 0] >= 0);
            //    Debug.Assert(Fix[iEdge, 1] >= 0);
            //}

            // sealing test
            // =================

            TestSealing(this.GridData);


            // cell volume and edge area check, if possible
            // ===============================================


            if (this.CellVolume > 0) {
                double err = 0;
                double Treshold = 1.0e-10;


                for(int j = 0; j < J; j++) {
                    err += Math.Abs(this.GridData.Cells.GetCellVolume(j) - this.CellVolume);
                }

                bool passed = (err < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine("Cell volume error: " + err + " passed? " + passed);
                Console.WriteLine("--------------------------------------------");

            }

            if(this.EdgeArea > 0) {
                double err = 0;
                double Treshold = 1.0e-10;

                int E = this.GridData.Edges.Count;

                for(int e = 0; e < E; e++) {
                    err += Math.Abs(this.GridData.Edges.GetEdgeArea(e) - this.EdgeArea);
                }

                bool passed = (err < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine("Edge area error: " + err + " passed? " + passed);
                Console.WriteLine("--------------------------------------------");

            }

            // Orthonormality of basis in physical coords
            // ==========================================

            {
                Basis Bs = this.f1.Basis;
                int N = Bs.Length;
                int degQuad = this.GridData.Cells.GetInterpolationDegree(0) * D + Bs.Degree + 3;

                // mass matrix: should be identity!
                MultidimensionalArray MassMatrix = MultidimensionalArray.Create(J, N, N);

                // compute mass matrix by quadrature.
                var quad = CellQuadrature.GetQuadrature(new int[] { N, N }, base.GridData,
                    (new CellQuadratureScheme()).Compile(base.GridData, degQuad),
                    delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        NodeSet QuadNodes = QR.Nodes;
                        MultidimensionalArray BasisVals = Bs.CellEval(QuadNodes, i0, Length);
                        EvalResult.Multiply(1.0, BasisVals, BasisVals, 0.0, "jknm", "jkn", "jkm");
                    },
                    delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        MassMatrix.SetSubArray(ResultsOfIntegration, new int[] { i0, 0, 0 }, new int[] { i0 + Length - 1, N - 1, N - 1 });
                    },
                    cs: CoordinateSystem.Physical);
                quad.Execute();

                // check that mass matrix is Id.
                int MaxErrorCell = -1;
                double MaxError = -1;
                for(int j = 0; j < J; j++) {
                    MultidimensionalArray MassMatrix_j = MassMatrix.ExtractSubArrayShallow(j, -1, -1);
                    MassMatrix_j.AccEye(-1.0);

                    double Norm_j = MassMatrix_j.InfNorm();
                    if(Norm_j > MaxError) {
                        MaxError = Norm_j;
                        MaxErrorCell = j;
                    }

                }

                bool passed = (MaxError < 1.0e-8);
                m_passed = m_passed && passed;
                Console.WriteLine("Mass Matrix, maximum error in Cell #" + MaxErrorCell + ", mass matrix error norm: " + MaxError + " passed? " + passed);
            }
            
            // Broken Derivatives
            // =================

            double totalVolume = (new SubGrid(CellMask.GetFullMask(this.GridData))).Volume;
            
            for (int d = 0; d < D; d++) {

                // compute
                f1Gradient_Numerical[d].Clear();
                f1Gradient_Numerical[d].Derivative(1.0, f1, d);
                f2Gradient_Numerical[d].Clear();
                f2Gradient_Numerical[d].Derivative(1.0, f2, d);

                // subtract analytical
                var Errfield1 = f1Gradient_Numerical[d].CloneAs();
                Errfield1.Acc(-1, f1Gradient_Analytical[d]);

                var Errfield2 = f2Gradient_Numerical[d].CloneAs();
                Errfield2.Acc(-1, f2Gradient_Analytical[d]);

                Console.WriteLine("Broken Derivatives: ");

                double Treshold = 1.0e-10;
                if(AltRefSol)
                    Treshold = 1.0e-4; // not exactly polynomial, therefore a higher threshold

                double err1_dx = Errfield1.L2Norm() / totalVolume;
                bool passed = (err1_dx < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| df1/dx{0}_Numerical - df1/dx{0}_Analytical ||_2 = {1}, passed? {2}", d, err1_dx, passed));

                double err2_dx = Errfield2.L2Norm() / totalVolume;
                passed = (err2_dx < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| df2/dx{0}_Numerical - df2/dx{0}_Analytical ||_2 = {1}, passed? {2}", d, err2_dx, passed));

                Console.WriteLine("--------------------------------------------");
            }
            

            // Flux Derivatives
            // =================
            for(int d = 0; d < D; d++) {
                // compute
                f1Gradient_Numerical[d].Clear();
                f1Gradient_Numerical[d].DerivativeByFlux(1.0, f1, d);
                f2Gradient_Numerical[d].Clear();
                f2Gradient_Numerical[d].DerivativeByFlux(1.0, f2, d);

                f1Gradient_Numerical[d].CheckForNanOrInf(true, true, true);
                f2Gradient_Numerical[d].CheckForNanOrInf(true, true, true);

                // subtract analytical
                var Errfield1 = f1Gradient_Numerical[d].CloneAs();
                Errfield1.Acc(-1, f1Gradient_Analytical[d]);

                var Errfield2 = f2Gradient_Numerical[d].CloneAs();
                Errfield2.Acc(-1, f2Gradient_Analytical[d]);

                Console.WriteLine("Flux Derivatives: ");

                double Treshold = 1.0e-10;
                if(AltRefSol)
                    Treshold = 1.0e-4; // not exactly polynomial, therefore a higher threshold

                double err1_dx = Errfield1.L2Norm() / totalVolume;
                bool passed = (err1_dx < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| df1/dx{0}_Numerical - df1/dx{0}_Analytical ||_2 = {1}, passed? {2}", d, err1_dx, passed));

                double err2_dx = Errfield2.L2Norm() / totalVolume;
                passed = (err2_dx < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| df2/dx{0}_Numerical - df2/dx{0}_Analytical ||_2 = {1}, passed? {2}", d, err2_dx, passed));

                Console.WriteLine("--------------------------------------------");
            }
            
            
            // Linear flux Derivatives
            // =======================
            for(int d = 0; d < D; d++) {

                double[] korrekto = f1Gradient_Numerical[d].CoordinateVector.ToArray();

                // compute
                DerivativeByFluxLinear(f1, f1Gradient_Numerical[d], d, f1);
                DerivativeByFluxLinear(f2, f2Gradient_Numerical[d], d, f2);

                // subtract analytical
                var Errfield1 = f1Gradient_Numerical[d].CloneAs();
                Errfield1.Acc(-1, f1Gradient_Analytical[d]);

                var Errfield2 = f2Gradient_Numerical[d].CloneAs();
                Errfield2.Acc(-1, f2Gradient_Analytical[d]);

                Console.WriteLine("Linear Flux Derivatives: ");

                double Treshold = 1.0e-10;
                if(AltRefSol)
                    Treshold = 1.0e-4; // not exactly polynomial, therefore a higher threshold

                double err1_dx = Errfield1.L2Norm() / totalVolume;
                bool passed = (err1_dx < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| df1/dx{0}_Numerical - df1/dx{0}_Analytical ||_2 = {1}, passed? {2}", d, err1_dx, passed));

                double err2_dx = Errfield2.L2Norm() / totalVolume;
                passed = (err2_dx < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| df2/dx{0}_Numerical - df2/dx{0}_Analytical ||_2 = {1}, passed? {2}", d, err2_dx, passed));

                Console.WriteLine("--------------------------------------------");
            }
            
            // Laplacian, nonlinear
            // ====================

            if(!AltRefSol) {
                var Laplace = (new ipLaplace()).Operator(1);

                Laplace.Evaluate(new DGField[] { this.f1 }, new DGField[] { this.Laplace_f1_Numerical });
                Laplace.Evaluate(new DGField[] { this.f2 }, new DGField[] { this.Laplace_f2_Numerical });

                double Treshold = 1.0e-8;

                // subtract analytical
                var Errfield1 = Laplace_f1_Numerical.CloneAs();
                Errfield1.Acc(-1, Laplace_f1_Analytical);

                var Errfield2 = Laplace_f2_Numerical.CloneAs();
                Errfield2.Acc(-1, Laplace_f2_Analytical);

                double err_Lf1 = Errfield1.L2Norm() / totalVolume;
                bool passed = (err_Lf1 < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| /\\f1 Numerical - /\\f1 Analytical ||_2 = {0} (nonlinear evaluation), passed? {1}", err_Lf1, passed));

                double err_Lf2 = Errfield2.L2Norm() / totalVolume;
                passed = (err_Lf2 < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| /\\f2 Numerical - /\\f2 Analytical ||_2 = {0} (nonlinear evaluation), passed? {1}", err_Lf2, passed));

                Console.WriteLine("--------------------------------------------");
            }


            // Laplacian, linear
            // ====================

            if(!AltRefSol) {
                var Laplace = (new ipLaplace()).Operator(1);

                var LaplaceMtx = new MsrMatrix(this.f1.Mapping, this.Laplace_f1_Numerical.Mapping);
                var LaplaceAffine = new double[LaplaceMtx.RowPartitioning.LocalLength];
                
                Laplace.ComputeMatrix(this.f1.Mapping, null, this.Laplace_f1_Numerical.Mapping,
                    LaplaceMtx, LaplaceAffine, false);

                this.Laplace_f1_Numerical.CoordinateVector.SetV(LaplaceAffine);
                LaplaceMtx.SpMVpara(1.0, this.f1.CoordinateVector, 1.0, this.Laplace_f1_Numerical.CoordinateVector);

                this.Laplace_f2_Numerical.CoordinateVector.SetV(LaplaceAffine);
                LaplaceMtx.SpMVpara(1.0, this.f2.CoordinateVector, 1.0, this.Laplace_f2_Numerical.CoordinateVector);

                // subtract analytical
                var Errfield1 = Laplace_f1_Numerical.CloneAs();
                Errfield1.Acc(-1, Laplace_f1_Analytical);

                var Errfield2 = Laplace_f2_Numerical.CloneAs();
                Errfield2.Acc(-1, Laplace_f2_Analytical);


                double Treshold = 1.0e-8;

                double err_Lf1 = Errfield1.L2Norm() / totalVolume;
                bool passed = (err_Lf1 < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| /\\f1 Numerical - /\\f1 Analytical ||_2 = {0} (linear evaluation), passed? {1}", err_Lf1, passed));

                double err_Lf2 = Errfield2.L2Norm() / totalVolume;
                passed = (err_Lf2 < Treshold);
                m_passed = m_passed && passed;
                Console.WriteLine(string.Format("|| /\\f2 Numerical - /\\f2 Analytical ||_2 = {0} (linear evaluation), passed? {1}", err_Lf2, passed));

                               
                Console.WriteLine("--------------------------------------------");
            }

    
            // finally...
            // =================

            if(m_passed)
                Console.WriteLine("All tests passed. *****************************");
            else
                Console.WriteLine("Some error above threshold. *******************");



            return 0.0; // return some artificial timestep
        }





        static void Plot2dGridGnuplot(GridCommons grd) {

            using(var gp = new Gnuplot()) {
                int J = grd.Cells.Length;


                double xmin = double.MaxValue, xmax = double.MinValue, ymin = double.MaxValue, ymax = double.MinValue;
                for(int j = 0; j < J; j++) {

                    var Cell = grd.Cells[j];
                    var Kref = grd.GetRefElement(Cell.Type);

                    int[] R;
                    if(Kref.GetType() == typeof(Triangle)) {
                        R = new int[] { 0, 1, 2, 0 };

                    } else if(Kref.GetType() == typeof(Square)) {
                        R = new int[] { 0, 1, 3, 2, 0 };


                    } else {
                        throw new NotSupportedException();
                    }

                    int I = 4;
                    int K = 20;
                    NodeSet LocNodes = new NodeSet(Kref, I* K * (R.Length - 1), 2);
                    var vtx = Kref.Vertices;
                    double alpha = 1.0 / (K - 1);
                    for(int iFace = 0; iFace < R.Length - 1; iFace++) {
                        for(int k = 0; k < K; k++) {
                            double a = alpha * k;

                            for(int d = 0; d < 2; d++)
                                LocNodes[iFace * K + k, d] = vtx[R[iFace], d] * (1 - a) + vtx[R[iFace + 1], d] * a;
                        }
                    }

                    for(int i = 0; i < I; i++) {
                        int ind0 = K * (R.Length - 1) * i;
                        int indE = K * (R.Length - 1) * (i + 1) - 1;
                        int indp0 = K * (R.Length - 1) * (i - 1);
                        int indpE = K * (R.Length - 1) * i - 1;

                        var LocNodes_i = LocNodes.ExtractSubArrayShallow(new int[] { ind0, 0 }, new int[] { indE, 1 });

                        if(i > 0) {
                            var LocNodes_iP = LocNodes.ExtractSubArrayShallow(new int[] { indp0, 0 }, new int[] { indpE, 1 });
                            LocNodes_i.Set(LocNodes_iP);
                            LocNodes_i.Scale(0.65);
                        } else {
                            LocNodes_i.Scale(0.9);
                        }
                    }


                    LocNodes.LockForever();

                    MultidimensionalArray GlobalNodes = Transform(Kref, Cell, LocNodes);

                    xmin = Math.Min(xmin, GlobalNodes.ExtractSubArrayShallow(-1, 0).Min());
                    xmax = Math.Max(xmax, GlobalNodes.ExtractSubArrayShallow(-1, 0).Max());
                    ymin = Math.Min(ymin, GlobalNodes.ExtractSubArrayShallow(-1, 1).Min());
                    ymax = Math.Max(ymax, GlobalNodes.ExtractSubArrayShallow(-1, 1).Max());

                    gp.PlotXY(GlobalNodes.GetColumn(0), GlobalNodes.GetColumn(1), title:("tri" + j), 
                        format: new PlotFormat( lineColor: ((LineColors)j)));
                }

                gp.SetXRange(xmin - 0.1, xmax + 0.1);
                gp.SetYRange(ymin - 0.1, ymax + 0.1);

                gp.Execute();


                Console.WriteLine("press any key to continue...");
                Console.ReadKey();
            }
        }

        static MultidimensionalArray Transform(RefElement Kref, Cell Cl, NodeSet Nodes) {

            int D = Kref.SpatialDimension;
            PolynomialList polys = Kref.GetInterpolationPolynomials(Cl.Type);
            MultidimensionalArray polyVals = polys.Values.GetValues(Nodes);
            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.Create(Nodes.NoOfNodes, D);

            for(int d = 0; d < D; d++) {
                GlobalVerticesOut.ExtractSubArrayShallow(-1, d)
                    .Multiply(1.0, polyVals, Cl.TransformationParams.ExtractSubArrayShallow(-1, d), 0.0, "k", "kn", "n");
            }

            return GlobalVerticesOut;
        }


        /// <summary>
        /// Compares the cell surface and the boundary integral.
        /// </summary>
        public static void TestSealing(GridData gdat) {
            int J = gdat.Cells.NoOfLocalUpdatedCells;
            int[,] E2C = gdat.Edges.CellIndices;

            //
            // compute cell surface via edge integrals
            //
            MultidimensionalArray CellSurf1 = MultidimensionalArray.Create(J);
            EdgeQuadrature.GetQuadrature(new int[] { 1 },
                gdat, new EdgeQuadratureScheme().Compile(gdat, 2),
                delegate (int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) { // evaluate
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // save results
                    for (int i = 0; i < Length; i++) {
                        int iEdge = i + i0;
                        int jCellIn = E2C[iEdge, 0];
                        int jCellOt = E2C[iEdge, 1];

                        CellSurf1[jCellIn] += ResultsOfIntegration[i, 0];
                        if (jCellOt >= 0 && jCellOt < J)
                            CellSurf1[jCellOt] += ResultsOfIntegration[i, 0];
                    }
                }).Execute();

            //
            // compute cell surface via cell boundary integrals
            //
            var cbqs = new CellBoundaryQuadratureScheme(false, null);
            foreach (var Kref in gdat.Cells.RefElements)
                cbqs.AddFactory(new StandardCellBoundaryQuadRuleFactory(Kref));

            MultidimensionalArray CellSurf2 = MultidimensionalArray.Create(J);
            CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(new int[] { 1 },
                gdat, cbqs.Compile(gdat, 2),
                delegate (int i0, int Length, CellBoundaryQuadRule rule, MultidimensionalArray EvalResult) { // evaluate
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // save results
                    int NoOfFaces = ResultsOfIntegration.GetLength(1);
                    for (int i = 0; i < Length; i++) {
                        int jCell = i + i0;
                        for (int iFace = 0; iFace < NoOfFaces; iFace++) {
                            CellSurf2[jCell] += ResultsOfIntegration[i, iFace, 0];
                        }
                    }
                }, cs: CoordinateSystem.Physical).Execute();

            //
            // compare
            //
            MultidimensionalArray Err = CellSurf1.CloneAs();
            Err.Acc(-1.0, CellSurf2);
            double ErrNorm = Err.L2Norm();
            double TotSurf = CellSurf1.L2Norm();
            Console.WriteLine("Area Check " + Err.L2Norm());
            Assert.LessOrEqual(ErrNorm / TotSurf, 1.0e-8);

            //for (int j = 0; j < J; j++) {
            //    if (Err[j].Abs() >= 1.0e-6) {
            //        Console.WriteLine("Mismatch between edge area and cell surface area in cell {0}, GlobalId {1}, No. of neighbors {3},  {2:0.####E-00}", j, gdat.CurrentGlobalIdPermutation.Values[j], Err[j], gdat.Cells.CellNeighbours[j].Length);
            //        schas.SetMeanValue(j, 1);
            //    }
            //}


         }



    }
}
