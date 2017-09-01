using System;
using System.Collections.Generic;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using ilPSP.LinSolvers;
using NUnit.Framework;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Quadrature;
using ilPSP.Utils;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Gnuplot;
using System.Linq;
using System.Diagnostics;
using BoSSS.Foundation.Comm;
using BoSSS.Solution.GridImport;

namespace BoSSS.Application.HangingNodeTest {


    class HangingNodeTestMain : BoSSS.Solution.Application {
        static void Main(string[] args) {
            BoSSS.Solution.Application._Main(args, true, null, delegate() {
				Console.WriteLine("Pointer size: " + IntPtr.Size);
                return new HangingNodeTestMain();
            });
        }

        protected override void CreateEquationsAndSolvers() {
        }

        Field f1;
        Field f2;
        Field df1_dx_Analytical;
        Field df1_dy_Analytical;
        Field df1_dx_Numerical;
        Field df1_dy_Numerical;

        Field df2_dx_Analytical;
        Field df2_dy_Analytical;
        Field df2_dx_Numerical;
        Field df2_dy_Numerical;



        protected override void CreateFields() {
            Basis b = new Basis(this.GridDat, 1);
            f1 = new SinglePhaseField(b, "f1");
            f2 = new SinglePhaseField(b, "f2");
            df1_dx_Analytical = new SinglePhaseField(b, "df1_dx_Analytical");
            df1_dy_Analytical = new SinglePhaseField(b, "df1_dy_Analytical");
            df1_dx_Numerical = new SinglePhaseField(b, "df1_dx_Numerical");
            df1_dy_Numerical = new SinglePhaseField(b, "df1_dy_Numerical");

            df2_dx_Analytical = new SinglePhaseField(b, "df2_dx_Analytical");
            df2_dy_Analytical = new SinglePhaseField(b, "df2_dy_Analytical");
            df2_dx_Numerical = new SinglePhaseField(b, "df2_dx_Numerical");
            df2_dy_Numerical = new SinglePhaseField(b, "df2_dy_Numerical");
        }

        //protected override IIODriver GetIODriver() {
        //    return new StandartFsDriver(new string[] { "c:\\BoSSS\\bosss_db" });
        //}




        protected override GridCommons CreateOrLoadGrid() {
            //var grd = Grid2D.Cartesian2DGrid(Grid1D.Linspace(-1, 3, 3), Grid1D.Linspace(-1, 1, 3), type:CellType.Square_9, periodicX: false, periodicY: false);

            //var Center = MultidimensionalArray.Create (1, 2);
            //var Glob = MultidimensionalArray.Create (1, 1, 2);
            //grd.RefElements [0].TransformLocal2Global (Center, Glob, 0, grd.Cells [0].Type, grd.Cells [0].TransformationParams);
            
            //grd.EdgeTagsNames.Add(1, "Wall");
            //grd.DefineEdgeTags(delegate(double[] X) {
            //    return (byte)1;
            //});
            //var grd = Grid2D.CurvedSquareGrid(Grid1D.Linspace(1, 2, 2), Grid1D.Linspace(0, 1.0/2.0, 4), CellType.Square_Linear, PeriodicS:true);

            //var grd = Grid2D.CurvedSquareGrid(Grid1D.Linspace(1.6, 2, 3), Grid1D.Linspace(0, 1, 9), CellType.Square_9);


            Gmesh import = new Gmesh("untitled.msh");
            var grd = import.GenerateBoSSSGrid();

            int J = grd.Cells.Length;
            int cnt = 0;
            for (int j = 0; j < J; j++) {
                int[] NodeIdx = new int[4];
                //for (int k = 0; k < NodeIdx.Length; k++) {
                //    NodeIdx[k] = cnt;
                //    cnt++;
                //}
                //grd.Cells[j].NodeIndices = NodeIdx;
                grd.Cells[j].CellFaceTags = null;
            }

            PlotGrid(grd);

            return grd; 
        }

        protected override void SetInitial() {

            if (this.GridDat.SpatialDimension == 3) {

                f1.ProjectField((x, y, z) => (3*y));
                df1_dx_Analytical.ProjectField((x, y, z) => 0);
                df1_dy_Analytical.ProjectField((x, y, z) => 3);

                // ----------------------------------------------------------------------------------------------------

                f2.ProjectField((x, y, z) => x + 2 * y);
                df2_dx_Analytical.ProjectField((x, y, z) => 1.0);
                df2_dy_Analytical.ProjectField((x, y, z) => 2.0);

            } else if (this.GridDat.SpatialDimension == 2) {
                f1.ProjectField((x, y) => (3*y));
                df1_dx_Analytical.ProjectField((x, y) => 0);
                df1_dy_Analytical.ProjectField((x, y) => 3);

                // ----------------------------------------------------------------------------------------------------

                f2.ProjectField((x, y) => x + 2 * y);
                df2_dx_Analytical.ProjectField((x, y) => 1.0);
                df2_dy_Analytical.ProjectField((x, y) => 2.0);

            } else
                throw new NotImplementedException();
        }



        protected override void PlotCurrentState(double physTime, int timestepNo) {
            //var fields = new Field[] { f1, df1_dx_Analytical, df1_dx_Numerical, df1_dy_Analytical, df1_dy_Numerical };
            var fields = new Field[] { f1, df1_dx_Numerical, df1_dy_Numerical };
            Tecplot.PlotFields(fields, GridDat, "derivatives", "derivatives", 0.0, 4);
        }

        internal bool passed = true;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            base.EndTime = 0.0;
            base.NoOfTimesteps = 0;

            // Flux Derivatives
            // =================
            
            {

                // compute 
                df1_dx_Numerical.Clear();
                df1_dx_Numerical.DerivativeByFlux(1.0, f1, 0);
                df1_dy_Numerical.DerivativeByFlux(1.0, f1, 1);


                Console.WriteLine("Local Derivatives: ");

                // print diagnostic:
                df1_dx_Analytical.Acc(-1.0, df1_dx_Numerical);
                df1_dy_Analytical.Acc(-1.0, df1_dy_Numerical);

                double err_dx = df1_dx_Analytical.L2Norm();
                Console.WriteLine("|| df/dx_Numerical - df/dx_Analytical ||_2 = " + err_dx);
                passed = passed && (err_dx < 1.0e-10);

                double err_dy = df1_dy_Analytical.L2Norm();
                Console.WriteLine("|| df/dy_Numerical - df/dy_Analytical ||_2 = " + err_dy);
                passed = passed && (err_dy < 1.0e-10);

                Console.WriteLine("--------------------------------------------");
            }
            





            //PlotCurrentState(0, 0);


            //OrthonormalityTest();

            //PlotCell0();

            return 0.0; // return some artificial timestep
        }


        void PlotCell0() {
            BoundingBox BB = new BoundingBox(2);
            GridDat.Cells.GetCellBoundingBox(0, BB);

            int Nx = 50;
            int Ny = 55;
            var xNodes = Grid1D.Linspace(BB.Min[0], BB.Max[0], Nx);
            var yNodes = Grid1D.Linspace(BB.Min[1], BB.Max[1], Ny);

            var GlobalNodes = MultidimensionalArray.Create(Nx, Ny, 2);
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    GlobalNodes[i, j, 0] = xNodes[i];
                    GlobalNodes[i, j, 1] = yNodes[j];
                }
            }

            var GlobalNodes_C = GlobalNodes.ResizeShallow(Nx*Ny, 2);
            var LocalNodes = MultidimensionalArray.Create(1,Nx*Ny,2);
            GridDat.TransformGlobal2Local(GlobalNodes_C, LocalNodes, 0, 1, 0);

            var Values = MultidimensionalArray.Create(1, Nx, Ny);
            var Values_ = Values.ResizeShallow(1, Nx*Ny);
            
            var Gradients = MultidimensionalArray.Create(1, Nx, Ny, 2);
            var Gradients_ = Gradients.ResizeShallow(1, Nx*Ny, 2);


            var NodeSet = GridDat.NSC.CreateContainer(LocalNodes.ExtractSubArrayShallow(0, -1, -1), 0);
            var lh = GridDat.NSC.LockNodeSetFamily(NodeSet);
            f1.Evaluate(0, 1, 0, Values_);
            f1.EvaluateGradient(0, 1, 0, Gradients_);
            GridDat.NSC.UnlockNodeSetFamily(lh);


            xNodes.SaveToTextFile("C:\\tmp\\xNodes.txt");
            yNodes.SaveToTextFile("C:\\tmp\\yNodes.txt");

            FullMatrix output1 = new FullMatrix(Values.ExtractSubArrayShallow(0, -1, -1));
            output1.ToTxtFile("C:\\tmp\\f.txt");

            FullMatrix output2 = new FullMatrix(Gradients.ExtractSubArrayShallow(0, -1, -1, 0));
            output2.ToTxtFile("C:\\tmp\\df.txt");

            
        }



        void OrthonormalityTest() {
            Basis B = f1.Basis;
            int N = B.Length;
            var Mtx = new FullMatrix(N, N);

            double TotErrSum = 0.0;

            CellQuadrature.GetQuadrature(new int[] { N, N }, this.GridDat,
                (new CellQuadratureScheme(true)).Compile(GridDat, Math.Min(B.Degree*2, 16)),
                delegate(MultidimensionalArray NodesUntransformed, int iKref) {
                    var ret = new NodeSetController.NodeSetContainer[] {
                        GridDat.NSC.CreateContainer(NodesUntransformed, iKref)
                    };
                    return ret;
                },
                delegate(int i0, int Length, int NoOfNodes, MultidimensionalArray EvalResult) { // void Del_Evaluate
                    var BasisVal = B.CellEval(0, i0, Length);

                    EvalResult.Multiply(1.0, BasisVal, BasisVal, 0.0, "jknm", "jkn", "jkm");
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        var MassMtx = ResultsOfIntegration.ExtractSubArrayShallow(i, -1, -1);

                        double errsum = 0;
                        for (int n = 0; n < N; n++) {
                            for (int m = 0; m < N; m++) {
                                double soll = (n == m) ? 1.0 : 0.0;
                                errsum += Math.Abs(MassMtx[m, n] - soll);
                            }
                        }

                        TotErrSum += errsum;
                    }
                }).Execute();


            Console.WriteLine("orthonormality error sum:" + TotErrSum);

        }




    }
}
