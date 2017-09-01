//using System;
//using System.Collections.Generic;
//using BoSSS.Foundation;
//using BoSSS.Foundation.Grid;
//using BoSSS.Platform;
//using BoSSS.Foundation.IO;
//using ilPSP.LinSolvers;
//using MPI.Wrappers;
//using System.Diagnostics;
//using System.IO;
//using log4net;
//using BoSSS.Solution.Tecplot;
//using BoSSS.Solution.Utils;
//using BoSSS.Platform.LinAlg;

//namespace CurvedElementsTest {

//    class MyGrid : GridCommons {

//        public MyGrid() {

//            // define reference element
//            m_GridSimplices = new RefElement[1];
//            m_GridSimplices[0] = new Square();    // simplex 0 is a square
//            //m_GridSimplices[0] = new Triangle();

//            // define cells
//            base.Cells = new Cell[2];


//            base.Cells[0] = new Cell();
//            base.Cells[0].GlobalID = 0;
//            base.Cells[0].MajorCellTypeIndex = 0; // Simplex is a square
//            base.Cells[0].MinorCellTypeIndex = 0; // Simplex is linear
//            base.Cells[0].TransformationParams = new double[2, 4];
//            //base.Cells[0].TransformationParams = new double[2, 3];
//            base.Cells[0].TransformationParams[0, 0] = 0.0;


//        }



//    }


//    class CurvedElementsTest_Main {

//        static ILog Logger = LogManager.GetLogger(typeof(CurvedElementsTest_Main));


//        public static void Main(string[] args) {
//            bool dummy;
//            ilPSP.Enviroment.Bootstrap(args, BoSSS.Solution.Application.GetBoSSSInstallDir(), out dummy);

//            //var grd = new Grid1D(Grid1D.Linspace(-5,5,11)); // 1D grid
//            ///var grd = new Cartesian2DGrid(Grid1D.Linspace(-5, 5, 3), Grid1D.Linspace(-2, 2, 2), null, null, false, false);
//            //var grd = new UnstructuredTriangleGrid(Grid1D.Linspace(-5, 5, 2), Grid1D.Linspace(-2, 2, 3));
//            //var grd = new BilinearSquareGrid(Grid1D.Linspace(-4, 5, 4), Grid1D.Linspace(-2, 2, 3));//vertices with zero values
//            var grd = new BilinearSquareGrid(Grid1D.Linspace(-5, 5, 3), Grid1D.Linspace(-2, 2, 2));//cells- neighbors only in one direction
//            //var grd = new BilinearSquareGrid(Grid1D.Linspace(-4, 5, 4), Grid1D.Linspace(-4, 5, 4));//cells- neighbors only in two directions
//            //var grd = new BilinearSquareGrid(Grid1D.Linspace(-4, 5, 2), Grid1D.Linspace(-2, 2, 2)); //grid with only one cell

//            GridData grdData = new GridData(grd);


//            Basis b = new Basis(grdData, 2);

//            Field f = new SinglePhaseField(b, "f");

//            //f.ProjectField(delegate(MultidimensionalArray x, MultidimensionalArray r) {
//            //    int I = x.GetLength(0);
//            //    for (int i = 0; i < I; i++) {
//            //        r[i] = x[i, 0] + x[i, 1];
//            //    }
//            //});

//            f.ProjectField((x, y, z) => (x + y));


//            //var vtxLoc = grd.RefElements[0].Vertices;
//            //var vtxGlb = MultidimensionalArray.Create(1,vtxLoc.GetLength(0),vtxLoc.GetLength(1));
//            //var vtxLocBack = MultidimensionalArray.Create(1,vtxLoc.GetLength(0),vtxLoc.GetLength(1));
//            //grdData.TransformLocal2Global(vtxLoc, vtxGlb, 0, 1, 0);
//            //grdData.TransformGlobal2Local(vtxGlb.ExtractSubArrayShallow(0, -1, -1),vtxLocBack, 0, 1, 0);


//            int J = grdData.Cells.NoOfLocalUpdatedCells;

//            for (int j = 0; j < J; j++) {
//                double v = f.GetMeanValue(j);
//                f.SetMeanValue(j, v );
//               Console.WriteLine(v);
//            }


//            Tecplot.PlotFields(new Field[] { f }, null, grdData, "test", "test", 0, 0);


//            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
//            csMPI.Raw.mpiFinalize();
//        }
//    }
//}

using System.Collections.Generic;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using System.Diagnostics;
using System.IO;
using log4net;

using BoSSS.Solution.Utils;
using BoSSS.Platform.LinAlg;
using System;
using BoSSS.Solution;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
//using BoSSS.Solution.CgnsExport;
using ilPSP.Tracing;
using BoSSS.Solution.GridImport;

namespace CurvedElementsTest {

    /// <summary>
    /// Application class of my own BoSSS Application
    /// </summary>
    class Program : Application {

        /// <summary>
        /// Entry point of my program
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            Application._Main(args, true, "ScalarTransport,BoSSS.Solution,BoSSS.Foundation", delegate() { return new Program(); });
        }

        /// <summary>
        /// creates a simple 2d/3d bilinear grid;
        /// </summary>


        protected override GridCommons CreateOrLoadGrid() {
            // 2d Grid 
            // =======     
            var xNodes = Grid1D.Linspace(-1, 1, 3);
            var yNodes = Grid1D.Linspace(-1, 1, 3);
            var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, CellType.Square_8, false, false);
            //var grd = Grid2D.BilinearSquareGrid(Grid1D.Linspace(-100, 100, 8), Grid1D.Linspace(-100, 100, 8), 0.8, 0, 0, null);
            //var grd = Grid2D.CurvedSquareGrid(Grid1D.Linspace(1, 2, 7), Grid1D.Linspace(0, 1, 16), 9);
            // var grd = Grid2D.CurvedSquareGrid(Grid1D.Linspace(1, 2, 3), Grid1D.Linspace(0, 1, 101), CellType.Square_9, true);
            //var grd = Grid3D.Cartesian3DGrid(Grid1D.Linspace(0, 1, 3), Grid1D.Linspace(0, 1, 3), Grid1D.Linspace(0, 1, 3),false, false, false,CellType.Cube_8);
            //var grd = Grid3D.Cartesian3DGrid(Grid1D.Linspace(-1, 1, 4), Grid1D.Linspace(-1, 1, 4), Grid1D.Linspace(-1, 1, 4),true, true, true, CellType.Cube_Linear);
            //var grd = Grid3D.CylinderGrid(Grid1D.Linspace(1, 2, 4), Grid1D.Linspace(0, 0.25, 4), Grid1D.Linspace(0, 2, 4), CellType.Cube_27, true, true);
            return grd;
        }


        /// <summary>
        /// the one and only field in this application
        /// </summary>
        Field u;


        Field dx_u;
        Field dy_u;

        Field dx_u_byflux;
        Field dy_u_byflux;
        Field dz_u_byflux;

        Field du_dx_Analytical;
        Field du_dy_Analytical;
        Field du_dz_Analytical;


        /// <summary>
        /// creates the field <see cref="u"/>;
        /// </summary>
        protected override void CreateFields() {

            Console.Write("Creating Fields ... \t");
            Basis uBasis = new Basis(base.GridDat,4);
            u = new SinglePhaseField(uBasis, "u");
            m_IOFields.Add(u);
            dx_u = new SinglePhaseField(uBasis, "dx_u");
            dy_u = new SinglePhaseField(uBasis, "dy_u");
            dx_u_byflux = new SinglePhaseField(uBasis, "dx_u_byflux");
            dy_u_byflux = new SinglePhaseField(uBasis, "dy_u_byflux");
            dz_u_byflux = new SinglePhaseField(uBasis, "dz_u_byflux");
            du_dx_Analytical = new SinglePhaseField(uBasis, "du_dx_Analytical");
            du_dy_Analytical = new SinglePhaseField(uBasis, "du_dy_Analytical");
            du_dz_Analytical = new SinglePhaseField(uBasis, "du_dz_Analytical");
            Console.WriteLine(" done");
        }

        /// <summary>
        /// timestepper for this method;
        /// </summary>
        ExplicitEuler eEule;


        /// <summary>
        /// 
        /// </summary>
        protected override void CreateEquationsAndSolvers() {
            SpatialOperator diffOp = new SpatialOperator(1, 0, 1, "u", "codom1");


            switch (base.GridDat.Grid.SpatialDimension) {
                case 2: diffOp.EquationComponents["codom1"].Add(new ScalarTransportFlux()); break;
                case 3: diffOp.EquationComponents["codom1"].Add(new ScalarTransportFlux3D()); break;
                default:
                    throw new NotImplementedException("spatial dim. not supported.");
            }
            diffOp.Commit();

            eEule = new RungeKutta(
                                   RungeKutta.RungeKuttaScheme.TVD3,
                                   diffOp,
                                   new CoordinateMapping(u), null);
        }

        /// <summary>
        /// performs some timesteps
        /// </summary>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                /*
                u.ProjectField((x, y, z) => (x*x*3));

                //dx_u.Derivative(1.0, u, 0);
                //dy_u.Derivative(1.0, u, 1);

                dx_u_byflux.DerivativeByFlux(1.0, u, 0);
                //dy_u_byflux.DerivativeByFlux(1.0, u, 1);

                double min, max;
                dx_u_byflux.GetExtremalValues(out min, out max);
                Console.WriteLine(min + " < dx_u_byflux < " + max);

                base.TerminationKey = true;

                base.NoOfTimesteps = -1;
                return dt;

                */

                ScalarFunction sf = ((_2D)((x, y) => 1-x*x - y*y)).Vectorize();


                u.ProjectField(sf);
                var err = u.L2Error(sf, 4);


                               
                double min_x, max_x, min_y, max_y, min_u, max_u;
                dx_u_byflux.GetExtremalValues(out min_x, out max_x);
                dy_u_byflux.GetExtremalValues(out min_y, out max_y);
                u.GetExtremalValues(out min_u, out max_u);

                Console.WriteLine(min_u + " < u < " + max_u);
                Console.WriteLine(min_x + " < dx_u_byflux < " + max_x);
                Console.WriteLine(min_y + " < dy_u_byflux < " + max_y);

                //double Apl, Ami;
                //EvaluateLevelSetVolume(u, out Apl, out Ami);

                if (dt <= 0) {
                    // if dt <= 0, we're free to set the timestep on our own
                    //base.NoOfTimesteps = -1;
                    //dt = 1;

                    base.NoOfTimesteps = 1;
                    base.EndTime = Math.PI / 200.0;
                    dt = base.EndTime / ((double)base.NoOfTimesteps);


                    base.NoOfTimesteps *= 1; // run till' infinity
                    base.EndTime = 1.0e20;

                   
                }
                
                if (base.GridDat.MyRank  == 0)
                    Console.Write("Timestp. #" + TimestepNo + " of " + base.NoOfTimesteps + " ... \t");
                //eEule.Perform(dt);

                //if ((TimestepNo) % 1 == 0)
                //    PlotCurrentState(phystime, (int)TimestepNo);
                if (base.GridDat.MyRank == 0)
                    Console.WriteLine("finished!");
                

                // print diagnostic:
                if (TimestepNo == 1) {
                    du_dx_Analytical.Acc(-1.0, dx_u_byflux);
                    du_dy_Analytical.Acc(-1.0, dy_u_byflux);

                    double err_dx = du_dx_Analytical.L2Norm();
                    Console.WriteLine("|| df/dx_Numerical - df/dx_Analytical ||_2 = " + err_dx);
                    double err_dy = du_dy_Analytical.L2Norm();
                    Console.WriteLine("|| df/dy_Numerical - df/dy_Analytical ||_2 = " + err_dy);
                    Console.WriteLine("--------------------------------------------");
                }
                return dt;           
                 

            }
        }


        /// <summary>
        /// performs tecplot output of field <see cref="u"/>
        /// </summary>
        /// <param name="timestepNo"></param>
        /// <param name="phystime"></param>
        protected override void PlotCurrentState(double phystime, int timestepNo) {
            Tecplot plt2 = new Tecplot(base.GridDat, false, false, 3);
            switch (GridDat.SpatialDimension) {
                case 2:
                    plt2.PlotFields("transport_Recon_" + timestepNo, "Scalar Transport", phystime, new Field[] { u, dx_u_byflux, dy_u_byflux, du_dx_Analytical, du_dy_Analytical });
                    break;
                case 3:
                    plt2.PlotFields("transport_Recon_" + timestepNo, "Scalar Transport", phystime, new Field[] { u, dx_u_byflux, dy_u_byflux, dz_u_byflux, du_dx_Analytical, du_dy_Analytical, du_dz_Analytical });
                    break;
            }

}

        /// <summary>
        /// sets some initial value for field <see cref="u"/>;
        /// </summary>
        protected override void SetInitial() {
            Console.Write("Projecting Fields ... \t");
            switch (GridDat.SpatialDimension) {
                case 3:
                    u.ProjectField((x, y, z) => (3 * x * x + 3 * y * y + 3 * z * z));
                    dx_u_byflux.DerivativeByFlux(1.0, u, 0);
                    dy_u_byflux.DerivativeByFlux(1.0, u, 1);
                    dz_u_byflux.DerivativeByFlux(1.0, u, 2);
                    du_dx_Analytical.ProjectField((x, y, z) => (6 * x));
                    du_dy_Analytical.ProjectField((x, y, z) => (6 * y));
                    du_dz_Analytical.ProjectField((x, y, z) => (6 * z));
                    break;
                case 2:
                    u.ProjectField((x, y) => (3 * x * x + 3 * y * y));
                    dx_u_byflux.DerivativeByFlux(1.0, u, 0);
                    dy_u_byflux.DerivativeByFlux(1.0, u, 1);
                    du_dx_Analytical.ProjectField((x, y) => (6 * x));
                    du_dy_Analytical.ProjectField((x, y) => (6 * y));
                    break;
                default:
                    throw new InvalidDataException("Spatial Dimension invalid");
            }
            Console.WriteLine(" done");
            //u.ProjectField(Jump);
            //double rank = base.GridDat.MyRank;
            //int J = base.GridDat.Grid.NoOfUpdateCells;
            

            PlotCurrentState(0, 0);
        }

        /// <summary>
        /// (one possible) initial value for field <see cref="u"/>
        /// </summary>
        /// <param name="input">positions in space at which the function should be evaluated;
        /// 1st index: point index;
        /// 2nd index: spatial coordinate vector
        /// </param>
        /// <param name="output">result of function evaluation;
        /// 1st index: point index, corresponds with 1st index of <paramref name="input"/>
        /// </param>
        void Jump(MultidimensionalArray input, MultidimensionalArray output) {
            int D = input.GetLength(1);

            double max = 1;
            double min = 0;

            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];
                double z = double.NaN; if (D > 2) z = input[i, 2] - 5;

                if (Math.Abs(x) <= 3.0 && Math.Abs(y) <= 3.0) {
                    if (D > 2) {
                        // 3D - brach
                        if (Math.Abs(z) <= 3)
                            output[i] = max;
                        else
                            output[i] = min;
                    } else {
                        output[i] = max;
                    }
                } else {
                    output[i] = min;
                }
            }
        }

        /// <summary>
        /// function identical to 1.0
        /// </summary>
        /// <param name="input"></param>
        /// <param name="output"></param>
        void One(MultidimensionalArray input, MultidimensionalArray output) {

            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];

                output[i] = 1.0;
            }
        }
    }
}
