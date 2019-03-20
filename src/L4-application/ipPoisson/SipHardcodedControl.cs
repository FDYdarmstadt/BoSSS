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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Connectors.Matlab;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Gnuplot;
using BoSSS.Foundation.Grid;

namespace BoSSS.Application.SipPoisson {

    /// <summary>
    /// predefined control-objects
    /// </summary>
    static public class SipHardcodedControl {

        /// <summary>
        /// Test on a curved grid.
        /// </summary>
        public static SipControl TestCurved() {
            var R = new SipControl();
            R.ProjectName = "ipPoison/curved";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 15 });
            R.InitialValues_Evaluators.Add("RHS", X => 0.0);
            R.InitialValues_Evaluators.Add("Tex", X => (Math.Log(X[0].Pow2() + X[1].Pow2()) / Math.Log(4.0)) + 1.0);
            R.ExactSolution_provided = true;

            R.GridFunc = delegate () {
                var grd = Grid2D.CurvedSquareGrid(GenericBlas.Linspace(1, 2, 3), GenericBlas.Linspace(0, 1, 11), CellType.Square_9, true);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(X => 1);
                return grd;
            };

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     double x = X[0], y = X[1];
                     return Math.Sqrt(x * x + y * y);
                 });


            return R;
        }

        /// <summary>
        /// Creates Nodes, yes it really does!
        /// </summary>
        /// <param name="res"></param>
        /// <param name="stetch">
        /// Factor which determines how much the intervals in the output grow; 1.0 is no stretching.
        /// </param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        static double[] CreateNodes(int res, double stetch, double min, double max) {
            if (stetch == 1.0)
                return GenericBlas.Linspace(min, max, res + 1);
            else
                return Grid1D.ExponentialSpaceing(min, max, res + 1, stetch); // without proper preconditioning,
            // a stretched grid is much more expensive than
            // an equidistant grid !!!
        }

        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// </summary>
        public static SipControl TestCartesian1(int xRes = 32, double xStretch = 1.0, int yRes = 16, double yStretch = 1.01, int pDG = 2) {
            var RR = new SipControl();
            RR.ProjectName = "ipPoison/cartesian";
            RR.savetodb = false;

            RR.FieldOptions.Add("T", new FieldOpts() { Degree = pDG, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            RR.FieldOptions.Add("Tex", new FieldOpts() { Degree = pDG * 2 });
            RR.InitialValues_Evaluators.Add("RHS", X => 1.0);
            RR.InitialValues_Evaluators.Add("Tex", X => (0.5 * X[0].Pow2() - 10 * X[0]));
            RR.ExactSolution_provided = true;

            RR.GridFunc = delegate () {
                double[] xNodes = CreateNodes(xRes, xStretch, 0, 10);
                double[] yNodes = CreateNodes(yRes, yStretch, -1, +1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if (Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };


            RR.AddBoundaryValue(BoundaryType.Dirichlet.ToString());
            RR.AddBoundaryValue(BoundaryType.Neumann.ToString());


            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;


            return RR;
        }

        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// </summary>
        public static SipControl TestCartesian3D(int xRes = 32, double xStretch = 1.0, int yRes = 16, double yStretch = 1.0, int zRes = 16, double zStretch = 1.0) {
            var R = new SipControl();
            R.ProjectName = "ipPoison/cartesian";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = 6, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 6 });
            R.InitialValues_Evaluators.Add("RHS", X => 1.0);
            R.InitialValues_Evaluators.Add("Tex", X => (0.5 * X[0].Pow2() - 10 * X[0]));
            R.ExactSolution_provided = true;

            R.GridFunc = delegate () {
                double[] xNodes = CreateNodes(xRes, xStretch, 0, 10);
                double[] yNodes = CreateNodes(yRes, yStretch, -1, +1);
                double[] zNodes = CreateNodes(zRes, zStretch, -1, +1);

                var grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if (Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString());
            R.AddBoundaryValue(BoundaryType.Neumann.ToString());


            return R;
        }





        /// <summary>
        /// Test on a Cartesian grid, with a sinusodial solution.
        /// </summary>
        /// <param name="Res">
        /// Grid resolution
        /// </param>
        /// <param name="Dim">
        /// spatial dimension
        /// </param>
        /// <param name="deg">
        /// polynomial degree
        /// </param>
        /// <param name="solver_name">
        /// Name of solver to use.
        /// </param>
        public static SipControl TestCartesian2(int Res, int Dim, LinearSolverConfig.Code solver_name = LinearSolverConfig.Code.exp_softpcg_schwarz_directcoarse, int deg = 3) {
            if (Dim != 2 && Dim != 3)
                throw new ArgumentOutOfRangeException();

            var R = new SipControl();
            R.ProjectName = "ipPoison/cartesian";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg * 2 });
            R.InitialValues_Evaluators.Add("RHS", X => -Math.Sin(X[0]));
            R.InitialValues_Evaluators.Add("Tex", X => Math.Sin(X[0]));
            R.ExactSolution_provided = true;
            R.LinearSolver.NoOfMultigridLevels = int.MaxValue;
            R.LinearSolver.SolverCode = solver_name;
            //R.TargetBlockSize = 100;

            R.TracingNamespaces = "BoSSS,ilPSP";


            R.GridFunc = delegate () {
                GridCommons grd = null;
                if (Dim == 2) {
                    double[] xNodes = GenericBlas.Linspace(0, 10, Res * 5 + 1);
                    double[] yNodes = GenericBlas.SinLinSpacing(-1, +1, 0.6, Res + 1);

                    grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                } else if (Dim == 3) {
                    double[] xNodes = GenericBlas.Linspace(0, 10, Res * 5 + 1);
                    double[] yNodes = GenericBlas.SinLinSpacing(-1, +1, 0.6, Res + 1);
                    double[] zNodes = GenericBlas.SinLinSpacing(-1, +1, 0.6, Res + 1);

                    grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                } else {
                    throw new NotSupportedException();
                }
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if (Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     double x = X[0], y = X[1];

                     if (Math.Abs(X[0] - (0.0)) < 1.0e-8)
                         return 0.0;

                     throw new ArgumentOutOfRangeException();
                 });

            R.AddBoundaryValue(BoundaryType.Neumann.ToString(), "T",
                 delegate (double[] X) {
                     if (Math.Abs(X[1] - 1.0) < 1.0e-8 || Math.Abs(X[1] + 1.0) < 1.0e-8) // y = -1, y = +1
                         return 0;

                     if (X.Length > 2 && (Math.Abs(X[2] - 1.0) < 1.0e-8 || Math.Abs(X[2] + 1.0) < 1.0e-8)) // z = -1, z = +1
                         return 0;

                     if (Math.Abs(X[0] - (+10.0)) < 1.0e-8)
                         return Math.Cos(10.0);

                     throw new ArgumentOutOfRangeException();
                 });
            return R;
        }


        /// <summary>
        /// Poisson Equation on a (-1,1)x(-1,1), Dirichlet everywhere
        /// </summary>
        public static SipControl Square(int xRes = 5, int yRes = 5, int deg = 5) {

            //Func<double[], double> exRhs = X => 2 * X[0] * X[0] + 2 * X[1] * X[1] - 4;
            //Func<double[], double> exSol = X => (1.0 - X[0] * X[0]) * (1.0 - X[1] * X[1]);

            //Func<double[], double> exSol = X => (1.0 - X[1]);
            //Func<double[], double> exRhs = X => 0.0;

            Func<double[], double> exSol = X => -Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5);
            Func<double[], double> exRhs = X => (Math.PI * Math.PI * 0.5 * Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5)); // == - /\ exSol


            var R = new SipControl();
            R.ProjectName = "ipPoison/square";
            R.savetodb = false;
            //R.DbPath = "D:\\BoSSS-db";

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 4 });
            R.InitialValues_Evaluators.Add("RHS", exRhs);
            R.InitialValues_Evaluators.Add("Tex", exSol);
            R.ExactSolution_provided = true;
            //R.LinearSolver.NoOfMultigridLevels = 2;
            //R.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_softpcg_mg;
            R.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_softpcg_schwarz_directcoarse;
            R.SuppressExceptionPrompt = true;
            //R.LinearSolver.SolverCode = LinearSolverConfig.Code.classic_mumps;

            R.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(-1, 1, xRes);
                double[] yNodes = GenericBlas.Linspace(-1, 1, yRes);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret = 1;
                    return ret;
                });


                return grd;
            };

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T", exSol);

            R.NoOfSolverRuns = 1;

            R.AdaptiveMeshRefinement = true;
            R.NoOfTimesteps = 5;


            return R;
        }

        /// <summary>
        /// Test on a square 2D Voronoi mesh
        /// </summary>
        /// <param name="Res">
        /// number of randomly chosen Delaunay vertices
        /// </param>
        /// <param name="deg">
        /// polynomial degree
        /// </param>
        /// <param name="NoOfLlyodsIter">
        /// Number of Llyods iterations. 
        /// </param>
        /// <param name="mirror">
        /// Mirror vertices of boundary cells along boundary to approximate boundary
        /// </param>
        /// <param name="solver_name">
        /// Name of solver to use.
        /// </param>
        /// <returns></returns>
        public static SipControl TestVoronoi_Square(
            int Res,
            int deg = 1,
            int NoOfLlyodsIter = 10,
            bool mirror = false,
            LinearSolverConfig.Code solver_name = LinearSolverConfig.Code.classic_pardiso,
            Foundation.IO.IDatabaseInfo db = null)
        {
            return TestGrid(new VoronoiGrids.Square(Res, NoOfLlyodsIter), deg, solver_name, db);
        }

        /// <summary>
        /// Test on a L-shaped 2D Voronoi mesh
        /// </summary>
        /// <param name="Res">
        /// number of randomly chosen Delaunay vertices
        /// </param>
        /// <param name="deg">
        /// polynomial degree
        /// </param>
        /// <param name="NoOfLlyodsIter">
        /// Number of Llyods iterations. 
        /// </param>
        /// <param name="mirror">
        /// Mirror vertices of boundary cells along boundary to approximate boundary
        /// </param>
        /// <param name="solver_name">
        /// Name of solver to use.
        /// </param>
        /// <returns></returns>
        public static SipControl TestVoronoi_LDomain(
            int Res,
            int deg = 1,
            int NoOfLlyodsIter = 10,
            bool mirror = false,
            LinearSolverConfig.Code solver_name = LinearSolverConfig.Code.classic_pardiso,
            Foundation.IO.IDatabaseInfo db = null)
        {
            return TestGrid(new VoronoiGrids.LDomain(Res, NoOfLlyodsIter), deg, solver_name, db);
        }

        /// <summary>
        /// Refines a voronoi mesh on a L-shaped domain by increasing the number of Lloyd-iterations.  
        /// </summary>
        /// <param name="res">
        /// number of randomly chosen Delaunay vertices
        /// </param>
        /// <returns></returns>
        public static SipControl[] VoronoiPStudy(int res) {
            double[] NoOfLli = GenericBlas.Linspace(0, 100, 101);
            List<SipControl> R = new List<SipControl>();
            for(int i = 0; i < NoOfLli.Length; i++) {
                R.Add(TestVoronoi_LDomain(res, deg: 1, NoOfLlyodsIter: (int)NoOfLli[i]));
            }

            return R.ToArray();
        }

        //Base case for Voronoi Testing
        static SipControl TestGrid(
            VoronoiGrid grid,
            int deg = 1,
            LinearSolverConfig.Code solver_name = LinearSolverConfig.Code.classic_pardiso,
            Foundation.IO.IDatabaseInfo db = null)
        {
            var R = new SipControl
            {
                ProjectName = "SipPoisson-Voronoi",
                SessionName = "testrun"
            };

            if (db != null)
            {
                R.savetodb = true;
                R.SetDatabase(db);
            }
            R.ImmediatePlotPeriod = 1;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg * 2 });
            R.InitialValues_Evaluators.Add("RHS", X => 0.0);// -1.0 + X[0] * X[0]);
            R.InitialValues_Evaluators.Add("Tex", X => X[0]);
            R.ExactSolution_provided = false;
            R.LinearSolver.NoOfMultigridLevels = int.MaxValue;
            R.LinearSolver.SolverCode = solver_name;
            R.LinearSolver.NoOfMultigridLevels = 1;
            //R.TargetBlockSize = 100;

            grid.SetGridAndBoundaries(R);
            return R;
        }

        abstract class VoronoiGrid
        {
            protected readonly int NoOfLlyodsIter;
            protected readonly int Res;
            public VoronoiGrid(int res, int noOfLlyodsIter)
            {
                Res = res;
                NoOfLlyodsIter = noOfLlyodsIter;
            }

            public void SetGridAndBoundaries(AppControl R)
            {
                R.GridFunc = GridFunc;
                SetBoundaryValues(R);
            }

            protected abstract IGrid GridFunc();

            protected abstract void SetBoundaryValues(AppControl R);
        }

        static class VoronoiGrids
        {
            public class LDomain : VoronoiGrid
            {
                public LDomain(int res, int noOfLlyodsIter) : base(res, noOfLlyodsIter) { }
               
                protected override IGrid GridFunc()
                {
                    Vector[] DomainBndyPolygon = new[] {
                        new Vector(-1,1),
                        new Vector(1,1),
                        new Vector(1,-1),
                        new Vector(0,-1),
                        new Vector(0,0),
                        new Vector(-1,0)
                    };
                    AggregationGrid grid;
                    grid = BoSSS.Foundation.Grid.Voronoi.VoronoiGrid2D.FromPolygonalDomain(DomainBndyPolygon, NoOfLlyodsIter, Res);
                    grid.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                    grid.DefineEdgeTags(X => (byte)1);
                    return grid;
                }

                protected override void SetBoundaryValues(AppControl R)
                {
                    R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                        X => Math.Pow(X[0],2) + Math.Pow(X[1], 2));
                }
            }

            public class Square : VoronoiGrid
            {
                public Square(int res, int noOfLlyodsIter) : base(res, noOfLlyodsIter) { }

                protected override IGrid GridFunc()
                {
                    Vector[] DomainBndyPolygon = new[] {
                        new Vector(-1,1),
                        new Vector(1,1),
                        new Vector(1,-1),
                        new Vector(-1,-1)
                    };
                    AggregationGrid grid;
                    grid = BoSSS.Foundation.Grid.Voronoi.VoronoiGrid2D.FromPolygonalDomain(DomainBndyPolygon, NoOfLlyodsIter, Res);
                    grid.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                    grid.DefineEdgeTags(X => (byte)1);
                    return grid;
                }

                protected override void SetBoundaryValues(AppControl R)
                {
                    R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                        X => Math.Pow(X[0], 2) + Math.Pow(X[1], 2));
                }
            }
        }
    }
}
