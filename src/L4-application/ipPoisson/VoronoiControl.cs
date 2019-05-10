#define LOGGING
using System;
using System.Collections.Generic;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid;



namespace BoSSS.Application.SipPoisson{

    static public class VoronoiControl {
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

            VoronoiGrid grid = new VoronoiGrids.LDomain(Res, NoOfLlyodsIter);

            return TestGrid(grid, deg, solver_name, db);
        }

        /// <summary>
        /// Refines a voronoi mesh on a L-shaped domain by increasing the number of Lloyd-iterations.  
        /// </summary>
        /// <param name="res">
        /// number of randomly chosen Delaunay vertices
        /// </param>
        /// <returns></returns>
        public static SipControl[] VoronoiPStudy(int res)
        {
            double[] NoOfLli = GenericBlas.Linspace(0, 100, 101);
            List<SipControl> R = new List<SipControl>();
            for (int i = 0; i < NoOfLli.Length; i++)
            {
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

            R.ImmediatePlotPeriod = 1;
            if (db != null)
            {
                R.savetodb = true;
                R.SetDatabase(db);
            }

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg * 2 });
            R.InitialValues_Evaluators.Add("RHS", X => 1.0);
            R.InitialValues_Evaluators.Add("Tex", X => X[0]);
            R.ExactSolution_provided = false;
            R.LinearSolver.NoOfMultigridLevels = int.MaxValue;
            R.LinearSolver.SolverCode = solver_name;
            R.LinearSolver.NoOfMultigridLevels = 1;
            //R.TargetBlockSize = 100;

            grid.SetGridAndBoundaries(R);
            return R;
        }
    }

    abstract class VoronoiGrid
    {
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
        public class LDomainLine : VoronoiGrid
        {
            protected override IGrid GridFunc()
            {
                Vector[] DomainBndyPolygon = new[] {
                    new Vector(0,0),
                    new Vector(-1,0),
                    new Vector(-1,1),
                    new Vector(1,1),
                    new Vector(1,-1),
                    new Vector(0,-1),
                };
                AggregationGrid grid;
                MultidimensionalArray nodes = GetNodes();
                grid = BoSSS.Foundation.Grid.Voronoi.VoronoiGrid2D.FromPolygonalDomain(nodes, DomainBndyPolygon, 0, 20);
                grid.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grid.DefineEdgeTags(X => (byte)1);
                return grid;
            }

            MultidimensionalArray GetNodes()
            {
                MultidimensionalArray nodes = MultidimensionalArray.Create(40, 2);
                double[] x = ilPSP.Utils.GenericBlas.Linspace(-1, 1, 40);
                nodes.SetColumn(0, x);
                nodes.SetColumn(1, x);
                return nodes;
            }

            protected override void SetBoundaryValues(AppControl R)
            {
                R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                    X => Math.Pow(X[0], 2) + Math.Pow(X[1], 2));
            }
        }

        public class LDomain : VoronoiGrid
        {
            readonly int res;
            readonly int noOfLlyodsIter;
            public LDomain(int res, int noOfLlyodsIter)
            {
                this.res = res;
                this.noOfLlyodsIter = noOfLlyodsIter;
            }

            Vector[] LShape()
            {
                double a = 10000;
                Vector[] LShapedPolygon = new[] {
                    new Vector(-a,a),
                    new Vector(a,a),
                    new Vector(a,-a),
                    new Vector(0,-a),
                    new Vector(0,0),
                    new Vector(-a,0)
                };
                return LShapedPolygon;
            }

            protected override IGrid GridFunc()
            {
#if LOGGING
                Console.WriteLine("Calculating Grid...");
#endif
                Vector[] DomainBndyPolygon = LShape();
                AggregationGrid grid;
                grid = BoSSS.Foundation.Grid.Voronoi.VoronoiGrid2D.FromPolygonalDomain(DomainBndyPolygon, noOfLlyodsIter, res);
                grid.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grid.DefineEdgeTags(X => (byte)1);
#if LOGGING
                Console.WriteLine("Done.");
#endif
                return grid;
            }

            protected override void SetBoundaryValues(AppControl R)
            {
                R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                    X => Math.Pow(X[0], 2) + Math.Pow(X[1], 2));
            }
        }

        public class Square : VoronoiGrid
        {
            readonly int res;
            readonly int noOfLlyodsIter;
            public Square(int res, int noOfLlyodsIter)
            {
                this.res = res;
                this.noOfLlyodsIter = noOfLlyodsIter;
            }

            protected override IGrid GridFunc()
            {
                Vector[] DomainBndyPolygon = new[] {
                    new Vector(-1,1),
                    new Vector(1,1),
                    new Vector(1,-1),
                    new Vector(-1,-1)
                };
                AggregationGrid grid;
                grid = BoSSS.Foundation.Grid.Voronoi.VoronoiGrid2D.FromPolygonalDomain(DomainBndyPolygon, noOfLlyodsIter, res);
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
