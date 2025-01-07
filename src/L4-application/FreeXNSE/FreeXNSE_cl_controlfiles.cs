using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Platform.Utils.Geom;
using MathNet.Numerics.Interpolation;

namespace FreeXNSE {

    [TestFixture]
    public static class FreeXNSE_Contactline_test {

        [Test]
        public static void FreeXNSE_SlugInChannel_Static() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.SlugInChannel_Static();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_SlugInChannel_Equilibrium() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.SlugInChannel_Equilibrium();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_SlugInChannel_Couette() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.SlugInChannel_Couette();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_StaticDroplet_FixedInterface() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.StaticDroplet_FixedInterface();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_StaticDroplet_DynamicContactAngle() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.StaticDroplet_DynamicContactAngle();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_PlanarInterface_DynamicContactAngle() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.PlanarInterface_DynamicContactAngle();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_StaticDroplet_FreeContactAngle() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.StaticDroplet_FreeContactAngle();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_SlidingDroplet_TiltedPlane(double alpha = 0.0) {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.SlidingDroplet_TiltedPlane(alpha: alpha);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_SlidingDroplet_ContactAngleHysteresis(double theta_adv, double theta_rec) { 
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.SlidingDroplet_ContactAngleHysteresis(3, 4, theta_adv, theta_rec);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void AggregationFail() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_Contactline.SlidingDroplet_TiltedPlane();

            C.activeAMRlevelIndicators.Clear();
            int level = 3;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = level });

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

    }


    internal class FreeXNSE_Contactline {

        internal static FreeXNSE_Control SlugInChannel_Static(int GridRes = 3, int k = 4) {
            FreeXNSE_Control C = new FreeXNSE_Control(false);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;

            C.dtFixed = 0.01;
            C.NoOfTimesteps = 1;


            // Re
            C.DimensionlessNumbers.Oh = 1.0;


            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Off;
            C.ActiveTerms.Temporal = Temporal.Off;


            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-5, 5, 5 * GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Robin_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Robin_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Neumann_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            double U = 0.0;
            C.AddBoundaryValue("Robin_bottom", "VelocityX#A", $"X => {-U}", false);
            C.AddBoundaryValue("Robin_top", "VelocityX#A", $"X => {U}", false);


            C.AddInitialValue("Phi", "X => (X[0] - 1.1)*(X[0] + 1.1)", false);
            //C.AddInitialValue("Phi", "X => -1.0", false);

            C.DimensionlessNumbers.beta = 10.0; //double.PositiveInfinity; // freeslip (0.0)
            C.DimensionlessNumbers.Theta = Math.PI/2.0; //

            C.ReInitPeriod = 0;

            int level = 20;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });

            return C;
        }

        internal static FreeXNSE_Control SlugInChannel_Equilibrium(int GridRes = 3, int k = 4) {
            FreeXNSE_Control C = new FreeXNSE_Control(false);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;

            C.dtFixed = 0.01;
            C.NoOfTimesteps = 300;


            // Re
            C.DimensionlessNumbers.Oh = 1.0;


            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Off;
            C.ActiveTerms.Temporal = Temporal.Off;


            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-5, 5, 5 * GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Robin_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Robin_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Neumann_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            double U = 0.0;
            C.AddBoundaryValue("Robin_bottom", "VelocityX#A", $"X => {-U}", false);
            C.AddBoundaryValue("Robin_top", "VelocityX#A", $"X => {U}", false);


            C.AddInitialValue("Phi", "X => (X[0] - 1.1)*(X[0] + 1.1)", false);
            //C.AddInitialValue("Phi", "X => -1.0", false);

            C.SlipScaling = 0.0;
            C.DimensionlessNumbers.alpha = 1.0; //
            C.DimensionlessNumbers.beta = 0.0; // freeslip
            C.DimensionlessNumbers.Theta = Math.PI / 4.0; // 

            int level = 4;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });

            C.PostprocessingModules.Add(new ContactLineLogging());

            return C;
        }

        internal static FreeXNSE_Control SlugInChannel_Couette(int GridRes = 3, int k = 4) {
            FreeXNSE_Control C = new FreeXNSE_Control(false);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.CustomLevelSet;

            C.dtFixed = 0.01;
            C.NoOfTimesteps = 300;


            // Re
            C.DimensionlessNumbers.Oh = 1.0;


            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Off;
            C.ActiveTerms.Temporal = Temporal.Off;


            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-5, 5, 5 * GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Robin_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Robin_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Neumann_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            double U = 1.0;
            C.AddBoundaryValue("Robin_bottom", "VelocityX#A", $"X => {-U}", false);
            C.AddBoundaryValue("Robin_top", "VelocityX#A", $"X => {U}", false);

            if(C.Option_LevelSetEvolution != BoSSS.Solution.LevelSetTools.LevelSetEvolution.CustomLevelSet) {
                C.AddInitialValue("Phi", $"X => (X[0]-1.1)*(X[0]+1.1)", false);
            } else if(C.Option_LevelSetEvolution == BoSSS.Solution.LevelSetTools.LevelSetEvolution.CustomLevelSet) {
                C.DualSplinePhi0Initial = new Tuple<double[], double[]>[2];
                C.DualSplinePhi0Initial[0] = Tuple.Create(new double[] { -1, 1 }, new double[] {-1.1, -1.1 }); ;
                C.DualSplinePhi0Initial[1] = Tuple.Create(new double[] { -1, 1 }, new double[] { 1.1, 1.1 }); ;
            }

            C.SlipScaling = 0.0;
            C.ContactAngleScaling= 1.0;
            //C.ActiveTerms.SurfaceTension = SurfaceTension.Off;
            C.DimensionlessNumbers.alpha = 0.5; // Math.Sqrt(2) / 2.0;
            C.DimensionlessNumbers.beta = 0.0; // freeslip
            C.DimensionlessNumbers.Theta = Math.PI / 2.0; // 

            int level = 2;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = level });

            C.PostprocessingModules.Add(new ContactLineLogging());           

            return C;
        }

        internal static FreeXNSE_Control PlanarInterface_DynamicContactAngle(int GridRes = 3, int k = 5) {
            FreeXNSE_Control C = new FreeXNSE_Control(true);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;

            C.dtFixed = 0.01;
            C.NoOfTimesteps = 1;


            // Re
            C.DimensionlessNumbers.Oh = 1.0;

            C.ActiveTerms.SurfaceTension = SurfaceTension.LaplaceBeltrami;
            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Off;
            C.ActiveTerms.Temporal = Temporal.Off;

            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-5, 5, 5 * GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Robin_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Robin_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Neumann_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            double U = 0.0;
            C.AddBoundaryValue("Robin_bottom", "VelocityX#A", $"X => {-U}", false);
            C.AddBoundaryValue("Robin_top", "VelocityX#A", $"X => {U}", false);

            C.AddInitialValue("Phi", "X => X[0]", false);

            C.DimensionlessNumbers.alpha = 1.0;
            C.DimensionlessNumbers.beta = 1.0;
            C.DimensionlessNumbers.Theta = 0.0;

            int level = 10;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });

            return C;
        }

        private static FreeXNSE_Control StaticDroplet_Base(int GridRes, int k, double theta, double R) {
            FreeXNSE_Control C = new FreeXNSE_Control(true);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "StaticDroplet";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;

            C.dtFixed = 0.01;
            C.NoOfTimesteps = 1;

            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Off;
            C.ActiveTerms.Temporal = Temporal.Off;

            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1.1*R, 1.1*R, GridRes + 1);
                var _yNodes = GenericBlas.Linspace(0, 2.2*R, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Robin_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Neumann_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Neumann_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };            

            double yOffset = -R * Math.Cos(theta);
            C.AddInitialValue("Phi", $"X => Math.Sqrt(Math.Pow(X[0], 2) + Math.Pow(X[1] - {yOffset},2)) - {R}", false);

            int level = 20;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });

            return C;
        }

        internal static FreeXNSE_Control StaticDroplet_FixedInterface(int GridRes = 3, int k = 5) {
            FreeXNSE_Control C = StaticDroplet_Base(GridRes, k, Math.PI / 4.0, 10.0);

            // Re
            C.DimensionlessNumbers.Oh = 1.0;

            C.ActiveTerms.SurfaceTension = SurfaceTension.Off;
            C.FixedInterface = true;            

            double U = -1.0;
            C.AddBoundaryValue("Robin_bottom", "VelocityX#A", $"X => {-U}", false);

            C.SlipScaling = -1.0;
            C.DimensionlessNumbers.beta = 1.0; // no-slip (double.PositiveInfinity); // freeslip (0.0)
            C.DimensionlessNumbers.Theta = 0; // ignore surface tension for now

            return C;
        }

        internal static FreeXNSE_Control StaticDroplet_DynamicContactAngle(int GridRes = 3, int k = 5) {
            FreeXNSE_Control C = StaticDroplet_Base(GridRes, k, Math.PI / 4.0, 10.0);

            // Re
            C.DimensionlessNumbers.Oh = 1.0;

            C.ActiveTerms.SurfaceTension = SurfaceTension.LaplaceBeltrami;
            C.FixedInterface = false;

            C.SlipScaling = 0.0;
            C.DimensionlessNumbers.beta = 0.0; // no-slip (double.PositiveInfinity); // freeslip (0.0)
            C.DimensionlessNumbers.Theta = 0;
            C.DimensionlessNumbers.alpha = 1.0;

            return C;
        }

        internal static FreeXNSE_Control StaticDroplet_FreeContactAngle(int GridRes = 3, int k = 5) {
            FreeXNSE_Control C = StaticDroplet_Base(GridRes, k, Math.PI / 2.0, 10.0);

            // Re
            C.DimensionlessNumbers.Oh = 1.0;

            C.ActiveTerms.SurfaceTension = SurfaceTension.LaplaceBeltrami;
            C.FixedInterface = false;

            C.SlipScaling = 1.0;
            C.DimensionlessNumbers.beta = 1.0; // no-slip (double.PositiveInfinity); // freeslip (0.0)
            C.DimensionlessNumbers.Theta = 0; // ignore surface tension for now
            C.DimensionlessNumbers.alpha = -1.0;

            C.NoOfTimesteps = 200;

            C.activeAMRlevelIndicators.Clear();
            int level = 4;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = level-2 });

            return C;
        }

        private static FreeXNSE_Control SlidingDroplet_Base(int GridRes, int k, double theta, double R) {
            FreeXNSE_Control C = new FreeXNSE_Control(true);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "SlidingDroplet";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;

            C.dtFixed = 0.1;
            C.NoOfTimesteps = 1;

            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Off;
            C.ActiveTerms.Temporal = Temporal.Off;

            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-2.2 * R, 6.6 * R, 4 * GridRes + 1);
                var _yNodes = GenericBlas.Linspace(0, 2.2 * R, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Robin_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Neumann_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Neumann_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            double yOffset = -R * Math.Cos(theta);
            C.AddInitialValue("Phi", $"X => Math.Sqrt(Math.Pow(X[0], 2) + Math.Pow(X[1] - {yOffset},2)) - {R}", false);

            int level = 20;
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });

            return C;
        }

        
        internal static FreeXNSE_Control SlidingDroplet_TiltedPlane(int GridRes = 3, int k = 4, double alpha = 0.0) {
            FreeXNSE_Control C = SlidingDroplet_Base(GridRes, k, Math.PI / 2.0, 10.0);

            // Re
            C.DimensionlessNumbers.Oh = 1.0;
            C.DimensionlessNumbers.Fr = 100.0;
            C.VolumeForce = new ConstantGravity(new double[] { Math.Sin(alpha), -Math.Cos(alpha) });
            C.ActiveTerms.SurfaceTension = SurfaceTension.LaplaceBeltrami;
            C.FixedInterface = false;

            C.SlipScaling = 0.0;
            C.DimensionlessNumbers.beta = 0.1; // no-slip (double.PositiveInfinity); // freeslip (0.0)
            C.DimensionlessNumbers.Theta = Math.PI / 2.0; // ignore surface tension for now
            C.DimensionlessNumbers.ThetaRec = 0.0; // static value should be Math.PI/3
            C.DimensionlessNumbers.ThetaAdv = 106.584 / 180.0 * Math.PI;
            C.DimensionlessNumbers.alpha = 1.0;            

            C.NoOfTimesteps = 500;

            // fail with k=5, level=4, around ts 60, not 100% reproducible ?!
            C.activeAMRlevelIndicators.Clear();
            int level = 2; // fail for level = 3
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = level });

            //C.InitialValues.Clear();
            //C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.CustomLevelSet;
            //C.SemiCircleSplinePhi0Initial = a => new[] { 9 * Math.Cos(a), 9 * Math.Sin(a) };
            //C.NoOfNodes = 20;

            C.PostprocessingModules.Add(new ContactLineLogging());

            return C;
        }

        internal static FreeXNSE_Control SlidingDroplet_ContactAngleHysteresis(int GridRes = 3, int k = 4, double theta_adv = 0.0, double theta_rec = Math.PI) {
            FreeXNSE_Control C = new FreeXNSE_Control(true);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "SlidingDroplet";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;

            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Off;
            C.ActiveTerms.Temporal = Temporal.Off;

            // degree
            C.SetDGdegree(k);

            // Create Grid
            double R = 1;
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 2.2 * R, GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-2.2 * R, 2.2 * R, 2 * GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Neumann_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Neumann_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Robin_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            C.AddInitialValue("Phi", $"X => Math.Sqrt(Math.Pow(X[0], 2) + Math.Pow(X[1],2)) - {R}", false);

            // Re
            C.DimensionlessNumbers.Oh = 1.0;
            Console.WriteLine("Weber number: {0}", C.DimensionlessNumbers.We);
            Console.WriteLine("Laplace number: {0}", C.DimensionlessNumbers.La);
            C.DimensionlessNumbers.Fr = Math.PI/2.0*R*R* C.DimensionlessNumbers.We;
            C.VolumeForce = new ConstantGravity(new double[] { 0.0, -1.0 });
            C.ActiveTerms.SurfaceTension = SurfaceTension.LaplaceBeltrami;
            C.FixedInterface = false;

            C.SlipScaling = 0.0;
            C.DimensionlessNumbers.beta = 0.0; // no-slip (double.PositiveInfinity); // freeslip (0.0)
            C.DimensionlessNumbers.Theta = Math.PI / 2.0;
            C.DimensionlessNumbers.ThetaRec = theta_rec;
            C.DimensionlessNumbers.ThetaAdv = theta_adv;
            C.DimensionlessNumbers.alpha = 0.0;

            C.dtFixed = 0.01;
            C.NoOfTimesteps = 500;

            // fail with k=5, level=4, around ts 60, not 100% reproducible ?!
            int level = 2; // fail for level = 3
            C.AdaptiveMeshRefinement = level > 0;
            C.AMR_startUpSweeps = level;
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = level });
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = level });

            //C.InitialValues.Clear();
            //C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.CustomLevelSet;
            //C.SemiCircleSplinePhi0Initial = a => new[] { 9 * Math.Cos(a), 9 * Math.Sin(a) };
            //C.NoOfNodes = 20;

            C.PostprocessingModules.Add(new ContactLineLogging());

            return C;
        }

    }
}
