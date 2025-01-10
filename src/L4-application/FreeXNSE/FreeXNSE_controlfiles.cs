using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FreeXNSE {
    internal static class FreeXNSE_controlfiles {
        internal static FreeXNSE_Control ChannelFlow(int GridRes = 8, int k = 2) {
            FreeXNSE_Control C = new FreeXNSE_Control();

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            

            // Re
            C.DimensionlessNumbers.Oh = 0.1;

            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 10, GridRes * 5 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Dirichlet_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Dirichlet_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Dirichlet_inlet";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Neumann_Outlet";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            // set boundary conditions
            C.AddBoundaryValue("Dirichlet_inlet", "VelocityX", "X => 1 - X[1] * X[1]", TimeDependent: false);

            return C;
        }

        internal static FreeXNSE_Control Ellipse(int GridRes = 7, int k = 3, double ecc = 1.0) {
            FreeXNSE_Control C = new FreeXNSE_Control(false);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.01;
            C.NoOfTimesteps = 100;               


            // Re
            C.DimensionlessNumbers.Oh = 1e-3;

            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Temam;


            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Dirichlet_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Dirichlet_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Dirichlet_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Dirichlet_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            C.InitialValues.Add("Phi", new Formula($"X => {ecc} * X[0] * X[0] + 1 / {ecc} * X[1] * X[1] - 0.25"));

            return C;
        }

        internal static FreeXNSE_Control EllipseParameterized(int GridRes = 5, int k = 4, double ecc = 1.0) {
            FreeXNSE_Control C = new FreeXNSE_Control(true);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.01;
            C.NoOfTimesteps = 250;


            // Re
            C.DimensionlessNumbers.Oh = 1e-3;


            C.ActiveTerms.Viscous = Viscous.OBB;
            C.ActiveTerms.Convective = Convective.ConservativeTemam;


            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Dirichlet_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Dirichlet_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Dirichlet_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Dirichlet_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            //{ecc} * X[0] * X[0] + 1 / {ecc} * X[1] * X[1] - 0.25
            List<ParameterFunctionPair> LevelSet = new List<ParameterFunctionPair>();
            LevelSet.Add(new ParameterFunctionPair(-0.25, X => 1.0, new Func<double[], double>[] { X => 0.0, X => 0.0 }));
            LevelSet.Add(new ParameterFunctionPair(ecc, X => X[0] * X[0], new Func<double[], double>[] { X => 2 * X[0], X => 0.0 }));
            LevelSet.Add(new ParameterFunctionPair(1.0/ecc, X => X[1] * X[1], new Func<double[], double>[] { X => 0.0, X => 2 * X[1] }));
            LevelSet.Add(new ParameterFunctionPair(0.0, X => X[0], new Func<double[], double>[] { X => 0.0, X => 2 * X[1] }));
            LevelSet.Add(new ParameterFunctionPair(0.0, X => X[1], new Func<double[], double>[] { X => 0.0, X => 2 * X[1] }));
            LevelSet.Add(new ParameterFunctionPair(0.0, X => X[0] * X[1], new Func<double[], double>[] { X => 0.0, X => 2 * X[1] }));

            C.ParameterLevelSet = LevelSet.ToArray();

            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.None;

            return C;
        }

        internal static FreeXNSE_Control ParameterizedAdvection(int GridRes = 5, int k = 2) {
            FreeXNSE_Control C = new FreeXNSE_Control(false);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.01;
            C.NoOfTimesteps = 10;


            // Re
            C.DimensionlessNumbers.Oh = 1e-3;


            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Temam;


            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, BoSSS.Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Dirichlet_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Dirichlet_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Dirichlet_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Dirichlet_right";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            // X[0] * X[0] +  X[1] * X[1] - 1.0
            List<ParameterFunctionPair> LevelSet = new List<ParameterFunctionPair>();
            //LevelSet.Add(new ParameterFunctionPair(-1.0, X => 1.0, new Func<double[], double>[] { X => 0.0, X => 0.0 }));
            //LevelSet.Add(new ParameterFunctionPair(0.0, X => X[0], new Func<double[], double>[] { X => 1.0, X => 0.0 }));
            ////LevelSet.Add(new ParameterFunctionPair(0.0, X => X[1], new Func<double[], double>[] { X => 0.0, X => 1.0 }));
            //LevelSet.Add(new ParameterFunctionPair(1.0, X => X[0] * X[0], new Func<double[], double>[] { X => 2 * X[0], X => 0.0 }));
            ////LevelSet.Add(new ParameterFunctionPair(0.0, X => X[0] * X[1], new Func<double[], double>[] { X => X[1], X => X[0] }));
            //LevelSet.Add(new ParameterFunctionPair(1.0, X => X[1] * X[1], new Func<double[], double>[] { X => 0.0, X => 2 * X[1] }));

            LevelSet.Add(new ParameterFunctionPair(0.0, X => 1.0, new Func<double[], double>[] { X => 0.0, X => 0.0 }));
            LevelSet.Add(new ParameterFunctionPair(1.0, X => X[0], new Func<double[], double>[] { X => 1.0, X => 0.0 }));


            C.ParameterLevelSet = LevelSet.ToArray();

            C.InitialValues.Add("VelocityX", new Formula($"X => 1.0"));
            C.AddBoundaryValue("Dirichlet_left", "VelocityX", new Formula($"X => 1.0"));
            C.AddBoundaryValue("Dirichlet_top", "VelocityX", new Formula($"X => 1.0"));
            C.AddBoundaryValue("Dirichlet_bottom", "VelocityX", new Formula($"X => 1.0"));


            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.None;
            

            return C;
        }

        internal static FreeXNSE_Control Ellipsoid(int GridRes = 5, int k = 3, double ecc = 1.0) {
            FreeXNSE_Control C = new FreeXNSE_Control(true);

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "Ellipsoid";

            // Solver Options
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.01;
            C.NoOfTimesteps = 100;


            // Re
            C.DimensionlessNumbers.Oh = 1e-1;


            C.ActiveTerms.Viscous = Viscous.SIP;
            C.ActiveTerms.Convective = Convective.Temam;


            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);
                var _zNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, BoSSS.Foundation.Grid.RefElements.CellType.Cube_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];
                    double z = _X[2];

                    if(Math.Abs(y - _yNodes.First()) < 1.0e-8)
                        // bottom
                        return "Dirichlet_bottom";
                    if(Math.Abs(y - _yNodes.Last()) < 1.0e-8)
                        // top
                        return "Dirichlet_top";
                    if(Math.Abs(x - _xNodes.First()) < 1.0e-8)
                        // left
                        return "Dirichlet_left";
                    if(Math.Abs(x - _xNodes.Last()) < 1.0e-8)
                        // right
                        return "Dirichlet_right";
                    if(Math.Abs(z - _zNodes.First()) < 1.0e-8)
                        // front
                        return "Dirichlet_front";
                    if(Math.Abs(z - _zNodes.Last()) < 1.0e-8)
                        // back
                        return "Dirichlet_back";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            List<ParameterFunctionPair> LevelSet = new List<ParameterFunctionPair>();
            LevelSet.Add(new ParameterFunctionPair(-0.25, X => 1.0, new Func<double[], double>[] { X => 0.0, X => 0.0, X => 0.0 }));
            LevelSet.Add(new ParameterFunctionPair(ecc, X => X[0] * X[0], new Func<double[], double>[] { X => 2 * X[0], X => 0.0, X => 0.0 }));
            LevelSet.Add(new ParameterFunctionPair(1.0 / ecc, X => X[1] * X[1], new Func<double[], double>[] { X => 0.0, X => 2 * X[1], X => 0.0 }));
            LevelSet.Add(new ParameterFunctionPair(1.0, X => X[2] * X[2], new Func<double[], double>[] { X => 0.0, X => 0.0, X => 2 * X[2] }));


            C.ParameterLevelSet = LevelSet.ToArray();

            //C.InitialValues.Add("Phi", new Formula($"X => {ecc} * X[0] * X[0] + 1 / {ecc} * X[1] * X[1] + X[2] * X[2] - 0.25"));

            return C;
        }
    }
}
