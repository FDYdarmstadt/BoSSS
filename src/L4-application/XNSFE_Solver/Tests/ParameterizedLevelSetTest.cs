using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.XNSFE_Solver.Tests {

    static class ParameterizedLevelSetTest_Elemental {



        public static void Test() {
            var C = Control();
            using(var solver = new XNSFE() ){
                solver.Init(C);
                solver.RunSolverMode();
            }
        }
        static XNSFE_Control Control() {
            XNSFE_Control C = new XNSFE_Control();

            //C.GridFunc = () => {
            //    var xNodes = GenericBlas.Linspace(-1.0, 1.0, 20 + 1);
            //    var yNodes = GenericBlas.Linspace(-5.0, 5.0, 100 + 1);

            //    var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: true);
            //    grd.EdgeTagNames.Add(1, "wall_ZeroGradient");
            //    grd.DefineEdgeTags((Vector X) => 1);
            //    return grd;
            //};

            double Lx = 1.0;
            double Ly = 5 * Lx;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-Lx, Lx, 20 + 1);
                double[] Ynodes = GenericBlas.Linspace(- Ly,Ly, 100 + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "NavierSlip_Linear_ZeroGradient");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_ZeroGradient");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet_ZeroGradient");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if ((Math.Abs(X[0] - Xnodes.First()) < 1e-8) || (Math.Abs(X[0] - Xnodes.Last()) < 1e-8))
                        return 1; // walls
                    else if ((Math.Abs(X[1] - Ynodes.Last()) < 1e-8))
                        return 2; // upper border
                    else if ((Math.Abs(X[1] - Ynodes.First()) < 1e-8))
                        return 3; // bottom

                    return et;
                });

                return grd;
            };

            double xSemiAxis0 = 1.5;
            double ySemiAxis0 = 1.0;
            double yCenter0 = 0.0;
            //C.AddBoundaryValue("wall_ZeroGradient");

            C.AddBoundaryValue("NavierSlip_Linear_ZeroGradient");
            C.AddBoundaryValue("Pressure_Outlet_ZeroGradient");
            C.AddBoundaryValue("Velocity_Inlet_ZeroGradient");

            //C.AddInitialValue("VelocityY#A", new Formula("X => 2.0"));
            //C.AddInitialValue("VelocityY#B", new Formula("X => 2.0"));
            double ka = 0.5;
            double kb = -0.5;
            double kc = 0.0;
            Func<double[], double> VelFunc = (X => kc - (xSemiAxis0 * ySemiAxis0 * kb * (xSemiAxis0.Pow2() - X[0].Pow2()) + ySemiAxis0.Pow2() * ka * X[0].Pow2()) / ( (ySemiAxis0.Pow2() *(1 - X[0].Pow2() / xSemiAxis0.Pow2()) ).Sqrt() * xSemiAxis0.Pow2() * xSemiAxis0));
            C.InitialValues_Evaluators.Add("VelocityY#A", VelFunc);
            C.InitialValues_Evaluators.Add("VelocityY#B", VelFunc);

            Func<double[], double> PhiFunc = (X => X[1] - yCenter0 + (ySemiAxis0.Pow2() * (1 - X[0].Pow2() / xSemiAxis0.Pow2())).Sqrt());
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);
            //C.AddInitialValue("Phi", new Formula($"X => X[0]*X[0]/({xSemiAxis0} * {xSemiAxis0}) +  (X[1] - {yCenter0}).Pow2()/({ySemiAxis0} * {ySemiAxis0}) - 1"));
            //C.AddInitialValue("Phi", new Formula($"X => X[1] - {yCenter0} + Math.Sqrt({ySemiAxis0} * {ySemiAxis0}* (1 - X[0].Pow2() / ({xSemiAxis0} * {xSemiAxis0} ) ))"));

            C.SkipSolveAndEvaluateResidual = true;


            C.SetDGdegree(2);

            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.LieSplitting;
            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.ParameterizedLevelSet;

            if (C.Option_LevelSetEvolution == LevelSetEvolution.ParameterizedLevelSet) {
                C.ParameterizedLevelSetControl = new BoSSS.Solution.LevelSetTools.ParameterizedLevelSet.ParameterizedLevelSetControlEllipse(xSemiAxis0, ySemiAxis0, yCenter0);
            }

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.01;
            C.NoOfTimesteps = 50;
            C.Endtime = 1.0;

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            return C; 
        }
    }

    public static class ParameterizedLevelSet_Translation {

        public static XNSFE_Control Translation() {
            XNSFE_Control C = new XNSFE_Control();

            double Lx = 1.0;
            double Ly = 5 * Lx;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-Lx, Lx, 20 + 1);
                double[] Ynodes = GenericBlas.Linspace(-Ly, Ly, 100 + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "NavierSlip_Linear_ZeroGradient");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_ZeroGradient");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet_ZeroGradient");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if ((Math.Abs(X[0] - Xnodes.First()) < 1e-8) || (Math.Abs(X[0] - Xnodes.Last()) < 1e-8))
                        return 1; // walls
                    else if ((Math.Abs(X[1] - Ynodes.Last()) < 1e-8))
                        return 2; // upper border
                    else if ((Math.Abs(X[1] - Ynodes.First()) < 1e-8))
                        return 3; // bottom

                    return et;
                });

                return grd;
            };

            double xSemiAxis0 = 1.5;
            double ySemiAxis0 = 1.0;
            double yCenter0 = 0.0;

            C.AddBoundaryValue("NavierSlip_Linear_ZeroGradient");
            C.AddBoundaryValue("Pressure_Outlet_ZeroGradient");
            C.AddBoundaryValue("Velocity_Inlet_ZeroGradient");

            double ka = 0.0;
            double kb = 0.0;
            double kc = 1.0;
            Func<double[], double> VelFunc = (X => kc - (xSemiAxis0 * ySemiAxis0 * kb * (xSemiAxis0.Pow2() - X[0].Pow2()) + ySemiAxis0.Pow2() * ka * X[0].Pow2()) / ((ySemiAxis0.Pow2() * (1 - X[0].Pow2() / xSemiAxis0.Pow2())).Sqrt() * xSemiAxis0.Pow2() * xSemiAxis0));
            C.InitialValues_Evaluators.Add("VelocityY#A", VelFunc);
            C.InitialValues_Evaluators.Add("VelocityY#B", VelFunc);

            Func<double[], double> PhiFunc = (X => X[1] - yCenter0 + (ySemiAxis0.Pow2() * (1 - X[0].Pow2() / xSemiAxis0.Pow2())).Sqrt());
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.SkipSolveAndEvaluateResidual = true;


            C.SetDGdegree(2);

            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.LieSplitting;
            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.ParameterizedLevelSet;

            if (C.Option_LevelSetEvolution == LevelSetEvolution.ParameterizedLevelSet) {
                C.ParameterizedLevelSetControl = new BoSSS.Solution.LevelSetTools.ParameterizedLevelSet.ParameterizedLevelSetControlEllipse(xSemiAxis0, ySemiAxis0, yCenter0);
            }

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.01;
            C.NoOfTimesteps = 5;
            C.Endtime = 1.0;

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            return C;
        }
    }

    public static class ParameterizedLevelSet_ShapeChange {

        public static XNSFE_Control ShapeChange() {
            XNSFE_Control C = new XNSFE_Control();

            double Lx = 1.0;
            double Ly = 5 * Lx;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-Lx, Lx, 20 + 1);
                double[] Ynodes = GenericBlas.Linspace(-Ly, Ly, 100 + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "NavierSlip_Linear_ZeroGradient");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_ZeroGradient");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet_ZeroGradient");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if ((Math.Abs(X[0] - Xnodes.First()) < 1e-8) || (Math.Abs(X[0] - Xnodes.Last()) < 1e-8))
                        return 1; // walls
                    else if ((Math.Abs(X[1] - Ynodes.Last()) < 1e-8))
                        return 2; // upper border
                    else if ((Math.Abs(X[1] - Ynodes.First()) < 1e-8))
                        return 3; // bottom

                    return et;
                });

                return grd;
            };

            double xSemiAxis0 = 1.5;
            double ySemiAxis0 = 1.0;
            double yCenter0 = 0.0;

            C.AddBoundaryValue("NavierSlip_Linear_ZeroGradient");
            C.AddBoundaryValue("Pressure_Outlet_ZeroGradient");
            C.AddBoundaryValue("Velocity_Inlet_ZeroGradient");

            double ka = 0.0;
            double kb = -0.5;
            double kc = 0.0;
            Func<double[], double> VelFunc = (X => kc - (xSemiAxis0 * ySemiAxis0 * kb * (xSemiAxis0.Pow2() - X[0].Pow2()) + ySemiAxis0.Pow2() * ka * X[0].Pow2()) / ((ySemiAxis0.Pow2() * (1 - X[0].Pow2() / xSemiAxis0.Pow2())).Sqrt() * xSemiAxis0.Pow2() * xSemiAxis0));
            C.InitialValues_Evaluators.Add("VelocityY#A", VelFunc);
            C.InitialValues_Evaluators.Add("VelocityY#B", VelFunc);

            Func<double[], double> PhiFunc = (X => X[1] - yCenter0 + (ySemiAxis0.Pow2() * (1 - X[0].Pow2() / xSemiAxis0.Pow2())).Sqrt());
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.SkipSolveAndEvaluateResidual = true;


            C.SetDGdegree(2);

            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.LieSplitting;
            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.ParameterizedLevelSet;

            if (C.Option_LevelSetEvolution == LevelSetEvolution.ParameterizedLevelSet) {
                C.ParameterizedLevelSetControl = new BoSSS.Solution.LevelSetTools.ParameterizedLevelSet.ParameterizedLevelSetControlEllipse(xSemiAxis0, ySemiAxis0, yCenter0);
            }

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.01;
            C.NoOfTimesteps = 1;
            C.Endtime = 1.0;

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            return C;
        }
    }

    class ParameterizedLevelSetTest : IXNSFETest {

        public ParameterizedLevelSetTest() {

        }

        public bool TestImmersedBoundary => false;

        public Func<double[], double, double> GetPhi2() {
            throw new NotImplementedException(); // will never be called, as long as 'TestImmersedBoundary' == false;
        }

        public Func<double[], double, double> GetPhi2U(int d) {
            throw new NotImplementedException();
        }

        public int SpatialDimension => 2;

        const double DomainLength = 1e-04;
        const double Twall = 37.5;


        public double xSemiAxis0 => (2 * DomainLength / (3.0).Sqrt());

        public double ySemiAxis0 => (2 * DomainLength / (3.0).Sqrt());
            
        public double yCenter0 => 3e-04;

        public double k_x0 => 0.0;
       
        public double k_y0 => 0.0;

        public double k_yc0 => 1.0;

        public Func<double[], double, double> GetPhi() {
            //return (X, t) => X[1] - (k_yc0 * t + yCenter0) + (k_y0 * t + ySemiAxis0) * (1 - (X[0] / (k_x0 * t + xSemiAxis0)).Pow2()).Sqrt();
            return (X, t) => X[1] - (k_yc0 * t + yCenter0) +  ((k_y0 * t + ySemiAxis0) * (1 - (X[0] / (k_x0 * t + xSemiAxis0)).Pow2())).Sqrt();
        }

        public GridCommons CreateGrid(int Resolution) {

            if (Resolution < 1)
                throw new ArgumentException();

            double[] Xnodes = GenericBlas.Linspace(-DomainLength, DomainLength, 2 * Resolution + 1); 
            double[] Ynodes = GenericBlas.Linspace(0, 5 * DomainLength, 5 * Resolution + 1); 

            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes); // creation of 2D grid 


            grd.EdgeTagNames.Add(1, "NavierSlip_Linear_ZeroGradient");
            grd.EdgeTagNames.Add(2, "Pressure_Outlet_ZeroGradient");
            grd.EdgeTagNames.Add(3, "Velocity_Inlet_ConstantTemperature");

            grd.DefineEdgeTags(delegate (double[] X) {
                byte et = 0;
                if ((Math.Abs(X[0] - Xnodes.First()) < 1e-8) || (Math.Abs(X[0] - Xnodes.Last()) < 1e-8))
                    return 1; // walls
                else if ((Math.Abs(X[1] - Ynodes.Last()) < 1e-8))
                    return 2; // upper border
                else if ((Math.Abs(X[1] - Ynodes.First()) < 1e-8))
                    return 3; // bottom

                return et;
            });

            return grd;
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("NavierSlip_Linear_ZeroGradient", new AppControl.BoundaryValueCollection());
            config.Add("Pressure_Outlet_ZeroGradient", new AppControl.BoundaryValueCollection());
            config.Add("Velocity_Inlet_ConstantTemperature", new AppControl.BoundaryValueCollection());

            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add("Temperature#A", (X, t) => Twall);
            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add("VelocityX#A", (X, t) => 0.0);
            // config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add("VelocityY#A", (X, t) => k_yc0 - (ySemiAxis0 * k_y0 * (1 - X[0].Pow2() / xSemiAxis0.Pow2()) + ySemiAxis0.Pow2() * X[0].Pow2() * k_x0 / (xSemiAxis0.Pow2() * xSemiAxis0)) / Math.Sqrt(ySemiAxis0.Pow2()* (1 - X[0].Pow2() / xSemiAxis0.Pow2())));
            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add("VelocityY#A", (X, t) => 1 - 1 * X[0].Pow2() / DomainLength.Pow2()) ;
            //config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add("VelocityY#A", (X, t) => 1) ;

            return config;
        }

        public Func<double, double> GetE() => null;

        public Func<double[], double> GetF(string species, int d) {
            return (X => 0.0);
        }

        public Func<double[], double> GetQ(string species) {
            return (X => 0.0);
        }

        public Func<double[], double, double> GetU(string species, int d) {
            switch (species) {
                case "A": { return (X, t) => d == 0 ? 0.0 : d == 1 ?  1 - 1 * X[0].Pow2() / DomainLength.Pow2() : throw new ArgumentException(); } 
                case "B": { return (X, t) => 0.0; }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetPress(string species) {
            switch (species) {
                case "A": { return (X, t) => 0.0; }
                case "B": { return (X, t) => 0.0; }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetT(string species) {
            switch (species) {
                case "A": { return (X, t) => Twall; }
                case "B": { return (X, t) => 0.0; }
                default: { throw new ArgumentException(); }
            }
        }


        public double rho_A => 600.0;

        public double rho_B => 9.0;
  
        public double mu_A => 0.13;

        public double mu_B => 0.01019;

        public double Sigma => 0.02;

        public double c_A => 100;

        public double c_B => 37;

        public double k_A => 4.8;

        public double k_B => 0.0251;
    
        public double h_vap => 10000;
  
        public double T_sat => 0;

        public double theta_e => Math.PI / 6.0;

        public bool Material => false;

        public bool steady => false;

        public double dt => 1e-2;


        public bool IncludeConvection => true;

        public double[] AcceptableL2Error =>  new double[] { 1.0e-5, 1.0e-5, 0.02, 1.0e-5 };

        public double[] AcceptableResidual => new double[] { 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5 };

        public int LevelsetPolynomialDegree => 2;

        public bool CheckT => true;

        public bool CheckE => false;
    }
}