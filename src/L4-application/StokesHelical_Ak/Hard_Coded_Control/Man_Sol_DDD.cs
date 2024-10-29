using ilPSP.Utils;
using System;
using System.IO;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.IO;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Application.XNSE_Solver;


namespace StokesHelical_Ak.Hard_Coded_Control
{
    public static class Man_Sol_DDD
    {
        //#region Steady         
        ////===========
        //public static HelicalControl ManSol_Steady_DDD_Paper(int degree = 3, int noOfCellsR = 32, int noOfCellsXi = 32, bool _periodicXi = true, string _prjnm = "Helical", double rMin = 0.1)
        //{
        //    // DDD Paper Equation: 4.4 (a,b,c,d)
        //    // Convergance is confirmed!
        //    // Test with and without R0_fix
        //    HelicalControl Ctrl = new HelicalControl();
        //    #region db
        //    Ctrl.DbPath = null;
        //    Ctrl.savetodb = Ctrl.DbPath != null;
        //    #endregion
        //    #region Grid
        //    Ctrl.Resolution_R = noOfCellsR;
        //    Ctrl.Resolution_Xi = noOfCellsXi;
        //    Ctrl.rMin = rMin;
        //    Ctrl.rMax = 1;
        //    Ctrl.GridFunc = delegate
        //    {
        //        double[] xnodes = GenericBlas.Linspace(rMin, Ctrl.rMax, noOfCellsR + 1);
        //        double[] ynodes = GenericBlas.Linspace(0, 2 * Math.PI, noOfCellsXi + 1);
        //        GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, periodicY: true);

        //        grd.EdgeTagNames.Add(1, "Dirichlet");
        //        grd.EdgeTagNames.Add(2, "Stuff");

        //        grd.DefineEdgeTags(delegate (double[] _X)
        //        {
        //            var X = _X;
        //            double r, xi;
        //            r = X[0];
        //            xi = X[1];
        //            if (Math.Abs(r - Ctrl.rMax) < 1E-8 || Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8)
        //                return 1;
        //            else
        //                return 2;
        //        });
        //        return grd;
        //    };
        //    #endregion
        //    // ===================================================== //
        //    // for Timestepping
        //    #region Timestepping
        //    double dt = 1E20;
        //    Ctrl.dtMax = dt;
        //    Ctrl.dtMin = dt;
        //    Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
        //    Ctrl.NoOfTimesteps = 1;
        //    //###########
        //    // Need this to distinguish in the forcing terms
        //    Ctrl.steady = true;
        //    Globals.steady = Ctrl.steady;
        //    // Need this to distinguish in the forcing terms
        //    //###########
        //    #endregion
        //    // DG degree
        //    // =========
        //    Ctrl.dg_degree = degree;
        //    Ctrl.SetDGdegree(degree);
        //    // Initial Values
        //    // ==============
        //    #region Initial Values
        //    double a = Globals.a;
        //    double b = Globals.b;
        //    Ctrl.AddInitialValue("Pressure", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
        //    Ctrl.AddInitialValue("ur", new Formula("(X) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
        //    Ctrl.AddInitialValue("ueta", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(0)"));
        //    Ctrl.AddInitialValue("uxi", new Formula("(X) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(0)"));
        //    #endregion
        //    // Boundary Conditions
        //    // ==============
        //    #region Boundary Conditions
        //    Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
        //    Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
        //    Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(t)", true));
        //    Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(t)", true));
        //    #endregion
        //    return Ctrl;
        //}
        //#endregion

        // #region transient         
        // DDD Paper Equation: 4.4 (a,b,c,d)
        // Convergance is confirmed!
        // Test with and without R0_fix
        public static HelicalControl ManSol_DDD_Paper(int degree = 3, int noOfCellsR = 32, int noOfCellsXi = 32, int dtRefining = 32, string bdfOrder = "BDF1", double rMin = 0, bool steady = false)
        {

            HelicalControl Ctrl = new HelicalControl();
            #region db

            Ctrl.DbPath = null;
            Ctrl.savetodb = Ctrl.DbPath != null;
            Ctrl.ProjectName = "NStransient";
            Ctrl.SessionName = "degree= " + degree + " " + "noOfCellsR= " + noOfCellsR + " " + "noOfCellsXi= " + noOfCellsXi;
            Ctrl.rMin = rMin;
            Ctrl.rMax = 1;
            #endregion
            // Grid
            // ==============
            #region Grid
            Ctrl.Resolution_R = noOfCellsR;
            Ctrl.Resolution_Xi = noOfCellsXi;

            Ctrl.GridFunc = delegate
            {

                double[] xnodes = GenericBlas.Linspace(rMin, Ctrl.rMax, noOfCellsR + 1);
                double[] ynodes = GenericBlas.Linspace(0, 2 * Math.PI, noOfCellsXi + 1);
                GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, periodicY: true);

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Stuff");

                grd.DefineEdgeTags(delegate (double[] _X)
                {
                    var X = _X;
                    double r, xi;
                    r = X[0];
                    xi = X[1];
                    if (Math.Abs(r - Ctrl.rMax) < 1E-8 || Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8)
                        return 1;
                    else
                        return 2;
                });

                return grd;
            };
            #endregion
            // ===================================================== //
            // for Timestepping
            #region Timestepping
            if (steady == true) {
                double dt = 1E20;
                Ctrl.dtMax = dt;
                Ctrl.dtMin = dt;
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
                Ctrl.NoOfTimesteps = 1;
                //###########
                // Need this to distinguish in the forcing terms
                // And Boundary Conditions
                Ctrl.steady = true;
                Globals.steady = Ctrl.steady;
            } else {
                Ctrl.NoOfTimesteps = dtRefining;
                Ctrl.steady = false;
                double dt = 2 * Math.PI / dtRefining;
                Ctrl.dtMax = dt;
                Ctrl.dtMin = dt;
                if (bdfOrder == "BDF3") {
                    Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
                } else if (bdfOrder == "BDF1") {
                    Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
                }
            }
            #endregion

            // DG degree
            // =========
            Ctrl.dg_degree = degree;
            Ctrl.SetDGdegree(degree);

            // Initial Values
            // ==============
            double a = Globals.a;
            double b = Globals.b;
            #region Initial Values
            Ctrl.AddInitialValue("Pressure", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("ur", new Formula("(X) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("ueta", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("uxi", new Formula("(X) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(0)"));
            #endregion
            #region Boundary Values
            // Boundary Conditions
            // ==============
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(t)", true));
            #endregion
            return Ctrl;
        }
        // #endregion
    }
}
