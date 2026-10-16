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


namespace StokesHelical_Ak {
    public static class HardcodedControl {
        #region Steady         
        //===========
        public static HelicalControl ManSol_Steady_DDD_Paper(string _DbPath = @"", int degree = 3, int noOfCellsR = 32, int noOfCellsXi = 32, bool _periodicXi = true, string _prjnm = "Helical", double rMin = 0.1) {
            // DDD Paper Equation: 4.4 (a,b,c,d)
            // Convergance is confirmed!
            // Test with and without R0_fix
            HelicalControl Ctrl = new HelicalControl();

            #region db
            Ctrl.DbPath = _DbPath;
            //Ctrl.DbPath = @"P:\BoSSSpostprocessing\Akbari";
            Ctrl.DbPath = null;
            Ctrl.savetodb = Ctrl.DbPath != null;
            Ctrl.ProjectName = _prjnm;
            Ctrl.SessionName = "degree = " + degree + " " + "noOfCellsR= " + noOfCellsR + " " + "noOfCellsXi= " + noOfCellsXi;
            #endregion
            //Ctrl.savetodb = true;
            #region Grid
            Ctrl.Resolution_R = noOfCellsR;
            Ctrl.Resolution_Xi = noOfCellsXi;
            Ctrl.rMin = rMin;
            Ctrl.rMax = 1;
            Ctrl.GridFunc = delegate {
                double[] xnodes = GenericBlas.Linspace(rMin, Ctrl.rMax, noOfCellsR + 1);
                double[] ynodes = GenericBlas.Linspace(0, 2 * Math.PI, noOfCellsXi + 1);
                GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, periodicY: true);

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Stuff");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double r, xi;
                    r = X[0];
                    xi = X[1];
                    if(Math.Abs(r - Ctrl.rMax) < 1E-8 || (Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8))
                        // Question: Should I include Ctrl.rMin >= 1E-6 ?? R0_fix is ruling.
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
            double dt = 1E20;
            Ctrl.dtMax = dt;
            Ctrl.dtMin = dt;
            Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            Ctrl.NoOfTimesteps = 1;
            //###########
            // Need this to distinguish in the forcing terms
            Ctrl.steady = true;
            Globals.steady = Ctrl.steady;
            // Need this to distinguish in the forcing terms
            //###########
            #endregion
            // DG degree
            // =========
            Ctrl.dg_degree = degree;
            Ctrl.SetDGdegree(degree);
            // Initial Values
            // ==============
            #region Initial Values
            double a = Globals.a;
            double b = Globals.b;
            Ctrl.AddInitialValue("Pressure", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("ur", new Formula("(X) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("ueta", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("uxi", new Formula("(X) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(0)"));
            #endregion
            // Boundary Conditions
            // ==============
            #region Boundary Conditions
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(t)", true));
            #endregion
            return Ctrl;
        }


        public static HelicalControl ManSol_Steady_DDD_Hanging_Nodes(string _DbPath = @"", int degree = 3, int noOfCellsR = 32, int noOfCellsXi = 32, bool _periodicXi = true, string _prjnm = "Helical", double rMin = 0.1) {
            // DDD Paper Equation: 4.4 (a,b,c,d)
            // Convergance is confirmed!
            // Test with and without R0_fix
            HelicalControl Ctrl = new HelicalControl();

            #region db
            Ctrl.DbPath = _DbPath;
            //Ctrl.DbPath = @"P:\BoSSSpostprocessing\Akbari";
            Ctrl.DbPath = null;
            Ctrl.savetodb = Ctrl.DbPath != null;
            Ctrl.ProjectName = _prjnm;
            Ctrl.SessionName = "degree = " + degree + " " + "noOfCellsR= " + noOfCellsR + " " + "noOfCellsXi= " + noOfCellsXi;
            #endregion
            //Ctrl.savetodb = true;
            #region Grid
            Ctrl.Resolution_R = noOfCellsR;
            Ctrl.Resolution_Xi = noOfCellsXi;
            Ctrl.rMin = rMin;
            Ctrl.rMax = 1;

            Ctrl.GridFunc = delegate {
                int numberOfSegments = 4;
                // Initialisierung einer Liste für die Gitter
                GridCommons[] grids = new GridCommons[numberOfSegments];
                // Basisanzahl der Zellen im ersten Segment
                int N = (noOfCellsR + 1) / numberOfSegments;
                // Startpunkt für xNodes
                double xStart = rMin;
                // Die Schrittweite für xNodes-Bereiche erhöhen
                double xStep = (Ctrl.rMax / numberOfSegments) - xStart;

                for(int i = 0; i < numberOfSegments; i++) {
                    // Berechnung der Knotenanzahl für das aktuelle Segment
                    int nodesInSegment = N * (int)Math.Pow(2, i) + 1;

                    // Bestimmung der x-Bereichsgrenzen für das aktuelle Segment
                    double xEnd = xStart + xStep;

                    // Erzeugung der x-Knoten und y-Knoten für das aktuelle Segment
                    double[] xNodes = GenericBlas.Linspace(xStart, xEnd, nodesInSegment);
                    double[] yNodes = GenericBlas.Linspace(0, 2 * Math.PI, nodesInSegment);

                    // Erstellung des Gitters für das aktuelle Segment und Hinzufügen in die Liste
                    grids[i] = Grid2D.Cartesian2DGrid(xNodes, yNodes, type: CellType.Square_Linear);
                    // Aktualisierung des Startpunkts für den nächsten Durchlauf
                    xStart = xEnd;
                }

                // Zusammenführen der Gitter aus der Liste
                var grdJ = grids[0];
                for(int i = 1; i < grids.Length; i++) {
                    grdJ = GridCommons.MergeLogically(grids[i], grdJ);
                }
                // Versiegelung des zusammengeführten Gitters
                GridCommons grd = GridCommons.Seal(grdJ, 4);

                // Hinzufügen der Edge Tag Names und Definition der Edge Tags...
                // (Füge hier den restlichen Teil deines Codes ein)

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Stuff");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double r, xi;
                    r = X[0];
                    xi = X[1];
                    if(Math.Abs(r - Ctrl.rMax) < 1E-8 || (Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8))
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
            double dt = 1E20;
            Ctrl.dtMax = dt;
            Ctrl.dtMin = dt;
            Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            Ctrl.NoOfTimesteps = 1;
            //###########
            // Need this to distinguish in the forcing terms
            Ctrl.steady = true;
            Globals.steady = Ctrl.steady;
            // Need this to distinguish in the forcing terms
            //###########
            #endregion
            // DG degree
            // =========
            Ctrl.dg_degree = degree;
            Ctrl.SetDGdegree(degree);
            // Initial Values
            // ==============
            #region Initial Values
            double a = Globals.a;
            double b = Globals.b;
            Ctrl.AddInitialValue("Pressure", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("ur", new Formula("(X) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("ueta", new Formula("(X) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(0)"));
            Ctrl.AddInitialValue("uxi", new Formula("(X) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(0)"));
            #endregion
            // Boundary Conditions
            // ==============
            #region Boundary Conditions
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(t)", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) => (0.2e1 * X[0] * X[0] * Math.Pow(" + a * a + " * X[0] * X[0] +" + b * b + ", -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(" + a * a + " * X[0] * X[0] + " + b * b + ", -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(t)", true));
            #endregion
            return Ctrl;
        }
        #endregion

        #region transient         
        // DDD Paper Equation: 4.4 (a,b,c,d)
        // Convergance is confirmed!
        // Test with and without R0_fix
        public static HelicalControl ManSol_Transient_DDD_Paper(string _DbPath = @"", int degree = 3, int noOfCellsR = 32, int noOfCellsXi = 32, int dtRefining = 32, string bdfOrder = "BDF1", double rMin = 0) {

            HelicalControl Ctrl = new HelicalControl();
            #region db

            //Ctrl.DbPath = @"P:\BoSSSpostprocessing\Akbari"; // _DbPath;

            Ctrl.DbPath = null;
            Ctrl.savetodb = Ctrl.DbPath != null;
            Ctrl.ProjectName = "NStransient";
            Ctrl.SessionName = "degree= " + degree + " " + "noOfCellsR= " + noOfCellsR + " " + "noOfCellsXi= " + noOfCellsXi;
            Ctrl.rMin = rMin;
            Ctrl.rMax = 1;
            #endregion

            //Ctrl.savetodb = true;
            Ctrl.Resolution_R = noOfCellsR;
            Ctrl.Resolution_Xi = noOfCellsXi;

            Ctrl.GridFunc = delegate {

                double[] xnodes = GenericBlas.Linspace(rMin, Ctrl.rMax, noOfCellsR + 1);
                double[] ynodes = GenericBlas.Linspace(0, 2 * Math.PI, noOfCellsXi + 1);
                GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, periodicY: true);

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Stuff");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double r, xi;
                    r = X[0];
                    xi = X[1];
                    if(Math.Abs(r - Ctrl.rMax) < 1E-8 || (Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8))
                        return 1;
                    else
                        return 2;
                });

                return grd;
            };

            double dt = 2 * Math.PI / dtRefining;

            Ctrl.dtMax = dt;
            Ctrl.dtMin = dt;
            if(bdfOrder == "BDF3") {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            } else if(bdfOrder == "BDF1") {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            } else {
                throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
            }
            Ctrl.NoOfTimesteps = dtRefining;
            Ctrl.steady = false;


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

        #endregion
        public static HelicalControl HelicalControlRestart(string SID, int ts) {
            // @"\\dc3\userspace\akbari\cluster\Helical_DNS\1cb79396-cda0-49b7-851f-5c7fea5e8731"
            //StokesHelical_Ak.DNS.SimpleFlow(degree: 2, noOfCellsR: 8, noOfCellsXi: 8, rMin: 0, dtRefining: 32);
            HelicalControl C = StokesHelical_Ak.DNS.SimpleFlow();
            //var C = ManSol_Transient_DDD_Paper(degree: 4, noOfCellsR: 32, noOfCellsXi: 32, rMin: 0);
            C.InitialValues.Clear();
            C.InitialValues_Evaluators.Clear();
            C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid(SID), new TimestepNumber(ts));
            // C.GridGuid = tsi.GridID;
            C.GridFunc = null;
            return C;
        }

        public static void TSFP() {
            int[] degrees = { 2 };
            //int[] degrees = {5}; 
            int[] noOfCellsRs = { 5 };
            //int[] noOfCellsRs = {256};
            double rMin_ = 0;
            int[] dtRefinings = {5};
            //int[] dtRefinings = { 1 };
            string[] bdfOrders = { "BDF1" };
            //string[] bdfOrders = {"BDF3"};
            //double[] MaxAmps = {0,1,0.01,0.001};
            double[] MaxAmps = { 5 };
            //string[] grids = {"regular","hangingNodes"};
            string[] grids = { "regular" };

            foreach(string grid in grids) {
                foreach(double MaxAmp in MaxAmps) {
                    foreach(int degree in degrees) {
                        foreach(int noOfCellsR in noOfCellsRs) {
                            foreach(int dtRefining in dtRefinings) {
                                foreach(string bdfOrder in bdfOrders) {
                                    HelicalControl Ctrl = new HelicalControl();

                                    int noOfCellsXi = noOfCellsR; // Aquidistnant
                                    Ctrl.rMin = rMin_;
                                    Ctrl.rMax = 1;
                                    Ctrl.grid = grid;
                                    Ctrl.Resolution_R = noOfCellsR;
                                    Ctrl.Resolution_Xi = noOfCellsXi;
                                    if(grid == "regular") {
                                        Ctrl.GridFunc = delegate {
                                            double[] xnodes = GenericBlas.Linspace(Ctrl.rMin, Ctrl.rMax, noOfCellsR + 1);
                                            double[] ynodes = GenericBlas.Linspace(0, 2 * Math.PI, noOfCellsXi + 1);
                                            GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, periodicY: true);

                                            grd.EdgeTagNames.Add(1, "Dirichlet_rmax");
                                            grd.EdgeTagNames.Add(2, "Dirichlet_rmin");
                                            grd.EdgeTagNames.Add(3, "Stuff");

                                            grd.DefineEdgeTags(delegate (double[] _X) {
                                                var X = _X;
                                                double r, xi;
                                                r = X[0];
                                                xi = X[1];
                                                if(Math.Abs(r - Ctrl.rMax) < 1E-8)
                                                    return 1;
                                                else if(Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8)
                                                    return 2;
                                                else
                                                    return 3;
                                            });

                                            return grd;
                                        };
                                    } else {
                                        Ctrl.GridFunc = delegate {
                                            int numberOfSegments = noOfCellsR / 16;
                                            // Initialisierung einer Liste für die Gitter
                                            GridCommons[] grids = new GridCommons[numberOfSegments];
                                            // Basisanzahl der Zellen im ersten Segment
                                            int N = (noOfCellsR + 1) / numberOfSegments;
                                            // Startpunkt für xNodes
                                            double xStart = rMin_;
                                            // Die Schrittweite für xNodes-Bereiche erhöhen
                                            double xStep = (Ctrl.rMax / numberOfSegments) - xStart;

                                            for(int i = 0; i < numberOfSegments; i++) {
                                                // Berechnung der Knotenanzahl für das aktuelle Segment
                                                int nodesInSegment = N * (int)Math.Pow(2, i) + 1;

                                                // Bestimmung der x-Bereichsgrenzen für das aktuelle Segment
                                                double xEnd = xStart + xStep;

                                                // Erzeugung der x-Knoten und y-Knoten für das aktuelle Segment
                                                double[] xNodes = GenericBlas.Linspace(xStart, xEnd, nodesInSegment);
                                                double[] yNodes = GenericBlas.Linspace(0, 2 * Math.PI, nodesInSegment);

                                                // Erstellung des Gitters für das aktuelle Segment und Hinzufügen in die Liste
                                                grids[i] = Grid2D.Cartesian2DGrid(xNodes, yNodes, type: CellType.Square_Linear);
                                                // Aktualisierung des Startpunkts für den nächsten Durchlauf
                                                xStart = xEnd;
                                            };

                                            // Zusammenführen der Gitter aus der Liste
                                            var grdJ = grids[0];
                                            for(int i = 1; i < grids.Length; i++) {
                                                grdJ = GridCommons.MergeLogically(grids[i], grdJ);
                                            }
                                            // Versiegelung des zusammengeführten Gitters
                                            GridCommons grd = GridCommons.Seal(grdJ, 4);

                                            // Hinzufügen der Edge Tag Names und Definition der Edge Tags...
                                            // (Füge hier den restlichen Teil deines Codes ein)

                                            grd.EdgeTagNames.Add(1, "Dirichlet");
                                            grd.EdgeTagNames.Add(2, "Stuff");

                                            grd.DefineEdgeTags(delegate (double[] _X) {
                                                var X = _X;
                                                double r, xi;
                                                r = X[0];
                                                xi = X[1];
                                                if(Math.Abs(r - Ctrl.rMax) < 1E-8 || (Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8))
                                                    return 1;
                                                else
                                                    return 2;
                                            });
                                            return grd;
                                        };
                                    }
                                    double dt = 2 * Math.PI / (dtRefining * 100);
                                    //double dt = 1E20;
                                    Ctrl.dtMax = dt;
                                    Ctrl.dtMin = dt;
                                    if(bdfOrder == "BDF3") {
                                        Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
                                    } else if(bdfOrder == "BDF1") {
                                        Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
                                    } else {
                                        throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
                                    }
                                    //Ctrl.NoOfTimesteps = 1* dtRefining*10000;
                                    Ctrl.NoOfTimesteps = 1;

                                    // Solver Properties
                                    //=============
                                    Ctrl.steady = false;
                                    Ctrl.ExactResidual = false;
                                    Ctrl.HagenPoisseulle = true;
                                    // DG degree
                                    // =========
                                    Ctrl.dg_degree = degree;
                                    Ctrl.SetDGdegree(degree);

                                    double a = Globals.a;
                                    double b = Globals.b;
                                    double nu = Globals.nu;

                                    // Initial Values
                                    // ==============
                                    Ctrl.AddInitialValue("Pressure", new Formula("(X) =>0"));
                                    Ctrl.AddInitialValue("ur", new Formula("(X) => 0"));
                                    Ctrl.AddInitialValue("ueta", new Formula($"(X) => {MaxAmp} *(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} ))) * ({a * b} * ({Ctrl.rMax * Ctrl.rMax} - X[0]*X[0]) )/(X[0]*4* {nu})"));
                                    Ctrl.AddInitialValue("uxi", new Formula($"(X) => -{MaxAmp} *(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} ))) * ({a * a} * ({Ctrl.rMax * Ctrl.rMax} - X[0]*X[0]) )/(4* {nu})"));

                                    //Ctrl.AddInitialValue("Pressure", new Formula("(X) =>0"));
                                    //Ctrl.AddInitialValue("ur", new Formula("(X) => 0"));
                                    //Ctrl.AddInitialValue("ueta", new Formula($"(X) => 0"));
                                    //Ctrl.AddInitialValue("uxi", new Formula($"(X) => 0"));
                                    // Boundary Conditions
                                    // ==============
                                    Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) => 0", true));
                                    Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) =>0", true));
                                    Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) =>0", true));
                                    Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) =>0", true));

                                    // Solver Propertiees 
                                    // ==============
                                    Ctrl.maxAmpli = MaxAmp;

                                    //Ctrl.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { ConvergenceCriterion = 1e-13 , CoarseKickIn = 1000000};
                                    //Ctrl.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { ConvergenceCriterion = 1e-13 , TargetBlockSize =1000000};

                                    // Plotting
                                    // ==============
                                    //Ctrl.ImmediatePlotPeriod = 1;
                                    Ctrl.ImmediatePlotPeriod = 1;
                                    Ctrl.SuperSampling = 1;
                                    Ctrl.TracingNamespaces = "*";


                                    string sessionName = $"Steady_Stokes_{Ctrl.TimeSteppingScheme}_-Amp={Ctrl.maxAmpli}_Grid{Ctrl.grid}_degree={Ctrl.dg_degree}_noOfCellsR={Ctrl.Resolution_R}_noOfCellsXi={Ctrl.Resolution_Xi}_{DateTime.Now.ToString("MMMdd_HHmm")}";
                                    //string sessionName = $"Stationary_Stokes_{Ctrl.TimeSteppingScheme}_-Amp={Ctrl.maxAmpli}_Grid{Ctrl.grid}_degree={Ctrl.dg_degree}_noOfCellsR={Ctrl.Resolution_R}_noOfCellsXi={Ctrl.Resolution_Xi}_{DateTime.Now.ToString("MMMdd_HHmm")}";
                                    string projectName = "Stability_Check";

                                    var myBatch = BoSSS.Application.BoSSSpad.BoSSSshell.ExecutionQueues[2];
                                    BoSSS.Application.BoSSSpad.BoSSSshell.WorkflowMgm.SetNameBasedSessionJobControlCorrelation();
                                    // var dB = BoSSS.Application.BoSSSpad.BoSSSshell.WorkflowMgm.DefaultDatabase;
                                    Ctrl.SessionName = sessionName;
                                    Ctrl.ProjectName = projectName;

                                     string dB = "\\\\dc3\\userspace\\akbari\\Local\\Stability_Check";

                                    //string restartDBfullPath = Path.Combine(Directory.GetCurrentDirectory(), dB);
                                    //string restartDBfullPath = Path.Combine(Directory.GetCurrentDirectory(), dB);
                                    // Database
                                    // ==============
                                    //DatabaseUtils.CreateDatabase(dB);
                                    Ctrl.DbPath = dB; // _DbPath;
                                    Ctrl.savetodb = Ctrl.DbPath != null;
                                    // Solver Init
                                    // ==============
                                    var solver = new HelicalMain();
                                    solver.Init(Ctrl);
                                    solver.RunSolverMode();
                                }
                            }
                        }
                    }
                }
            }
        }


        static public void HagenPoiseulle_Check_Conv_Terms(int pOrder = 2, bool NavierStokes = true,  int BDForder = 1) {
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");   
            //ilPSP.Environment.NumThreads = 1;
            var Ctrl = StokesHelical_Ak.DNS_Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 5, noOfCellsXi: 5, dtRefining: 1, Tend: 5, _DbPath: tempDB.Path, bdfOrder: 1, rMin:0.1);

            Ctrl.ImmediatePlotPeriod = -1;
            Ctrl.NavierStokes = NavierStokes;
            Ctrl.InitialValues.Clear();
            Ctrl.InitialValues_Evaluators.Clear();

            Ctrl.AddInitialValue("Pressure", new Formula("(X) =>0"));
            Ctrl.AddInitialValue("ur", new Formula("(X) => X[0]"));
            Ctrl.AddInitialValue("ueta", new Formula($"(X) => X[0]"));
            Ctrl.AddInitialValue("uxi", new Formula($"(X) => 1"));

            Ctrl.savetodb = true;
            Ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;

            using(var solverStat = new HelicalMain()) {
                solverStat.Init(Ctrl);
                solverStat.RunSolverMode();
            }

        }
    }

}
