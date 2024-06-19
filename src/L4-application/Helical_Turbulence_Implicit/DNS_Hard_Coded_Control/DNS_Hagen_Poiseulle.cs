using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.IO;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.Utils;
using MathNet.Numerics;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using MPI.Wrappers;
using NUnit.Framework;
using NUnitLite;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Security.Policy;
using System.Text;
using System.Threading.Tasks;
using System.Web;
using BoSSS.Solution.AdvancedSolvers;

namespace BoSSS.Application.IncompressibleNSE {

    public static class DNS_Hagen_Poiseulle {


        /// <summary>
        /// laminar solution for a Hagen Poiseulle flow (aka. flow in a circular pipe)
        /// </summary>
        public static IncompressibleControl HagenPoiseulle(string _DbPath = null, int degree = 3, int noOfCellsR = 32, int noOfCellsXi = 32, int bdfOrder = 1, double rMin = 0.3) {

            IncompressibleControl Ctrl = new IncompressibleControl();
            #region db
            //Ctrl.DbPath = @"P:\BoSSSpostprocessing\Akbari"; // _DbPath;
            // Ctrl.DbPath = @"\\dc3\userspace\akbari\cluster\Helical_DNS";
            //Ctrl.DbPath = null;
            Ctrl.DbPath = _DbPath;

            const double MaxAmp = 5;

            if (rMin != 0) {
                for (int i = 0; i < 9; i++)
                    Console.WriteLine($"Remember: r min = {rMin} !!!!!!");
            }
            Ctrl.maxAmpli = MaxAmp;

            bool transient = true;
            Ctrl.savetodb = Ctrl.DbPath != null;
            Ctrl.ProjectName = "NStransient";
            Ctrl.SessionName = "degree= " + degree + " " + "noOfCellsR= " + noOfCellsR + " " + "noOfCellsXi= " + noOfCellsXi;
            Ctrl.rMin = rMin;
            Ctrl.rMax = 1;
            #endregion

            // Solver Options
            if (transient) {
                Ctrl.NoOfTimesteps = 1000;
                Ctrl.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                Ctrl.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
                Ctrl.dtFixed = 0.01;
            } else {
                Ctrl.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                Ctrl.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            }

            Ctrl.Resolution_R = noOfCellsR;
            Ctrl.Resolution_Xi = noOfCellsXi;

            Ctrl.GridFunc = delegate {

                double[] xnodes = GenericBlas.Linspace(rMin, Ctrl.rMax, noOfCellsR + 1);
                double[] ynodes = GenericBlas.Linspace(0, 2 * Math.PI, noOfCellsXi + 1);

                GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear,
                    periodicX: false,
                    periodicY: true);

                grd.EdgeTagNames.Add(1, "Dirichlet_rmax");
                grd.EdgeTagNames.Add(2, "Dirichlet_rmin");
                grd.EdgeTagNames.Add(3, "Stuff");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double r, xi;
                    r = X[0];
                    xi = X[1];
                    if (Math.Abs(r - Ctrl.rMax) < 1E-8)
                        return 1;
                    else if (Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8)
                        return 2;
                    else
                        return 3;
                });


                grd.DefineEdgeTags(delegate (double[] _X) {
                    double r = _X[0];
                    double xi = _X[1];

                    if (Math.Abs(r - Ctrl.rMax) < 1E-8)
                        // right
                        return "Dirichlet_outer_wall";
                    else if (Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8)
                        // right
                        return "Dirichlet_inner_wall";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };
            if (bdfOrder == 3) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            } else if (bdfOrder == 1) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            } else {
                throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
            }
            // DG degree
            // =========
            Ctrl.dg_degree = degree;
             Ctrl.SetDGdegree(degree);
            // Initial Values
            // ==============
            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;
            Ctrl.HagenPoisseulle = true;

            var random = 0.1 * MaxAmp / 4;
            // Initial Values
            // ==============
            string InitialValue_ur_p =
            "static class MyInitialValue_ur_p {" // class must be static
                                                 // Warning: static constants are allowed,
                                                 // but any changes outside of the current text box in BoSSSpad
                                                 // will not be recorded for the code that is passed to the solver.
                                                 // A method, which should be used for an initial value,
                                                 // must be static!
            + " public static double GenerateRandomValue(double[] X, double t) {"
            + "    var random = new Random();"
            + "    double randomValue = random.NextDouble() * " + random + " - " + (random / 2) + ";"
            + "   return randomValue;"
            + " }"
            + "}";
            string InitialValue_uxi =
             "static class MyInitialValue_uxi {" // class must be static
                                                 // Warning: static constants are allowed,
                                                 // but any changes outside of the current text box in BoSSSpad
                                                 // will not be recorded for the code that is passed to the solver.
                                                 // A method, which should be used for an initial value,
                                                 // must be static!
             + " public static double GenerateRandomValue(double[] X, double t) {"
             + "    var random = new Random();"
             + "    double randomValue_and_lami = - " + MaxAmp + " * (X[0] / (Math.Sqrt(" + a * a + " * X[0] * X[0] + " + b * b + "))) * (" + a * a + " * (" + Ctrl.rMax * Ctrl.rMax + " - X[0] * X[0])) / (4 * " + nu + ") + random.NextDouble() * " + random + " - " + (random / 2) + ";"
             + "   return randomValue_and_lami;"
             + " }"
             + "}";


            string InitialValue_ueta =
             "static class MyInitialValue_ueta {" // class must be static
                                                  // Warning: static constants are allowed,
                                                  // but any changes outside of the current text box in BoSSSpad
                                                  // will not be recorded for the code that is passed to the solver.
                                                  // A method, which should be used for an initial value,
                                                  // must be static!
             + " public static double GenerateRandomValue(double[] X, double t) {"
             + "    var random = new Random();"
             + "    double randomValue_and_lami = " + MaxAmp + " * (X[0] / (Math.Sqrt(" + a * a + " * X[0] * X[0] +" + b * b + "))) * (" + a * b + " * (" + Ctrl.rMax * Ctrl.rMax + " - X[0] * X[0])) / (X[0] *4 * " + nu + ") + random.NextDouble() * " + random + " - " + (random / 2) + ";"
             + "   return randomValue_and_lami;"
             + " }"
             + "}";


            var Velocity_R_and_Pressure_0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ur_p.GenerateRandomValue", true, InitialValue_ur_p);
            var Velocity_XI_0 = new BoSSS.Solution.Control.Formula("MyInitialValue_uxi.GenerateRandomValue", true, InitialValue_uxi);
            var Velocity_ETA_0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ueta.GenerateRandomValue", true, InitialValue_ueta);


            Ctrl.AddInitialValue("Pressure", Velocity_R_and_Pressure_0);
            Ctrl.AddInitialValue("Velocity_R", Velocity_R_and_Pressure_0);
            Ctrl.AddInitialValue("Velocity_ETA", Velocity_ETA_0);
            Ctrl.AddInitialValue("Velocity_XI", Velocity_XI_0);



            Ctrl.AddInitialValue("Pressure", new Formula($"(X) =>0 "));
            Ctrl.AddInitialValue("Velocity_R", new Formula($"(X) => 0"));
            Ctrl.AddInitialValue("Velocity_ETA", new Formula($"(X) => {MaxAmp} *(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} ))) * ({a * b} * ({Ctrl.rMax * Ctrl.rMax} - X[0]*X[0]) )/(X[0]*4* {nu})"));
            Ctrl.AddInitialValue("Velocity_XI", new Formula($"(X) => -{MaxAmp} *(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} ))) * ({a * a} * ({Ctrl.rMax * Ctrl.rMax} - X[0]*X[0]) )/(4* {nu})"));
            // Boundary Conditions
            // ==============
            Ctrl.AddBoundaryValue("Dirichlet_outer_wall", "Velocity_R", new Formula("(X,t) =>  0", true));
            Ctrl.AddBoundaryValue("Dirichlet_outer_wall", "Velocity_ETA", new Formula("(X,t) =>0", true));
            Ctrl.AddBoundaryValue("Dirichlet_outer_wall", "Velocity_XI", new Formula("(X,t) => 0", true));
          
            Ctrl.AddBoundaryValue("Dirichlet_inner_wall", "Velocity_R", new Formula("(X,t) =>  0", true));
            Ctrl.AddBoundaryValue("Dirichlet_inner_wall", "Velocity_ETA", new Formula($"(X,t) => {MaxAmp} * (X[0] / (Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))) * ({a * b} * ({Ctrl.rMax} * {Ctrl.rMax} - X[0] * X[0])) / (X[0] * 4 * {nu})", true));
            Ctrl.AddBoundaryValue("Dirichlet_inner_wall", "Velocity_XI", new Formula($"(X,t) =>- {MaxAmp} * (X[0] / (Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))) * ({a * a} * ({Ctrl.rMax} * {Ctrl.rMax} - X[0] * X[0])) / (4 * {nu})", true));

            //Ctrl.AddBoundaryValue("Dirichlet_outer_wall", "Pressure", new Formula("(X,t) =>0", true));

            if (rMin < 10e-6) {
                Globals.activeMult = Globals.Multiplier.Bsq;
            } else {
                Globals.activeMult = Globals.Multiplier.one;
            }
            return Ctrl;
        }

    }
}