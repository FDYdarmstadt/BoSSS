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

namespace StokesHelical_Ak {

    public static class DNS_Hagen_Poiseulle {


        /// <summary>
        /// laminar solution for a Hagen Poiseulle flow (aka. flow in a circular pipe)
        /// </summary>
        public static HelicalControl HagenPoiseulle(string _DbPath = null, int degree = 5, int noOfCellsR = 512, int noOfCellsXi = 512, int dtRefining = 4, int bdfOrder = 3, double Tend = 2 * Math.PI, double rMin = 0) {

            HelicalControl Ctrl = new HelicalControl();
            #region db
            //Ctrl.DbPath = @"P:\BoSSSpostprocessing\Akbari"; // _DbPath;
            // Ctrl.DbPath = @"\\dc3\userspace\akbari\cluster\Helical_DNS";
            //Ctrl.DbPath = null;
            Ctrl.DbPath = _DbPath;

            const double MaxAmp = 5;

            if(rMin != 0) {
                for(int i = 0; i < 9; i++)
                    Console.WriteLine($"Remember: r min = {rMin} !!!!!!");
            }
            Ctrl.maxAmpli = MaxAmp;


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

                GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear, 
                    periodicX: false,
                    periodicY: true);
                //for(int i = 0; i < 9; i++)
                //    Console.WriteLine("Remember: turn periodic on again !!!!!! All Dirichlet!!");

                //GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear,
                //periodicY: true);

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

                //grd.DefineEdgeTags(delegate (double[] _X) {
                //    var X = _X;
                //    double r, xi;
                //    r = X[0];
                //    xi = X[1];
                //    if(Math.Abs(r - Ctrl.rMax) < 1E-8) {
                //        return "Dirichlet_rmax";
                //    }
                //if(Math.Abs(r - Ctrl.rMin) < 1E-8) {
                //    return "Dirichlet_rmin";
                //}

                //Math.Abs(xi - 2 * Math.PI) < 1E-8 || Math.Abs(xi - 0) < 1E-8)
                //        return 1;
                //    else
                //        return 2;
                //    throw new ArgumentException("unknown bndy coordinate: " + new Vector(r, xi));
                //});


                return grd;
            };
            double dt = Tend / (dtRefining * 10);

            Ctrl.dtFixed = dt;
            if(bdfOrder == 3) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            } else if(bdfOrder == 1) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            } else {
                throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
            }
            //Ctrl.NoOfTimesteps = dtRefining * 200*4;
            Ctrl.NoOfTimesteps = dtRefining;
            Ctrl.steady = false;
            Ctrl.ExactResidual = false;

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

            string InitialValue =
            "static class MyInitialValue {" // class must be static
            // Warning: static constants are allowed,
            // but any changes outside of the current text box in BoSSSpad
            // will not be recorded for the code that is passed to the solver.
            // A method, which should be used for an initial value,
            // must be static!
            + " public static double GenerateRandomValue(double[] X, double t) {"
            + "    var random = new Random();"
            + "    double randomValue = random.NextDouble() * 200 - 100;"
            + "   return randomValue;"
            + " }"
            + "}";
            var fo = new BoSSS.Solution.Control.Formula("MyInitialValue.GenerateRandomValue", true, InitialValue);
            //Ctrl.AddInitialValue("Pressure", fo);
            //Ctrl.AddInitialValue("ur", fo);
            //Ctrl.AddInitialValue("ueta", fo);
            //Ctrl.AddInitialValue("uxi", fo);

            Ctrl.AddInitialValue("Pressure", new Formula($"(X) =>0 "));
            Ctrl.AddInitialValue("ur", new Formula($"(X) => 0"));
            Ctrl.AddInitialValue("ueta", new Formula($"(X) => {MaxAmp} *(X[0]/(Math.Sqrt({ a * a} * X[0] * X[0] + { b * b} ))) * ({a * b} * ({Ctrl.rMax* Ctrl.rMax} - X[0]*X[0]) )/(X[0]*4* {nu})"));
            Ctrl.AddInitialValue("uxi", new Formula($"(X) => -{MaxAmp} *(X[0]/(Math.Sqrt({ a * a} * X[0] * X[0] + { b * b} ))) * ({a * a} * ({Ctrl.rMax* Ctrl.rMax} - X[0]*X[0]) )/(4* {nu})"));
            // Boundary Conditions
            // ==============
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) =>  0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) =>0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) => 0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) =>0", true));

            //Ctrl.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { ConvergenceCriterion = 1e-13 , TargetBlockSize =1000000};
            return Ctrl;

        }

    }
}
