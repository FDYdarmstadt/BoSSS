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
    internal class Centrifuge {

        /// <summary>
        /// Centrifuge flow (aka. flow in a rotating pipe)
        /// </summary>
        public static HelicalControl Centrifuge_Flow(string _DbPath = null, int degree = 3, int noOfCellsR = 64, int noOfCellsXi = 64, int bdfOrder = 3, int numOfTimesteps = 1, double deltaT = 1, double rMin = 0, double MaxAmp = 5) {

            HelicalControl Ctrl = new HelicalControl();

            // Data Base
            // ==============
            #region db
            Ctrl.DbPath = _DbPath;
            #endregion

            // Settings
            // ==============
            #region Settings
            Ctrl.maxAmpli = MaxAmp;
            Ctrl.rMin = rMin;
            Ctrl.rMax = 1;
            Ctrl.dg_degree = degree;
            Ctrl.SetDGdegree(degree);
            #endregion

            // Grid
            // ==============
            #region Grid            
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
                    if(Math.Abs(r - Ctrl.rMax) < 1E-8)
                        return 1;
                    else if(Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8)
                        return 2;
                    else
                        return 3;
                });
                return grd;
            };
            #endregion

            // Time Stepping
            // ==============
            #region Timestepping
            if (deltaT >= 10E10) {
                Ctrl.dtFixed = deltaT;
                Ctrl.NoOfTimesteps = 1;
                Ctrl.steady = true;
            } else {
                Ctrl.dtFixed = deltaT;
                Ctrl.NoOfTimesteps = numOfTimesteps;
                Ctrl.steady = false;
            }

            if(bdfOrder == 3) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            } else if(bdfOrder == 1) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            } else {
                throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
            }


            #endregion

            // Initial Values
            // ==============
            #region InitialValue
            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;

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

            Ctrl.AddInitialValue("Pressure", new Formula("(X) =>0"));
            Ctrl.AddInitialValue("ur", new Formula("(X) => 0"));
            Ctrl.AddInitialValue("ueta", new Formula($"(X) => 0"));
            Ctrl.AddInitialValue("uxi", new Formula($"(X) => 0"));
            #endregion

            // Boundary Conditions
            // ==============
            #region BoundaryCOnditions
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) =>  0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula($"(X,t) =>(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} )))*{a}*{MaxAmp}* X[0]", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula($"(X,t) => (X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} )))*{b}*{MaxAmp}", true));
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) =>0", true));
            #endregion

            //Ctrl.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { ConvergenceCriterion = 1e-13 , TargetBlockSize =1000000};
            return Ctrl;
        }

    }
}
