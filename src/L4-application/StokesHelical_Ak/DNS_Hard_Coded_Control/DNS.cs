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

    public static class DNS {

        public static HelicalControl SimpleFlow(string _DbPath = @"", int degree = 5, int noOfCellsR = 512, int noOfCellsXi = 512, int dtRefining = 4, string bdfOrder = "BDF3", double rMin = 0.1) {

            HelicalControl Ctrl = new HelicalControl();
            #region db
            //Ctrl.DbPath = @"P:\BoSSSpostprocessing\Akbari"; // _DbPath;
             Ctrl.DbPath = @"\\dc3\userspace\akbari\cluster\Helical_DNS";
            //Ctrl.DbPath = null;

            double MaxAmp = 0.1;
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
            double dt = 2 * Math.PI / (dtRefining * 10);
            //double dt = Math.Pow(10,20)*2 * Math.PI / dtRefining;

            Ctrl.dtMax = dt;
            Ctrl.dtMin = dt;
            if(bdfOrder == "BDF3") {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            } else if(bdfOrder == "BDF1") {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            } else {
                throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
            }
            Ctrl.NoOfTimesteps = dtRefining * 200*4;
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

            Ctrl.AddInitialValue("Pressure", new Formula("(X) => 0"));
            Ctrl.AddInitialValue("ur", new Formula("(X) => 0"));
            Ctrl.AddInitialValue("ueta", new Formula("(X) => " + MaxAmp + " "));
            Ctrl.AddInitialValue("uxi", new Formula("(X) => 0"));
            // Boundary Conditions
            // ==============
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) => 0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) => 0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) => " + MaxAmp + " ", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) =>" + MaxAmp + "*(1-Math.Exp(-t))", true));

            //Ctrl.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { ConvergenceCriterion = 1e-13 , TargetBlockSize =1000000};
            return Ctrl;
        }

    }
}
