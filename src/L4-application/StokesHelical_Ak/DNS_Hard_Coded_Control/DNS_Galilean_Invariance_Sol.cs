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
    internal class DNS_Galilean_Invariance_Sol {
        /// <summary>
        /// ATTENTION!!!!!!! Singularity at r =0;
        /// Exact Solution for helically symetric Navier Stokes Equation. (Dierkes et. al. 2020, Equation )
        /// </summary>
        public static HelicalControl Gali_Invar_Sol(string _DbPath = null, int degree = 2, int noOfCellsR = 4, int noOfCellsXi = 4, int dtRefining = 128, int bdfOrder = 1, double Tend = 2 * Math.PI, double rMin = 0) {

            HelicalControl Ctrl = new HelicalControl();
            #region db
            //Ctrl.DbPath = @"P:\BoSSSpostprocessing\Akbari"; // _DbPath;
            // Ctrl.DbPath = @"\\dc3\userspace\akbari\cluster\Helical_DNS";
            //Ctrl.DbPath = null;
            Ctrl.DbPath = _DbPath;

            const double MaxAmp = 0.5;
            const double t0 = 2;

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
            double dt = Tend / (dtRefining * 1000);

            Ctrl.dtFixed = dt;
            if(bdfOrder == 3) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            } else if(bdfOrder == 1) {
                Ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            } else {
                throw new ArgumentException("Unsupported BDF scheme: " + bdfOrder);
            }
            //Ctrl.NoOfTimesteps = dtRefining * 200*4;
            Ctrl.NoOfTimesteps = dtRefining * 100;
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
            Ctrl.HagenPoisseulle = false;
            Ctrl.AddInitialValue("Pressure", new Formula($"(X) =>(({MaxAmp}*{MaxAmp})/(2*X[0]*X[0]))*Math.Exp(-(X[0] * X[0])/(2*{nu}*(0+{t0})))"));
            Ctrl.AddInitialValue("ur"      , new Formula($"(X) => ({MaxAmp}/X[0])*Math.Exp(-(X[0] * X[0])/(4*{nu}*(0+{t0})))"));
            Ctrl.AddInitialValue("ueta"    , new Formula($"(X) =>(({MaxAmp}*{b}*X[1]*(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))))/(2*{nu}*{a}*X[0]*(0+{t0})))*Math.Exp(-(X[0] * X[0])/(4*{nu}*(0+{t0})))"));
            Ctrl.AddInitialValue("uxi"     , new Formula($"(X) =>(({MaxAmp}    *X[1]*(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))))/(2*{nu}         *(0+{t0})))*Math.Exp(-(X[0] * X[0])/(4*{nu}*(0+{t0})))"));
            // Boundary Conditions
            // ==============
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula($"(X,t)       => ({MaxAmp}/X[0])*Math.Exp(-(X[0] * X[0])/(4*{nu}*(t+{t0})))", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula($"(X,t)     =>(({MaxAmp}*{b}*X[1]*(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))))/(2*{nu}*{a}*X[0]*(t+{t0})))*Math.Exp(-(X[0] * X[0])/(4*{nu}*(t+{t0})))", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula($"(X,t)      =>(({MaxAmp}    *X[1]*(X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))))/(2*{nu}*         (t+{t0})))*Math.Exp(-(X[0] * X[0])/(4*{nu}*(t+{t0})))", true));
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula($"(X,t) =>(({MaxAmp}*{MaxAmp})/(2*X[0]*X[0]))*Math.Exp(-(X[0] * X[0])/(2*{nu}*(t+{t0})))", true));
            Ctrl.ImmediatePlotPeriod = 1;
            //Ctrl.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { ConvergenceCriterion = 1e-13 , TargetBlockSize =1000000};
            return Ctrl;
        }

    }
}
