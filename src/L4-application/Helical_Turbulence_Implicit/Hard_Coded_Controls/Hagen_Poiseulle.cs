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

namespace BoSSS.Application.IncompressibleNSE.Hard_Coded_Controls {
    internal class Hagen_Poiseulle {
        /// <summary>
        /// laminar solution for a Hagen Poiseulle flow (aka. flow in a circular pipe)
        /// </summary>
        public static IncompressibleControl HagenPoiseulle(string _DbPath = null, int degree = 3, int noOfCellsR = 32, int noOfCellsXi = 32, int bdfOrder = 1, double rMin = 0.4, double MaxAmp = 1, double rMax =2, double deltaT = 0.01,int numOfTimesteps = 1) {

            IncompressibleControl Ctrl = new IncompressibleControl();

            // Database
            // ==============
            #region db
            Ctrl.DbPath = _DbPath;
            #endregion
            // Settings
            // ==============
            #region Settings
            Ctrl.maxAmpli = MaxAmp;
            Ctrl.rMin = rMin;
            Ctrl.rMax = rMax;
            Ctrl.HagenPoisseulle = true;
            Ctrl.dg_degree = degree;
            Ctrl.SetDGdegree(degree);
            if(rMin != 0) {
                for(int i = 0; i < 9; i++)
                    Console.WriteLine($"Remember: r min = {rMin} !!!!!!");
            }
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


                grd.DefineEdgeTags(delegate (double[] _X) {
                    double r = _X[0];
                    double xi = _X[1];

                    if(Math.Abs(r - Ctrl.rMax) < 1E-8)
                        // right
                        return "Dirichlet_outer_wall";
                    else if(Ctrl.rMin >= 1E-6 && Math.Abs(r - Ctrl.rMin) < 1E-8)
                        // right
                        return "Dirichlet_inner_wall";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };
            #endregion
            // Time Stepping
            // ==============
            #region Timestepping
            if(deltaT >= 10E10) {
                Ctrl.dtFixed = deltaT;
                Ctrl.NoOfTimesteps = 1;
                Ctrl.steady = true;
                Ctrl.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                Ctrl.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            } else {
                Ctrl.dtFixed = deltaT;
                Ctrl.NoOfTimesteps = numOfTimesteps;
                Ctrl.steady = false;
                Ctrl.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                Ctrl.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            }
            #endregion
            #region Solveroptions
            Ctrl.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            #endregion

            // Initial Values
            // ==============
            #region InitialValues
            Ctrl.AddInitialValue("Pressure", new Formula($"(X) =>   0"));
            Ctrl.AddInitialValue("ur", new Formula($"(X) =>         0"));
            Ctrl.AddInitialValue("ueta", new Formula($"(X) =>       0"));
            Ctrl.AddInitialValue("uxi", new Formula($"(X) =>        0"));
            #endregion
            // Boundary Conditions
            // ==============
            #region BoundaryConditions
            Ctrl.AddBoundaryValue("Dirichlet", "ur", new Formula("(X,t) =>  0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "ueta", new Formula("(X,t) =>0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "uxi", new Formula("(X,t) => 0", true));
            Ctrl.AddBoundaryValue("Dirichlet", "Pressure", new Formula("(X,t) =>0", true));
            #endregion

            return Ctrl;
        }

    }
}
