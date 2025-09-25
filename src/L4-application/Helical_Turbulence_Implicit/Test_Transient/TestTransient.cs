using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.XdgTimestepping;
using ilPSP.Utils;
using MathNet.Numerics;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using MPI.Wrappers;
using NUnitLite;
using System.Collections;
using System.Diagnostics;
using System.IO;
using System.Security.Policy;
using System.Web;
using BoSSS.Solution.AdvancedSolvers;



namespace BoSSS.Application.IncompressibleNSE.Test_Transient {
    [TestFixture]
    static public class TestTransient {
        #region Hagen Poiseulle
        /// <summary>
        /// Full Navier Stokes Hagen Poiseulle flow (aka. Pipe flow),
        /// With 0,1% White Noise over laminar Solutin
        /// Reynolds 10
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void Transient_HP_Re_10_White_Noise_1_ProMil_with_R0fix([Values(3)] int pOrder = 3) {

            double maxAmpitude = 40;
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");
            var ctrl = BoSSS.Application.IncompressibleNSE.Hard_Coded_Controls.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1000, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, MaxAmp: maxAmpitude);


            ctrl.savetodb = true;
            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();
            ctrl.ImmediatePlotPeriod = -1;
            // Initial Values
            // ==============
            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;
            ctrl.HagenPoisseulle = true;

            var random = 0.001 * maxAmpitude;
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
             + "    double randomValue_and_lami = - " + maxAmpitude + " * (X[0] / (Math.Sqrt(" + a * a + " * X[0] * X[0] + " + b * b + "))) * (" + a * a + " * (" + ctrl.rMax * ctrl.rMax + " - X[0] * X[0])) / (4 * " + nu + ") + random.NextDouble() * " + random + " - " + (random / 2) + ";"
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
             + "    double randomValue_and_lami = " + maxAmpitude + " * (X[0] / (Math.Sqrt(" + a * a + " * X[0] * X[0] +" + b * b + "))) * (" + a * b + " * (" + ctrl.rMax * ctrl.rMax + " - X[0] * X[0])) / (X[0] *4 * " + nu + ") + random.NextDouble() * " + random + " - " + (random / 2) + ";"
             + "   return randomValue_and_lami;"
             + " }"
             + "}";


            var Velocity_R_and_Pressure_0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ur_p.GenerateRandomValue", true, InitialValue_ur_p);
            var Velocity_XI_0 = new BoSSS.Solution.Control.Formula("MyInitialValue_uxi.GenerateRandomValue", true, InitialValue_uxi);
            var Velocity_ETA_0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ueta.GenerateRandomValue", true, InitialValue_ueta);


            ctrl.AddInitialValue("Pressure", Velocity_R_and_Pressure_0);
            ctrl.AddInitialValue("Velocity_R", Velocity_R_and_Pressure_0);
            ctrl.AddInitialValue("Velocity_ETA", Velocity_ETA_0);
            ctrl.AddInitialValue("Velocity_XI", Velocity_XI_0);

            // Boundary Conditions
            // ==============
            ctrl.AddBoundaryValue("Dirichlet_outer_wall", "Velocity_R", new Formula("(X,t) =>  0", true));
            ctrl.AddBoundaryValue("Dirichlet_outer_wall", "Velocity_ETA", new Formula("(X,t) =>0", true));
                ctrl.AddBoundaryValue("Dirichlet_outer_wall", "Velocity_XI", new Formula("(X,t) => 0", true));

            ctrl.AddBoundaryValue("Dirichlet_inner_wall", "Velocity_R", new Formula("(X,t) =>  0", true));
            ctrl.AddBoundaryValue("Dirichlet_inner_wall", "Velocity_ETA", new Formula($"(X,t) => {maxAmpitude} * (X[0] / (Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))) * ({a * b} * ({ctrl.rMax} * {ctrl.rMax} - X[0] * X[0])) / (X[0] * 4 * {nu})", true));
            ctrl.AddBoundaryValue("Dirichlet_inner_wall", "Velocity_XI", new Formula($"(X,t) =>- {maxAmpitude} * (X[0] / (Math.Sqrt({a * a} * X[0] * X[0] + {b * b}))) * ({a * a} * ({ctrl.rMax} * {ctrl.rMax} - X[0] * X[0])) / (4 * {nu})", true));
            

            var ur0_p0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ur_p.GenerateRandomValue", true, InitialValue_ur_p);
            var uxi0 = new BoSSS.Solution.Control.Formula("MyInitialValue_uxi.GenerateRandomValue", true, InitialValue_uxi);
            var ueta0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ueta.GenerateRandomValue", true, InitialValue_ueta);


            ctrl.AddInitialValue("Pressure", ur0_p0);
            ctrl.AddInitialValue("ur", ur0_p0);
            ctrl.AddInitialValue("ueta", ueta0);
            ctrl.AddInitialValue("uxi", uxi0);

            using(var solverStat = new Helical_Turbulence_Implicit_Main()) {
                solverStat.Init(ctrl);
                solverStat.RunSolverMode();


                // CloneAs zuvor
                SinglePhaseField exSol_ur = solverStat.Velocity[0].CloneAs();
                SinglePhaseField exSol_uxi = solverStat.Velocity[1].CloneAs();
                SinglePhaseField exSol_ueta = solverStat.Velocity[2].CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();
                SinglePhaseField pressureError = solverStat.Pressure.CloneAs();

                // Exact Solution
                // Func<double[], double> uxi = (X) => -maxAmpitude * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (4 * nu);
                // Func<double[], double> ueta = (X) => maxAmpitude * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                // Func<double[], double> ur = (X) => 0;
                // Func<double[], double> pressure = (X) => 0;

                ScalarFunction ur = (MultidimensionalArray inp, MultidimensionalArray outp) => {
                    int N = inp.GetLength(0); // number of evaluation points
                    for(int i = 0; i < N; i++) {
                        outp[i] = 0.0;
                    }
                };

                ScalarFunction pressure = (MultidimensionalArray inp, MultidimensionalArray outp) => {
                    int N = inp.GetLength(0); // number of evaluation points
                    for(int i = 0; i < N; i++) {
                        outp[i] = 0.0;
                    }
                };

                ScalarFunction uxi = (MultidimensionalArray inp, MultidimensionalArray outp) => {
                    int N = inp.GetLength(0);          // Anzahl Punkte
                    double rMax = solverStat.Control.rMax;
                    for(int i = 0; i < N; i++) {
                        double x = inp[i, 0];          // 1. Koordinate
                        double denom = Math.Sqrt(a * a * x * x + b * b);
                        double val = -maxAmpitude * (x / denom)
                                     * (a * a * (rMax * rMax - x * x)) / (4.0 * nu);
                        outp[i] = val;
                    }
                };


                ScalarFunction ueta = (MultidimensionalArray inp, MultidimensionalArray outp) => {
                    int N = inp.GetLength(0);
                    double rMax = solverStat.Control.rMax;
                    for(int i = 0; i < N; i++) {
                        double x = inp[i, 0];
                        double denom = Math.Sqrt(a * a * x * x + b * b);

                        // ursprüngliche Formel hatte:  (...) / (x * 4*nu)
                        // → numerisch instabil bei x≈0; hier stabilisiert:
                        double val = maxAmpitude * (a * b * (rMax * rMax - x * x))
                                     / (4.0 * nu * denom);

                        outp[i] = val;
                    }
                };

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                pressureError.ProjectField(pressure);
                pressureError.Acc(-1, solverStat.Pressure);

                double ur_L2 = solverStat.Velocity[0].L2Error(ur);
                double uxi_L2 = solverStat.Velocity[1].L2Error(exSol_uxi);
                double ueta_L2 = solverStat.Velocity[2].L2Error(exSol_ueta);
                double psiErrorMean = pressureError.GetMeanValueTotal(null);
                pressureError.AccConstant(-psiErrorMean);

                double pressure_L2 = pressureError.L2Norm();
                Console.WriteLine($"ur       L2 Error/maxAmpitude: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error/maxAmpitude: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error/maxAmpitude: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error/maxAmpitude: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2 / maxAmpitude, 1.0e-5, $"ur L2 Error/maxAmpitude out of range: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(uxi_L2 / maxAmpitude, 1.0e-5, $"uxi L2 Error/maxAmpitude out of range: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(ueta_L2 / maxAmpitude, 1.0e-4, $"ueta L2 Error/maxAmpitude out of range: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2 / maxAmpitude, 1.0e-4, $"pressure L2 Error/maxAmpitude out of range: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                // Exact solutions als ScalarFunction (input: [Npts, D], output: [Npts])

            }

        }

        #endregion
    }
}
