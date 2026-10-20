/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Globalization;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.XdgTimestepping;
using ilPSP.Utils;
using NUnit.Framework;

namespace BoSSS.Application.XNSFE_Solver.Tests {

    /// <summary>
    /// Consistency of linearization and evaluation, M(u0)*u0 + b(u0) == f(u0), for XNSFE
    /// on a curved (<see cref="CellType.Cube_27"/>), rotationally periodic cylinder wedge.
    /// </summary>
    /// <remarks>
    /// The check is performed at a non-trivial initial state, without any Newton iteration,
    /// see <see cref="OperatorSelfConsistencyCheck"/>.
    /// </remarks>
    [TestFixture]
    public static class CurvedRotationalPeriodicRegressionTest {

        /// <summary>
        /// Coarse Taylor-Couette setup on a curved, rotationally periodic cylinder wedge: single phase, rotating inner wall, convection and heat equation.
        /// </summary>
        public static XNSFE_Control BuildControl(int nSCells = 3, int nZCells = 1, int dgDegree = 2) {
            var C = new XNSFE_Control();
            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "CurvedRotationalPeriodicRegressionTest";

            // CylinderGrid requires at least 3 cells in a periodic direction
            double rMin = 75, rMax = 76;
            double sMin = 0, sMax = 0.0085;
            double zMin = 0, zMax = 2;
            int nRnodes = 2, nSnodes = nSCells + 1, nZnodes = nZCells + 1;

            C.GridFunc = delegate () {
                var rNodes = GenericBlas.Linspace(rMin, rMax, nRnodes);
                var sNodes = GenericBlas.Linspace(sMin, sMax, nSnodes);
                var zNodes = GenericBlas.Linspace(zMin, zMax, nZnodes);
                GridCommons grid = Grid3D.CylinderGrid(rNodes, sNodes, zNodes,
                    CellType.Cube_27, PeriodicZ: false, PeriodicS: true);

                grid.DefineEdgeTags(delegate (double[] X) {
                    double r = Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
                    // velocity and temperature boundary condition types
                    if (Math.Abs(r - rMin) <= 1.0e-8) return "wall_inner_ConstantTemperature";
                    if (Math.Abs(r - rMax) <= 1.0e-8) return "wall_outer_ConstantTemperature";
                    if (Math.Abs(X[2] - zMin) <= 1.0e-8 || Math.Abs(X[2] - zMax) <= 1.0e-8) return "wall_z_ConstantTemperature";
                    throw new ArgumentOutOfRangeException($"unexpected boundary point X=({X[0]:G17},{X[1]:G17},{X[2]:G17})");
                });
                return grid;
            };

            // rotating inner wall, outer wall at rest
            string Phi_ = "Math.Atan2(X[1], X[0])";
            double wallSpeed = 100.0;
            foreach (string species in new[] { "A", "B" }) {
                C.AddBoundaryValue("wall_inner_ConstantTemperature", $"VelocityX#{species}", $"X => {wallSpeed.ToString(CultureInfo.InvariantCulture)}*-Math.Sin({Phi_})", false);
                C.AddBoundaryValue("wall_inner_ConstantTemperature", $"VelocityY#{species}", $"X => {wallSpeed.ToString(CultureInfo.InvariantCulture)}*Math.Cos({Phi_})", false);
                C.AddBoundaryValue("wall_inner_ConstantTemperature", $"VelocityZ#{species}", "X => 0.0", false);
                C.AddBoundaryValue("wall_outer_ConstantTemperature", $"VelocityX#{species}", "X => 0.0", false);
                C.AddBoundaryValue("wall_outer_ConstantTemperature", $"VelocityY#{species}", "X => 0.0", false);
                C.AddBoundaryValue("wall_outer_ConstantTemperature", $"VelocityZ#{species}", "X => 0.0", false);
                C.AddBoundaryValue("wall_z_ConstantTemperature", $"VelocityX#{species}", "X => 0.0", false);
                C.AddBoundaryValue("wall_z_ConstantTemperature", $"VelocityY#{species}", "X => 0.0", false);
                C.AddBoundaryValue("wall_z_ConstantTemperature", $"VelocityZ#{species}", "X => 0.0", false);

                C.AddBoundaryValue("wall_inner_ConstantTemperature", $"Temperature#{species}", "X => 1.0", false);
                C.AddBoundaryValue("wall_outer_ConstantTemperature", $"Temperature#{species}", "X => 0.0", false);
                C.AddBoundaryValue("wall_z_ConstantTemperature", $"Temperature#{species}", "X => 0.0", false);
            }

            // single phase, no cut cells
            C.AddInitialValue("Phi", 1.0);

            // non-trivial initial velocity (linear profile across the gap), at which the check is performed
            foreach (string species in new[] { "A", "B" }) {
                C.InitialValues_Evaluators.Add($"VelocityX#{species}", X => {
                    double r0 = Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
                    double phi0 = Math.Atan2(X[1], X[0]);
                    double frac = (rMax - r0) / (rMax - rMin);
                    return wallSpeed * frac * -Math.Sin(phi0);
                });
                C.InitialValues_Evaluators.Add($"VelocityY#{species}", X => {
                    double r0 = Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
                    double phi0 = Math.Atan2(X[1], X[0]);
                    double frac = (rMax - r0) / (rMax - rMin);
                    return wallSpeed * frac * Math.Cos(phi0);
                });
            }
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.PhysicalParameters.rho_A = 1.0; C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_A = 1.0; C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;
            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            C.ThermalParameters.rho_A = 1.0; C.ThermalParameters.rho_B = 1.0;
            C.ThermalParameters.k_A = 1.0; C.ThermalParameters.k_B = 1.0;
            C.ThermalParameters.c_A = 1.0; C.ThermalParameters.c_B = 1.0;
            C.ThermalParameters.IncludeConvection = true;
            C.ThermalParameters.hVap = 0.0;
            C.ThermalParameters.T_sat = 0.0;

            C.SetDGdegree(dgDegree);

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            // no Newton iterations: only the initial residual evaluation is needed
            C.NonLinearSolver.MinSolverIterations = 0;
            C.NonLinearSolver.MaxSolverIterations = 0;
            C.NonLinearSolver.ConvergenceCriterion = 1e-30;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            return C;
        }

        /// <summary>
        /// Runs the XNSFE solver for <see cref="BuildControl"/> and returns the values recorded by <see cref="OperatorSelfConsistencyCheck"/>.
        /// </summary>
        public static (double LastAbsoluteError, double LastComparisonDenominator, double LastRelativeError, double MaxRelativeError) RunAndGetSelfConsistency(int nSCells = 3, int nZCells = 1, int dgDegree = 2, bool doMpiInit = true) {
            OperatorSelfConsistencyCheck.Reset();
            bool oldForce = OperatorSelfConsistencyCheck.ForceCheckEveryCall;
            OperatorSelfConsistencyCheck.ForceCheckEveryCall = true;

            var C = BuildControl(nSCells, nZCells, dgDegree);

            bool mustFinalizeMPI = doMpiInit && BoSSS.Solution.Application.InitMPI(new string[0]);
            try {
                using (var solver = new XNSFE()) {
                    solver.Init(C);
                    try {
                        solver.RunSolverMode();
                    } catch (ArithmeticException) {
                        // expected, since no Newton iterations are performed
                    }
                }
            } finally {
                OperatorSelfConsistencyCheck.ForceCheckEveryCall = oldForce;
                if (mustFinalizeMPI)
                    BoSSS.Solution.Application.FinalizeMPI();
            }

            if (!OperatorSelfConsistencyCheck.LastAbsoluteError.HasValue || !OperatorSelfConsistencyCheck.LastComparisonDenominator.HasValue || !OperatorSelfConsistencyCheck.MaxRelativeError.HasValue)
                throw new InvalidOperationException("Operator self-consistency check did not run -- test setup is broken (expected it to fire on the initial residual evaluation).");

            return (OperatorSelfConsistencyCheck.LastAbsoluteError.Value, OperatorSelfConsistencyCheck.LastComparisonDenominator.Value, OperatorSelfConsistencyCheck.LastRelativeError.Value, OperatorSelfConsistencyCheck.MaxRelativeError.Value);
        }

        /// <summary>
        /// Checks that the linearization M(u0)*u0 + b(u0) is consistent with the evaluation f(u0).
        /// </summary>
        [Test]
        public static void CurvedElementRotationalPeriodicSelfConsistency() {
            var (lastAbsErr, lastComp, lastRelErr, maxRelErr) = RunAndGetSelfConsistency();

            Console.WriteLine($"CurvedElementRotationalPeriodicSelfConsistency: lastAbsoluteError={lastAbsErr:G17}, lastComparisonDenominator={lastComp:G17}, lastRelativeError={lastRelErr:G17}, maxRelativeError={maxRelErr:G17}");

            const double RelativeErrorThreshold = 1.0e4;
            Assert.Less(maxRelErr, RelativeErrorThreshold,
                $"M(u0)*u0+b(u0) vs. f(u0) self-consistency mismatch on the curved rotationally-periodic regression configuration is too large (worst relative error {maxRelErr:G6}, threshold {RelativeErrorThreshold:G6}) -- one of the curved-element/rotational-periodic-boundary fixes (LECVolumeQuadrature2.cs, LECEdgeQuadrature2.cs, DifferentialOperator.cs, NECQuadratureEdge.cs, JacobianParamUpdate.cs, PeriodicBoundaryUtils.cs) may have regressed.");
        }
    }
}
