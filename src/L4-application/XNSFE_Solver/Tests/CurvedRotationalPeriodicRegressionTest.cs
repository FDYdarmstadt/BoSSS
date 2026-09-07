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
    /// Permanent regression test for curved-element + rotational-periodic-boundary support in
    /// the XDG/XNSFE discretization pipeline: a minimal cylindrical-wedge Taylor-Couette-style
    /// configuration (curved <see cref="CellType.Cube_27"/> elements, rotationally-periodic in
    /// the azimuthal direction) that is known to have been broken -- by a combination of
    /// six distinct bugs, in <c>LECVolumeQuadrature2.cs</c>, <c>LECEdgeQuadrature2.cs</c>,
    /// <c>DifferentialOperator.cs</c>, <c>NECQuadratureEdge.cs</c>, <c>JacobianParamUpdate.cs</c>
    /// and <c>PeriodicBoundaryUtils.cs</c> -- prior to a dedicated investigation that fixed all
    /// six. This does not re-run the full historical investigation; it distills it into a single,
    /// fast, deterministic assertion.
    /// </summary>
    /// <remarks>
    /// Rather than checking for exact convergence or an exactly-zero residual, this test checks
    /// the operator's own <c>M(u0)*u0+b(u0) == f(u0)</c> self-consistency identity (linearization
    /// vs. direct nonlinear evaluation of the same operator, at the same state <c>u0</c>) via
    /// <see cref="OperatorSelfConsistencyCheck"/> -- the same identity, and the same numbers,
    /// that <see cref="NonlinearSolver.EvaluateOperator"/> already computes internally as part of
    /// normal production solves (there, only reported via a <see cref="Console.Error"/> warning).
    /// This identity must hold to within floating-point/quadrature noise for ANY correctly
    /// assembled Jacobian, independent of whether the underlying nonlinear solve ever converges,
    /// or even runs at all -- <see cref="BuildControl"/> seeds a non-trivial (already-rotating)
    /// initial velocity field for exactly this reason: the check is performed at that initial
    /// state, with the nonlinear solver bounded to zero iterations, so this test needs only ONE
    /// assembled residual/Jacobian pair, never an actual (expensive) linear solve, keeping it fast.
    ///
    /// A genuinely exact ("must be zero") threshold is not appropriate here: a separate,
    /// structurally distinct, still-open defect in momentum's periodic-edge convection Jacobian
    /// (a `JacobianParamUpdate` linearization-snapshot projection error, unrelated to the six
    /// fixes this test guards) is known to leave a small residual on this exact configuration even
    /// with all six fixes applied. The threshold below was determined empirically -- measured
    /// relative errors (<see cref="OperatorSelfConsistencyCheck.MaxRelativeError"/>) at the exact
    /// configuration <see cref="CurvedElementRotationalPeriodicSelfConsistency"/> uses:
    /// current (all six fixes applied): ~1.8; any ONE of the six fixes reverted to its
    /// pre-investigation state: from ~3.9e6 (<c>PeriodicBoundaryUtils.cs</c>) up to ~1.5e9
    /// (<c>NECQuadratureEdge.cs</c>) -- reverting <c>JacobianParamUpdate.cs</c> alone does not even
    /// compile, since <c>DifferentialOperator.cs</c>'s fix depends on a property it adds. The
    /// chosen threshold sits about 3 orders of magnitude above the current value and at least 2
    /// orders of magnitude below the smallest reverted-fix value -- ample margin on both sides.
    /// </remarks>
    [TestFixture]
    public static class CurvedRotationalPeriodicRegressionTest {

        /// <summary>
        /// Builds a minimal, deliberately coarse cylindrical-wedge Taylor-Couette-style control
        /// object: curved (<see cref="CellType.Cube_27"/>) elements, rotationally-periodic in the
        /// azimuthal ("S") direction, single phase via a spatially-constant level-set (so the XDG
        /// machinery is genuinely engaged but no cell is actually cut), rotating no-slip wall at
        /// the inner cylinder, convection and heat conduction both active (both are needed to
        /// engage the fixed code paths -- see the investigation notes).
        /// </summary>
        public static XNSFE_Control BuildControl(int nSCells = 3, int nZCells = 1, int dgDegree = 2) {
            var C = new XNSFE_Control();
            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "CurvedRotationalPeriodicRegressionTest";

            // Minimal legal CylinderGrid resolution: at least 1 cell in R (never periodic), and
            // at least 3 cells in S (azimuthal) since Grid3D.CylinderGrid requires >= 3 cells in
            // any direction that is periodic.
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
                    // the tag must contain both the velocity ("wall") and thermal
                    // ("ConstantTemperature") BC-type substrings BoundaryCondMap matches on.
                    if (Math.Abs(r - rMin) <= 1.0e-8) return "wall_inner_ConstantTemperature";
                    if (Math.Abs(r - rMax) <= 1.0e-8) return "wall_outer_ConstantTemperature";
                    // only reachable since PeriodicZ is off -- the z=zMin/zMax faces need a tag.
                    if (Math.Abs(X[2] - zMin) <= 1.0e-8 || Math.Abs(X[2] - zMax) <= 1.0e-8) return "wall_z_ConstantTemperature";
                    throw new ArgumentOutOfRangeException($"unexpected boundary point X=({X[0]:G17},{X[1]:G17},{X[2]:G17})");
                });
                return grid;
            };

            // Rotating no-slip wall at the inner cylinder (r=rMin), stationary at the outer
            // (r=rMax) -- tangential velocity v = speed*(-sin(phi),cos(phi),0).
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

            // Single phase via a trivial (never-cutting) level-set: the XDG machinery is
            // genuinely engaged, but zero cells are actually cut -- keeps the test focused on the
            // curved-element/periodic-rotation machinery, not on cut-cell/agglomeration behavior.
            C.AddInitialValue("Phi", 1.0);

            // A non-trivial (already-rotating) initial velocity, matching the wall speed profile
            // linearly across the gap: this is deliberately NOT the true steady solution, just
            // something with genuine rotation/gradients from the start, so the self-consistency
            // check (evaluated at THIS initial state) has something non-zero to be inconsistent
            // about -- avoiding the need for an actual (expensive) Newton solve just to move the
            // state away from the trivial all-zero default.
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
            // We only need ONE assembled residual/Jacobian pair, at the non-trivial initial state
            // seeded above -- see remarks on the class -- not an actual solve, so the solver is
            // bounded to zero iterations. The very first residual evaluation (needed regardless,
            // to check the initial convergence criterion) still happens, which is all this test
            // needs.
            C.NonLinearSolver.MinSolverIterations = 0;
            C.NonLinearSolver.MaxSolverIterations = 0;
            C.NonLinearSolver.ConvergenceCriterion = 1e-30;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            return C;
        }

        /// <summary>
        /// Runs <see cref="BuildControl"/> through the real XNSFE solver (zero Newton iterations
        /// -- see remarks on the class) and returns the resulting M(u0)*u0+b(u0) vs. f(u0)
        /// self-consistency numbers recorded by <see cref="OperatorSelfConsistencyCheck"/>. Forces
        /// the check to run on this single residual evaluation (it would otherwise only run every
        /// 10th/1000th call in DEBUG/RELEASE builds, respectively).
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
                        // Expected: with MaxSolverIterations=0 and a deliberately tight
                        // ConvergenceCriterion, the solver reports non-convergence -- irrelevant
                        // here, since we only need the internal self-consistency check (forced to
                        // run, see above) to have fired on the initial residual evaluation.
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
        /// Regression test: on a curved, rotationally-periodic wedge, the linearization
        /// M(u0)*u0+b(u0) must stay close to a direct nonlinear evaluation f(u0) of the same
        /// operator, at the same state -- to within a bound well above the small, still-open,
        /// separately-tracked residual left by this configuration, but far below what any of the
        /// six historical bugs (in curved-cell orthonormalization, periodic-rotation parameter
        /// handling, quadrature-order selection, or vector-tuple grouping) would produce if
        /// reintroduced.
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
