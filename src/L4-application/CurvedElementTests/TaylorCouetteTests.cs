using System;
using System.Linq;
using System.Reflection;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.CurvedElementTests {

    /// <summary>
    /// Steady laminar Taylor-Couette flow (Stokes and Navier-Stokes) on the curved grid of <see cref="CurvedWedgeHelpers.CreateGrid"/> (rotationally periodic wedge or full ring, see <see cref="FullRing"/>),
    /// solved with the (single-phase) XNSE solver.
    /// </summary>
    [TestFixture]
    public static class TaylorCouetteTests {

        /// <summary>
        /// MPI init.
        /// </summary>
        [OneTimeSetUp]
        public static void SetUp() {
            BoSSS.Solution.Application.InitMPI();
        }

        /// <summary>
        /// false: rotationally periodic wedge; true: full ring, without periodic boundaries,
        /// e.g. to check whether an effect is caused by the periodic boundary
        /// </summary>
        internal const bool FullRing = false;

        // must match the radial extent of the grid, see CurvedWedgeHelpers.CreateGrid()
        const double r_i = 1.0, r_o = 2.0;

        // inner cylinder: tangential velocity 1 (i.e. angular velocity 1/r_i), outer cylinder at rest
        const double Omega_i = 1.0 / r_i, Omega_o = 0.0;

        // viscosity 1, gap width 1, wall velocity 1, i.e. the Reynolds number is equal to the density
        const double mu = 1.0;
        const double rho = 1.0;

        static readonly string innerWallTag = IncompressibleBcType.Velocity_Inlet.ToString() + "_inner";
        static readonly string outerWallTag = IncompressibleBcType.Wall.ToString() + "_outer";

        /// <summary>
        /// Control object: single phase (species A everywhere), steady.
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
        /// <param name="resolution">grid refinement level, see <see cref="CurvedWedgeHelpers.CreateGrid"/></param>
        /// <param name="convection">false: Stokes, true: Navier-Stokes</param>
        static XNSE_Control BuildControl(int degree, int resolution, bool convection) {
            var C = new XNSE_Control();
            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "CurvedElementTests/TaylorCouette";

            C.SetDGdegree(degree);

            C.GridFunc = delegate () {
                var grd = CurvedWedgeHelpers.CreateGrid(resolution, FullRing);
                grd.DefineEdgeTags(delegate (double[] X) {
                    // X is the edge midpoint, which is in general not an interpolation node (e.g. Square_16):
                    // it lies on the polynomial approximation of the arc, not exactly on the circle,
                    // so assign each wall edge to the closer cylinder instead of testing r == r_i, r == r_o.
                    double r = Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
                    double tol = 1.0e-3 * (r_o - r_i);
                    if(Math.Abs(r - r_i) <= tol)
                        return innerWallTag;
                    if(Math.Abs(r - r_o) <= tol)
                        return outerWallTag;
                    throw new ArgumentOutOfRangeException($"unexpected boundary point ({X[0]}, {X[1]}), r = {r}");
                });
                return grd;
            };

            // rotating inner cylinder: u = u_theta * (-y/r, x/r), with u_theta = Omega_i * r_i
            double uWall = Omega_i * r_i;
            foreach(string spc in new[] { "A", "B" }) {
                C.AddBoundaryValue(innerWallTag, $"VelocityX#{spc}", $"X => -{uWall} * X[1] / Math.Sqrt(X[0]*X[0] + X[1]*X[1])", false);
                C.AddBoundaryValue(innerWallTag, $"VelocityY#{spc}", $"X => +{uWall} * X[0] / Math.Sqrt(X[0]*X[0] + X[1]*X[1])", false);
                C.AddBoundaryValue(outerWallTag);
            }

            // single phase: level set negative everywhere, i.e. only species A, no cut cells
            C.AddInitialValue("Phi", -1.0);
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.PhysicalParameters.rho_A = rho;
            C.PhysicalParameters.rho_B = rho;
            C.PhysicalParameters.mu_A = mu;
            C.PhysicalParameters.mu_B = mu;
            C.PhysicalParameters.Sigma = 0.0;
            C.PhysicalParameters.IncludeConvection = convection;
            C.PhysicalParameters.Material = true;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;

            return C;
        }

        /// <summary>
        /// Converts the species-A part of an XDG field of the solver into a standard DG field.
        /// </summary>
        static SinglePhaseField SpeciesA(XNSE solver, string name) {
            var xf = (XDGField)solver.CurrentState.Fields.Single(f => f.Identification == name);
            var u = new SinglePhaseField(xf.Basis.NonX_Basis, name);
            u.AccLaidBack(1.0, xf.GetSpeciesShadowField("A"));
            return u;
        }

        /// <summary>
        /// L2 error of <paramref name="p"/> against <paramref name="exact"/>, ignoring a constant offset (e.g. for a pressure which is only defined up to a constant).
        /// </summary>
        /// <remarks>
        /// Note: <see cref="DGField.L2ErrorNoMean"/> subtracts the integral of the difference, not its mean value, i.e. it is only correct for a domain of measure 1.
        /// </remarks>
        static double L2ErrorIgnoringMean(DGField p, ScalarFunction exact, int order) {
            var gdat = p.GridDat;
            var one = new SinglePhaseField(new Basis(gdat, 0));
            one.AccConstant(1.0);

            // the L2 projection preserves the integral, since constants are contained in the DG space
            var pExact = new SinglePhaseField(p.Basis);
            pExact.ProjectField(exact);
            double meanDiff = (p.IntegralOver(null) - pExact.IntegralOver(null)) / one.IntegralOver(null);

            var q = p.CloneAs();
            q.AccConstant(-meanDiff);
            return q.L2Error(exact, order);
        }

        /// <summary>
        /// Same self-consistency check as in <see cref="NonlinearSolver"/> (see <see cref="OperatorSelfConsistencyCheck"/>):
        /// \f$ \| M(u_0) u_0 + b(u_0) - f(u_0) \| \f$, with the linearization (M, b) and the evaluation f
        /// assembled by the timestepper of the solver, at the state \f$ u_0 = \textrm{scale} \cdot u \f$.
        /// The check is done here explicitly, since for a linear operator (Stokes),
        /// the linear solver is used and the check in <see cref="NonlinearSolver"/> is never executed.
        /// </summary>
        static (double Absolute, double Denominator, double NormF) LinearizationMismatch(XNSE solver, double scale) {
            // `AssembleMatrixCallback` is `internal protected`, therefore reflection
            var tsb = solver.Timestepping.TimesteppingBase;
            var mi = typeof(XdgTimesteppingBase).GetMethod("AssembleMatrixCallback", BindingFlags.Instance | BindingFlags.NonPublic);
            var assemble = (OperatorEvalOrLin)Delegate.CreateDelegate(typeof(OperatorEvalOrLin), tsb, mi);

            DGField[] u0 = solver.CurrentState.Fields.Select(f => {
                var g = f.CloneAs();
                g.Scale(scale);
                return g;
            }).ToArray();

            // note: the callback copies its argument into the current state of the solver, therefore backup and restore
            double[] backup = solver.CurrentStateVector.ToArray();
            BlockMsrMatrix M;
            double[] b, f;
            try {
                assemble(out M, out b, out _, u0, true, out _);
                assemble(out _, out f, out _, u0, false, out _);
            } finally {
                solver.CurrentStateVector.SetV(backup);
            }

            var check = b.CloneAs();
            M.SpMV(1.0, new CoordinateVector(u0), 1.0, check);
            var err = check.CloneAs();
            err.AccV(-1.0, f);

            // same comparison denominator as in NonlinearSolver.EvaluateOperator
            double normF = f.MPI_L2Norm();
            double comp = Math.Sqrt(Math.Max(normF, check.MPI_L2Norm()) * BLAS.MachineEps + BLAS.MachineEps);
            return (err.MPI_L2Norm(), comp, normF);
        }

        /// <summary>
        /// Steady Stokes flow between a rotating inner cylinder (tangential velocity 1) and an outer cylinder at rest,
        /// on the curved, rotationally periodic wedge.
        /// The exact velocity is the laminar Taylor-Couette profile; for Stokes flow (no convection), the exact pressure is constant.
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
        /// <param name="resolution">grid refinement level, see <see cref="CurvedWedgeHelpers.CreateGrid"/></param>
        [Test]
        public static void SteadyStokesOnPeriodicWedge([Values(1, 2, 3)] int degree, [Values(1, 2)] int resolution) {
            SteadyTaylorCouette(degree, resolution, convection: false);
        }

        /// <summary>
        /// Same as <see cref="SteadyStokesOnPeriodicWedge"/>, but for the Navier-Stokes equations:
        /// the velocity profile is the same, the exact pressure is \f$ p(r) \f$ with \f$ \partial_r p = \rho u_\theta^2 / r \f$.
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
        /// <param name="resolution">grid refinement level, see <see cref="CurvedWedgeHelpers.CreateGrid"/></param>
        [Test]
        public static void SteadyNavierStokesOnPeriodicWedge([Values(1, 2, 3)] int degree, [Values(1, 2)] int resolution) {
            SteadyTaylorCouette(degree, resolution, convection: true);
        }

        /// <summary>
        /// Common implementation of <see cref="SteadyStokesOnPeriodicWedge"/> and <see cref="SteadyNavierStokesOnPeriodicWedge"/>.
        /// </summary>
        static void SteadyTaylorCouette(int degree, int resolution, bool convection) {
            var C = BuildControl(degree, resolution, convection);
            string tag = $"{(convection ? "Navier-Stokes" : "Stokes")}, degree {degree}, resolution {resolution}";

            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();

                // linearization vs. evaluation, at the solution u and at 2*u (where f is clearly non-zero);
                // in BoSSS, the linearization around u0 is defined such that M(u0)*u0 + b(u0) == f(u0) for any state u0
                // number of periodic edges whose two cells live on different MPI processes (relevant for MPI-parallel runs)
                {
                    var gdat = solver.GridData;
                    int Jup = gdat.iGeomCells.NoOfLocalUpdatedCells;
                    int periodicEdges = 0, periodicEdgesAcrossProcs = 0;
                    for(int e = 0; e < gdat.iGeomEdges.Count; e++) {
                        if(gdat.iGeomEdges.EdgeTags[e] < GridCommons.FIRST_PERIODIC_BC_TAG)
                            continue;
                        periodicEdges++;
                        if(gdat.iGeomEdges.CellIndices[e, 1] >= Jup)
                            periodicEdgesAcrossProcs++;
                    }
                    Console.WriteLine($"{tag}: {gdat.MpiSize} MPI process(es), periodic edges (local, summed over processes) = {periodicEdges.MPISum()}, of which across process boundaries = {periodicEdgesAcrossProcs.MPISum()}");
                }

                double mismatchAtU = double.NaN, mismatchAt2u = double.NaN;
                foreach(double scale in new[] { 1.0, 2.0 }) {
                    var m = LinearizationMismatch(solver, scale);
                    Console.WriteLine($"{tag}, linearization mismatch at {scale}*u: absolute = {m.Absolute:0.####e-00}, relative = {m.Absolute / m.Denominator:0.####e-00} (comparison denominator {m.Denominator:0.####e-00}, |f(u)| = {m.NormF:0.####e-00})");
                    if(scale == 1.0)
                        mismatchAtU = m.Absolute / m.Denominator; // |f(u)| ~ 0 at the solution, so use the NonlinearSolver's own measure
                    else
                        mismatchAt2u = m.Absolute / m.NormF;
                }

                var u = new SinglePhaseField[2];
                var velErr = new double[2];
                var projErr = new double[2];
                for(int d = 0; d < 2; d++) {
                    u[d] = SpeciesA(solver, d == 0 ? VariableNames.VelocityX : VariableNames.VelocityY);
                    var exact = CurvedWedgeHelpers.TaylorCouetteVelocity(d, r_i, r_o, Omega_i, Omega_o);
                    velErr[d] = u[d].L2Error(exact, 2 * degree + 4);

                    // best approximation in the same DG space (lower bound for the error of any discrete solution)
                    var uProj = new SinglePhaseField(u[d].Basis);
                    uProj.ProjectField(exact);
                    projErr[d] = uProj.L2Error(exact, 2 * degree + 4);

                    Console.WriteLine($"{tag}, {u[d].Identification}: L2 error solution = {velErr[d]:0.####e-00}, L2 error projection = {projErr[d]:0.####e-00}, ratio = {velErr[d] / projErr[d]:0.##}");
                }

                // pressure: only defined up to a constant (walls everywhere), i.e. the mean is ignored
                var p = SpeciesA(solver, VariableNames.Pressure);
                var pExact = CurvedWedgeHelpers.TaylorCouettePressure(convection ? rho : 0.0, r_i, r_o, Omega_i, Omega_o);
                double presErr = L2ErrorIgnoringMean(p, pExact, 2 * degree + 4);
                Console.WriteLine($"{tag}, Pressure: L2 error (ignoring the mean) = {presErr:0.####e-00}");

                // continuity of the solution across the periodic boundary (if there is one)
                var j = CurvedWedgeHelpers.VelocityJumps(u, 2 * degree + 2);
                bool hasPeriodicEdges = j.LengthPeriodic > 0;
                double rmsPeriodic = hasPeriodicEdges ? j.Periodic / Math.Sqrt(j.LengthPeriodic) : double.NaN;
                double rmsInterior = j.Interior / Math.Sqrt(j.LengthInterior);
                Console.WriteLine($"{tag}: RMS velocity jump per edge length: periodic (rotated) = {(hasPeriodicEdges ? rmsPeriodic.ToString("0.####e-00") : "none")}, interior = {rmsInterior:0.####e-00}");

                // TODO: thresholds are roughly 2x the errors measured on the resolution-1 grid; replace e.g. by a convergence-rate check
                double velThreshold = degree switch { 1 => 7.0e-2, 2 => 3.0e-3, 3 => 1.0e-4, _ => throw new NotSupportedException() };
                for(int d = 0; d < 2; d++) {
                    Assert.Less(velErr[d], velThreshold, $"L2 error of {u[d].Identification} too large.");

                    // best approximation; measured: at most ~3.2e-3
                    Assert.Less(projErr[d], 1.0e-2, $"L2 error of the projection of the exact {u[d].Identification} too large.");

                    // the discrete solution should stay within a fixed factor of the best approximation; measured: ~3.5 to ~13.8
                    Assert.LessOrEqual(velErr[d] / projErr[d], 20.0, $"L2 error of {u[d].Identification} too large compared to the best approximation.");
                }
                Assert.Less(presErr, 1.0e-1, "L2 error of pressure too large.");

                // linearization at the solution, with the NonlinearSolver's own warning threshold (measure < 1);
                // measured: at most ~8e-3 (before the FixOutSelfDiagonalBlock fix: up to ~690 for Navier-Stokes)
                Assert.Less(mismatchAtU, 1.0, "linearization M(u0)*u0 + b(u0) is inconsistent with the evaluation f(u0) at the solution.");
                if(hasPeriodicEdges)
                    Assert.LessOrEqual(rmsPeriodic, 10 * rmsInterior, "velocity jump across periodic edges is too large compared to interior edges.");

                // Regression check for the consistency of the linearization at rotated periodic edges
                // (LECEdgeQuadrature2.FixOutSelfDiagonalBlock): at 2*u, |M*u0 + b - f(u0)| / |f(u0)| is ~1e-15 (round-off);
                // before that fix, it was between ~6e-13 and ~7e-9 for Navier-Stokes on this wedge (all degrees and resolutions).
                Assert.Less(mismatchAt2u, 1.0e-13, "linearization M(u0)*u0 + b(u0) is inconsistent with the evaluation f(u0) at 2*u.");
            }
        }
    }
}
