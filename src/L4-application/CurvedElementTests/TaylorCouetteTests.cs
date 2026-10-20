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
    /// Steady laminar Taylor-Couette flow on the curved grids of <see cref="CurvedWedgeHelpers.CreateGrid"/> and <see cref="CurvedWedgeHelpers.CreateGrid3D"/>
    /// </summary>
    [TestFixture]
    public static class TaylorCouetteTests {

        [OneTimeSetUp]
        public static void SetUp() {
            BoSSS.Solution.Application.InitMPI();
        }

        /// <summary>
        /// false: rotationally periodic wedge; true: full ring
        /// </summary>
        internal const bool FullRing = false;

        // must match the radial extent of the grid
        const double r_i = 1.0, r_o = 2.0;

        const double Omega_i = 1.0 / r_i, Omega_o = 0.0;

        const double mu = 1.0;
        const double rho = 1.0;

        static readonly string innerWallTag = IncompressibleBcType.Velocity_Inlet.ToString() + "_inner";
        static readonly string outerWallTag = IncompressibleBcType.Wall.ToString() + "_outer";

        /// <summary>
        /// Control object: single phase, steady.
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
        /// <param name="resolution">grid refinement level</param>
        /// <param name="convection">false: Stokes, true: Navier-Stokes</param>
        /// <param name="D">spatial dimension: 2 (wedge) or 3 (cylinder wedge, periodic in z)</param>
        static XNSE_Control BuildControl(int degree, int resolution, bool convection, int D) {
            var C = new XNSE_Control();
            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "CurvedElementTests/TaylorCouette";

            C.SetDGdegree(degree);

            C.GridFunc = delegate () {
                var grd = D == 3 ? CurvedWedgeHelpers.CreateGrid3D(resolution) : CurvedWedgeHelpers.CreateGrid(resolution, FullRing);
                grd.DefineEdgeTags(delegate (double[] X) {
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

            double uWall = Omega_i * r_i;
            foreach(string spc in new[] { "A", "B" }) {
                C.AddBoundaryValue(innerWallTag, $"VelocityX#{spc}", $"X => -{uWall} * X[1] / Math.Sqrt(X[0]*X[0] + X[1]*X[1])", false);
                C.AddBoundaryValue(innerWallTag, $"VelocityY#{spc}", $"X => +{uWall} * X[0] / Math.Sqrt(X[0]*X[0] + X[1]*X[1])", false);
                if(D == 3)
                    C.AddBoundaryValue(innerWallTag, $"VelocityZ#{spc}", "X => 0.0", false);
                C.AddBoundaryValue(outerWallTag);
            }

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
        /// L2 error of <paramref name="p"/> against <paramref name="exact"/>, ignoring a constant offset
        /// </summary>
        /// <remarks>
        /// <see cref="DGField.L2ErrorNoMean"/> is only correct for a domain of measure 1.
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
        /// \f$ \| M(u_0) u_0 + b(u_0) - f(u_0) \| \f$ at \f$ u_0 = \textrm{scale} \cdot u \f$, assembled by the timestepper
        /// (same check as in <see cref="NonlinearSolver"/>, see <see cref="OperatorSelfConsistencyCheck"/>)
        /// </summary>
        static (double Absolute, double Denominator, double NormF) LinearizationMismatch(XNSE solver, double scale) {
            var tsb = solver.Timestepping.TimesteppingBase;
            var mi = typeof(XdgTimesteppingBase).GetMethod("AssembleMatrixCallback", BindingFlags.Instance | BindingFlags.NonPublic);
            var assemble = (OperatorEvalOrLin)Delegate.CreateDelegate(typeof(OperatorEvalOrLin), tsb, mi);

            DGField[] u0 = solver.CurrentState.Fields.Select(f => {
                var g = f.CloneAs();
                g.Scale(scale);
                return g;
            }).ToArray();

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

            double normF = f.MPI_L2Norm();
            double comp = Math.Sqrt(Math.Max(normF, check.MPI_L2Norm()) * BLAS.MachineEps + BLAS.MachineEps);
            return (err.MPI_L2Norm(), comp, normF);
        }

        /// <summary>
        /// Steady Stokes flow
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
        /// <param name="resolution">grid refinement level, see <see cref="CurvedWedgeHelpers.CreateGrid"/></param>
        [Test]
        public static void SteadyStokesOnPeriodicWedge([Values(1, 2, 3)] int degree, [Values(1, 2)] int resolution) {
            SteadyTaylorCouette(degree, resolution, convection: false, D: 2);
        }

        /// <summary>
        /// Steady Navier-Stokes flow
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
        /// <param name="resolution">grid refinement level, see <see cref="CurvedWedgeHelpers.CreateGrid"/></param>
        [Test]
        public static void SteadyNavierStokesOnPeriodicWedge([Values(1, 2, 3)] int degree, [Values(1, 2)] int resolution) {
            SteadyTaylorCouette(degree, resolution, convection: true, D: 2);
        }

        /// <summary>
        /// Steady Stokes flow on the 3D cylinder wedge (periodic in z)
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
#if !DEBUG
        [Test]
#endif
        public static void SteadyStokesOnPeriodicCylinderWedge([Values(1, 2)] int degree) {
            SteadyTaylorCouette(degree, 1, convection: false, D: 3);
        }

        /// <summary>
        /// Steady Navier-Stokes flow on the 3D cylinder wedge (periodic in z)
        /// </summary>
        /// <param name="degree">DG degree of the velocity</param>
#if !DEBUG
        [Test]
#endif
        public static void SteadyNavierStokesOnPeriodicCylinderWedge([Values(1, 2)] int degree) {
            SteadyTaylorCouette(degree, 1, convection: true, D: 3);
        }

        /// <summary>
        /// Common implementation of the Stokes and Navier-Stokes tests, in 2D and 3D.
        /// </summary>
        static void SteadyTaylorCouette(int degree, int resolution, bool convection, int D) {
            var C = BuildControl(degree, resolution, convection, D);
            string tag = $"{(convection ? "Navier-Stokes" : "Stokes")}, {D}D, degree {degree}, resolution {resolution}";

            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();

                // number of periodic edges across MPI process boundaries
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

                // linearization vs. evaluation, at the solution u and at 2*u
                double mismatchAtU = double.NaN, mismatchAt2u = double.NaN;
                foreach(double scale in new[] { 1.0, 2.0 }) {
                    var m = LinearizationMismatch(solver, scale);
                    Console.WriteLine($"{tag}, linearization mismatch at {scale}*u: absolute = {m.Absolute:0.####e-00}, relative = {m.Absolute / m.Denominator:0.####e-00} (comparison denominator {m.Denominator:0.####e-00}, |f(u)| = {m.NormF:0.####e-00})");
                    if(scale == 1.0)
                        mismatchAtU = m.Absolute / m.Denominator; // |f(u)| ~ 0, therefore the NonlinearSolver's measure
                    else
                        mismatchAt2u = m.Absolute / m.NormF;
                }

                var u = new SinglePhaseField[D];
                var velErr = new double[D];
                var projErr = new double[D];
                for(int d = 0; d < D; d++) {
                    u[d] = SpeciesA(solver, VariableNames.Velocity_d(d));
                    var exact = CurvedWedgeHelpers.TaylorCouetteVelocity(d, r_i, r_o, Omega_i, Omega_o);
                    velErr[d] = u[d].L2Error(exact, 2 * degree + 4);

                    // best approximation in the same DG space
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

                double velThreshold = degree switch { 1 => 1.0e-1, 2 => 1.0e-2, 3 => 1.0e-4, _ => throw new NotSupportedException() };
                for(int d = 0; d < D; d++) {
                    Assert.Less(velErr[d], velThreshold, $"L2 error of {u[d].Identification} too large.");

                    Assert.Less(projErr[d], 1.0e-2, $"L2 error of the projection of the exact {u[d].Identification} too large.");

                    // ratio only for components with a non-trivial exact solution (u_z = 0 is represented exactly)
                    if(projErr[d] > 0)
                        Assert.LessOrEqual(velErr[d] / projErr[d], 20.0, $"L2 error of {u[d].Identification} too large compared to the best approximation.");
                }
                Assert.Less(presErr, 1.0e-1, "L2 error of pressure too large.");

                // linearization at the solution: NonlinearSolver's warning threshold
                Assert.Less(mismatchAtU, 1.0, "linearization M(u0)*u0 + b(u0) is inconsistent with the evaluation f(u0) at the solution.");
                if(hasPeriodicEdges)
                    Assert.LessOrEqual(rmsPeriodic, 10 * rmsInterior, "velocity jump across periodic edges is too large compared to interior edges.");

                // linearization at 2*u, relative to |f|
                Assert.Less(mismatchAt2u, 1.0e-13, "linearization M(u0)*u0 + b(u0) is inconsistent with the evaluation f(u0) at 2*u.");
            }
        }
    }
}
