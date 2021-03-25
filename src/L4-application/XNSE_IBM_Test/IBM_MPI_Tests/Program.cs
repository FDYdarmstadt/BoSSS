using System;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.XNSECommon;
using NUnit.Framework;
using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Application.XNSE_Solver;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.XdgTimestepping;


namespace IBM_MPI_Tests {
    public static class Program {
        static void Main(string[] args) {
        }

        [Test]
        public static void IBMChannelSolverTest(
            [Values(1, 2, 3)] int FlowSolverDegree = 2,
            [Values(0)] double angle = 0.0,
            [Values(LinearSolverCode.exp_Kcycle_schwarz, LinearSolverCode.exp_gmres_levelpmg)] LinearSolverCode solvercode = LinearSolverCode.exp_Kcycle_schwarz
            ) {
            double AgglomerationTreshold = 0.3;

            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;


            int GridResolution = 1;

            var Tst = new IBMChannel(30 * Math.PI / 180, true);

            var C = ASUnitTest.TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, ViscosityMode.Standard, SurfTensionMode: SurfaceStressTensor_IsotropicMode.Curvature_Projected, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolution, solvercode: solvercode);
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.ConvergenceCriterion = 1e-11;


            XNSESolverTest(Tst, C);
        }

        public static void RotatingCube(
            [Values(1, 2, 3)] int FlowSolverDegree = 2,
            [Values(0)] double angle = 0.0,
            [Values(LinearSolverCode.exp_Kcycle_schwarz, LinearSolverCode.exp_gmres_levelpmg)] LinearSolverCode solvercode = LinearSolverCode.exp_Kcycle_schwarz
            ) {
            double AgglomerationTreshold = 0.3;

            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;


            int GridResolution = 1;

            var Tst = new IBMChannel(30 * Math.PI / 180, true);

            var C = ASUnitTest.TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, ViscosityMode.Standard, SurfTensionMode: SurfaceStressTensor_IsotropicMode.Curvature_Projected, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolution, solvercode: solvercode);
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.ConvergenceCriterion = 1e-11;


            XNSESolverTest(Tst, C);
        }

        private static void XNSESolverTest(IXNSETest Tst, XNSE_Control C) {

            if (Tst.SpatialDimension == 3) {
                Console.WriteLine($"Reminder: skipping 3D test for now...");
                return;
            }

            using (var solver = new XNSE()) {

                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;

                solver.Init(C);
                solver.RunSolverMode();
                //if(C.TimesteppingMode == AppControl._TimesteppingMode.Steady) // deavtivated by FK; has only value for a series of meshes, but not for a single calc.
                //    solver.OperatorAnalysis();

                //-------------------Evaluate Error ---------------------------------------- 
                var evaluator = new XNSEErrorEvaluator<XNSE_Control>(solver);
                double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length; i++) {
                    bool ok = LastErrors[i] <= ErrThresh[i];
                    Console.Write("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);

                    if (ok)
                        Console.WriteLine("   (ok)");
                    else
                        Console.WriteLine("   Above Threshold (" + ErrThresh[i] + ")");
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    bool ok = ResNorms[i] <= ResThresh[i];
                    Console.Write("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);

                    if (ok)
                        Console.WriteLine("   (ok)");
                    else
                        Console.WriteLine("   Above Threshold (" + ResThresh[i] + ")");
                }

                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i], $"Error {solver.CurrentState.Fields[i].Identification} above threshold.");

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i], $"Residual {solver.CurrentResidual.Fields[i].Identification} above threshold.");
            }
        }

        public static XNSE_Control Rotating_Cube(int k = 2, int Res = 10, int SpaceDim = 2, bool useAMR = false) {
            XNSE_Control C = new XNSE_Control();
            // basic database options
            // ======================

            C.savetodb = false;
            C.ProjectName = "XNSE/IBM_benchmark";
            C.ProjectDescription = "rotating cube";
            C.Tags.Add("rotating");
            C.Tags.Add("tracing");

            // DG degrees
            // ==========

            C.SetFieldOptions(k, Math.Max(6, k * 2));
            C.GridPartType = GridPartType.Hilbert;
            C.SessionName = "XNSE_rotsphere";
            C.saveperiod = 1;
            //C.TracingNamespaces = "*";


            // grid and boundary conditions
            // ============================

            //// Create Grid
            Console.WriteLine("...generating grid");
            C.GridFunc = delegate {

                // x-direction
                double xMin = -1, yMin = -1, zMin = -1;
                double xMax = 1, yMax = 1, zMax = 1;
                var _xNodes = GenericBlas.Linspace(xMin, xMax, Res + 1);
                //var _xNodes = GenericBlas.Logspace(0, 3, cells_x + 1);
                // y-direction
                var _yNodes = GenericBlas.Linspace(yMin, yMax, Res + 1);
                // z-direction
                var _zNodes = GenericBlas.Linspace(zMin, zMax, Res + 1);


                GridCommons grd;
                switch (SpaceDim) {
                    case 2:
                    grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                    break;

                    case 3:
                    grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, BoSSS.Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                    break;

                    default:
                    throw new ArgumentOutOfRangeException();
                }

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x, y, z;
                    x = X[0];
                    y = X[1];
                    if (SpaceDim == 3)
                        z = X[2];

                    return 2;
                });

                return grd;

            };

            //// Set Initial Conditions
            //C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            //C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //if (SpaceDim == 3)
            //    C.InitialValues_Evaluators.Add("VelocityZ", X => 0);


            // Phi (X,t): p-norm cube with forced rotation

            // Physical Parameters
            // ===================
            const double rhoA = 1;
            const double muA = 1;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;
            double anglev = 10;
            double[] pos = new double[SpaceDim];
            double particleRad = 0.26;



            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                double power = 10;
                //anglev *= t < 0.005 ? Math.Sin(2000 * Math.PI * t - Math.PI / 2) / 2 + 0.5 : 1;
                double angle = -(anglev * t) % (2 * Math.PI);
                switch (SpaceDim) {
                    case 2:
                    //return -Math.Pow((Math.Pow((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle), power)
                    //+ Math.Pow((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle), power)), 1.0/power)
                    //+ particleRad; // 1e6

                    return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                        Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
                        + particleRad;

                    //return -Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle))
                    //- Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle))
                    //+ Math.Abs(particleRad);
                    //return -X[0] * X[0] - X[1] * X[1] + particleRad * particleRad;

                    case 3:
                    return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                                            Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)), Math.Abs(X[2] - pos[2])))
                                            + particleRad;

                    //return -Math.Pow(Math.Pow((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle), power)
                    //+ Math.Pow((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle), power)
                    //+ Math.Pow(X[2] - pos[2], power),1.0/power)
                    //+ particleRad;

                    //return -Math.Max(Math.Pow((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle), power))
                    //+ Math.Pow((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle), power)
                    //+ Math.Pow(X[2] - pos[2], power), 1.0 / power)
                    //+ particleRad;

                    //return -Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle))
                    //- Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle))
                    //- Math.Abs(X[2] - pos[2])
                    //+ Math.Abs(particleRad);

                    //return -X[0] * X[0] - X[1] * X[1] - X[2] * X[2] + particleRad * particleRad;
                    default:
                    throw new NotImplementedException();
                }
            };

            Func<double[], double, double[]> VelocityAtIB = delegate (double[] X, double time) {

                if (pos.Length != X.Length)
                    throw new ArgumentException("check dimension of center of mass");

                Vector angVelo = new Vector(new double[] { 0, 0, anglev });
                Vector CenterofMass = new Vector(pos);
                Vector radialVector = new Vector(X) - CenterofMass;
                Vector transVelocity = new Vector(new double[SpaceDim]);
                Vector pointVelocity;

                switch (SpaceDim) {
                    case 2:
                    pointVelocity = new Vector(transVelocity[0] - angVelo[2] * radialVector[1], transVelocity[1] + angVelo[2] * radialVector[0]);
                    break;
                    case 3:
                    pointVelocity = transVelocity + angVelo.CrossProduct(radialVector);
                    break;
                    default:
                    throw new NotImplementedException("this number of dimensions is not supported");
                }

                return pointVelocity;
            };

            Func<double[], double, double> VelocityX = delegate (double[] X, double time) { return VelocityAtIB(X, time)[0]; };
            Func<double[], double, double> VelocityY = delegate (double[] X, double time) { return VelocityAtIB(X, time)[1]; };
            Func<double[], double, double> VelocityZ = delegate (double[] X, double time) { return VelocityAtIB(X, time)[2]; };

            var PhiFuncDelegate = BoSSS.Solution.Utils.NonVectorizedScalarFunction.Vectorize(PhiFunc);

            C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), X => -1);
            C.UseImmersedBoundary = true;
            if (C.UseImmersedBoundary) {
                //C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), PhiFunc);
                C.InitialValues_EvaluatorsVec.Add(VariableNames.LevelSetCGidx(1), PhiFuncDelegate);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityX@Phi2", VelocityX);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityY@Phi2", VelocityY);
                if (SpaceDim == 3)
                    C.InitialValues_Evaluators_TimeDep.Add("VelocityZ@Phi2", VelocityZ);
            }
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            C.AddBoundaryValue("Wall");

            //C.OperatorMatrixAnalysis = false;

            // misc. solver options
            // ====================

            //C.EqualOrder = false;
            //C.PressureStabilizationFactor = 1;
            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            C.UseSchurBlockPrec = true;
            //C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            //C.PressureBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LinearSolver.NoOfMultigridLevels = 5;
            C.LinearSolver.ConvergenceCriterion = 1E-8;
            C.LinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MaxKrylovDim = 30;
            C.LinearSolver.TargetBlockSize = 100;
            C.LinearSolver.verbose = true;
            C.LinearSolver.SolverCode = LinearSolverCode.exp_Kcycle_schwarz;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
                C.AMR_startUpSweeps = 1;
            }

            // Timestepping
            // ============

            //C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            double dt = 0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.dtFixed = dt;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;

        }

    }
}
