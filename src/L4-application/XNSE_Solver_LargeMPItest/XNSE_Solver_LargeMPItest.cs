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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using ilPSP;
using System.Diagnostics;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Application.XNSE_Solver.Legacy;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Tests whether the XNSE solver (<see cref="XNSE_SolverMain"/>) also works MPI-parallel for larger cases (8 MPI cores)
    /// </summary>
    [TestFixture]
    public static class XNSE_Solver_LargeMPItest {

        [Test]
        static public void ParallelRotatingSphere() {
            var C = RotatingSphere(k:1, Res:20, SpaceDim:3, useAMR:true, loadbalancing:true);
            //C.TracingNamespaces = "*";

            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
            //Assert.IsTrue(false, "this should fail!");
        }



        
        public static XNSE_Control RotatingSphere(int k = 1, int Res = 20, int SpaceDim = 3, bool useAMR = true, bool loadbalancing = true) {
            XNSE_Control C = new XNSE_Control();
            // basic database options
            // ======================

            
            C.savetodb = false;
            C.ProjectName = "XNSE/IBM_benchmark";
            C.ProjectDescription = "rotating sphere";
            C.Tags.Add("rotating");
            C.Tags.Add("tracing");

            // DG degrees
            // ==========

            C.SetFieldOptions(k, Math.Max(6, k * 2));
            C.SessionName = "XNSE_rotsphere";
            C.saveperiod = 1;
            //C.TracingNamespaces = "*";
            
            // grid and boundary conditions
            // ============================

            //// Create Grid
            Console.WriteLine("...generating grid");
            double xMin = -1, yMin = -1, zMin = -1;
            double xMax = 1, yMax = 1, zMax = 1;

            Func<double[], int> MakeDebugPart = delegate (double[] X) {
                double x = X[0];

                double range = xMax - xMin;
                double interval = range / ilPSP.Environment.MPIEnv.MPI_Size;

                return (int)((x - xMin) / interval);
            };

            C.GridFunc = delegate {

                // x-direction
                
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
                    grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                    break;

                    default:
                    throw new ArgumentOutOfRangeException();
                }

                //grd.AddPredefinedPartitioning("debug", MakeDebugPart);

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
            //C.GridPartType = GridPartType.Predefined;
            //C.GridPartOptions = "debug";
            C.GridPartType = GridPartType.clusterHilbert;

            C.DynamicLoadBalancing_On = loadbalancing;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_CellCostEstimatorFactories = Loadbalancing.XNSECellCostEstimator.Factory().ToList();
            C.DynamicLoadBalancing_ImbalanceThreshold = -0.1;

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
            double particleRad = 0.261;



            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                double power = 10;
                //anglev *= t < 0.005 ? Math.Sin(2000 * Math.PI * t - Math.PI / 2) / 2 + 0.5 : 1;
                double angle = -(anglev * t) % (2 * Math.PI);
                switch (SpaceDim) {
                    case 2:
                    // Inf-Norm square
                    return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                        Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
                        + particleRad;

                    // p-Norm square
                    //return -Math.Pow((Math.Pow((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle), power)
                    //+ Math.Pow((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle), power)), 1.0/power)
                    //+ particleRad; // 1e6

                    // 0-Norm square
                    //return -Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle))
                    //- Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle))
                    //+ Math.Abs(particleRad);

                    // circle
                    //return -X[0] * X[0] - X[1] * X[1] + particleRad * particleRad;

                    case 3:
                    // Inf-Norm cube
                    //return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                    //                        Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)),
                    //                        Math.Abs(X[2] - pos[2])))
                    //                        + particleRad;

                    // p-Norm cube
                    //return -Math.Pow(Math.Pow((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle), power)
                    //+ Math.Pow((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle), power)
                    //+ Math.Pow(X[2] - pos[2], power),1.0/power)
                    //+ particleRad;

                    // 0-Norm cube
                    //return -Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle))
                    //- Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle))
                    //- Math.Abs(X[2] - pos[2])
                    //+ Math.Abs(particleRad);

                    // sphere
                    return -X[0] * X[0] - X[1] * X[1] - X[2] * X[2] + particleRad * particleRad;

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
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            C.UseSchurBlockPrec = true;
            //C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            //C.PressureBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LinearSolver.NoOfMultigridLevels = 5;
            C.LinearSolver.ConvergenceCriterion = 1E-8;
            C.LinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MaxKrylovDim = 30;
            C.LinearSolver.TargetBlockSize = 10000;
            C.LinearSolver.verbose = true;
            C.LinearSolver.SolverCode = LinearSolverCode.exp_Kcycle_schwarz;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.verbose = true;
            
            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 2 });
                C.AMR_startUpSweeps = 1;
            }

            // Timestepping
            // ============


            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            double dt = 0.01;
            //C.dtMax = dt;
            //C.dtMin = dt*1E-2;
            C.dtFixed = dt;
            C.NoOfTimesteps = 3;

            // haben fertig...
            // ===============

            return C;

        }



  
        /// <summary>
        /// 
        /// </summary>
        static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI();
            ParallelRotatingSphere();
            BoSSS.Solution.Application.FinalizeMPI();
        }




    }
}
