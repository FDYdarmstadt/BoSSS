﻿/* =======================================================================
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
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System.Linq;
using System.Collections.Generic;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Tests whether the XNSE solver (<see cref="XNSE_SolverMain"/>) also works MPI-parallel for 
    /// non-trivial cases.
    /// </summary>
    [TestFixture]
    public class XNSE_Solver_MPItest {

        [Test]
        static public void ParallelRisingDroplet([Values(3)] int p) {
            var C = RisingBubble(p: p);
            C.TracingNamespaces = "*";

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;

            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        [Test]
        static public void RotCube_GetSpeciesIDError() {
            // Tritt nur mit 4 cores auf !!!
            // Fixed: Bei AMR wird LevelsetTracker "genullt", dieser wurde bis dato noch vollständig an die Flüsse übergeben
            var C = Rotating_Cube(4, 30, 2, true);

            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        [Test]
        static public void RotCube_DomainDecompoitionError() {
            // error occurs only with 3 processors and AMR true

            var C = HardcodedControl.RotCubeDomainDecompoitionError();

            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        [Test]
        static public void RotCube_OrderNotSupportedInHMF() {
            // in 3D partially HMF was executed for surface integration, although Saye was selected 
            // For Saye the order is raised internally, to get same accuracy as with HMF
            // This higher order may not be supported by HMF
            // Fixed: in 3D now Saye is used everywhere

            var C = Rotating_Cube(4, 10, 3, false);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        [Test]
        public static void BadInitiallyDistributionTest(
            [Values(true,false)] bool useAMR) {
            var C = Rotating_Cube(k: 1, Res: 10, SpaceDim: 3, useAMR, useLoadBal: true , UsePredefPartitioning: true);
            // dbg_launch();
            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        [Test]
        public static void SayeBug() {
            // fallback method (which was not implemented at this point), would have been triggered in Saye quadrature
            var C = Rotating_Cube(1, 20, 3, true, false);
            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        //[Test]
        //public static void PardisoFailsInProjection() {
        //    // 4 cores
        //    var C = Rotating_Sphere(1, 10, 3, false, false, false);
        //    using (var solver = new XNSE()) {
        //        solver.Init(C);
        //        solver.RunSolverMode();
        //    }
        //}

        [Test]
        public static void EmptyMaskInSchwarz() {
            //--test=BoSSS.Application.XNSE_Solver.XNSE_Solver_MPItest.EmptyMaskInSchwarz

            // This test simulates bad initial distribution of void cells over ranks
            // which would lead to an error within Schwarz solver
            // because of voidcells Schwarzblocks would be empty
            // Remedy: force repartitioning at startup and fallback in schwarz if only some blocks are empty 
            var C = PartlyCoverdDomain(2, 50, 2, false, true, false);
            C.LinearSolver = new Solution.AdvancedSolvers.OrthoMGSchwarzConfig() {
                TargetBlockSize = 1000
            };
            C.GridPartType = GridPartType.clusterHilbert;
            C.DynamicLoadBalancing_On = true;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0;

            //this test takes too much time with 3 procs and exceed the 4 hr limit.
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;

            C.NoOfTimesteps = 50;
            using (var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        static void Main(string[] args) {
            
            BoSSS.Solution.Application.InitMPI();
            //Debugger.Launch();
            //ParallelRisingDroplet(1);
            //ParallelRisingDroplet(2);
            //ParallelRisingDroplet(3);
            //BoSSS.Application.XNSE_Solver.XNSE_Solver_MPItest.BadInitiallyDistributionTest(true);
            //BoSSS.Application.XNSE_Solver.XNSE_Solver_MPItest.RotCube_OrderNotSupportedInHMF();
            BoSSS.Application.XNSE_Solver.XNSE_Solver_MPItest.EmptyMaskInSchwarz();
            BoSSS.Solution.Application.FinalizeMPI();            

        }




        /// <summary>
        /// Configuration which performs three timesteps of the rising droplet;
        /// uses a pre-defined partitioning for 4 processors.
        /// </summary>
        public static XNSE_Control RisingBubble(int p = 3, int kelem = 20, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath =  @" C.NoOfTimesteps = 10;";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;//@"\\fdyprime\userspace\nietz\databases\big_test";
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Bubble";
            C.ProjectDescription = "rising bubble";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            /*C.FieldOptions.Add("FilteredVelocityX", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("FilteredVelocityY", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });*/
            C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Testcase 1");
            C.PhysicalParameters.rho_A = 100;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 10;
            C.PhysicalParameters.Sigma = 24.5; //2.0*




            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // grid generation
            // ===============
            #region grid


            double xSize = 1.0;
            double ySize = 2.0;


            //int kelem = 40;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1); //ohne 2*
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                //grd.EdgeTagNames.Add(3, "wall_left");
                //grd.EdgeTagNames.Add(4, "wall_right");
                //grd.EdgeTagNames.Add(3, "freeslip_left");
                //grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                grd.AddPredefinedPartitioning("ZwoProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if(x < 0.5)
                        rank = 0;
                    else
                        rank = 1;

                    return rank;
                });

                grd.AddPredefinedPartitioning("VierProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if(x < 0.35)
                        rank = 0;
                    else if(x < 0.5)
                        rank = 1;
                    else if(x < 0.75)
                        rank = 2;
                    else
                        rank = 3;

                    return rank;
                });


                return grd;
            };

            int MpiSize;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MpiSize);
            switch(MpiSize) {
                case 4:
                C.GridPartType = GridPartType.Predefined;
                C.GridPartOptions = "VierProcSplit";
                break;

                case 2:
                C.GridPartType = GridPartType.Predefined;
                C.GridPartOptions = "ZwoProcSplit";
                break;

                default:
                C.GridPartType = GridPartType.METIS;
                break;
            }

            #endregion



            // Initial Values
            // ==============
            #region init

            double[] center = new double[] { 0.5, 0.5 }; //0.5,0.5
            double radius = 0.25;


            //Func<double[], double> PhiFunc = (X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2()); // quadratic form
            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius); // signed-distance form
            C.InitSignedDistance = false;

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double> PeriodicFunc = x => radius;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("cd1a6c18-2659-4405-bf56-1e461441c0a0");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");


            C.AddBoundaryValue("wall_lower", VariableNames.LevelSet, PhiFunc);

            #endregion

            // Level-Set
            // =================
            #region Fourier

            int numSp = 1024;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for(int sp = 0; sp < numSp; sp++) {
                FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                samplP[sp] = radius;
            }

            //double circum = 2.0 * Math.PI * radius;
            //double filter = (circum * 20.0) / ((double)numSp / 2.0);
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem)
            //{
            //    center = center,
            //    FourierEvolve = Fourier_Evolution.FourierPoints,
            //    //Timestepper = FourierLevelSet_Timestepper.AdamsBashforth2,
            //    //UnderRelax = underrelax
            //    centerMove = CenterMovement.Reconstructed,
            //    //curvComp_extended = false
            //};


            #endregion

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching; // Test with old XNSE_SolverMain was using ScalarConvection

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;


            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            //C.AdvancedDiscretizationOptions.surfTensionMode = SurfaceTensionMode.Curvature_Fourier;
            //C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            //C.Timestepper_MassMatrix = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            //C.TimeStepper = XNSE_Control._Timestepper.BDF2;
            //double dt = 75e-4; // (1.0 / (double)kelem) / 16.0;
            //CFL condition > capillary time-step constraint
            double dt = 0.25 * Math.Sqrt((C.PhysicalParameters.rho_A + C.PhysicalParameters.rho_B) * Math.Pow((Math.Min(xSize, ySize) / kelem), 3) / (4 * Math.PI * C.PhysicalParameters.Sigma));

            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 3;
            C.saveperiod = 1;



            #endregion

            return C;

        }

        public static XNSE_Control PartlyCoverdDomain(int k = 4, int Res = 30, int SpaceDim = 2, bool useAMR = true, bool useLoadBal = false, bool UsePredefPartitioning = false) {
            return Rotating_Something(k, Res, SpaceDim, useAMR, useLoadBal, Geometry.Parted, UsePredefPartitioning);
        }

        public static XNSE_Control Rotating_Cube(int k = 4, int Res = 30, int SpaceDim = 2, bool useAMR = true, bool useLoadBal = false, bool UsePredefPartitioning = false) {
            return Rotating_Something(k, Res, SpaceDim, useAMR, useLoadBal, Geometry.Cube, UsePredefPartitioning);
        }

        public static XNSE_Control Rotating_Sphere(int k = 4, int Res = 30, int SpaceDim = 2, bool useAMR = true, bool useLoadBal = false, bool UsePredefPartitioning = false) {
            return Rotating_Something(k, Res, SpaceDim, useAMR, useLoadBal, Geometry.Sphere, UsePredefPartitioning);
        }

        enum Geometry {
            Cube = 0,
            Sphere = 1,
            Parted = 2
        }

        public static Func<IGrid> GridFuncFactory(int SpaceDim, int Res, bool UsePredefPartitioning) {
            double xMin = -1, yMin = -1, zMin = -1;
            double xMax = 1, yMax = 1, zMax = 1;

            // Predefined Partitioning
            Func<double[], int> MakeDebugPart = delegate (double[] X) {
                double x = X[0];
                double range = xMax - xMin;
                double interval = range * 0.9 / (ilPSP.Environment.MPIEnv.MPI_Size-1); // last rank gets only a 10% stripe of the domain
                return (int)((x - xMin) / interval);
            };

            // The Grid Function
            return delegate {
                
                var _xNodes = GenericBlas.Linspace(xMin, xMax, Res + 1);
                var _yNodes = GenericBlas.Linspace(yMin, yMax, Res + 1);
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
                if(UsePredefPartitioning) grd.AddPredefinedPartitioning("testgrid", MakeDebugPart);
                grd.EdgeTagNames.Add(2, "Wall");
                grd.DefineEdgeTags(delegate (double[] _X) {
                    return 2;
                });
                return grd;
            };
        }

        private static Func<double[], double, double> GenPhiFunc(double anglev, Geometry Gshape, int SpaceDim, double particleRad, double[] pos) {
            return delegate (double[] X, double t) {
                double angle = -(anglev * t) % (2 * Math.PI);
                switch (Gshape) {
                    case Geometry.Cube:
                    switch (SpaceDim) {
                        case 2:
                        return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                            Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
                            + particleRad;
                        case 3:
                        return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                                                Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)), Math.Abs(X[2] - pos[2])))
                                                + particleRad;
                        default:
                        throw new NotImplementedException("Dimension not supported");
                    };
                    case Geometry.Sphere:
                    switch (SpaceDim) {
                        case 2:
                        return -X[0] * X[0] - X[1] * X[1] + particleRad * particleRad;
                        case 3:
                        return -X[0] * X[0] - X[1] * X[1] - X[2] * X[2] + particleRad * particleRad;
                        default:
                        throw new NotImplementedException("Dimension not supported");
                    }
                    case Geometry.Parted:
                        return -X[0] + 0.7;
                    default:
                    throw new NotImplementedException("Shape unknown");
                }
            };
        }

        private static XNSE_Control Rotating_Something(int k, int Res, int SpaceDim, bool useAMR, bool useLoadBal, Geometry Gshape, bool UsePredefPartitioning) {
            XNSE_Control C = new XNSE_Control();
            // basic database options
            // ======================

            C.savetodb = false;
            //C.DbPath = @"D:\trash_db";
            C.ProjectName = "XNSE/IBM_test";
            C.ProjectDescription = "rotating cube";
            C.Tags.Add("rotating");
            C.Tags.Add("level set");
            C.Tags.Add(String.Format("{0}D", SpaceDim));

            // DG degrees
            // ==========

            C.SetFieldOptions(k, Math.Max(2, k * 2));
            if(UsePredefPartitioning) {
                C.GridPartType = GridPartType.Predefined;
                C.GridPartOptions = "testgrid";
            } else
                C.GridPartType = GridPartType.clusterHilbert;

            C.SessionName = "XNSE_rotcube_test";
            C.saveperiod = 1;


            // grid and boundary conditions
            // ============================

            //// Create Grid
            C.GridFunc = GridFuncFactory(SpaceDim, Res, UsePredefPartitioning);

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

            var PhiFunc = GenPhiFunc(anglev, Gshape, SpaceDim, particleRad, pos);

            Func<double[], double, double[]> VelocityAtIB = delegate (double[] X, double time) {

                if(pos.Length != X.Length)
                    throw new ArgumentException("check dimension of center of mass");

                Vector angVelo = new Vector(new double[] { 0, 0, anglev });
                Vector CenterofMass = new Vector(pos);
                Vector radialVector = new Vector(X) - CenterofMass;
                Vector transVelocity = new Vector(new double[SpaceDim]);
                Vector pointVelocity;

                switch(SpaceDim) {
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
            if(C.UseImmersedBoundary) {
                //C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), PhiFunc);
                C.InitialValues_EvaluatorsVec.Add(VariableNames.LevelSetCGidx(1), PhiFuncDelegate);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityX@Phi2", VelocityX);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityY@Phi2", VelocityY);
                if(SpaceDim == 3)
                    C.InitialValues_Evaluators_TimeDep.Add("VelocityZ@Phi2", VelocityZ);
            }
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            C.AddBoundaryValue("Wall");

            // misc. solver options
            // ====================

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LinearSolver = new Solution.AdvancedSolvers.OrthoMGSchwarzConfig() {
                NoOfMultigridLevels = 5,
                ConvergenceCriterion = 1E-8,
                MaxSolverIterations = 100,
                //MaxKrylovDim = 30,
                TargetBlockSize = 10000,
                //verbose = true
            };
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
                C.AMR_startUpSweeps = 1;
            }

            C.DynamicLoadBalancing_On = useLoadBal;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            //C.DynamicLoadBalancing_CellCostEstimators = Loadbalancing.XNSECellCostEstimator.Factory().ToList();
            C.DynamicLoadBalancing_ImbalanceThreshold = 0;

            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            double dt = 0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.dtFixed = dt;
            C.NoOfTimesteps = 1;

            return C;
        }

        public Func<double[], double, double> GetPhi() {
            throw new NotImplementedException();
        }

        public GridCommons CreateGrid(int Resolution) {
            throw new NotImplementedException();
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            throw new NotImplementedException();
        }
    }
}
