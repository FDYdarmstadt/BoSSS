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

using BoSSS.Application.IBM_Solver;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.Multigrid;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IBM_Solver {
    public class HardcodedPrecTest {

        static public IBM_Control PrecTest3DChannel(int k, int cells_x, int cells_yz, int Procs) {
            IBM_Control C = new IBM_Control();

            // basic database options
            // ======================
            C.savetodb = true;
            C.DbPath = @"/home/oe11okuz/BoSSS_DB/Lichtenberg_DB";

            //string restartSession = "727da287-1b6a-463e-b7c9-7cc19093b5b3";
            //string restartGrid = "3f8f3445-46f1-47ed-ac0e-8f0260f64d8f";

            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication<AppControl> app, int noOfPerformanceClasses) {
                Console.WriteLine("i was called");
                int[] map = new int[] { 1, 5, 100 };
                return new StaticCellCostEstimator(map);
            });

            // Assign correct names
            C.SessionName = "Channel_" + k + "_" + cells_x + "x" + cells_yz + "_" + Procs + "Procs_4Deeep";

            C.saveperiod = 1;
            //C.SessionName = "Sphere_k" + k + "_h" + h+"Re100";
            C.ProjectName = "3DChannel";
            C.ProjectDescription = "Sphere_k" + k + cells_x + "x" + cells_yz + "x" + cells_yz;
            C.Tags.Add("Prec Test");

            // Create Fields
            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #region Creates grid () and sets BC
            //// Create Grid
            Console.WriteLine("...generating grid");
            C.GridFunc = delegate {

                // x-direction
                var _xNodes = GenericBlas.Linspace(-0.5, 1.5, cells_x + 1);

                // y-direction
                var _yNodes = GenericBlas.Linspace(-0.5, 0.5, cells_yz + 1);

                // z-direction
                var _zNodes = GenericBlas.Linspace(-0.5, 0.5, cells_yz + 1);

                // Cut Out
                var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, true, false, CellType.Cube_Linear);

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];
                    double z = X[2];

                    if (Math.Abs(x - (-0.5)) < 1.0e-6)
                        // inlet
                        return 1;

                    if (Math.Abs(x - (1.5)) < 1.0e-6)
                        // outlet
                        return 3;

                    if (Math.Abs(y - (-0.5)) < 1.0e-6)
                        // left
                        return 2;

                    if (Math.Abs(y - (0.5)) < 1.0e-6)
                        // right
                        return 2;

                    if (Math.Abs(z - (-0.5)) < 1.0e-6)
                        // top left
                        return 2;

                    if (Math.Abs(z - (0.5)) < 1.0e-6)
                        // top right
                        return 2;

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            #endregion



            // Set Initial Conditions
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            C.InitialValues_Evaluators.Add("VelocityZ", X => 0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0);

            // Because its only channeö
            C.InitialValues_Evaluators.Add("Phi", X => -1);

            Console.WriteLine("...starting calculation of Preconditioning test with 3D Channel");

            // Physical values
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1.0 / 10.0;

            // Boundary conditions
            C.AddBoundaryCondition("Velocity_inlet", "VelocityX", (X, t) => 1 - 4 * (X[2] * X[2]));
            C.AddBoundaryCondition("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryCondition("Wall");
            C.AddBoundaryCondition("Pressure_Outlet");


            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 1000;
            C.MaxSolverIterations = 1;
            C.Solver_ConvergenceCriterion = 1E-5;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NonlinearMethod = Solution.XdgTimestepping.NonlinearSolverMethod.Newton;


            // Choosing the Preconditioner
            ISolverSmootherTemplate Prec;

            //Prec = new SchurPrecond()
            //{
            //    SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
            //};

            //Prec = new SchurPrecond()
            //{
            //    SchurOpt = SchurPrecond.SchurOptions.SIMPLE
            //};

            //C.LinearSolver = new Schwarz() {
            //    m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
            //        NoOfParts = Procs,
            //    },
            //    Overlap = 1
            //};

            //C.LinearSolver = new ClassicMultigrid() {
            //    CoarserLevelSolver = new DirectSolver() {
            //        WhichSolver = DirectSolver._whichSolver.MUMPS

            //    },
            //};

            C.LinearSolver = new Schwarz() {
                m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                    Depth = 4,
                },
                CoarseSolver = new ClassicMultigrid() {
                    CoarserLevelSolver = new ClassicMultigrid() {
                        CoarserLevelSolver = new ClassicMultigrid() {
                            CoarserLevelSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS },
                        },
                    },
                },
                Overlap = 1,
            };




            //C.LinearSolver = new SoftGMRES()
            //{
            //    MaxKrylovDim = C.MaxKrylovDim,
            //    Precond = Prec,
            //    m_Tolerance = 1E-6,
            //    m_MaxIterations = 50
            //};       

            // Timestepping
            // ============
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1E20;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10000000;
            C.NoOfTimesteps = 1;
            C.NoOfMultigridLevels = 5;

            return C;
        }


        static public IBM_Control PrecTest3DChannelDegenhardt(int precNo, int k, int cells_x, int cells_yz, int re, int ASparts = 5, int ASDepth = 2, int MGLevels = 3, int maxKrDim = 1000, int saveToDB = 1) {
            IBM_Control C = new IBM_Control();
            bool name_newton = true;

            // basic database options
            // ======================
            if (saveToDB == 1) {
                C.savetodb = true;
            } else {
                C.savetodb = false;
            }
            //C.savetodb = saveToDB;
            //C.savetodb = true;
            C.DbPath = @"/work/scratch/ws35kire/performance_db";
            //C.DbPath = @"\\dc1\scratch\Stange\\HiWi\performance_db";


            //string restartSession = "727da287-1b6a-463e-b7c9-7cc19093b5b3";
            //string restartGrid = "3f8f3445-46f1-47ed-ac0e-8f0260f64d8f";

            C.DynamicLoadBalancing_Period = 1;
            //C.DynamicLoadBalancing_CellCostEstimatorFactory = delegate (IApplication<AppControl> app, int noOfPerformanceClasses) {
            //    Console.WriteLine("i was called");
            //    int[] map = new int[] { 1, 5, 100 };
            //    return new StaticCellCostEstimator(map);
            //};

            // Assign correct names
            if (name_newton) {
                C.SessionName = "Newton_Channel_prec" + precNo + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;
                C.ProjectDescription = "Newton_Sphere_k_prec" + precNo + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;

            } else {
                C.SessionName = "Channel_prec" + precNo + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;
                C.ProjectDescription = "Sphere_k_prec" + precNo + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;

            }

            C.saveperiod = 1;
            //C.SessionName = "Sphere_k" + k + "_h" + h+"Re100";
            C.ProjectName = "3DChannel";
            C.Tags.Add("Prec param study");

            // Create Fields
            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #region Creates grid () and sets BC
            //// Create Grid
            Console.WriteLine("...generating grid");
            C.GridFunc = delegate {

                // x-direction
                var _xNodes = GenericBlas.Linspace(-0.5, 1.5, cells_x + 1);

                // y-direction
                var _yNodes = GenericBlas.Linspace(-0.5, 0.5, cells_yz + 1);

                // z-direction
                var _zNodes = GenericBlas.Linspace(-0.5, 0.5, cells_yz + 1);

                // Cut Out
                var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, true, false, CellType.Cube_Linear);

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X)
                {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];
                    double z = X[2];

                    if (Math.Abs(x - (-0.5)) < 1.0e-6)
                        // inlet
                        return 1;

                    if (Math.Abs(x - (1.5)) < 1.0e-6)
                        // outlet
                        return 3;

                    if (Math.Abs(y - (-0.5)) < 1.0e-6)
                        // left
                        return 2;

                    if (Math.Abs(y - (0.5)) < 1.0e-6)
                        // right
                        return 2;

                    if (Math.Abs(z - (-0.5)) < 1.0e-6)
                        // top left
                        return 2;

                    if (Math.Abs(z - (0.5)) < 1.0e-6)
                        // top right
                        return 2;

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            #endregion



            // Set Initial Conditions
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            C.InitialValues_Evaluators.Add("VelocityZ", X => 0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0);

            // Because its a sphere

            // C.particleRadius = 0.1;
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + -(X[2]).Pow2() + C.particleRadius.Pow2());
            C.InitialValues_Evaluators.Add("Phi", X => -1);

            Console.WriteLine("...starting calculation of Preconditioning test with 3D Channel");
            if (name_newton) {
                Console.WriteLine("Newton_Channel_k_prec" + precNo + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim);

            } else {
                Console.WriteLine("Sphere_k_prec" + precNo + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim);

            }

            // Physical values
            C.PhysicalParameters.rho_A = 1;
            // 1/Re
            //C.PhysicalParameters.mu_A = 1.0 / 10.0;
            //C.PhysicalParameters.mu_A = 0.2 / re;

            C.PhysicalParameters.mu_A = 1.0 / re;

            // Boundary conditions
            C.AddBoundaryCondition("Velocity_inlet", "VelocityX", (X, t) => 1 - 4 * (X[2] * X[2]));
            C.AddBoundaryCondition("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryCondition("Wall");
            C.AddBoundaryCondition("Pressure_Outlet");


            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            //C.MaxKrylovDim = 1000;
            C.MaxKrylovDim = maxKrDim;
            C.MaxSolverIterations = 50;
            // C.MaxSolverIterations = 10000;
            C.Solver_ConvergenceCriterion = 1E-5;
            //C.Solver_ConvergenceCriterion = 1E-6;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NonlinearMethod = Solution.XdgTimestepping.NonlinearSolverMethod.Newton;

            // Choosing the Preconditioner
            ISolverSmootherTemplate Prec;



            switch (precNo) {
                case 0: {
                        Prec = null;
                        break;
                    }
                case 1: {
                        Prec = new SchurPrecond() {
                            SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
                        };
                        break;
                    }
                case 2: {
                        Prec = new SchurPrecond() {
                            SchurOpt = SchurPrecond.SchurOptions.SIMPLE
                        };
                        break;
                    }
                case 3: {
                        Prec = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //NoOfParts = 5,
                                NoOfParts = ASparts,
                            },
                            CoarseSolver = new DirectSolver() {
                                WhichSolver = DirectSolver._whichSolver.MUMPS
                            },
                            Overlap = 1
                        };
                        C.NoOfMultigridLevels = 2;
                        break;
                    }
                case 4: {
                        //C.NoOfMultigridLevels = 3;
                        C.NoOfMultigridLevels = MGLevels;
                        Prec = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                //Depth = 2,
                                Depth = ASDepth,
                            },
                            CoarseSolver = new DirectSolver() {
                                WhichSolver = DirectSolver._whichSolver.MUMPS
                            },
                            Overlap = 1
                        };
                        break;
                    }
                default: {
                        Prec = new SchurPrecond() {
                            SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
                        };
                        break;
                    }

            }

            //Prec = new SchurPrecond()
            //{
            //    SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
            //};

            //Prec = new SchurPrecond()
            //{
            //    SchurOpt = SchurPrecond.SchurOptions.SIMPLE
            //};

            //Prec = new Schwarz()
            //{
            //    m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
            //    {
            //        NoOfParts = 5,
            //    },
            //    CoarseSolver = new DirectSolver()
            //    {
            //        WhichSolver = DirectSolver._whichSolver.MUMPS
            //    },
            //    overlap = 1
            //};


            //C.NoOfMultigridLevels = 3;
            //Prec = new Schwarz()
            //{
            //    m_BlockingStrategy = new Schwarz.MultigridBlocks()
            //    {
            //        Depth = 2,
            //    },
            //    CoarseSolver = new DirectSolver()
            //    {
            //        WhichSolver = DirectSolver._whichSolver.MUMPS
            //    },
            //    overlap = 1
            //};

            // For Newton
            C.LinearSolver = Prec;

            ////For Picard
            //C.LinearSolver = new SoftGMRES()
            //{
            //    MaxKrylovDim = C.MaxKrylovDim,
            //    Precond = Prec,
            //    m_Tolerance = 1E-6,
            //    m_MaxIterations = 50
            //};

         

            // Timestepping
            // ============
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1E20;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10000000;
            C.NoOfTimesteps = 1;


            return C;
        }
    }
}
