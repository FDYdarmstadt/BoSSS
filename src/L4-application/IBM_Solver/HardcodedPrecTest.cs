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
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IBM_Solver {
    public class HardcodedPrecTest {

        static public IBM_Control PrecTest3DChannel(int k, int cells_x, int cells_yz) {
            IBM_Control C = new IBM_Control();

            // basic database options
            // ======================
            C.savetodb = false;
            //C.DbPath = @"/home/oe11okuz/BoSSS_DB/Lichtenberg_DB";
            //C.DbPath = @"P:\BoSSS_DBs\Bug";

            //string restartSession = "727da287-1b6a-463e-b7c9-7cc19093b5b3";
            //string restartGrid = "3f8f3445-46f1-47ed-ac0e-8f0260f64d8f";

            //C.DynamicLoadBalancing_Period = 1;
            //C.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication app, int noOfPerformanceClasses) {
            //    Console.WriteLine("i was called");
            //    int[] map = new int[] { 1, 5, 100 };
            //    return new StaticCellCostEstimator(map);
            //});
            C.DynamicLoadBalancing_RedistributeAtStartup = false;

            //c.DynamicLoadBalancing_CellClassifier = new IndifferentCellClassifier();
            C.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 1}));
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 10, 1 }));
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(ArtificialViscosityCellCostEstimator.GetStaticCostBasedEstimator());
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(ArtificialViscosityCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());

            // Assign correct names
            C.SessionName = "Channel_" + k + "_" + cells_x + "x" + cells_yz ;

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
                var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, CellType.Cube_Linear, false, true, false);

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
            C.InitialValues_Evaluators.Add("VelocityX", X => 1 - 4 * (X[2] * X[2]));
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
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", (X, t) => 1 - 4 * (X[2] * X[2]));
            C.AddBoundaryValue("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryValue("Wall");
            C.AddBoundaryValue("Pressure_Outlet");


            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxKrylovDim = 1000;
            C.LinearSolver.MaxSolverIterations = 1;
            C.NonLinearSolver.MaxSolverIterations = 1;
            C.LinearSolver.ConvergenceCriterion = 1E-5;
            C.NonLinearSolver.ConvergenceCriterion = 1E-5;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;                    

            // Solver configuration
            C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.NewtonGMRES;
            C.LinearSolver.SolverCode = LinearSolverConfig.Code.classic_mumps;
     

            // Timestepping
            // ============
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1E20;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10000000;
            C.NoOfTimesteps = 1;
            C.LinearSolver.NoOfMultigridLevels = 7;

            return C;
        }


        static public IBM_Control PrecTest3dDegenhardt(int precNo = 4, int channel = 1, int name_newton =1, int k =3, int cells_x = 4, int cells_yz =5, int re = 100, int ASparts = 3, int ASDepth = 2, int MGLevels = 3, int maxKrDim = 1000, int saveToDB = 1)
        {
            IBM_Control C = new IBM_Control();

            //in SolverFactory die DoF parts ändern

            //Possibilities:
            //channel = 0 --> channel 3D with sphere
            //channel = 1 --> channel 3D empty
            //channel = 2 --> channel 2D with cylinder
            //channel = 3 --> channel 2D empty

            string sessName = "";
            if (channel == 0)
                sessName = "Channel_3D_Sphere";
            else if (channel == 1)
                sessName = "Channel_3D_empty";
            else if (channel == 2)
                sessName = "Channel_2D_Cylinder";
            else if (channel == 3)
                sessName = "Channel_2D_empty";

            string precString = "";
            if (precNo == 0)
                precString = "_noPrec";
            if (precNo == 1)
                precString = "_Schur";
            if (precNo == 2)
                precString = "_Simple";
            if (precNo == 3)
                precString = "_AS-1000";
            if (precNo == 4)
                precString = "_AS-5000";
            if (precNo == 5)
                precString = "_AS-10000";
            if (precNo == 6)
                precString = "_AS-MG";
            if (precNo == 7)
                precString = "_localPrec";



            // basic database options
            // ======================
           // if (saveToDB == 1)
               // C.savetodb = true;
            //else
                //C.savetodb = false;

            C.savetodb = true;

            C.DbPath = @" \\dc1\scratch\Krause\Datenbank_Louis\degenhardt_final";
            //C.DbPath = @"\\hpccluster\hpccluster-scratch\krause\cluster_db";
            //C.DbPath = @"/home/oe11okuz/BoSSS_DB/Lichtenberg_DB";


            //string restartSession = "727da287-1b6a-463e-b7c9-7cc19093b5b3";
            //string restartGrid = "3f8f3445-46f1-47ed-ac0e-8f0260f64d8f";

            C.DynamicLoadBalancing_Period = 1;
            //C.DynamicLoadBalancing_CellCostEstimatorFactory = delegate (IApplication<AppControl> app, int noOfPerformanceClasses) {
            //    Console.WriteLine("i was called");
            //    int[] map = new int[] { 1, 5, 100 };
            //    return new StaticCellCostEstimator(map);
            //};

            

            if (name_newton == 1)
            {
                C.SessionName = "Newton_" + sessName + precString + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;
                C.ProjectDescription = "Newton_" + sessName + precString + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;
            }
            else
            {
                C.SessionName = "Picard_" + sessName + precString + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;
                C.ProjectDescription = "Picard_" + sessName + precString + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim;
            }

            C.saveperiod = 1;
            //C.SessionName = "Sphere_k" + k + "_h" + h+"Re100";
            C.ProjectName = "iteration-study";
            C.Tags.Add("Prec param study");

            // Create Fields
            C.FieldOptions.Add("VelocityX", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts()
            {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            if (channel == 0 || channel == 1) //3D
            {
                C.FieldOptions.Add("VelocityZ", new FieldOpts()
                {
                    Degree = k,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
            }



            #region Creates grid () and sets BC
            //// Create Grid
            Console.WriteLine("...generating grid");
            if (channel == 0 || channel == 1) //3D
            {
                #region grid 3D
                C.GridFunc = delegate
                {
                    // x-direction
                    var _xNodes = GenericBlas.Linspace(-0.5, 1.5, cells_x + 1);

                    // y-direction
                    var _yNodes = GenericBlas.Linspace(-0.5, 0.5, cells_yz + 1);

                    // z-direction
                    var _zNodes = GenericBlas.Linspace(-0.5, 0.5, cells_yz + 1);

                    // Cut Out
                    var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, CellType.Cube_Linear, false, true, false);

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
            }else
            {
                #region grid 2D
                C.GridFunc = delegate {

                    // x-direction
                    var _xnodes = GenericBlas.Linspace(-0.5, 1.5, cells_x + 1);
                    // y-direction
                    var _ynodes = GenericBlas.Linspace(-0.5, 0.5, cells_yz + 1);

                    var grd = Grid2D.Cartesian2DGrid(_xnodes, _ynodes, CellType.Square_Linear, false, false);

                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(2, "Wall");
                    grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                    grd.DefineEdgeTags(delegate (double[] _X)
                    {
                        var X = _X;
                        double x = X[0];
                        double y = X[1];

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

                        throw new ArgumentOutOfRangeException();
                    });

                    return grd;
                                       
                };
                #endregion 
            }
            #endregion


            // set initial conditions
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);

            if (channel == 0 | channel == 1)  //3D
            {
                //C.InitialValues_Evaluators.Add("VelocityX", X => 1 - 4 * (X[2] * X[2]));
                C.InitialValues_Evaluators.Add("VelocityX", X => 0);

                C.InitialValues_Evaluators.Add("VelocityZ", X => 0);
            }
            else
            {
                //C.InitialValues_Evaluators.Add("VelocityX", X => 1 - 4 * (X[1] * X[1]));
                C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            }



            // Because its a sphere

            if (channel == 0)  //3D channel sphere
            {
                C.particleRadius = 0.1;
                C.InitialValues_Evaluators.Add("Phi", x => -(x[0]).Pow2() + -(x[1]).Pow2() + -(x[2]).Pow2() + C.particleRadius.Pow2());
            }
            else if (channel == 1 || channel == 3)  //3D channel empty or 2D channel empty
            {
                C.InitialValues_Evaluators.Add("Phi", x => -1);
            }
            else if (channel == 2)   //2D channel cylinder
            {
                var radius = 0.1;
                C.particleRadius = radius;
                C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + radius.Pow2());
            }




            Console.WriteLine("...starting calculation of Preconditioning test with 3D Channel");
            if (name_newton == 1)
            {
                Console.WriteLine("newton_" + sessName + precString + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim);

            }
            else
            {
                Console.WriteLine("picard_" + sessName + precString + "_k" + k + "_x" + cells_x + "_yz" + cells_yz + "_re" + re + "_asp" + ASparts + "_asd" + ASDepth + "_mgl" + MGLevels + "_kr" + maxKrDim);
            }

           
            // Physical values
            C.PhysicalParameters.rho_A = 1;
            // 1/Re
            //C.PhysicalParameters.mu_A = 1.0 / 10.0;
            //C.PhysicalParameters.mu_A = 0.2 / re;

            C.PhysicalParameters.mu_A = 1.0 / re;

            // Boundary conditions
            C.AddBoundaryValue("Velocity_inlet", "VelocityY", (x, t) => 0);
            C.AddBoundaryValue("Wall");
            C.AddBoundaryValue("Pressure_Outlet");
        
            if (channel == 0 || channel == 1) //3D
            {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", (x, t) => 1 - 4 * (x[2] * x[2]));
            }
            else
            {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", (x, t) => 1 - 4 * (x[1] * x[1]));
            }



            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            //C.LinearSolver.MaxKrylovDim = 1000;
            C.LinearSolver.MaxKrylovDim = maxKrDim;
            C.LinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.ConvergenceCriterion = 1E-5;
            C.NonLinearSolver.ConvergenceCriterion = 1E-5;
            //C.LinearSolver.ConvergenceCriterion = 1E-6;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;

            // Choosing the Preconditioner
            ISolverSmootherTemplate Prec;

            if (name_newton == 1)
                C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.NewtonGMRES;
            else
                C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.PicardGMRES;


            switch (precNo)
            {
                case 0:
                    {
                        Prec = null;
                        break;
                    }
                case 1:
                    {
                        C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_Schur;
                        break;
                    }
                case 2:
                    {

                        C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_Simple;
                        break;
                    }
                case 3:
                    {
                        C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_AS_1000;
                        C.LinearSolver.NoOfMultigridLevels = MGLevels;  // 3 // --> grobes MG am Ende nochmal
                        break;
                    }
                case 4:
                    {
                        C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_AS_5000;
                        C.LinearSolver.NoOfMultigridLevels = MGLevels;
                        break;
                    }
                case 5:
                    {
                        C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_AS_10000;
                        C.LinearSolver.NoOfMultigridLevels = MGLevels;
                        break;
                    }
                case 6:
                    {
                        //depth = 2,
                        //   Depth = ASDepth,  //--> MG bei der Blockzerlegung --> Resultat ergibt die Blöcke zur Berechnung (kleine Blöcke--> schlecht)
                        C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_AS_MG;
                        C.LinearSolver.NoOfMultigridLevels = MGLevels;
                        break;
                    }
                case 7:
                    {
                        C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_localPrec; ;
                        C.LinearSolver.NoOfMultigridLevels = MGLevels;
                        break;
                    }
                case 8:
                    {
                        C.LinearSolver.NoOfMultigridLevels = 5;
                        Prec = new Schwarz()
                        {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy()
                            {
                                //noofparts = 5,
                                NoOfPartsPerProcess = ASparts,
                            },
                            CoarseSolver = new ClassicMultigrid()
                            {
                                CoarserLevelSolver = new ClassicMultigrid()
                                {
                                    CoarserLevelSolver = new ClassicMultigrid()
                                    {
                                        CoarserLevelSolver = new SparseSolver()
                                        {
                                            WhichSolver = SparseSolver._whichSolver.MUMPS
                                        },
                                    },
                                },
                            },
                            Overlap = 1
                        };
                        break;
                    }
                default:
                    {
                        Prec = new SchurPrecond()
                        {
                            SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
                        };
                        break;
                    }
            }


            // For Newton
            //  C.LinearSolver.SolverCoder = Prec;

            ////For Picard
            //C.LinearSolver.SolverCoder = new SoftGMRES()
            //{
            //    MaxKrylovDim = C.LinearSolver.MaxKrylovDim,
            //    Precond_solver = Prec,
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
