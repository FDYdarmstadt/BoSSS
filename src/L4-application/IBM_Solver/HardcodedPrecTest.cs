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
using BoSSS.Foundation.Grid;

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
            C.SetDGdegree(k);

            #region Creates grid () and sets BC
            //// Create Grid
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
            C.NonLinearSolver.MaxSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1E-5;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;                    

            // Solver configuration
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton; // Newton GMRES will be executed if a GMRES linsolver is chosen
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
     

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
