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
using System.Collections.Generic;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedControl_straightChannel : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control ActiveRod_noBackroundFlow(int k = 2, int angle = 20, double aspectRatio = 3, double activeStress = 10) {
            FSI_Control C = new FSI_Control(k, "activeRod_noBackroundFlow", "active Particles");
            //C.SetSaveOptions(dataBasePath: @"/home/ij83requ/default_bosss_db", savePeriod: 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            string ID = "3debf19b-821a-4a76-8b06-e1b945b4393e";
            C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), -1);
            C.IsRestart = true;
            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Wall_lower",
                "Wall_upper",
                "Pressure_Dirichlet_right",
                "Pressure_Dirichlet_left"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 2.5, lengthY: 0.5, cellsPerUnitLength: 25, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(3);
            C.minDistanceThreshold = 0.005;

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-1;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 500;

            // Particle Properties
            // =============================   
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, false, false, false, 1);
            double particleRadius = 0.05;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, aspectRatio * particleRadius, particleRadius, new double[] { 0, 0.0 }, angle, activeStress)
            };   

            // misc. solver options
            // =============================  
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: int.MaxValue);

            return C;
        }

        public static FSI_Control ActiveRod_noWalls(int k = 3) {
            FSI_Control C = new FSI_Control(k, "activeRod_noBackroundFlow", "active Particles");
            C.SetSaveOptions(dataBasePath: @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\Channel", savePeriod: 1);

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 20, lengthY: 8, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 6);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.maxIterationsFullyCoupled = 100000;
            C.hydrodynamicsConvergenceCriterion = 1e-1;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = true;
            double particleDensity = 1;
            // Particle Properties
            // =============================   
            
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, false, false, false, 1);
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, 0.5, 0.05, new double[] { 0.0, 0.0 }, startAngl: 0, activeStress: 1) 
            };

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            //Initial Values
            // =============================   
            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);

            // For restart
            // =============================  
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("42c82f3c-bdf1-4531-8472-b65feb713326"), 400);
            //C.GridGuid = new Guid("f1659eb6 -b249-47dc-9384-7ee9452d05df");

            // misc. solver options
            // =============================  
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            
            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-3, noOfTimesteps: 1000000000);

            return C;
        }
    }
}