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
        public static FSI_Control ActiveRod_noBackroundFlow(int k = 3, int amrLevel = 2, double aspectRatio = 3, double relaxationFactor = 0.1, bool addaptiveUnderrelaxation = false) {
            FSI_Control C = new FSI_Control(k, "activeRod_noBackroundFlow", "active Particles");
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet_left",
                "Pressure_Outlet_right",
                "Pressure_Outlet_lower",
                "Pressure_Outlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 10, lengthY: 4, cellsPerUnitLength: 3, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-2;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 2;
            C.gravity = new double[] { 0, 0 };

            // Particle Properties
            // =============================   
            C.underrelaxationParam = new ParticleUnderrelaxationParam(convergenceLimit: C.hydrodynamicsConvergenceCriterion, relaxationFactor, addaptiveUnderrelaxation);
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, C.underrelaxationParam, 1.5);
            double particleRadius = 0.25;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, aspectRatio * particleRadius, particleRadius, new double[] { 0.0, 0.0 }, startAngl: 0, activeStress: 1)
            };   

            // misc. solver options
            // =============================  
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            

            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-3, noOfTimesteps: 100000);

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
            C.gravity = new double[] { 0, 0 };
            double particleDensity = 1;
            // Particle Properties
            // =============================   
            C.underrelaxationParam = new ParticleUnderrelaxationParam(convergenceLimit: C.hydrodynamicsConvergenceCriterion, relaxationFactor: 3.0, useAddaptiveUnderrelaxation: true);
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, C.underrelaxationParam, 1);
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