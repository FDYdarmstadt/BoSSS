/* =======================================================================
Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using BoSSS.Solution.XdgTimestepping;
using ilPSP;

namespace BoSSS.Application.FSI_Solver {
    public class HC_CollisionRuntime : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control TwoSpheres(int k = 2) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "TwoSpheres");

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 4, lengthY: 4, cellsPerUnitLength: 20, periodicX: false, periodicY: false);
            
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.pureDryCollisions = true;
            C.gravity = new Vector( 0, -9.81 );

            // Particle Properties
            // =============================
            double particleDensity = 1;
            InitializeMotion motion1 = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false);
            InitializeMotion motion2 = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, true, true);
            C.Particles.Add(new Particle_Ellipsoid(motion1, 0.25, 0.1, new double[] { 0 , 0.5 }, startAngl: 45, activeStress: 0));
            C.Particles.Add(new Particle_Ellipsoid(motion2, 0.25, 0.1, new double[] { 0, 0 }, startAngl: 0, activeStress: 0));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000000;
            C.NoOfTimesteps = 100;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LSunderrelax = 1.0;

            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 1000000;

            return C;
        }
    }
}