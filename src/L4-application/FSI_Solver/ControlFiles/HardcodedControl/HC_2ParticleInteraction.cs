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

using System.Collections.Generic;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HC_2ParticleInteraction : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control Main(double angle = 45, double distance = 4) {
            FSI_Control C = new FSI_Control(2, "2particleInteractions", "active Particles");
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 20, lengthY: 20, cellsPerUnitLength: 4, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(2);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-6;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 1;

            // Particle Properties
            // =============================   
            C.fixPosition = true;
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, false, false, false, 1.5, true);
            double particleRadius = 0.5;
            double aspectRatio = 2;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, aspectRatio * particleRadius, particleRadius, new double[] { -distance / 2, 0.0 }, angle, 1),
                new Particle_Ellipsoid(motion, aspectRatio * particleRadius, particleRadius, new double[] { distance / 2, 0.0 }, 180 - angle, 1)
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
            C.SetTimesteps(dt: 1e-1, noOfTimesteps: 25);

            return C;
        }

        public static FSI_Control Single(double angle = 45, double distance = 6) {
            FSI_Control C = new FSI_Control(2, "2particleInteractions", "active Particles");
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 40, lengthY: 40, cellsPerUnitLength: 2, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(3);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-6;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 1;

            // Particle Properties
            // =============================   
            C.fixPosition = true;
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, false, false, false, 1.5, true);
            double particleRadius = 0.5;
            double aspectRatio = 2;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, aspectRatio * particleRadius, particleRadius, new double[] { -distance / 2, 0.0 }, angle, 1)
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
            C.SetTimesteps(dt: 1e-1, noOfTimesteps: 25);

            return C;
        }

    }
}