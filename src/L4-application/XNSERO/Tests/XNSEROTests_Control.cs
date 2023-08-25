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

using ilPSP;
using System.Collections.Generic;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.XNSERO_Solver {
    public class XNSEROTest_Control {

        public static XNSERO_Control LevelSetTest(int k = 2) {
            XNSERO_Control C = new XNSERO_Control(k, "ParticleParameterTest");

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(lengthX: 4, lengthY: 4, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(MaxRefinementLevel: 0);

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            double dt = 1e-2;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 0;
            double particleDensity1 = 1;
            Motion motion1 = new(particleDensity1);
            List<Particle> particles = new List<Particle> {
                new ParticleDisk(motion1, 0.5, new double[] { 0.0, 0.0 }),
            };
            C.InitialiseParticles(particles);

            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.2;
            C.NonLinearSolver.MaxSolverIterations = 10;

            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e-2;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }

        public static XNSERO_Control TestParticleProperties(Particle TestParticle) {
            XNSERO_Control C = new XNSERO_Control(2, "ParticleParameterTest");

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(lengthX: 10, lengthY: 10, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(MaxRefinementLevel: 1);

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 1.0;
            C.CoefficientOfRestitution = 0;
            double dt = 1e-2;
            C.InitialiseParticles(new List<Particle>() { TestParticle });

            C.PhysicalParameters.IncludeConvection = false;

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.2;
            C.NonLinearSolver.MaxSolverIterations = 10;

            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e-2;
            C.NoOfTimesteps = 1;

            return C;
        }

        public static XNSERO_Control TestParticleParameter(int k = 2) {
            XNSERO_Control C = new XNSERO_Control(k, "ParticleParameterTest");

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(lengthX: 8, lengthY: 4, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(MaxRefinementLevel: 1);

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 0;
            double particleDensity1 = 2; 
            Motion motion1 = new(particleDensity1);
            double particleDensity2 = 1;
            Motion motion2 = new(particleDensity2);
            List<Particle> particles = new List<Particle> {
                new ParticleDisk(motion1, 1, new double[] { -2.0, 0.0 }),
                new ParticleEllipse(motion2, 1, 1, new double[] { 2.0, 0.0 }, startAngl: 0)
            };
            double dt = 1e-2;
            C.InitialiseParticles(particles);

            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.2;
            C.NonLinearSolver.MaxSolverIterations = 10;

            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e-2;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }

        public static XNSERO_Control TestParticleInShearFlow(int k = 2) {
            XNSERO_Control C = new XNSERO_Control(k, "ParticleInShearFlow");
            // grid and boundary conditions
            // ============================
            List<string> boundaryValues = new List<string> {
                "Velocity_Inlet_left",
                "Velocity_Inlet_right",
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(lengthX: 4, lengthY: 6, cellsPerUnitLength: 5, periodicX: false, periodicY: true);
            C.SetAddaptiveMeshRefinement(MaxRefinementLevel: 1);

            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityY#A", X => 0.02);
            C.AddBoundaryValue("Velocity_Inlet_right", "VelocityY#A", X => -0.02);
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            double dt = 0.1;

            double particleDensity = 1;
            MotionWetNoTranslation motion = new(particleDensity);
            C.InitialiseParticles(new List<Particle> { new ParticleDisk(motion, 0.4, new double[] { 0.0, 0.0 }) });
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.25;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================
            C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            C.AgglomerationThreshold = 0.1;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MinSolverIterations = 5;
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            // Timestepping
            // ============

            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 120;
            C.NoOfTimesteps = 10;


            // haben fertig...
            // ===============

            return C;
        }

        public static XNSERO_Control TestParticleInShearFlow_Phoretic(int k = 2) {
            //BoSSS.Application.XNSERO_Solver.XNSEROTest_Control.TestParticleInShearFlow_Phoretic();

            var C = TestParticleInShearFlow(k);
            C.UsePhoreticField = true;
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
            C.NonLinearSolver.ConvergenceCriterion = 1e-12;
            return C;
        }

        public static XNSERO_Control TestStickyTrap(int k = 2) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "ParticleCollisionTest");
            
            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            
            // grid and boundary conditions
            // ============================ 
            List<string> boundaryValues = new List<string> {
            "Wall_left",
            "Wall_right",
            "Pressure_Dirichlet_lower",
            "Pressure_Dirichlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(lengthX: 3, lengthY: 3, cellsPerUnitLength: 20, periodicX: false, periodicY: false);
            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.CoefficientOfRestitution = 0;
            C.SetGravity(new Vector(0, -9.81));
            C.SkipSolveAndEvaluateResidual = true;

            // Particle Properties
            double particleDensity1 = 100.0;
            double particleDensity2 = 1.0;
            List<Particle> particles = new List<Particle> {
                new ParticleDisk(new Motion(particleDensity1), 0.18, new double[] { 0.0, 0.6 }),
                new ParticleSuperEllipsoidFlat(new MotionFixed(particleDensity2), 0.4, 0.2, 4, new double[] { 0.45, 0 }, startAngl: 45),
                new ParticleSuperEllipsoidFlat(new MotionFixed(particleDensity2), 0.4, 0.2, 4, new double[] { -0.45, 0 }, startAngl: -45),
            };
            double dt = 1e-1;
            C.InitialiseParticles(particles);
            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.1;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MinSolverIterations = 3;

            // ============

            // Timestepping
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 300;

            // haben fertig...
            // ===============

            return C;
        }
    }
}
