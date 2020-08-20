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
using System;
using System.Collections.Generic;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedTestExamples {

        public static FSI_Control TestParticleParameter(int k = 2)
        {
            FSI_Control C = new FSI_Control(k, "ParticleParameterTest");

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 8, lengthY: 4, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 1);
            C.pureDryCollisions = true;


            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 0;
            double particleDensity1 = 2;
            InitializeMotion motion1 = new InitializeMotion(C.gravity, particleDensity1, C.pureDryCollisions);
            double particleDensity2 = 1;
            InitializeMotion motion2 = new InitializeMotion(C.gravity, particleDensity2, C.pureDryCollisions);
            C.Particles.Add(new Particle_Sphere(motion1, 1, new double[] { -2.0, 0.0 }));
            C.Particles.Add(new Particle_Ellipsoid(motion2, 1, 1, new double[] { 2.0, 0.0 }, startAngl: 0));

            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = true;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.hydrodynamicsConvergenceCriterion = 1e-2;


            // Timestepping
            // ============

            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e-2;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }


        public static FSI_Control TestParticleInShearFlow(int k = 2) {
            FSI_Control C = new FSI_Control(k, "ParticleInShearFlow");
            // grid and boundary conditions
            // ============================
            List<string> boundaryValues = new List<string> {
                "Velocity_Inlet_left",
                "Velocity_Inlet_right",
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 4, lengthY: 6, cellsPerUnitLength: 5, periodicX: false, periodicY: true);
            C.SetAddaptiveMeshRefinement(amrLevel: 1);

            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityY", X => 0.02);
            C.AddBoundaryValue("Velocity_Inlet_right", "VelocityY", X => -0.02);
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.hydrodynamicsConvergenceCriterion = 1e-1;
            double particleDensity = 1;
            C.gravity = new Vector(0, 0);
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, true);
            C.Particles.Add(new Particle_Sphere(motion, 0.4, new double[] { 0.0, 0.0 }));
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.25;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================
            C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.LevelSetSmoothing = false;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.LinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;

            // Timestepping
            // ============

            C.Timestepper_Scheme = FSI_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 120;
            C.NoOfTimesteps = 25;

           

            // haben fertig...
            // ===============

            return C;
        }

        /// <summary>
        /// Testing of particle/wall interactions using a single particle
        /// </summary>
        public static FSI_Control TestSingleDryParticleAgainstWall(bool meshRefine = false) {
            FSI_Control C = new FSI_Control(degree: 2, projectName: "DryParticleWallCollision");

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 2, lengthY: 2, cellsPerUnitLength: 10, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(meshRefine ? 1 : 0);
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;
            C.pureDryCollisions = true;
            C.PhysicalParameters.IncludeConvection = true;
            // Particle Properties
            // =============================
            double particleDensity = 1;
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions);
            C.fixPosition = true;
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -0.5, -0.5 }, startAngl: 90.0, activeStress: 0, startTransVelocity: new double[] { 1, -1 }));

            double V = 0;
            foreach (var p in C.Particles) {
                V = Math.Max(V, p.Motion.GetTranslationalVelocity(0).L2Norm());
            }
            if (V <= 0)
                throw new ArithmeticException();

            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = (1 / (14 * V)) * (meshRefine ? 0.5 * 0.5 * 0.5 * 0.2 : 0.1);
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 5;
            C.NoOfTimesteps = 500;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LSunderrelax = 1.0;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 100;

            return C;
        }

        public static FSI_Control TestStickyTrap(int k = 2)
        {
            FSI_Control C = new FSI_Control(degree: k, projectName: "ParticleCollisionTest") {
                pureDryCollisions = true
            };
            // grid and boundary conditions
            // ============================ 
            List<string> boundaryValues = new List<string> {
            "Wall_left",
            "Wall_right",
            "Pressure_Outlet_lower",
            "Pressure_Outlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 3, lengthY: 3, cellsPerUnitLength: 10, periodicX: false, periodicY: false);
            //C.SetAddaptiveMeshRefinement(amrLevel: 2);
            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 0;
            C.gravity = new Vector( 0, -9.81 );

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleMass = 1;
            double particleDensity1 = 4.0;
            InitializeMotion motion1 = new InitializeMotion(C.gravity, particleDensity1, C.pureDryCollisions, true);
            double particleDensity2 = 1.0;
            InitializeMotion motion2 = new InitializeMotion(C.gravity, particleDensity2, C.pureDryCollisions, true, true);

            C.Particles.Add(new Particle_Sphere(motion1, 0.18, new double[] { 0.0, 0.6 }));
            C.Particles.Add(new Particle_superEllipsoid(motion2, 0.4, 0.2, 4, new double[] { 0.45, 0 }, startAngl: 45));
            C.Particles.Add(new Particle_superEllipsoid(motion2, 0.4, 0.2, 4, new double[] { -0.45, 0 }, startAngl: -45));

            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.hydrodynamicsConvergenceCriterion = 1e-2;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 35;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control TestHydrodynamicForces(int k = 2)
        {
            FSI_Control C = new FSI_Control(degree: k, projectName: "HydrodynamicForces");
            //C.SetSaveOptions(@"/home/ij83requ/default_bosss_db", 1);

            List<string> boundaryValues = new List<string> {
                "Velocity_Inlet_left",
                "Pressure_Outlet_right",
                "Wall_lower",
                "Wall_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 8, lengthY: 6, cellsPerUnitLength: 2, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 1);
            C.hydrodynamicsConvergenceCriterion = 1e-4;
            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityX", X => 1.0);

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 0.1;
            C.PhysicalParameters.mu_A = 1e-1;
            C.PhysicalParameters.Material = true;
            double particleDensity = 1.0;
            C.gravity = new Vector(0, 0);
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false, 1);
            // Particle Properties
            // =============================   
            C.Particles = new List<Particle>();
            int numOfParticles = 1;
            for (int d = 0; d < numOfParticles; d++)
            {
                C.Particles.Add(new Particle_Sphere(motion, 0.5, new double[] { 0.0, 0.0 }, startAngl: 0));
            }

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = true;

            // misc. solver options
            // =============================  
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
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 1000;



            // Timestepping
            // =============================  
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e-3;
            C.NoOfTimesteps = 2;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control TestPeriodicBoundaries(int k = 2, double aspectRatio = 2) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "2_active_Rods");
            C.SetGrid(lengthX: 2, lengthY: 2, cellsPerUnitLength: 16, periodicX: true, periodicY: true);

            // Fluid Properties
            // =============================
            C.CoefficientOfRestitution = 1;
            C.pureDryCollisions = true;

            // Particle Properties
            // =============================
            double particleDensity = 1;
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, true, false, false, 1);
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.2, 0.2 * aspectRatio, new double[] { 0, 0 }, startAngl: -45, activeStress: 0, startTransVelocity: new double[] { 1, -1 }));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000000;
            C.NoOfTimesteps = 80;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LSunderrelax = 1.0;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 1000000;


            return C;
        }
    }
}
