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

using ilPSP;
using System;
using System.Collections.Generic;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.XNSERO_Solver {
    public class XNSEROTest_Control {

        public static XNSERO_Control LevelSetTest(int k = 2) {
            XNSERO_Control C = new XNSERO_Control(k, "ParticleParameterTest");

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 4, lengthY: 4, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 1);

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 0;
            double particleDensity1 = 1;
            InitializeMotion motion1 = new InitializeMotion(particleDensity1);
            List<Particle> particles = new List<Particle> {
                new Particle_Sphere(motion1, 0.5, new double[] { 0.0, 0.0 }),
            };
            C.SetParticles(particles);

            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.2;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;

            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e-2;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }

        public static XNSERO_Control TestParticleParameter(int k = 2) {
            XNSERO_Control C = new XNSERO_Control(k, "ParticleParameterTest");

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 8, lengthY: 4, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 1);

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 0;
            double particleDensity1 = 2;
            InitializeMotion motion1 = new InitializeMotion(particleDensity1);
            double particleDensity2 = 1;
            InitializeMotion motion2 = new InitializeMotion(particleDensity2);
            List<Particle> particles = new List<Particle> {
                new Particle_Sphere(motion1, 1, new double[] { -2.0, 0.0 }),
                new Particle_Ellipsoid(motion2, 1, 1, new double[] { 2.0, 0.0 }, startAngl: 0)
            };
            C.SetParticles(particles);

            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.2;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;

            double dt = 1e-2;
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
            C.SetGrid(lengthX: 4, lengthY: 6, cellsPerUnitLength: 5, periodicX: false, periodicY: true);
            C.SetAddaptiveMeshRefinement(amrLevel: 1);

            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityY#A", X => 0.02);
            C.AddBoundaryValue("Velocity_Inlet_right", "VelocityY#A", X => -0.02);
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            double particleDensity = 1;
            InitializeMotion motion = new InitializeMotion(particleDensity, false, false, true);
            C.SetParticles(new List<Particle> { new Particle_Sphere(motion, 0.4, new double[] { 0.0, 0.0 }) });
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.25;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================
            C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            C.AgglomerationThreshold = 0.1;
            C.LinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MinSolverIterations = 5;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_mumps;
            C.NonLinearSolver.ConvergenceCriterion = 1e-12;

            // Timestepping
            // ============

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

        public static XNSERO_Control TestParticleInShearFlow_Phoretic(int k = 2) {
            //BoSSS.Application.XNSERO_Solver.XNSEROTest_Control.TestParticleInShearFlow_Phoretic();

            var C = TestParticleInShearFlow(k);
            C.UsePhoreticField = true;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_mumps;
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
            C.SetGrid(lengthX: 3, lengthY: 3, cellsPerUnitLength: 10, periodicX: false, periodicY: false);
            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.CoefficientOfRestitution = 0;
            C.SetGravity(new Vector(0, -9.81));

            // Particle Properties
            double particleDensity1 = 100.0;
            InitializeMotion motion1 = new InitializeMotion(particleDensity1, false, true);
            double particleDensity2 = 1.0;
            InitializeMotion motion2 = new InitializeMotion(particleDensity2, false, true, true);
            List<Particle> particles = new List<Particle> {
                new Particle_Sphere(motion1, 0.18, new double[] { 0.0, 0.6 }),
                new Particle_superEllipsoid(motion2, 0.4, 0.2, 4, new double[] { 0.45, 0 }, startAngl: 45),
                new Particle_superEllipsoid(motion2, 0.4, 0.2, 4, new double[] { -0.45, 0 }, startAngl: -45),
            };
            C.SetParticles(particles);
            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.LinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.LinearSolver.NoOfMultigridLevels = 1;

            // ============

            // Timestepping
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 25;

            // haben fertig...
            // ===============

            return C;
        }
    }
}
