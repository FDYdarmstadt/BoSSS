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
    public class HC_Dry : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control SingleParticleFalling(int k = 2, int amrLevel = 1) {
            FSI_Control C = new FSI_Control(k, "activeRod_noBackroundFlow", "active Particles");
            //C.SetSaveOptions(dataBasePath: @"/home/ij83requ/default_bosss_db", savePeriod: 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            //C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 1, lengthY: 1, cellsPerUnitLength: 4, periodicX: true, periodicY: true);
            C.SetAddaptiveMeshRefinement(amrLevel);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-4;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 100;
            C.pureDryCollisions = true;

            // Particle Properties
            // =============================   
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, true, false, false);
            double particleRadius = 0.1;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, particleRadius, particleRadius, new double[] {0,0 },0, 0, new double[] {0,-0.1 })
            };   

            // misc. solver options
            // =============================  
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();
            C.LevelSetSmoothing = false;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: int.MaxValue);

            return C;
        }
        public static FSI_Control TwoParticleCollision(int k = 2, int amrLevel = 2) {
            FSI_Control C = new FSI_Control(k, "activeRod_noBackroundFlow", "active Particles");
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 3, lengthY: 3, cellsPerUnitLength: 5, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-4;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 100;
            C.pureDryCollisions = true;

            // Particle Properties
            // =============================   
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, true, false, false);
            double particleRadius = 0.1;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, particleRadius, 0.25 * particleRadius, new double[] { -0.1, 0.0 },0, 0, new double[] {0,-0.05 }, -0.05),
                new Particle_Ellipsoid(motion, particleRadius, 0.25 * particleRadius, new double[] {-0,0.125 },0, 0, new double[] {0,-0.1 }, 1)
            };

            // misc. solver options
            // =============================  
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();
            C.LevelSetSmoothing = false;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: 300);

            return C;
        }

        public static FSI_Control ThreeParticleCollision(int k = 2, int amrLevel = 0) {
            FSI_Control C = new FSI_Control(k, "activeRod_noBackroundFlow", "active Particles");
            //C.SetSaveOptions(dataBasePath: @"/home/ij83requ/default_bosss_db", savePeriod: 1);
            C.SetSaveOptions(@"D:\BoSSS_databases\wetParticleCollision", 1);
            // Domain
            // =============================
            C.SetGrid(lengthX: 1, lengthY: 1, cellsPerUnitLength: 40, periodicX: true, periodicY: true);
            C.SetAddaptiveMeshRefinement(amrLevel);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-4;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1e-3;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 1;
            C.pureDryCollisions = true;

            // Particle Properties
            // =============================   
            InitializeMotion motion1 = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false);
            double particleRadius = 0.1;
            C.Particles = new List<Particle>();
            //for (int i = 0; i < 9; i++) {
            //    for (int j = 0; j < 9; j++) {
            //        C.Particles.Add(new Particle_Sphere(motion1, particleRadius, new double[] { -2 + i * 0.5+0.1*Math.Pow(-1,j), 2 - j * 0.5 }, 0, 0, new double[] { 0.1 * Math.Pow(-1, i), -0.1 * Math.Pow(-1, j) }));
            //    }
            //}
            C.Particles.Add(new Particle_Sphere(motion1, particleRadius, new double[] { 0, 0 }, 0, 0));
            C.Particles.Add(new Particle_Sphere(motion1, particleRadius, new double[] { -0.36, 0 }, 0, 0, new double[] { 1, 0 }));
            //C.Particles.Add(new Particle_Sphere(motion1, particleRadius, new double[] { 0.36, 0 }, 180, 1));

            // misc. solver options
            // =============================  
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();
            C.LevelSetSmoothing = false;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: 5000);

            return C;
        }
    }
}