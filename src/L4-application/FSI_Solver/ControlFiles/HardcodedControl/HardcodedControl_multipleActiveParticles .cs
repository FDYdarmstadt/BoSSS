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
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedControl_multipleActiveParticles : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control ActiveRods_noBackroundFlow(int k = 2) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "41_active_Rods");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 8, lengthY: 8, cellsPerUnitLength: 8, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 2);
            C.hydrodynamicsConvergenceCriterion = 1e-2;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.CoefficientOfRestitution = 1;

            // Particle Properties
            // =============================
            double particleDensity = 2;
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, 1);
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -1, 0 }, startAngl: 12, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 1, 0 }, startAngl: 39, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -2, 2 }, startAngl: -1, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 2.5, 2 }, startAngl: 192, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -2, -1.5 }, startAngl: 20, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 2.5, -1.5 }, startAngl: 42, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0 , 1 }, startAngl: 51, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -3, -3 }, startAngl: 180, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0, -1 }, startAngl: -33, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -0.5, -3 }, startAngl: -48, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 1, -3 }, startAngl: -1, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -2.5, 0 }, startAngl: 15, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 2.5, 0 }, startAngl: 115, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0, 3 }, startAngl: -15, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -2.5, 3 }, startAngl: 5, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 2.5, 3 }, startAngl: 111, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 3, -3 }, startAngl: 187, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -2, 1 }, startAngl: 1, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0.5, -2 }, startAngl: 8, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -3.5, -1.5 }, startAngl: 92, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -3.5, 1.5 }, startAngl: -90, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -1.8, -3 }, startAngl: 90, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 1.8, 2 }, startAngl: -90, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 3, 1 }, startAngl: -180, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -0.7, 2.2 }, startAngl: -180, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -1, 3.1 }, startAngl: -90, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 1.8, 3.1 }, startAngl: -96, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 1.8, -1 }, startAngl: 45, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 3.5, 3.1 }, startAngl: -96, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0.6, 2 }, startAngl: 172, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -0.9, -2.0 }, startAngl: -162, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 1.5, 1 }, startAngl: 181, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -1, 1 }, startAngl: -91, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 3.2, -2 }, startAngl: -179, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -3.5, 3 }, startAngl: -72, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0, 0 }, startAngl: 213, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -1.2, -1 }, startAngl: -43, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -2.4, -0.7 }, startAngl: -3, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -3.2, 0.7 }, startAngl: 43, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 2, -3.4 }, startAngl: 91, activeStress: 10));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1e-3;
            C.dtMax = dt;   
            C.dtMin = dt;
            C.Endtime = 100000000;
            C.NoOfTimesteps = int.MaxValue;
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

        public static FSI_Control FourParticles(int k = 2) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "2_active_Rods");
            C.SetSaveOptions(@"/home/ij83requ/default_bosss_db", 1);
            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet_left",
                "Pressure_Outlet_right",
                "Wall_lower",
                "Wall_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 6, lengthY: 6, cellsPerUnitLength: 4, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 4);
            C.hydrodynamicsConvergenceCriterion = 1e-2;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.01;
            C.PhysicalParameters.IncludeConvection = false;

            // Particle Properties
            // =============================
            double particleDensity = 10000;
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, 1.5);
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.2, 0.1, new double[] { -2, 0 }, startAngl: 6, activeStress: 10000));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.2, 0.1, new double[] { 2, 0 }, startAngl: 178, activeStress: 10000));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.2, 0.1, new double[] { 0, 1.8 }, startAngl: 96, activeStress: 10000));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.2, 0.1, new double[] { 0.2, -1.4 }, startAngl: -56, activeStress: 10000));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000000;
            C.NoOfTimesteps = int.MaxValue;
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

        public static FSI_Control SingleParticles(int k = 2, double distance = -1, double angle = 30) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 10, lengthY: 10, cellsPerUnitLength: 5, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 3);
            C.hydrodynamicsConvergenceCriterion = 1e-3;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.CoefficientOfRestitution = 1;

            // Particle Properties
            // =============================
            double particleDensity = 2;
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, 1);
            C.fixPosition = true;
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.5, 0.25, new double[] { 0.5 + distance / 2, 0 }, startAngl: 180 + angle, activeStress: 1));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1e-1;
            C.SetTimesteps(dt, 50, false);
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
            C.maxIterationsFullyCoupled = 100;


            return C;
        }

        public static FSI_Control TwoParticles(int k = 2, double distance = 3, double angle = 140) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 10, lengthY: 10, cellsPerUnitLength: 5, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 3);
            C.hydrodynamicsConvergenceCriterion = 1e-3;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.CoefficientOfRestitution = 1;

            // Particle Properties
            // =============================
            double particleDensity = 2;
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, 1);
            C.fixPosition = true;
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.5, 0.25, new double[] { -0.5 - distance / 2, 0 }, startAngl: 30, activeStress: 1));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.5, 0.25, new double[] { 0.5 + distance / 2, 0 }, startAngl: 180 + angle, activeStress: 1));
                        
            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1e-1;
            C.SetTimesteps(dt, 50, false);
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
            C.maxIterationsFullyCoupled = 100;
            

            return C;
        }

        public static FSI_Control PackedParticles(int k = 2, double aspectRatio = 0.5) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/home/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);// Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.CoefficientOfRestitution = 1;

            // Particle Properties
            // =============================
            double particleDensity = 1.1;
            double activeStress = 1;
            double particleLength = 0.1;
            double position = 1;
            double distance = position / 2;

            C.SetGrid(lengthX: 2 * position + distance, lengthY: 2 * position + distance, cellsPerUnitLength: 30, periodicX: true, periodicY: true);
            C.SetAddaptiveMeshRefinement(0);
            C.hydrodynamicsConvergenceCriterion = 1e-1;
            
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, 1.5);

            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position, position - 0 * distance }, -181, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position, position - 0.5 * distance }, 55, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 0.5 * distance, position - 2.5 * distance }, -168, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 0.5 * distance, position - 3 * distance }, 134, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 0.5 * distance, position - 3.5 * distance }, -14, activeStress));

            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 0.5 * distance, position - 0 * distance }, 7, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 0.5 * distance, position - 0.5 * distance }, 123, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 1 * distance, position - 2 * distance }, -89, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 1 * distance, position - 3.5 * distance }, 86, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 1 * distance, position - 4 * distance }, -123, activeStress));

            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 2 * distance, position - 0 * distance }, 63, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 2 * distance, position - 1 * distance }, 145, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 1.5 * distance, position - 2.3 * distance }, -91, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 1.5 * distance, position - 3 * distance }, 12, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 1.5 * distance, position - 3.9 * distance }, -65, activeStress));

            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 3.1 * distance, position - 0 * distance }, 321, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 3.0 * distance, position - 1 * distance }, -4, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 2.6 * distance, position - 2.1 * distance }, 78, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 3.2 * distance, position - 2.9 * distance }, -41, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 2.9 * distance, position - 3.6 * distance }, 5, activeStress));

            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 4 * distance, position - 0 * distance }, -156, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 3.6 * distance, position - 1.2 * distance }, 184, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 3.5 * distance, position - 2 * distance }, -274, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 3.9 * distance, position - 3.1 * distance }, 364, activeStress));
            C.Particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { position - 3.5 * distance, position - 4 * distance }, -45, activeStress));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(1e-1, int.MaxValue, false);
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

            string ID = "e302ab4e-04fe-4848-9cb4-8db8128bbbdc";
            C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), -1);

            return C;
        }

        public static FSI_Control FixedParticle(int k = 2, int amrLevel = 7, double aspectRatio = 0.2) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/home/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 200, lengthY: 200, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: amrLevel);
            C.hydrodynamicsConvergenceCriterion = 1e-4;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 10;
            C.PhysicalParameters.IncludeConvection = false;
            C.CoefficientOfRestitution = 1;

            // Particle Properties
            // =============================
            double particleDensity = 1;
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, 1);
            C.Particles.Add(new Particle_Ellipsoid(motion, 1, 1 * aspectRatio, new double[] { 0, 0 }, startAngl: 0, activeStress: 1));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
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