﻿/* =======================================================================
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

using System.Collections.Generic;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using BoSSS.Solution.Control;
using System;

namespace BoSSS.Application.FSI_Solver {
    public class ParticleStokesFlow : IBM_Solver.HardcodedTestExamples {

        public static FSI_Control StokesFlow(int k = 2, int amrLevel = 1) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\2particleInteractions", 1);

            List<string> boundaryValues = new List<string> {
                "Wall",
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 5, lengthY: 5, cellsPerUnitLength: 8, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel);
            C.hydrodynamicsConvergenceCriterion = 1e-2;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector(0, -0.0 );
            // Particle Properties
            // =============================   
            double particleDensity = 10;
            C.Particles = new List<Particle>();
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, false, false, false );
            C.Particles.Add(new Particle_superEllipsoid(motion, 0.4,0.2,4,new double[] { 0.0, 0.0 },45, 0, new double[] { 0, 0 }, 1));
            //C.AddBoundaryValue("Velocity_Inlet_left", "VelocityY", X => 0.1);
            //C.AddBoundaryValue("Velocity_Inlet_right", "VelocityY", X => -0.1);
            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = false;
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-1, noOfTimesteps: 100000);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = true;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.LinearSolver.TargetBlockSize = 10000;

            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSILieSplittingFullyCoupled;
            C.fullyCoupledSplittingMaxIterations = 100;
            C.FullOutputToConsole = true;

            return C;
        }

        public static FSI_Control WetParticleWallCollision(double DensityFactor = 400) {
            FSI_Control C = new FSI_Control(degree: 3, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\wetParticleCollision", 1);
            //C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\WetParticleCollision", 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);

            List<string> boundaryValues = new List<string> {
                "Wall_left",
                "Wall_right",
                "Wall_lower",
                "Pressure_Outlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetAddaptiveMeshRefinement(4);
            C.SetGrid(lengthX: 5, lengthY: 1, cellsPerUnitLength: 6, periodicX: false, periodicY: false);
            C.hydrodynamicsConvergenceCriterion = 1e-3;
            C.pureDryCollisions = false;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector( 0, -10 );
            double particleDensity = 1 * DensityFactor;
            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false, 0);
            C.Particles.Add(new Particle_Sphere(motion, 0.125, new double[] { 0.0, -0.0002 }, 0, 0, new double[] { 0, 0 }));

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = false;

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
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSILieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.fullyCoupledSplittingMaxIterations = 2000;

            
            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(1e-3, 500, false);

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control WetParticleWallCollision2(double DensityFactor = 220, double kRes = 0.8, double minDis = 3e-3) {
            FSI_Control C = new FSI_Control(degree: 3, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\wetParticleCollisionWOGRavitiy", 1);
            //C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\WetParticleCollision", 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);

            List<string> boundaryValues = new List<string> {
                "Wall_left",
                "Wall_right",
                "Wall_lower",
                "Pressure_Outlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetAddaptiveMeshRefinement(4);
            C.SetGrid(lengthX: 5, lengthY: 1, cellsPerUnitLength: 6, periodicX: false, periodicY: false);
            C.hydrodynamicsConvergenceCriterion = 1e-3;
            C.pureDryCollisions = false;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector(0, 0);
            double particleDensity = 1 * DensityFactor;
            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            C.minDistanceThreshold = minDis;
            C.CoefficientOfRestitution = kRes;
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false, 0);
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.0, -0.3 }, 0, 0, new double[] { 0, -2 }));

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = false;

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
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSILieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.fullyCoupledSplittingMaxIterations = 2000;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(1e-5, 5200, true);

            // haben fertig...
            // ===============

            return C;
        }
        public static FSI_Control WetParticlParticleCollision(double DensityFactor = 50) {
            FSI_Control C = new FSI_Control(degree: 3, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\wetParticleParticleCollision", 1);
            //C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\WetParticleCollision", 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetAddaptiveMeshRefinement(2);
            C.SetGrid(lengthX: 2, lengthY: 2, cellsPerUnitLength: 4, periodicX: false, periodicY: false);
            C.hydrodynamicsConvergenceCriterion = 1e-3;
            C.pureDryCollisions = true;
            C.FullOutputToConsole = true;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector(0, 0);
            double particleDensity = 1 * DensityFactor;
            C.minDistanceThreshold = 1e-2;
            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false, 0);
            C.Particles.Add(new Particle_Ellipsoid(motion,0.2, 0.1, new double[] { -0.3, 0.0 }, 0, 0, new double[] { 1, 0 }));
            C.Particles.Add(new Particle_Ellipsoid(motion,0.2, 0.1, new double[] { 0.3, 0.0 }, 90, 0, new double[] { -1, 0 }, 1));

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = false;

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
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSILieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.fullyCoupledSplittingMaxIterations = 2000;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(1e-2, 5000, true);

            

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control TripleCollision(double DensityFactor = 500) {
            FSI_Control C = new FSI_Control(degree: 3, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\wetParticleParticleCollision", 1);
            //C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\WetParticleCollision", 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);

            List<string> boundaryValues = new List<string> {
                "Wall_left",
                "Wall_right",
                "Wall_lower",
                "Pressure_Outlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetAddaptiveMeshRefinement(3);
            C.SetGrid(lengthX: 2, lengthY: 2, cellsPerUnitLength: 5, periodicX: false, periodicY: false);
            C.hydrodynamicsConvergenceCriterion = 1e-4;
            C.pureDryCollisions = false;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector(0, 0);
            double particleDensity = 1 * DensityFactor;
            C.minDistanceThreshold = 3e-3;
            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false, 0);
            C.Particles.Add(new Particle_Sphere(motion, 0.075, new double[] { -0.3, 0.0 }, 0, 0, new double[] { 1, 0 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.075, new double[] { 0.3, 0.0 }, 180, 0, new double[] { -1, 0 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.075, new double[] { 0.0001, 0.3 }, 180, 0, new double[] { 0, -1 }));

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = false;

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
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSILieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.fullyCoupledSplittingMaxIterations = 2000;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(1e-3, 200, true);



            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control AgglomerationCollision(double DensityFactor = 5000) {
            FSI_Control C = new FSI_Control(degree: 3, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\agglomerationCollision", 10);
            //C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\WetParticleCollision", 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);

            //List<string> boundaryValues = new List<string> {
            //    "Wall_lower",
            //    "Pressure_Outlet_upper"
            //};
            //C.SetBoundaries(boundaryValues);

            C.SetGrid(lengthX: 3, lengthY: 3, cellsPerUnitLength: 30, periodicX: true, periodicY: true);
            C.hydrodynamicsConvergenceCriterion = 1e-3;
            C.pureDryCollisions = false;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector(0, 0);
            double particleDensity = 1 * DensityFactor;
            C.minDistanceThreshold = 1/20;
            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, C.pureDryCollisions, false, false, 0);
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 1, 1 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 1, 0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 1, 0 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 1, -0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 1, -1 }, 0, 0, new double[] { 0, -1 }));

            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.5, 1 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.5, 0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.5, 0 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.5, -0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.5, -1 }, 0, 0, new double[] { 0, -1 }));

            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.0, 1 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.0, 0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.0, 0.0 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.0, -0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { 0.0, -1 }, 0, 0, new double[] { 0, -1 }));

            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -0.5, -1 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -0.5, -0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -0.5, 0.0 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -0.5, 0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -0.5, 1 }, 0, 0, new double[] { 0, -1 }));

            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -1, -1 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -1, -0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -1, 0.0 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -1, 0.5 }, 0, 0, new double[] { 0, -1 }));
            C.Particles.Add(new Particle_Sphere(motion, 0.1, new double[] { -1, 1 }, 0, 0, new double[] { 0, -1 }));

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = false;

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
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSILieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.fullyCoupledSplittingMaxIterations = 2000;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(1e-4, int.MaxValue, true);

            // haben fertig...
            // ===============

            return C;
        }
    }
}
