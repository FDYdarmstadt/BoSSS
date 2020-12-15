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
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HC_2ParticleInteraction : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control Main(double angle = 180, double verticalDistance = 0, double distance = 1.1, double halfAxis = 0.5, double aspectRatio = 0.5) {
            FSI_Control C = new FSI_Control(2, "2particleInteractions", "active Particles");

            C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\2PI", 1);
            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            //C.AlternateDbPaths = new[] { new ValueTuple<string, string>(@"/work/scratch/ij83requ/default_bosss_db", ""), new ValueTuple<string, string>(@"U:\default_bosss_db", "") };
            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            double length = 10;
            C.SetGrid(lengthX: length, lengthY: length, cellsPerUnitLength: 4, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(0);

            // Coupling Properties
            // =============================
            C.hydrodynamicsConvergenceCriterion = 1e-6;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.IsStationary = false;
            double particleDensity = 10;

            // Particle Properties
            // =============================   
            C.fixPosition = true;
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, false, false, false, 0, false);
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, halfAxis, aspectRatio * halfAxis, new double[] { -distance / 2, -verticalDistance / 2 }, 0, 1),
                new Particle_Ellipsoid(motion, halfAxis, aspectRatio * halfAxis, new double[] { distance / 2, verticalDistance / 2 }, angle, 1)
            };


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-1, noOfTimesteps: 10);
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
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.fullyCoupledSplittingMaxIterations = 100;
            C.FullOutputToConsole = true;

            return C;
        }

        public static FSI_Control Single(double angle = 0, double distance = 0, double aspectRatio = 0.5, double activeStress = 1) {
            FSI_Control C = new FSI_Control(2, "2particleInteractions", "active Particles");
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            //C.AlternateDbPaths = new[] { new ValueTuple<string, string>(@"/work/scratch/ij83requ/default_bosss_db", ""), new ValueTuple<string, string>(@"U:\default_bosss_db", "") };
            //string ID = "b09f7099-871a-49b3-8abb-509c1c2b58f1";
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), -1);
            //C.IsRestart = true;

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Pressure_Dirichlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 100, lengthY: 100, cellsPerUnitLength: 1, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(2);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LevelSetSmoothing = true;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-3;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 10;
            C.PhysicalParameters.IncludeConvection = false;
            C.IsStationary = false;
            double particleDensity = 1;

            // Particle Properties
            // =============================   
            C.fixPosition = true;
            InitializeMotion motion = new InitializeMotion(C.gravity, particleDensity, false, false, false, 0, false);
            double particleRadius = 2.5;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, particleRadius, aspectRatio * particleRadius, new double[] { -distance / 2, -0.0 }, angle, activeStress)
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
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: 100000);

            return C;
        }

    }
}