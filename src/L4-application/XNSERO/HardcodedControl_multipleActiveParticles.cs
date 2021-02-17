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
using MPI.Wrappers;

namespace BoSSS.Application.XNSERO_Solver {

    public static class MultiplePacticles {
        public static XNSERO_Control Main(int k = 3, double particleLength = 0.5, double aspectRatio = 0.5, int cellsPerUnitLength = 6, double noOfParticles = 5) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            //C.SetSaveOptions(dataBasePath: @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\packedParticles", savePeriod: 1);
            //string ID = "ccf040a2-e5db-47a7-b226-b642676636cd";
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), -1);
            //C.IsRestart = true;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.IncludeConvection = true;

            // Particle Properties
            // =============================
            double particleDensity = C.PhysicalParameters.rho_A * 1000;
            double activeStress = 1;
            double nextParticleDistance = particleLength * 2;
            double domainLength = nextParticleDistance * noOfParticles;
            C.SetGrid(domainLength, domainLength, cellsPerUnitLength, true, true);
            C.SetAddaptiveMeshRefinement(0);
            C.minDistanceThreshold = 1 / cellsPerUnitLength;
            C.CoefficientOfRestitution = 1;
            InitializeMotion motion = new InitializeMotion(particleDensity, false, false, false, 0);
            double leftCorner = -domainLength / 2 + nextParticleDistance / 2;
            Random angle = new Random();
            int j = 0;
            List<Particle> particles = new List<Particle>();
            while(leftCorner + j * nextParticleDistance < domainLength / 2) {
                int i = 0;
                while(leftCorner + i * nextParticleDistance < domainLength / 2) {
                    double angle2 = (double)angle.Next(0, 6) * 180 + angle.Next(0, 361) * Math.Pow(-1, i * j);
                    angle2 = angle2.MPIBroadcast(0);
                    particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { leftCorner + i * nextParticleDistance, leftCorner + j * nextParticleDistance }, angle2, activeStress, new double[] { 0, 0 }));
                    i += 1;
                }
                j += 1;
            }
            C.SetParticles(particles);
            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.SetTimesteps(1e-1, 50000000);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.0;
            C.LinearSolver.NoOfMultigridLevels = 10;
            C.LinearSolver.MaxSolverIterations = 1;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.verbose = false;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.UseSchurBlockPrec = false;
            C.LinearSolver.pMaxOfCoarseSolver = k;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.TargetBlockSize = 10000;
            return C;
        }

        public static XNSERO_Control TWOPI(double angle = 0, double angleXAxis = 150, double distance = 5, double halfAxis = 0.5, double aspectRatio = 0.5) {
            XNSERO_Control C = new XNSERO_Control(2, "2particleInteractions", "active Particles");

            C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\2PI4", 1);

            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            //C.AlternateDbPaths = new[] { new ValueTuple<string, string>(@"/work/scratch/ij83requ/default_bosss_db", ""), new ValueTuple<string, string>(@"U:\default_bosss_db", "") };
            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Pressure_Dirichlet"
            };
            C.SetBoundaries(boundaryValues);
            double length = 10;
            C.SetGrid(lengthX: length, lengthY: length, cellsPerUnitLength: 14, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(3);

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = true;
            double particleDensity = 10;

            // Particle Properties
            // =============================   
            C.fixPosition = true;
            InitializeMotion motion = new InitializeMotion(particleDensity, false, false, false, 0);
            List<Particle> particles = new List<Particle>();
            particles = new List<Particle> {
                new Particle_Ellipsoid(motion, halfAxis, aspectRatio * halfAxis, new double[] { -distance * Math.Cos(angleXAxis * Math.PI / 180) / 2, -distance * Math.Sin(angleXAxis * Math.PI / 180) / 2 }, 0, 10),
                new Particle_Ellipsoid(motion, halfAxis, aspectRatio * halfAxis, new double[] { distance * Math.Cos(angleXAxis * Math.PI / 180) / 2, distance * Math.Sin(angleXAxis * Math.PI / 180) / 2 }, angle, 10)
            };
            C.SetParticles(particles);
            C.NonLinearSolver.ConvergenceCriterion = 1e-4;

            // Timestepping
            // =============================  
            C.SetTimesteps(dt: 1e-1, noOfTimesteps: 1);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.LinearSolver.TargetBlockSize = 10000;

            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            return C;
        }

    }

}