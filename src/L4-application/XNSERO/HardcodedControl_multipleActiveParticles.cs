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
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using MPI.Wrappers;

namespace BoSSS.Application.XNSERO_Solver {

    public static class MultiplePacticles {
        public static XNSERO_Control Main(int k = 2, double particleLength = 0.5, double aspectRatio = 0.5, int cellsPerUnitLength = 10, double noOfParticles = 7) {
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
            C.PhysicalParameters.IncludeConvection = false;

            // Particle Properties
            // =============================
            double particleDensity = C.PhysicalParameters.rho_A * 1000;
            double activeStress = 1;
            double nextParticleDistance = particleLength * 2.25;
            double domainLength = nextParticleDistance * noOfParticles;
            C.SetGrid(domainLength, domainLength, cellsPerUnitLength, true, true);
            C.minDistanceThreshold = 2 / cellsPerUnitLength;
            C.CoefficientOfRestitution = 1;
            InitializeMotion motion = new InitializeMotion(particleDensity, false, false, false, 0);
            double leftCorner = -domainLength / 2 + nextParticleDistance / 2;
            Random angle = new Random();
            int j = 0;
            List<Particle> particles = new List<Particle>();
            while (leftCorner + j * nextParticleDistance < domainLength / 2) {
                int i = 0;
                while (leftCorner + i * nextParticleDistance < domainLength / 2) {
                    double angle2 = (double)angle.Next(0, 6) * 180 + angle.Next(0, 361) * Math.Pow(-1, i * j);
                    angle2 = angle2.MPIBroadcast(0);
                    particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { leftCorner + i * nextParticleDistance, leftCorner + j * nextParticleDistance }, angle2, activeStress, new double[] { 0, 0 }));

                    i += 1;
                }
                j += 1;
            }
            C.SetParticles(particles);
            C.SetTimesteps(dt: 1e-1, noOfTimesteps: int.MaxValue);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.LinearSolver.TargetBlockSize = 10000;
            C.LinearSolver.verbose = false;
            C.UseSchurBlockPrec = false;
            C.LinearSolver.pMaxOfCoarseSolver = 1;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.ConvergenceCriterion = 1e-4;

            return C;
        }

        public static XNSERO_Control TWOPI(double angle = 180, double angleXAxis = 0, double distance = 2, double halfAxis = 0.5, double aspectRatio = 0.5, double deviation = 0) {
            XNSERO_Control C = new XNSERO_Control(2, "2particleInteractions", "active Particles");

            C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\2PI6", 1);

            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/2PI", 1);
            //C.AlternateDbPaths = new[] { new ValueTuple<string, string>(@"/work/scratch/ij83requ/default_bosss_db", ""), new ValueTuple<string, string>(@"U:\default_bosss_db", "") };
            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 20, lengthY: 20, cellsPerUnitLength: 2, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(3);
            //List<string> boundaryValues = new List<string> {
            //    "Pressure_Outlet"
            //};
            //C.SetBoundaries(boundaryValues);
            //double length = 6;
            //double cellsPerUnitLength = 8;
            //C.SetGrid(lengthX: length, lengthY: length, cellsPerUnitLength, periodicX: false, periodicY: false);
            C.AdaptiveMeshRefinement = true;
            C.AMR_startUpSweeps = 2;
            C.LS_TrackerWidth = 2;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            
            double particleDensity = 100;
            // Particle Properties
            // =============================   
            C.fixPosition = true;
            InitializeMotion motion = new InitializeMotion(particleDensity, false, false, false, 1.5);
            List<Particle> particles = new List<Particle>();
            angleXAxis = angleXAxis * Math.PI / 180;
            particles = new List<Particle> {
                new Particle_Ellipsoid(motion, halfAxis, aspectRatio * halfAxis, new double[] { -distance * Math.Cos(angleXAxis) / 2 + deviation, -distance * Math.Sin(angleXAxis) / 2+deviation }, 0, 1),
                new Particle_Ellipsoid(motion, halfAxis, aspectRatio * halfAxis, new double[] { distance * Math.Cos(angleXAxis) / 2 + deviation, distance * Math.Sin(angleXAxis) / 2+deviation }, angle, 1)
            };
            C.SetParticles(particles);
            C.NonLinearSolver.ConvergenceCriterion = 1e-12;
            List<Func<Vector,double, bool>> Contains = new List<Func<Vector, double, bool>>();
            for(int p = 0; p < particles.Count; p++) {
                Contains.Add(particles[p].Contains);
            }
            C.activeAMRlevelIndicators.Add(new AMRForRigidObject(Contains, 2) { maxRefinementLevel = 3 });
            // Timestepping
            // =============================  
            C.SetTimesteps(dt: 1e-1, noOfTimesteps: 10);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LinearSolver.NoOfMultigridLevels = 10;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.LinearSolver.TargetBlockSize = 10000;
            C.LinearSolver.verbose = false;
            C.UseSchurBlockPrec = false;
            C.LinearSolver.pMaxOfCoarseSolver = 1;

            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            return C;
        }

    }

}