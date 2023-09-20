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
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using MPI.Wrappers;

namespace BoSSS.Application.XNSERO_Solver {

    public static class MultiplePacticles {
        public static XNSERO_Control Main(int k = 2, double particleLength = 0.5, double aspectRatio = 0.3, int cellsPerUnitLength =25, double noOfParticles = 3) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            //C.SetSaveOptions(dataBasePath: @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\packedParticles", savePeriod: 1);
            string ID = "e5a0398f-1cc2-45cf-85ef-ab16c76b538d";
            Guid SID = new Guid(ID);
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));
            int timestep = 699;
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), timestep);
            //C.IsRestart = true;
            string PathToOldSessionDir = @"D:\BoSSS_databases\Channel\sessions\e5a0398f-1cc2-45cf-85ef-ab16c76b538d";
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.LS_TrackerWidth = 4;

            // Particle Properties
            // =============================
            double particleDensity = C.PhysicalParameters.rho_A * 1000;
            double activeStress = 1e1;
            double nextParticleDistance = particleLength * 3;
            double domainLength = nextParticleDistance * noOfParticles;
            List<string> boundaryValues = new List<string> {
                "Wall_upper","Wall_lower"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(5*domainLength, domainLength, cellsPerUnitLength, true, false);
            C.CoefficientOfRestitution = 0.5;
            Motion motion = new(particleDensity);
            Motion noMotion = new(particleDensity * 1000);
            double leftCorner = -domainLength / 2 + nextParticleDistance / 2;
            Random angle = new Random();
            int j = 0;
            List<Particle> particles = new List<Particle>();
            while (-2.5 + j * nextParticleDistance * aspectRatio * 1.6 < domainLength / 2) {
                int i = 0;
                while (-domainLength * 5 / 2 + nextParticleDistance / 2 + i * nextParticleDistance < domainLength*5/2){//-domainLength /3) {
                    double angle2 = (double)angle.Next(0, 6) * 180 + angle.Next(0, 361) * Math.Pow(-1, i * j);
                    angle2 = 0;// angle2.MPIBroadcast(0);
                    particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { -domainLength * 5 / 2 + nextParticleDistance / 2 + i * nextParticleDistance, -2.5 + 1.6 * j * nextParticleDistance * aspectRatio }, angle2, activeStress, new double[] { 1, 0 }));

                    i += 1;
                }
                j += 1;
            }
            //particles.Add(new Particle_Ellipsoid(noMotion, domainLength / 5, domainLength / 5, new double[] { 0, 0 }, 0, 0, new double[] { 0, 0 }));
            //j = 0;
            //while (-4.5 + j * nextParticleDistance * aspectRatio * 1.6 < domainLength / 2) {
            //    int i = 0;
            //    while (domainLength / 3 + nextParticleDistance / 2 + i * nextParticleDistance < domainLength * 5 / 2) {
            //        double angle2 = (double)angle.Next(0, 6) * 180 + angle.Next(0, 361) * Math.Pow(-1, i * j);
            //        angle2 = 0;// angle2.MPIBroadcast(0);
            //        particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { domainLength / 3 + nextParticleDistance / 2 + i * nextParticleDistance, -4.5 + 1.6 * j * nextParticleDistance * aspectRatio }, angle2, activeStress, new double[] { 1, 0 }));

            //        i += 1;
            //    }
            //    j += 1;
            //}
            C.AddBoundaryValue("Wall_upper", "VelocityX#A", (X, t) => 1);
            C.AddBoundaryValue("Wall_lower", "VelocityX#A", (X, t) => 0);
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: int.MaxValue);
            C.InitialiseParticles(particles);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.UseSchurBlockPrec = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.ConvergenceCriterion = 1e-4;

            return C;
        }

        public static XNSERO_Control Closed(int k = 2, double particleLength = 0.5, double aspectRatio = 0.4, int cellsPerUnitLength = 10, double lengthParameter = 4) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            //C.SetSaveOptions(dataBasePath: @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\packedParticles", savePeriod: 1);
            string ID = "386a16bc-b0f4-4031-a300-f00f310e98c3";
            Guid SID = new Guid(ID);
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));
            int timestep = 5553;
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), timestep);
            string PathToOldSessionDir = @"D:\BoSSS_databases\Channel\sessions\386a16bc-b0f4-4031-a300-f00f310e98c3";
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.LS_TrackerWidth = 2;
            C.TracingNamespaces = "BoSSS.Foundation";

            // Particle Properties
            // =============================
            double particleDensity = C.PhysicalParameters.rho_A * 1000;
            double activeStress = 1e1;
            double distanceParam = 2.45;
            double nextParticleDistance = particleLength * distanceParam;
            double domainLength = nextParticleDistance * lengthParameter;
            List<string> boundaryValues = new() {  "wall" };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(domainLength, domainLength, cellsPerUnitLength, false, false);
            C.CoefficientOfRestitution = 1;
            double leftCorner = -domainLength / 2 + nextParticleDistance / 2 ;
            int j = 0;
            List<Particle> particles = new List<Particle>();
            double densityParam = 6;
            while (leftCorner + j * nextParticleDistance * aspectRatio * densityParam < domainLength / 2) {//nextParticleDistance * aspectRatio * 1.3
                int i = 0;
                while (leftCorner + i * nextParticleDistance < domainLength / 2) {
                    double angle2 = 180 * (-i + j);
                    Motion motion = new(particleDensity);
                    particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { leftCorner + i * nextParticleDistance, leftCorner + j * nextParticleDistance * aspectRatio * densityParam }, angle2, activeStress));

                    i += 1;
                }
                j += 1;
            }
            //particles.Clear();
            //particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { 0, 0 }, 0, -activeStress));
            //for (int i = 0; i < 211; i++) {
            //   particles.Add(new Particle_Ellipsoid(motion, particleLength, particleLength * aspectRatio, new double[] { 0, 0 }, 0, activeStress));
            //}
            //particles.Add(new Particle_Ellipsoid(motion, 1, 1, new double[] { 0, 0 }, 0, 0));
            //particles = C.LoadParticlesOnRestart(PathToOldSessionDir, particles, timestep);
            //Console.WriteLine("particle length " + particles.Count);
            //for (int i = particles.Count() - 1; i >= 0; i--) {
            //    if (particles[i].Motion.GetPosition().Abs() < 1.5)
            //        particles.RemoveAt(i);
            //}
            C.InitialiseParticles(particles);
            C.WallPositionPerDimension[0][0] = -domainLength / 2;
            C.WallPositionPerDimension[0][1] = domainLength / 2;
            C.WallPositionPerDimension[1][0] = -domainLength / 2;
            C.WallPositionPerDimension[1][1] = domainLength / 2;
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: int.MaxValue);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.UseSchurBlockPrec = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.ConvergenceCriterion = 1e-4;
            C.NonLinearSolver.MinSolverIterations = 2;

            return C;
        }

        public static XNSERO_Control SaffmanForce(int k = 4, double particleLength = 0.5, double aspectRatio = 0.7, double cellsPerUnitLength = 0.5) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "2_active_Rods", IsRestart: false, PathToOldSessionDir: @"D:\BoSSS_databases\Channel\sessions\a5553a60-868e-44bf-834d-ede0561fd7c6", TimestepNoForRestart: 4960);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            //C.SetSaveOptions(dataBasePath: @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\packedParticles", savePeriod: 1);
            string ID = "a5553a60-868e-44bf-834d-ede0561fd7c6";
            Guid SID = new Guid(ID);
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));
            int timestep = 4960;
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), timestep);
            string PathToOldSessionDir = @"D:\BoSSS_databases\Channel\sessions\a5553a60-868e-44bf-834d-ede0561fd7c6";
           
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 10;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.LS_TrackerWidth = 2;
            //C.TracingNamespaces = "BoSSS.Foundation";

            // Particle Properties
            // =============================
            double particleDensity = C.PhysicalParameters.rho_A * 10;
            double activeStress = 10;
            List<string> boundaryValues = new() {
                "Wall_upper",
                "Wall_lower",
                "Pressure_Outlet_right",
                "Pressure_Outlet_left"
            };
            C.AddBoundaryValue("Wall_lower", "VelocityX#A", X => -0.1);
            C.AddBoundaryValue("Wall_upper", "VelocityX#A", X => 0.1);
            C.AddInitialValue("VelocityX", new Formula(@"X=>0.1*X[1]/6"));
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(12, 12, cellsPerUnitLength, false, false);
            C.SetAddaptiveMeshRefinement(4);
            C.CoefficientOfRestitution = 1;
            MotionFixedWithVelocity motion = new(particleDensity);
            List<Particle> particles = new List<Particle> {
                new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { 0, 0 }, 0, activeStress)
            };
            C.InitialiseParticles(particles);
            C.SetTimesteps(dt: 1e-4, noOfTimesteps: 10000);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.UseSchurBlockPrec = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.NonLinearSolver.MinSolverIterations = 3;

            return C;
        }

        public static XNSERO_Control Channel(int k = 2, double particleLength = 0.3, double aspectRatio = 0.4, int cellsPerUnitLength = 20) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "2_active_Rods", IsRestart: false, PathToOldSessionDir: @"D:\BoSSS_databases\Channel\sessions\85c5edac-706a-4f5d-b0a6-6aa1f489efac");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            //C.SetSaveOptions(dataBasePath: @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\packedParticles", savePeriod: 1);
            string ID = "85c5edac-706a-4f5d-b0a6-6aa1f489efac";
            Guid SID = new Guid(ID);
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));
            int timestep = 220;
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), timestep);
            string PathToOldSessionDir = @"D:\BoSSS_databases\Channel\sessions\85c5edac-706a-4f5d-b0a6-6aa1f489efac";
           
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.LS_TrackerWidth = 2;
            C.TracingNamespaces = "BoSSS.Foundation";
            //C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;

            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 3 });

            // Particle Properties
            // =============================
            double particleDensity = C.PhysicalParameters.rho_A * 1000;
            double activeStress = 100;
            double nextParticleDistance = 0.5 * 2.45;
            double domainLength =3;
            List<string> boundaryValues = new() {
                "Wall_lower",
                "Wall_upper",
                //"Pressure_Outlet_right",
               // "Pressure_Outlet_left"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(domainLength, 4, cellsPerUnitLength, true, false);
            C.CoefficientOfRestitution = 1;
            Motion motion = new(particleDensity);
            double leftCorner = -domainLength / 2 + nextParticleDistance / 2;
            double j = -domainLength / 2 + 0.5;
            List<Particle> particles = new List<Particle>();
            while (j * 1 < domainLength / 2)
            {
                double posX = j * 1;
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, 2.333 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, 2.000 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, 1.666 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, 1.333 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, 1.5 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, 1 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, 0.5 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, 0.000 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, -0.5 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, -1 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, -1.5 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, -1.333 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, -1.666 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { posX, -2.000 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, -2.333 }, 0, activeStress));
                j += 1;
            }
            C.InitialiseParticles(particles);
           

            //Console.WriteLine(particles[0].Motion.GetPosition(0));
            //Console.WriteLine(particles[1].Motion.GetPosition(0));
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: int.MaxValue);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.UseSchurBlockPrec = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.NonLinearSolver.MinSolverIterations = 2;
            C.SetGravity(new Vector(1, 0));
            C.AgglomerationThreshold = 0.1;
            return C;
        }

        public static XNSERO_Control TWOPI(double angle = 78.75, double angleXAxis = 0, double distance = 2, double halfAxis = 0.5, double aspectRatio = 0.4, double deviation = 0) {
            XNSERO_Control C = new XNSERO_Control(2, "2particleInteractions", "active Particles");

            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/2PI", 1);
            //C.AlternateDbPaths = new[] { new ValueTuple<string, string>(@"/work/scratch/ij83requ/default_bosss_db", ""), new ValueTuple<string, string>(@"U:\default_bosss_db", "") };
            // Domain
            // =============================
            //List<string> boundaryValues = new List<string> {
            // "Wall"
            //};
            //C.SetBoundaries(boundaryValues);
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));
            string ID = "c9f879d9-07da-4f1f-860b-4b2f9fb447c5";
            int timestep = 940;
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), timestep);
            string PathToOldSessionDir = @"D:\BoSSS_databases\Channel\sessions\c9f879d9-07da-4f1f-860b-4b2f9fb447c5";
            
            C.SetGrid2D(lengthX: 10, lengthY: 10, cellsPerUnitLength: 10 , periodicX: true, periodicY: true);
            C.SetAddaptiveMeshRefinement(5);
            double dt = 1e-1;
            //List<string> boundaryValues = new List<string> {
            //    "Pressure_Outlet"
            //};
            //C.SetBoundaries(boundaryValues);
            //double length = 6;
            //double cellsPerUnitLength = 8;
            //C.SetGrid(lengthX: length, lengthY: length, cellsPerUnitLength, periodicX: false, periodicY: false);
            C.AdaptiveMeshRefinement = false;
            C.AMR_startUpSweeps = 5;
            C.LS_TrackerWidth = 3;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            
            double particleDensity = 1000;
            // Particle Properties
            // =============================   
            Motion motion = new(particleDensity);
            List<Particle> particles = new List<Particle>();
            angleXAxis = angleXAxis * Math.PI / 180;
            particles = new List<Particle> {
                new ParticleEllipse(motion, halfAxis, aspectRatio * halfAxis, new double[] { -1, -1 }, 45, 1),
                new ParticleEllipse(motion, halfAxis, aspectRatio * halfAxis, new double[] { 1, 1 }, -90, 1)
            };
            C.InitialiseParticles(particles);
            C.NonLinearSolver.ConvergenceCriterion = 1e-4;
            List<Func<Vector,double, bool>> Contains = new List<Func<Vector, double, bool>>();
            for(int p = 0; p < particles.Count; p++) {
                Contains.Add(particles[p].Contains);
            }
            //tC.SetBoundariesactiveAMRlevelIndicators.Add(new AMRForRigidObject(Contains, 4) { maxRefinementLevel = 3 });
            // Timestepping
            // =============================  
            C.SetTimesteps(dt, noOfTimesteps: int.MaxValue);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            return C;
        }

        public static XNSERO_Control newtonPendel(double angle = 78.75, double angleXAxis = 0, double distance = 2, double halfAxis = 0.39, double aspectRatio = 1, double deviation = 0) {
            XNSERO_Control C = new XNSERO_Control(2, "2particleInteractions", "active Particles");

            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\2particleInteractions", savePeriod: 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/2PI", 1);
            //C.AlternateDbPaths = new[] { new ValueTuple<string, string>(@"/work/scratch/ij83requ/default_bosss_db", ""), new ValueTuple<string, string>(@"U:\default_bosss_db", "") };
            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
             "Pressure_Outlet"
            };
            C.SetBoundaries(boundaryValues);
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));
            string ID = "c9f879d9-07da-4f1f-860b-4b2f9fb447c5";
            int timestep = 940;
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), timestep);
            string PathToOldSessionDir = @"D:\BoSSS_databases\Channel\sessions\c9f879d9-07da-4f1f-860b-4b2f9fb447c5";

            C.SetGrid2D(lengthX: 10, lengthY:10, cellsPerUnitLength: 10, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(5);
            double dt = 1e-1;
            //List<string> boundaryValues = new List<string> {
            //    "Pressure_Outlet"
            //};
            //C.SetBoundaries(boundaryValues);
            //double length = 6;
            //double cellsPerUnitLength = 8;
            //C.SetGrid(lengthX: length, lengthY: length, cellsPerUnitLength, periodicX: false, periodicY: false);
            C.AdaptiveMeshRefinement = false;
            C.AMR_startUpSweeps = 5;
            C.LS_TrackerWidth = 3;
            C.SkipSolveAndEvaluateResidual = false;
            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = true;

            double particleDensity = 1000;
            // Particle Properties
            // =============================   
            Motion motion = new(particleDensity);
            List<Particle> particles = new List<Particle>();
            angleXAxis = angleXAxis * Math.PI / 180;
            particles = new List<Particle> {
                new ParticleEllipse(motion, halfAxis, aspectRatio * halfAxis, new double[] { -1.5, 0 }, 0, 0, new double[] {1,0}),
                new ParticleEllipse(motion, halfAxis, aspectRatio * halfAxis, new double[] { 0, 0 }, 0, 0, new double[] {0,0 }),
                new ParticleEllipse(motion, halfAxis, aspectRatio * halfAxis, new double[] { 1, 0 }, -180, 0)
            };
            C.InitialiseParticles(particles);
            C.NonLinearSolver.ConvergenceCriterion = 1e-4;
            List<Func<Vector, double, bool>> Contains = new List<Func<Vector, double, bool>>();
            for (int p = 0; p < particles.Count; p++) {
                Contains.Add(particles[p].Contains);
            }
            //tC.SetBoundariesactiveAMRlevelIndicators.Add(new AMRForRigidObject(Contains, 4) { maxRefinementLevel = 3 });
            // Timestepping
            // =============================  
            C.SetTimesteps(dt, noOfTimesteps: int.MaxValue);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            return C;
        }

        public static XNSERO_Control Wirbelbett(int k = 3, double particleLength = 0.25, double aspectRatio = 1, int cellsPerUnitLength = 5, double noOfRows = 13) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "2_active_Rods", IsRestart: true, PathToOldSessionDir: @"D:\BoSSS_databases\Channel\sessions\cdfe23b1-6f88-4f4f-83ee-39168c6893ad");
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);
            C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);
            //C.SetSaveOptions(dataBasePath: @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\packedParticles", savePeriod: 1);
            string ID = "cdfe23b1-6f88-4f4f-83ee-39168c6893ad";
            Guid SID = new Guid(ID);
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));
            int timestep = 320;
            C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(ID), timestep);
            string PathToOldSessionDir = @"D:\BoSSS_databases\Channel\sessions\cdfe23b1-6f88-4f4f-83ee-39168c6893ad";

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.IncludeConvection = false;
            C.LS_TrackerWidth = 3;
            C.TracingNamespaces = "BoSSS.Foundation";

            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 3 });

            // Particle Properties
            // =============================
            double particleDensity = C.PhysicalParameters.rho_A * 1000;
            double activeStress = 10;
            double nextParticleDistance = 0.5 * 2.45;
            double domainLength = 1;
            List<string> boundaryValues = new() {
                "Wall_lower",
                "Wall_upper",
                //"Pressure_Dirichlet_right",
                //"Pressure_Dirichlet_left"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(domainLength, 3, cellsPerUnitLength, true, false);
            C.CoefficientOfRestitution = 1;
            Motion motion = new(particleDensity);
            double leftCorner = -domainLength / 2 + nextParticleDistance / 2;
            double j = -domainLength / 2 + 0.5;
            List<Particle> particles = new List<Particle>();
            while ((j) * 1 < domainLength / 2) {
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, 2.333 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, 2.000 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, 1.666 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, 1.333 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1 * 0, 1.000 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, 0.666 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, 0.333 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1 * 0, 0.000 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, -0.333 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, -0.666 }, 0, activeStress));
                particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1 * 0, -1.000 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, -1.333 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, -1.666 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, -2.000 }, 0, activeStress));
                //particles.Add(new ParticleEllipse(motion, particleLength, particleLength * aspectRatio, new double[] { j * 1, -2.333 }, 0, activeStress));
                j += 100000;
            }
            C.InitialiseParticles(particles);


            //Console.WriteLine(particles[0].Motion.GetPosition(0));
            //Console.WriteLine(particles[1].Motion.GetPosition(0));
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: int.MaxValue);
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.UseSchurBlockPrec = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.ConvergenceCriterion = 1e-4;
            C.NonLinearSolver.MinSolverIterations = 2;
            C.SetGravity(new Vector(1, 0));

            return C;
        }

    }

}