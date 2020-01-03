///* =======================================================================
//Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//*/

using System;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using System.Collections.Generic;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedControlDeriabina : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control RectangleTest(int k = 3) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "9_active_Rods");
            //C.SetSaveOptions(@"D:\BoSSS_databases\multipleActiveParticles", 1);

            List<string> boundaryValues = new List<string> {
                "Wall_left",
                "Wall_right",
                "Wall_lower",
                "Wall_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 3, lengthY: 3, cellsPerUnitLength: 30, periodicX: false, periodicY: false);
            //C.SetAddaptiveMeshRefinement(amrLevel: 2);
            C.hydrodynamicsConvergenceCriterion = 1e-2;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = true;
            C.pureDryCollisions = true;
            C.gravity = new Vector( 0, -9.81 );

            // Particle Properties
            // =============================
            double particleDensity = 20;
            ParticleMotionInit motion1 = new ParticleMotionInit(C.gravity, particleDensity, C.pureDryCollisions, false, false);
            ParticleMotionInit motion2 = new ParticleMotionInit(C.gravity, particleDensity, C.pureDryCollisions, true, true);
            C.Particles.Add(new Particle_Shell(motion2, 1, 0.5, 0.2, new double[] { 0, 0}, startAngl: 0));
            C.Particles.Add(new Particle_Sphere(motion1, 0.1, new double[] { 0, 1 }, startAngl: 0));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000000;
            C.NoOfTimesteps = 1000000;
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
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 1000000;

            return C;
        }
        public static FSI_Control DeriabinaPentagoneFalle(int k = 2) {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;


            // basic database options
            // ======================

            //C.DbPath = @"\\dc1\userspace\deriabina\bosss_db";
            C.savetodb = false;
            C.saveperiod = 1;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = false;
            C.SessionName = "fjkfjksdfhjk";
            C.RefinementLevel = 3;

            C.pureDryCollisions = true;
            C.SetDGdegree(k);

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {

                int q = new int();
                int r = new int();

                q = 10;
                r = 10;


                double[] Xnodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, q + 1);
                double[] Ynodes = GenericBlas.Linspace(3.5 * BaseSize, 6.0 * BaseSize, r + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Wall_left");
                grd.EdgeTagNames.Add(2, "Wall_right");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                grd.EdgeTagNames.Add(4, "Wall_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (3.5 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-6.0 * BaseSize)) <= 1.0e-8)
                        et = 4;


                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };

            C.GridPartType = GridPartType.Hilbert;

            C.AddBoundaryValue("Wall_left");
            C.AddBoundaryValue("Wall_right");
            C.AddBoundaryValue("Pressure_Outlet");
            C.AddBoundaryValue("Wall_upper");

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 1.0;
            C.gravity = new Vector( 0, -9.81 );
            double particleDensity = 2.01;

            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, C.pureDryCollisions);
            ParticleMotionInit fix = new ParticleMotionInit(C.gravity, particleDensity, C.pureDryCollisions, true, true);

            C.Particles.Add(new Particle_Sphere(motion, 0.35, new double[] { -1.0, 5.5 }) {
            });

            C.Particles.Add(new Particle_superEllipsoid(fix, 1, 0.3, 4, new double[] { 0.60, 4.0 }, startAngl: 45) {
            });

            C.Particles.Add(new Particle_superEllipsoid(fix, 1, 0.3, 4, new double[] { -0.60, 4.0 }, startAngl: -45) {
            });

            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);

            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("42c82f3c-bdf1-4531-8472-b65feb713326"), 400);
            //C.GridGuid = new Guid("f1659eb6-b249-47dc-9384-7ee9452d05df");


            // Physical Parameters
            // ===================

            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 30.0;
            C.NoOfTimesteps = 1000000;

            // haben fertig...
            // ===============

            return C;
        }
        //        public static FSI_Control DeriabinaHefezelleWORefinement(int k = 2)
        //        {
        //            FSI_Control C = new FSI_Control();


        //            const double BaseSize = 1.0;


        //            // basic database options
        //            // ======================

        //            C.DbPath = @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\Deriabina";
        //            C.savetodb = false;
        //            C.saveperiod = 1;
        //            C.ProjectName = "ParticleCollisionTest";
        //            C.ProjectDescription = "Gravity";
        //            C.SessionName = C.ProjectName;
        //            C.Tags.Add("with immersed boundary method");
        //            C.AdaptiveMeshRefinement = false;
        //            C.SessionName = "abc";
        //            C.RefinementLevel = 3;

        //            C.pureDryCollisions = false;
        //            C.SetDGdegree(k);

        //            // grid and boundary conditions
        //            // ============================

        //            C.GridFunc = delegate
        //            {

        //                int q = new int();
        //                int r = new int();

        //                q = 35;
        //                r = 140;

        //                double[] Xnodes = GenericBlas.Linspace(-1.0 * BaseSize, 1.0 * BaseSize, q + 1);
        //                double[] Ynodes = GenericBlas.Linspace(2 * BaseSize, 10 * BaseSize, r + 1);

        //                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

        //                grd.EdgeTagNames.Add(1, "Wall_left");
        //                grd.EdgeTagNames.Add(2, "Wall_right");
        //                grd.EdgeTagNames.Add(3, "Pressure_Outlet");
        //                grd.EdgeTagNames.Add(4, "Wall_upper");


        //                grd.DefineEdgeTags(delegate (double[] X)
        //                {
        //                    byte et = 0;
        //                    if (Math.Abs(X[0] - (-1.0 * BaseSize)) <= 1.0e-8)
        //                        et = 1;
        //                    if (Math.Abs(X[0] + (-1.0 * BaseSize)) <= 1.0e-8)
        //                        et = 2;
        //                    if (Math.Abs(X[1] - (2 * BaseSize)) <= 1.0e-8)
        //                        et = 3;
        //                    if (Math.Abs(X[1] + (-10 * BaseSize)) <= 1.0e-8)
        //                        et = 4;


        //                    return et;
        //                });

        //                Console.WriteLine("Cells:" + grd.NumberOfCells);

        //                return grd;
        //            };

        //            C.GridPartType = GridPartType.Hilbert;

        //            C.AddBoundaryValue("Wall_left");
        //            C.AddBoundaryValue("Wall_right");
        //            C.AddBoundaryValue("Pressure_Outlet");
        //            C.AddBoundaryValue("Wall_upper");



        //            // Initial Values
        //            // ==============

        //            // Coupling Properties
        //            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;

        //            // Fluid Properties
        //            C.PhysicalParameters.rho_A = 1.0;
        //            C.PhysicalParameters.mu_A = 0.01;
        //            C.CoefficientOfRestitution = 1;


        //            C.Particles.Add(new Particle_Sphere( 0.1, new double[] { 0.0, 9.5 })
        //            {
        //                particleDensity = 1.01,
        //                useAddaptiveUnderrelaxation = true,
        //                underrelaxation_factor = 3.0,
        //                clearSmallValues = true,
        //                UseAddedDamping = true
        //            });

        //            C.Particles.Add(new Particle_Sphere( 0.1, new double[] { 0.0, 9.1 })
        //            {
        //                particleDensity = 1.01,
        //                useAddaptiveUnderrelaxation = true,
        //                underrelaxation_factor = 3.0,
        //                clearSmallValues = true,
        //                UseAddedDamping = true
        //            });

        //            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
        //            C.InitialValues_Evaluators.Add("VelocityY", X => 0);


        //            // Physical Parameters
        //            // ===================

        //            C.PhysicalParameters.IncludeConvection = false;

        //            // misc. solver options
        //            // ====================

        //            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
        //            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
        //            C.LevelSetSmoothing = false;
        //            C.LinearSolver.MaxSolverIterations = 10;
        //            C.NonLinearSolver.MaxSolverIterations = 10;
        //            C.LinearSolver.NoOfMultigridLevels = 1;
        //            C.forceAndTorqueConvergenceCriterion = 5e-3;


        //            // Timestepping
        //            // ============

        //            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
        //            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
        //            double dt = 5e-4;
        //            C.dtMax = dt;
        //            C.dtMin = dt;
        //            C.Endtime = 1000000.0;
        //            C.NoOfTimesteps = 1000000000;

        //            // haben fertig...
        //            // ===============

        //            return C;
        //        }

        //        public static FSI_Control DeriabinaHefezelleWRefinement(int k = 2)
        //        {
        //            FSI_Control C = new FSI_Control();


        //            const double BaseSize = 1.0;


        //            // basic database options
        //            // ======================

        //            C.DbPath = @"D:\BoSSS_databases\wetParticleCollision";
        //            C.saveperiod = 1;
        //            C.ProjectName = "ParticleCollisionTest";
        //            C.ProjectDescription = "Gravity";
        //            C.SessionName = C.ProjectName;
        //            C.Tags.Add("with immersed boundary method");
        //            C.AdaptiveMeshRefinement = true;
        //            C.SessionName = "fjkfjksdfhjk";
        //            C.RefinementLevel = 3;

        //            C.pureDryCollisions = false;
        //            C.SetDGdegree(k);

        //            // grid and boundary conditions
        //            // ============================

        //            C.GridFunc = delegate
        //            {

        //                int q = new int();
        //                int r = new int();

        //                r = 80;
        //                q = r / 4;

        //                double[] Xnodes = GenericBlas.Linspace(-1.0 * BaseSize, 1.0 * BaseSize, q + 1);
        //                double[] Ynodes = GenericBlas.Linspace(2 * BaseSize, 10 * BaseSize, r + 1);

        //                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

        //                grd.EdgeTagNames.Add(1, "Wall_left");
        //                grd.EdgeTagNames.Add(2, "Wall_right");
        //                grd.EdgeTagNames.Add(3, "Pressure_Outlet");
        //                grd.EdgeTagNames.Add(4, "Wall_upper");


        //                grd.DefineEdgeTags(delegate (double[] X)
        //                {
        //                    byte et = 0;
        //                    if (Math.Abs(X[0] - (-1.0 * BaseSize)) <= 1.0e-8)
        //                        et = 1;
        //                    if (Math.Abs(X[0] + (-1.0 * BaseSize)) <= 1.0e-8)
        //                        et = 2;
        //                    if (Math.Abs(X[1] - (2 * BaseSize)) <= 1.0e-8)
        //                        et = 3;
        //                    if (Math.Abs(X[1] + (-10 * BaseSize)) <= 1.0e-8)
        //                        et = 4;


        //                    return et;
        //                });

        //                Console.WriteLine("Cells:" + grd.NumberOfCells);

        //                return grd;
        //            };

        //            C.GridPartType = GridPartType.Hilbert;

        //            C.AddBoundaryValue("Wall_left");
        //            C.AddBoundaryValue("Wall_right");
        //            C.AddBoundaryValue("Pressure_Outlet");
        //            C.AddBoundaryValue("Wall_upper");



        //            // Initial Values
        //            // ==============

        //            // Coupling Properties
        //            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;

        //            // Fluid Properties
        //            C.PhysicalParameters.rho_A = 1.0;
        //            C.PhysicalParameters.mu_A = 0.01;
        //            C.CoefficientOfRestitution = 1;


        //            C.Particles.Add(new Particle_Sphere( 0.1, new double[] { 0.0, 9.5 })
        //            {
        //                particleDensity = 1.01,
        //                useAddaptiveUnderrelaxation = true,
        //                underrelaxation_factor = 3.0,
        //                clearSmallValues = true,
        //                UseAddedDamping = true
        //            });

        //            C.Particles.Add(new Particle_Sphere( 0.1, new double[] { 0.0, 9.1 })
        //            {
        //                particleDensity = 1.01,
        //                useAddaptiveUnderrelaxation = true,
        //                underrelaxation_factor = 3.0,
        //                clearSmallValues = true,
        //                UseAddedDamping = true
        //            });

        //            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
        //            C.InitialValues_Evaluators.Add("VelocityY", X => 0);


        //            // Physical Parameters
        //            // ===================

        //            C.PhysicalParameters.IncludeConvection = false;

        //            // misc. solver options
        //            // ====================

        //            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
        //            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
        //            C.LevelSetSmoothing = false;
        //            C.LinearSolver.MaxSolverIterations = 10;
        //            C.NonLinearSolver.MaxSolverIterations = 10;
        //            C.LinearSolver.NoOfMultigridLevels = 1;
        //            C.forceAndTorqueConvergenceCriterion = 2e-3;


        //            // Timestepping
        //            // ============

        //            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
        //            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
        //            double dt = 1e-3;
        //            C.dtMax = dt;
        //            C.dtMin = dt;
        //            C.Endtime = 1000000.0;
        //            C.NoOfTimesteps = 1000000000;

        //            // haben fertig...
        //            // ===============

        //            return C;
        //        }


    }
}
