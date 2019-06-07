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
using System.Linq;
using BoSSS.Solution.NSECommon;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Utils;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedTestExamples {
        public static FSI_Control ParticleInShearFlow(string _DbPath = null, int k = 2, double VelXBase = 0.0) {
            FSI_Control C = new FSI_Control();

            const double BaseSize = 1.0;

            // basic database options
            // ======================

            C.savetodb = false;
            C.ProjectName = "ShearFlow_Test";
            C.ProjectDescription = "ShearFlow";
            C.Tags.Add("with immersed boundary method");

            // DG degrees
            // ==========

            C.SetDGdegree(k);

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {
                double[] Xnodes = GenericBlas.Linspace(-2 * BaseSize, 2 * BaseSize, 21);
                double[] Ynodes = GenericBlas.Linspace(-3 * BaseSize, 3 * BaseSize, 31);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: true);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-2 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-2 * BaseSize)) <= 1.0e-8)
                        et = 2;

                    Debug.Assert(et != 0);
                    return et;
                });

                return grd;
            };

            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityY", X => 0.02);
            C.AddBoundaryValue("Velocity_Inlet_right", "VelocityY", X => -0.02);

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============
            // Coupling Properties
            //C.LevelSetMovement = "coupled";
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            // Particle Properties
            C.Particles.Add(new Particle_Sphere(new double[] { 0.0, 0.0 }) {
                radius_P = 0.4,
                particleDensity = 1.0,
                IncludeTranslation = false,
                IncludeRotation = true
            });

            

            ////Define level-set
            //Func<double[], double, double> phiComplete = delegate (double[] X, double t) {
            //    int exp = C.Particles.Count - 1;
            //    double ret = Math.Pow(-1, exp);
            //    for (int i = 0; i < C.Particles.Count; i++) {
            //        ret *= C.Particles[i].Phi_P(X, t);
            //    }
            //    return ret;
            //};

            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;

            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            //C.InitialValues.Add("VelocityX#B", X => 1);
            //C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            //C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("fec14187-4e12-43b6-af1e-e9d535c78668"), -1);


            // Physical Parameters
            // ===================

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.25;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================
            C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.LevelSetSmoothing = false;
            C.LinearSolver.SolverCode = LinearSolverConfig.Code.classic_pardiso;
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
            C.NoOfTimesteps = 100;

           

            // haben fertig...
            // ===============

            return C;
        }

        /// <summary>
        /// Testing of particle/wall interactions using a single particle
        /// </summary>
        public static FSI_Control SingleDryParticleAgainstWall(string _DbPath = null, bool MeshRefine = false) {
            FSI_Control C = new FSI_Control();

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.saveperiod = 1;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


            // DG degrees
            // ==========

            C.SetDGdegree(1);

            // grid and boundary conditions
            // ============================

            double[] Xnodes = GenericBlas.Linspace(-1, 1, 31);
            double[] Ynodes = GenericBlas.Linspace(-1, 1, 31);
            double h = Math.Min((Xnodes[1] - Xnodes[0]), (Ynodes[1] - Ynodes[0]));

            C.GridFunc = delegate {
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);
                grd.EdgeTagNames.Add(1, "Wall");
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 1;
                    return et;
                });

                return grd;
            };

            C.AddBoundaryValue("Wall");

            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;


            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;

            // Particles
            // =========
            C.Particles.Add(new Particle_Sphere(new double[] { -0.5, -0.5 }, startAngl: 90.0)
            {
                particleDensity = 1.0,
                radius_P = 0.1,
            });
            C.Particles[0].TranslationalVelocity[0][0] = +1;
            C.Particles[0].TranslationalVelocity[0][1] = -1;
            C.Particles[0].RotationalVelocity[0] = 0;
            C.pureDryCollisions = true;
            C.collisionModel = FSI_Control.CollisionModel.MomentumConservation;

            double V = 0;
            foreach (var p in C.Particles) {
                V = Math.Max(V, p.TranslationalVelocity[0].L2Norm());
            }

            if (V <= 0)
                throw new ArithmeticException();


            // Physical Parameters
            // ===================

            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.AdaptiveMeshRefinement = MeshRefine;
            C.RefinementLevel = 1;

            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;

            double dt = (h / V) * (MeshRefine ? 0.5 * 0.5 * 0.5 * 0.2 : 0.1);
            C.dtMax = dt;
            C.dtMin = dt;

            C.Endtime = 100.0 / V;
            C.NoOfTimesteps = 500;

            // haben fertig...
            // ===============

            return C;
        }

        /// <summary>
        /// Testing of particle/wall interactions using a single particle
        /// </summary>
        public static FSI_Control DryParticleCollision(string _DbPath = null, bool MeshRefine = false) {
            FSI_Control C = new FSI_Control();

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.saveperiod = 1;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


            // DG degrees
            // ==========

            C.SetDGdegree(1);

            // grid and boundary conditions
            // ============================

            double[] Xnodes = GenericBlas.Linspace(-1, 1, 21);
            double[] Ynodes = GenericBlas.Linspace(-1, 1, 21);
            double h = Math.Min((Xnodes[1] - Xnodes[0]), (Ynodes[1] - Ynodes[0]));

            C.GridFunc = delegate {
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);
                grd.EdgeTagNames.Add(1, "Wall");
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 1;
                    return et;
                });

                return grd;
            };

            C.AddBoundaryValue("Wall");

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;


            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;


            // Particles
            // =========

            C.Particles.Add(new Particle_Sphere(new double[] { -0.6, +0.1 }, startAngl: 90.0) {
                particleDensity = 1.0,
                radius_P = 0.15
            });
            C.Particles[0].TranslationalVelocity[0][0] = +1;
            C.Particles[0].TranslationalVelocity[0][1] = 0;
            C.Particles[0].RotationalVelocity[0] = 0;

            C.Particles.Add(new Particle_Sphere(new double[] { +0.6, -0.1 }, startAngl: 90.0) {
                particleDensity = 1.0,
                radius_P = 0.15
            });
            C.Particles[1].TranslationalVelocity[0][0] = -1;
            C.Particles[1].TranslationalVelocity[0][1] = 0;
            C.Particles[1].RotationalVelocity[0] = 0;
            
            C.collisionModel = FSI_Control.CollisionModel.MomentumConservation;

            double V = 0;
            foreach (var p in C.Particles) {
                V = Math.Max(V, p.TranslationalVelocity[0].L2Norm());
            }

            if (V <= 0)
                throw new ArithmeticException();


            // Physical Parameters
            // ===================

            C.pureDryCollisions = true;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.AdaptiveMeshRefinement = MeshRefine;

            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;

            double dt = (h / V) * (MeshRefine ? 0.5 * 0.5 * 0.5 * 0.2 : 0.1);
            C.dtMax = dt;
            C.dtMin = dt;

            C.Endtime = 100.0 / V;
            C.NoOfTimesteps = 200;

            // haben fertig...
            // ===============

            return C;

        }

        /// <summary>
        /// Testing particle bouncing
        /// </summary>
        public static FSI_Control DryParticleBounce(string _DbPath = null)
        {
            FSI_Control C = new FSI_Control();

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.saveperiod = 1;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


            // DG degrees
            // ==========

            C.SetDGdegree(1);

            // grid and boundary conditions
            // ============================

            double[] Xnodes = GenericBlas.Linspace(-1, 1, 20);
            double[] Ynodes = GenericBlas.Linspace(-1, 2, 30);
            double h = Math.Min((Xnodes[1] - Xnodes[0]), (Ynodes[1] - Ynodes[0]));

            C.GridFunc = delegate {
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);
                grd.EdgeTagNames.Add(1, "Wall");
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 1;
                    return et;
                });

                return grd;
            };

            C.AddBoundaryValue("Wall");

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;


            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;


            // Particles
            // =========

            C.Particles.Add(new Particle_Sphere(new double[] { 0.0, 0.8 }, startAngl: 0.0)
            {
                particleDensity = 1.0,
                radius_P = 0.15,
                GravityVertical = -9.81
            });

            C.collisionModel = FSI_Control.CollisionModel.MomentumConservation;


            // Physical Parameters
            // ===================

            C.pureDryCollisions = true;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.AdaptiveMeshRefinement = false;

            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;

            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;

            C.Endtime = 2;
            C.NoOfTimesteps = 111;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control StickyTrap(string _DbPath = null, int k = 2, double VelXBase = 0.0, double angle = 0.0)
        {
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

            C.pureDryCollisions = true;
            C.SetDGdegree(k);

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate
            {

                int q, r;

                q = 50;
                r = 50;

                double[] Xnodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, q + 1);
                double[] Ynodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, r + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Wall_left");
                grd.EdgeTagNames.Add(2, "Wall_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet_upper");


                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 2;

                    if (Math.Abs(X[1] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 4;


                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };

            C.GridPartType = GridPartType.Hilbert;

            C.AddBoundaryValue("Wall_left");
            C.AddBoundaryValue("Wall_right");
            C.AddBoundaryValue("Wall_lower");
            C.AddBoundaryValue("Pressure_Outlet_upper");

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
            C.CoefficientOfRestitution = 0;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleMass = 1;


            C.Particles.Add(new Particle_Sphere(new double[] { 0.0, 0.7 })
            {
                radius_P = 0.13,
                //length_P = 0.2,
                //thickness_P = 0.1,
                particleDensity = 2.0,
                GravityVertical = -9.81,
                //AddaptiveUnderrelaxation = true,
                //underrelaxation_factor = 1,
                //ClearSmallValues = true,
                //neglectAddedDamping = false,
            });

            C.Particles.Add(new Particle_superEllipsoid(new double[] { 0.45, 0 }, startAngl: 45)
            {
                particleDensity = 1,
                thickness_P = 0.2,
                length_P = 0.4,
                //radius_P = 0.4,
                superEllipsoidExponent = 4,
                GravityVertical = -9.81,
                IncludeRotation = false,
                IncludeTranslation = false,
            });


            C.Particles.Add(new Particle_superEllipsoid(new double[] { -0.45, 0 }, startAngl: -45)
            {
                particleDensity = 1,
                thickness_P = 0.2,
                length_P = 0.4,
                //radius_P = 0.4,
                superEllipsoidExponent = 4,
                GravityVertical = -9.81,
                IncludeRotation = false,
                IncludeTranslation = false,
            });
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


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 38;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control ActiveParticle_ForceTest(int k = 2)
        {
            FSI_Control C = new FSI_Control();

            // basic database options
            // =============================
            //C.DbPath = @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\straightChannel";
            C.savetodb = false;
            C.saveperiod = 1;
            C.ProjectName = "Test_singleActiveParticle";
            C.ProjectDescription = "Test_singleActiveParticle";
            C.SessionName = C.ProjectName;
            C.Tags.Add("activeParticle");
            
            // DG degrees
            // =============================
            C.SetDGdegree(k);
            
            // Grid 
            // =============================
            //Generating grid
            C.GridFunc = delegate
            {
                int q = new int(); // #Cells in x-dircetion + 1
                int r = new int(); // #Cells in y-dircetion + 1

                q = 40;
                r = 30;

                double[] Xnodes = GenericBlas.Linspace(-4, 4, q);
                double[] Ynodes = GenericBlas.Linspace(-3, 3, r);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Pressure_Outlet_left");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Wall_upper");


                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-4)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-4)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-3)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-3)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };


            // Mesh refinement
            // =============================
            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 2;
            C.maxCurvature = 2;


            // Boundary conditions
            // =============================
            C.AddBoundaryValue("Pressure_Outlet_left");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Pressure_Outlet_right");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Wall_lower");
            C.AddBoundaryValue("Wall_upper");


            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;//pg/(mum^3)
            C.PhysicalParameters.mu_A = 1e4;//pg(mum*s)
            C.PhysicalParameters.Material = true;


            // Particle Properties
            // =============================   
            C.Particles = new List<Particle>();
            int numOfParticles = 1;
            for (int d = 0; d < numOfParticles; d++)
            {
                C.Particles.Add(new Particle_Ellipsoid(new double[] { 0.0, 0.0 }, startAngl: 0)
                {
                    particleDensity = 1,
                    ActiveParticle = true,
                    ActiveStress = 1e5,
                    thickness_P = 0.4,
                    length_P = 1,
                    AddaptiveUnderrelaxation = true,// set true if you want to define a constant underrelaxation (not recommended)
                    underrelaxation_factor = 1,// underrelaxation with [factor * 10^exponent]
                    ClearSmallValues = true,
                    neglectAddedDamping = false,
                    IncludeRotation = false,
                    IncludeTranslation = false
                });
            }
            //Define level-set
            double phiComplete(double[] X, double t)
            {
                //Generating the correct sign
                int exp = C.Particles.Count - 1;
                double ret = Math.Pow(-1, exp);
                //Level-set function depending on # of particles
                for (int i = 0; i < C.Particles.Count; i++)
                {
                    ret *= C.Particles[i].Phi_P(X);
                }
                return ret;
            }
            
            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            
            //Initial Values
            // =============================   
            C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);

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
            C.ForceAndTorque_ConvergenceCriterion = 1e-2;
            C.LSunderrelax = 1.0;
            
            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.max_iterations_fully_coupled = 250;



            // Timestepping
            // =============================  
            C.instationarySolver = true;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e-3;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }
    }
}
