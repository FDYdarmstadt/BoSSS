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

namespace BoSSS.Application.FSI_Solver
{
    public class HardcodedControlDeussen : IBM_Solver.HardcodedTestExamples
    {
        public static FSI_Control TestActiveParticle(string _DbPath = null, int k = 2, double VelXBase = 0.0, double stressM = 1e0, double cellAgg = 0.2, int maxCurv = 20, double muA = 1e0, double timestepX = 1e-3)
        {
            FSI_Control C = new FSI_Control();


            // General scaling parameter
            // =============================
            const double BaseSize = 1.0e0;


            // basic database options
            // =============================
            //C.DbPath = @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\active_particle_test";
            C.savetodb = false;
            C.saveperiod = 1;
            C.ProjectName = "ActiveParticleTest";
            C.ProjectDescription = "Active";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");


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

                q = 20;
                r = 8;

                double[] Xnodes = GenericBlas.Linspace(-5 * BaseSize, 5 * BaseSize, q);
                double[] Ynodes = GenericBlas.Linspace(-2 * BaseSize, 2 * BaseSize, r);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Wall_left");
                grd.EdgeTagNames.Add(2, "Wall_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Wall_upper");


                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-2 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-2 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };


            // Mesh refinement
            // =============================
            C.AdaptiveMeshRefinement = true;
            C.RefinementLevel =2;
            C.maxCurvature = maxCurv;


            // Boundary conditions
            // =============================
            C.AddBoundaryValue("Wall_left");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Wall_right");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Wall_lower");
            C.AddBoundaryValue("Wall_upper");
            

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 0.9982e0;//pg/(mum^3)
            C.PhysicalParameters.mu_A = muA;//pg(mum*s)
            C.PhysicalParameters.Material = true;


            // Particle Properties
            // =============================   
            // Defining particles
            int numOfParticles = 1;
            for (int d = 0; d < numOfParticles; d++)
            {
                C.Particles.Add(new Particle_Sphere(new double[] { 0 + 14 * d, 0.0 }, startAngl: 180 * d)
                {
                    radius_P = 1,
                    particleDensity = 1.5,//pg/(mum^3)
                    GravityVertical = 0,
                    ActiveParticle = true,
                    ActiveStress = stressM,
                    //thickness_P = 0.1 * BaseSize,  Sphere kann nur einen radius haben! fk.
                    //length_P = 2 * BaseSize,       Sphere kann nur einen radius haben! fk.
                    //superEllipsoidExponent = 4, // only even numbers are supported
                    useAddaptiveUnderrelaxation = true,// set true if you want to define a constant underrelaxation (not recommended)
                    underrelaxation_factor = 9,// underrelaxation with [factor * 10^exponent]
                });
            }
            //Define level-set
            //Func<double[], double, double> phiComplete = delegate (double[] X, double t)
            //{
            //    //Generating the correct sign
            //    int exp = C.Particles.Count - 1;
            //    double ret = Math.Pow(-1, exp);
            //    //Level-set function depending on #particles
            //    for (int i = 0; i < C.Particles.Count; i++)
            //    {
            //        ret *= C.Particles[i].Phi_P(X);
            //    }
            //    return ret;
            //};


            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;


            //Initial Values
            // =============================   
            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);


            // For restart
            // =============================  
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("42c82f3c-bdf1-4531-8472-b65feb713326"), 400);
            //C.GridGuid = new Guid("f1659eb6 -b249-47dc-9384-7ee9452d05df");


            // Physical Parameters
            // =============================  
            C.PhysicalParameters.IncludeConvection = false;


            // misc. solver options
            // =============================  
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = cellAgg;
            C.LevelSetSmoothing = false;
            //C.MaxSolverIterations = 1000;
            //C.MinSolverIterations = 1;
            //C.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;

            C.ForceAndTorque_ConvergenceCriterion = 1e-6;
            C.LSunderrelax = 1.0;
            

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.max_iterations_fully_coupled = 10000;
            
            // Timestepping
            // =============================  
            //switch (C.Timestepper_LevelSetHandling)
            //{
            //    case LevelSetHandling.Coupled_Once:
            //        C.Timestepper_Mode = FSI_Control.TimesteppingMode.MovingMesh;
            //        break;

            //    case LevelSetHandling.Coupled_Iterative:
            //        C.Timestepper_Mode = FSI_Control.TimesteppingMode.MovingMesh;
            //        break;

            //    case LevelSetHandling.LieSplitting:
            //        C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            //        break;

            //    case LevelSetHandling.StrangSplitting:
            //        C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            //        break;

            //    case LevelSetHandling.None:
            //        C.Timestepper_Mode = FSI_Control.TimesteppingMode.None;
            //        break;
            //}
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = timestepX;//s
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000000;
            C.NoOfTimesteps = 10000;
            
            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control RefinementTest(string _DbPath = null, bool MeshRefine = true)
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

            double[] Xnodes = GenericBlas.Linspace(-1, 1, 5);
            double[] Ynodes = GenericBlas.Linspace(-1, 1, 5);

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
            C.Particles.Add(new Particle_Sphere(new double[] { -0.5, -0.5 }, startAngl: 90.0)
            {
                particleDensity = 1.0,
                radius_P = 0.4
            });

            C.pureDryCollisions = true;
            C.collisionModel = FSI_Control.CollisionModel.MomentumConservation;


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

            double dt = 0.01;
            C.dtMax = dt;
            C.dtMin = dt;

            C.Endtime = 0.02;
            C.NoOfTimesteps = 2;

            // haben fertig...
            // ===============

            return C;

        }
    }
}