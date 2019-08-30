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
    public class HardcodedControl_multipleActiveParticles : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control ActiveRods_noBackroundFlow(int k = 3) {
            FSI_Control C = new FSI_Control();
            // General scaling parameter
            // =============================
            const double BaseSize = 1.0;


            // basic database options
            // =============================
            C.DbPath = @"D:\BoSSS_databases\multipleActiveParticles";
            C.savetodb = true;
            C.saveperiod = 1;
            C.ProjectName = "5activeRods_noBackroundFlow";
            C.ProjectDescription = "Active";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");


            // DG degrees
            // =============================
            C.SetDGdegree(k);


            // Grid 
            // =============================
            //Generating grid
            C.GridFunc = delegate {

                int q = new int(); // #Cells in x-dircetion + 1
                int r = new int(); // #Cells in y-dircetion + 1

                q = 12;
                r = 12;

                double[] Xnodes = GenericBlas.Linspace(-5 * BaseSize, 5 * BaseSize, q + 1);
                double[] Ynodes = GenericBlas.Linspace(-5 * BaseSize, 5 * BaseSize, r + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Wall_left");
                grd.EdgeTagNames.Add(2, "Wall_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Wall_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-5 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-5 * BaseSize)) <= 1.0e-8)
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
            C.RefinementLevel = 4;


            // Boundary conditions
            // =============================
            C.AddBoundaryValue("Wall_left");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Wall_right");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Wall_lower");
            C.AddBoundaryValue("Wall_upper");


            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;//pg/(mum^3)
            C.PhysicalParameters.mu_A = 1;//pg(mum*s)
            C.PhysicalParameters.Material = true;


            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles.Add(new Particle_Ellipsoid(new double[] { 0.0, 0.0 }, startAngl: 0) {
                particleDensity = 1,
                activeStress = 1e3,
                thickness_P = 0.1 * BaseSize,
                length_P = 0.2 * BaseSize,
                useAddaptiveUnderrelaxation = true,
                underrelaxation_factor = 3.0,
                clearSmallValues = true,
                UseAddedDamping = true,
            });

            C.Particles.Add(new Particle_Ellipsoid(new double[] { 2.5, 2.5 }, startAngl: 10) {
                particleDensity = 1,
                activeStress = 1e3,
                thickness_P = 0.1 * BaseSize,
                length_P = 0.2 * BaseSize,
                useAddaptiveUnderrelaxation = true,
                underrelaxation_factor = 3.0,
                clearSmallValues = true,
                UseAddedDamping = true
            });

            C.Particles.Add(new Particle_Ellipsoid(new double[] { 2.5, -2.5 }, startAngl: 50) {
                particleDensity = 1,
                activeStress = 1e3,
                thickness_P = 0.1 * BaseSize,
                length_P = 0.2 * BaseSize,
                useAddaptiveUnderrelaxation = true,
                underrelaxation_factor = 3.0,
                clearSmallValues = true,
                UseAddedDamping = true
            });

            C.Particles.Add(new Particle_Ellipsoid(new double[] { -2.5, 2.5 }, startAngl: 80) {
                particleDensity = 1,
                activeStress = 1e3,
                thickness_P = 0.1 * BaseSize,
                length_P = 0.2 * BaseSize,
                useAddaptiveUnderrelaxation = true,
                underrelaxation_factor = 3.0,
                clearSmallValues = true,
                UseAddedDamping = true
            });

            C.Particles.Add(new Particle_Ellipsoid(new double[] { -2.5, -2.5 }, startAngl: 160) {
                particleDensity = 1,
                activeStress = 1e3,
                thickness_P = 0.1 * BaseSize,
                length_P = 0.2 * BaseSize,
                useAddaptiveUnderrelaxation = true,
                underrelaxation_factor = 3.0,
                clearSmallValues = true,
                UseAddedDamping = true
            });


            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;


            //Initial Values
            // =============================   
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
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.forceAndTorqueConvergenceCriterion = 10;
            C.LSunderrelax = 1.0;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.max_iterations_fully_coupled = 10000;


            // Timestepping
            // =============================  
            C.instationarySolver = true;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000000;
            C.NoOfTimesteps = 1000000;

            // haben fertig...
            // ===============

            return C;
        }
    }
}