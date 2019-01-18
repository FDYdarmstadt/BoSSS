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
using BoSSS.Solution.Multigrid;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver
{
    public class HardcodedControlDeussen : IBM_Solver.HardcodedTestExamples
    {
        public static FSI_Control TestActiveParticle(string _DbPath = null, int k = 2, double VelXBase = 0.0, double stressM = 1e0, double cellAgg = 0.2, int maxCurv = 20, double muA = 1e6, double timestepX = 1e-3)
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

                q = 28;
                r = 28;

                double[] Xnodes = GenericBlas.Linspace(-8 * BaseSize, 8 * BaseSize, q);
                double[] Ynodes = GenericBlas.Linspace(-8 * BaseSize, 8 * BaseSize, r);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Pressure_Outlet_left");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Wall_upper");


                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-8 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-8 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-8 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-8 * BaseSize)) <= 1.0e-8)
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
            C.AddBoundaryValue("Pressure_Outlet_left");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Pressure_Outlet_right");//, "VelocityX", X => 0.0);
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
            C.Particles = new List<Particle>();
            int numOfParticles = 1;
            for (int d = 0; d < numOfParticles; d++)
            {
                C.Particles.Add(new Particle_Ellipsoid(2, 4, new double[] { 0 + 4 * d, 0.0 }, startAngl: 180 * d)
                {
                    radius_P = 1,
                    rho_P = 1.02e0,//pg/(mum^3)
                    includeGravity = false,
                    active_P = true,
                    stress_magnitude_P = stressM,
                    thickness_P = 0.4 * BaseSize,
                    length_P = 1 * BaseSize,
                    underrelaxationFT_constant = false,// set true if you want to define a constant underrelaxation (not recommended)
                    underrelaxation_factor = 0.75,// underrelaxation with [factor * 10^exponent]
                    underrelaxationFT_exponent_min = -1,
                });
            }
            //Define level-set
            Func<double[], double, double> phiComplete = delegate (double[] X, double t)
            {
                //Generating the correct sign
                int exp = C.Particles.Count - 1;
                double ret = Math.Pow(-1, exp);
                //Level-set function depending on #particles
                for (int i = 0; i < C.Particles.Count; i++)
                {
                    ret *= C.Particles[i].phi_P(X, t);
                }
                return ret;
            };


            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;


            //Initial Values
            // =============================   
            C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
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
            C.MaxSolverIterations = 1000;
            C.MinSolverIterations = 1;
            C.NoOfMultigridLevels = 1;
            C.ForceAndTorque_ConvergenceCriterion = 1e-5;
            C.LSunderrelax = 1.0;
            

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.splitting_fully_coupled = true;
            C.max_iterations_fully_coupled = 10000;
            C.includeRotation = true;
            C.includeTranslation = true;


            // Timestepping
            // =============================  
            switch (C.Timestepper_LevelSetHandling)
            {
                case LevelSetHandling.Coupled_Once:
                    C.Timestepper_Mode = FSI_Control.TimesteppingMode.MovingMesh;
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    C.Timestepper_Mode = FSI_Control.TimesteppingMode.MovingMesh;
                    break;

                case LevelSetHandling.LieSplitting:
                    C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
                    break;

                case LevelSetHandling.StrangSplitting:
                    C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
                    break;

                case LevelSetHandling.None:
                    C.Timestepper_Mode = FSI_Control.TimesteppingMode.None;
                    break;
            }
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
    }
}