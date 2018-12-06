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
        public static FSI_Control TestActiveParticle(string _DbPath = null, int k = 2, double VelXBase = 0.0, double stressM = 1 , double cellAgg = 0.2, int maxCurv = 20, double muA = 1e3)
        {
            FSI_Control C = new FSI_Control();


            // General scaling parameter
            // =============================
            const double BaseSize = 1.0;


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
            C.FieldOptions.Add("VelocityX", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts()
            {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });


            // Grid 
            // =============================
            //Generating grid
            C.GridFunc = delegate {

                int q = new int();
                int r = new int();

                switch (99)
                {
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;

                    case 99:
                        q = 20;
                        r = 100;
                        break;

                    default:

                        throw new ApplicationException();
                }
                
                double[] Xnodes = GenericBlas.Linspace(-1 * BaseSize, 1 * BaseSize, q); 
                double[] Ynodes = GenericBlas.Linspace(-10 * BaseSize, 0 * BaseSize, r); 

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Pressure_Outlet_left");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-10 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-0 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };


            //Mesh refinement
            // =============================
            C.AdaptiveMeshRefinement = true;
            C.RefinementLevel = 2;
            C.maxCurvature = maxCurv;


            // Boundary conditions
            // =============================
            C.AddBoundaryValue("Pressure_Outlet_left");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Pressure_Outlet_right");//, "VelocityX", X => 0.0);
            C.AddBoundaryValue("Wall_lower");
            C.AddBoundaryValue("Pressure_Outlet_upper");
            

            // Coupling Properties
            // =============================
            //C.LevelSetMovement = "coupled";
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.includeRotation = true;
            C.includeTranslation = true;


            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 0.9982;//pg/(mum^3)
            C.PhysicalParameters.mu_A = muA;//pg(mum*mus)
            C.PhysicalParameters.Material = true;


            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            int numOfParticles = 1;
            for (int d = 0; d < numOfParticles; d++)
            {
                C.Particles.Add(new Particle(2, 9, new double[] { 0 + 14.0 * d, -1.0 }, startAngl: 180.0*d, shape: Particle.ParticleShape.elliptic)
                //Generates a series of opposing particles
                {
                    radius_P = 1,
                    rho_P = 1.1,//pg/(mum^3)
                    includeGravity = true,
                    active_P = false,
                    stress_magnitude_P = stressM,
                    thickness_P = 0.5,
                    length_P = 0.5,
                    C_v = 10000,
                    velResidual_ConvergenceCriterion = 1e-15
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
            C.LevelSet_ConvergenceCriterion = 1e-8;
            C.LSunderrelax = 1.0;


            // Timestepping
            // =============================  
            C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-6;//ms
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10000;
            C.NoOfTimesteps = 200;
            
            

            // haben fertig...
            // ===============

            return C;

        }

        public static FSI_Control Test2ActiveParticle(string _DbPath = null, int k = 2, double VelXBase = 0.0, double stressM = 100)
        {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;

            //C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
            //                    new Tuple<string,object>("k", k),
            //                };

            // k = i;

            // basic database options
            // ======================
            //C.DbPath = @"\\hpccluster\hpccluster-scratch\deussen\cluster_db\active_particle_test";
            C.savetodb = false;
            C.saveperiod = 5;
            C.ProjectName = "Active2ParticleTest";
            C.ProjectDescription = "Active2";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;
            C.RefinementLevel = 2;


        // DG degrees
        // ==========

        C.FieldOptions.Add("VelocityX", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts()
            {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });


            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {

                int q = new int();
                int r = new int();

                switch (99)
                {
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;

                    case 99:
                        q = 80;
                        r = 25;
                        break;

                    default:

                        throw new ApplicationException();
                }

                //q = 16;
                //r = 46;

                double[] Xnodes = GenericBlas.Linspace(-5 * BaseSize, 15 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(-3 * BaseSize, 3 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-15 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-3 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-3 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };




            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityX", X => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_right", "VelocityX", X => 0.0);
            C.AddBoundaryValue("Pressure_Outlet_lower");
            C.AddBoundaryValue("Pressure_Outlet_upper");

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeRotation = true;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.25;
            C.PhysicalParameters.Material = true;

            //Defining particles
            C.Particles = new List<Particle>();
            int numOfParticles = 2;            
            for (int d = 0; d < numOfParticles; d++)
            {
                C.Particles.Add(new Particle(2, 4, new double[] { 10.0*d, 0.0 }, startAngl: d*180.0, shape: Particle.ParticleShape.elliptic)
                {
                  radius_P = 1,
                    rho_P = 1.0,
                   active_P = true,
                   stress_magnitude_P = stressM,
                   thickness_P = 1.0,
                   length_P = 3.0
                 });
            }
            //Define level-set
            Func<double[], double, double> phiComplete = delegate (double[] X, double t)
            {
                int exp = C.Particles.Count - 1;
                double ret = Math.Pow(-1, exp);
                for (int i = 0; i < C.Particles.Count; i++)
                {
                    ret *= C.Particles[i].phi_P(X, t);
                }
                return ret;
            };

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            

            //Func<double[], double, double> phiComplete = (X, t) => -1 * (C.Particles[0].phi_P(X, t)) * (C.Particles[1].phi_P(X, t));
            //Func<double[], double, double> phiComplete = (X, t) => -1 * (C.Particles[0].phi_P(X, t) * C.Particles[1].phi_P(X, t));
            //for (int i = 0;i<C.Particles.Count; i++) {
            //    phiComplete = (X,t) => phiComplete(X,t)*C.Particles[i].phi_P(X,t);
            //}


            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;         

            C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            //C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("42c82f3c-bdf1-4531-8472-b65feb713326"), 400);
            //C.GridGuid = new Guid("f1659eb6 -b249-47dc-9384-7ee9452d05df");


            // Physical Parameters
            // ===================

            C.PhysicalParameters.IncludeConvection = true;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxSolverIterations = 100;
            C.NoOfMultigridLevels = 1;


            // Timestepping
            // ============

            C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 100000;

            // haben fertig...
            // ===============

            return C;

        }
    }
}