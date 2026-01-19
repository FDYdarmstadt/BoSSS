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

using System.Collections.Generic;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using BoSSS.Solution.Control;
using System;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Application.XNSERO_Solver {
    public static class ParticleStokesFlow {
        public static XNSERO_Control StokesFlow(int k = 2, int amrLevel = 1) {
            XNSERO_Control C = new XNSERO_Control(degree: k, projectName: "Test");
            //C.SetSaveOptions(@"D:\BoSSS_databases\2particleInteractions", 1);
            //string _DbPath = null;
            C.DbPath = null;
            List<string> boundaryValues = new List<string> {
                "Pressure_Dirichlet_Left",
                "Pressure_Dirichlet_right",
                "Pressure_Dirichlet_upper",
                "Wall_lower",
            };
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.UseImmersedBoundary = true;
            C.dtFixed = 1e-3;
            C.NoOfTimesteps = 50000;
            C.SetBoundaries(boundaryValues);
            C.SetGrid2D(lengthX: 4, lengthY: 4, cellsPerUnitLength: 4, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(0);

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 10;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81);
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.IncludeConvection = true;
            C.SetGravity(new Vector(0, -9.81 ));
            // Particle Properties
            // =============================   
            double particleDensity = 150;
            List<Particle> particles = new List<Particle>();
            Motion motion = new(particleDensity);
            particles.Add(new ParticleEllipse(motion, 0.4, 0.4, new double[] { 1.0, 0.0 }, 0, 0, new double[] { 0, 0 }, 0));
            particles.Add(new ParticleEllipse(motion, 0.4, 0.4, new double[] { -1.0, 0.0 }, 0, 0, new double[] { 0, 0 }, 0));
            C.InitialiseParticles(particles);
            C.CutCellQuadratureType = Foundation.XDG.CutCellQuadratureMethod.Saye;
            double levelSet0(double[] X) => X[0];
            C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), levelSet0);
            C.Option_LevelSetEvolution2 = Solution.LevelSetTools.LevelSetEvolution.RigidObject;
            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.FastMarching;


            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.CutCellQuadratureMethod.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            return C;
        }
    }

    public static class MicroFluidChannel {

        public static XNSERO_Control ParticleSeparation_Steady(int k = 3) {

            XNSERO_Control C = new XNSERO_Control();

            C.SetDGdegree(k);

            // basic database options
            // ======================
            #region db

            string _DbPath = null; // @"D:\local\local_test_db";
            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSERO/MicroChannel";
            C.ProjectDescription = "Micro Channel flow for Particle Separation";

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1000;
            C.PhysicalParameters.mu_A = 1e-3;

            C.PhysicalParameters.IncludeConvection = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double SizeScale = 1.0e-6; 

            double H = 2.0 * SizeScale;
            double H_center = 0.5 * SizeScale;  // defines resolution
            double H_wedge = 0.75 * SizeScale; 
            double L_inflow = 2.0 * SizeScale;
            double L_wedge = 0.75 * SizeScale;
            double L_outflow = 3.0 * SizeScale;
            double W = H / 4.0;

            int resolution = 1;
            int kelemH = 8 * resolution;
            int kelemL = 23 * resolution;
            int kelemW = 2 * resolution;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L_inflow, L_wedge + L_outflow, kelemL + 1);
                double[] Ynodes = GenericBlas.Linspace(-H/2.0, H/2.0, kelemH + 1);
                //var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                double[] Znodes = GenericBlas.Linspace(-W/2.0, W/2.0, kelemW + 1);
                var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                grd.EdgeTagNames.Add(5, "wall_front");
                grd.EdgeTagNames.Add(6, "wall_back");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (H / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (H / 2.0)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] + L_inflow) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - (L_wedge + L_outflow)) <= 1.0e-8)
                        et = 4;
                    if (Math.Abs(X[2] + (W / 2.0)) <= 1.0e-8)
                        et = 5;
                    if (Math.Abs(X[2] - (W / 2.0)) <= 1.0e-8)
                        et = 6;

                    return et;
                });

                return grd;
            };

            //C.SetPeriodicity2D();
            //double[][] bndPos = new double[2][];
            //C.SetBoundaryPosition(bndPos);
            //double[][] wallPos = new double[2][];
            //C.SetWallPosition(wallPos);

            #endregion


            // Particle Properties
            // ===================
            #region particle

            double particle_Density = 1.0;
            double particel_Radius = 0.1 * SizeScale;
            double[] particle_InitPosition = new double[] { -0.2, 0.0 };
            double[] particele_InitTransVelocity = new double[] { 0.0, 0.0 };
            double particele_InitRotlocity = 0.0;

            List<Particle> particles = new List<Particle>();

            MotionFixed noMotion = new (density: particle_Density, spatDim: 3);

            // lower wegde
            MultidimensionalArray pointsL = MultidimensionalArray.Create(3, 2);
            pointsL.SetSubVector(new double[] { 0.0 * SizeScale, -1.0 * SizeScale }, new int[] { 0, -1 });
            pointsL.SetSubVector(new double[] { 0.0 * SizeScale, -0.25 * SizeScale }, new int[] { 1, -1 });
            pointsL.SetSubVector(new double[] { 0.75 * SizeScale, -1.0 * SizeScale }, new int[] { 2, -1 });
            particles.Add(new ParticleTriangle(noMotion, pointsL, true));

            // upper wegde
            MultidimensionalArray pointsU = MultidimensionalArray.Create(3, 2);
            pointsU.SetSubVector(new double[] { 0.75 * SizeScale, 1.0 * SizeScale }, new int[] { 0, -1 });
            pointsU.SetSubVector(new double[] { 0.0 * SizeScale, 0.25 * SizeScale }, new int[] { 1, -1 });
            pointsU.SetSubVector(new double[] { 0.0 * SizeScale, 1.0 * SizeScale }, new int[] { 2, -1 });
            particles.Add(new ParticleTriangle(noMotion, pointsU, true));

            // particle 
            //particles.Add(new ParticleDisk(noMotion, particel_Radius, particle_InitPosition));

            //C.InitialiseParticles(particles);
            C.Particles = particles;

            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.None;
            C.Option_LevelSetEvolution2 = Solution.LevelSetTools.LevelSetEvolution.None;

            //C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new AMRaroundRigidObject(particles, H / kelemH ) { maxRefinementLevel = 1 });
            //C.AMR_startUpSweeps = 1;

            #endregion


            // Initial Values
            // ==============
            #region init


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");

            double U = 0.6;
            //C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => ((-4.0 * U / H.Pow2()) * X[1].Pow2() + U)); // * Math.Sin(2.0*Math.PI*(t/T)));
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => U); // * Math.Sin(2.0*Math.PI*(t/T)));

            C.AddBoundaryValue("pressure_outlet_right");

            C.AddBoundaryValue("wall_front");
            C.AddBoundaryValue("wall_back");

            #endregion


            // misc. solver options
            // ====================
            #region solver


            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            //double dt = 0.138; // 5e-2; // 5e-2;
            //C.dtMax = dt;
            //C.dtMin = dt;
            //C.Endtime = 1000;
            //C.NoOfTimesteps = 100; // 500;
            //C.saveperiod = 10;

            #endregion


            return C;

        }
    }
}
