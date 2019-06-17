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
    public class DryParticles : IBM_Solver.HardcodedTestExamples
    {
        /// <summary>
        /// Testing of particle/wall interactions using a single particle
        /// </summary>
        public static FSI_Control HundredDryActiveParticles(string _DbPath = null)
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
            C.AdaptiveMeshRefinement = false;
            C.pureDryCollisions = true;

            // DG degrees
            // ==========

            C.SetDGdegree(2);

            // grid and boundary conditions
            // ============================

            double[] Xnodes = GenericBlas.Linspace(-30, 30, 300);
            double[] Ynodes = GenericBlas.Linspace(-10, 10, 100);
            double h = Math.Min((Xnodes[1] - Xnodes[0]), (Ynodes[1] - Ynodes[0]));

            C.GridFunc = delegate {
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true, periodicY: true);
                grd.EdgeTagNames.Add(1, "Velocity_Inlet");
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 1;
                    return et;
                });

                return grd;
            };

            C.AddBoundaryValue("Velocity_Inlet");

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;

            // Particles
            // =========
            for (int i = 0; i < 29; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    double StartAngle = 10 * i - 10 * i * j;
                    C.Particles.Add(new Particle_Ellipsoid(new double[] { -28 + 2 * i, 8 - 2 * j }, StartAngle)
                    {
                        particleDensity = 1.0,
                        length_P = 0.5,
                        thickness_P = 0.2,
                        ActiveVelocity = 1,
                    });
                }
            }

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
            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 1;

            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;

            double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;

            C.Endtime = 100000.0;
            C.NoOfTimesteps = 50000;

            // haben fertig...
            // ===============

            return C;

        }
    }


}
