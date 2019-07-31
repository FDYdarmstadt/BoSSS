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
    public class FourCellsTest : IBM_Solver.HardcodedTestExamples
    {
        public static FSI_Control Main()
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

            C.pureDryCollisions = false;
            C.SetDGdegree(2);

            // grid and boundary conditions
            // ============================

            C.AdaptiveMeshRefinement = true;
            C.RefinementLevel = 4;

            C.GridFunc = delegate
            {
                double[] Xnodes = GenericBlas.Linspace(-1 * BaseSize, 1 * BaseSize, 14);
                double[] Ynodes = GenericBlas.Linspace(-1 * BaseSize, 1 * BaseSize, 14);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Wall_left");
                grd.EdgeTagNames.Add(2, "Wall_right");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet_upper");


                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1 * BaseSize)) <= 1.0e-8)
                        et = 2;

                    if (Math.Abs(X[1] - (-1 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-1 * BaseSize)) <= 1.0e-8)
                        et = 4;


                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };

            C.Particles.Add(new Particle_Sphere(new double[] { 0, 0 }, startAngl: 0)
            {
                particleDensity = 7.8,
                radius_P = 0.1,
                GravityVertical = -9.81,
                useAddaptiveUnderrelaxation = true,
                underrelaxation_factor = 5,
                clearSmallValues = true,
                UseAddedDamping = true
            });

            C.GridPartType = GridPartType.Hilbert;

            C.AddBoundaryValue("Wall_left");
            C.AddBoundaryValue("Wall_right");
            C.AddBoundaryValue("Pressure_Outlet_lower");
            C.AddBoundaryValue("Pressure_Outlet_upper");

            // Initial Values
            // ==============

            // Coupling Properties
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;
            C.CoefficientOfRestitution = 0;

            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            C.PhysicalParameters.IncludeConvection = false;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.ForceAndTorque_ConvergenceCriterion = 1e-2;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 5000;

            // haben fertig...
            // ===============

            return C;
        }
    }
}
