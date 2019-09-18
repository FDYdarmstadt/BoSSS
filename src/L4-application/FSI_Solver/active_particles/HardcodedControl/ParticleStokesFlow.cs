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
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class ParticleStokesFlow : IBM_Solver.HardcodedTestExamples {

        public static FSI_Control WetParticleWallCollision(int k = 3, double DensityFactor = 1000) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\Channel", 1);

            List<string> boundaryValues = new List<string> {
                "Pressure_Outlet_left",
                "Pressure_Outlet_right",
                "Wall_lower",
                "Pressure_Outlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 10, lengthY: 10, cellsPerUnitLength: 0.5, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 6);
            C.hydrodynamicsConvergenceCriterion = 1e-2;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new double[] { 0, -9.81 };

            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            C.underrelaxationParam = new ParticleUnderrelaxationParam(convergenceLimit: C.hydrodynamicsConvergenceCriterion, underrelaxationFactorIn: 1.0, useAddaptiveUnderrelaxationIn: true);
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, false, false, false, C.underrelaxationParam, 1);
            C.Particles.Add(new Particle_Sphere(motion, 0.5, new double[] { 0.0, 0.0 }, startAngl: 0) {
                particleDensity = 1 * DensityFactor,
            });

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

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
            C.LSunderrelax = 1.0;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 1000000;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-3, noOfTimesteps: 50000);

            // haben fertig...
            // ===============

            return C;
        }
    }
}
