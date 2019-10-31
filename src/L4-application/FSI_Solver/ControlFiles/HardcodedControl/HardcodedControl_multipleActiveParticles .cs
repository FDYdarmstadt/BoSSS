﻿/* =======================================================================
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
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedControl_multipleActiveParticles : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control ActiveRods_noBackroundFlow(int k = 2) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "9_active_Rods");
            C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\25_particles", 1);

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            int sqrtPart = 2;
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: sqrtPart + 1, lengthY: sqrtPart + 1, cellsPerUnitLength: 4, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 1);
            C.hydrodynamicsConvergenceCriterion = 1e-3;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.01;
            C.PhysicalParameters.IncludeConvection = false;

            // Particle Properties
            // =============================
            double particleDensity = 2;
            C.underrelaxationParam = new ParticleUnderrelaxationParam(C.hydrodynamicsConvergenceCriterion, ParticleUnderrelaxationParam.UnderrelaxationMethod.ProcentualRelaxation, relaxationFactor: 3.0, useAddaptiveUnderrelaxation: true);
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, C.underrelaxationParam, 1);
            for (int x = 0; x < sqrtPart; x++) {
                for (int y = 0; y < sqrtPart; y++) {
                    C.Particles.Add(new Particle_Ellipsoid(motion, 0.25, 0.1, new double[] { -0.5 + 1 * x, 0.5 - 1 * y }, startAngl: 180 -30 - 90 * (x - y) + 180 * (1 - x * y), activeStress: 50));
                }
            }

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000000;
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
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 1000000;

            return C;
        }
    }
}