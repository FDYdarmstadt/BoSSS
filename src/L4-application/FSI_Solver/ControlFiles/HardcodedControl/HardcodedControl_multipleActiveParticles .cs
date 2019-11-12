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
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedControl_multipleActiveParticles : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control ActiveRods_noBackroundFlow(int k = 2) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "7_active_Rods");
            C.SetSaveOptions(@"/home/ij83requ/default_bosss_db", 1);

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 10, lengthY: 10, cellsPerUnitLength: 3, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 4);
            C.hydrodynamicsConvergenceCriterion = 1e-2;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;

            // Particle Properties
            // =============================
            double particleDensity = 1000;
            C.underrelaxationParam = new ParticleUnderrelaxationParam(C.hydrodynamicsConvergenceCriterion, ParticleUnderrelaxationParam.UnderrelaxationMethod.Jacobian, relaxationFactor: 3.0, useAddaptiveUnderrelaxation: true);
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, C.underrelaxationParam, 1.5);
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -3, 3 }, startAngl: -32, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -2.8, 0 }, startAngl: 49, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0.2, -3.1 }, startAngl: -98, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 2.8, -0.5 }, startAngl: 182, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 0.6, 1.5 }, startAngl: 99, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { -3.0, -3.0 }, startAngl: 56, activeStress: 10));
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.4, 0.2, new double[] { 1.0, 3.0 }, startAngl: 180, activeStress: 10));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000000;
            C.NoOfTimesteps = int.MaxValue;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LSunderrelax = 1.0;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

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