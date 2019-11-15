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
    public class HC_Dry : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control SingleParticleFalling(int k = 4, int amrLevel = 4) {
            FSI_Control C = new FSI_Control(k, "activeRod_noBackroundFlow", "active Particles");
            //C.SetSaveOptions(dataBasePath: @"/home/ij83requ/default_bosss_db", savePeriod: 1);

            // Domain
            // =============================
            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 0.5, lengthY: 0.5, cellsPerUnitLength: 40, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel);

            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LevelSetSmoothing = false;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.hydrodynamicsConvergenceCriterion = 1e-4;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.IncludeConvection = false;
            double particleDensity = 100;
            C.gravity = new double[] { 0, 0 };
            C.pureDryCollisions = true;

            // Particle Properties
            // =============================   
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, true, false, false);
            double particleRadius = 0.1;
            C.Particles = new List<Particle> {
                new Particle_Ellipsoid(motion, particleRadius, particleRadius, new double[] {0,0 },0, 0, new double[] {0,-0.1 })
            };   

            // misc. solver options
            // =============================  
            C.NonLinearSolver.MaxSolverIterations = 1000;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 1000;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.LevelSetSmoothing = false;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-2, noOfTimesteps: int.MaxValue);

            return C;
        }
    }
}