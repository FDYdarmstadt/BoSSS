﻿//Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//*/

using System.Collections.Generic;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class DryParticles : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control periodicTest(int k = 2, int amrLevel = 1, double aspectRatio = 2, double angle = -10) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "2_active_Rods");
            //C.SetSaveOptions(@"/home/ij83requ/default_bosss_db", 1);
            //C.SetSaveOptions(dataBasePath: @"D:\BoSSS_databases\Channel", savePeriod: 1);

            List<string> boundaryValues = new List<string> {
                "Pressure_Dirichlet"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 2, lengthY: 2, cellsPerUnitLength: 4, periodicX: false, periodicY: true);
            //C.SetAddaptiveMeshRefinement(amrLevel: amrLevel);

            // Fluid Properties
            // =============================
            C.CoefficientOfRestitution = 1;
            C.pureDryCollisions = true;

            // Particle Properties
            // =============================
            double particleDensity = 1;
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, true, false, false, 1);
            C.Particles.Add(new Particle_Ellipsoid(motion, 0.2, 0.2 * aspectRatio, new double[] { 0, -2 }, startAngl: 0, activeStress: 0, startTransVelocity: new double[] { 0, 1 }));

            // misc. solver options
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
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
