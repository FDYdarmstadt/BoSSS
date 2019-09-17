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
    public class HardcodedControl_multipleActiveParticles : IBM_Solver.HardcodedTestExamples {
        public static FSI_Control ActiveRods_noBackroundFlow(int k = 3) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "9_active_Rods");
            C.SetSaveOptions(@"D:\BoSSS_databases\multipleActiveParticles", 1);

            List<string> boundaryValues = new List<string> {
                "Wall_left",
                "Wall_right",
                "Wall_lower",
                "Wall_upper"
            };

            int sqrtPart = 3;
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 2 * sqrtPart + 1, lengthY: 2 * sqrtPart + 1, cellsPerUnitLength: 1.5, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel: 4);
            C.hydrodynamicsConvergenceCriterion = 0.1;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;

            // Particle Properties
            // =============================
            C.underrelaxationParam = new ParticleUnderrelaxationParam(convergenceLimit: C.hydrodynamicsConvergenceCriterion, underrelaxationFactorIn: 5.0, useAddaptiveUnderrelaxationIn: true);
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, false, false, false, C.underrelaxationParam, 1);
            for (int x = 0; x < sqrtPart; x++) {
                for (int y = 0; y < sqrtPart; y++) {
                    C.Particles.Add(new Particle_Ellipsoid(motion, 0.2, 0.1, new double[] { -sqrtPart + 1 + 2 * x, sqrtPart - 1 - 2 * y }, startAngl: x * y + 26 * y) {
                        particleDensity = 1,
                        activeStress = 1,
                    });
                }
            }
            //

            // Quadrature rules
            // =============================   
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;


            //Initial Values
            // =============================   
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
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000000;
            C.NoOfTimesteps = 1000000;

            // haben fertig...
            // ===============

            return C;
        }
    }
}