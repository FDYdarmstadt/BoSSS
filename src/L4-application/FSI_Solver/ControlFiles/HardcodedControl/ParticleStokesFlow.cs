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

namespace BoSSS.Application.FSI_Solver {
    public class ParticleStokesFlow : IBM_Solver.HardcodedTestExamples {

        public static FSI_Control StokesFlow(int k = 2, int amrLevel = 1) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\wetParticleCollision", 1);

            List<string> boundaryValues = new List<string> {
                "Wall"
            };
            C.SetBoundaries(boundaryValues);
            C.SetGrid(lengthX: 2, lengthY: 2, cellsPerUnitLength: 12, periodicX: false, periodicY: false);
            C.SetAddaptiveMeshRefinement(amrLevel, true);
            C.hydrodynamicsConvergenceCriterion = 1e-2;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1e-2;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector(0, -0.01 );
            // Particle Properties
            // =============================   
            double particleDensity = 2;
            C.Particles = new List<Particle>();
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, false, false, false, 0);
            C.Particles.Add(new Particle_Sphere(motion, 0.25, new double[] { 0.0, 0.0 },-90, 0, new double[] { 0, 0 }));

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
            C.maxIterationsFullyCoupled = 20000;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: 1e-3, noOfTimesteps: int.MaxValue);

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control WetParticleWallCollision(int k = 3, double DensityFactor = 500, double dtMax = 1e-3, double minDist = 1e-3) {
            FSI_Control C = new FSI_Control(degree: k, projectName: "wetParticleWallCollision");
            C.SetSaveOptions(@"D:\BoSSS_databases\wetParticleCollision", 1);
            //C.SetSaveOptions(@"\\hpccluster\hpccluster-scratch\deussen\cluster_db\WetParticleCollision", 1);
            //C.SetSaveOptions(@"/work/scratch/ij83requ/default_bosss_db", 1);

            List<string> boundaryValues = new List<string> {
                "Wall_left",
                "Wall_right",
                "Wall_lower",
                "Pressure_Outlet_upper"
            };
            C.SetBoundaries(boundaryValues);
            C.SetAddaptiveMeshRefinement(4);
            C.SetGrid(lengthX: 5, lengthY: 1, cellsPerUnitLength: 6, periodicX: false, periodicY: false);
            C.hydrodynamicsConvergenceCriterion = 1e-3;
            C.pureDryCollisions = false;
            C.minDistanceThreshold = minDist;

            // Fluid Properties
            // =============================
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.Material = true;
            C.gravity = new Vector( 0, -5 );
            double particleDensity = 1 * DensityFactor;
            // Particle Properties
            // =============================   
            // Defining particles
            C.Particles = new List<Particle>();
            ParticleMotionInit motion = new ParticleMotionInit(C.gravity, particleDensity, C.pureDryCollisions, false, false, 0);
            C.Particles.Add(new Particle_Sphere(motion, 0.125, new double[] { 0.0, 0 }, 0, 0, new double[] { 0, 0 }));

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
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;


            // Coupling Properties
            // =============================
            C.Timestepper_LevelSetHandling = LevelSetHandling.FSI_LieSplittingFullyCoupled;
            C.LSunderrelax = 1;
            C.maxIterationsFullyCoupled = 2000;


            // Timestepping
            // =============================  
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.SetTimesteps(dt: dtMax, noOfTimesteps: 10000);

            // haben fertig...
            // ===============

            return C;
        }
    }
}
