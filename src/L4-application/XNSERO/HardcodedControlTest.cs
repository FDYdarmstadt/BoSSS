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
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            double levelSet0(double[] X) => X[0];
            C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), levelSet0);
            C.Option_LevelSetEvolution2 = Solution.LevelSetTools.LevelSetEvolution.RigidObject;
            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.FastMarching;


            // Coupling Properties
            // =============================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            return C;
        }
    }
}
