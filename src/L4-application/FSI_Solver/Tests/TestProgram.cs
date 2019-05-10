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
using System.Linq;
using System.Text;
using System.IO;
using BoSSS.Foundation;
using BoSSS.Solution;
using MPI.Wrappers;
using NUnit.Framework;
using BoSSS.Platform.LinAlg;
using ilPSP;

namespace BoSSS.Application.FSI_Solver {

    [TestFixture]
    static class TestProgram {

        [TestFixtureSetUp]
        public static void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application<FSI_Control>.GetBoSSSInstallDir(),
                out dummy);
        }

        [TestFixtureTearDown]
        public static void Cleanup() {
            //Console.Out.Dispose();
            csMPI.Raw.mpiFinalize();
        }

        [Test]
        public static void TestFlowRotationalCoupling() {
            using (FSI_SolverMain p = new FSI_SolverMain()) {

                var ctrl = BoSSS.Application.FSI_Solver.HardcodedTestExamples.ParticleInShearFlow(k: 1);
                //ctrl.ImmediatePlotPeriod = 1;
                //ctrl.SuperSampling = 2;
                p.Init(ctrl);
                p.RunSolverMode();

                double angularVelocity_Sol = 0.00487;

                double angularVelocity = (double)p.QueryHandler.QueryResults["Angular_Velocity"];

                double diff_Velocity = Math.Abs(angularVelocity - angularVelocity_Sol);
                Console.WriteLine("   angular velocity is " + angularVelocity);
                Console.WriteLine("         should be     " + angularVelocity_Sol);
                Console.WriteLine("         difference is " + diff_Velocity);


                Assert.LessOrEqual(diff_Velocity, 0.00025, "Error in expected angular velocity is to high");

            }
        }

        [Test]
        public static void SingleDryParticleAgainstWall([Values(false, true)]  bool MeshRefine) { 
            using (FSI_SolverMain p = new FSI_SolverMain()) {

                var ctrl = BoSSS.Application.FSI_Solver.HardcodedTestExamples.SingleDryParticleAgainstWall(MeshRefine:MeshRefine);
                p.Init(ctrl);
                p.RunSolverMode();


                Vector Dest_Should;
                if (MeshRefine)
                    Dest_Should = new Vector(0.0695097474610063, -0.594908028831844); 
                else
                    Dest_Should = new Vector(-0.238381401305482, 0.341721015345088); 

                Vector Dest_Is = new Vector(p.Particles[0].Position[0]);

                double dist = (Dest_Should - Dest_Is).L2Norm();

                Console.WriteLine("Particle reached position " + Dest_Is + ", expected at " + Dest_Should + ", distance is " + dist);

                Assert.Less(dist, 0.1, "Particle to far from expected position");
            }
        }

        [Test]
        public static void DryParticleBounce()
        {
            using (FSI_SolverMain p = new FSI_SolverMain())
            {

                var ctrl = BoSSS.Application.FSI_Solver.HardcodedTestExamples.DryParticleBounce();
                p.Init(ctrl);
                p.RunSolverMode();

                Vector Dest_Should;
                Dest_Should = new Vector(0.263026905796573, 0.788180688520332);

                Vector Dest_Is = new Vector(p.Particles[0].Position[0]);

                double dist = (Dest_Should - Dest_Is).L2Norm();

                Console.WriteLine("Particle reached position " + Dest_Is + ", expected at " + Dest_Should + ", distance is " + dist);

                Assert.Less(dist, 0.1, "Particle to far from expected position");
            }
        }

        [Test]
        public static void StickyTrap()
        {
            using (FSI_SolverMain p = new FSI_SolverMain())
            {

                var ctrl = BoSSS.Application.FSI_Solver.HardcodedTestExamples.StickyTrap();
                p.Init(ctrl);
                p.RunSolverMode();

                Vector Dest_Should;
                Dest_Should = new Vector(0.0, 0.115138000000001);
                double VelY_Should = 0;

                Vector Dest_Is = new Vector(p.Particles[0].Position[0]);
                double VelY_Is = p.Particles[0].TranslationalVelocity[0][0];

                double dist = (Dest_Should - Dest_Is).L2Norm();
                double Vel_Div = Math.Abs(VelY_Should - VelY_Is);

                Console.WriteLine("Particle reached position " + Dest_Is + ", expected at " + Dest_Should + ", distance is " + dist);
                Console.WriteLine("Particle end velocitiy " + VelY_Is + ", expected velocity: " + VelY_Should + ", difference is " + Vel_Div);

                Assert.Less(dist, 0.1, "Particle to far from expected position.");
                Assert.Less(Vel_Div, 0.01, "Particle is moving.");
            }
        }

        [Test]
        public static void ActiveParticle_ForceTest()
        {
            using (FSI_SolverMain p = new FSI_SolverMain())
            {
                var ctrl = HardcodedTestExamples.ActiveParticle_ForceTest();
                p.Init(ctrl);
                p.RunSolverMode();

                double ForcesSoll = 15376.4338960998;

                double Forces = p.Particles[0].HydrodynamicForces[0][0];

                double DiffForces = Math.Abs(ForcesSoll - Forces); 

                Assert.LessOrEqual(DiffForces, 20);
            }
        }
    }
}
