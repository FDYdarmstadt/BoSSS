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
                p.Init(ctrl);
                p.RunSolverMode();

                double angularVelocity_Sol = 0.00487;

                double angularVelocity = (double)p.QueryHandler.QueryResults["Angular_Velocity"];

                double diff_Velocity = Math.Abs(angularVelocity - angularVelocity_Sol);

                Assert.LessOrEqual(diff_Velocity, 0.00025);

            }
        }

        /// <summary>
        /// Note: this test is fucked; the results are nowhere near where you would expext.
        /// </summary>
        [Test]
        public static void SingleDryParticleAgainstWall([Values(false, true)]  bool MeshRefine) {
            using (FSI_SolverMain p = new FSI_SolverMain()) {

                var ctrl = BoSSS.Application.FSI_Solver.HardcodedTestExamples.SingleDryParticleAgainstWall(MeshRefine:MeshRefine);
                p.Init(ctrl);
                p.RunSolverMode();


                Vector Dest_Should;
                if (MeshRefine)
                    Dest_Should = new Vector(0.089255650988794, -1.08925565098877); //new Vector(0.420719299693095, -0.907165088781989);
                else
                    Dest_Should = new Vector(1.80535999455424, -0.785548829055413); //new Vector(0.748512025578859, -0.578342794422653);

                Vector Dest_Is = new Vector(p.Particles[0].Position[0]);

                double dist = (Dest_Should - Dest_Is).L2Norm();

                Console.WriteLine("Particle reached position " + Dest_Is + ", expected at " + Dest_Should + ", distance is " + dist);

                Assert.Less(dist, 0.1, "Particle to far from expected position");
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

                double ForcesSoll = 1481.4254921133;

                double Forces = p.Particles[0].HydrodynamicForces[0][0];

                double DiffForces = Math.Abs(ForcesSoll - Forces); 

                Assert.LessOrEqual(DiffForces, 20);
            }
        }
    }
}
