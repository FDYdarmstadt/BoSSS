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
using MPI.Wrappers;
using NUnit.Framework;
using ilPSP;
using System.Collections;

namespace BoSSS.Application.XNSERO_Solver {

    [TestFixture]
    static class TestProgram {

        [Test]
        public static void TestRigidLevelSetProjection() {
            using(XNSERO p = new XNSERO()) {

                var ctrl = XNSEROTest_Control.LevelSetTest();
                p.Init(ctrl);
                p.RunSolverMode();
                BitArray cutCellMask = p.LsTrk.Regions.GetCutCellMask().GetBitMask();
                BitArray cutCellMaskSoll = new BitArray(cutCellMask.Length);
                cutCellMaskSoll[5] = true;
                cutCellMaskSoll[6] = true;
                cutCellMaskSoll[8] = true;
                cutCellMaskSoll[10] = true;
                for(int i = 0; i < cutCellMask.Length; i++) {
                    Assert.IsTrue(cutCellMask[i] == cutCellMaskSoll[i], "Error in level set projection, cell " + i);
                }
            }
        }

        [Test]
        public static void TestParticleParameter() {
            using(XNSERO p = new XNSERO()) {

                var ctrl = XNSEROTest_Control.TestParticleParameter();
                p.Init(ctrl);
                p.RunSolverMode();

                double p0_area = p.Particles[0].Area;
                double p0_area_soll = Math.PI;
                double p0_Mass = p.Particles[0].Motion.Mass;
                double p1_area = p.Particles[1].Area;


                double diff_Area1 = Math.Abs(p0_area - p0_area_soll);
                double diff_Area2 = Math.Abs(p0_area - p1_area);
                double diff_Mass = Math.Abs(2 * p0_area - p0_Mass);
                Assert.LessOrEqual(diff_Area1, 1e-12, "Error in calculation of particle area");
                Assert.LessOrEqual(diff_Area2, 1e-12, "Error in calculation of particle area");
                Assert.LessOrEqual(diff_Mass, 1e-12, "Error in calculation of particle mass");

            }
        }

        [Test]
        public static void TestParticleInShearFlow() {
            using(XNSERO p = new XNSERO()) {

                var ctrl = XNSEROTest_Control.TestParticleInShearFlow(k: 2);
                p.Init(ctrl);
                p.RunSolverMode();

                double angularVelocitySol = -0.00283955477369256;
                double angularVelocityIs = p.Particles[0].Motion.GetRotationalVelocity(0);

                double diff_Velocity = Math.Abs(angularVelocityIs - angularVelocitySol);
                Console.WriteLine("   angular velocity is " + angularVelocityIs);
                Console.WriteLine("         should be     " + angularVelocitySol);
                Console.WriteLine("         difference is " + diff_Velocity);


                Assert.LessOrEqual(diff_Velocity, 0.00025, "Error in expected angular velocity is to high");

            }
        }

        [Test]
        public static void TestStickyTrap() {
            using(XNSERO p = new XNSERO()) {

                var ctrl = XNSEROTest_Control.TestStickyTrap();
                p.Init(ctrl);
                p.RunSolverMode();

                Vector Dest_Should;
                Dest_Should = new Vector(0.0, 0.0792658735354393);
                double VelY_Should = 0;

                Vector Dest_Is = new Vector((double[])p.Particles[0].Motion.GetPosition(0));
                double VelY_Is = p.Particles[0].Motion.GetTranslationalVelocity(0)[0];

                double dist = (Dest_Should - Dest_Is).L2Norm();
                double Vel_Div = Math.Abs(VelY_Should - VelY_Is);

                Console.WriteLine("Particle reached position " + Dest_Is + ", expected at " + Dest_Should + ", distance is " + dist);
                Console.WriteLine("Particle end velocitiy " + VelY_Is + ", expected velocity: " + VelY_Should + ", difference is " + Vel_Div);

                Assert.Less(dist, 0.1, "Particle to far from expected position.");
                Assert.Less(Vel_Div, 0.05, "Particle is moving.");
            }
        }
    }
}
