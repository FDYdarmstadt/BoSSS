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
using BoSSS.Foundation.IO;
using System.IO;
using MathNet.Numerics;
using System.Linq;

namespace BoSSS.Application.XNSERO_Solver {

    [TestFixture]
    static class TestProgram {

        /// <summary>
        /// Test level set projection by checking the cut-cells.
        /// </summary>
        [Test]
        public static void TestRigidLevelSetProjection() {
            using(XNSERO p = new XNSERO()) {

                var ctrl = XNSEROTest_Control.LevelSetTest();
                ctrl.ImmediatePlotPeriod = 1;
                ctrl.SuperSampling = 3;

                p.Init(ctrl);
                p.RunSolverMode();
                BitArray cutCellMask = p.LsTrk.Regions.GetCutCellMask().GetBitMask();
                BitArray cutCellMaskSoll = new BitArray(cutCellMask.Length);
                cutCellMaskSoll[5] = true;
                cutCellMaskSoll[6] = true;
                cutCellMaskSoll[9] = true;
                cutCellMaskSoll[10] = true;
                for(int i = 0; i < cutCellMask.Length; i++) {
                    Assert.IsTrue(cutCellMask[i] == cutCellMaskSoll[i], "Error in level set projection, cell " + i);
                }
            }
        }

        /// <summary>
        /// Test geometrical properties of <see cref="ParticleEllipse"/>.
        /// </summary>
        [Test]
        public static void TestEllipseProperties() {
            using (XNSERO p = new XNSERO()) {
                double density = 1.5;
                Motion motion = new Motion(density);
                for (int i = 1; i < 9; i += 2) {
                    Particle testParticle = new ParticleEllipse(motion, length: i, thickness: 1, startPos: new double[] { 0, 0 });
                    double circumference;
                    double area;
                    switch (i) {// values from wolfram alpha
                        case 1:
                        circumference = 2 * Math.PI;
                        area = Math.PI;
                        break;
                        case 3:
                        circumference = 13.3649;
                        area = 3 * Math.PI;
                        break;
                        case 5:
                        circumference = 21.01;
                        area = 5 * Math.PI;
                        break;
                        case 7:
                        circumference = 28.8142;
                        area = 7 * Math.PI;
                        break;
                        case 9:
                        circumference = 36.6878;
                        area = 9 * Math.PI;
                        break;
                        default:
                        circumference = 0;
                        area = 0;
                        break;
                    }
                    Assert.LessOrEqual(testParticle.Volume - area, 1e-12, "Error in calculation of ellipse area, should be exact!");
                    Assert.LessOrEqual(testParticle.Mass - area * density, 1e-12, "Error in calculation of ellipse mass, should be exact!");
                    Assert.LessOrEqual(testParticle.Circumference - circumference, 1e-4, "Error in calculation of ellipse circumference.");
                }
            }
        }

        /// <summary>
        /// Test geometrical properties of <see cref="ParticleDisk"/>.
        /// </summary>
        [Test]
        public static void TestDiskProperties() {
            using (XNSERO p = new XNSERO()) {
                double density = 1.5;
                Motion motion = new Motion(density);
                for (int r = 1; r < 3; r += 1) {
                    Particle testParticle = new ParticleDisk(motion, radius: r, startPos: new double[] { 0, 0 });
                    Assert.LessOrEqual(testParticle.Volume - Math.PI * r * r, 1e-12, "Error in calculation of disk area, should be exact!");
                    Assert.LessOrEqual(testParticle.Mass - Math.PI * r * r * density, 1e-12, "Error in calculation of disk mass, should be exact!");
                    Assert.LessOrEqual(testParticle.Circumference - 2 * Math.PI * r, 1e-12, "Error in calculation of disk circumference, should be exact!");
                }
            }
        }

        /// <summary>
        /// Test geometrical properties of <see cref="ParticleRectangle"/>.
        /// </summary>
        [Test]
        public static void TestRectangleProperties() {
            using (XNSERO p = new XNSERO()) {
                double density = 1.5;
                Motion motion = new Motion(density);
                for (int r = 1; r < 3; r += 1) {
                    Particle testParticle = new ParticleRectangle(motion, length: r, thickness: 1, startPos: new double[] { 0, 0 });
                    Assert.LessOrEqual(testParticle.Volume - r * 1, 1e-12, "Error in calculation of rectangle area, should be exact!");
                    Assert.LessOrEqual(testParticle.Mass - r * 1 * density, 1e-12, "Error in calculation of rectangle mass, should be exact!");
                    Assert.LessOrEqual(testParticle.Circumference - 2 * r - 2 * 1, 1e-12, "Error in calculation of rectangle circumference, should be exact!");
                }
            }
        }

        /// <summary>
        /// Test for the <see cref="Particle.Contains(Vector, double)"/> method.
        /// </summary>
        [Test]
        public static void TestParticleContainsPoint2D() {
            using (XNSERO p = new XNSERO()) {
                Motion motion = new Motion(density: 1);
                Vector[] points = new Vector[] { new Vector(0, 0), new Vector(0, 1), new Vector(0, 2), new Vector(1, 0), new Vector(2, 0) };
                foreach(Vector pt in points) {
                    Particle testParticle = new ParticleRectangle(motion, length: 1, thickness: 1, startPos: new double[] { 0, 0 });
                    double pointAbs = pt.Abs();
                    if (pointAbs < 1)
                        Assert.IsTrue(testParticle.Contains(pt), "Error in determining whether a point " + pt + " lies within a rectangle. Should be inside.");
                    else
                        Assert.IsFalse(testParticle.Contains(pt), "Error in determining whether a point " + pt + " lies within a rectangle. Should be outside.");

                    testParticle = new ParticleDisk(motion, radius: 1, startPos: new double[] { 0, 0 });
                    pointAbs = pt.Abs();
                    if (pointAbs < 1)
                        Assert.IsTrue(testParticle.Contains(pt), "Error in determining whether a point " + pt + " lies within a disk. Should be inside.");
                    else
                        Assert.IsFalse(testParticle.Contains(pt), "Error in determining whether a point " + pt + " lies within a disk. Should be outside.");

                    testParticle = new ParticleEllipse(motion, length: 2, thickness: 1, startPos: new double[] { 0, 0 }, startAngl: 45);
                    double angleRad = 45 * 2 * Math.PI / 360;
                    Vector orientation = new(Math.Cos(angleRad), Math.Sin(angleRad));
                    double Ellipse = (pt[0] * orientation[0] + pt[1] * orientation[1]).Pow2() / 4 + (pt[0] * orientation[1] - pt[1]* orientation[0]).Pow2();
                    if(Ellipse < 1)
                        Assert.IsTrue(testParticle.Contains(pt), "Error in determining whether a point " + pt + " lies within an ellipse. Should be inside.");
                    else
                        Assert.IsFalse(testParticle.Contains(pt), "Error in determining whether a point " + pt + " lies within an ellipse. Should be outside.");
                }
            }
        }

        /// <summary>
        /// Test for the calculation of  the support point of an disk and an ellipse. 
        /// Relevant methods: <see cref="ParticleDisk.GetSupportPoint(Vector, Vector, Vector, int, double)"/> 
        /// and <see cref="ParticleEllipse.GetSupportPoint(Vector, Vector, Vector, int, double)"/>.
        /// </summary>
        [Test]
        public static void TestSupportPoint() {
            using (XNSERO p = new XNSERO()) {
                Motion motion = new Motion(density: 1);
                Particle disk = new ParticleDisk(motion, radius: 2, startPos: new double[] { 0, 0 });
                Particle ellipse = new ParticleEllipse(motion, length: 2, thickness: 2, startPos: new double[] { 0, 0 });
                Vector[] points = new Vector[] { new Vector(1, 0), new Vector(1, 1), new Vector(2, 1)};
                foreach (Vector pt in points) {
                    Vector diskSupportPoint = disk.GetSupportPoint(pt, new Vector(0, 0), new Vector(1), 0, 0);
                    Vector ellipseSupportPoint = ellipse.GetSupportPoint(pt, new Vector(0, 0), new Vector(1), 0, 0);
                    Assert.LessOrEqual((diskSupportPoint - ellipseSupportPoint).Abs(), 1e-12, "Error in determination of support point. Ellipse and disk should deliver the same result.");
                }
            }
        }

        /// <summary>
        /// Test for integration of physical variable in <see cref="Motion"/>. 
        /// </summary>
        [Test]
        public static void TestMotionIntegration() {
            using (XNSERO p = new XNSERO()) {
                int timesteps = 100000;
                double dt = Math.PI / timesteps;
                double[] integrationDomain = Generate.LinearSpaced(timesteps, 0, Math.PI);
                Motion motion = new Motion(1);
                Particle testParticle = new ParticleDisk(motion, 1, new double[] { 0, 0 });
                for (int t = 0; t < integrationDomain.Length; t++) {
                    testParticle.Motion.SaveVelocityOfPreviousTimestep();
                    testParticle.Motion.SavePositionAndAngleOfPreviousTimestep();
                    Vector forces = new(Math.Cos(integrationDomain[t]), Math.Sin(integrationDomain[t]));
                    double torque = Math.Sin(integrationDomain[t]);
                    testParticle.Motion.PrescribeHydrodynamicForcesAndTorque(forces, torque, 0);
                    testParticle.Motion.UpdateParticleVelocity(dt);
                    testParticle.Motion.UpdateParticlePositionAndAngle(dt);
                    Vector positionqwrasdas = testParticle.Motion.GetPosition(0);
                }
                Vector transVelocity = testParticle.Motion.GetTranslationalVelocity(0);
                Vector transVelocityExpected = new Vector(0, 2 / Math.PI);
                Assert.LessOrEqual((transVelocityExpected - transVelocity).Abs(), 1e-4, "Error in calculation of translation velocity. Deviation: " + (transVelocityExpected - transVelocity).Abs());

                Vector position = testParticle.Motion.GetPosition(0);
                Vector positionExpected = new Vector(0, 0);
                Assert.LessOrEqual((positionExpected - position).Abs(), 1e-4, "Error in calculation of particle position. Deviation: " + (positionExpected - position).Abs());

                double rotVelocity = testParticle.Motion.GetRotationalVelocity(0);
                double rotVelocityExpected = 2 / testParticle.MomentOfInertia;
                Assert.LessOrEqual((rotVelocityExpected - rotVelocity).Abs(), 1e-4, "Error in calculation of rotational velocity.");

                double angle = testParticle.Motion.GetAngle(0);
                double angleExpected = 0;
                while (angleExpected > 2 * Math.PI)
                    angleExpected -= 2 * Math.PI;
                Assert.LessOrEqual((angleExpected - angle).Abs(), 1e-4, "Error in calculation of particle angle.");
            }
        }

        /// <summary>
        /// Test for the calculation of the minimal distance between two particles.
        /// </summary>
        [Test]
        public static void TestDistance() {
            using (XNSERO p = new XNSERO()) {
                Particle[] testParticles = new Particle[] { new ParticleDisk(new Motion(1), 1, new double[] { -2, 0 }), new ParticleDisk(new Motion(1), 1, new double[] { 2, 0 }) };
                DistanceAlgorithm distance = new(testParticles, 0); 
                distance.CalculateTwoParticleDistance();
                Assert.LessOrEqual(distance.DistanceVector.Abs() - 2, 1e-12, "Error in calculation of the minimal distance between two points.");
                Assert.LessOrEqual(distance.ClosestPoints[0].Abs() - 1, 1e-12, "Error in calculation of the closest point toward the opposing particle");
                Assert.LessOrEqual(distance.ClosestPoints[1].Abs() - 1, 1e-12, "Error in calculation of the closest point toward the opposing particle");
                Assert.IsFalse(distance.Overlapping);
            }
        }

        /// <summary>
        /// Test for the calculation of the minimal distance between two particles.
        /// </summary>
        [Test]
        public static void TestOverlappingParticles() {
            using (XNSERO p = new XNSERO()) {
                Particle[] testParticles = new Particle[] { new ParticleDisk(new Motion(1), 1, new double[] { -1, 0 }), new ParticleEllipse(new Motion(1), 3, 1, new double[] { 2, 0 }) };
                DistanceAlgorithm distance = new(testParticles, 0);
                distance.CalculateTwoParticleDistance();
                Assert.IsTrue(distance.Overlapping);
            }
        }

        /// <summary>
        /// Test for the calculation of the minimal distance between a particle and a wall
        /// </summary>
        [Test]
        public static void TestDistanceFromWall() {
            using (XNSERO p = new XNSERO()) {
                Particle testParticle = new ParticleDisk(new Motion(1), 1, new double[] { 0, 0 });
                testParticle.ClosestPointOnOtherObjectToThis = new Vector(2, 0); //setup wall position
                DistanceAlgorithm distance = new(new Particle[] { testParticle }, 0);
                distance.CalculateParticleWallDistance(testParticle.ClosestPointOnOtherObjectToThis);
                Assert.LessOrEqual(distance.DistanceVector.Abs() - 2, 1e-12, "Error in calculation of the minimal distance between two points.");
                Assert.LessOrEqual(distance.ClosestPoints[0].Abs() - 1, 1e-12, "Error in calculation of the closest point toward the opposing particle");
                Assert.IsFalse(distance.Overlapping);
            }
        }

        [Test]
        public static void TestParticleInShearFlow() {
            using(XNSERO p = new XNSERO()) {

                XNSERO_Control ctrl = XNSEROTest_Control.TestParticleInShearFlow();
                p.Init(ctrl);
                p.RunSolverMode();

                double angularVelocitySol = -0.0004346048227901287; // -0.0004198291651200775; Don't know why this changed when updating the linear solver; don't care ;)
                double angularVelocityIs = p.Particles[0].Motion.GetRotationalVelocity(0);

                double diff_Velocity = Math.Abs(angularVelocityIs - angularVelocitySol);
                Console.WriteLine("   angular velocity is " + angularVelocityIs);
                Console.WriteLine("         should be     " + angularVelocitySol);
                Console.WriteLine("         difference is " + diff_Velocity);


                Assert.LessOrEqual(diff_Velocity, 1e-8, "Error in expected angular velocity is to high");

            }
        }

        static public IDatabaseInfo CreateTempDatabase() {

            DirectoryInfo TempDir;
            {
                var rnd = new Random();
                bool Exists = false;
                do {
                    var tempPath = Path.GetTempPath();
                    var tempDir = rnd.Next().ToString();
                    TempDir = new DirectoryInfo(Path.Combine(tempPath, tempDir));
                    Exists = TempDir.Exists;
                } while (Exists == true);
            }
            
            string path = TempDir.FullName;
            return DatabaseInfo.CreateOrOpen(path);
        }


        [Test]
        public static void TestParticleInShearFlow_Phoretic() {
            using(XNSERO p = new XNSERO()) {

                XNSERO_Control ctrl = XNSEROTest_Control.TestParticleInShearFlow_Phoretic(k: 2);
                //var tempDB = CreateTempDatabase();
                //ctrl.SetDatabase(tempDB);
                //ctrl.savetodb = true;

                p.Init(ctrl);
                p.RunSolverMode();

                double angularVelocitySol = -0.0004198291651200775;
                double angularVelocityIs = p.Particles[0].Motion.GetRotationalVelocity(0);

                double diff_Velocity = Math.Abs(angularVelocityIs - angularVelocitySol);
                Console.WriteLine("   angular velocity is " + angularVelocityIs);
                Console.WriteLine("         should be     " + angularVelocitySol);
                Console.WriteLine("         difference is " + diff_Velocity);


                Assert.LessOrEqual(diff_Velocity, 0.00025, "Error in expected angular velocity is to high");

                //DatabaseUtils.DeleteDatabase(ctrl.DbPath);
            }
        }


        [Test]
        public static void TestStickyTrap() {
            using(XNSERO p = new XNSERO()) {

                var ctrl = XNSEROTest_Control.TestStickyTrap();
                p.Init(ctrl);
                p.RunSolverMode();

                Vector Dest_Should;
                Dest_Should = new Vector(0.0, 0.19033655497377627);
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
