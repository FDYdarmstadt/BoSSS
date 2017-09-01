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

namespace BoSSS.Application.IBM_Solver {

    [TestFixture]
    static class TestProgram {

        [TestFixtureSetUp]
        public static void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application<IBM_Control>.GetBoSSSInstallDir(),
                out dummy);
        }

        [TestFixtureTearDown]
        public static void Cleanup() {
            //Console.Out.Dispose();
            csMPI.Raw.mpiFinalize();
        }

        [Test]
        public static void TestNSECylinder_stationary() {

            IBM_SolverMain p = null;
            Application<IBM_Control>._Main(new string[] {
                "--control","cs:BoSSS.Application.IBM_Solver.HardcodedTestExamples.IBMCylinderFlow(k:2)","--delplt"
            }, false, "", delegate () {
                p = new IBM_SolverMain();
                return p;
            });


            double C_Drag_Sol = 5.58;
            double C_Lift_Sol = 0.0107;
            double C_Drag = (double)p.QueryHandler.QueryResults["C_Drag"];
            double C_Lift = (double)p.QueryHandler.QueryResults["C_Lift"];

            double diff_Drag = Math.Abs(C_Drag - C_Drag_Sol);
            double diff_Lift = Math.Abs(C_Lift - C_Lift_Sol);

            Console.WriteLine(" Drag difference: {0:0.####E+00}, ok ? {1}, less or equal {2:0.####E+00}.", diff_Drag, (bool)(diff_Drag < 0.01), 0.01);
            Console.WriteLine(" Lift difference: {0:0.####E+00}, ok ? {1}, less or equal {2:0.####E+00}.", diff_Lift, (bool)(diff_Lift < 0.001), 0.001);

            Assert.LessOrEqual(diff_Drag, 0.01, "Drag coefficient seems wrong.");
            Assert.LessOrEqual(diff_Lift, 0.001, "Lift coefficient seems wrong.");
        }

        [Test]
        public static void TestChannel() {


            IBM_SolverMain p = null;
            Application<IBM_Control>._Main(new string[] {
                "--control","cs:BoSSS.Application.IBM_Solver.HardcodedTestExamples.ChannelFlow(k:2,pardiso:true)","--delplt"
            }, false, "", delegate () {
                p = new IBM_SolverMain();
                return p;
            });

            double L2VelX = (double)p.QueryHandler.QueryResults["L2err_VelocityX"];

            Assert.LessOrEqual(L2VelX, 1E-8, "L2 Error of VelocityX seems wrong");

        }

        [Test]
        public static void TestChannel_MUMPS() {

            IBM_SolverMain p = null;

            Application<IBM_Control>._Main(new string[] {
                "--control","cs:BoSSS.Application.IBM_Solver.HardcodedTestExamples.ChannelFlow(k:2,pardiso:false)","--delplt"
            }, false, "", delegate () {
                p = new IBM_SolverMain();
                return p;
            });



            double L2VelX = (double)p.QueryHandler.QueryResults["L2err_VelocityX"];

            Assert.LessOrEqual(L2VelX, 1E-8, "L2 Error of VelocityX seems wrong");

        }

        [Test]
        public static void TestChannel_periodic() {

            IBM_SolverMain p = null;
            Application<IBM_Control>._Main(new string[] {
                "--control","cs:BoSSS.Application.IBM_Solver.HardcodedTestExamples.ChannelFlow(k:2,periodic:true,pardiso:true)","--delplt"
            }, false, "", delegate () {
                p = new IBM_SolverMain();
                return p;
            });

            double L2VelX = (double)p.QueryHandler.QueryResults["L2err_VelocityX"];

            Assert.LessOrEqual(L2VelX, 1E-8, "L2 Error of VelocityX seems wrong");
        }
    }

}

