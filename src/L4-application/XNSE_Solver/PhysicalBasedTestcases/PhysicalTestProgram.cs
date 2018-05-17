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
using System.Threading.Tasks;

using NUnit.Framework;
using MPI.Wrappers;
using BoSSS.Solution;


namespace BoSSS.Application.XNSE_Solver.Tests {

    [TestFixture]
    static public partial class UnitTest {

#if !DEBUG

        /// <summary>
        /// See <see cref="PhysicalBasedTestcases.CapillaryWave.CW_Test"/>.
        /// </summary>
        [Test]
        public static void TestCapillaryWave() {
            var C = BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases.CapillaryWave.CW_Test();
            using (var solver = new XNSE_SolverMain()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        /// <summary>
        /// See <see cref="PhysicalBasedTestcases.RayleighTaylorInstability.RT_Test"/>.
        /// </summary>
        [Test]
        public static void TestRayleighTaylorInstability() {

            var C = BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases.RayleighTaylorInstability.RT_Test();
            using (var solver = new XNSE_SolverMain()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        /*
        [Test]
        public static void TestRisingBubble() {
                        
            var C = BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases.RisingBubble.RB_Test();
            using (var p = new XNSE_SolverMain()) {
                p.Init(C);
                p.RunSolverMode();


                double yCM_end = p.QueryResultTable["yCM"].Select(kvp => (double)kvp.Value).ToArray().Last();
                double[] circ = p.QueryResultTable["circ"].Select(kvp => (double)kvp.Value).ToArray();
                double[] riseV = p.QueryResultTable["riseV"].Select(kvp => (double)kvp.Value).ToArray();
                double[] time = p.QueryResultTable["time"].Select(kvp => (double)kvp.Value).ToArray();

                double circ_min = circ[0];
                double t_circ_min = time[0];
                double riseV_max = riseV[0];
                double t_riseV_max = time[0];

                int Timesteps = time.Length;
                for (int i = 1; i < Timesteps; i++) {
                    if (circ[i] < circ_min) {
                        circ_min = circ[i];
                        t_circ_min = time[i];
                    }
                    if (riseV[i] > riseV_max) {
                        riseV_max = riseV[i];
                        t_riseV_max = time[i];
                    }
                }

                Console.WriteLine("Center of mass y-Position = {0} at t_end = {1}", yCM_end, time[Timesteps - 1]);
                Console.WriteLine("minimum circularity = {0} at t = {1}", circ_min, t_circ_min);
                Console.WriteLine("maximum rise Velocity = {0} at t = {1}", riseV_max, t_riseV_max);

                Assert.LessOrEqual(circ_min, 1.0);
                Assert.GreaterOrEqual(riseV_max, 0.0);
            }
        }
        */
#endif

    }
}
