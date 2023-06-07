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
using System.IO;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Statistic;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.CahnHilliard.Tests {

    [TestFixture]
    static class TestProgram {

        [Test]
        public static void TestEllipticDroplet([Values(2, 3, 4)] int pDeg = 2) {

            using(var p = new CahnHilliardMain()) {
                var ctrl = Examples.EllipticDroplet(xRes: 20, yRes: 20, pDG: pDeg);
                ctrl.ImmediatePlotPeriod = 1;
                ctrl.SuperSampling = 3;
                p.Init(ctrl);
                p.RunSolverMode();


                // 01feb23: reference values for total concentration and total nuber of iterations, in dependence of degree:
                double TotalConcentration_Ref;
                int TotIter_ref;
                switch (pDeg) {
                    case 2:
                        TotalConcentration_Ref =-715.6359672572; TotIter_ref = 62; break;
                    case 3:
                        TotalConcentration_Ref =-715.63651730895; TotIter_ref = 62; break;
                    case 4:
                        TotalConcentration_Ref =-715.63648317368; TotIter_ref = 62; break;
                    default:
                        throw new ArgumentOutOfRangeException();
                }

                int i = 0;
                foreach (var tt in p.LogHistory) {
                    Console.WriteLine($"c[{i}] = {tt.bq.concentration} @ t = {tt.time}");

                    double concErr = Math.Abs(tt.bq.concentration - TotalConcentration_Ref);
                    Assert.LessOrEqual(concErr, 1.0e-9, "Conservation of concentration is violated");
                    i++;
                }

                p.c.GetExtremalValues(out double c_min, out double c_max);
                Assert.LessOrEqual(c_max, 1.9, "Maximum concentration out-of-range");
                Assert.GreaterOrEqual(c_max, 0.9, "Maximum concentration out-of-range");
                Assert.LessOrEqual(c_min, -0.9, "Minimum concentration out-of-range");
                Assert.GreaterOrEqual(c_min, -1.9, "Minimum concentration out-of-range");

                Assert.IsTrue(p.ResLogger.TimeStep == 11, "Number of timesteps have been modified");
                Assert.GreaterOrEqual(p.ResLogger.LineCounter, 5, "Number of total iterations dropped beyond plausibility");
                Assert.LessOrEqual(p.ResLogger.LineCounter, TotIter_ref, "Number of total iterations (over all timesteps) got out of hand.");

                //Console.WriteLine("No Of Iter = " + p.QueryHandler.QueryResults[QueryHandler.AccNonLinIter]);
            }
        }


        [Test]
        public static void TestEllipticDropletConvergence([Values(2)] int pDeg = 2) {

            int[] Res = new int[] { 10, 20, 40, 80, 160 };

            var Controls = new List<CahnHilliardControl>();
            foreach (int res in Res) {
                var ctrl = Examples.EllipticDroplet(xRes: res, yRes: res, pDG: pDeg);
                Controls.Add(ctrl);
            }


            // 01feb23:
            // Convergence slope for Error of 'c':     2.821947945144496       Intercept:      -1.1137827966318687     (Expecting: 2/-1 to 1 in norm L2_embedded)
            //  Convergence slope for Error of 'mu':   2.8593703784657984      Intercept:      -1.0080261719610017(Expecting: 2/-1 to 1 in norm L2_embedded)


            if (pDeg == 2) {
                ConvergenceTest.SolverConvergenceTest_Experimental(
                     Controls,
                     "CahnHilliard",
                     ("c", NormType.L2_embedded, 2.5, -1.1137827966318687, 0.2),
                     ("mu", NormType.L2_embedded, 2.5, -1.0080261719610017, 0.2)
                     );
            } else {
                throw new NotImplementedException("unknown polynomila degree");
            }
            /*
            using (var p = new CahnHilliardMain()) {
                var ctrl = Examples.EllipticDroplet(xRes: 80, yRes: 80, pDG: 2);
                ctrl.ImmediatePlotPeriod = 1;
                ctrl.SuperSampling = 3;
                p.Init(ctrl);
                p.RunSolverMode();


                //Application<SipControl>._Main(new string[] {
                //        "--control", "cs:ipPoisson.ippHardcodedControl.TestCartesian1()"
                //    },
                //    false,
                //    "",
                //    delegate() {
                //        p = new SipPoissonMain();
                //        return p;
                //    });


                double err = (double)p.QueryHandler.QueryResults["SolL2err"];
                double thres = 5.0e-9;

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }
            */
        }
        

    }
}
