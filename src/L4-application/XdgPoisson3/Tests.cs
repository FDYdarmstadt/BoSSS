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

using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using Code = BoSSS.Solution.Control.LinearSolverCode;

namespace BoSSS.Application.XdgPoisson3 {
    [TestFixture]
    static public class Tests {


        /// <summary>
        /// MPI finalize.
        /// </summary>
        [TestFixtureTearDown]
        static public void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// MPI init.
        /// </summary>
        [TestFixtureSetUp]
        static public void TestFixtureSetUp() {
            BoSSS.Solution.Application.InitMPI(new string[0]);
            XQuadFactoryHelper.CheckQuadRules = true;
        }

        /// <summary>
        /// Test for iterative solver
        /// </summary>
        /// <param name="SolverName"></param>
        [Test]
        public static void SolverTest([Values(Code.exp_Kcycle_schwarz, Code.exp_gmres_levelpmg)] Code SolverName) {
            using (var solver = new XdgPoisson3Main()) {

                int Res, p;
#if DEBUG
                Res = 6;
                p = 2;
#else
                Res = 12;
                p = 3;
#endif
                var C = HardCodedControl.Ball3D(pDeg: p,Res: Res, solverCode: SolverName);

                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        private static double Regression(double[] xValues, double[] yValues) {
            double xAvg = xValues.Average();
            double yAvg = yValues.Average();

            double v1 = 0.0;
            double v2 = 0.0;

            for (int i = 0; i < yValues.Length; i++) {
                v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
                v2 += Math.Pow(xValues[i] - xAvg, 2);
            }

            double a = v1 / v2;
            double b = yAvg - a * xAvg;

            return a;
        }


        /// <summary>
        /// Grid scale tests for condition numbers
        /// </summary>
        public static void DiscretizationScalingTest(
            [Values(1,2)] int dgDegree
            ) {

            var Controls = new List<XdgPoisson3Control>();
            foreach (int res in new int[] { 8, 9, 16, 17, 32, 33, 64, 65 }) {
                var C = HardCodedControl.Circle(Resolution: res, p: dgDegree);
                C.LinearSolver.SolverCode = Code.classic_pardiso;
                C.savetodb = false;

                Controls.Add(C);
            }
            
            
            var data = RunAndLog(Controls);

            string[] xKeys = data.Keys.Where(name => name.StartsWith("Grid:")).ToArray();
            string[] yKeys = data.Keys.Where(name => !name.StartsWith("Grid:")).ToArray();

            foreach(var yKey in yKeys) {
                double[] xVals = data["Grid:1Dres"];
                double[] yVals = data[yKey];

                double Slope = Regression(xVals.Select(x => Math.Log10(x)).ToArray(), yVals.Select(y => Math.Log10(y)).ToArray());

                Console.WriteLine($"Slope for {yKey}: {Slope:0.###e-00}");
            }
        }

        private static IDictionary<string, double[]> RunAndLog(List<XdgPoisson3Control> Controls) {
            var ret = new Dictionary<string, List<double>>();

            foreach (var C in Controls) {
                using (var solver = new XdgPoisson3Main()) {
                    solver.Init(C);
                    solver.RunSolverMode();

                    int J = Convert.ToInt32(solver.CurrentSessionInfo.KeysAndQueries["Grid:NoOfCells"]);
                    double hMin = Convert.ToDouble(solver.CurrentSessionInfo.KeysAndQueries["Grid:hMin"]);
                    double hMax = Convert.ToDouble(solver.CurrentSessionInfo.KeysAndQueries["Grid:hMax"]);
                    int D = Convert.ToInt32(solver.CurrentSessionInfo.KeysAndQueries["Grid:SpatialDimension"]);
                    double J1d = Math.Pow(J, 1.0 / D);

                    var prop = solver.OperatorAnalysis();

                    if (ret.Count == 0) {
                        ret.Add("Grid:NoOfCells", new List<double>());
                        ret.Add("Grid:hMin", new List<double>());
                        ret.Add("Grid:hMax", new List<double>());
                        ret.Add("Grid:1Dres", new List<double>());

                        foreach (var kv in prop) {
                            ret.Add(kv.Key, new List<double>());
                        }
                    }

                    {
                        ret["Grid:NoOfCells"].Add(J);
                        ret["Grid:hMin"].Add(hMin);
                        ret["Grid:hMax"].Add(hMax);
                        ret["Grid:1Dres"].Add(J1d);

                        foreach (var kv in prop) {
                            ret[kv.Key].Add(kv.Value);
                        }
                    }
                }
            }

            // data conversion & return 
            // ========================
            {
                var realRet = new Dictionary<string, double[]>();
                foreach (var kv in ret) {
                    realRet.Add(kv.Key, kv.Value.ToArray());
                }

                return realRet;
            }
        }
    }
}
