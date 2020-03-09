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
        public static void IterativeSolverTest([Values(Code.exp_Kcycle_schwarz, Code.exp_gmres_levelpmg)] Code SolverName) {
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

        private static double LogLogRegression(IEnumerable<double> _xValues, IEnumerable<double> _yValues) {
            double[] xValues = _xValues.Select(x => Math.Log10(x)).ToArray();
            double[] yValues = _yValues.Select(y => Math.Log10(y)).ToArray();

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
        [Test]
        public static void DiscretizationScalingTest(
#if DEBUG
            [Values(1)] 
#else
            [Values(1,2,3,4)] 
#endif
            int dgDegree
            ) //
        {

            var Controls = new List<XdgPoisson3Control>();

            int[] ResS = null;
            switch(dgDegree) {
                case 1: ResS = new int[] { 8, 9, 16, 17, 32, 33, 64, 65 }; break;
                case 2: ResS = new int[] { 8, 9, 16, 17, 32, 33, 64, 65 }; break;
                case 3: ResS = new int[] { 8, 9, 16, 17, 32, 33, 64 }; break;
                case 4: ResS = new int[] { 8, 9, 16, 17, 32, 33 }; break;
                default: throw new NotImplementedException();
            }
            
            foreach (int res in ResS) {
                var C = HardCodedControl.Circle(Resolution: res, p: dgDegree);
                C.LinearSolver.SolverCode = Code.classic_pardiso;
                C.savetodb = false;

                Controls.Add(C);
            }
            
            
            var data = RunAndLog(Controls);

            //string[] xKeys = data.Keys.Where(name => name.StartsWith("Grid:")).ToArray();
            //string[] yKeys = data.Keys.Where(name => !name.StartsWith("Grid:")).ToArray();


            /*
            Für p = 1:
            Slope for TotCondNo-Var0: 2.153e00
            Slope for InnerCondNo-Var0: 2.59e00
            Slope for StencilCondNo-innerUncut-Var0: 1.486e-01
            Slope for StencilCondNo-innerCut-Var0: 9.957e-02
            Slope for StencilCondNo-bndyUncut-Var0: 1.226e-04
            Slope for StencilCondNo-bndyCut-Var0: 0e00             
            */

            var ExpectedSlopes = new List<ValueTuple<string, string, double>>();
            ExpectedSlopes.Add(("Grid:1Dres", "TotCondNo-Var0", 2.5));
            ExpectedSlopes.Add(("Grid:1Dres", "StencilCondNo-innerUncut-Var0", 0.5));
            ExpectedSlopes.Add(("Grid:1Dres", "StencilCondNo-innerCut-Var0", 0.5));


            foreach (var ttt in ExpectedSlopes) {
                double[] xVals = data[ttt.Item1];
                double[] yVals = data[ttt.Item2];

                double Slope = LogLogRegression(xVals, yVals);

                Console.WriteLine($"Slope for {ttt.Item2}: {Slope:0.###e-00}");
            }

            foreach (var ttt in ExpectedSlopes) {
                double[] xVals = data[ttt.Item1];
                double[] yVals = data[ttt.Item2];

                double Slope = LogLogRegression(xVals, yVals);

                Assert.LessOrEqual(Slope, ttt.Item3, $"Condition number slope for {ttt.Item2} to high; at max. {ttt.Item3}");
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
