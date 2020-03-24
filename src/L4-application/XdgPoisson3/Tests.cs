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
using BoSSS.Solution.AdvancedSolvers.Testing;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using Code = BoSSS.Solution.Control.LinearSolverCode;

namespace BoSSS.Application.XdgPoisson3 {
    
    /// <summary>
    /// NUnit test routines
    /// </summary>
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
                case 1: ResS = new int[] { 8, 9, 16, 17, 32, 33, 64, 65, 128 }; break;
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

            /*
            Für p = 1:
            Slope for TotCondNo-Var0: 2.153e00
            Slope for InnerCondNo-Var0: 2.59e00
            Slope for StencilCondNo-innerUncut-Var0: 1.486e-01
            Slope for StencilCondNo-innerCut-Var0: 9.957e-02
            Slope for StencilCondNo-bndyUncut-Var0: 1.226e-04
            Slope for StencilCondNo-bndyCut-Var0: 0e00             
            */

            //var ExpectedSlopes = new List<ValueTuple<Solution.AdvancedSolvers.Testing.OpAnalysisBase.XAxisDesignation, string, double>>();

            //ExpectedSlopes.Add((Solution.AdvancedSolvers.Testing.OpAnalysisBase.XAxisDesignation.Grid_1Dres, "TotCondNo-Var0", 2.5));
            //ExpectedSlopes.Add((Solution.AdvancedSolvers.Testing.OpAnalysisBase.XAxisDesignation.Grid_1Dres, "StencilCondNo-innerUncut-Var0", 0.5));
            //ExpectedSlopes.Add((Solution.AdvancedSolvers.Testing.OpAnalysisBase.XAxisDesignation.Grid_1Dres, "StencilCondNo-innerCut-Var0", 0.5));

            //Solution.AdvancedSolvers.Testing.OpAnalysisBase.TestSlopes(Controls, ExpectedSlopes);

            var t = new ConditionNumberScalingTest();
            t.SetControls(Controls);
            t.CheckResults();


        }

       
    }
}
