using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Tests {
    
    //[TestFixture]
    static class LevelSetEvolverTest {

        //[Test]
        public static void Test(
            [Values(1, 2, 3, 4)] int degree,
            [Values(0.0, 0.1)] double agglomerationTreshold
            ) {
            ZLS_Control control = ZwoLevelSetSolver.ControlFiles.HardCodedControl.QuasiStationaryDroplet(degree, agglomerationTreshold);
            
            Test<ZLS_Control> levelSetEvolverTest = new Test<ZLS_Control>(control);
            levelSetEvolverTest.AddSolutionTestField(
                new SolutionTestField {

            });

            TestRunner<ZLS, ZLS_Control>.SolverTest(levelSetEvolverTest);
        }

        
    }
}
