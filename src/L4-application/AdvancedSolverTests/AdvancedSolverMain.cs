using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution;

namespace AdvancedSolverTests {
    class AdvancedSolverMain {
        public static void Main() {
            TestBench selectedTest = availableTests["SubBlockTests"];
            RunTest(selectedTest);
        }

        static void RunTest(TestBench test) {
            test.SetUp();
            test.Run();
            test.TearDown();
        }

        static readonly Dictionary<string, TestBench> availableTests = new Dictionary<string, TestBench>()
        {
            {"SubBlockTests", new AdvancedSolverTests.SubBlocking.SubBlockTests() },
        };

    }
}
