using System.Collections.Generic;

namespace VoronoiTests
{
    class Program
    {
        public static void Main() {
            TestBench selectedTest = availableTests["Solver_IpPoisson"];
            RunTest(selectedTest);
        }

        static void RunTest(TestBench test)
        {
            test.SetUp();
            test.Run();
            test.TearDown();
        }

        static readonly Dictionary<string, TestBench> availableTests = new Dictionary<string, TestBench>()
        {
            {"ArithmeticTests", new Database.Session.DGFieldArithmeticTests() },
            {"SessionIOTests", new Database.Session.SessionIOTests() },
            {"GridIOTests", new Database.GridIOTests()},
            {"BoSSSpadTests", new Database.Session.BoSSSpadTests()},
            {"GridTests", new GridCreation.GridTests()},
            {"Solver_IpPoisson", new Solver.IpPoissonTests()}
        };
    }
}
