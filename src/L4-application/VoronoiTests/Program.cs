using System.Collections.Generic;

namespace VoronoiTests {
    class Program {
        public static void Main() {

            TestBench selectedTest = availableTests["Boundary Tests"];
            RunTest(selectedTest);
        }

        static void RunTest(TestBench test) {
            test.SetUp();
            test.Run();
            test.TearDown();
        }

        static readonly Dictionary<string, TestBench> availableTests = new Dictionary<string, TestBench>()
        {
            {"Arithmetic Tests", new Database.Session.DGFieldArithmeticTests() },
            {"Session IO Tests", new Database.Session.SessionIOTests() },
            {"Grid IO Tests", new Database.GridIOTests()},
            {"BoSSSpad Tests", new Database.Session.BoSSSpadTests()},
            {"Grid Creation", new Grid.GridCreationTests()},
            {"Solver IpPoisson", new Solver.IpPoissonTests()},
            {"Grid Movement", new Grid.MovementTests()},
            {"Boundary Tests", new Grid.BoundaryTests()},
            {"Mapping", new Grid.NodeMappingTests()},
        };
    }
}
