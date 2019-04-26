using System.Collections.Generic;

namespace VoronoiTests
{
    class Program
    {
        public static void Main() {
            TestBench selectedTest = availableTests[4];
            RunTest(selectedTest);
        }

        static void RunTest(TestBench test)
        {
            test.SetUp();
            test.Run();
            test.TearDown();
        }

        static readonly List< TestBench> availableTests = new List<TestBench>
        {
            new Database.Session.DGFieldArithmeticTests(),
            new Database.Session.SessionIOTests(),
            new Database.GridIOTests(),
            new Database.Session.BoSSSpadTests(),
            new GridCreation.GridTests()
        };
    }
}
