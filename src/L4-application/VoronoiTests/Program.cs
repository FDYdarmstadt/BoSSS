using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using VoronoiTests.Database;
using VoronoiTests.Database.Session;

namespace VoronoiTests
{
    class Program
    {
        public static void Main() {
            TestBench selectedTest = availableTests[3];
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
            {new DGFieldArithmeticTests() },
            {new SessionIOTests() },
            {new GridIOTests() },
            {new BoSSSpadTests() }
        };
        
    }
}
