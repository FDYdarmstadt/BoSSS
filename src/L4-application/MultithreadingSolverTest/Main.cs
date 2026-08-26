using System;
using System.Linq;

namespace BoSSS.Application.MultithreadingSolverTest {
    class Program {
        public static void Main() {
            BoSSS.Solution.Application.InitMPI();

            // get one test from the list and run it, to check if the test framework is working at all; the actual tests are run in the NUnit test class
            var test = MultithreadingSolverTestcases.All.First();
            var c = test.CreateControl();
            var s = test.RunControl(c);

            Console.WriteLine($"Test {test.Name} completed with time {s}");
        }
    }
}
