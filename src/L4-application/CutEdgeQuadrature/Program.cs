using BoSSS.Solution;
using CutEdgeQuadrature.QuadratureTests;
using System;

namespace CutEdgeQuadrature {
    class Program {

        static IEdgeQuadratureTest3D[] tests3D = new IEdgeQuadratureTest3D[] {
            new SkewLinearLevelSet(),
            new CubicLevelSet(),
            new QuadraticLevelSet(),
            new LinearLevelSet(),
        };

        static void Main(string[] args) {
            Application.InitMPI();
            Tester3D tester = new Tester3D();
            foreach(var test in tests3D) {
                int order = 3;
                double error = tester.Test(test, order);
                Console.WriteLine($"Order: {order}, Error: {error}");
            }
            Console.WriteLine("Press any key...");
            Console.ReadKey();
        }
    }
}
