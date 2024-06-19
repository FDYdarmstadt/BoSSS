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
                test.MomentFittingVariant = BoSSS.Foundation.XDG.XQuadFactoryHelperBase.MomentFittingVariants.Algoim;
                Console.WriteLine($"### Testing: {test.GetType()} ###");
                int[] orders = new int[] { 1,2,3,4};

                foreach (int order in orders) { 
                    double error = tester.Test(test, order);
                    Console.WriteLine($"Order: {order}, Error: {error}");
                }
            }

            Console.WriteLine("Press any key...");
            Console.ReadKey();
        }
    }
}
