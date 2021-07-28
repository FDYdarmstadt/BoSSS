using System;
using NUnit.Framework;
using ilPSP.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Utils;

namespace CDG_Projection_MPI {

    [TestFixture]
    class ConstrainedDGField_Tests : ConstrainedDGField {
        static void Main(string[] args) {
            Application.InitMPI();
            InterProcProjectionTest();
            Application.FinalizeMPI();
        }

        private ConstrainedDGField_Tests(Basis b) : base(b) { }

        [Test]
        static public void InterProcProjectionTest() {
            // Arrange -- Grid
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            int numOfCells = 4 * NoOfCores;
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCells + 1);
            var grid = Grid3D.Cartesian3DGrid(xNodes, xNodes, xNodes);

            // Arrange -- target field
            var Basis = new Basis(grid, 2);
            var CDGTestField = new ConstrainedDGField_Tests(Basis);
            var mask = CellMask.GetFullMask(grid.GridData);

            // Arrange -- source field
            var field = new SinglePhaseField(Basis);
            var idontcarefield = new SinglePhaseField(Basis,"bla");
            var quadrature = new CellQuadratureScheme(true, mask);
            double particleRad = 0.261;
            double angle = 0;
            Func<double[], double> projFunc = (double[] X) =>
                -Math.Max(Math.Abs((X[0]) * Math.Cos(angle) - (X[1]) * Math.Sin(angle)),
                Math.Max(Math.Abs((X[0]) * Math.Sin(angle) + (X[1]) * Math.Cos(angle)),
                Math.Abs(X[2]))) + particleRad;
            var ProjF = NonVectorizedScalarFunction.Vectorize(projFunc);
            field.ProjectField(ProjF);

            // Act -- Do the CG-projection
            CDGTestField.SetDGCoordinatesOnce(mask, field);
            double jumpNorm = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump norm of initial: "+jumpNorm);
            
            CDGTestField.ProjectDGField_patchwise(field, mask, 1);
            double jumpNorm_before = CDGTestField.CheckLocalProjection(mask, false);
            
            //CDGTestField.ProjectDGFieldOnPatch(field, mask);
            //double jumpNorm_total = CDGTestField.CheckLocalProjection(mask, false);
            //Console.WriteLine("jump Norm after total projection: "+ jumpNorm_total);

            CDGTestField.InterProcessProjectionBranch(mask, field);
            double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, true);

            // Assert -- effective inter process projection
            Assert.IsTrue(jumpNorm_after < jumpNorm_before);
        }
    }
}
