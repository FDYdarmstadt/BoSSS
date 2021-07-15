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
            CellMask mask = null; // CellMask.GetFullMask(grid.GridData);

            // Arrange -- source field
            var field = new SinglePhaseField(Basis);
            var idontcarefield = new SinglePhaseField(Basis,"bla");
            var quadrature = new CellQuadratureScheme(true, mask);
            double particleRad = 0.2;

            Func<double[], double> projFunc = (double[] X) => -X[0] * X[0] - X[1] * X[1] - X[2] * X[2] + particleRad * particleRad;
            var ProjF = NonVectorizedScalarFunction.Vectorize(projFunc);
            field.ProjectField(ProjF);

            // Act -- Do the CG-projection
            CDGTestField.SetDGCoordinatesOnce(mask, field);
            CDGTestField.ProjectDGField_patchwise(field,mask,1);
            double jumpNorm_before = CDGTestField.CheckLocalProjection(mask,false);
            CDGTestField.InterProcessProjectionBranch(mask, field);
            double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, true);

            // Assert -- effective inter process projection
            Assert.IsTrue(jumpNorm_after < jumpNorm_before * 0.01);
        }
    }
}
