using System;
using NUnit.Framework;
using ilPSP.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.ConstrainedDGprojection;
using System.Collections.Generic;
using MPI.Wrappers;
using System.Collections;
using System.Linq;

namespace CDG_Projection_MPI {

    [TestFixture]
    class ConstrainedDGField_Tests : ConstrainedDGField {
        static void Main(string[] args) {
            Application.InitMPI();
            ProjectionIn3D(ProjectionStrategy.patchwiseOnly, 10);
            Application.FinalizeMPI();
        }

        private ConstrainedDGField_Tests(Basis b) : base(b) { }

        private static Grid3D Create3Dgrid(int Res) {
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            int numOfCells = Res * NoOfCores;
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCells + 1);
            var grid = Grid3D.Cartesian3DGrid(xNodes, xNodes, xNodes);
            return grid;
        }

        private static Grid2D Create2DGrid(int Res) {
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            int numOfCells = Res * NoOfCores;
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCells + 1);
            var grid = Grid2D.Cartesian2DGrid(xNodes, xNodes);
            return grid;
        }

        private static CellMask ComputeNarrowband(GridCommons grid, double radius) {
            var cells = grid.GridData.Cells;
            long i0 = grid.GridData.CellPartitioning.i0;
            long iE = grid.GridData.CellPartitioning.iE;
            int L = grid.CellPartitioning.LocalLength;
            double CellLength = grid.GridData.Cells.h_minGlobal;
            int D = grid.SpatialDimension;
            var barray = new BitArray(L);
            double tol = CellLength * 2; // Narrowband of CLength *2 width
            var UpperBnd = GenFunc(D, radius + tol);
            var LowerBnd = GenFunc(D, radius - tol);

            for (int iCell = 0; iCell < L; iCell++) {
                var center = cells.GetCenter(iCell);
                double[] tmp = center.GetSubVector(0, D);
                barray[iCell] = LowerBnd(tmp) < 0 && UpperBnd(tmp) > 0;
            }
            return new CellMask(grid.GridData, barray);
        }

        private static Func<double[], double> GenFunc(int Dim, double particleRad) {
            double angle = 1;
            int power = 10;
            switch (Dim) {
                case 2:
                    return (double[] X) =>
                        -Math.Max(Math.Abs((X[0]) * Math.Cos(angle) - (X[1]) * Math.Sin(angle)),
                        Math.Abs((X[0]) * Math.Sin(angle) + (X[1]) * Math.Cos(angle))) + particleRad;
                case 3:
                    return (double[] X) =>
                    -Math.Pow(Math.Pow((X[0]) * Math.Cos(angle) - (X[1]) * Math.Sin(angle), power)
                    + Math.Pow((X[0]) * Math.Sin(angle) + (X[1]) * Math.Cos(angle), power)
                    + Math.Pow(X[2], power), 1.0 / power)
                    + particleRad;
                    //Func<double[], double> projFunc = (double[] X) =>
                    //    -Math.Max(Math.Abs((X[0]) * Math.Cos(angle) - (X[1]) * Math.Sin(angle)),
                    //    Math.Max(Math.Abs((X[0]) * Math.Sin(angle) + (X[1]) * Math.Cos(angle)),
                    //    Math.Abs(X[2]))) + particleRad;
                default:
                throw new NotSupportedException();
            }
            
        }

        private static void PlotDGFields(DGField[] input, GridData gridData) {
            var Basis = new Basis(gridData,0);
            var listofinputs = input.Select(s => new SinglePhaseField(Basis, "prj_" + s));
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(listofinputs, "PrjMasks.plt", 0.0, 0);
        }

        static public void ProjectionIn3D(ProjectionStrategy PStrategy, int GridRes) {
            // Arrange -- Grid
            var grid = Create3Dgrid(GridRes);
            // Arrange -- target field
            double particleRad = 0.6;
            var Basis = new Basis(grid, 2);
            var CDGTestField = new ConstrainedDGField_Tests(Basis);
            //var mask = CellMask.GetFullMask(grid.GridData);
            var mask = ComputeNarrowband(grid, particleRad);

            // Arrange -- source field
            var projFunc = GenFunc(grid.SpatialDimension, particleRad);
            var field = new SinglePhaseField(Basis, "DG");
            var ProjMask = new SinglePhaseField(Basis, "ProjectionMask");
            var ProjF = NonVectorizedScalarFunction.Vectorize(projFunc);
            field.ProjectField(ProjF);
            var list = new List<SinglePhaseField>();
            ProjMask.AccConstant(1.0, mask);
            list.Add(ProjMask);
            list.Add(field.CloneAs());

            // Act -- Do the CG-projection
            CDGTestField.Setup();
            CDGTestField.SetDGCoordinatesOnce(mask, field);
            double jumpNorm_before = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump norm of initial: " + jumpNorm_before);
            switch (PStrategy) {
                case ProjectionStrategy.globalOnly:
                    CDGTestField.ProjectDGFieldGlobal(mask);
                break;
                case ProjectionStrategy.patchwiseOnly:                  
                    var PrjMasks = CDGTestField.ProjectDGField_patchwise(mask,1);
                    PlotDGFields(PrjMasks,grid.GridData);
                break;
            }
            double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, false);

            //CDGTestField.ProjectDGFieldOnPatch(field, mask);
            //double jumpNorm_total = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump Norm after total projection: "+ jumpNorm_after);

            field.Clear();
            CDGTestField.AccToDGField(1.0, field);
            field.Identification = "CG";
            list.Add(field);
            //BoSSS.Solution.Tecplot.Tecplot.PlotFields(list, "BLargh.plt", 0.0, 1);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(CDGTestField.ProjectionSnapshots, "PrjSnap.plt", 0.0, 1);
            // Assert -- effective inter process projection
            // Assert.IsTrue(jumpNorm_after < jumpNorm_before * 1E-8);
        }


        static public void ProjectionIn2D(ProjectionStrategy PStrategy, int GridRes) {
            // Arrange -- Grid
            var grid = Create2DGrid(GridRes);
            double particleRad = 0.5;
            // Arrange -- target field
            var Basis = new Basis(grid, 2);
            var CDGTestField = new ConstrainedDGField_Tests(Basis);
            var mask = CellMask.GetFullMask(grid.GridData);

            // Arrange -- source field
            var field = new SinglePhaseField(Basis, "DG");
            var projFunc = GenFunc(grid.SpatialDimension, particleRad);

            var ProjF = NonVectorizedScalarFunction.Vectorize(projFunc);
            field.ProjectField(ProjF);
            var list = new List<SinglePhaseField>();
            list.Add(field.CloneAs());

            // Act -- Do the CG-projection
            CDGTestField.Setup();
            CDGTestField.SetDGCoordinatesOnce(mask, field);
            double jumpNorm_before = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump norm of initial: " + jumpNorm_before);

            switch (PStrategy) {
                case ProjectionStrategy.globalOnly:
                CDGTestField.ProjectDGFieldGlobal(mask);
                break;
                case ProjectionStrategy.patchwiseOnly:
                CDGTestField.ProjectDGField_patchwise(mask);
                break;
            }

            double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump Norm after total projection: " + jumpNorm_after);

            field.Clear();
            CDGTestField.AccToDGField(1.0, field);
            field.Identification = "CG";
            list.Add(field);
            //BoSSS.Solution.Tecplot.Tecplot.PlotFields(list, "BLargh.plt", 0.0, 0);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(CDGTestField.ProjectionSnapshots, "PrjSnap.plt", 0.0, 1);
            // Assert -- effective inter process projection
            // Assert.IsTrue(jumpNorm_after < jumpNorm_before * 1E-8);
        }
    }
}
