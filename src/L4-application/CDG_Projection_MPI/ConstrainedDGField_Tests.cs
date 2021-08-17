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
using BoSSS.Foundation.IO;
using ilPSP;
using BoSSS.Solution.Statistic;

namespace CDG_Projection_MPI {

    [TestFixture]
    class ConstrainedDGField_Tests : ConstrainedDGField {
        static void Main(string[] args) {
            Application.InitMPI();
            //Application.DeleteOldPlotFiles();
            //Projection(ProjectionStrategy.patchwiseOnly, 3, 12, GeomShape.Cube);
            ProjectionPseudo1D(ProjectionStrategy.patchwiseOnly, 10);
            Application.FinalizeMPI();
        }

        private static void ZeroPivotWithPardiso() {
            Projection(ProjectionStrategy.patchwiseOnly, 3, 10, GeomShape.Sphere);
        }

        public enum GeomShape {
            Sphere = 0,
            Cube = 1
        }

        private ConstrainedDGField_Tests(Basis b) : base(b) { }


        private static GridCommons CreateGrid(int Res, int Dim) {
            switch(Dim) {
                case 1: return Create1DGrid(Res);
                case 2: return Create2DGrid(Res);
                case 3: return Create3Dgrid(Res);
                default: throw new ArgumentOutOfRangeException("Un-supported spatial dimension: " + Dim);
            }
        }


        private static Grid3D Create3Dgrid(int Res) {
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            //int numOfCells = Res * NoOfCores;
            int numOfCells = Res;
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

        private static Grid1D Create1DGrid(int Res) {
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            int numOfCells = Res * NoOfCores;
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCells + 1);
            var grid = Grid1D.LineGrid(xNodes);
            return grid;
        }

        private static CellMask ComputeNarrowband(GridCommons grid, double radius, Func<double, Func<double[], double>> FGen) {
            var cells = grid.GridData.Cells;
            long i0 = grid.GridData.CellPartitioning.i0;
            long iE = grid.GridData.CellPartitioning.iE;
            int L = grid.CellPartitioning.LocalLength;
            double CellLength = grid.GridData.Cells.h_minGlobal;
            int D = grid.SpatialDimension;
            var barray = new BitArray(L);
            double tol = CellLength * 2; // Narrowband of CLength *2 width
            var UpperBnd = FGen(radius + tol);
            var LowerBnd = FGen(radius - tol);

            for (int iCell = 0; iCell < L; iCell++) {
                var center = cells.GetCenter(iCell);
                double[] tmp = center.GetSubVector(0, D);
                barray[iCell] = LowerBnd(tmp) < 0 && UpperBnd(tmp) > 0;
            }
            return new CellMask(grid.GridData, barray);
        }

        private static Func<double[], double> GenCube(int Dim, double particleRad) {
            int power = 10;
            double angle = 0.7;
            switch(Dim) {
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
                //return (double[] X) =>
                //        -Math.Max(Math.Abs((X[0]) * Math.Cos(angle) - (X[1]) * Math.Sin(angle)),
                //        Math.Max(Math.Abs((X[0]) * Math.Sin(angle) + (X[1]) * Math.Cos(angle)),
                //        Math.Abs(X[2]))) + particleRad;
                default:
                throw new NotSupportedException();
            }

        }

        private static Func<double[], double> GenSphere(int Dim, double particleRad) {
            switch (Dim) {
                case 2:
                return (double[] X) => -X[0] * X[0] - X[1] * X[1] + particleRad * particleRad;
                case 3:
                return (double[] X) => -X[0] * X[0] - X[1] * X[1] - X[2] * X[2] + particleRad * particleRad;
                default:
                throw new NotSupportedException();
            }
        }

        private static Func<double, Func<double[], double>> FuncGenerator(GeomShape shape, int Dim) {
            switch (shape) {
                case GeomShape.Sphere:
                return (double Radius) => GenSphere(Dim,Radius);
                case GeomShape.Cube:
                return (double Radius) => GenCube(Dim, Radius);
                default:
                throw new NotImplementedException();
            }
        }

        static public void Projection(ProjectionStrategy PStrategy, int Dim, int GridRes, GeomShape shape) {
            // Arrange -- Grid
            var grid = CreateGrid(GridRes, Dim);
            // Arrange -- target field
            double particleRad = 0.625;
            var DGBasis = new Basis(grid, 2);
            var CDGBasis = new Basis(grid, 2);
            var CDGTestField = new ConstrainedDGField_Tests(CDGBasis);
            //var mask = CellMask.GetFullMask(grid.GridData);
            var FGenerator = FuncGenerator(shape, grid.GridData.SpatialDimension);
            var mask = ComputeNarrowband(grid, particleRad, FGenerator);
            Console.WriteLine("Cells in mask: " + mask.NoOfItemsLocally.MPISum());
            Console.WriteLine("edges in mask: " + mask.AllEdges().NoOfItemsLocally.MPISum());

            // Arrange -- source field
            var projFunc = FGenerator(particleRad);
            var field = new SinglePhaseField(DGBasis, "DG");
            var ProjF = NonVectorizedScalarFunction.Vectorize(projFunc);
            field.ProjectField(ProjF);

            // Act -- Do the CG-projection
            CDGTestField.Setup();
            CDGTestField.SetDGCoordinatesOnce(field, CellMask.GetFullMask(grid.GridData));
            double jumpNorm_before = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump norm of initial: " + jumpNorm_before);
            switch (PStrategy) {
                case ProjectionStrategy.globalOnly:
                    CDGTestField.ProjectDGFieldGlobal(mask);
                break;
                case ProjectionStrategy.patchwiseOnly:                  
                    var PrjMasks = CDGTestField.ProjectDGField_patchwise(mask, 2); 
                    //BoSSS.Solution.Tecplot.Tecplot.PlotFields(PrjMasks, "PrjMasks.plt", 0.0, 0);
                break;
            }
            double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump Norm after total projection: "+ jumpNorm_after);

            if (CDGTestField.ProjectionSnapshots != null) 
                BoSSS.Solution.Tecplot.Tecplot.PlotFields(CDGTestField.ProjectionSnapshots, "PrjSnap.plt", 0.0, 2);
            
            // Assert -- effective inter process projection
            // Assert.IsTrue(jumpNorm_after < jumpNorm_before * 1E-8);
        }

        /*
        static public void ProjectionIn2D(ProjectionStrategy PStrategy, int GridRes, GeomShape shape) {
            // Arrange -- Grid
            var grid = Create2DGrid(GridRes);
            double particleRad = 0.625;
            // Arrange -- target field
            var Basis = new Basis(grid, 2);
            var CDGTestField = new ConstrainedDGField_Tests(Basis);
            var FGenerator = FuncGenerator(shape, grid.GridData.SpatialDimension);
            var mask = ComputeNarrowband(grid,particleRad,FGenerator);

            // Arrange -- source field
            var field = new SinglePhaseField(Basis, "DG");
            
            var projFunc = FGenerator(particleRad);
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
                CDGTestField.ProjectDGField_patchwise(mask,1);
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
        */


        static public void ProjectionPseudo1D(ProjectionStrategy PStrategy, int GridRes) {
            // Arrange -- Grid
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            int numOfCells = GridRes * NoOfCores;
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCells);
            var grid = Grid2D.Cartesian2DGrid(xNodes, new double[] { -1, 1 });

            // Arrange -- initial DG projection
            var BasisDG = new Basis(grid, 0);
            var BasisCDG = new Basis(grid, 2);
            var CDGTestField = new ConstrainedDGField_Tests(BasisCDG);
            //double slope = 0.6;
            //Func<double[],double> ProjFunc = (double[] X) => X[0]* slope * (X[0]<=0.0?1:-1) + slope;
            Func<double[], double> ProjFunc = (double[] X) => X[0];
            var mask = CellMask.GetFullMask(grid.GridData);
            var field = new SinglePhaseField(BasisDG, "DG");
            var ProjF = NonVectorizedScalarFunction.Vectorize(ProjFunc);
            field.ProjectField(ProjF);

            // Act -- Do the CG-projection
            CDGTestField.Setup();
            CDGTestField.SetDGCoordinatesOnce(field, mask);
            double jumpNorm_before = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump norm of initial: " + jumpNorm_before);

            switch (PStrategy) {
                case ProjectionStrategy.globalOnly:
                CDGTestField.ProjectDGFieldGlobal(mask);
                break;
                case ProjectionStrategy.patchwiseOnly:
                var dbgFields = CDGTestField.ProjectDGField_patchwise(mask, 2);
                BoSSS.Solution.Tecplot.Tecplot.PlotFields(dbgFields, "debug", 0.0, 1);
                break;
            }

            double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, false);
            Console.WriteLine("jump Norm after total projection: " + jumpNorm_after);

            //var Aux1DGrid = Create1DGrid(GridRes);
            //var field1D = new SinglePhaseField(new Basis(Aux1DGrid, 2));
            double[] XCoor = GenericBlas.Linspace(xMin+1E-4, xMax - 1E-4, GridRes * 10);
            var inital = MultidimensionalArray.Create(XCoor.Length, 2);
            var interproc = MultidimensionalArray.Create(XCoor.Length, 2);
            var proclocal = MultidimensionalArray.Create(XCoor.Length, 2);
            var FieldDingen = new FieldEvaluation(grid.GridData);
            double LineplotAt = 0.0;
            for (int i=0;i<XCoor.Length;i++) {
                inital[i, 1] = XCoor[i];
                interproc[i, 1] = XCoor[i];
                proclocal[i, 1] = XCoor[i];
                //Results[i, 0]=FieldDingen.Evaluate(new Vector(XCoor[i],0.0),field);
                inital[i, 0] = field.ProbeAt(XCoor[i], LineplotAt);
                interproc[i, 0] = CDGTestField.ProjectionSnapshots[2].ProbeAt(XCoor[i], LineplotAt);
                proclocal[i, 0] = CDGTestField.ProjectionSnapshots.Last().ProbeAt(XCoor[i], LineplotAt);
            }

            //double alpha, IEnumerable< DGField > Flds, MultidimensionalArray Points, double beta, MultidimensionalArray Result, BitArray UnlocatedPoints = null, int[] LocalCellIndices = null

            inital.SaveToTextFile("initial");
            interproc.SaveToTextFile("interproc");
            proclocal.SaveToTextFile("proclocal");

            BoSSS.Solution.Tecplot.Tecplot.PlotFields(CDGTestField.ProjectionSnapshots, "PrjSnap.plt", 0.0, 2);
            // Assert -- effective inter process projection
            // Assert.IsTrue(jumpNorm_after < jumpNorm_before * 1E-8);
        }
    }
}
