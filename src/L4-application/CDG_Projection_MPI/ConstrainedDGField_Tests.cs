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
using System.Diagnostics;

namespace CDG_Projection_MPI {

    [TestFixture]
    public static class ConstrainedDGField_Tests {
        static void Main(string[] args) {
            Application.InitMPI();
            Application.DeleteOldPlotFiles();
            NoOfPatches = 20;
            Projection(ProjectionStrategy.patchwiseOnly, DGDeg: 2, incDeg: 1, Dim: 3, GridRes: 15, shape: GeomShape.Cube, narrowBand: true);
            //CDG_Projection_MPI.ConstrainedDGField_Tests.ProjectionTest_Global(2, 1, 2, 10, GeomShape.Sphere);
            //ProjectionPseudo1D(ProjectionStrategy.globalOnly, 10);
            //Test4Idempotence(3, 10);
            Application.FinalizeMPI();
        }

        //private static void ZeroPivotWithPardiso() {
        //    Projection(ProjectionStrategy.patchwiseOnly, 3, 10, GeomShape.Sphere);
        //}

        public enum GeomShape {
            Sphere = 0,
            Cube = 1
        }

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

        private static CellMask ComputeNarrowband(GridCommons grid, double radius, Func<double[], double> FGen) {
            var cells = grid.GridData.Cells;
            long i0 = grid.GridData.CellPartitioning.i0;
            long iE = grid.GridData.CellPartitioning.iE;
            int L = grid.CellPartitioning.LocalLength;
            double CellLength = grid.GridData.Cells.h_minGlobal;
            int D = grid.SpatialDimension;
            var barray = new BitArray(L);
            double tol = CellLength * 2; // Narrowband of CLength *2 width

            for (int iCell = 0; iCell < L; iCell++) {
                var Nodes = grid.GridData.GlobalNodes.GetValue_Cell(grid.GridData.Cells.GetRefElement(iCell).Vertices, iCell, 1).ExtractSubArrayShallow(0, -1, -1);

                bool isNegtive = false, isPositive = false;
                for(int k = 0; k <= Nodes.NoOfRows; k++) {
                    Vector pt;
                    if(k < Nodes.NoOfRows)
                        pt = Nodes.GetRowPt(k);
                    else
                        pt = cells.GetCenter(iCell);

                    double val = FGen(pt);
                    isNegtive = val <= tol;
                    isPositive = val >= -tol;
                }
                barray[iCell] = isNegtive && isPositive;
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

#if DEBUG
        [Test]
        static public void ProjectionTest_Patchwise( 
            [Values(2)] int DGDeg, 
            [Values(1)] int incDeg, 
            [Values(2,3)] int Dim,
            [Values(10,15)] int GridRes, 
            [Values(GeomShape.Cube)] GeomShape shape, 
            [Values(2, 3)] int __NoOfPatches) {

            NoOfPatches = __NoOfPatches;
            Projection(ProjectionStrategy.patchwiseOnly, DGDeg, incDeg, Dim, GridRes, shape, false);
        }


        [Test]
        static public void ProjectionTest_Global( 
            [Values(2)] int DGDeg, 
            [Values(1)] int incDeg, 
            [Values(2)] int Dim,
            [Values(10,15,20)] int GridRes, 
            [Values(GeomShape.Sphere, GeomShape.Cube)] GeomShape shape) {

            Projection(ProjectionStrategy.globalOnly, DGDeg, incDeg, Dim, GridRes, shape, false);
        }
#else
        [Test]
        static public void ProjectionTest_Patchwise( 
            [Values(2)] int DGDeg, 
            [Values(1)] int incDeg, 
            [Values(3)] int Dim,
            [Values(15)] int GridRes, 
            [Values(GeomShape.Cube)] GeomShape shape, 
            [Values(-1, 4)] int __NoOfPatches) {
            NoOfPatches = __NoOfPatches;
            Projection(ProjectionStrategy.patchwiseOnly, DGDeg, incDeg, Dim, GridRes, shape, false);
        }


        [Test]
        static public void ProjectionTest_Global( 
            [Values(2)] int DGDeg, 
            [Values(1)] int incDeg, 
            [Values(3)] int Dim,
            [Values(10, 15, 20)] int GridRes, 
            [Values(GeomShape.Sphere, GeomShape.Cube)] GeomShape shape) {

            Projection(ProjectionStrategy.globalOnly, DGDeg, incDeg, Dim, GridRes, shape, false);
        }

#endif



        static public void Projection(ProjectionStrategy PStrategy, int DGDeg, int incDeg, int Dim, int GridRes, GeomShape shape, bool narrowBand) {
            // Arrange -- Grid
            var grid = CreateGrid(GridRes, Dim);
            // Arrange -- target field
            double particleRad = 0.625;
            var DGBasis = new Basis(grid, DGDeg);
            var CDGBasis = new Basis(grid, DGDeg + incDeg);
            
            var FGenerator = FuncGenerator(shape, grid.GridData.SpatialDimension);
            var projFunc = FGenerator(particleRad);

            var mask = narrowBand ? ComputeNarrowband(grid, particleRad, projFunc) : CellMask.GetFullMask(grid.GridData);
            Console.WriteLine("Cells in mask: " + mask.NoOfItemsLocally.MPISum());
            Console.WriteLine("edges in mask: " + mask.AllEdges().NoOfItemsLocally.MPISum());

            // Arrange -- source field
            var field = new SinglePhaseField(DGBasis, "DG");
            var ProjF = NonVectorizedScalarFunction.Vectorize(projFunc);
            field.ProjectField(ProjF);

            // Act -- Do the CG-projection
            SinglePhaseField BestApprox = new SinglePhaseField(CDGBasis, "CG");
            using(var CDGTestField = Factory(CDGBasis, mask, PStrategy)) {

                CDGTestField.ScheissDrauf(field);
                double jumpNorm_before = CDGTestField.CheckLocalProjection(mask, false);
                Console.WriteLine("jump norm of initial: " + jumpNorm_before);

                CDGTestField.ProjectDGField(field);
                CDGTestField.AccToDGField(1.0, BestApprox, mask);

                //BoSSS.Solution.Tecplot.Tecplot.PlotFields(new DGField[] { field, BestApprox }, "PrjSnap", 0.0, 3);
                
                double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, false);
                Console.WriteLine("jump Norm after total projection: " + jumpNorm_after);
                Assert.LessOrEqual(jumpNorm_after, Math.Max(jumpNorm_before * 1E-8, 1e-10), // in case of sphere, we start wit a low jump norm.
                    "Jump norm not sufficiently reduced.");


            }
            
        }

        /// <summary>
        /// Tests that the global projection operator is idempotent,
        /// i.e. executing it twice gives the same result as executing once.
        /// (A projector must be, per definition, linear and idempotent.)
        /// </summary>
        [Test]
        static public void Test4Idempotence([Values(2, 3)] int Dim, [Values(10)] int GridRes) {
            GeomShape shape = GeomShape.Cube; // use something with jumps on the initial value
            ProjectionStrategy PStrategy = ProjectionStrategy.globalOnly; // of course, only global projection is idempotent (local is only an approx)

            // Arrange -- Grid
            var grid = CreateGrid(GridRes, Dim);
            // Arrange -- target field
            double particleRad = 0.625;
            var DGBasis = new Basis(grid, 2);
            var CDGBasis = new Basis(grid, 3);

            var mask = CellMask.GetFullMask(grid.GridData);

            // Arrange -- source field
            var FGenerator = FuncGenerator(shape, grid.GridData.SpatialDimension);
            var projFunc = FGenerator(particleRad);
            var field = new SinglePhaseField(DGBasis, "DG");
            var ProjF = NonVectorizedScalarFunction.Vectorize(projFunc);
            field.ProjectField(ProjF);

            // Act -- Do the CG-projection
            SinglePhaseField BestApprox = null;
            double FirstCorrection = 0;
            using(var CDGTestField = Factory(CDGBasis, mask, PStrategy)) {
                for(int i = 0; i < 3; i++) {
                    SinglePhaseField Rest;
                    if(BestApprox == null) {
                        Rest = field.CloneAs();
                    } else {
                        Rest = BestApprox.CloneAs();
                        Rest.AccLaidBack(-1.0, field);
                        Rest.Scale(-1.0);
                    }

                    
                    CDGTestField.ProjectDGField(Rest);

                    if(BestApprox == null) {
                        BestApprox = new SinglePhaseField(CDGTestField.Basis, "C0approx");
                        CDGTestField.AccToDGField(1.0, BestApprox);
                        double l2Correction = BestApprox.L2Norm();
                        FirstCorrection = l2Correction;
                        Console.WriteLine("===================== L2 norm correction " + i + ": " + l2Correction);
                    } else {
                        var next = new SinglePhaseField(CDGTestField.Basis, "C0approx");
                        CDGTestField.AccToDGField(1.0, next);
                        double l2Correction = next.L2Norm();
                        next.Acc(1.0, BestApprox);
                        BestApprox = next;
                        Console.WriteLine("===================== L2 norm correction " + i + ": " + l2Correction);
                        Assert.LessOrEqual(l2Correction, FirstCorrection * 1e-8, "Idempotence test failed; after the first projection, no more change should occur");
                    }
                }


                //BoSSS.Solution.Tecplot.Tecplot.PlotFields(new DGField[] { field, BestApprox }, "PrjSnap-" + i, 0.0, 3);
            }
        }

        static int NoOfPatches = -1;

        static public ConstrainedDGFieldMk3 Factory(Basis b, CellMask __domainLimit, ProjectionStrategy s) {
            switch(s) {
                case ProjectionStrategy.globalOnly: return new ConstrainedDGField_Global(b, __domainLimit);
                case ProjectionStrategy.patchwiseOnly: return new ConstrainedDgField_Patchwise(b, __domainLimit, NoOfPatches);
                default: throw new NotImplementedException();
            }
        }

       
        static public void ProjectionPseudo1D(ProjectionStrategy PStrategy, int GridRes) {
            // Arrange -- Grid
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            int numOfCells = GridRes * NoOfCores;
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCells);
            var grid = Grid2D.Cartesian2DGrid(xNodes, new double[] { -1, 1 });
            var mask = CellMask.GetFullMask(grid.GridData);

            // Arrange -- initial DG projection
            var BasisDG = new Basis(grid, 0);
            var BasisCDG = new Basis(grid, 2);

            
            //var CDGTestField =  Factory(BasisCDG, mask, PStrategy); 
            Func<double[], double> ProjFunc = (double[] X) => X[0];
            var field = new SinglePhaseField(BasisDG, "DG");
            var ProjF = NonVectorizedScalarFunction.Vectorize(ProjFunc);
            field.ProjectField(ProjF);

           
            SinglePhaseField BestApprox = null;
            for(int i = 0; i < 5; i++) {
                SinglePhaseField Rest;
                if(BestApprox == null) {
                    Rest = field.CloneAs();
                } else {
                    Rest = BestApprox.CloneAs();
                    Rest.AccLaidBack(-1.0, field);
                    Rest.Scale(-1.0);
                }


                var CDGTestField = Factory(BasisCDG, mask, PStrategy); 
                CDGTestField.ScheissDrauf(Rest);
                double jumpNorm_before = CDGTestField.CheckLocalProjection(mask, false);
                Console.WriteLine("jump norm of initial: " + jumpNorm_before);
                
               
                CDGTestField.ProjectDGField(Rest);
                
                
                double jumpNorm_after = CDGTestField.CheckLocalProjection(mask, false);
                Console.WriteLine("jump Norm after total projection: " + jumpNorm_after);

                if(BestApprox == null) {
                    BestApprox = new SinglePhaseField(CDGTestField.Basis, "C0approx");
                    CDGTestField.AccToDGField(1.0, BestApprox);
                    double l2Correction = BestApprox.L2Norm();
                    double dist = BestApprox.L2Distance(field);
                    Console.WriteLine("===================== L2 norm correction " + i + ": " + l2Correction + "\t" + dist);
                } else {
                    var next = new SinglePhaseField(CDGTestField.Basis, "C0approx");
                    CDGTestField.AccToDGField(1.0, next);
                    double l2Correction = next.L2Norm();
                    next.Acc(1.0, BestApprox);
                    BestApprox = next;
                    double dist = BestApprox.L2Distance(field);

                    Console.WriteLine("===================== L2 norm correction " + i + ": " + l2Correction + "\t" + dist);
                }



                BoSSS.Solution.Tecplot.Tecplot.PlotFields(new DGField[] { field, BestApprox }, "PrjSnap-" + i, 0.0, 2);
            }

        }
    
    
    
         
    }
}
