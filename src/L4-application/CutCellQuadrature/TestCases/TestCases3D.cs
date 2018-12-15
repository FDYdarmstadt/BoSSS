/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;

namespace CutCellQuadrature.TestCases {

    abstract class Generic3DTestCase : TestCase<Shift3D> {

        private IEnumerable<Shift3D> shifts = null;

        public Generic3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift3D> AllShifts {
            get {
                if (shifts == null) {
                    List<Shift3D> myShifts = new List<Shift3D>();
                    Random rng = new Random(1);
                    //for (int i = 1; i < 11; i++) {
                    //for (int i = 1; i < 21; i++) {
                    //for (int i = 1; i < 26; i++) {
                    //for (int i = 1; i < 51; i++) {
                    for (int i = 1; i < 101; i++) {
                        double x = rng.NextDouble();
                        double y = rng.NextDouble();
                        double z = rng.NextDouble();

                        myShifts.Add(new Shift3D() {
                            OffsetX = x,
                            OffsetY = y,
                            OffsetZ = z
                        });
                    }
                    shifts = myShifts;
                }

                return shifts;
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            IGrid grid;
            switch (GridType) {
                case GridTypes.Structured:
                    int noOfCellsPerDirection;
                    switch (GridSize) {
                        case GridSizes.Tiny:
                            noOfCellsPerDirection = 5;
                            break;

                        case GridSizes.Small:
                            noOfCellsPerDirection = 10;
                            break;

                        case GridSizes.Normal:
                            noOfCellsPerDirection = 20;
                            break;

                        case GridSizes.Large:
                            noOfCellsPerDirection = 40;
                            break;

                        case GridSizes.Huge:
                            noOfCellsPerDirection = 80;
                            break;

                        default:
                            throw new Exception();
                    }
                    double[] nodes = GenericBlas.Linspace(-2.0, 2.0, noOfCellsPerDirection + 1);
                    grid = Grid3D.Cartesian3DGrid(nodes, nodes, nodes);
                    break;

                case GridTypes.Unstructured:
                    switch (GridSize) {
                        case GridSizes.Tiny:
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("523a0225-2c4f-46dd-ac75-86da5d547ad4"), db);
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("ec92c7eb-890e-468b-9728-5e5a3235822a"), db);
                            break;

                        case GridSizes.Small:
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("b55f328e-49dc-4e0b-aa0c-a8f3529adddc"), db);
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("45e51c76-e636-4bf1-b577-966d23aa744f"), db);
                            break;

                        case GridSizes.Normal:
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("189ac579-f9b9-4e5f-b847-4bc32b77c812"), db);
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("990cd7df-fe82-4c36-86b3-ea67d24057c7"), db);
                            break;

                        case GridSizes.Large:
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("51b0b6d8-f52a-45e3-a785-258fee6f573c"), db);
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("989ea7bb-86d4-4429-98a7-27471aab4bb8"), db);
                            break;

                        case GridSizes.Huge:
                            throw new NotImplementedException();

                        default:
                            throw new Exception();
                    }
                    break;

                default:
                    throw new Exception();
            }

            return grid;
        }

        public override double GridSpacing {
            get {
                return 0.2;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        public abstract int LevelSetDegree {
            get;
        }

        public abstract void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output);
    }

    abstract class Sphere3DTestCase : Generic3DTestCase {

        public Sphere3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                //return 5;
                return 2;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                //output[i] = 1.0 - Math.Sqrt(x * x + y * y + z * z);
                output[i] = 1.0 - (x * x + y * y + z * z);
            }
        }
    }

    class SphereVolume3DTestCase : Sphere3DTestCase, IVolumeTestCase {

        public SphereVolume3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 4.0 / 3.0 * Math.PI;
            }
        }
    }

    class SphereVolume3DNonSignedDistanceTestCase : SphereVolume3DTestCase {

        public SphereVolume3DNonSignedDistanceTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                output[i] = 1.0 - (x * x + y * y + z * z);
            }
        }
    }

    class ConeVolume3DTestCase : Sphere3DTestCase, IVolumeTestCase {

        public ConeVolume3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override void JumpingFieldSpeciesBInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                output[i] = 1.0 - Math.Sqrt(x * x + y * y + z * z);
            }
        }

        public override int IntegrandDegree {
            get {
                return 5;
            }
        }

        public override double Solution {
            get {
                return Math.PI / 3.0;
            }
        }
    }

    class InverseConeVolume3DTestCase : Sphere3DTestCase, IVolumeTestCase {

        public InverseConeVolume3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override void JumpingFieldSpeciesBInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                output[i] = Math.Sqrt(x * x + y * y + z * z);
            }
        }

        public override void JumpingFieldSpeciesAInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                output[i] = 1.0;
            }
        }

        public override int IntegrandDegree {
            get {
                return 5;
            }
        }

        public override double Solution {
            get {
                return Math.PI;
            }
        }
    }

    abstract class GenericTorus3DTestCase : Generic3DTestCase {

        //protected double R = 1.0;
        protected double R = 2.0;

        //protected double r = 0.25;
        protected double r = 1.0;

        public GenericTorus3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 10;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                output[i] = r * r - Math.Pow(Math.Sqrt(x * x + y * y) - R, 2.0) - z * z;
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            GridCommons grid;

            switch (GridType) {
                case GridTypes.Structured:
                    double[] xNodes = GenericBlas.Linspace(-3.2, 3.2, 33);
                    double[] yNodes = GenericBlas.Linspace(-3.2, 3.2, 33);
                    double[] zNodes = GenericBlas.Linspace(-1.4, 1.4, 15);
                    grid = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                    break;

                case GridTypes.Unstructured:
                    throw new NotImplementedException("To do: Create grid");

                default:
                    throw new NotImplementedException();
            }

            return grid;
        }
    }

    class TorusVolume3DTestCase : GenericTorus3DTestCase, IVolumeTestCase {

        public TorusVolume3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 2.0 * Math.PI * Math.PI * R * r * r;
            }
        }
    }

    class ConstantIntgreandSphereSurfaceIntegral3DTestCase : Sphere3DTestCase, ISurfaceTestCase {

        public ConstantIntgreandSphereSurfaceIntegral3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 4.0 * Math.PI;
            }
        }
    }

    class ConstantIntegrandSingleCellSphereCutout3DTestCase : Generic3DTestCase, ISurfaceTestCase {

        public ConstantIntegrandSingleCellSphereCutout3DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift3D> AllShifts {
            get {
                yield return new Shift3D() {
                    OffsetX = 0.0,
                    OffsetY = 0.0,
                    OffsetZ = 0.0
                };
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            if (GridType != GridTypes.Structured) {
                throw new NotImplementedException();
            }

            return Grid3D.Cartesian3DGrid(
                GenericBlas.Linspace(0.0, 0.5, 2),
                GenericBlas.Linspace(0.0, 0.5, 2),
                GenericBlas.Linspace(0.0, 0.5, 2));
        }

        public override double Solution {
            get {
                return 0.0;
            }
        }

        public override int LevelSetDegree {
            get {
                return 10;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                output[i] = 0.25 - Math.Sqrt(x * x + y * y + z * z);
            }
        }
    }

    class MinGibou3TorusQuadraticIntegrand : GenericTorus3DTestCase, ISurfaceTestCase {

        public MinGibou3TorusQuadraticIntegrand(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                //return 22.0 * Math.PI * Math.PI;
                return 11.0 * Math.PI * Math.PI;
            }
        }

        public override int IntegrandDegree {
            get {
                return 2;
            }
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                output[i] = x * x;
            }
        }
    }

    abstract class SingleCubeTestCase : Generic3DTestCase {

        public SingleCubeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift3D> AllShifts {
            get {
                yield return new Shift3D();
            }
        }

        public override int LevelSetDegree {
            get {
                return 1;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                //output[i] = x - 0.5;
                output[i] = x + 0.5 * z;
                //output[i] = x + y + z;


            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            switch (GridType) {
                case GridTypes.Structured:
                    return Grid3D.Cartesian3DGrid(
                        GenericBlas.Linspace(-2.0, 2.0, 2),
                        GenericBlas.Linspace(-1.0, 1.0, 2),
                        GenericBlas.Linspace(-1.0, 1.0, 2));
                //return Grid3D.Cartesian3DGrid(
                //    GenericBlas.Linspace(0.0, 2.0, 2),
                //    GenericBlas.Linspace(-1.0, 1.0, 2),
                //    GenericBlas.Linspace(-1.0, 1.0, 2));
                //return Grid3D.Cartesian3DGrid(
                //    GenericBlas.Linspace(0.0, 4.0, 2),
                //    GenericBlas.Linspace(-1.0, 1.0, 3),
                //    GenericBlas.Linspace(-1.0, 1.0, 4));

                case GridTypes.Unstructured:
                    throw new NotImplementedException();

                default:
                    throw new ArgumentException();
            }
        }
    }

    class SingleCubeSurfaceTestCase : SingleCubeTestCase, ISurfaceTestCase {

        public SingleCubeSurfaceTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                //return 4.0;
                //return 8.0;
                //return 16.0;
                return 2.0 * Math.Sqrt(5.0);
                //return 2.0 * Math.Sqrt(3.0);
                //return 16.0;
                //return 8.0 * Math.Sqrt(5.0);
                //return Math.Sqrt(5.0);
            }
        }
    }

    class SingleCubeVolumeTestCase : SingleCubeTestCase, IVolumeTestCase {

        public SingleCubeVolumeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                //return 2.0;
                //return 4.0;
                //return 14.0;
                return 16.0;
            }
        }
    }

    abstract class SingleCubeParaboloidTestCase : SingleCubeTestCase {

        public SingleCubeParaboloidTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];
                double z = input[i, 2];
                output[i] = 0.5 * x * x + x - 0.5 - y;
                //output[i] = 3.0 / 2.0 * x - y - 0.5;
            }
        }
    }


    class SingleCubeParaboloidVolumeTestCase : SingleCubeParaboloidTestCase, IVolumeTestCase {

        public SingleCubeParaboloidVolumeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 4.0 * 4.0 / 3.0 + 1.0;
            }
        }
    }

    class SingleCubeParaboloidSurfaceAreaTestCase : SingleCubeParaboloidTestCase, ISurfaceTestCase {

        public SingleCubeParaboloidSurfaceAreaTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                // 2 * int sqrt(t^2 + 2t + 2) from 0 to 1
                return 2.0 * 1.8100921403928758;
                //return Math.Sqrt(13.0);
            }
        }
    }

    abstract class SingleTetraTestCase : Generic3DTestCase {

        public SingleTetraTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift3D> AllShifts {
            get {
                yield return new Shift3D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            throw new NotImplementedException("Not ported yet");

            //if (GridType != GridTypes.Unstructured) {
            //    throw new NotImplementedException();
            //}

            //Tetra tetra = new Tetra();

            //long[] globalId = new long[] { 0 };
            //double[, ,] vertices = new double[1, 4, 3];
            ////vertices[0, 0, 0] = tetra.Vertices[0, 0];
            ////vertices[0, 0, 1] = tetra.Vertices[0, 1];
            ////vertices[0, 0, 2] = tetra.Vertices[0, 2];
            ////vertices[0, 1, 0] = tetra.Vertices[1, 0];
            ////vertices[0, 1, 1] = tetra.Vertices[1, 1];
            ////vertices[0, 1, 2] = tetra.Vertices[1, 2];
            ////vertices[0, 2, 0] = tetra.Vertices[2, 0];
            ////vertices[0, 2, 1] = tetra.Vertices[2, 1];
            ////vertices[0, 2, 2] = tetra.Vertices[2, 2];
            ////vertices[0, 3, 0] = tetra.Vertices[3, 0];
            ////vertices[0, 3, 1] = tetra.Vertices[3, 1];
            ////vertices[0, 3, 2] = tetra.Vertices[3, 2];

            ////vertices[0, 0, 0] = 2.0 * tetra.Vertices[0, 0];
            ////vertices[0, 0, 1] = 2.0 * tetra.Vertices[0, 1];
            ////vertices[0, 0, 2] = 2.0 * tetra.Vertices[0, 2];
            ////vertices[0, 1, 0] = 2.0 * tetra.Vertices[1, 0];
            ////vertices[0, 1, 1] = 2.0 * tetra.Vertices[1, 1];
            ////vertices[0, 1, 2] = 2.0 * tetra.Vertices[1, 2];
            ////vertices[0, 2, 0] = 2.0 * tetra.Vertices[2, 0];
            ////vertices[0, 2, 1] = 2.0 * tetra.Vertices[2, 1];
            ////vertices[0, 2, 2] = 2.0 * tetra.Vertices[2, 2];
            ////vertices[0, 3, 0] = 2.0 * tetra.Vertices[3, 0];
            ////vertices[0, 3, 1] = 2.0 * tetra.Vertices[3, 1];
            ////vertices[0, 3, 2] = 2.0 * tetra.Vertices[3, 2];

            ////vertices[0, 0, 0] = 2.0 * tetra.Vertices[0, 0];
            ////vertices[0, 0, 1] = 2.0 * tetra.Vertices[0, 1];
            ////vertices[0, 0, 2] = 2.0 * tetra.Vertices[0, 2];
            ////vertices[0, 1, 0] = tetra.Vertices[1, 0];
            ////vertices[0, 1, 1] = tetra.Vertices[1, 1];
            ////vertices[0, 1, 2] = tetra.Vertices[1, 2];
            ////vertices[0, 2, 0] = tetra.Vertices[2, 0];
            ////vertices[0, 2, 1] = tetra.Vertices[2, 1];
            ////vertices[0, 2, 2] = tetra.Vertices[2, 2];
            ////vertices[0, 3, 0] = tetra.Vertices[3, 0];
            ////vertices[0, 3, 1] = tetra.Vertices[3, 1];
            ////vertices[0, 3, 2] = tetra.Vertices[3, 2];

            //vertices[0, 0, 0] = 0.0;
            //vertices[0, 0, 1] = 0.0;
            //vertices[0, 0, 2] = 0.0;
            //vertices[0, 1, 0] = 1.0;
            //vertices[0, 1, 1] = 0.0;
            //vertices[0, 1, 2] = 0.0;
            //vertices[0, 2, 0] = 0.0;
            //vertices[0, 2, 1] = 1.0;
            //vertices[0, 2, 2] = 0.0;
            //vertices[0, 3, 0] = 0.0;
            //vertices[0, 3, 1] = 0.0;
            //vertices[0, 3, 2] = 1.0;

            //long[,] cellNeighbours = new long[,] {
            //    { -1, -1, -1, -1 }
            //};
            //byte[,] edgeTags = new byte[,] {
            //    { 0, 0, 0, 0 }
            //};

            //return new UnstructuredTetraGrid(globalId, vertices, cellNeighbours, edgeTags);
        }

        public override int LevelSetDegree {
            get {
                return 1;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new ArgumentException();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                //output[i] = input[i, 2] - Math.Sqrt(2.0) / 3.0;
                //output[i] = input[i, 2] - 2.0 * Math.Sqrt(2.0) / 3.0;
                //output[i] = input[i, 2] - 4.0 * Math.Sqrt(2.0) / 3.0;
                //output[i] = input[i, 2] - 5.0 * Math.Sqrt(2.0) / 6.0;
                output[i] = input[i, 0] - 0.5;
            }
        }
    }

    class SingleTetraSurfaceTestCase : SingleTetraTestCase, ISurfaceTestCase {

        public SingleTetraSurfaceTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                //return Math.Sqrt(3.0) / 3.0;
                //return Math.Sqrt(3.0) / 3.0 / 4.0;
                //return 4.0 * Math.Sqrt(3.0) / 3.0;
                return 1.0 / 8.0;
            }
        }
    }

    class SingleTetraVolumeTestCase : SingleTetraTestCase, IVolumeTestCase {

        public SingleTetraVolumeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                throw new NotImplementedException();
            }
        }
    }
}