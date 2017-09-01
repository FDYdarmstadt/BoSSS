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
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;

namespace CutCellQuadrature.TestCases {

    class Smereka1SquareArcLength : Generic2DTestCase, ISurfaceTestCase {

        private static double L = Math.Sqrt(2.0);

        public Smereka1SquareArcLength(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D() {
                    OffsetX = 0.0,
                    OffsetY = 0.0
                };
            }
        }

        public override double Solution {
            get {
                return 4.0 * L;
            }
        }

        public override int LevelSetDegree {
            get {
                return 1;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double offset = 0.0;

                double x = input[i, 0] - offset;
                double y = input[i, 1] - offset;

                output[i] = (Math.Abs(x) + Math.Abs(y) - 1.0) / Math.Sqrt(2.0);
            }
        }

        public override double GridSpacing {
            get {
                return 0.2;
            }
        }
    }


    class Smereka2EllipseArcLength : Generic2DTestCase, ISurfaceTestCase {

        private const double xMajor = 3.0;

        private const double yMajor = 1.5;

        public Smereka2EllipseArcLength(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 7.26633616541076;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            //return new AnalyticEllipseLevelSet(context);
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            //AnalyticEllipseLevelSet analyticeLevelSet = levelSet as AnalyticEllipseLevelSet;
            //if (analyticeLevelSet == null) {
            //    throw new Exception();
            //}

            //analyticeLevelSet.SetParameters(3.0, 1.5, CurrentShift.OffsetX, CurrentShift.OffsetY);
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            double aSquare = xMajor * xMajor / 4.0;
            double bSquare = yMajor * yMajor / 4.0;

            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                //output[i] = 1.0 - Math.Sqrt(x * x / aSquare + y * y / bSquare);
                output[i] = 1.0 - (x * x / aSquare + y * y / bSquare);
            }
        }

        public override double GridSpacing {
            get {
                return 0.4;
            }
        }
    }
    
    class Smereka3CircleQuadraticIntegrand : Generic2DTestCase, ISurfaceTestCase {

        public Smereka3CircleQuadraticIntegrand(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 2.0 * Math.PI;
            }
        }

        public override int LevelSetDegree {
            get {
                //return 10;
                //return 6;
                //return 4;
                return 2;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
            //return new AnalyticCircleLevelSet(context);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);

            //AnalyticCircleLevelSet analyticeLevelSet = levelSet as AnalyticCircleLevelSet;
            //if (analyticeLevelSet == null) {
            //    throw new Exception();
            //}

            //analyticeLevelSet.SetParameters(1.0, CurrentShift.OffsetX, CurrentShift.OffsetY);
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                //output[i] = 1.0 - Math.Sqrt(x * x + y * y);
                output[i] = 1.0 - (x * x + y * y);
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
                double y = input[i, 1] - CurrentShift.OffsetY;
                output[i] = 3.0 * x * x - y * y;
            }
        }
    }

    class Smereka4EllipsoidSurface : Generic3DTestCase, ISurfaceTestCase {

        private const double xMajor = 3.0;

        private const double yMajor = 1.5;

        private const double zMajor = 1.0;

        public Smereka4EllipsoidSurface(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 9.90182151329315;
                //return 9.90182151329315 / 8.0;
            }
        }

        public override int LevelSetDegree {
            get {
                //return 5;
                return 2;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
            //return new AnalyticEllipsoidLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            ((LevelSet)levelSet).ProjectField(LevelSetInitialValue);
            //AnalyticEllipsoidLevelSet analyticeLevelSet = levelSet as AnalyticEllipsoidLevelSet;
            //if (analyticeLevelSet == null) {
            //    throw new Exception();
            //}

            //analyticeLevelSet.SetOffset(CurrentShift.OffsetX, CurrentShift.OffsetY, CurrentShift.OffsetZ);
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            double aSquare = xMajor * xMajor / 4.0;
            double bSquare = yMajor * yMajor / 4.0;
            double cSquare = zMajor * zMajor / 4.0;

            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                //output[i] = 1.0 - Math.Sqrt(x * x / aSquare + y * y / bSquare + z * z / cSquare);
                output[i] = 1.0 - (x * x / aSquare + y * y / bSquare + z * z / cSquare);
            }
        }

        //public override GridCommons GetGrid(Context context) {
        //    GridCommons grid;
        //    switch (GridType) {
        //        case GridTypes.Structured:
        //            double[] xNodes;
        //            double[] yNodes;
        //            double[] zNodes;

        //            switch (GridSize) {
        //                case GridSizes.Tiny:
        //                    GenericBlas.Linspace(-2.0, 2.0, 6);
        //                    GenericBlas.Linspace(-1.2, 1.2, 4);
        //                    GenericBlas.Linspace(-0.8, 0.8, 3);
        //                    break;

        //                case GridSizes.Small:
        //                    GenericBlas.Linspace(-2.0, 2.0, 11);
        //                    GenericBlas.Linspace(-1.2, 1.2, 7);
        //                    GenericBlas.Linspace(-0.8, 0.8, 5);
        //                    break;

        //                case GridSizes.Normal:
        //                    GenericBlas.Linspace(-2.0, 2.0, 21);
        //                    GenericBlas.Linspace(-1.2, 1.2, 13);
        //                    GenericBlas.Linspace(-0.8, 0.8, 9);
        //                    break;

        //                case GridSizes.Large:
        //                    GenericBlas.Linspace(-2.0, 2.0, 21);
        //                    GenericBlas.Linspace(-1.2, 1.2, 13);
        //                    GenericBlas.Linspace(-0.8, 0.8, 9);
        //                    break;

        //                case GridSizes.Huge:
        //                    GenericBlas.Linspace(-2.0, 2.0, 21);
        //                    GenericBlas.Linspace(-1.2, 1.2, 13);
        //                    GenericBlas.Linspace(-0.8, 0.8, 9);
        //                    break;

        //                default:
        //                    throw new Exception();
        //            }
        //            grid = new Cartesian3DGrid(context.CommMaster, xNodes, yNodes, zNodes, null, null, false, false, false);
        //            break;

        //        case GridTypes.Unstructured:
        //            grid = context.IOMaster.LoadGrid(new Guid("8d36703a-58ca-41a3-96a5-2859357bc60b"));
        //            break;

        //        default:
        //            throw new Exception();
        //    }

        //    return grid;
        //}

        public override double GridSpacing {
            get {
                return 0.2;
                //return 0.4;
            }
        }
    }

    class Smereka5SphereQuadraticIntegrand : Sphere3DTestCase, ISurfaceTestCase {

        public Smereka5SphereQuadraticIntegrand(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 40.0 * Math.PI / 3.0;
            }
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

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
            //return new AnalyticSphereLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            //AnalyticSphereLevelSet analyticLevelSet = levelSet as AnalyticSphereLevelSet;
            //if (analyticLevelSet == null) {
            //    throw new Exception();
            //}

            //analyticLevelSet.SetOffset(CurrentShift.OffsetX, CurrentShift.OffsetY, CurrentShift.OffsetZ);

            LevelSet field = levelSet as LevelSet;
            if (field == null) {
                throw new Exception();
            }

            field.ProjectField(LevelSetInitialValue);
        }

        public override int IntegrandDegree {
            get {
                return 2;
            }
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;
                output[i] = 4.0 - 3.0 * x * x + 2.0 * y * y - z * z;
            }
        }

        public override double GridSpacing {
            get {
                return 0.2;
            }
        }
    }

    class MinGibou1EllipseArea : Generic2DTestCase, IVolumeTestCase {

        private const double xMajor = 3.0;

        private const double yMajor = 1.5;

        public MinGibou1EllipseArea(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return xMajor * yMajor * Math.PI / 4.0;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            //return new AnalyticEllipseLevelSet(gridData);
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            //AnalyticEllipseLevelSet analyticeLevelSet = levelSet as AnalyticEllipseLevelSet;
            //if (analyticeLevelSet == null) {
            //    throw new Exception();
            //}

            //analyticeLevelSet.SetParameters(3.0, 1.5, CurrentShift.OffsetX, CurrentShift.OffsetY);

            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            double aSquare = xMajor * xMajor / 4.0;
            double bSquare = yMajor * yMajor / 4.0;

            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = 1.0 - (x * x / aSquare + y * y / bSquare);
            }
        }

        public override double GridSpacing {
            get {
                return 0.4;
            }
        }
    }

    class MinGibou2EllipsoidVolume : Generic3DTestCase, IVolumeTestCase {

        private const double xMajor = 3.0;

        private const double yMajor = 1.5;

        private const double zMajor = 1.0;

        public MinGibou2EllipsoidVolume(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 4.0 / 3.0 * Math.PI * xMajor * yMajor * zMajor / 8.0;
            }
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            ((LevelSet)levelSet).ProjectField(LevelSetInitialValue);
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            double aSquare = xMajor * xMajor / 4.0;
            double bSquare = yMajor * yMajor / 4.0;
            double cSquare = zMajor * zMajor / 4.0;

            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                double z = input[i, 2] - CurrentShift.OffsetZ;

                output[i] = 1.0 - (x * x / aSquare + y * y / bSquare + z * z / cSquare);
            }
        }

        public override double GridSpacing {
            get {
                return 0.2;
            }
        }
    }

    abstract class StarTestCase : Generic2DTestCase {

        public StarTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
            //return new AnalyticStarLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            //AnalyticStarLevelSet analyticLevelSet = levelSet as AnalyticStarLevelSet;
            //if (analyticLevelSet == null) {
            //    throw new Exception();
            //}
            //analyticLevelSet.SetParameters(CurrentShift.OffsetX, CurrentShift.OffsetY);

            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        public override int LevelSetDegree {
            get {
                return 7;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                //double norm = (x * x + y * y).Sqrt();
                //if (norm < 1e-4) {
                //    output[i] = -1.0;
                //} else {
                //    output[i] = norm - 0.5 -
                //        (Math.Pow(y, 5.0) + 5.0 * Math.Pow(x, 4.0) * y - 10.0 * Math.Pow(x, 2.0) * Math.Pow(y, 3.0)) /
                //        (3.0 * Math.Pow(norm, 4.0));
                //}

                //output[i] = 2.0 * normSquare * normSquare * normSquare - 2.0 * Math.Pow(normSquare, 2.5) -
                //    (Math.Pow(x, 5.0) + 5.0 * x * Math.Pow(y, 4.0) - 10.0 * x * x * x * y * y);

                double r = Math.Sqrt(x * x + y * y);
                double theta = Math.Atan2(y, x);
                double value = 1.0 + 0.5 * Math.Cos(5.0 * theta) - r;
                //if (value > 0.9 || r < 0.05) {
                //    value = 0.9;
                //}
                //output[i] = value.Sign() * Math.Min(value.Abs(), 1.0);

                output[i] = value;
            }
        }

        public override GridCommons GetGrid(IDatabaseInfo db) {
            GridCommons grid;
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

                    grid = Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-2.0, 2.0, noOfCellsPerDirection + 1),
                        GenericBlas.Linspace(-2.0, 2.0, noOfCellsPerDirection + 1),
                        CellType.Square_Linear,
                        false,
                        false,
                        new BoundingBox(new double[,] { {-0.2, -0.2}, { 0.2, 0.2 } }));
                    break;

                case GridTypes.Unstructured:
                    grid = db.Controller.DBDriver.LoadGrid(
                        new Guid("7a4cf525-76e0-44fc-add2-7ce683a082c3"), db);
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
    }

    class MinGibou3StarArcLength : StarTestCase, ISurfaceTestCase {

        public MinGibou3StarArcLength(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 12.329044714372;
            }
        }
    }

    class MinGibou4StarArea : StarTestCase, IVolumeTestCase {

        public MinGibou4StarArea(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 9.0 / 8.0 * Math.PI;
            }
        }
    }

    class ChengFriesSinusoidalArea : Generic2DTestCase, IVolumeTestCase {

        public ChengFriesSinusoidalArea(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 5;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];
                output[i] = Math.Cos(2.0 * x) / 3.0 - y + 0.4;
            }
        }

        public override double Solution {
            get {
                return 2.0 + 2.0 * Math.Sin(2.0) / 6.0 + 0.8;
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override GridCommons GetGrid(IDatabaseInfo db) {
            GridCommons grid;
            switch (GridType) {
                case GridTypes.Structured:
                    int noOfCellsPerDirection;
                    switch (GridSize) {
                        case GridSizes.Tiny:
                            noOfCellsPerDirection = 1;
                            break;

                        case GridSizes.Small:
                            noOfCellsPerDirection = 2;
                            break;

                        case GridSizes.Normal:
                            noOfCellsPerDirection = 4;
                            break;

                        case GridSizes.Large:
                            noOfCellsPerDirection = 8;
                            break;

                        case GridSizes.Huge:
                            noOfCellsPerDirection = 16;
                            break;

                        default:
                            throw new Exception();
                    }
                    double[] nodes = GenericBlas.Linspace(-1.0, 1.0, noOfCellsPerDirection + 1);
                    grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                    break;

                case GridTypes.Unstructured:
                    grid = null;
                    switch (GridSize) {
                        case GridSizes.Tiny:
                            //grid = context.IOMaster.LoadGrid(new Guid("90bfc279-6b6a-45c4-914e-3dfddbcbff45"));
                            break;

                        case GridSizes.Small:
                            //grid = context.IOMaster.LoadGrid(new Guid("f0ba754c-49a5-45dd-8ab8-b49a514efe6d"));
                            break;

                        case GridSizes.Normal:
                            //grid = context.IOMaster.LoadGrid(new Guid("943c3b1a-f254-4904-be9e-2124e2322739"));
                            break;

                        case GridSizes.Large:
                            //grid = context.IOMaster.LoadGrid(new Guid("db230415-78b4-4fd2-bc81-1be4b89dfa75"));
                            break;

                        case GridSizes.Huge:
                            //grid = context.IOMaster.LoadGrid(new Guid("2b6d74ed-5b91-4328-9449-5c1a9d12ba93"));
                            break;

                        default:
                            throw new Exception();
                    }
                    break;

                default:
                    throw new Exception();
            }

            return grid;
        }
    }

    class ChengFriesSinusoidalLength : Generic2DTestCase, ISurfaceTestCase {

        public ChengFriesSinusoidalLength(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 10;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];
                output[i] = Math.Cos(2.0 * x) / 3.0 - y + 0.4;
            }
        }

        public override double Solution {
            get {
                return 2.0 * 1.1220338643772864;
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override GridCommons GetGrid(IDatabaseInfo db) {
            GridCommons grid;
            switch (GridType) {
                case GridTypes.Structured:
                    int noOfCellsPerDirection;
                    switch (GridSize) {
                        case GridSizes.Tiny:
                            noOfCellsPerDirection = 1;
                            break;

                        case GridSizes.Small:
                            noOfCellsPerDirection = 2;
                            break;

                        case GridSizes.Normal:
                            noOfCellsPerDirection = 4;
                            break;

                        case GridSizes.Large:
                            noOfCellsPerDirection = 8;
                            break;

                        case GridSizes.Huge:
                            noOfCellsPerDirection = 16;
                            break;

                        default:
                            throw new Exception();
                    }
                    double[] nodes = GenericBlas.Linspace(-1.0, 1.0, noOfCellsPerDirection + 1);
                    grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                    break;

                case GridTypes.Unstructured:
                    grid = null;
                    switch (GridSize) {
                        case GridSizes.Tiny:
                            //grid = context.IOMaster.LoadGrid(new Guid("90bfc279-6b6a-45c4-914e-3dfddbcbff45"));
                            break;

                        case GridSizes.Small:
                            //grid = context.IOMaster.LoadGrid(new Guid("f0ba754c-49a5-45dd-8ab8-b49a514efe6d"));
                            break;

                        case GridSizes.Normal:
                            //grid = context.IOMaster.LoadGrid(new Guid("943c3b1a-f254-4904-be9e-2124e2322739"));
                            break;

                        case GridSizes.Large:
                            //grid = context.IOMaster.LoadGrid(new Guid("db230415-78b4-4fd2-bc81-1be4b89dfa75"));
                            break;

                        case GridSizes.Huge:
                            //grid = context.IOMaster.LoadGrid(new Guid("2b6d74ed-5b91-4328-9449-5c1a9d12ba93"));
                            break;

                        default:
                            throw new Exception();
                    }
                    break;

                default:
                    throw new Exception();
            }

            return grid;
        }
    }
}
