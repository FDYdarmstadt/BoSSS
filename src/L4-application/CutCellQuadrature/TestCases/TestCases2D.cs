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
using System.Linq;
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

    /// <summary>
    /// Generic base class for all two-dimensional test cases.
    /// </summary>
    abstract class Generic2DTestCase : TestCase<Shift2D> {

        /// <summary>
        /// Cache for <see cref="AllShifts"/>
        /// </summary>
        private IEnumerable<Shift2D> shifts = null;

        /// <summary>
        /// Constructs a new test case
        /// </summary>
        /// <param name="gridSize">
        /// <see cref="ITestCase.GridSize"/>
        /// </param>
        /// <param name="gridType">
        /// <see cref="ITestCase.GridType"/>
        /// </param>
        public Generic2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        /// <summary>
        /// Generates a list of 100 random shifts.
        /// </summary>
        protected override IEnumerable<Shift2D> AllShifts {
            get {
                if (shifts == null) {
                    List<Shift2D> myShifts = new List<Shift2D>();
                    Random rng = new Random(1);
                    for (int i = 1; i < 101; i++) {
                        double x = rng.NextDouble();
                        double y = rng.NextDouble();

                        myShifts.Add(new Shift2D() {
                            OffsetX = x,
                            OffsetY = y
                        });
                    }
                    shifts = myShifts;
                }

                return shifts;
            }
        }

        /// <summary>
        /// Constructs an appropriate depending on the current settings.
        /// </summary>
        /// <param name="db">
        /// The (optional) db the grid may be stored in
        /// </param>
        /// <returns>
        /// A grid for the standard domain [-2,2]x[-2,2]
        /// </returns>
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
                    grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                    break;

                case GridTypes.PseudoStructured:
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
                    nodes = GenericBlas.Linspace(-2.0, 2.0, noOfCellsPerDirection + 1);
                    grid = Grid2D.UnstructuredTriangleGrid(nodes, nodes);
                    break;

                case GridTypes.Unstructured:
                    switch (GridSize) {
                        case GridSizes.Tiny:
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("90bfc279-6b6a-45c4-914e-3dfddbcbff45"), db);
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("8c374ad2-d196-4447-9dda-e56557b8a747"), db);
                            break;

                        case GridSizes.Small:
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("f0ba754c-49a5-45dd-8ab8-b49a514efe6d"), db);
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("5d2b78ff-238a-41d4-a4e8-76a25380f4a4"), db);
                            break;

                        case GridSizes.Normal:
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("943c3b1a-f254-4904-be9e-2124e2322739"), db);
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("12c57f4d-f9e2-4f0b-bc7a-2ac0499341c9"), db);
                            break;

                        case GridSizes.Large:
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("db230415-78b4-4fd2-bc81-1be4b89dfa75"), db);
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("6a771d0e-c351-4437-965e-abbf43b15497"), db);
                            break;

                        case GridSizes.Huge:
                            //grid = db.DBController.DBDriver.LoadGrid(
                            //    new Guid("2b6d74ed-5b91-4328-9449-5c1a9d12ba93"), db);
                            grid = db.Controller.DBDriver.LoadGrid(
                                new Guid("0dbffb53-2bb1-4ec8-95bb-8aec417c78d4"), db);
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

        /// <summary>
        /// The reference grid spacing.
        /// </summary>
        public override double GridSpacing {
            get {
                return 0.4;
            }
        }

        /// <summary>
        /// Returns an empty level set field, <see cref="LevelSet"/>.
        /// </summary>
        /// <param name="gridData">
        /// <see cref="ITestCase.GetLevelSet"/>
        /// </param>
        /// <returns><see cref="ITestCase.GetLevelSet"/></returns>
        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        /// <summary>
        /// <see cref="ITestCase.UpdateLevelSet"/>
        /// </summary>
        /// <param name="levelSet">
        /// <see cref="ITestCase.UpdateLevelSet"/>
        /// </param>
        public override void UpdateLevelSet(ILevelSet levelSet) {
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        /// <summary>
        /// The degree of the level set function defined by
        /// <see cref="LevelSetInitialValue"/>
        /// </summary>
        public abstract int LevelSetDegree {
            get;
        }

        /// <summary>
        /// A <see cref="ScalarFunction"/> defining the level set field and
        /// thus the domain of integration
        /// </summary>
        /// <param name="input"><see cref="ScalarFunction"/></param>
        /// <param name="output"><see cref="ScalarFunction"/></param>
        public abstract void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output);
    }

    abstract class Circle2DTestCase : Generic2DTestCase {

        public Circle2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            //return new AnalyticCircleLevelSet(gridData);
            return new LevelSet(new Basis(gridData, 2), "levelSet");
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                //yield return new Shift2D() {
                //    OffsetX = 0.049733716831418556,
                //    OffsetY = 0.022148795436205713
                //};
                //yield return new Shift2D() {
                //    OffsetX = 0.0,
                //    OffsetY = 0.1
                //};
                yield return new Shift2D();
            }
        }

        //public override void UpdateLevelSet(ILevelSet levelSet) {
        //    AnalyticCircleLevelSet analyticLevelSet =
        //        levelSet as AnalyticCircleLevelSet;
        //    if (analyticLevelSet == null) {
        //        throw new ArgumentException();
        //    }

        //    analyticLevelSet.SetParameters(1.0, CurrentShift.OffsetX, CurrentShift.OffsetY);
        //    //((LevelSet)levelSet).ProjectField(LevelSetInitialValue);
        //}

        public override int LevelSetDegree {
            get {
                //return 10;
                return 2;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                //output[i] = 1.0 - Math.Sqrt(x * x + y * y);
                output[i] = 1.0 - (x * x + y * y);
            }
        }
    }

    class QuarterCircleSignedDistance :
        Generic2DTestCase,
        IVolumeTestCase
    {

        public QuarterCircleSignedDistance(GridSizes gridSize, GridTypes gridType) :
                base(gridSize, gridType)
        {
            
        }

        public override double Solution {
            get {
                return Math.PI / 4.0;
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override int LevelSetDegree {
            get {
                return 10;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output)
        {
            for (int i = 0; i < output.GetLength(0); i++)
            {
                double x = input[i, 0] + 1;
                double y = input[i, 1] - 1;
                output[i] = 1.0 - Math.Sqrt(x * x + y * y);
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db)
        {
            switch (GridType)
            {
                case GridTypes.Structured:
                    int noOfCells = (int)this.GridSize + 1;
                    return Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-1.0, 1.0, noOfCells + 1),
                        GenericBlas.Linspace(-1.0, 1.0, noOfCells + 1));

                case GridTypes.Unstructured:
                    throw new NotImplementedException();

                default:
                    throw new Exception();
            }
        }
    }


    class CircleVolume2DTestCase : Circle2DTestCase, IVolumeTestCase {

        public CircleVolume2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override double Solution {
            get {
                return Math.PI;
            }
        }
    }

    class ConeVolume2DTestCase : Circle2DTestCase, IVolumeTestCase {

        public ConeVolume2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return Math.PI / 3.0;
            }
        }

        public override int IntegrandDegree {
            get {
                return 5;
            }
        }

        public override void JumpingFieldSpeciesBInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = 1.0 - Math.Sqrt(x * x + y * y);
            }
        }
    }

    class InverseConeVolume2DTestCase : Circle2DTestCase, IVolumeTestCase {

        public InverseConeVolume2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 0.5 * Math.PI;
                //return 2.0 * Math.PI / 3.0;
            }
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = 1.0 - (x * x + y * y);
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        public override int IntegrandDegree {
            get {
                return 2;
            }
        }

        public override void JumpingFieldSpeciesBInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                //output[i] = Math.Sqrt(x * x + y * y);
                output[i] = x * x + y * y;
            }
        }
    }

    class BowlVolume2DTestCase : Circle2DTestCase, IVolumeTestCase {

        public BowlVolume2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return Math.PI / 2.0;
            }
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = 1.0 - (x * x + y * y);
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");
        }

        public override int IntegrandDegree {
            get {
                return 2;
            }
        }

        public override void JumpingFieldSpeciesBInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = x * x + y * y;
            }
        }

        public override void JumpingFieldSpeciesAInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                output[i] = 1.0;
            }
        }
    }

    abstract class Ellipse2DTestCase : Generic2DTestCase {

        public Ellipse2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public abstract double xMajor {
            get;
        }

        public abstract double yMajor {
            get;
        }

        public override int LevelSetDegree {
            get {
                //return 5;
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
    }

    class EllipseVolume2DTestCase : Ellipse2DTestCase, IVolumeTestCase {

        public EllipseVolume2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double xMajor {
            get {
                return 3.0;
            }
        }

        public override double yMajor {
            get {
                return 1.5;
            }
        }

        public override double Solution {
            get {
                return Math.PI * xMajor * yMajor / 4.0;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new AnalyticEllipseLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            AnalyticEllipseLevelSet analyticeLevelSet = levelSet as AnalyticEllipseLevelSet;
            if (analyticeLevelSet == null) {
                throw new Exception();
            }

            analyticeLevelSet.SetParameters(3.0, 1.5, CurrentShift.OffsetX, CurrentShift.OffsetY);
        }
    }

    class SingleSquareStraightLineLengthTestCase : Generic2DTestCase, ISurfaceTestCase {

        public SingleSquareStraightLineLengthTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            if (GridType != GridTypes.Structured) {
                throw new NotImplementedException();
            }

            return Grid2D.Cartesian2DGrid(
                GenericBlas.Linspace(-1.0, 1.0, 2), GenericBlas.Linspace(-1.0, 1.0, 2));
        }

        public override double Solution {
            get {
                return 2.0;
            }
        }

        public override int LevelSetDegree {
            get {
                return 1;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                output[i] = input[i, 1];
            }
        }

        public override int IntegrandDegree {
            get {
                return 0;
            }
        }
    }

    class SingleSquareStraightLineVolumeTestCase : Generic2DTestCase, IVolumeTestCase {

        public SingleSquareStraightLineVolumeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            if (GridType != GridTypes.Structured) {
                throw new NotImplementedException();
            }

            return Grid2D.Cartesian2DGrid(
                GenericBlas.Linspace(-1.0, 1.0, 2), GenericBlas.Linspace(-1.0, 1.0, 2));
            //return Grid2D.Cartesian2DGrid(
            //    GenericBlas.Linspace(-1.0, 1.0, 2), GenericBlas.Linspace(-1.0, 1.0, 2));
            //return Grid2D.Cartesian2DGrid(
            //    GenericBlas.Linspace(-2.0, 2.0, 2), GenericBlas.Linspace(-2.0, 2.0, 2));
        }

        public override double Solution {
            get {
                //return 3.5;
                return 8.0;
            }
        }

        public override int LevelSetDegree {
            get {
                return 1;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                output[i] = input[i, 0] + input[i, 1] - 1.0;
                //output[i] = input[i, 0] - 0.5;
                //output[i] = input[i, 0] + input[i, 1];
            }
        }

        public override int IntegrandDegree {
            get {
                return 0;
            }
        }
    }

    abstract class SingleTriangleTestCase : Generic2DTestCase {

        public SingleTriangleTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            if (GridType != GridTypes.Unstructured) {
                throw new NotImplementedException();
            }

            throw new NotImplementedException("Not ported yet");
            //return Grid2D.UnstructuredTriangleGrid(
            //    GenericBlas.Linspace(0.0, 1.0, 2),
            //    GenericBlas.Linspace(0.0, 1.0, 2));
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                output[i] = input[i, 0] - 1e-13;
                //output[i] = input[i, 0] - 0.5;
                //output[i] = input[i, 0] - 0.75;
                //output[i] = input[i, 0] - 1.0 / Math.Sqrt(3.0);
                //output[i] = input[i, 0] - 2.0 / Math.Sqrt(3.0);
            }
        }

        public override int IntegrandDegree {
            get {
                return 1;
            }
        }
    }

    class SingleTriangleSurfaceTestCase : SingleTriangleTestCase, ISurfaceTestCase {

        public SingleTriangleSurfaceTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 0.5;
                //return 1.0;
                //return 2.0;
                //return 4.0;
            }
        }
    }

    class SingleTriangleVolumeTestCase : SingleTriangleTestCase, IVolumeTestCase {

        public SingleTriangleVolumeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 2.0 / 3.0 * Math.Sqrt(3.0);
                //return 2.0 / 3.0 * Math.Sqrt(3.0) * 4.0;
                //return 0.5;
                //return 1.0 / 8.0;
            }
        }
    }

    class CircleArcLength : Circle2DTestCase, ISurfaceTestCase {

        public CircleArcLength(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 2.0 * Math.PI;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, 2), "levelSet");
        }
    }

    class LinearIntgreandCircleSurfaceIntegral2DTestCase : Circle2DTestCase, ISurfaceTestCase {

        public LinearIntgreandCircleSurfaceIntegral2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 2.0 * Math.PI;
            }
        }

        public override int IntegrandDegree {
            get {
                return 1;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new AnalyticCircleLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            AnalyticCircleLevelSet analyticeLevelSet = levelSet as AnalyticCircleLevelSet;
            if (analyticeLevelSet == null) {
                throw new Exception();
            }

            analyticeLevelSet.SetParameters(1.0, CurrentShift.OffsetX, CurrentShift.OffsetY);
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                // We want to integrate over x + 1 but the circle has an offset
                double x = input[i, 0] - CurrentShift.OffsetX;
                output[i] = x + 1.0;
            }
        }
    }

    class QuadraticIntegrandCircleSurfaceIntegral2DTestCase : Circle2DTestCase, ISurfaceTestCase {

        public QuadraticIntegrandCircleSurfaceIntegral2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return Math.PI;
            }
        }

        public override int IntegrandDegree {
            get {
                return 2;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new AnalyticCircleLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            AnalyticCircleLevelSet analyticeLevelSet = levelSet as AnalyticCircleLevelSet;
            if (analyticeLevelSet == null) {
                throw new Exception();
            }

            analyticeLevelSet.SetParameters(1.0, CurrentShift.OffsetX, CurrentShift.OffsetY);
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                // We want to integrate over x^2 but the circle has an offset
                double x = input[i, 0] - CurrentShift.OffsetX;
                output[i] = x * x;
            }
        }
    }

    class ExponentialIntegrandCircleSurfaceIntegral2DTestCase : Circle2DTestCase, ISurfaceTestCase {

        public ExponentialIntegrandCircleSurfaceIntegral2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                // 2 * Pi * I(0, 1) where I is the Bessel function of the first kind
                return 7.9549265210128453;
            }
        }

        public override int IntegrandDegree {
            get {
                return 10;
            }
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                output[i] = Math.Exp(x);
            }
        }
    }


    class SquareVolumeTestCase : Generic2DTestCase, IVolumeTestCase {

        public SquareVolumeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
        }

        public override double Solution {
            get {
                return 1.0;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new AnalyticSquareLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            AnalyticSquareLevelSet analyticeLevelSet = levelSet as AnalyticSquareLevelSet;
            if (analyticeLevelSet == null) {
                throw new Exception();
            }

            analyticeLevelSet.SetOffset(CurrentShift.OffsetX, CurrentShift.OffsetY);
        }
    }


    class SquareArcLengthTestCase : Generic2DTestCase, ISurfaceTestCase {

        public SquareArcLengthTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 10;
            }
        }

        private static double L2Distance(double x1, double y1, double x2, double y2) {
            return Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        }

        private static double L1Norm(double x, double y) {
            return Math.Abs(x) + Math.Abs(y);
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            double c = Math.Sqrt(2.0) / 2.0;
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];

                if (y > x + c && y > -x + c) {
                    output[i] = -L2Distance(x, y, 0, c);
                } else if (y > x + c && y < -x - c) {
                    output[i] = -L2Distance(x, y, -c, 0);
                } else if (y < -x - c && y < x - c) {
                    output[i] = -L2Distance(x, y, 0, -c);
                } else if (y < x - c && y > -x + c) {
                    output[i] = -L2Distance(x, y, c, 0);
                } else {
                    output[i] = (c - L1Norm(x, y)) * c;
                }
            }
        }

        public override double Solution {
            get {
                return 4.0;
            }
        }

        public override ILevelSet GetLevelSet(GridData gridData) {
            return new LevelSet(new Basis(gridData, LevelSetDegree), "levelSet");

            //return new AnalyticSquareLevelSet(gridData);
        }

        public override void UpdateLevelSet(ILevelSet levelSet) {
            //AnalyticSquareLevelSet analyticeLevelSet = levelSet as AnalyticSquareLevelSet;
            //if (analyticeLevelSet == null) {
            //    throw new Exception();
            //}

            //analyticeLevelSet.SetOffset(CurrentShift.OffsetX, CurrentShift.OffsetY);
            LevelSet levelSetField = levelSet as LevelSet;
            if (levelSetField == null) {
                throw new Exception();
            }

            levelSetField.ProjectField(LevelSetInitialValue);
        }
    }

    abstract class ParabolaSingleCell2DTestCase : Generic2DTestCase {

        public ParabolaSingleCell2DTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
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
                //output[i] = 0.5 * x * x + x - 0.5 - y + 1e-13;
                output[i] = -x * x - 2.0 * x + 0.5 * y - 0.5;
                //output[i] = -x * x - 2.0 * x + y - 1.0;
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            switch (GridType) {
                case GridTypes.Structured:
                    int noOfCells = 2 << (int)this.GridSize;
                    return Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-1.0, 1.0, noOfCells),
                        GenericBlas.Linspace(-1.0, 1.0, noOfCells));

                case GridTypes.Unstructured:
                    throw new NotImplementedException();

                default:
                    throw new Exception();
            }
        }
    }


    class SingleSquareParabolaVolumeTestCase : ParabolaSingleCell2DTestCase, IVolumeTestCase {

        public SingleSquareParabolaVolumeTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                //return 2.0;
                return 4.0 / 3.0;
            }
        }
    }

    class SingleSquareParabolaLengthTestCase : ParabolaSingleCell2DTestCase, ISurfaceTestCase {

        public SingleSquareParabolaLengthTestCase(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 2.9578857150891949;
            }
        }
    }




    class FlorianSurface : Generic2DTestCase, ISurfaceTestCase {

        public FlorianSurface()
            : base(GridSizes.Tiny, GridTypes.Structured) {
        }

        public override double Solution {
            get {
                return Math.Sqrt(2.0) * Math.PI;
            }
        }

        public override int LevelSetDegree {
            get {
                return 10;
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            return Grid2D.Cartesian2DGrid(
                GenericBlas.Linspace(-1.0, 1.0, 10),
                GenericBlas.Linspace(-1.0, 1.0, 10));
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

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                //output[i] = 1.0 / 2.0 - (x * x + y * y);
                output[i] = Math.Sqrt(2.0) / 2.0 - Math.Sqrt(x * x + y * y);
            }
        }

        public override int IntegrandDegree {
            get {
                return 0;
            }
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;
                output[i] = 1.0;
            }
        }
    }

    class FlorianVolume : Generic2DTestCase, IVolumeTestCase {

        public FlorianVolume()
            : base(GridSizes.Tiny, GridTypes.Structured) {
        }

        public override double Solution {
            get {
                return 0.5 * Math.PI;
            }
        }

        public override int LevelSetDegree {
            get {
                return 10;
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            return Grid2D.Cartesian2DGrid(
                GenericBlas.Linspace(-1.0, 1.0, 10),
                GenericBlas.Linspace(-1.0, 1.0, 10));
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

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                //output[i] = 1.0 / 2.0 - (x * x + y * y);
                output[i] = Math.Sqrt(2.0) / 2.0 - Math.Sqrt(x * x + y * y);
            }
        }

        public override int IntegrandDegree {
            get {
                return 0;
            }
        }
    }

    class DuesterSurface : Generic2DTestCase, ISurfaceTestCase {

        private double xM;

        private double r;

        public DuesterSurface(GridSizes gridSize, double xM)
            : base(gridSize, GridTypes.Structured) {
            this.xM = xM;
            r = Math.Sqrt(Math.Pow(0.5 - xM, 2.0) + 1.0);
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
                output[i] = (x - xM) * (x - xM) + y * y - r * r;
            }
        }

        public override double Solution {
            get {
                return 2.0 * r * Math.Asin(1.0 / r);
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            switch (this.GridSize) {
                case GridSizes.Tiny:
                    return Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-1.0, 1.0, 2),
                        GenericBlas.Linspace(-1.0, 1.0, 2));

                case GridSizes.Small:
                    return Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-1.0, 1.0, 3),
                        GenericBlas.Linspace(-1.0, 1.0, 3));

                case GridSizes.Normal:
                    return Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-1.0, 1.0, 5),
                        GenericBlas.Linspace(-1.0, 1.0, 5));

                case GridSizes.Large:
                    return Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-1.0, 1.0, 9),
                        GenericBlas.Linspace(-1.0, 1.0, 9));

                case GridSizes.Huge:
                    return Grid2D.Cartesian2DGrid(
                        GenericBlas.Linspace(-1.0, 1.0, 17),
                        GenericBlas.Linspace(-1.0, 1.0, 17));

                default:
                    throw new Exception();
            }
        }
    }


    abstract class TestEgger : Generic2DTestCase {

        public readonly double Offset = 0.0;

        private readonly double a = 2.0;

        public TestEgger(GridSizes gridSize)
            : base(gridSize, GridTypes.Structured) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
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

                    double[] xnodes = GenericBlas.Linspace(-0.125, 0.125, noOfCellsPerDirection + 1);
                    double[] ynodes = GenericBlas.Linspace(-0.5, 0.5, noOfCellsPerDirection + 1);
                    grid = Grid2D.Cartesian2DGrid(xnodes, ynodes);
                    break;

                default:
                    throw new NotImplementedException();
            }

            return grid;
        }

        public override int LevelSetDegree {
            get {
                return 3;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                double x = input[i, 0] - Offset;
                double y = input[i, 1];
                output[i] = x * (1.0 + a * x * y);
            }
        }
    }

    class TestEggerSurface : TestEgger, ISurfaceTestCase {

        public TestEggerSurface(GridSizes gridSize)
            : base(gridSize) {
        }

        public override double Solution {
            get {
                return 1.0;
            }
        }
    }

    class TestEggerVolume : TestEgger, IVolumeTestCase {

        public TestEggerVolume(GridSizes gridSize)
            : base(gridSize) {
        }

        public override double Solution {
            get {
                return (0.125 - Offset) * 2.0;
            }
        }
    }


    abstract class KarakusCircle : Generic2DTestCase {


        public KarakusCircle(GridSizes gridSize, GridTypes gridTypes)
            : base(gridSize, gridTypes) {
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = -1 * ((x - 1).Pow2() + (y - 1).Pow2() + 10) * (x * x + y * y - 1);
            }
        }

        public override int LevelSetDegree {
            get {
                return 4;
            }
        }
    }

    class KarakusCircleVolume : KarakusCircle, IVolumeTestCase {
        public KarakusCircleVolume(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return Math.PI;
            }
        }

    }

    class KarakusCircleArcLength : KarakusCircle, ISurfaceTestCase {
        public KarakusCircleArcLength(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        public override double Solution {
            get {
                return 2 * Math.PI;
            }
        }

    }

    class IBMCircleVolume : Generic2DTestCase, IVolumeTestCase {

        public IBMCircleVolume() : base(GridSizes.Normal, GridTypes.Structured) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            int numberOfCells = 32 << 1;
            double stretching = 2.0;
            double D = 40.0;

            double[] nodesA = Grid1D.TanhSpacing(-D, 0, numberOfCells / 2 + 1, stretching, false);
            double[] nodesB = Grid1D.TanhSpacing(0, D, numberOfCells / 2 + 1, stretching, true);
            double[] nodes = nodesA.Take(nodesA.Length - 1).Concat(nodesB).ToArray();
            GridCommons grid = Grid2D.Cartesian2DGrid(nodes, nodesB);

            return grid;
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override double Solution {
            get {
                double D = 40.0;
                return (2.0 * D) * D - Math.PI / 2.0;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = x * x + y * y - 1.0;
            }
        }
    }

    class IBMCircleSurface : Generic2DTestCase, ISurfaceTestCase {

        public IBMCircleSurface() : base(GridSizes.Normal, GridTypes.Structured) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            int numberOfCells = 32 << 1;
            double stretching = 2.0;
            double D = 40.0;

            double[] nodesA = Grid1D.TanhSpacing(-D, 0, numberOfCells / 2 + 1, stretching, false);
            double[] nodesB = Grid1D.TanhSpacing(0, D, numberOfCells / 2 + 1, stretching, true);
            double[] nodes = nodesA.Take(nodesA.Length - 1).Concat(nodesB).ToArray();
            GridCommons grid = Grid2D.Cartesian2DGrid(nodes, nodesB);

            return grid;
        }

        public override int LevelSetDegree {
            get {
                return 2;
            }
        }

        public override double Solution {
            get {
                return Math.PI;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0] - CurrentShift.OffsetX;
                double y = input[i, 1] - CurrentShift.OffsetY;

                output[i] = x * x + y * y - 1.0;
            }
        }
    }

    class NACA0012Base : Generic2DTestCase {

        public NACA0012Base(GridSizes gridSizes)
            : base(gridSizes, GridTypes.Structured) {
        }

        public override int LevelSetDegree {
            get {
                return 8;
            }
        }

        public override double Solution {
            get { // Maple arc length solution
                return 2.039548931156011151948822860131849293000;
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                // OffsetX is used as angle of attack
                List<Shift2D> shifts = new List<Shift2D>();
                for (int i = 0; i < 2; i++) {
                    shifts.Add(new Shift2D() {
                        OffsetX = 2.0 * i,
                        OffsetY = 0.0
                    });
                }

                return shifts;
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {

            int MeshPara;
            switch (GridSize) {
                case GridSizes.Tiny:
                    MeshPara = 4;
                    break;

                case GridSizes.Small:
                    MeshPara = 8;
                    break;

                case GridSizes.Normal:
                    MeshPara = 16;
                    break;

                case GridSizes.Large:
                    MeshPara = 32;
                    break;

                case GridSizes.Huge:
                    MeshPara = 64;
                    break;

                default:
                    throw new Exception();
            }

            double xBegin = -0.025;
            double xEnd = 1.025;
            int chords = 20;

            int xleft = -chords;
            int xRight = chords;
            int yBottom = -chords;
            int yTop = chords;

            double spacingFactor = 3.9;

            int spacingElements = (int)(0.5 * MeshPara) + 1;

            double[] xnodes1 = Grid1D.TanhSpacing(xleft, xBegin, spacingElements, spacingFactor, false);
            double[] xnodes2 = Grid1D.TanhSpacing_DoubleSided(xBegin, xEnd, MeshPara + 1, 2.0);
            double[] xnodes3 = Grid1D.TanhSpacing(xEnd, xRight, spacingElements, spacingFactor, true);
            double[] xComplete = new double[xnodes1.Length + xnodes2.Length + xnodes3.Length - 2];

            for (int i = 0; i < xnodes1.Length; i++) {
                xComplete[i] = xnodes1[i];
            }
            for (int i = 1; i < xnodes2.Length; i++) {
                xComplete[i + xnodes1.Length - 1] = xnodes2[i];
            }
            for (int i = 1; i < xnodes3.Length; i++) {
                xComplete[i + xnodes1.Length + xnodes2.Length - 2] = xnodes3[i];
            }

            double yrefinedTop = 0.2;
            double yrefinedBottom = -0.2;

            double[] ynodes1 = Grid1D.TanhSpacing(yBottom, yrefinedBottom, MeshPara + 1, spacingFactor, false);
            double[] ynodes2 = GenericBlas.Linspace(yrefinedBottom, yrefinedTop, (int)(0.75 * MeshPara) + 1);
            double[] ynodes3 = Grid1D.TanhSpacing(yrefinedTop, yTop, MeshPara + 1, spacingFactor, true);

            double[] yComplete = new double[ynodes1.Length + ynodes2.Length + ynodes3.Length - 2];
            for (int i = 0; i < ynodes1.Length; i++) {
                yComplete[i] = ynodes1[i];
            }
            for (int i = 1; i < ynodes2.Length; i++) {
                yComplete[i + ynodes1.Length - 1] = ynodes2[i];
            }
            for (int i = 1; i < ynodes3.Length; i++) {
                yComplete[i + ynodes1.Length + ynodes2.Length - 2] = ynodes3[i];
            }

            GridCommons grid = Grid2D.Cartesian2DGrid(
                xComplete,
                yComplete
                );

            return grid;
        }


        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {

                double x = input[i, 0];
                double y = input[i, 1];

                double alpha = CurrentShift.OffsetX;

                double radian = alpha * Math.PI / 180;

                double xRotated = 1 + Math.Cos(radian) * (x - 1) - Math.Sin(radian) * y;
                double yRotated = Math.Sin(radian) * (x - 1) + Math.Cos(radian) * y;

                double a = 0.6;
                //double b = 0.2969;
                double c1 = 0.126;
                double d = 0.3516;
                double e = 0.2843;
                double f = 0.1036;


                if ((yRotated >= 0.0 && x > 0.062165) || (yRotated >= 0.0 && x < -0.025) || (x > 0.562875 && y > 0) || (y > 0.03333385 && x > -0.025 && x < 0.062165)) {
                    output[i] = Math.Pow((yRotated + a * (c1 * xRotated + d * Math.Pow(xRotated, 2) - e * Math.Pow(xRotated, 3) + f * Math.Pow(xRotated, 4))), 2) - 0.0317338596 * xRotated;
                } else {
                    output[i] = Math.Pow((-yRotated + a * (c1 * xRotated + d * Math.Pow(xRotated, 2) - e * Math.Pow(xRotated, 3) + f * Math.Pow(xRotated, 4))), 2) - 0.0317338596 * xRotated;
                }
            }

        }
    }

    class NACA0012ArcLength : NACA0012Base, ISurfaceTestCase {
        public NACA0012ArcLength(GridSizes sizes) :
            base(sizes) {
        }

        public override double Solution {
            get { // Maple arc length solution
                return 2.039548931156011151948822860131849293000;
            }
        }
    }

    class NACA0012QudraticIntegrand : NACA0012Base, ISurfaceTestCase {
        public NACA0012QudraticIntegrand(GridSizes sizes) :
            base(sizes) {

        }

        public override double Solution {
            get { // Maple solution for line integral with f(x,y)=x^2
                return 0.6746469715584554876004040497116201630508;
            }
        }

        public override int IntegrandDegree {
            get {
                return 2;
            }
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];

                double alpha = CurrentShift.OffsetX;

                double radian = alpha * Math.PI / 180;

                double xRotated = 1 + Math.Cos(radian) * (x - 1) - Math.Sin(radian) * y;
                double yRotated = Math.Sin(radian) * (x - 1) + Math.Cos(radian) * y;
                output[i] = xRotated * xRotated + yRotated * yRotated;
            }
        }
    }

    class Olshanskii : Generic2DTestCase, IVolumeTestCase {

        public Olshanskii(GridSizes gridSize, GridTypes gridType) : base(gridSize, gridType) {
        }

        public override int LevelSetDegree {
            get {
                return 8;
            }
        }

        public override double Solution {
            get {
                return 0.747519761188567;
            }
        }

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                double x = Math.Sqrt(input[i, 0] * input[i, 0] + input[i, 1] * input[i, 1]);
                output[i] = Math.Abs(Math.Abs(x) - 1.0) - 0.1;
            }
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < output.GetLength(0); i++) {
                double x = input[i, 0];
                double y = input[i, 1];

                double r = Math.Sqrt(x * x + y * y);
                double theta = Math.Atan2(y, x);

                output[i] = 1e5 * Math.Sin(21.0 * theta) * Math.Sin(5.0 * Math.PI * r);
            }
        }

        public override int IntegrandDegree {
            get {
                return 5;
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
            int noOfCellsPerDirection;
            switch (GridSize) {
                case GridSizes.Tiny:
                    noOfCellsPerDirection = 10;
                    break;

                case GridSizes.Small:
                    noOfCellsPerDirection = 20;
                    break;

                case GridSizes.Normal:
                    noOfCellsPerDirection = 40;
                    break;

                case GridSizes.Large:
                    noOfCellsPerDirection = 80;
                    break;

                case GridSizes.Huge:
                    noOfCellsPerDirection = 160;
                    break;

                default:
                    throw new Exception();
            }
            double[] nodes = GenericBlas.Linspace(0.0, 1.0, noOfCellsPerDirection + 1);

            GridCommons grid;
            switch (GridType) {
                case GridTypes.Structured:
                    grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                    break;

                case GridTypes.PseudoStructured:
                    grid = Grid2D.UnstructuredTriangleGrid(nodes, nodes);
                    break;

                default:
                    throw new Exception();
            }

            return grid;
        }
    }





    class ProductOfLinears : Generic2DTestCase, ISurfaceTestCase {

        public ProductOfLinears(GridSizes gridSize, GridTypes gridType)
            : base(gridSize, gridType) {
        }

        protected override IEnumerable<Shift2D> AllShifts {
            get {
                yield return new Shift2D();
            }
        }

        public override IGrid GetGrid(IDatabaseInfo db) {
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
            GridCommons grid;
            switch (GridType) {
                case GridTypes.Structured:
                    grid = Grid2D.Cartesian2DGrid(nodes, nodes);
                    break;

                case GridTypes.PseudoStructured:
                    grid = Grid2D.UnstructuredTriangleGrid(nodes, nodes);
                    break;

                case GridTypes.Unstructured:
                    throw new NotImplementedException();

                default:
                    throw new Exception();
            }

            return grid;
        }

        public override int LevelSetDegree => 2;

        public override void LevelSetInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i=0; i<output.Length; i++) {
                double x = input[i, 0] - this.CurrentShift.OffsetX;
                double y = input[i, 1] - this.CurrentShift.OffsetY;
                output[i] = x * x - y * y;
            }
        }

        public override double Solution => 4.0 * Math.Sqrt(2);
    }
}
