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
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Foundation.XDG.Quadrature.Subdivision;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Tecplot;
using CutCellQuadrature.TestCases;
using ilPSP;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;

namespace CutCellQuadrature {

    enum Modes {

        Standard,

        BruteForce,

        Adaptive,

        Regularized,

        HMFClassic,

        HMFOneStepGauss,

        HMFOneStepGaussAndStokes,

        EquivalentPolynomials,

        SayeGaussRules
    }

    public partial class Program : Application {

        private static ITestCase[] testCases = new ITestCase[] {
            //new SingleSquareStraightLineLengthTestCase(GridSizes.Tiny, GridTypes.Structured),
            new SingleSquareStraightLineVolumeTestCase(GridSizes.Tiny, GridTypes.Structured),
            //new SingleSquareParabolaLengthTestCase(GridSizes.Tiny, GridTypes.Structured),
            //new SingleSquareParabolaVolumeTestCase(GridSizes.Tiny, GridTypes.Structured),

            //new SingleCubeParaboloidVolumeTestCase(GridSizes.Tiny, GridTypes.Structured),
            //new SphereVolume3DTestCase_NoShifts(GridSizes.Tiny, GridTypes.Structured),
            //new SphereVolume3DTestCase_NoShifts(GridSizes.Small, GridTypes.Structured),
            //new SphereVolume3DTestCase_NoShifts(GridSizes.Normal, GridTypes.Structured),
            //new SphereVolume3DTestCase_NoShifts(GridSizes.Large, GridTypes.Structured),
            //new SphereVolume3DTestCase_NoShifts(GridSizes.Huge, GridTypes.Structured)

            //new ConstantIntgreandSphereSurfaceIntegral3DTestCase(GridSizes.Tiny, GridTypes.Structured),
            //new ConstantIntgreandSphereSurfaceIntegral3DTestCase(GridSizes.Small, GridTypes.Structured),
            //new ConstantIntgreandSphereSurfaceIntegral3DTestCase(GridSizes.Normal, GridTypes.Structured),
            
            /*
            new SingleSquareParabolaLengthTestCase(GridSizes.Tiny, GridTypes.Structured),
            new SingleSquareParabolaLengthTestCase(GridSizes.Small, GridTypes.Structured),
            new SingleSquareParabolaLengthTestCase(GridSizes.Normal, GridTypes.Structured),
            new SingleSquareParabolaLengthTestCase(GridSizes.Large, GridTypes.Structured),
            new SingleSquareParabolaLengthTestCase(GridSizes.Huge, GridTypes.Structured),
            */

            //new Smereka2EllipseArcLength(GridSizes.Tiny, GridTypes.Structured),
            //new Smereka2EllipseArcLength(GridSizes.Small, GridTypes.Structured),
            //new Smereka2EllipseArcLength(GridSizes.Normal, GridTypes.Structured),
            //new Smereka2EllipseArcLength(GridSizes.Large, GridTypes.Structured),
            //new Smereka2EllipseArcLength(GridSizes.Huge, GridTypes.Structured),

            //new Smereka3CircleQuadraticIntegrand(GridSizes.Tiny, GridTypes.Structured),
            //new Smereka3CircleQuadraticIntegrand(GridSizes.Small, GridTypes.Structured),
            //new Smereka3CircleQuadraticIntegrand(GridSizes.Normal, GridTypes.Structured),
            //new Smereka3CircleQuadraticIntegrand(GridSizes.Large, GridTypes.Structured),
            //new Smereka3CircleQuadraticIntegrand(GridSizes.Huge, GridTypes.Structured),


            //new MinGibou1EllipseArea(GridSizes.Tiny, GridTypes.Structured),
            //new MinGibou1EllipseArea(GridSizes.Small, GridTypes.Structured),
            //new MinGibou1EllipseArea(GridSizes.Normal, GridTypes.Structured),
            //new MinGibou1EllipseArea(GridSizes.Large, GridTypes.Structured),
            //new MinGibou1EllipseArea(GridSizes.Huge, GridTypes.Structured),


            
            //new Smereka4EllipsoidSurface(GridSizes.Tiny, GridTypes.Structured),
            //new Smereka4EllipsoidSurface(GridSizes.Small, GridTypes.Structured),
            //new Smereka4EllipsoidSurface(GridSizes.Normal, GridTypes.Structured),
            //new Smereka4EllipsoidSurface(GridSizes.Large, GridTypes.Structured),
            //new MinGibou2EllipsoidVolume(GridSizes.Large, GridTypes.Structured),
            //new MinGibou2EllipsoidVolume(GridSizes.Tiny, GridTypes.Structured),
            //new MinGibou2EllipsoidVolume(GridSizes.Normal, GridTypes.Structured),
            //new MinGibou2EllipsoidVolume(GridSizes.Large, GridTypes.Structured),

            //new KarakusCircleVolume(GridSizes.Tiny, GridTypes.Structured),
            //new KarakusCircleVolume(GridSizes.Small, GridTypes.Structured),
            //new KarakusCircleVolume(GridSizes.Normal, GridTypes.Structured),
            //new KarakusCircleVolume(GridSizes.Large, GridTypes.Structured),

            //new IBMCircleSurface()
            //new IBMCircleVolume()
            //new NACA0012ArcLength(GridSizes.Normal),
            //new NACA0012ArcLength(GridSizes.Large),
            //new NACA0012ArcLength(GridSizes.Huge),

            
            //new CircleVolume2DTestCase(GridSizes.Tiny, GridTypes.Structured),
            //new CircleVolume2DTestCase(GridSizes.Small, GridTypes.Structured),
            //new CircleVolume2DTestCase(GridSizes.Normal, GridTypes.Structured),
            //new CircleVolume2DTestCase(GridSizes.Large, GridTypes.Structured),
            //new CircleVolume2DTestCase(GridSizes.Huge, GridTypes.Structured),
            

            //new CircleArcLength(GridSizes.Small, GridTypes.Structured)
            //new Olshanskii(GridSizes.Tiny, GridTypes.PseudoStructured),
            //new Olshanskii(GridSizes.Small, GridTypes.PseudoStructured),
            //new Olshanskii(GridSizes.Normal, GridTypes.PseudoStructured),
            //new Olshanskii(GridSizes.Large, GridTypes.PseudoStructured),
            //new Olshanskii(GridSizes.Huge, GridTypes.PseudoStructured)
        };

        private ITestCase testCase;

        static void Main(string[] args) {
            InitMPI(args);

            foreach (var testCase in testCases) {

                Program app = new Program(testCase);
                app.Init(null);
                //app.Init(null, opt, "BoSSS.Platform, BoSSS.Foundation, BoSSS.Foundation.Grid, BoSSS.Foundation.XDG, BoSSS.Solution, CutCellQuadrature");
                app.RunSolverMode();
                app.ProfilingLog();
            }

            FinalizeMPI();
        }

        public Program(ITestCase testCase) {
            this.testCase = testCase;
        }

        public Program() {
        }

        private ILevelSet levelSet;

        private LevelSetTracker levelSetTracker;

        private XDGField XDGField;

        private SinglePhaseField SinglePhaseField;

        private static NumberFormatInfo formatInfo = NumberFormatInfo.InvariantInfo;

        protected override GridCommons CreateOrLoadGrid() {
            return testCase.GetGrid(CurrentSessionInfo.Database);
        }

        protected override IDatabaseInfo GetDatabase() {
            //return new StandardDBInfo("c:\\bosss_db\\GridOfTomorrow");
            return NullDatabaseInfo.Instance;
        }

        protected override void CreateFields() {
            levelSet = testCase.GetLevelSet(GridData);
            levelSetTracker = new LevelSetTracker(GridData, 
                XQuadFactoryHelper.MomentFittingVariants.Classic, // should have no effect, this app creates its own quad-rules independent of the tracker
                1, new string[] { "A", "B" }, levelSet);

            XDGField = new XDGField(
                new XDGBasis(levelSetTracker, testCase.IntegrandDegree), "XDGField");
            SinglePhaseField = new SinglePhaseField(
                new Basis(GridData, testCase.IntegrandDegree), "SinglePhaseField");

            if (levelSet is LevelSet) {
                m_IOFields.Add((LevelSet)levelSet);
            } else {
                LevelSet projectedLevelSet = new LevelSet(new Basis(GridData, 4), "projectedAnalyticLevelSet");
                projectedLevelSet.ProjectField(testCase.GetLevelSet(GridData).Evaluate);
                m_IOFields.Add(projectedLevelSet);
            }
            m_IOFields.Add(XDGField);
            m_IOFields.Add(SinglePhaseField);
        }

        protected override void SetInitial() {
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            Stopwatch globalWatch = new Stopwatch();
            globalWatch.Start();

            /*
             * Configuration options
             */
            // Only applies to subdivision-based quadrature (brute force and adaptive)
            int[] divisions = new int[] { (int)testCase.GridSize };

            // Only applies to adaptive quadrature
            int leafDivisions = -1;

            // Only applies to regularized quadrature
            double[] widths = new double[] { double.NaN };
            int[] vanishingMonents = new int[] { int.MinValue };
            int[] continuousDerivatives = new int[] { int.MinValue };

            // Only applies to HMF quadrature
            LineSegment.IRootFindingAlgorithm rootFindingAlgorithm;


            /*
             * ENTER CONFIGURATION HERE
             */

            // Export options
            int plotSuperSampling = 3;
            bool logVolumeNodes = true;
            int logVolumeNodes_selectedCell = -1;
            bool logSurfaceNodes = false;
            bool logConsole = true;
            int selectedShift = -1;

            // Quadrature variant

            Modes mode = Modes.SayeGaussRules;
            int[] orders = new int[] { 3};

            //Modes mode = Modes.HMFClassic;
            //int[] orders = new int[] { 3, 4, 5, 6, 7, 8 };
            //LevelSetSurfaceQuadRuleFactory.RestrictNodes = false;
            //LevelSetSurfaceQuadRuleFactory.UseGaussNodes = true;
            //LevelSetVolumeQuadRuleFactory.RestrictNodes = false;
            //LevelSetVolumeQuadRuleFactory.UseGaussNodes = true;
            //LevelSetVolumeQuadRuleFactory.NodeCountSafetyFactor = 1.0;

            //Modes mode = Modes.HMFClassic;
            //int[] orders = new int[] { 1 };
            //LevelSetSurfaceQuadRuleFactory.RestrictNodes = false;
            //LevelSetSurfaceQuadRuleFactory.UseGaussNodes = true;
            //LevelSetVolumeQuadRuleFactory.RestrictNodes = false;
            //LevelSetVolumeQuadRuleFactory.UseGaussNodes = true;
            //LevelSetVolumeQuadRuleFactory.NodeCountSafetyFactor = 1.0;

            //Modes mode = Modes.HMFOneStepGaussAndStokes;
            //int[] orders = new int[] { 2, 4, 6, 8, 10 };

            //Modes mode = Modes.Adaptive;
            //int[] orders = new int[] { 1 };
            //int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
            //int leafDivisions = 1;

            //Modes mode = Modes.Regularized;
            //int[] orders = Enumerable.Range(1, 20).Select(k => 2 * k - 1).ToArray();
            //widths = new double[] { 0.1, 0.2 };
            //vanishingMonents = new int[] { 0, 2, 4, 6 };
            //continuousDerivatives = new int[] { 1, 3, 5 };

            //Modes mode = Modes.BruteForce;
            //int[] orders = new int[] { 1 };
            //divisions = new int[] { 0, 1, 2, 3, 4 };

            //Modes mode = Modes.Standard;
            //int[] orders = Enumerable.Range(1, 20).Select(k => 2 * k - 1).ToArray();

            if (levelSet is IAnalyticLevelSet) {
                rootFindingAlgorithm = new AnalyticLevelSetRootFindingAlgorithm();
            } else {
                rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);
                //rootFindingAlgorithm = new LineSegment.GSLRootFindingAlgorithm(1e-14);
            }
            /*
             * END OF CONFIGURATION
             */


            SubGrid cutCellGrid = levelSetTracker.Regions.GetCutCellSubGrid();
            CellMask cellMask = CellMask.GetFullMask(GridData);
            SubGrid selectedSubGrid = new SubGrid(cellMask);

            testCase.ScaleShifts(0.5 * testCase.GridSpacing);
            double hBase = GridData.Cells.h_maxGlobal;

            string logName = ""
                + testCase.GetType().Name
                + "_" + Grid.RefElements[0].GetType().Name
                + "_" + testCase.GridSize
                + "_" + mode.ToString();

            if (logConsole) {
                string filename = logName + "_stdout.txt";
                ilPSP.Environment.StdOut.WriterS.Add(new StreamWriter(filename));
            }

            Console.WriteLine("Test case: " + testCase.GetType().Name);

            var errorMap = new Dictionary<Tuple<int, int, double, int, int>, List<Tuple<double, double, int>>>();
            foreach (int division in divisions) {
                foreach (int order in orders) {
                    foreach (double width in widths) {
                        foreach (int vanishingMoment in vanishingMonents) {
                            foreach (int continuousDerivative in continuousDerivatives) {
                                errorMap[Tuple.Create(division, order, width, vanishingMoment, continuousDerivative)] =
                                    new List<Tuple<double, double, int>>();
                            }
                        }
                    }
                }
            }

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                Console.WriteLine("Processing shift " + i + " of " + testCase.NumberOfShifts);

                if (selectedShift > 0 && i != selectedShift) {
                    i++;
                    continue;
                }

                cutCellGrid = new SubGrid(
                    levelSetTracker.Regions.GetCutCellSubGrid().VolumeMask.Intersect(
                        selectedSubGrid.VolumeMask));

                double referenceValue = SetUpConfiguration();

                StreamWriter volumeNodesLog = null;
                StreamWriter surfaceNodesLog = null;
                if (plotSuperSampling >= 0 && i == 1) {
                    PlotCurrentState(0.0, 0, plotSuperSampling, cutCellGrid);
                }

                Console.WriteLine();
                foreach (int division in divisions) {
                    Console.WriteLine("Number of divisions: " + division);

                    if (logVolumeNodes) {
                        string filename = Path.Combine(
                            Path.GetFullPath("."),
                            "volumeNodes_" + testCase.GetType().Name + "_" + i + "_" + division);
                        volumeNodesLog = new StreamWriter(filename + ".txt");
                        if (GridData.SpatialDimension == 2) {
                            volumeNodesLog.WriteLine("Cell\tNode\tx\ty\tweight");
                        } else {
                            volumeNodesLog.WriteLine("Cell\tNode\tx\ty\tz\tweight");
                        }
                    }

                    if (logSurfaceNodes) {
                        string filename = Path.Combine(
                            Path.GetFullPath("."),
                            "surfaceNodes_" + testCase.GetType().Name + "_" + i + "_" + division);
                        surfaceNodesLog = new StreamWriter(filename + ".txt");
                        if (Grid.SpatialDimension == 2) {
                            surfaceNodesLog.WriteLine("Cell\tNode\tx\ty\tweight");
                        } else {
                            surfaceNodesLog.WriteLine("Cell\tNode\tx\ty\tz\tweight");
                        }
                    }

                    foreach (int order in orders) {
                        Console.WriteLine("Order: " + order);

                        foreach (double width in widths) {
                            foreach (int vanishingMoment in vanishingMonents) {
                                foreach (int continuousDerivative in continuousDerivatives) {
                                    var result = PerformConfiguration(
                                        mode,
                                        order,
                                        division,
                                        volumeNodesLog,
                                        surfaceNodesLog,
                                        leafDivisions,
                                        vanishingMoment,
                                        continuousDerivative,
                                        width,
                                        hBase,
                                        rootFindingAlgorithm,
                                        logVolumeNodes_selectedCell);
                                    double error = Math.Abs(result.Item1 - referenceValue);

                                    var key = Tuple.Create(
                                        division, order, width, vanishingMoment, continuousDerivative);
                                    errorMap[key].Add(Tuple.Create(
                                        result.Item1 + testCase.Solution - referenceValue,
                                        error,
                                        result.Item2));

                                }
                            }
                        }
                    }

                    if (volumeNodesLog != null) {
                        volumeNodesLog.Close();
                    }
                    if (surfaceNodesLog != null) {
                        surfaceNodesLog.Close();
                    }
                }
                Console.WriteLine();

                i++;
            }
            testCase.ResetShifts();

            using (StreamWriter log = new StreamWriter(logName + ".txt")) {
                log.WriteLine("Divisions\tOrder\tWidth\tMoments\tDerivatives\tResult\tMeanError\tMeanRelError\tStdDev\tMinError\tMaxError\tMaxErrorCase\tTime");

                foreach (var entry in errorMap) {
                    if (entry.Value.Count == 0) {
                        continue;
                    }

                    IEnumerable<double> errors = entry.Value.Select((tuple) => tuple.Item2);
                    double meanResult = entry.Value.Select(tuple => tuple.Item1).Average();
                    double meanError = errors.Average();
                    double stdDev = Math.Sqrt(errors.Select((error) => error * error).Average() - meanError * meanError);
                    double time = entry.Value.Select((tuple) => tuple.Item3).Average();

                    string line = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}",
                        entry.Key.Item1, // division
                        entry.Key.Item2, // order
                        entry.Key.Item3, // width
                        entry.Key.Item4, // vanishing moments
                        entry.Key.Item5, // continuous derivative
                        meanResult.ToString(formatInfo),
                        meanError.ToString(formatInfo),
                        (meanError / testCase.Solution).ToString(formatInfo),
                        stdDev.ToString(formatInfo),
                        errors.Min().ToString(formatInfo),
                        errors.Max().ToString(formatInfo),
                        errors.IndexOfMax((d) => d) + 1,
                        time.ToString(formatInfo));


                    log.WriteLine(line);
                }
            }

            globalWatch.Stop();
            Console.WriteLine("Finished case " + testCase.GetType().Name + " after " + globalWatch.ElapsedMilliseconds + "ms");

            return dt;
        }

        private double SetUpConfiguration() {
            testCase.UpdateLevelSet(levelSet);
            levelSetTracker.UpdateTracker(__NearRegionWith: 0, incremental: false, __LevSetAllowedMovement: 2);

            XDGField.Clear();
            XDGField.GetSpeciesShadowField("A").ProjectField(
                1.0, testCase.JumpingFieldSpeciesAInitialValue, default(CellQuadratureScheme));
            XDGField.GetSpeciesShadowField("B").ProjectField(
                1.0, testCase.JumpingFieldSpeciesBInitialValue, default(CellQuadratureScheme));

            SinglePhaseField.Clear();
            SinglePhaseField.ProjectField(testCase.ContinuousFieldInitialValue);

            double referenceValue = testCase.Solution;
            if (testCase is IVolumeTestCase) {
                QuadRule standardRule = Grid.RefElements[0].GetQuadratureRule(2 * XDGField.Basis.Degree + 1);
                ScalarFieldQuadrature uncutQuadrature = new ScalarFieldQuadrature(
                    GridData,
                    XDGField,
                    new CellQuadratureScheme(
                        new FixedRuleFactory<QuadRule>(standardRule),
                        levelSetTracker.Regions.GetCutCellSubGrid().Complement().VolumeMask),
                    standardRule.OrderOfPrecision);
                uncutQuadrature.Execute();

                referenceValue -= uncutQuadrature.Result;
            }

            return referenceValue;
        }

        private Tuple<double, int, int> PerformConfiguration(
            Modes mode, 
            int order, 
            int division = 0, 
            StreamWriter volumeNodesLog = null, 
            StreamWriter surfaceNodesLog = null, 
            int leafDivisions = -1, 
            int vanishingMoment = -1, 
            int continuousDerivative = -1, 
            double width = double.NaN, 
            double hBase = double.NaN, 
            LineSegment.IRootFindingAlgorithm rootFindingAlgorithm = null, 
            int logVolumeNodes_selectedCell = -1)
        {
            SubGrid cutCellGrid = levelSetTracker.Regions.GetCutCellSubGrid();

            IQuadRuleFactory<QuadRule> volumeFactory = null;
            IQuadRuleFactory<QuadRule> edgeFactory = null;
            switch (mode) {
                case Modes.Standard: //
                    {
                        volumeFactory = new StandardQuadRuleFactory(Grid.RefElements[0]);
                        edgeFactory = new StandardQuadRuleFactory(
                            Grid.RefElements[0].FaceRefElement);
                        break;
                    }

                case Modes.BruteForce: //
                    {
                        volumeFactory = (IQuadRuleFactory<QuadRule>)new CutCellQuadRuleFactory(
                            new BruteForceSubdivisionStrategy(
                                Grid.RefElements[0], division),
                            order);
                        edgeFactory = new CutCellQuadRuleFactory(
                            new BruteForceSubdivisionStrategy(
                                Grid.RefElements[0].FaceRefElement, division),
                            order);
                        break;
                    }

                case Modes.Adaptive: //
                    {
                        volumeFactory = (IQuadRuleFactory<QuadRule>)new CutCellQuadRuleFactory(
                            new AdaptiveSubdivisionStrategy(
                                Grid.RefElements[0], levelSetTracker.DataHistories[0].Current, division),
                            leafDivisions);
                        edgeFactory = new CutCellQuadRuleFactory(
                            new AdaptiveSubdivisionStrategy(
                                Grid.RefElements[0].FaceRefElement, levelSetTracker.DataHistories[0].Current, division),
                            leafDivisions);
                        break;
                    }

                case Modes.Regularized: //
                    {
                        volumeFactory = new RegularizedQuadRuleFactory(
                            new StandardQuadRuleFactory(Grid.RefElements[0]),
                            levelSetTracker,
                            testCase.GetPolynomial(vanishingMoment, continuousDerivative),
                            0.5 * width * hBase);
                        edgeFactory = null;
                        break;
                    }

                case Modes.HMFClassic: //
                    {
                        IQuadRuleFactory<CellBoundaryQuadRule> volumeRuleFactoryEdge;
                        if (Grid.SpatialDimension == 2) {
                            LineAndPointQuadratureFactory bndrule = new LineAndPointQuadratureFactory(
                                this.Grid.RefElements[0],
                                levelSetTracker.DataHistories[0].Current,
                                true,
                                rootFindingAlgorithm);

                            volumeRuleFactoryEdge = bndrule.GetLineFactory();

                            //volumeRuleFactoryEdge = new CutLineQuadRuleFactory(
                            //    levelSetTracker,
                            //    Grid.RefElements[0],
                            //    rootFindingAlgorithm: rootFindingAlgorithm);
                        } else {
                            volumeRuleFactoryEdge = new LevelSetEdgeVolumeQuadRuleFactory(
                                levelSetTracker.DataHistories[0].Current,
                                rootFindingAlgorithm,
                                JumpTypes.Heaviside);
                        }

                        LevelSetSurfaceQuadRuleFactory surfaceFactory =
                            new LevelSetSurfaceQuadRuleFactory(
                                levelSetTracker.DataHistories[0].Current, volumeRuleFactoryEdge);


                        if (testCase is ISurfaceTestCase) {
                            volumeFactory = surfaceFactory;
                        } else {
                            volumeFactory = new LevelSetVolumeQuadRuleFactory(
                            levelSetTracker.DataHistories[0].Current,
                            volumeRuleFactoryEdge,
                            surfaceFactory,
                            JumpTypes.Heaviside);
                        }
                        edgeFactory = null;
                        break;
                    }

                case Modes.HMFOneStepGauss: //
                    {
                        if (Grid.SpatialDimension != 2)
                            throw new NotImplementedException();

                        LineAndPointQuadratureFactory bndrule = new LineAndPointQuadratureFactory(
                            this.Grid.RefElements[0],
                            levelSetTracker.DataHistories[0].Current,
                            true,
                            rootFindingAlgorithm);

                        LevelSetComboRuleFactory2 Factory = new LevelSetComboRuleFactory2(
                            levelSetTracker.DataHistories[0].Current,
                            bndrule.GetLineFactory(),
                            null,
                            _SurfaceNodesOnZeroLevset: false,
                            _UseAlsoStokes: false,
                            _DoCheck: false);

                        if (testCase is ISurfaceTestCase) {
                            volumeFactory = Factory.GetSurfaceFactory();
                        } else {
                            volumeFactory = Factory.GetVolumeFactory();
                        }
                        edgeFactory = null;

                        break;

                    }
                case Modes.HMFOneStepGaussAndStokes: //
                    {
                        if (Grid.SpatialDimension != 2)
                            throw new NotImplementedException();

                        LineAndPointQuadratureFactory bndrule = new LineAndPointQuadratureFactory(
                            this.Grid.RefElements[0],
                            levelSetTracker.DataHistories[0].Current,
                            true,
                            rootFindingAlgorithm);

                        LevelSetComboRuleFactory2 Factory = new LevelSetComboRuleFactory2(
                            levelSetTracker.DataHistories[0].Current,
                            bndrule.GetLineFactory(), bndrule.GetPointFactory(),
                            _SurfaceNodesOnZeroLevset: false,
                            _UseAlsoStokes: true,
                            _DoCheck: false);

                        if (testCase is ISurfaceTestCase) {
                            volumeFactory = Factory.GetSurfaceFactory();
                        } else {
                            volumeFactory = Factory.GetVolumeFactory();
                        }
                        edgeFactory = null;

                        break;

                    }
                case Modes.SayeGaussRules: //
                    {
                        if (testCase is ISurfaceTestCase)
                        {
                            volumeFactory = new SayeGaussRule_LevelSet2D(
                                levelSetTracker.DataHistories[0].Current,
                                rootFindingAlgorithm);
                        }
                        else
                        {
                            volumeFactory = new SayeGaussRule_Volume2D(
                                levelSetTracker.DataHistories[0].Current,
                                rootFindingAlgorithm);
                        }

                        edgeFactory = null;
                        break;
                    }

                case Modes.EquivalentPolynomials: //
                    {
                        var lineAndPointFactory = new LineAndPointQuadratureFactory(
                                    Grid.RefElements[0],
                                    levelSetTracker.DataHistories[0].Current,
                                    true,
                                    rootFindingAlgorithm);
                        if (testCase is ISurfaceTestCase) {
                            volumeFactory = new LinearReconstructionQuadRuleFactory(
                                levelSetTracker, lineAndPointFactory);
                        } else {
                            volumeFactory = new EquivalentPolynomialQuadRuleFactory(
                                new StandardQuadRuleFactory(Grid.RefElements[0]),
                                levelSetTracker,
                                lineAndPointFactory);
                        }

                        edgeFactory = null;
                        break;
                    }
            }

            if (volumeNodesLog != null) {
                WriteVolumeNodes(volumeNodesLog, volumeFactory, order, cutCellGrid, testCase, logVolumeNodes_selectedCell);
            }

            if (surfaceNodesLog != null) {
                WriteSurfaceNodes(surfaceNodesLog, edgeFactory, order, cutCellGrid);
            }

            Stopwatch timer = new Stopwatch();
            Stopwatch localTimer = new Stopwatch();
            double result = double.NaN;

            timer.Start();
            if (testCase is IVolumeTestCase
                || mode == Modes.Regularized
                || mode == Modes.HMFClassic
                || mode == Modes.HMFOneStepGauss
                || mode == Modes.HMFOneStepGaussAndStokes
                || mode == Modes.EquivalentPolynomials
                || mode == Modes.SayeGaussRules) {

                result = PerformVolumeQuadrature(
                    mode, volumeFactory, cutCellGrid, order, localTimer);
            } else {
                result = PerformSurfaceQuadrature(
                    mode, volumeFactory, edgeFactory, cutCellGrid, order, localTimer);
            }
            timer.Stop();

            return new Tuple<double, int, int>(
                result,
                (int)timer.ElapsedMilliseconds,
                (int)localTimer.ElapsedMilliseconds);
        }

        private void WriteSurfaceNodes(StreamWriter log, IQuadRuleFactory<QuadRule> ruleFactory, int order, SubGrid subGrid) {
            var edgeRules = ruleFactory.GetQuadRuleSet(
                levelSetTracker.Regions.GetCutCellSubGrid().AllEdgesMask, order);

            foreach (var chunkRulePair in edgeRules) {
                foreach (int edge in chunkRulePair.Chunk.Elements) {
                    QuadRule rule = chunkRulePair.Rule;

                    int cell = GridData.Edges.CellIndices[edge, 0];

                    NodeSet volumeVertices = new NodeSet(
                        GridData.Cells.GetRefElement(cell),
                        rule.NoOfNodes, Grid.SpatialDimension);
                    Grid.RefElements[0].TransformFaceCoordinates(
                        GridData.Edges.FaceIndices[edge, 0], rule.Nodes, volumeVertices);
                    volumeVertices.LockForever();

                    MultidimensionalArray globalVertices = MultidimensionalArray.Create(
                        1, rule.NoOfNodes, Grid.SpatialDimension);
                    GridData.TransformLocal2Global(volumeVertices, cell, 1, globalVertices, 0);
                    for (int k = 0; k < rule.NoOfNodes; k++) {
                        if (Grid.SpatialDimension == 2) {
                            log.WriteLine(
                                "{0}\t{1}\t{2}\t{3}\t{4}",
                                cell,
                                k,
                                globalVertices[0, k, 0].ToString(formatInfo),
                                globalVertices[0, k, 1].ToString(formatInfo),
                                rule.Weights[k].ToString(formatInfo));
                        } else {
                            log.WriteLine(
                                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                                cell,
                                k,
                                globalVertices[0, k, 0].ToString(formatInfo),
                                globalVertices[0, k, 1].ToString(formatInfo),
                                globalVertices[0, k, 2].ToString(formatInfo),
                                rule.Weights[k].ToString(formatInfo));
                        }
                    }
                }
            }
        }

        private void WriteVolumeNodes(StreamWriter log, IQuadRuleFactory<QuadRule> ruleFactory, int order, SubGrid subGrid, ITestCase testCase, int selectedCell = -1) {
            foreach (var chunkRulePair in ruleFactory.GetQuadRuleSet(subGrid.VolumeMask, order)) {
                foreach (int cell in chunkRulePair.Chunk.Elements) {
                    QuadRule rule = chunkRulePair.Rule;

                    MultidimensionalArray globalVertices = MultidimensionalArray.Create(
                        1, rule.NoOfNodes, Grid.SpatialDimension);
                    MultidimensionalArray metrics = levelSetTracker.DataHistories[0].Current.GetLevelSetNormalReferenceToPhysicalMetrics(
                        rule.Nodes, cell, 1);
                    GridData.TransformLocal2Global(rule.Nodes, cell, 1, globalVertices, 0);
                    
                    if (selectedCell >= 0 && cell != selectedCell) {
                        continue;
                    }

                    for (int k = 0; k < rule.NoOfNodes; k++) {
                        double weight = rule.Weights[k];

                        if (testCase is ISurfaceTestCase) {
                            // Use to get correct HMF weights in reference coordinate
                            // system (surface weights are already divided by $metrics
                            // to save this step in actual quadrature)
                            //weight *= metrics[0, k];
                        }

                        if (Grid.SpatialDimension == 2) {
                            log.WriteLine(
                                "{0}\t{1}\t{2}\t{3}\t{4}",
                                cell,
                                k,
                                globalVertices[0, k, 0].ToString(formatInfo),
                                globalVertices[0, k, 1].ToString(formatInfo),
                                Math.Round(weight, 2).ToString(formatInfo));
                        } else {
                            log.WriteLine(
                                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                                cell,
                                k,
                                globalVertices[0, k, 0].ToString(formatInfo),
                                globalVertices[0, k, 1].ToString(formatInfo),
                                globalVertices[0, k, 2].ToString(formatInfo),
                                weight.ToString(formatInfo));
                        }
                    }
                }
            }
        }

        private double PerformVolumeQuadrature(Modes mode, IQuadRuleFactory<QuadRule> factory, SubGrid cutCellGrid, int order, Stopwatch timer) {
            using (new FuncTrace()) {
                ScalarFieldQuadrature quadrature;
                CellQuadratureScheme quadInstr = new CellQuadratureScheme(
                    factory, cutCellGrid.VolumeMask);
                if (testCase is ISurfaceTestCase) {
                    quadrature = new ScalarFieldQuadrature(GridData, SinglePhaseField, quadInstr, order);
                } else {
                    quadrature = new ScalarFieldQuadrature(GridData, XDGField, quadInstr, order);
                }

                timer.Start();
                quadrature.Execute();
                timer.Stop();

                return quadrature.Result;
            }
        }

        private double PerformSurfaceQuadrature(Modes mode, IQuadRuleFactory<QuadRule> volumeFactory, IQuadRuleFactory<QuadRule> edgeFactory, SubGrid cutCellGrid, int order, Stopwatch timer) {
            using (new FuncTrace()) {
                CellQuadratureScheme volInstr = new CellQuadratureScheme(
                    volumeFactory, cutCellGrid.VolumeMask);
                CellBoundaryQuadratureScheme edgeInstr = new CellBoundaryQuadratureScheme(
                    new CellBoundaryFromEdgeRuleFactory<CellBoundaryQuadRule>(
                        GridData, Grid.RefElements[0], edgeFactory),
                    cutCellGrid.VolumeMask);

                ScalarFieldLevelSetIntegrator quadrature = new ScalarFieldLevelSetIntegrator(
                    levelSetTracker,
                    SinglePhaseField,
                    volInstr.Compile(GridData, order),
                    edgeInstr.Compile(GridData, order),
                    order,
                    cutCellGrid,
                    0);

                timer.Start();
                double result = quadrature.ExecuteA().Storage.Sum();
                timer.Stop();

                return result;
            }
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            PlotCurrentState(physTime, timestepNo, 0, null);
        }

        private void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling, SubGrid subGrid) {
            Tecplot tecplot = new Tecplot(GridData, true, false, (uint)superSampling, subGrid);
            //Tecplot tecplot = new Tecplot(m_Context, true, false, (uint)superSampling, null);
            string path = Path.Combine(Path.GetFullPath("."), "plot_" + testCase.GetType().Name);
            tecplot.PlotFields(path, physTime, m_IOFields);
        }

        private static double Regression(double[] xValues, double[] yValues) {
            double xAvg = xValues.Average();
            double yAvg = yValues.Average();

            double v1 = 0.0;
            double v2 = 0.0;

            for (int i = 0; i < yValues.Length; i++) {
                v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
                v2 += Math.Pow(xValues[i] - xAvg, 2);
            }

            double a = v1 / v2;
            double b = yAvg - a * xAvg;

            return a;
        }

        public static void Create2DGridPlotFile(GridData gridData, string filename) {
            if (gridData.SpatialDimension != 2) {
                throw new ArgumentException("Only supported for 2D grids");
            }

            //MultidimensionalArray edgeVertices = gridData.Grid.RefElements[0].FaceRefElement.Vertices;
            //int noOfVerticesPerEdge = edgeVertices.GetLength(0);

            //NodeSet volumeVertices = MultidimensionalArray.Create(
            //    noOfVerticesPerEdge, gridData.SpatialDimension);

            int noOfVerticesPerEdge = 2;
            MultidimensionalArray globalVertices = MultidimensionalArray.Create(
                1,
                gridData.Edges.Count,
                noOfVerticesPerEdge,
                gridData.SpatialDimension);

            for (int i = 0; i < gridData.Edges.Count; i++) {
                int cell = gridData.Edges.CellIndices[i, 0];
                int e = gridData.Edges.FaceIndices[i, 0];
                //Cell cellInfo = gridData.Cells.GetCell(cell);

                NodeSet volumeVertices = gridData.Grid.GetRefElement(cell).GetFaceVertices(e);

                Debug.Assert(volumeVertices.NoOfNodes == 2);


                //gridData.Grid.RefElements[0].TransformFaceCoordinates(e, edgeVertices, volumeVertices);
                gridData.TransformLocal2Global(
                    volumeVertices,
                    cell, 1,
                    globalVertices.ExtractSubArrayShallow(-1, i, -1, -1),
                    0);
            }

            double xMax = globalVertices.ExtractSubArrayShallow(-1, -1, -1, 0).Max(x => x);
            double xMin = globalVertices.ExtractSubArrayShallow(-1, -1, -1, 0).Min(x => x);
            double yMax = globalVertices.ExtractSubArrayShallow(-1, -1, -1, 1).Max(y => y);
            double yMin = globalVertices.ExtractSubArrayShallow(-1, -1, -1, 1).Min(y => y);

            using (var file = new StreamWriter(filename + ".gp")) {
                file.WriteLine("reset");
                file.WriteLine("set terminal post enhanced 20 dl 2.25");
                file.WriteLine("set output \"" + filename + ".ps\"");
                file.WriteLine("set size square");
                file.WriteLine("set origin 0.0, 0.0");
                file.WriteLine("set xrange [{0}:{1}]", xMin.ToStringDot(), xMax.ToStringDot());
                file.WriteLine("set yrange [{0}:{1}]", yMin.ToStringDot(), yMax.ToStringDot());

                for (int i = 0; i < gridData.Edges.Count; i++) {
                    file.WriteLine(
                        "set arrow from {0},{1} to {2},{3} nohead",
                        globalVertices[0, i, 0, 0].ToStringDot(),
                        globalVertices[0, i, 0, 1].ToStringDot(),
                        globalVertices[0, i, 1, 0].ToStringDot(),
                        globalVertices[0, i, 1, 1].ToStringDot());
                }

                // Issue plot command without effect to force plot
                file.WriteLine("plot NaN notitle");
            }
        }
    }
}