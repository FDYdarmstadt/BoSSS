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

//using System;
//using System.Collections;
//using System.Collections.Generic;
//using System.Diagnostics;
//using System.Globalization;
//using System.IO;
//using System.Linq;
//using BoSSS.Foundation;
//using BoSSS.Foundation.Grid;
//using BoSSS.Foundation.IO;
//using BoSSS.Foundation.Quadrature;
//using BoSSS.Foundation.XDG;
//using BoSSS.Foundation.XDG.Quadrature;
//using BoSSS.Foundation.XDG.Quadrature.HMF;
//using BoSSS.Foundation.XDG.Quadrature.Subdivision;
//using BoSSS.Platform;
//using BoSSS.Solution;
//using BoSSS.Solution.Tecplot;
//using CutCellQuadrature.TestCases;
//using ilPSP.Tracing;

//namespace CutCellQuadrature {

//    enum Modes {

//        Standard,

//        BruteForce,

//        Adaptive,

//        RegularizedStandard,

//        RegularizedBruteForce,

//        RegularizedAdaptive,

//        HybridBruteForce,

//        MomentFitting
//    }

//    class ProgramMultipleLevelSets : Application {

//        private static ITestCase[] testCases = new ITestCase[] {
//           // new Smereka2EllipseArcLength(GridSizes.Normal, GridTypes.Structured),
//            //new Smereka3CircleQuadraticIntegrand(GridSizes.Tiny, GridTypes.Unstructured),
//            //new Smereka3CircleQuadraticIntegrand(GridSizes.Small, GridTypes.Unstructured),
//            //new Smereka3CircleQuadraticIntegrand(GridSizes.Normal, GridTypes.Unstructured),
//            //new Smereka3CircleQuadraticIntegrand(GridSizes.Large, GridTypes.Unstructured),
//            //new Smereka3CircleQuadraticIntegrand(GridSizes.Huge, GridTypes.Unstructured),
//            new Smereka4EllipsoidSurface(GridSizes.Tiny, GridTypes.Structured),
//            //new MinGibou1EllipseArea(GridSizes.Normal, GridTypes.Structured),
//            //new MinGibou2EllipsoidVolume(GridSizes.Tiny, GridTypes.Structured),
//            //new Smereka5SphereQuadraticIntegrand(GridSizes.Tiny, GridTypes.Structured),
//            //new Smereka5SphereQuadraticIntegrand(GridSizes.Tiny, GridTypes.Unstructured),
//            //new Smereka5SphereQuadraticIntegrand(GridSizes.Small, GridTypes.Unstructured),
//            //new Smereka5SphereQuadraticIntegrand(GridSizes.Normal, GridTypes.Unstructured),
//            //new Smereka5SphereQuadraticIntegrand(GridSizes.Large, GridTypes.Unstructured)
//        };

//        private ITestCase testCase;

//        static void Main(string[] args) {
//            bool dummy;
//            InitMPI(args, out dummy);

//            CommandLineOptions opt = new CommandLineOptions();
//            foreach (var testCase in testCases) {
//                ProgramMultipleLevelSets app = new ProgramMultipleLevelSets(testCase);
//                //app.Init(null, opt, "");
//                app.Init(null, opt, "BoSSS.Foundation.XDG,CutCellQuadrature");
//                app.RunSolverMode();
//            }

//            MPI.Wrappers.csMPI.Raw.mpiFinalize();
//        }

//        public ProgramMultipleLevelSets(ITestCase testCase) {
//            this.testCase = testCase;
//        }

//        // Dummy level set
//        private LevelSet levelSet1;

//        private ILevelSet levelSet2;

//        private LevelSetTracker levelSetTracker, levelSetTrackerHelp;

//        private XDGField XDGField;

//        private SinglePhaseField SinglePhaseField;

//        int levelSetIndex;

//        protected override void CreateOrLoadGrid() {
//            m_Context.SetGrid(testCase.GetGrid(m_Context));
//        }

//        protected override IIODriver GetIODriver() {
//            return new StandartFsDriver(new string[] { "\\\\fdyprime\\userspace\\kallendorf\\bosss-db" });
//        }

//        protected override void CreateFields() {
//            levelSetIndex = 0;

//            levelSet1 = new LevelSet(new Basis(m_Context, 2), "dummyLevelSet");
//            levelSet1.AccConstant(1.0);
//            levelSet2 = testCase.GetLevelSet(m_Context);
//            levelSetTracker = new LevelSetTracker(
//                m_Context, 1, new string[,] { { "A1", "B1" }, { "A2", "B2" } }, levelSet1, levelSet2);
//            levelSetTrackerHelp = new LevelSetTracker(
//                m_Context, 1, new string[] { "A", "B" }, levelSetTracker.LevelSets[levelSetIndex]);
//            XDGField = new XDGField(
//                new XDGBasis(m_Context, testCase.IntegrandDegree, levelSetTrackerHelp),
//                "XDGField");
//            SinglePhaseField = new SinglePhaseField(
//                new Basis(m_Context, testCase.IntegrandDegree),
//                "SinglePhaseField");
//            for (int k = 0; k < levelSetTracker.LevelSets.Count; k++) {
//                if (levelSetTracker.LevelSets[k] is LevelSet) {
//                    m_IOFields.Add((LevelSet)levelSetTracker.LevelSets[k]);
//                } else {
//                    LevelSet projectedLevelSet = new LevelSet(new Basis(m_Context, 4), "projectedAnalyticLevelSet");
//                    projectedLevelSet.ProjectField(testCase.GetLevelSet(m_Context).Evaluate);
//                    m_IOFields.Add(projectedLevelSet);
//                }
//            }
//            m_IOFields.Add(XDGField);
//            m_IOFields.Add(SinglePhaseField);
//        }

//        protected override void SetInitial() {
//        }

//        protected override void CreateEquationsAndSolvers() {
//        }

//        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
//            Console.WriteLine("Test case: " + testCase.GetType().Name);

//            Stopwatch globalWatch = new Stopwatch();
//            globalWatch.Start();

//            SubGrid cutCellGrid = levelSetTracker.GetCutCellSubgrid4LevSet(levelSetIndex);
//            CellMask uncutCellmask = cutCellGrid.VolumeMask.Complement();

//            SpeciesId idA = levelSetTrackerHelp.GetSpeciesId("A");
//            SpeciesId idB = levelSetTrackerHelp.GetSpeciesId("B");


//            //double[] xValues = new double[] { 0.8, 0.4, 0.2, 0.1, 0.05 };
//            //double[] yValues = new double[] { 0.00604029688860881, 0.00150669487854926, 0.000552730954015259, 0.000147125750204156, 3.75884357728008E-05 };
//            //double bla = Regression(
//            //    xValues.Select(x => Math.Log(x)).ToArray(),
//            //    yValues.Select(x => Math.Log(x)).ToArray());
//            //Console.WriteLine("EOC: " + bla);
//            //Environment.Exit(-1);



//            // Modes mode = Modes.MomentFitting;
//            //int[] orders = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
//            //int[] orders = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
//            //int[] orders = new int[] { 0, 1, 2, 3, 4, 5, 6 };
//            //int[] orders = new int[] { 0, 1, 2, 3, 4 };
//            //int[] orders = new int[] { 0, 1, 2, 3 };
//            //int[] orders = new int[] { 0, 1, 2 };
//            //int[] orders = new int[] { 0, 1 };
//            //int[] orders = new int[] { 1, 0 };
//            //int[] orders = new int[] { 1, 2 };
//            //int[] orders = new int[] { 0 };
//            //int[] orders = new int[] { 2 };
//            //int[] orders = new int[] { 5 };
//            //int[] orders = new int[] { 4, 5 };
//            //int[] divisions = new int[] { 0 };
//            //int[] divisions = new int[] { 1, 3, 5, 7, 9, 11, 13, 15, 17 };
//            //double[] widths = new double[] { -1 };
//            //int[] vanishingMonents = new int[] { int.MinValue };
//            //int[] continuousDerivatives = new int[] { int.MinValue };
//            //int leafDivisions = -1;

//            Modes mode = Modes.Adaptive;
//            //int[] orders = new int[] { 5 };
//            ////int[] orders = new int[] { 3 };
//            int[] orders = new int[] { 1 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
//            //int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
//            ////int[] divisions = new int[] {  7, 8 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6 };
//            int[] divisions = new int[] { 1, 6 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4 };
//            ////int[] divisions = new int[] { 4, 6, 7 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3 };
//            ////int[] divisions = new int[] { 0, 6 };
//            double[] widths = new double[] { double.NaN };
//            int[] vanishingMonents = new int[] { int.MinValue };
//            int[] continuousDerivatives = new int[] { int.MinValue };
//            int leafDivisions = 1;
//            ////int leafDivisions = 1;
//            ////int leafDivisions = 2;

//            //Modes mode = Modes.BruteForce;
//            //int[] orders = new int[] { 1 };
//            ////int[] orders = new int[] { 3 };
//            ////int[] orders = new int[] { 5 };
//            //int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4, 5, 6 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4, 5 };
//            ////int[] divisions = new int[] { 0, 1, 2, 3, 4 };
//            ////int[] divisions = new int[] { 0, 1, 2 };
//            ////int[] divisions = new int[] { 0 };
//            //double[] widths = new double[] { double.NaN };
//            //int[] vanishingMonents = new int[] { int.MinValue };
//            //int[] continuousDerivatives = new int[] { int.MinValue };
//            //int leafDivisions = -1;
//            //int edgeLeafDivisions = leafDivisions;

//            //Modes mode = Modes.Standard;
//            //int[] orders = Enumerable.Range(1, 20).Select(k => 2 * k - 1).ToArray();
//            ////int[] orders = Enumerable.Range(1, 41).Select(k => 2 * k - 1).ToArray();
//            ////int[] orders = new int[] { 81 };
//            //int[] divisions = new int[] { 0 };
//            //double[] widths = new double[] { double.NaN };
//            //int[] vanishingMonents = new int[] { int.MinValue };
//            //int[] continuousDerivatives = new int[] { int.MinValue };
//            //int leafDivisions = -1;
//            //int edgeLeafDivisions = leafDivisions;

//            int plotSuperSampling = 0;
//            bool logVolumeNodes = false;
//            bool logSurfaceNodes = false;
//            //int selectedShift = 62;
//            //int selectedShift = -1;
//            int selectedShift = 1;

//            CellMask cellMask = CellMask.GetFullMask(m_Context.GridDat);
//            //CellMask cellMask = new CellMask(
//            //    m_Context.GridDat,
//            //    Chunk.GetSingleElementChunk(268),
//            //    Chunk.GetSingleElementChunk(286));
//            //CellMask cellMask = new CellMask(
//            //    m_Context.GridDat,
//            //    Chunk.GetSingleElementChunk(268));
//            SubGrid selectedSubGrid = new SubGrid(m_Context.GridDat, cellMask);

//            LineSegment.IRootFindingAlgorithm rootFindingAlgorithm;

//            if (levelSetTracker.LevelSets[levelSetIndex] is IAnalyticLevelSet) {

//                rootFindingAlgorithm = new AnalyticLevelSetRootFindingAlgorithm();
//            } else {
//                rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);
//                //rootFindingAlgorithm = new LineSegment.GSLRootFindingAlgorithm(1e-14);
//            }

//            testCase.ScaleShifts(0.5 * testCase.GridSpacing);
//            double hBase = m_Context.GridDat.h_max.Max();

//            var errorMap = new Dictionary<Tuple<int, int, double, int, int>, List<Tuple<double, double, int>>>();
//            foreach (int division in divisions) {
//                foreach (int order in orders) {
//                    foreach (double width in widths) {
//                        foreach (int vanishingMoment in vanishingMonents) {
//                            foreach (int continuousDerivative in continuousDerivatives) {
//                                errorMap[Tuple.Create(division, order, width, vanishingMoment, continuousDerivative)] =
//                                    new List<Tuple<double, double, int>>();
//                            }
//                        }
//                    }
//                }
//            }

//            int i = 1;
//            while (testCase.ProceedToNextShift()) {
//                Console.WriteLine("Processing shift " + i + " of " + testCase.NumberOfShifts);

//                if (selectedShift > 0 && i != selectedShift) {
//                    i++;
//                    continue;
//                }

//                // testCase.UpdateLevelSet(levelSetTracker.LevelSets[levelSetIndex]);
//                testCase.UpdateLevelSet(levelSet2);
//                //levelSetTracker.UpdateTracker(0, 2);
//                levelSetTracker.UpdateTracker(0, new int[] { 1, 1 });

//                XDGField.Clear();
//                XDGField.GetSpeciesShadowField("A").ProjectField(1.0, testCase.JumpingFieldSpeciesAInitialValue, null);
//                XDGField.GetSpeciesShadowField("B").ProjectField(1.0, testCase.JumpingFieldSpeciesBInitialValue, null);

//                SinglePhaseField.Clear();
//                SinglePhaseField.ProjectField(testCase.ContinuousFieldInitialValue);

//                cutCellGrid = levelSetTracker.GetCutCellSubgrid4LevSet(levelSetIndex);
//                uncutCellmask = cutCellGrid.VolumeMask.Complement();

//                cutCellGrid = new SubGrid(m_Context.GridDat, cutCellGrid.VolumeMask.Intersect(selectedSubGrid.VolumeMask));

//                StreamWriter volumeNodesLog = null;
//                StreamWriter surfaceNodesLog = null;
//                if (plotSuperSampling >= 0) {
//                    PlotCurrentState(0.0, 0, plotSuperSampling, cutCellGrid);
//                }

//                double referenceValue = testCase.Solution;
//                if (testCase is IVolumeTestCase) {
//                    QuadRule standardRule = m_Context.Grid.GridSimplex.GetQuadratureRule(2 * XDGField.Basis.Degree + 1);
//                    ScalarFieldQuadrature uncutQuadrature = new ScalarFieldQuadrature(
//                        m_Context,
//                        XDGField,

//                        new CellQuadratureScheme(
//                            new FixedRuleFactory<QuadRule>(
//                                m_Context.Grid.GridSimplex, standardRule),
//                            uncutCellmask),
//                        standardRule.OrderOfPrecision);
//                    uncutQuadrature.Execute();
//                    referenceValue -= uncutQuadrature.Result;
//                }

//                Console.WriteLine();
//                foreach (int division in divisions) {
//                    Console.WriteLine("Number of divisions: " + division);

//                    //  levelSetTracker.SetBruteForceQuadratureRules(5, 2);
//                    //levelSetTracker.SetBruteForceQuadratureRules(7, 1);

//                    if (logVolumeNodes) {
//                        string filename = Path.Combine(
//                            Path.GetFullPath("."),
//                            "volumeNodes_" + testCase.GetType().Name + "_" + i + "_" + division);
//                        volumeNodesLog = new StreamWriter(filename + ".txt");
//                        if (m_Context.Grid.SpatialDimension == 2) {
//                            volumeNodesLog.WriteLine("Cell\tNode\tx\ty\tweight");
//                        } else {
//                            volumeNodesLog.WriteLine("Cell\tNode\tx\ty\tz\tweight");
//                        }
//                    }

//                    if (logSurfaceNodes) {
//                        string filename = Path.Combine(
//                            Path.GetFullPath("."),
//                            "surfaceNodes_" + testCase.GetType().Name + "_" + i + "_" + division);
//                        surfaceNodesLog = new StreamWriter(filename + ".txt");
//                        if (m_Context.Grid.SpatialDimension == 2) {
//                            surfaceNodesLog.WriteLine("Cell\tNode\tx\ty\tweight");
//                        } else {
//                            surfaceNodesLog.WriteLine("Cell\tNode\tx\ty\tz\tweight");
//                        }
//                    }

//                    foreach (int order in orders) {
//                        Console.WriteLine("Order: " + order);

//                        foreach (double width in widths) {
//                            foreach (int vanishingMoment in vanishingMonents) {
//                                foreach (int continuousDerivative in continuousDerivatives) {
//                                    IQuadRuleFactory<QuadRule> volumeFactory = null;
//                                    IQuadRuleFactory<QuadRule> edgeFactory = null;
//                                    switch (mode) {
//                                        case Modes.Standard:
//                                            volumeFactory = new StandardQuadRuleFactory(
//                                                m_Context.Grid.GridSimplex);
//                                            edgeFactory = new StandardQuadRuleFactory(
//                                                m_Context.Grid.GridSimplex.EdgeSimplex);
//                                            break;

//                                        case Modes.BruteForce:
//                                            volumeFactory = (IQuadRuleFactory<QuadRule>)new CutCellQuadRuleFactory(
//                                                new BruteForceSubdivisionStrategy(
//                                                    m_Context.Grid.GridSimplex, division),
//                                                order);
//                                            edgeFactory = new CutCellQuadRuleFactory(
//                                                new BruteForceSubdivisionStrategy(
//                                                    m_Context.Grid.GridSimplex.EdgeSimplex, division),
//                                                order);
//                                            break;

//                                        case Modes.Adaptive:
//                                            volumeFactory = (IQuadRuleFactory<QuadRule>)new CutCellQuadRuleFactory(
//                                                new AdaptiveSubdivisionStrategy(
//                                                    m_Context.Grid.GridSimplex, levelSetTracker, levelSetIndex,
//                                                    division),
//                                                leafDivisions);
//                                            edgeFactory = new CutCellQuadRuleFactory(
//                                                new AdaptiveSubdivisionStrategy(
//                                                    m_Context.Grid.GridSimplex.EdgeSimplex, levelSetTracker, levelSetIndex, division),
//                                                leafDivisions);
//                                            break;

//                                        case Modes.RegularizedStandard:
//                                            volumeFactory = new RegularizedQuadRuleFactory(
//                                                new StandardQuadRuleFactory(m_Context.Grid.GridSimplex),
//                                                levelSetTracker,
//                                                testCase.GetPolynomial(vanishingMoment, continuousDerivative),
//                                                0.5 * width * hBase);
//                                            edgeFactory = null;
//                                            break;

//                                        case Modes.RegularizedBruteForce:
//                                            volumeFactory = new RegularizedQuadRuleFactory(
//                                                new CutCellQuadRuleFactory(
//                                                    new BruteForceSubdivisionStrategy(m_Context.Grid.GridSimplex, division),
//                                                    leafDivisions),
//                                                levelSetTracker,
//                                                testCase.GetPolynomial(vanishingMoment, continuousDerivative),
//                                                0.5 * width * hBase / Math.Pow(2.0, division));
//                                            edgeFactory = null;
//                                            break;

//                                        case Modes.RegularizedAdaptive:
//                                            volumeFactory = new RegularizedQuadRuleFactory(
//                                                new CutCellQuadRuleFactory(
//                                                    new AdaptiveSubdivisionStrategy(
//                                                        m_Context.Grid.GridSimplex, levelSetTracker, division),
//                                                        leafDivisions),
//                                                    levelSetTracker,
//                                                    testCase.GetPolynomial(vanishingMoment, continuousDerivative),
//                                                    0.5 * width * hBase / Math.Pow(2.0, division));
//                                            edgeFactory = null;
//                                            break;

//                                        case Modes.MomentFitting:
//                                            IQuadRuleFactory<CellBoundaryQuadRule> volumeRuleFactoryEdge;
//                                            if (m_Context.Grid.SpatialDimension == 2) {
//                                                volumeRuleFactoryEdge = new CutLineQuadRuleFactory(
//                                                    levelSetTracker, levelSetIndex, rootFindingAlgorithm, JumpType.Heaviside);
//                                            } else {
//                                                volumeRuleFactoryEdge = new LevelSetEdgeVolumeQuadRuleFactory(
//                                                    levelSetTracker, levelSetIndex, rootFindingAlgorithm, JumpType.Heaviside);
//                                            }

//                                            LevelSetSurfaceQuadRuleFactory surfaceFactory = new LevelSetSurfaceQuadRuleFactory(levelSetTracker, levelSetIndex, volumeRuleFactoryEdge);

//                                            if (testCase is ISurfaceTestCase) {
//                                                volumeFactory = surfaceFactory;
//                                            } else {
//                                                volumeFactory = new LevelSetVolumeQuadRuleFactory(
//                                                    m_Context.Grid.GridSimplex,
//                                                    levelSetTracker, levelSetIndex,
//                                                    volumeRuleFactoryEdge,
//                                                    surfaceFactory,
//                                                    JumpType.Heaviside);
//                                            }
//                                            edgeFactory = null;




//                                            //volumeFactory = volumeRuleFactoryEdge;

//                                            //volumeFactory = new LevelSetEdgeSurfaceQuadRuleFactory(
//                                            //    levelSetTracker, new CutLineOnEdgeQuadRuleFactory(levelSetTracker, rootFindingAlgorithm, JumpType.Heaviside), JumpType.Heaviside);

//                                            break;
//                                    }

//                                    var key = Tuple.Create(division, order, width, vanishingMoment, continuousDerivative);
//                                    Stopwatch timer = new Stopwatch();
//                                    Stopwatch localTimer = new Stopwatch();
//                                    double result = double.NaN;
//                                    double error = 0.0;

//                                    timer.Start();
//                                    if (testCase is IVolumeTestCase
//                                        || mode == Modes.RegularizedStandard
//                                        || mode == Modes.RegularizedAdaptive
//                                        || mode == Modes.MomentFitting) {

//                                        result = PerformVolumeQuadrature(
//                                            mode, volumeFactory, cutCellGrid, order, localTimer);
//                                    } else {
//                                        result = PerformSurfaceQuadrature(
//                                            mode, volumeFactory, edgeFactory, cutCellGrid, order, localTimer);
//                                    }
//                                    timer.Stop();

//                                    error = Math.Abs(result - referenceValue);

//                                    if (volumeNodesLog != null) {
//                                        foreach (var chunkRulePair in volumeFactory.GetQuadRuleSet(cutCellGrid.VolumeMask, order)) {
//                                            foreach (int cell in chunkRulePair.Chunk.Elements) {
//                                                QuadRule rule = chunkRulePair.Rule;

//                                                MultidimensionalArray globalVertices = MultidimensionalArray.Create(
//                                                    1, rule.NoOfNodes, m_Context.Grid.SpatialDimension);
//                                                m_Context.GridDat.TransformLocal2Global(
//                                                    rule.Nodes, globalVertices, cell, 1, 0);
//                                                for (int k = 0; k < rule.NoOfNodes; k++) {
//                                                    if (m_Context.Grid.SpatialDimension == 2) {
//                                                        volumeNodesLog.WriteLine(
//                                                            "{0}\t{1}\t{2}\t{3}\t{4}",
//                                                            cell,
//                                                            k,
//                                                            globalVertices[0, k, 0].ToString(NumberFormatInfo.InvariantInfo),
//                                                            globalVertices[0, k, 1].ToString(NumberFormatInfo.InvariantInfo),
//                                                            rule.Weights[k].ToString(NumberFormatInfo.InvariantInfo));
//                                                    } else {
//                                                        volumeNodesLog.WriteLine(
//                                                            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
//                                                            cell,
//                                                            k,
//                                                            globalVertices[0, k, 0].ToString(NumberFormatInfo.InvariantInfo),
//                                                            globalVertices[0, k, 1].ToString(NumberFormatInfo.InvariantInfo),
//                                                            globalVertices[0, k, 2].ToString(NumberFormatInfo.InvariantInfo),
//                                                            rule.Weights[k].ToString(NumberFormatInfo.InvariantInfo));
//                                                    }
//                                                }
//                                            }
//                                        }
//                                    }

//                                    if (surfaceNodesLog != null) {
//                                        foreach (var chunkRulePair in edgeFactory.GetQuadRuleSet(levelSetTracker.GetCutCellSubgrid4LevSet(levelSetIndex).AllEdgesMask, order)) {
//                                            foreach (int edge in chunkRulePair.Chunk.Elements) {
//                                                QuadRule rule = chunkRulePair.Rule;

//                                                int cell = m_Context.GridDat.Edges[edge, 0];
//                                                int localEdge = -1;
//                                                for (int e = 0; e < m_Context.Grid.GridSimplex.NoOfEdges; e++) {
//                                                    if (Math.Abs(m_Context.GridDat.LocalCellIndexToEdges[cell, e]) - 1 == edge) {
//                                                        localEdge = e;
//                                                        break;
//                                                    }
//                                                }
//                                                if (localEdge < 0) {
//                                                    throw new Exception("Internal error");
//                                                }

//                                                MultidimensionalArray volumeVertices = MultidimensionalArray.Create(
//                                                    rule.NoOfNodes, m_Context.Grid.SpatialDimension);
//                                                m_Context.Grid.GridSimplex.EdgeToVolumeCoordinates(
//                                                    localEdge, rule.Nodes, volumeVertices);

//                                                MultidimensionalArray globalVertices = MultidimensionalArray.Create(
//                                                    1, rule.NoOfNodes, m_Context.Grid.SpatialDimension);
//                                                m_Context.GridDat.TransformLocal2Global(
//                                                    volumeVertices, globalVertices, cell, 1, 0);
//                                                for (int k = 0; k < rule.NoOfNodes; k++) {
//                                                    if (m_Context.Grid.SpatialDimension == 2) {
//                                                        surfaceNodesLog.WriteLine(
//                                                            "{0}\t{1}\t{2}\t{3}\t{4}",
//                                                            edge,
//                                                            k,
//                                                            globalVertices[0, k, 0].ToString(NumberFormatInfo.InvariantInfo),
//                                                            globalVertices[0, k, 1].ToString(NumberFormatInfo.InvariantInfo),
//                                                            rule.Weights[k].ToString(NumberFormatInfo.InvariantInfo));
//                                                    } else {
//                                                        surfaceNodesLog.WriteLine(
//                                                            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
//                                                            edge,
//                                                            k,
//                                                            globalVertices[0, k, 0].ToString(NumberFormatInfo.InvariantInfo),
//                                                            globalVertices[0, k, 1].ToString(NumberFormatInfo.InvariantInfo),
//                                                            globalVertices[0, k, 2].ToString(NumberFormatInfo.InvariantInfo),
//                                                            rule.Weights[k].ToString(NumberFormatInfo.InvariantInfo));
//                                                    }
//                                                }
//                                            }
//                                        }
//                                    }

//                                    errorMap[key].Add(Tuple.Create(result + testCase.Solution - referenceValue, error, (int)timer.ElapsedMilliseconds));
//                                }
//                            }
//                        }
//                    }

//                    if (volumeNodesLog != null) {
//                        volumeNodesLog.Close();
//                    }
//                    if (surfaceNodesLog != null) {
//                        surfaceNodesLog.Close();
//                    }
//                }
//                Console.WriteLine();

//                i++;
//            }
//            testCase.ResetShifts();

//            string logName = "log_"
//                + testCase.GetType().Name
//                + "_" + m_Context.Grid.GridSimplex.GetType().Name
//                + "_" + testCase.GridSize
//                + "_" + mode.ToString();
//            using (StreamWriter log = new StreamWriter(logName + ".txt")) {
//                log.WriteLine("Divisions\tOrder\tWidth\tMoments\tDerivatives\tResult\tMeanError\tMeanRelError\tStdDev\tMinError\tMaxError\tMaxErrorCase\tTime");

//                foreach (var entry in errorMap) {
//                    if (entry.Value.Count == 0) {
//                        continue;
//                    }

//                    IEnumerable<double> errors = entry.Value.Select((tuple) => tuple.Item2);
//                    double meanResult = entry.Value.Select(tuple => tuple.Item1).Average();
//                    double meanError = errors.Average();
//                    double stdDev = Math.Sqrt(errors.Select((error) => error * error).Average() - meanError * meanError);
//                    double time = entry.Value.Select((tuple) => tuple.Item3).Average();

//                    log.WriteLine(
//                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}",
//                        entry.Key.Item1, // division
//                        entry.Key.Item2, // order
//                        entry.Key.Item3, // width
//                        entry.Key.Item4, // vanishing moments
//                        entry.Key.Item5, // continuous derivative
//                        meanResult.ToString(NumberFormatInfo.InvariantInfo),
//                        meanError.ToString(NumberFormatInfo.InvariantInfo),
//                        (meanError / testCase.Solution).ToString(NumberFormatInfo.InvariantInfo),
//                        stdDev.ToString(NumberFormatInfo.InvariantInfo),
//                        errors.Min().ToString(NumberFormatInfo.InvariantInfo),
//                        errors.Max().ToString(NumberFormatInfo.InvariantInfo),
//                        errors.IndexOfMax((d) => d) + 1,
//                        time.ToString(NumberFormatInfo.InvariantInfo));
//                }
//            }

//            using (StreamWriter log = new StreamWriter(logName + "_paper.txt")) {
//                log.WriteLine("Divisions\tMeanRelError\tTime\tEOC\tRegression");

//                foreach (var group in errorMap.GroupBy(e => e.Key.Item2)) {
//                    log.WriteLine("Order: " + group.Key);

//                    double lastMeanError = -1.0;
//                    int lastDivision = -1;
//                    int index = 1;
//                    foreach (var entry in group) {
//                        if (entry.Value.Count == 0) {
//                            continue;
//                        }

//                        IEnumerable<double> errors = entry.Value.Select((tuple) => tuple.Item2);
//                        double meanError = errors.Average();
//                        double stdDev = Math.Sqrt(errors.Select((error) => error * error).Average() - meanError * meanError);
//                        double time = entry.Value.Select((tuple) => tuple.Item3).Average();

//                        double eoc = 0.0;
//                        if (lastMeanError >= 0.0) {
//                            eoc = (Math.Log(meanError) - Math.Log(lastMeanError)) / (Math.Log(entry.Key.Item1) - Math.Log(lastDivision));
//                        }
//                        lastMeanError = meanError;
//                        lastDivision = entry.Key.Item1;

//                        double slope = Regression(
//                            group.Take(index).Select(e => Math.Log(0.4 / e.Key.Item1)).ToArray(),
//                            //group.Take(index).Select(e => Math.Log(0.4 * Math.Pow(2.0, e.Key.Item1))).ToArray(),
//                            group.Take(index).Select(e => Math.Log(e.Value.Average(tuple => tuple.Item2))).ToArray());

//                        log.WriteLine(
//                            "{0}\t{1}\t{2}\t{3}\t{4}",
//                            entry.Key.Item1, // division
//                            (meanError / testCase.Solution).ToString("0.00e-00", NumberFormatInfo.InvariantInfo),
//                            Math.Ceiling(time),
//                            eoc.ToString("0.00", NumberFormatInfo.InvariantInfo),
//                            slope.ToString("0.00", NumberFormatInfo.InvariantInfo));
//                        index++;
//                    }

//                    log.WriteLine();
//                }
//            }

//            globalWatch.Stop();
//            Console.WriteLine("Finished case " + testCase.GetType().Name + " after " + globalWatch.ElapsedMilliseconds + "ms");

//            return dt;
//        }

//        private double PerformVolumeQuadrature(Modes mode, IQuadRuleFactory<QuadRule> factory, SubGrid cutCellGrid, int order, Stopwatch timer) {
//            using (new FuncTrace()) {
//                ScalarFieldQuadrature quadrature;
//                CellQuadratureScheme quadInstr = new CellQuadratureScheme(
//                    factory, cutCellGrid.VolumeMask);
//                if (testCase is ISurfaceTestCase) {
//                    quadrature = new ScalarFieldQuadrature(m_Context, SinglePhaseField, quadInstr, order);
//                } else {
//                    quadrature = new ScalarFieldQuadrature(m_Context, XDGField, quadInstr, order);
//                }

//                timer.Start();
//                quadrature.Execute();
//                timer.Stop();

//                return quadrature.Result;
//            }
//        }

//        private double PerformSurfaceQuadrature(Modes mode, IQuadRuleFactory<QuadRule> volumeFactory, IQuadRuleFactory<QuadRule> edgeFactory, SubGrid cutCellGrid, int order, Stopwatch timer) {
//            using (new FuncTrace()) {
//                CellQuadratureScheme volInstr = new CellQuadratureScheme(
//                    volumeFactory, cutCellGrid.VolumeMask);
//                CellBoundaryQuadratureScheme edgeInstr = new CellBoundaryQuadratureScheme(
//                    new CellBoundaryFromEdgeRuleFactory<CellBoundaryQuadRule>(m_Context, edgeFactory),
//                    cutCellGrid.VolumeMask);

//                ScalarFieldLevelSetIntegrator quadrature = new ScalarFieldLevelSetIntegrator(
//                    levelSetTracker,
//                    SinglePhaseField,
//                    volInstr.Compile(m_Context.GridDat, order),
//                    edgeInstr.Compile(m_Context.GridDat, order),
//                    order,
//                    cutCellGrid,
//                    levelSetIndex);

//                timer.Start();
//                double result = quadrature.ExecuteA().Storage.Sum();
//                timer.Stop();

//                return result;
//            }
//        }

//        private static double RelativeError(double value, double referenceValue) {
//            return Math.Abs(referenceValue - value) / referenceValue;
//        }

//        protected override void PlotCurrentState(double physTime, int timestepNo) {
//            PlotCurrentState(physTime, timestepNo, 0, null);
//        }

//        private void PlotCurrentState(double physTime, int timestepNo, int superSampling, SubGrid subGrid) {
//            Tecplot tecplot = new Tecplot(m_Context, true, false, (uint)superSampling, subGrid);
//            string path = Path.Combine(Path.GetFullPath("."), "plot_" + testCase.GetType().Name);
//            tecplot.PlotFields(path, "bla", physTime, m_IOFields);
//        }

//        private static void ConstantOne(MultidimensionalArray input, MultidimensionalArray output) {
//            for (int i = 0; i < output.Length; i++) {
//                output[i] = 1.0;
//            }
//        }

//        private double Regression(double[] xValues, double[] yValues) {
//            double xAvg = xValues.Average();
//            double yAvg = yValues.Average();

//            double v1 = 0.0;
//            double v2 = 0.0;

//            for (int i = 0; i < yValues.Length; i++) {
//                v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
//                v2 += Math.Pow(xValues[i] - xAvg, 2);
//            }

//            double a = v1 / v2;
//            double b = yAvg - a * xAvg;

//            return a;
//        }
//    }
//}