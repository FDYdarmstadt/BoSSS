using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers.Testing {
    
    /// <summary>
    /// Utility class for executing a series of solver runs and studying the condition number slope over mesh resolution;
    /// Works only for solvers which implemented <see cref="Application{T}.OperatorAnalysis()"/>,
    /// see also <see cref="OpAnalysisBase.GetNamedProperties"/>
    /// </summary>
    public class ConditionNumberScalingTest {

        /// <summary>
        /// Easy-to-use driver routine
        /// </summary>
        /// <param name="controls">
        /// a set of control object over which the scaling is investigated
        /// </param>
        /// <param name="plot">
        /// if true, an interactive Gnuplot session is opened
        /// </param>
        /// <param name="title">
        /// Gnuplot title/output filename
        /// </param>
        /// <param name="ThrowAssertions">
        /// assertions are thrown if the slopes of the condition number are too high, c.f. <see cref="ExpectedSlopes"/>;
        /// should always be true for testing
        /// </param>
        /// <returns>
        /// <see cref="ResultData"/>
        /// </returns>
        static public IDictionary<string, double[]> Perform(IEnumerable<AppControl> controls, bool plot = false, string title = "", bool ThrowAssertions = true) {
            var t = new ConditionNumberScalingTest(title);
            t.SetControls(controls);
            t.RunAndLog();

            t.PrintResults(Console.Out);
            
            if(plot) {
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPIsize);

                string mpiS = "";
                if(MPIsize > 1) {
                    mpiS = "." + MPIrank + "of" + MPIsize;
                }

                using (var gp = t.Plot()) {

                    if(title.IsEmptyOrWhite()) {
                        gp.Execute();
                        Console.WriteLine("plotting in interactive gnuplot session - press any key to continue...");
                        Console.ReadKey();
                    } else {

                        // set terminal
                        int xRes = 1024;
                        int yRes = 768;
                        gp.Terminal = string.Format("pngcairo size {0},{1}", xRes, yRes);

                        string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");
                        gp.OutputFile = title + "-" + DateNtime + mpiS + ".png";

                        // call gnuplot
                        int exCode = gp.RunAndExit(); // run & close gnuplot
                        if(exCode != 0) {
                            Console.WriteLine("Gnuplot-internal error: exit code " + exCode);
                        }

                        // ----------------------------------------
                        using(var tw = new System.IO.StreamWriter(title + "-" + DateNtime + mpiS + ".txt")) {
                            t.PrintResults(tw);
                        }

                    }
                }
            }

            if(ThrowAssertions)
                t.CheckResults();

            return t.ResultData;
        }

        /// <summary>
        /// Interactive plotting using gnuplot.
        /// </summary>
        public Gnuplot.Gnuplot Plot() {
           
            
            var data = this.ResultData;
            if(data == null)
                throw new NotSupportedException("No data available: user must call 'RunAndLog()' first.");

            var gp = new Gnuplot.Gnuplot();

            gp.Cmd("set key t l");

            var fmt = new PlotFormat("rx-");
            fmt.LineWidth = 4;
            fmt.PointSize = 2.0;

            int Kount = 1;

            foreach (var ttt in ExpectedSlopes) {
                double[] xVals = data[ttt.Item1.ToString()];
                string[] allYNames = data.Keys.Where(name => ttt.Item2.WildcardMatch(name)).ToArray();

                foreach(string yName in allYNames) {
                    double[] yVals = data[yName];
                    double Slope = xVals.LogLogRegressionSlope(yVals);

                    gp.PlotXY(xVals, yVals, logX: true, logY: true, title:yName, format:(fmt.WithLineColor(Kount).WithPointType(Kount)));
                    gp.SetXLabel(ttt.Item1.ToString());

                    Kount++;
                }
            }

                       

            return gp;
        }

        string m_title;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="Title">
        /// Optional title used for plots, etc.
        /// </param>
        public ConditionNumberScalingTest(string Title) {
            m_title = Title;

            this.ExpectedSlopes = new List<ValueTuple<XAxisDesignation, string, double, double>>();

            ExpectedSlopes.Add((XAxisDesignation.Grid_1Dres, "TotCondNo-*", 2.4, 1.5));
            ExpectedSlopes.Add((XAxisDesignation.Grid_1Dres, "StencilCondNo-innerUncut-*", 0.5, -0.2));
            ExpectedSlopes.Add((XAxisDesignation.Grid_1Dres, "StencilCondNo-innerCut-*", 0.5, -0.2));
            ExpectedSlopes.Add((XAxisDesignation.Grid_1Dres, "StencilCondNo-bndyUncut-*", 0.5, -0.2));
            ExpectedSlopes.Add((XAxisDesignation.Grid_1Dres, "StencilCondNo-bndyCut-*", 0.5, -0.2));

            ExpectedMaximum.Add(("TotCondNo-*", 1e13));
            ExpectedMaximum.Add(("StencilCondNo-innerUncut-*", 1e6));
            ExpectedMaximum.Add(("StencilCondNo-innerCut-*", 1e6));
            ExpectedMaximum.Add(("StencilCondNo-bndyUncut-*", 1e6));
            ExpectedMaximum.Add(("StencilCondNo-bndyCut-*", 1e6));
        }


        /// <summary>
        /// A range of control objects over which the condition number scaling is performed.
        /// </summary>
        public void SetControls(IEnumerable<AppControl> controls) {
            this.Controls = controls.ToArray();
        }


        /// <summary>
        /// Entsetzlich viel code für was primitives
        /// </summary>
        class MyEnu : IEnumerable<AppControl> {

            Func<IGrid>[] GridFuncs;
            AppControl BaseControl;

            public MyEnu(AppControl __baseControl, Func<IGrid>[] __GridFuncs) {
                this.GridFuncs = __GridFuncs;
                this.BaseControl = __baseControl;
            }

            public IEnumerator<AppControl> GetEnumerator() {
                return new E() { GridFuncs = this.GridFuncs, BaseControl = this.BaseControl };
            }

            IEnumerator IEnumerable.GetEnumerator() {
                return new E() { GridFuncs = this.GridFuncs, BaseControl = this.BaseControl };
            }

            class E : IEnumerator<AppControl> {
                public AppControl BaseControl;
                int position = -1;
                public Func<IGrid>[] GridFuncs;

                public AppControl Current {
                    get {
                        return GetCurrent();
                    }
                }

                private AppControl GetCurrent() {
                    var c = BaseControl;
                    c.GridGuid = default(Guid);
                    c.GridFunc = GridFuncs[position];
                    return c;
                }

                object IEnumerator.Current {
                    get {
                        return GetCurrent();
                    }
                }

                public void Dispose() { }

                public bool MoveNext() {
                    position++;
                    return position < GridFuncs.Length;
                }

                public void Reset() {
                    position = -1;
                }
            }
        }


        /// <summary>
        /// A range of control objects over which the condition number scaling is performed.
        /// </summary>
        /// <param name="BaseControl">
        /// basic settings
        /// </param>
        /// <param name="GridFuncs">
        /// a sequence of grid-generating functions, (<see cref="AppControl.GridFunc"/>),
        /// which defines the grid for each run
        /// </param>
        public void SetControls(AppControl BaseControl, IEnumerable<Func<IGrid>> GridFuncs) {
            Controls = new MyEnu(BaseControl, GridFuncs.ToArray());
        }



        IEnumerable<AppControl> Controls;

        /// <summary>
        /// One tuple for each slope that should be tested
        /// - 1st item: name of x-axis
        /// - 2nd item: name of y-axis (wildcards accepted)
        /// - 3rd item: expected slope in the log-log-regression
        /// - 4th item: lower bound for expected slope in the log-log-regression
        /// </summary>
        public IList<(XAxisDesignation xDesign, string yName, double slopeMax, double slopeMin)> ExpectedSlopes;


        /// <summary>
        /// One tuple for maximum value that should be tested
        /// - 1st item: name of y-axis (wildcards accepted)
        /// - 2nd item: maximum for the respective condition number
        public IList<(string yName, double MaxCondNo)> ExpectedMaximum;


        /// <summary>
        /// Phase 2, Examination: prints slope thresholds to console output. 
        /// </summary>
        public virtual void PrintResults(System.IO.TextWriter tw) {

            var data = this.ResultData;
            if(data == null)
                throw new NotSupportedException("No data available: user must call 'RunAndLog()' first.");

            tw.WriteLine("Condition Number Scaling Test slopes:");

            IDictionary<string, IEnumerable<double>> testData = new Dictionary<string, IEnumerable<double>>();

            foreach (var ttt in ExpectedSlopes) {
                double[] xVals = data[ttt.Item1.ToString()];
                string[] allYNames = data.Keys.Where(name => ttt.yName.WildcardMatch(name)).ToArray();

                if (!testData.ContainsKey(ttt.Item1.ToString())) 
                    testData.Add(ttt.Item1.ToString(), xVals);

                foreach(string yName in allYNames) {
                    double[] yVals = data[yName];
                    double Slope = xVals.LogLogRegressionSlope(yVals);

                    testData.Add(yName, yVals);

                    string tstPasses = Slope <= ttt.slopeMax ? Slope >= ttt.slopeMin ? "passed" : $"FAILED (threshold is {ttt.slopeMin})" : $"FAILED (threshold is {ttt.slopeMax})";
                    tw.WriteLine($"    Slope for {yName}: {Slope:0.###e-00} -- {tstPasses}");
                }
            }

            foreach (var ttt in ExpectedMaximum) {
                string[] allYNames = data.Keys.Where(name => ttt.yName.WildcardMatch(name)).ToArray();

                foreach (string yName in allYNames) {
                    double[] yVals = data[yName];

                    double max = yVals.Max();

                    string tstPasses = max <= ttt.MaxCondNo ? "passed" : $"FAILED (threshold is {ttt.MaxCondNo})";
                    tw.WriteLine($"    Max. for {yName}: {max:0.###e-00} -- {tstPasses}");
                }
            }



            CSVFile.SaveToCSVFile<IEnumerable<double>>(testData, "ConditionNumberScalingTest_dataSet-" + DateTime.Now.ToString("yyyyMMMdd_HHmmss") + ".txt");
            //Console.WriteLine("warning no output-file - ToDo");

        }


        /// <summary>
        /// Phase 2, Examination: checks slope thresholds with NUnit assertions. 
        /// </summary>
        public virtual void CheckResults() {

            var data = this.ResultData;
            if(data == null)
                throw new NotSupportedException("No data available: user must call 'RunAndLog()' first.");

         
            foreach (var ttt in ExpectedSlopes) {
                double[] xVals = data[ttt.Item1.ToString()];
                string[] allYNames = data.Keys.Where(name => ttt.yName.WildcardMatch(name)).ToArray();

                foreach(string yName in allYNames) {
                    double[] yVals = data[yName];

                    double Slope = DoubleExtensions.LogLogRegressionSlope(xVals, yVals);

                    Assert.LessOrEqual(Slope, ttt.slopeMax, $"Condition number slope for {ttt.yName} to high; at max. {ttt.slopeMax}");
                    Assert.LessOrEqual(Slope, ttt.slopeMin, $"Condition number slope for {ttt.yName} to low; at min. {ttt.slopeMin}");
                }
            }

            foreach (var ttt in ExpectedMaximum) {
                string[] allYNames = data.Keys.Where(name => ttt.yName.WildcardMatch(name)).ToArray();

                foreach (string yName in allYNames) {
                    double[] yVals = data[yName];

                    double max = yVals.Max();

                    Assert.LessOrEqual(max, ttt.MaxCondNo, $"Condition number maximum for {ttt.yName} to high; at max. {ttt.MaxCondNo}");
                }
            }

        }


        /// <summary>
        /// Stores the result of <see cref="RunAndLog"/>; Contains
        /// a table, containing grid resolutions and measurements on condition number
        /// - keys: column names
        /// - values: measurements of each column
        /// </summary>
        public IDictionary<string, double[]> ResultData {
            get;
            private set;
        }

        /// <summary>
        /// 
        /// </summary>
        static public int RunNumber;

        /// <summary>
        /// Phase 1: runs the solvers and stores results in <see cref="ResultData"/>.
        /// Utility routine, performs an operator analysis on a sequence of control objects and returns a table of results.
        /// </summary>
        public void RunAndLog() {

            var ret = new Dictionary<string, List<double>>();
           

            int Counter = 0;
            foreach(var C in this.Controls) {
                var st = C.GetSolverType();

                Counter++;
                Console.WriteLine("================================================================");
                Console.WriteLine($"Condition Number Scaling Analysis:  Run {Counter} of {this.Controls.Count()}");
                Console.WriteLine("================================================================");
                RunNumber = Counter;

                using(var solver = (BoSSS.Solution.IApplication)Activator.CreateInstance(st)) {
                    Console.WriteLine("  Starting Solver...");
                    solver.Init(C);
                    solver.RunSolverMode();

                    Console.WriteLine("  Done solver; Now performing operator analysis...");

                    int J = Convert.ToInt32(solver.CurrentSessionInfo.KeysAndQueries["Grid:NoOfCells"]);
                    double hMin = Convert.ToDouble(solver.CurrentSessionInfo.KeysAndQueries["Grid:hMin"]);
                    double hMax = Convert.ToDouble(solver.CurrentSessionInfo.KeysAndQueries["Grid:hMax"]);
                    int D = Convert.ToInt32(solver.CurrentSessionInfo.KeysAndQueries["Grid:SpatialDimension"]);
                    double J1d = Math.Pow(J, 1.0 / D);

                    var prop = solver.OperatorAnalysis();
                    Console.WriteLine("  finished analysis.");


                    if(ret.Count == 0) {
                        ret.Add(XAxisDesignation.Grid_NoOfCells.ToString(), new List<double>());
                        ret.Add(XAxisDesignation.Grid_hMin.ToString(), new List<double>());
                        ret.Add(XAxisDesignation.Grid_hMax.ToString(), new List<double>());
                        ret.Add(XAxisDesignation.Grid_1Dres.ToString(), new List<double>());

                        foreach(var kv in prop) {
                            ret.Add(kv.Key, new List<double>());
                        }
                    }

                    {
                        ret[XAxisDesignation.Grid_NoOfCells.ToString()].Add(J);
                        ret[XAxisDesignation.Grid_hMin.ToString()].Add(hMin);
                        ret[XAxisDesignation.Grid_hMax.ToString()].Add(hMax);
                        ret[XAxisDesignation.Grid_1Dres.ToString()].Add(J1d);

                        foreach(var kv in prop) {
                            ret[kv.Key].Add(kv.Value);
                        }
                    }
                }
            }

            /*
            // write statistics
            // ================
            {
                var xdes = XAxisDesignation.Grid_1Dres.ToString();
                var xVals = ret[xdes];

                Console.WriteLine("Regression of condition number slopes:");
                foreach(string ydes in ret.Keys) {
                    if(!Enum.TryParse(ydes, out XAxisDesignation dummy)) {
                        var yVals = ret[ydes];

                        double slope = DoubleExtensions.LogLogRegression(xVals, yVals);
                        Console.WriteLine($"   slope of {ydes}: {slope}");

                    }
                }
                


            }
            */

            // data conversion & return 
            // ========================
            {
                var realRet = new Dictionary<string, double[]>();
                foreach(var kv in ret) {
                    realRet.Add(kv.Key, kv.Value.ToArray());
                }

                this.ResultData = realRet;
            }
        }


        

        /// <summary>
        /// Names for the x-axis, over which condition number scaling slopes are computed.
        /// </summary>
        public enum XAxisDesignation {

            /// <summary>
            /// total number of cells, <see cref="IGridData.CellPartitioning"/>
            /// </summary>
            Grid_NoOfCells,

            /// <summary>
            /// maximum cell size, <see cref="IGeometricalCellsData.h_min"/>
            /// </summary>
            Grid_hMin,

            /// <summary>
            /// maximum cell size, <see cref="IGeometricalCellsData.h_max"/>
            /// </summary>
            Grid_hMax,

            /// <summary>
            /// D-th root of number of cells, where D is the spatial dimension
            /// </summary>
            Grid_1Dres
        }

        /*
        private static double LogLogRegression(IEnumerable<double> _xValues, IEnumerable<double> _yValues) {
            double[] xValues = _xValues.Select(x => Math.Log10(x)).ToArray();
            double[] yValues = _yValues.Select(y => Math.Log10(y)).ToArray();

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
        */
        

    }
}
