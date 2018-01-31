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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace BoSSS.Solution.Gnuplot {

    /// <summary>
    /// Convenience interface to create gnuplot plots from BoSSS
    /// </summary>
    public class Gnuplot : IDisposable {

        private Process gnucmd;

        private string m_sGNUPlotFileName;

        private int MultiplotCols = 1;

        private int MultiplotRows = 1;

        private IFormatProvider nfoi = NumberFormatInfo.InvariantInfo;

        private List<string> m_TempFiles = new List<string>();

        /// <summary>
        /// Path of all temporary created data files
        /// </summary>
        public IEnumerable<string> TempFilesPath {
            get {
                return m_TempFiles.AsEnumerable();
            }
        }

        PlotFormat m_baseLineFormat;

        /// <summary>
        /// Gnuplot 'set output' option; 
        /// </summary>
        public string OutputFile {
            get;
            set;
        }

        /// <summary>
        /// Gnuplot 'set terminal' option; 
        /// </summary>
        public string Terminal {
            get;
            set;
        }

        StringWriter m_Commands = new StringWriter();

        /// <summary>
        /// constructor;
        /// </summary>
        /// <param name="AlternativeGnuplotPath">
        /// The path where the gnuplot executeable can be found; if null, the system PATH is searched for the gnuplot exe.
        /// </param>
        /// <param name="baseLineFormat"></param>
        public Gnuplot(string AlternativeGnuplotPath = null, PlotFormat baseLineFormat = default(PlotFormat)) {
            this.m_baseLineFormat = baseLineFormat;
            m_sGNUPlotFileName = GetGnuplotPath(AlternativeGnuplotPath);
            Console.WriteLine("Using gnuplot: " + m_sGNUPlotFileName);
        }

        /// <summary>
        /// Variant of <see cref="Cmd(string)"/> using a format string via
        /// <see cref="String.Format(string, object[])"/>
        /// </summary>
        /// <param name="fmt"></param>
        /// <param name="args"></param>
        public void Cmd(string fmt, params object[] args) {
            this.Cmd(string.Format(fmt, args));
        }

        /// <summary>
        /// Adds a custom gnuplot format (i.e., one which is not exposed by
        /// this interface) to the list of commands to be executed
        /// </summary>
        /// <param name="cmdstr"></param>
        public void Cmd(string cmdstr) {
            this.m_Commands.WriteLine(cmdstr);
            this.m_Commands.Flush();
        }

        private void CreateBoSSSGrid(IGridData c, int upsampling, out NodeSet localNodes, out MultidimensionalArray globalNodes) {
            RefElement simplex = null;
            int[] iArr1 = null;
            if (c.SpatialDimension != 1)
                throw new NotSupportedException("only 1D-domains plots are supported yet.");
            if (upsampling < 1)
                throw new ArgumentOutOfRangeException("upsampling must be at leas 1.");
            int NN = upsampling;
            simplex = c.iGeomCells.GetRefElement(0);

            localNodes = new NodeSet(simplex, NN, 1);
            double left = simplex.Vertices[0, 0];
            double right = simplex.Vertices[1, 0];
            if (NN >= 2) {
                double dx = (right - left) / ((double)NN - 1.0);
                if (dx <= 0.0)
                    throw new ApplicationException("unsupported grid: left >= right ???");
                for (int i = 0; i < NN; i++)
                    localNodes[i, 0] = left + dx * ((double)i);

            } else {
                localNodes[0, 0] = left + ((right - left) * 0.5);
            }
            localNodes.LockForever();
            int i3 = c.iLogicalCells.NoOfLocalUpdatedCells;
            iArr1 = new int[] { i3, NN, 1 };
            globalNodes = MultidimensionalArray.Create(iArr1);
            c.TransformLocal2Global(localNodes, 0, i3, globalNodes, 0);
        }

        private void DeleteTemporaryFiles() {
            foreach (string s in m_TempFiles) {
                try {
                    File.Delete(s);
#if DEBUG
                } catch (Exception e) {
                    Console.Error.WriteLine(e.Message);
                }
#else 
                } catch (Exception) {
                    // Swallow
                }
#endif

            }
            m_TempFiles.Clear();
        }

        /// <summary>
        /// Path to Gnuplot executable.
        /// </summary>
        /// <param name="alternativeGnuplotPath">Alternative search path, can be null.</param>
        /// <returns></returns>
        public static string GetGnuplotPath(string alternativeGnuplotPath = null) {
            string path = GetProgramPath("gnuplot.exe", alternativeGnuplotPath); // Windows gnuplot
            if (path == null)
                path = GetProgramPath("gnuplot", alternativeGnuplotPath); // Unix

            if (path == null) {
                throw new ApplicationException("Can't find gnuplot in your PATH");
            }

            return path;
        }

        /// <summary>
        /// A hack: on Windows, we
        /// prefer gnuplot which is shiped with BoSSS, e.g. in
        /// $BOSSS_INSTALL\bin\native\win\gnuplot-gp510-20160418-win32-mingw\gnuplot\bin
        /// </summary>
        /// <returns></returns>
        private static string GetBoSSSGnuplot() {
            if (System.Environment.OSVersion.Platform != PlatformID.Win32NT)
                return null;

            string bi = System.Environment.GetEnvironmentVariable("BOSSS_INSTALL");
            if (bi == null || bi.Length < 0)
                return null;

            string gpDir = Path.Combine(bi, "bin", "native", "win");
            gpDir = Directory.EnumerateDirectories(gpDir, "gnuplot*").FirstOrDefault();
            if (gpDir == null)
                return null;
            gpDir = Path.Combine(gpDir, "gnuplot", "bin");

            return gpDir;
        }

        private static string GetProgramPath(string pname, string alternativeGnuplotPath) {

            string path = System.Environment.GetEnvironmentVariable("PATH");
            if (path == null)
                path = "";

            List<string> pathDirs;
            if (alternativeGnuplotPath != null) {
                pathDirs = (new string[] { alternativeGnuplotPath }).ToList();
            } else {
                string WinPath = GetBoSSSGnuplot();
                pathDirs = path.Split(new[] { Path.PathSeparator }, StringSplitOptions.RemoveEmptyEntries).ToList();
                if (WinPath != null)
                    pathDirs.Insert(0, WinPath);
            }

            foreach (string s4 in pathDirs) {
                try {
                    string gnuplot_path = Path.Combine(s4, pname);
                    if (File.Exists(gnuplot_path))
                        return gnuplot_path;
                } catch (Exception) {
                }
            }

            return null;
        }

        /// <summary>
        /// Simplified variant of
        /// <see cref="Plot3D(double[], double[], double[,], string, bool, bool)"/>
        /// where x-axis and y-axis data are omitted
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="title"></param>
        /// <param name="Colored"></param>
        /// <param name="Hidden3D"></param>
        public void Plot3D(double[,] Z, string title, bool Colored, bool Hidden3D) {
            double[] x = new double[Z.GetLength(1)];
            for (int i = 0; i < x.Length; i++)
                x[i] = (double)i;

            double[] y = new double[Z.GetLength(0)];
            for (int i = 0; i < y.Length; i++)
                y[i] = (double)i;

            Plot3D(x, y, Z, title, Colored, Hidden3D);
        }

        /// <summary>
        /// Creates a 3D plot
        /// </summary>
        /// <param name="xScale"></param>
        /// <param name="yScale"></param>
        /// <param name="Z"></param>
        /// <param name="title"></param>
        /// <param name="Colored"></param>
        /// <param name="Hidden3D"></param>
        public void Plot3D(double[] xScale, double[] yScale, double[,] Z, string title, bool Colored, bool Hidden3D) {
            // checks
            if (xScale.Length != Z.GetLength(1))
                throw new ArgumentException("the xScale must have as much entries as Z coloums");
            if (yScale.Length != Z.GetLength(0))
                throw new ArgumentException("the yScale must have as much entries as Z rows");

            // create temp file
            String name = Path.GetTempFileName();

            FileStream f = new FileStream(name, FileMode.OpenOrCreate);
            StreamWriter s = new StreamWriter(f);

            // Save the temporary filename
            m_TempFiles.Add(name);

            for (int i = 0; i < yScale.Length; i++) {
                for (int j = 0; j < Z.GetLength(1); j++) {
                    s.Write(yScale[i].ToString(nfoi));
                    s.Write(" ");
                    s.Write(xScale[j].ToString(nfoi));
                    s.Write(" ");
                    s.Write(Z[i, j].ToString(nfoi));
                    s.Write(System.Environment.NewLine);
                }
            }

            s.Close();
            f.Close();

            // set grid
            int M = yScale.Length, N = xScale.Length;

            Cmd(
                String.Concat("set dgrid3d ",
                String.Concat(N.ToString(nfoi),
                String.Concat(",",
                String.Concat(M.ToString(nfoi), ",2")))));

            if (Colored)
                Cmd("set pm3d");
            else
                Cmd("unset pm3d");

            if (Hidden3D)
                Cmd("set hidden3d");
            else
                Cmd("unset hidden3d");

            if (title != null)
                Cmd(String.Concat("set title \"", title, "\""));

            Cmd("set nokey");

            String gpnmn = name.Replace(Path.DirectorySeparatorChar, '/');
            Cmd(String.Concat("splot \"",
                 String.Concat(gpnmn, "\" with lines")));
        }

        /// <summary>
        /// Variant of
        /// <see cref="PlotContour(double[], double[], double[,], string, bool, double[], bool)"/>
        /// where x-axis data and y-axis data is deduced from <paramref name="Z"/>
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="title"></param>
        /// <param name="Colored"></param>
        /// <param name="Levels"></param>
        /// <param name="hidden3D"></param>
        public void PlotContour(double[,] Z, string title, bool Colored, int Levels, bool hidden3D = false) {
            double[] x = new double[Z.GetLength(1)];
            for (int i = 0; i < x.Length; i++)
                x[i] = (double)i;

            double[] y = new double[Z.GetLength(0)];
            for (int i = 0; i < y.Length; i++)
                y[i] = (double)i;
            PlotContour(x, y, Z, title, Colored, Levels, hidden3D);
        }

        /// <summary>
        /// Plots <paramref name="Levels"/> 3D contour surfaces of the data
        /// provided in <paramref name="Z"/>
        /// </summary>
        /// <param name="xScale"></param>
        /// <param name="yScale"></param>
        /// <param name="Z"></param>
        /// <param name="title"></param>
        /// <param name="Colored"></param>
        /// <param name="Levels"></param>
        /// <param name="hidden3D"></param>
        public void PlotContour(double[] xScale, double[] yScale, double[,] Z, string title, bool Colored, int Levels, bool hidden3D = false) {
            Cmd("set contour");
            Cmd("set nosurface");
            Cmd("set view 0, 0");
            Cmd(String.Concat("set cntrparam levels ", Levels.ToString()));
            Cmd("set nokey");

            Plot3D(xScale, yScale, Z, title, Colored, hidden3D);
        }

        /// <summary>
        /// Plots the 3D contours of the data provided in <paramref name="Z"/>
        /// at the given <paramref name="Levels"/>s
        /// </summary>
        /// <param name="xScale"></param>
        /// <param name="yScale"></param>
        /// <param name="Z"></param>
        /// <param name="title"></param>
        /// <param name="Colored"></param>
        /// <param name="Levels"></param>
        /// <param name="hidden3D"></param>
        public void PlotContour(double[] xScale, double[] yScale, double[,] Z, string title, bool Colored, double[] Levels, bool hidden3D = false) {
            Cmd("set contour");
            Cmd("set nosurface");
            Cmd("set view map");
            Cmd("set contour base");
            Cmd(String.Concat("set cntrparam levels discrete ", Levels.Skip(1).Aggregate("" + Levels.First(), (a, b) => a + ", " + b)));
            Cmd("set nokey");

            Plot3D(xScale, yScale, Z, title, Colored, hidden3D);
        }

        /// <summary>
        /// Plots the 2D contours of <paramref name="function"/> at the given
        /// <paramref name="levels"/>
        /// </summary>
        /// <param name="xNodes"></param>
        /// <param name="yNodes"></param>
        /// <param name="function"></param>
        /// <param name="levels"></param>
        /// <param name="title"></param>
        public void PlotContour(double[] xNodes, double[] yNodes, Func<double, double, double> function, double[] levels, string title = "") {
            double[,] z = new double[xNodes.Length, yNodes.Length];
            for (int i = 0; i < xNodes.Length; i++) {
                for (int j = 0; j < yNodes.Length; j++) {
                    z[i, j] = function(xNodes[i], yNodes[j]);
                }
            }

            PlotContour(
                yNodes,
                xNodes,
                z,
                title,
                false,
                levels,
                true);
        }

        /// <summary>
        /// Plots the given <paramref name="equation"/> which has to be in a
        /// format understood by gnuplot
        /// </summary>
        /// <param name="equation"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        /// <param name="useY2"></param>
        /// <param name="useX2"></param>
        public void PlotEquation(String equation, String title = null, PlotFormat format = null, bool deferred = true, bool useX2 = false, bool useY2 = false) {
            using (StringWriter stw = new StringWriter()) {
                if (deferred) {
                    // nop
                } else {
                    stw.Write("plot ");
                }

                stw.Write("plot ");
                stw.Write(equation);
                stw.Write(" ");
                stw.Write(Format2D(title, format, useX2, useY2));

                if (deferred) {
                    this.m_DeferredPlotCommands.Add(stw.ToString());
                } else {
                    Cmd(stw.ToString());
                }
            }
        }

        /// <summary>
        /// Plots a one-dimensional <see cref="DGField"/> <paramref name="u"/>.
        /// </summary>
        /// <param name="u"></param>
        /// <param name="format"></param>
        /// <param name="upsampling"></param>
        /// <param name="deferred"></param>
        /// <param name="label"></param>
        public void PlotField(DGField u, PlotFormat format = null, int upsampling = 1, bool deferred = true, string label = null) {
            PlotField(u.GridDat, u.Evaluate, Label: label ?? u.Identification, upsampling: upsampling, format: format, deferred: deferred);
        }

        /// <summary>
        /// Plots a one-dimensional <see cref="ScalarFunctionEx"/> <paramref name="func"/>.
        /// </summary>
        /// <param name="c"></param>
        /// <param name="func"></param>
        /// <param name="Label"></param>
        /// <param name="upsampling"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotField(
            IGridData c, ScalarFunctionEx func, string Label, int upsampling = 10, PlotFormat format = null, bool deferred = true) {

            // create eval grid
            NodeSet nodesLocal;
            MultidimensionalArray nodesGlobal;
            CreateBoSSSGrid(c, upsampling,
                out nodesLocal,
                out nodesGlobal);
            int NN = nodesLocal.NoOfNodes;

            // evaluate
            int NoofCells = c.iLogicalCells.NoOfLocalUpdatedCells;
            MultidimensionalArray y = MultidimensionalArray.Create(new int[] { NoofCells, NN });
            func(0, NoofCells, nodesLocal, y);

            // plot
            this.PlotXY(nodesGlobal.Storage, y.Storage, title: Label, format: format, deferred: deferred);
        }

        /// <summary>
        /// Plots a one-dimensional <see cref="ScalarFunction"/> <paramref name="func"/>.
        /// </summary>
        /// <param name="func"></param>
        /// <param name="c"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotFunction(
            ScalarFunction func, GridData c, string title = null, PlotFormat format = null, bool deferred = true) {

            // create eval grid
            NodeSet nodes;
            MultidimensionalArray nodesGlobal;
            CreateBoSSSGrid(c, 2,
                out nodes,
                out nodesGlobal);

            int NN = nodes.GetLength(0);

            // Evaluate
            MultidimensionalArray y;
            y = MultidimensionalArray.Create(new int[] { nodesGlobal.GetLength(0) * NN });

            MultidimensionalArray globalNodesRsz;
            int[] sz = new int[2];
            sz[0] = y.GetLength(0);
            sz[1] = nodes.GetLength(1);
            globalNodesRsz = nodesGlobal.ResizeShallow(sz);

            func.Invoke(globalNodesRsz, y);

            // plot
            this.PlotXY(nodesGlobal.Storage, y.Storage, title ?? func.ToString(), format, deferred);
        }

        /// <summary>
        /// plots the affine-linear function y = x * <paramref name="slope"/> + <paramref name="intercept"/>
        /// </summary>
        public void PlotSlope(
            double slope, double intercept, string title = null, PlotFormat format = null, bool deferred = true, bool useX2 = false, bool useY2 = false) {

            using (StringWriter stringWriter = new StringWriter()) {
                if (deferred) {
                    // nop
                } else {
                    stringWriter.Write("plot ");
                }

                stringWriter.Write(" ( ");
                stringWriter.Write(slope.ToString(nfoi));
                stringWriter.Write("*x + ");
                stringWriter.Write((intercept).ToString(nfoi));
                stringWriter.Write(" ) ");
                stringWriter.Write(Format2D(title, format, useX2, useY2));

                if (deferred) {
                    this.m_DeferredPlotCommands.Add(stringWriter.ToString());
                } else {
                    Cmd(stringWriter.ToString());
                }
            }
        }

        /// <summary>
        /// plots the affine-linear function
        /// \f$ y = x ^ <paramref name="expo"/> * <paramref name="C"/> \f$
        /// </summary>
        public void PlotPow(double expo, double C, string title = null,
             PlotFormat format = null, bool deferred = true, bool useX2 = false, bool useY2 = false) {
            using (StringWriter stringWriter = new StringWriter()) {
                if (deferred) {
                    // nop
                } else {
                    stringWriter.Write("plot ");
                }

                stringWriter.Write(" ( (");
                stringWriter.Write(C.ToString(nfoi));
                stringWriter.Write(")*x**( ");
                stringWriter.Write((expo).ToString(nfoi));
                stringWriter.Write(" ) ) ");
                stringWriter.Write(Format2D(title, format, useX2, useY2));


                if (deferred) {
                    this.m_DeferredPlotCommands.Add(stringWriter.ToString());
                } else {
                    Cmd(stringWriter.ToString());
                }
            }
        }

        /// <summary>
        /// Plots the logarithm of the absolute difference between
        /// <paramref name="u"/> and the given
        /// <paramref name="exactSolution"/>
        /// </summary>
        /// <param name="u"></param>
        /// <param name="exactSolution"></param>
        /// <param name="title"></param>
        /// <param name="upsampling"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotLogError(
            DGField u, Func<double[], double> exactSolution, string title = null, int upsampling = 10, PlotFormat format = null, bool deferred = true) {

            // create eval grid
            NodeSet nodesLocal;
            MultidimensionalArray nodesGlobal;
            CreateBoSSSGrid(u.GridDat, upsampling, out nodesLocal, out nodesGlobal);

            // evaluate
            int noOfNodesPerCell = nodesLocal.NoOfNodes;
            int NoofCells = u.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            MultidimensionalArray y = MultidimensionalArray.Create(new int[] { NoofCells, noOfNodesPerCell });
            u.Evaluate(0, NoofCells, nodesLocal, y);

            for (int i = 0; i < NoofCells; i++) {
                for (int j = 0; j < noOfNodesPerCell; j++) {
                    double[] node = nodesGlobal.ExtractSubArrayShallow(i, j, -1).To1DArray();
                    y[i, j] = (Math.Abs(y[i, j] - exactSolution(node)));
                }
            }

            // plot
            this.PlotXY(
                nodesGlobal.Storage,
                y.Storage,
                title: title,
                format: format,
                deferred: deferred,
                logY: true);

        }

        private string Format2D(string title, PlotFormat format, bool useX2, bool useY2) {
            if (format == null)
                format = new PlotFormat(this.m_baseLineFormat);

            using (StringWriter stringWriter = new StringWriter()) {
                if ((title != null) && !title.Equals(String.Empty)) {
                    stringWriter.Write("title \"");
                    stringWriter.Write(title);
                    stringWriter.Write("\"");
                } else {
                    stringWriter.Write(" notitle");
                }
                if(useX2 || useY2) {
                    stringWriter.Write(" axes ");
                    stringWriter.Write(useX2 ? "x2" : "x1");
                    stringWriter.Write(useY2 ? "y2" : "y1");
                    stringWriter.Write(" ");
                }

                stringWriter.Write(" with ");
                stringWriter.Write(format.Style.ToString().ToLower());

                string LineColorString;
                if(Enum.IsDefined(typeof(LineColors), format.LineColor)) {
                    LineColorString = "\"" + format.LineColor.ToString().ToLowerInvariant() + "\"";
                } else {
                    LineColorString = format.LineColor.ToString();
                }

                stringWriter.Write(" linecolor  " + LineColorString);
                switch (format.Style) {
                    case Styles.Lines:
                    case Styles.LinesPoints:
                        stringWriter.Write(" dashtype " + (int)format.DashType);
                        if (format.DashType != DashTypes.Solid)
                            Cmd("set termoption dashed");
                        stringWriter.Write(" linewidth " + format.LineWidth.ToString(nfoi));
                        break;
                }

                switch (format.Style) {
                    case Styles.Errorbars:
                    case Styles.LinesPoints:
                    case Styles.Points:
                        stringWriter.Write(" pointtype " + (int)format.PointType);
                        stringWriter.Write(" pointsize " + format.PointSize.ToString(nfoi));
                        break;
                }

                string R = stringWriter.ToString();
                return R;
            }
        }

        /// <summary>
        /// Plots <paramref name="y"/> using
        /// <see cref="PlotXY(IEnumerable{double}, IEnumerable{double}, string, PlotFormat, bool, bool?, bool?, bool, bool)"/>
        /// where the x-axis data is deduced from <paramref name="y"/>
        /// </summary>
        /// <param name="y"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotY(IEnumerable<double> y, string title = null, PlotFormat format = null, bool deferred = true) {
            PlotXY(null, y, title, format, deferred);
        }

        /// <summary>
        /// Standard plot of <paramref name="x"/> vs. <paramref name="y"/>
        /// using linear scaling.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred">
        /// The actual plot command is deferred until <see cref="WriteDeferredPlotCommands"/> (or one of its callers, e.g. <see cref="RunAndExit"/>)
        /// is called.
        /// </param>
        /// <param name="logX">
        /// Tri-state:
        /// - true: 'logscale' for the x-axis is set.
        /// - false: 'logscale' for the x-axis is unset.
        /// - null: no change to the 'logscale' state is done.
        /// </param>
        /// <param name="logY">
        /// Analogous to <paramref name="logY"/>.
        /// </param>
        /// <param name="useX2">
        /// Use secondary x axis
        /// </param>
        /// <param name="useY2">
        /// Use secondary y axis
        /// </param>
        public void PlotXY(IEnumerable<double> x, IEnumerable<double> y, string title = null, PlotFormat format = null, bool deferred = true, bool? logX = null, bool? logY = null, bool useX2 = false, bool useY2 = false) {
            using (StringWriter stringWriter = new StringWriter()) {
                string s2 = null;

                if ((x != null) && (x.Count() != y.Count()))
                    throw new ArgumentException("x and y Vector must have the same lengths.");
                string s1 = Path.GetTempFileName();

                FileStream fileStream = new FileStream(s1, FileMode.Open);
                StreamWriter streamWriter = new StreamWriter(fileStream);
                m_TempFiles.Add(s1);
                if (x == null) {
                    var _y = y.ToArray();

                    for (int i = 0; i < _y.Length; i++) {
                        streamWriter.Write(i);
                        streamWriter.Write(" ");
                        streamWriter.WriteLine(_y[i].ToString(nfoi));
                    }
                } else {
                    var _x = x.ToArray();
                    var _y = y.ToArray();

                    for (int i = 0; i < _x.Length; i++) {
                        streamWriter.Write(_x[i].ToString(nfoi));
                        streamWriter.Write(" ");
                        streamWriter.WriteLine(_y[i].ToString(nfoi));
                    }
                }
                streamWriter.Close();
                fileStream.Close();

                if (deferred) {
                    // nop
                } else {
                    stringWriter.Write("plot ");
                }

                s2 = s1.Replace(Path.DirectorySeparatorChar, '/');
                stringWriter.Write("\"");
                stringWriter.Write(s2);
                stringWriter.Write("\" ");
                stringWriter.Write(Format2D(title, format, useX2, useY2));

                if (logX != null) {
                    if (logX.Value == true)
                        Cmd("set logscale x" + (useX2 ? "2" : ""));
                    else
                        Cmd("unset logscale x"  + (useX2 ? "2" : ""));
                }

                if (logY != null) {
                    if (logY.Value == true)
                        Cmd("set logscale y" + (useY2 ? "2" : ""));
                    else
                        Cmd("unset logscale y" + (useY2 ? "2" : ""));
                }


                if (deferred) {
                    m_DeferredPlotCommands.Add(stringWriter.ToString());
                } else {
                    Cmd(stringWriter.ToString());
                }
            }
        }

        /// <summary>
        /// Plots <paramref name="x"/> vs.
        /// <paramref name="yFunc"/>(<paramref name="x"/>)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="yFunc"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotXY(IEnumerable<double> x, Func<double[], double> yFunc, string title = null, PlotFormat format = null, bool deferred = true) {
            PlotXY(
                x,
                x.Select(x_i => yFunc(new double[] { x_i })),
                title,
                format,
                deferred);
        }

        /// <summary>
        /// Plots <paramref name="x"/> vs. <paramref name="y"/> using a
        /// logarithmic scaling for <paramref name="x"/>
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotLogXY(IEnumerable<double> x, IEnumerable<double> y, string title = null, PlotFormat format = null, bool deferred = true) {
            PlotXY(
                x,
                y,
                title,
                format,
                deferred,
                true,
                false);
        }

        /// <summary>
        /// Plots <paramref name="x"/> vs. <paramref name="y"/> using a
        /// logarithmic scaling for <paramref name="y"/>
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotXLogY(IEnumerable<double> x, IEnumerable<double> y, string title = null, PlotFormat format = null, bool deferred = true) {
            PlotXY(
                x,
                y,
                title,
                format,
                deferred,
                false,
                true);
        }

        /// <summary>
        /// Plots<paramref name="x"/> vs. <paramref name= "y" /> using a
        /// logarithmic scaling for <paramref name="x"/> and
        /// <paramref name="y"/>
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="title"></param>
        /// <param name="format"></param>
        /// <param name="deferred"></param>
        public void PlotLogXLogY(IEnumerable<double> x, IEnumerable<double> y, string title = null,
            PlotFormat format = null, bool deferred = true) {
            PlotXY(
                x,
                y,
                title,
                format,
                deferred,
                true,
                true);
        }

        /// <summary>
        /// Plots data stored in a file called <paramref name="fileName"/>
        /// </summary>
        public void PlotDataFile(string fileName, string title = null, PlotFormat format = null, bool deferred = true, bool useX2 = false, bool useY2 = false) {
            using (StringWriter stringWriter = new StringWriter()) {
                if (!deferred) {
                    stringWriter.Write("plot ");
                }

                stringWriter.Write("\"");
                stringWriter.Write(
                    fileName.Replace(Path.DirectorySeparatorChar, '/'));
                stringWriter.Write("\" ");
                stringWriter.Write(Format2D(title, format, useX2, useY2));

                if (deferred) {
                    m_DeferredPlotCommands.Add(stringWriter.ToString());
                } else {
                    Cmd(stringWriter.ToString());
                }
            }
        }

        /// <summary>
        /// Plots the matrix structure of the given <paramref name="matrix"/>.
        /// That is, plots a dot a for each non-zero entry in the matrix, where
        /// red dots are used for diagonal blocks and green dots are used for
        /// off-diagonal blocks. The blocking structure is defined by the user
        /// via the given <paramref name="blockSize"/>
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="blockSize"></param>
        public void PlotMatrixStructure(IMatrix matrix, int blockSize = 1) {
            if (matrix.NoOfRows % blockSize != 0) {
                throw new Exception();
            }
            int rowBlocks = matrix.NoOfRows / blockSize;

            if (matrix.NoOfCols % blockSize != 0) {
                throw new Exception();
            }
            int colBlocks = matrix.NoOfCols / blockSize;

            for (int rowBlock = 0; rowBlock < rowBlocks; rowBlock++) {
                for (int colBlock = 0; colBlock < colBlocks; colBlock++) {
                    LineColors lineColor = LineColors.Green;
                    if (rowBlock == colBlock) {
                        lineColor = LineColors.Red;
                    }

                    for (int i = 0; i < blockSize; i++) {
                        int row = rowBlock * blockSize + i;

                        for (int j = 0; j < blockSize; j++) {
                            int col = colBlock * blockSize + j;

                            if (Math.Abs(matrix[row, col]) > 1e-14) {
                                PlotFormat format = new PlotFormat(
                                    lineColor: lineColor,
                                    dashType: DashTypes.Dotted,
                                    pointSize: 1,
                                    pointType: PointTypes.Circle,
                                    Style: Styles.Points);
                                PlotXY(new double[] { row }, new double[] { col }, format: format);
                            }
                        }
                    }
                }
            }
        }

        List<string> m_DeferredPlotCommands = new List<string>();

        /// <summary>
        /// Sets up a multiplot environment that allows for including multiple
        /// plots (<paramref name="rows"/> rows and <paramref name="cols"/>
        /// columns) in a single graphic. Individual sub-plots are selected
        /// via <see cref="SetSubPlot(int, int, string)"/>.
        /// </summary>
        /// <param name="rows"></param>
        /// <param name="cols"></param>
        public void SetMultiplot(int rows, int cols) {
            Cmd("set multiplot");
            this.MultiplotRows = rows;
            this.MultiplotCols = cols;

            double xSubSize = (1.0 / this.MultiplotCols);
            double ySubSize = (1.0 / this.MultiplotRows);
            this.Cmd("set size {0},{1}", xSubSize.ToStringDot(), ySubSize.ToStringDot());
        }

        /// <summary>
        /// Changes the current sub-plot and optionally assigns a title to it
        /// </summary>
        /// <param name="iRow"></param>
        /// <param name="jCol"></param>
        /// <param name="title"></param>
        public void SetSubPlot(int iRow, int jCol, string title = null) {
            if (iRow < 0 || iRow >= this.MultiplotRows)
                throw new IndexOutOfRangeException();
            if (jCol < 0 || jCol >= this.MultiplotCols)
                throw new IndexOutOfRangeException();

            double xSubSize = (1.0 / this.MultiplotCols);
            double ySubSize = (1.0 / this.MultiplotRows);
            double XOrigin = (jCol) * xSubSize;
            double YOrigin = (this.MultiplotRows - iRow - 1) * ySubSize;

            this.Cmd("set size {0},{1}", xSubSize.ToStringDot(), ySubSize.ToStringDot());
            this.Cmd("set origin {0},{1}", XOrigin.ToStringDot(), YOrigin.ToStringDot());

            if (title != null) {
                this.Cmd("set title \"{0}\"", title);
            }
        }

        /// <summary>
        /// Guess what
        /// </summary>
        public void SetXAutorange() {
            Cmd("set autoscale x");
        }

        /// <summary>
        /// Guess what
        /// </summary>
        public void SetX2Autorange() {
            Cmd("set autoscale x2");
        }

        /// <summary>
        /// Guess what
        /// </summary>
        /// <param name="label"></param>
        public void SetXLabel(string label) {
            if (label != null) {
                Cmd("set xlabel \"" + label + "\"");
            } else {
                Cmd("unset xlabel");
            }
        }

        /// <summary>
        /// Guess what
        /// </summary>
        /// <param name="label"></param>
        public void SetX2Label(string label) {
            if (label != null) {
                Cmd("set x2label \"" + label + "\"");
            } else {
                Cmd("unset x2label");
            }
        }

        /// <summary>
        /// Defines the desired range of the x-axis
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        public void SetXRange(double min, double max) {
            using (StringWriter stw = new StringWriter()) {
                stw.Write("set xrange [");
                stw.Write(min.ToString(nfoi));
                stw.Write(":");
                stw.Write(max.ToString(nfoi));
                stw.Write("]");
                Cmd(stw.ToString());
            }
        }

        /// <summary>
        /// Defines the desired range of the x-axis
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        public void SetX2Range(double min, double max) {
            using (StringWriter stw = new StringWriter()) {
                stw.Write("set x2range [");
                stw.Write(min.ToString(nfoi));
                stw.Write(":");
                stw.Write(max.ToString(nfoi));
                stw.Write("]");
                Cmd(stw.ToString());
            }
        }

        /// <summary>
        /// Guess what
        /// </summary>
        public void SetYAutorange() {
            Cmd("set autoscale y");
        }

        /// <summary>
        /// Guess what
        /// </summary>
        public void SetY2Autorange() {
            Cmd("set autoscale y2");
        }

        /// <summary>
        /// Guess what
        /// </summary>
        /// <param name="label"></param>
        /// <param name="xOffset"></param>
        /// <param name="yOffset"></param>
        public void SetYLabel(string label, double xOffset = 0, double yOffset = 0) {
            if (label != null) {
                string offset;
                if (xOffset != 0.0 || yOffset != 0.0) {
                    offset = " offset " + xOffset.ToStringDot() + "," + yOffset.ToStringDot();
                } else {
                    offset = "";
                }

                Cmd("set ylabel \"" + label + "\"" + offset);
            } else {
                Cmd("unset ylabel");
            }
        }

        /// <summary>
        /// Guess what
        /// </summary>
        /// <param name="label"></param>
        /// <param name="xOffset"></param>
        /// <param name="yOffset"></param>
        public void SetY2Label(string label, double xOffset = 0, double yOffset = 0) {
            if (label != null) {
                string offset;
                if (xOffset != 0.0 || yOffset != 0.0) {
                    offset = " offset " + xOffset.ToStringDot() + "," + yOffset.ToStringDot();
                } else {
                    offset = "";
                }

                Cmd("set y2label \"" + label + "\"" + offset);
            } else {
                Cmd("unset y2label");
            }
        }

        /// <summary>
        /// Set minimum and maximum for y-axis
        /// </summary>
        public void SetYRange(double min, double max) {
            using (StringWriter stringWriter = new StringWriter()) {
                stringWriter.Write("set yrange [");
                stringWriter.Write(min.ToString(nfoi));
                stringWriter.Write(":");
                stringWriter.Write(max.ToString(nfoi));
                stringWriter.Write("]");
                Cmd(stringWriter.ToString());
            }
        }

        /// <summary>
        /// Set minimum and maximum for y2-axis (right side y-axis).
        /// </summary>
        public void SetY2Range(double min, double max) {
            using (StringWriter stringWriter = new StringWriter()) {
                stringWriter.Write("set y2range [");
                stringWriter.Write(min.ToString(nfoi));
                stringWriter.Write(":");
                stringWriter.Write(max.ToString(nfoi));
                stringWriter.Write("]");
                Cmd(stringWriter.ToString());
            }
        }

        /// <summary>
        /// Gnuplot command `set title`;
        /// </summary>
        /// <param name="title">
        /// `null` resets the title, i.e. means `unset`
        /// </param>
        public void SetTitle(string title) {
            if (title != null) {
                Cmd("set title \"{0}\"", title);
            } else {
                Cmd("unset title ");
            }
        }

        /// <summary>
        /// Deactivates <see cref="SetMultiplot(int, int)"/>
        /// </summary>
        public void UnsetMultiplot() {
            Cmd("unset multiplot");
        }

        /// <summary>
        /// Inserts af deferred plot commands at the current position of the gnuplot command stream;
        /// Needs to be called by the user only when e.g. a multiplot window is changed.
        /// </summary>
        public void WriteDeferredPlotCommands() {
            if (this.m_DeferredPlotCommands.Count > 0) {
                using (var stw = new StringWriter()) {
                    stw.Write("plot ");

                    for (int i = 0; i < this.m_DeferredPlotCommands.Count; i++) {
                        stw.Write(this.m_DeferredPlotCommands[i]);
                        if (i < this.m_DeferredPlotCommands.Count - 1)
                            stw.Write(", ");
                    }

                    stw.Flush();
                    this.Cmd(stw.ToString());

                    this.m_DeferredPlotCommands.Clear();
                }
            }
        }

        /// <summary>
        /// Writes the complete gnuplot script into a string.
        /// </summary>
        public string GetAllCommandsString() {
            this.m_Commands.Flush();
            string BeforeDeferred = this.m_Commands.ToString(); // backup commands

            var hack = m_DeferredPlotCommands.ToList();
            this.WriteDeferredPlotCommands();
            this.m_Commands.Flush();
            m_DeferredPlotCommands = hack;

            string AllCommands = this.m_Commands.ToString();
            this.m_Commands.Dispose();

            this.m_Commands = new StringWriter();

            this.m_Commands.Write(BeforeDeferred);

            using (var w = new StringWriter()) {

                if (this.OutputFile != null) {
                    w.WriteLine("set output '{0}'", this.OutputFile);
                }

                if (this.Terminal != null) {
                    w.WriteLine("set terminal {0}", this.Terminal);
                }

                w.WriteLine(AllCommands);

                return w.ToString();
            }
        }

        /// <summary>
        /// Starts the gnuplot process, executes the commands and exits.
        /// </summary>
        public int RunAndExit() {
            if (this.m_Commands == null)
                throw new NotSupportedException("Illegal call at this point.");

            string AllCommands = this.GetAllCommandsString();

            int ExitCode = 0;

            if (this.gnucmd == null) {
                ProcessStartInfo psi = new ProcessStartInfo() {
                    FileName = m_sGNUPlotFileName,
                    RedirectStandardOutput = true,
                    RedirectStandardInput = true,
                    RedirectStandardError = true,
                    UseShellExecute = false,
                    CreateNoWindow = true
                };
                gnucmd = Process.Start(psi);

                if (this.gnucmd == null) {
                    throw new ApplicationException("Couldn't open connection to gnuplot.");
                }

                this.gnucmd.StandardInput.WriteLine(AllCommands);
                this.gnucmd.StandardInput.Flush();


                this.gnucmd.StandardInput.WriteLine("exit");
                this.gnucmd.StandardInput.Flush();

                this.gnucmd.WaitForExit();

                string stderrString = this.gnucmd.StandardError.ReadToEnd();
                string stdoutString = this.gnucmd.StandardOutput.ReadToEnd();

                ExitCode = this.gnucmd.ExitCode;

                if (stderrString.Length > 0) {
                    Console.Error.WriteLine("Gnuplot Error: " + stderrString);
                }

                this.gnucmd = null;
            } else {
                throw new NotSupportedException("not possible if gnuplot is already running.");
            }

            return ExitCode;
        }

        /// <summary>
        /// Starts the gnuplot process and executes the commands. 
        /// After calling this function, no further gnuplot command can be called.
        /// The only further valid call is <see cref="Dispose"/> --or-- <see cref="Reset"/>.
        /// </summary>
        public void Execute() {
            if (this.m_Commands == null)
                throw new NotSupportedException("Illegal call at this point.");

            string AllCommands = this.GetAllCommandsString();
            this.m_Commands.Dispose();
            this.m_Commands = null;
            this.m_DeferredPlotCommands = null;

            if (this.gnucmd == null) {
                ProcessStartInfo psi = new ProcessStartInfo() {
                    FileName = m_sGNUPlotFileName,
                    RedirectStandardOutput = true,
                    RedirectStandardInput = true,
                    RedirectStandardError = true,
                    UseShellExecute = false,
                    CreateNoWindow = true
                };
                gnucmd = Process.Start(psi);

                if (this.gnucmd == null) {
                    throw new ApplicationException("Couldn't open connection to gnuplot.");
                }

                this.gnucmd.StandardInput.WriteLine(AllCommands);
                this.gnucmd.StandardInput.Flush();
            }
        }

        /// <summary>
        /// Writes the gnuplot script defined so far, and all auxiliary data
        /// files, to file <paramref name="filename"/>.
        /// </summary>
        public void WriteScriptAndData(string filename) {
            if (this.m_Commands == null)
                throw new NotSupportedException("Illegal call at this point.");

            using (var gpScript = new StreamWriter(filename)) {

                string AllCommands = this.GetAllCommandsString();

                string OutpDir = Path.GetDirectoryName(filename);
                string BaseName = Path.GetFileNameWithoutExtension(filename);

                int i = 0;
                foreach (var tmpFile in this.m_TempFiles) {
                    i++;
                    Debug.Assert(File.Exists(tmpFile));

                    string scriptName = tmpFile.Replace(Path.DirectorySeparatorChar, '/');
                    string newName = BaseName + "_data_" + i + ".csv";

                    AllCommands = AllCommands.Replace(scriptName, newName);

                    File.Copy(tmpFile, Path.Combine(OutpDir, newName), true);

                }

                AllCommands = AllCommands + System.Environment.NewLine + "exit" + System.Environment.NewLine;

                gpScript.Write(AllCommands);
                gpScript.Flush();
            }
        }

        /// <summary>
        /// Re-activates the plotting after <see cref="Execute"/> has been called.
        /// </summary>
        public void Reset() {
            if (this.m_Commands != null)
                m_Commands.Dispose();
            this.m_Commands = new StringWriter();
            this.m_DeferredPlotCommands = new List<string>();
            this.DeleteTemporaryFiles();
        }

        /// <summary>
        /// 
        /// </summary>
        public void Dispose() {
            if (gnucmd == null) {

            } else {
                if (this.gnucmd.StandardInput.BaseStream != null
                    && this.gnucmd.StandardInput.BaseStream.CanWrite) {
                    try {
                        this.gnucmd.StandardInput.WriteLine("exit");
                        this.gnucmd.StandardInput.Flush();
                    } catch (Exception) {
                        // Swallow
                    }
                }
                if (!gnucmd.WaitForExit(5000)) {
                    gnucmd.Kill();
                }
                gnucmd = null;
            }

            DeleteTemporaryFiles();
        }

        /// <summary>
        /// Destructor.
        /// </summary>
        ~Gnuplot() {
            Dispose();
        }

        /// <summary>
        /// clear output window
        /// </summary>
        public void Clear() {
            this.Cmd("clear");
        }
    }
}
