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
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using BoSSS.Platform;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using System.Runtime.Serialization;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// A type representing multiple sets of abscissas with the corresponding sets of values.
    /// </summary>
    [Serializable]
    [DataContract]
    public class Plot2Ddata {

        /// <summary>
        /// Represents a single set of abscissas with the corresponding set of values.
        /// </summary>
        [Serializable]
        [DataContract]
        public class XYvalues : ICloneable {

            /// <summary>
            /// Gnuplot Style to use.
            /// </summary>
            public BoSSS.Solution.Gnuplot.PlotFormat Format = new PlotFormat() {
                PointType = PointTypes.Circle,
                DashType = DashTypes.Solid,
                LineColor = LineColors.Black,
                LineWidth = 1,
                PointSize = 3,
                Style = Styles.LinesPoints
            };


            /// <summary>
            /// The name of the group, i.e. the name which may show up in the legend.
            /// </summary>
            [DataMember]
            public string Name;

            /// <summary>
            /// The points of evaluation, i.e. x-values
            /// </summary>
            [DataMember]
            public double[] Abscissas;

            /// <summary>
            /// The values at the <see cref="Abscissas"/>, i.e. y-values
            /// </summary>
            [DataMember]
            public double[] Values;

            /// <summary>
            /// Constructs a data group.
            /// </summary>
            /// <param name="name">
            /// The name of the group
            /// </param>
            /// <param name="abscissas">
            /// The points of evaluation, i.e. x-values
            /// </param>
            /// <param name="values">
            /// The values at the <see cref="Abscissas"/>, i.e. y-values. Note
            /// that the length must be equal to the length of
            /// <paramref name="abscissas"/>.
            /// </param>
            public XYvalues(string name, double[] abscissas, double[] values) {
                if (abscissas.Length != values.Length) {
                    throw new ArgumentException(
                        "Number of x and y values must be identical within each group");
                }

                this.Name = name;
                this.Abscissas = abscissas;
                this.Values = values;
            }

            /// <summary>
            /// Constructs a data group.
            /// </summary>
            /// <param name="name">
            /// The name of the group
            /// </param>
            public XYvalues(string name) {

                this.Name = name;
                this.Abscissas = new double[0];
                this.Values = new double[0];
            }

            #region ICloneable Members

            /// <summary>
            /// Clone
            /// </summary>
            public object Clone() {
                return new XYvalues(
                    this.Name.CloneAs(),
                    this.Abscissas.CloneAs(),
                    this.Values.CloneAs()) {
                    Format = this.Format.CloneAs()
                };
            }

            #endregion
        }

        /// <summary>
        /// For each group: A set of data consisting of the x- and the
        /// y-values.
        /// </summary>
        /// <remarks>
        /// The logarithm is not applied to this data, even if
        /// <see cref="LogX"/> or <see cref="LogY"/> is true.
        /// This will be done on the fly during
        /// the post-processing of the data
        /// </remarks>
        [DataMember]
        public XYvalues[] dataGroups;

        /// <summary>
        /// Indicates whether the abscissas should be scaled logarithmically.
        /// </summary>
        [DataMember]
        public bool LogX;

        /// <summary>
        /// Indicates whether the values should be scaled logarithmically.
        /// </summary>
        [DataMember]
        public bool LogY;


        /// <summary>
        /// y-range minimum, optional 
        /// </summary>
        [DataMember]
        public double? XrangeMin = null;

        /// <summary>
        /// x-range maximum, optional 
        /// </summary>
        [DataMember]
        public double? XrangeMax = null;

        /// <summary>
        /// y-range minimum, optional 
        /// </summary>
        [DataMember]
        public double? YrangeMin = null;

        /// <summary>
        /// y-range maximum, optional 
        /// </summary>
        [DataMember]
        public double? YrangeMax = null;

        /// <summary>
        /// Label for X-axis
        /// </summary>
        [DataMember]
        public string Xlabel = null;

        /// <summary>
        /// Label for secondary X-Axis
        /// </summary>
        [DataMember]
        public string X2label = null;

        /// <summary>
        /// Label for Y-Axis
        /// </summary>
        [DataMember]
        public string Ylabel = null;

        /// <summary>
        /// Label for secondary Y-Axis
        /// </summary>
        [DataMember]
        public string Y2label = null;
        

        /// <summary>
        /// Constructs a new, empty plot.
        /// </summary>
        public Plot2Ddata() {
            this.dataGroups = new XYvalues[0];
        }

        /// <summary>
        /// Constructs a new, lightweight <see cref="Plot2Ddata"/> for the given
        /// data.
        /// </summary>
        /// <param name="dataRows">
        /// For each (group) key: A data set of where the first index
        /// corresponds to the abscissas (a.k.a. the x-values) and the second
        /// index corresponds to the values at these abscissas (a.k.a. the
        /// y-values).
        /// </param>
        public Plot2Ddata(params KeyValuePair<string, double[][]>[] dataRows)
            : this() {
            this.dataGroups = dataRows.
                Select(p => new XYvalues(p.Key, p.Value[0], p.Value[1])).
                OrderBy(p => p.Name).
                ToArray();
        }

        /// <summary>
        /// Constructs a new, lightweight <see cref="Plot2Ddata"/> for a single set
        /// of values.
        /// </summary>
        /// <param name="xyPairs">
        /// The abscissas a.k.a. x-coordinates paired with the corresponding
        /// y-values
        /// the y-values.
        /// </param>
        /// <param name="groupKey">
        /// An optional name of the given row.
        /// </param>
        public Plot2Ddata(IEnumerable<KeyValuePair<double, double>> xyPairs, string groupKey = "unnamed")
            : this() {
            this.dataGroups = new XYvalues[] {
                new XYvalues(
                    groupKey,
                    xyPairs.Select(p => p.Key).ToArray(),
                    xyPairs.Select(p => p.Value).ToArray())
            };
        }

        /// <summary>
        /// Constructs a new, lightweight <see cref="Plot2Ddata"/> for a single set
        /// of values.
        /// </summary>
        /// <param name="abscissas">
        /// The abscissas a.k.a. x-coordinates.
        /// </param>
        /// <param name="values">
        /// A set of values at the given <paramref name="abscissas"/>; a.k.a.
        /// the y-values.
        /// </param>
        /// <param name="groupKey">
        /// An optional name of the given row.
        /// </param>
        public Plot2Ddata(IEnumerable<double> abscissas, IEnumerable<double> values, string groupKey = "unnamed")
            : this() {
            this.dataGroups = new XYvalues[] {
                new XYvalues(groupKey, abscissas.ToArray(), values.ToArray())
            };
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        /// <param name="originalSet">
        /// Object to be copied from.
        /// </param>
        private Plot2Ddata(Plot2Ddata originalSet)
            : this() {
            this.dataGroups = originalSet.dataGroups;
            this.LogX = originalSet.LogX;
            this.LogY = originalSet.LogY;
        }

        /// <summary>
        /// Another copy constructor used by <see cref="Merge"/> and <see cref="Extract"/>
        /// </summary>
        /// <param name="groups">
        /// The data groups (see <see cref="dataGroups"/>) to be included in
        /// the new object
        /// </param>
        private Plot2Ddata(params XYvalues[] groups) : this() {
            this.dataGroups = groups.OrderBy(p => p.Name).ToArray();
        }

        /// <summary>
        /// Creates a new data set where <see cref="LogX"/> is to true. Useful
        /// for command chaining on the console.
        /// </summary>
        /// <returns>
        /// A copy of this object where <see cref="LogX"/> equals true.
        /// </returns>
        public Plot2Ddata WithLogX() {
            var set = new Plot2Ddata(this);
            set.LogX = true;
            return set;
        }

        /// <summary>
        /// Creates a new data set where <see cref="LogY"/> is to true. Useful
        /// for command chaining on the console.
        /// </summary>
        /// <returns>
        /// A copy of this object where <see cref="LogY"/> equals true.
        /// </returns>
        public Plot2Ddata WithLogY() {
            var set = new Plot2Ddata(this);
            set.LogY = true;
            return set;
        }

        /// <summary>
        /// Merges the this object and the given data set
        /// <paramref name="other"/> into a new data set (while making a deep
        /// copy of the data)
        /// </summary>
        /// <param name="other">
        /// The data set to be merged into this one
        /// </param>
        /// <returns>
        /// A data set containing all data in this object and
        /// <paramref name="other"/>
        /// </returns>
        public Plot2Ddata Merge(Plot2Ddata other) {
            if (this.LogX != other.LogX || this.LogY != other.LogY) {
                throw new Exception("Data sets have incompatible logarithmic scaling options");
            }

            IList<XYvalues> mergedGroups = new List<XYvalues>(this.dataGroups.Length + other.dataGroups.Length);
            mergedGroups.AddRange(this.dataGroups.Select(g => g.CloneAs()));
            foreach (XYvalues otherGroup in other.dataGroups) {
                if (this.dataGroups.Any(g => g.Name == otherGroup.Name)) {
                    throw new NotSupportedException(String.Format(
                        "Group key '{0}' exists in both data sets. This is not supported.",
                        otherGroup.Name));
                }

                mergedGroups.Add(otherGroup.CloneAs());
            }

            Plot2Ddata result = new Plot2Ddata(mergedGroups.ToArray());
            result.LogX = this.LogX;
            result.LogY = this.LogY;
            return result;
        }

        /// <summary>
        /// Extracts the data rows whose the given
        /// <paramref name="groupNames"/> and creates new data set (while
        /// making a deep copy of the data)
        /// </summary>
        /// <param name="groupNames">
        /// The names of the groups to be selected
        /// </param>
        /// <returns>
        /// A data set containing all data corresponding to the groups with
        /// names in  <paramref name="groupNames"/>
        /// </returns>
        public Plot2Ddata Extract(params string[] groupNames) {
            return new Plot2Ddata(dataGroups.
                Where(g => groupNames.Contains(g.Name)).
                Select(g => g.CloneAs()).
                ToArray());
        }

        /// <summary>
        /// Subtracts data rows with the given <paramref name="groupNames"/>
        /// and creates new data set (while making a deep copy of the data)
        /// </summary>
        /// <param name="groupNames">
        /// The names of the groups to be excluded
        /// </param>
        /// <returns>
        /// A data set containing all data within this object, except of groups
        /// whose names are listed in <paramref name="groupNames"/>
        /// </returns>
        public Plot2Ddata Without(params string[] groupNames) {
            return new Plot2Ddata(dataGroups.
                Where(g => !groupNames.Contains(g.Name)).
                Select(g => g.CloneAs()).
                ToArray());
        }

        /// <summary>
        /// For each data row: Computes the regression coefficients for a
        /// linear fit to the data (or their respective logarithms depending on
        /// <see cref="LogX"/> and <see cref="LogY"/>).
        /// </summary>
        /// <returns>
        /// For each data row: A tuple where
        /// <see cref="Tuple{Double,Double}.Item1"/> and
        /// <see cref="Tuple{Double,Double}.Item2"/> denote the slope and the
        /// affine offset of the linear fit, respectively.
        /// </returns>
        public IEnumerable<KeyValuePair<string, double>> Regression() {
            foreach (var group in dataGroups) {
                double[] xValues;
                if (LogX) {
                    xValues = group.Abscissas.Select(x => x.Log10()).ToArray();
                } else {
                    xValues = group.Abscissas;
                }
                double xAvg = xValues.Average();

                double[] yValues;
                if (LogY) {
                    yValues = group.Values.Select(y => y.Log10()).ToArray();
                } else {
                    yValues = group.Values;
                }
                double yAvg = yValues.Average();

                double v1 = 0.0;
                double v2 = 0.0;

                for (int i = 0; i < yValues.Length; i++) {
                    v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
                    v2 += Math.Pow(xValues[i] - xAvg, 2);
                }

                yield return new KeyValuePair<string, double>(group.Name, v1 / v2);
            }
        }

        /// <summary>
        /// Visualizes the stored data using gnuplot.
        /// </summary>
        public void Plot() {
            Gnuplot gp = new Gnuplot();
            gp.SetXLabel("x");
            gp.SetYLabel("y");
            gp.Cmd("set terminal wxt noraise");
            if (LogX) {
                gp.Cmd("set logscale x");
                gp.Cmd("set format x \"10^{%L}\"");
            }
            if (LogY) {
                gp.Cmd("set logscale y");
                gp.Cmd("set format y \"10^{%L}\"");
            }
            gp.Cmd("set grid xtics ytics");

            int lineColor = 0;
            foreach (var group in dataGroups) {
                gp.PlotXY(group.Abscissas, group.Values, group.Name,  
                    new PlotFormat(lineColor:((LineColors)( ++lineColor)), pointType: ((PointTypes)4), pointSize:1.5, Style: Styles.LinesPoints));
            }
            gp.Execute();
        }

        /// <summary>
        /// Saves a tabular summary of the stored data to a text file
        /// </summary>
        /// <param name="path">
        /// Path to the file
        /// </param>
        public void SaveToTextFile(string path) {
            Console.WriteLine("Warning: Handling of logarithms still unclear in 'SaveToTextFile', Raw Data Written");

            // open file
            using (StreamWriter stw = new StreamWriter(path)) {
                stw.Write(this.ToString());
                stw.Close();
            }
        }

        /// <summary>
        /// Saves two text files: 
        /// 1) A tabular summary of the stored data, see <see cref="SaveTabular(string)"/>
        /// 2) A table of the linear regression values, see <see cref="Regression()"/>
        /// File name convention:
        /// 1) <paramref name="path"/>+Data.txt
        /// 1) <paramref name="path"/>+Rgrs.txt 
        /// </summary>
        /// <param name="path">
        /// Path to file
        /// </param>
        public void SaveTextFileToPublish(string path){
            // writing data
            string pathWithoutExt = System.IO.Path.ChangeExtension(path, null);
            string newPath = pathWithoutExt+"Data.txt";
            SaveTabular(newPath);
            // writing regression
            newPath = pathWithoutExt + "Rgrs.txt";
            var regressionData = this.Regression();
            using (StreamWriter stw = new StreamWriter(newPath)) {
                stw.WriteLine("\\$\\degree$ \t EOC");
                foreach (var item in regressionData) {
                    stw.WriteLine(item.Key + "\t" + item.Value);
                }
                stw.Close();
            }
        }

        /// <summary>
        /// Saves a gnuplot file that can be used to plot the data represented
        /// by this data set
        /// </summary>
        /// <param name="path">
        /// Path to the gnuplot file
        /// </param>
        public void SaveGnuplotFile(string path) {
            using (StreamWriter stw = new StreamWriter(path)) {
                if (LogX) {
                    stw.WriteLine("set logscale x");
                    stw.WriteLine("set format x \"10^{%L}\"");
                }
                if (LogY) {
                    stw.WriteLine("set logscale y");
                    stw.WriteLine("set format y \"10^{%L}\"");
                }
                stw.WriteLine("plot for [IDX=0:{0}] \"-\" using 1:2 with lines title columnheader(1)", dataGroups.Length - 1);
                stw.Write(this.ToString(true));
                stw.WriteLine("pause -1");
                stw.Close();
            }
        }

        /// <summary>
        /// Saves pgfplots code to a file which can in turn be used to plot
        /// the data represented by this data set
        /// </summary>
        /// <param name="path">
        /// Path to the gnuplot file
        /// </param>
        public void SavePgfplotsFile(string path) {
            using (StreamWriter s = new StreamWriter(path)) {
                string axisdefinition;
                if (LogX & !LogY) {
                    axisdefinition = "logxaxis";
                } else if (!LogX & LogY) {
                    axisdefinition = "logyaxis";
                } else if (LogX && LogY) {
                    axisdefinition = "loglogaxis";
                } else {
                    axisdefinition = "axis";
                }

                char RegressionCounter = 'a';

                s.Write(@"
%% Uncomment these lines to test the raw output.
%% The packages below are mandatory.
%\documentclass{scrartcl} 
%\usepackage{tikz} 
%\usepackage{pgfplots} 
%   \pgfplotsset{compat=1.3}
%\usepackage{pgfplotstable}

%\begin{document} 
%\begin{figure} 
%\centering ");
                s.WriteLine(@"\begin{tikzpicture}");
                s.WriteLine("\\begin{{{0}}}", axisdefinition);
                s.WriteLine(@"[legend style={at={(1.1,0.95)},anchor=north west,fill=none,draw=none,nodes=right}]");
                foreach (var group in dataGroups) {
                    //Write raw data points
                    s.WriteLine(@"\addplot[only marks] table{");
                    for (int i = 0; i < group.Abscissas.Length; i++) {
                        s.Write(group.Abscissas[i].ToString("E16", NumberFormatInfo.InvariantInfo));
                        s.Write("\t");
                        s.Write(group.Values[i].ToString("E16", NumberFormatInfo.InvariantInfo));
                        s.Write("\n");
                    }
                    s.WriteLine(@"};");
                    s.WriteLine("\\addlegendentry{{Group:{0}}};", group.Name);
                    //Write regression line
                    s.WriteLine(@"\addplot[mark=none, solid] table[x=X,y={create col/linear regression={y=Y}}]{");
                    s.WriteLine("X \t Y");
                    for (int i = 0; i < group.Abscissas.Length; i++) {
                        s.Write(group.Abscissas[i].ToString("E16", NumberFormatInfo.InvariantInfo));
                        s.Write("\t");
                        s.Write(group.Values[i].ToString("E16", NumberFormatInfo.InvariantInfo));
                        s.Write("\n");
                    }
                    s.WriteLine(@"};");
                    //add legend entry for slope of regression line
                    s.Write("\\makeatletter \n \\global\\let\\slope{0} = \\pgfplotstableregressiona \n	\\addlegendentry{{slope: $\\pgfmathprintnumber{{\\slope{0}}}$}}; \\makeatother", RegressionCounter);
                    s.Write("\n");
                    RegressionCounter++; //go to next letter
                }
                s.WriteLine("\\end{{{0}}}", axisdefinition);
                s.WriteLine(@"\end{tikzpicture}");
                s.Write(@"%\end{figure} 
%\end{document}
");
                s.Close();
            }
        }

        /// <summary>
        /// See <see cref="ToString(bool)"/>
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return this.ToString(false);
        }

        /// <summary>
        /// Assembles a tabular summary of the data.
        /// </summary>
        /// <param name="includeBlockEndStatement">
        /// If true "end" is written at the end of every block of grouped points.
        /// This is needed for Gnuplot Output.
        /// </param>
        /// <returns></returns>
        public string ToString(bool includeBlockEndStatement) {
            StringWriter s = new StringWriter(NumberFormatInfo.InvariantInfo);
            // Name of the Dataset
            foreach (var group in dataGroups) {
                s.WriteLine("Grouping:{0}", group.Name);
                s.WriteLine("# x \t y");
                for (int i = 0; i < group.Abscissas.Length; i++) {
                    s.Write(group.Abscissas[i].ToString("E16", NumberFormatInfo.InvariantInfo));
                    s.Write("\t");
                    s.Write(group.Values[i].ToString("E16", NumberFormatInfo.InvariantInfo));
                    s.Write("\n");
                }
                if (includeBlockEndStatement) {
                    s.WriteLine("end");
                };
                s.Write("\n \n");
            }
            s.Flush();
            s.Close();

            return s.ToString();
        }

        /// <summary>
        /// Saves the data set set in a purely tabular format with the
        /// following layout:
        /// <code>
        /// group   x   y
        /// groupName1  1   2
        /// groupName1  2   3
        /// groupName1  3   3
        /// groupName2  1   3
        /// groupName2  2   4
        /// groupName2  3   2
        /// </code>
        /// </summary>
        /// <param name="path">File path</param>
        public void SaveTabular(string path) {
            using (StreamWriter s = new StreamWriter(path)) {
                s.WriteLine("group\tx\ty");

                foreach (var group in dataGroups) {
                    for (int i = 0; i < group.Abscissas.Length; i++) {
                        s.Write(group.Name + "\t");
                        s.Write(group.Abscissas[i].ToString("E16", NumberFormatInfo.InvariantInfo) + "\t");
                        s.Write(group.Values[i].ToString("E16", NumberFormatInfo.InvariantInfo));
                        s.WriteLine();
                    }
                }
            }
        }
    }
}
