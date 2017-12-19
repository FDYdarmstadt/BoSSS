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

using BoSSS.Solution.Gnuplot;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Form for live plotting of residual curves
    /// </summary>
    public partial class ResidualForm : Form {

        ResidualLog log;

        /// <summary>
        /// Plots a logfile as a static graph.
        /// </summary>
        /// <param name="log">Logfile to be plotted.</param>
        public ResidualForm(ResidualLog log) {
            this.log = log;
            InitializeComponent();
            drawGnuplotGraph();
        }
        /// <summary>
        /// Draws a graph with Gnuplot and saves it as a gif. Subsequently sets Picture of Form to the plot. 
        /// </summary>
        private void drawGnuplotGraph() {
            using (Gnuplot gp = new Gnuplot()) {
                gp.SetXLabel("Iteration");
                gp.SetYLabel(log.Norm + " Norm");
                gp.Cmd("set title \"Residual Plot\"");
                gp.Cmd("set key outside under horizontal box");
                gp.Cmd("set logscale y");
                gp.Cmd("set format y \"10^{%L}\"");
                gp.Cmd("set grid xtics ytics");

                List<string> KeysToPlot = log.Values.Keys.ToList();
                KeysToPlot.Remove("#Line");
                KeysToPlot.Remove("#Time");
                KeysToPlot.Remove("#Iter");

                for (int i = 0; i < KeysToPlot.Count; i++) {
                    gp.PlotXY(log.Values["#Time"], log.Values[KeysToPlot[i]], KeysToPlot[i],
                        new PlotFormat(lineColor: ((LineColors)i), Style: Styles.Lines));
                }

                Image graph = gp.PlotGIF();
                pictureBox1.Image = graph;
                pictureBox1.Refresh();
            }
        }
    }
}
