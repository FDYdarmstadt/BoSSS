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
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using BoSSS.Solution.Gnuplot;
using BoSSS.Foundation.IO;


namespace BoSSS.Application.BoSSSpad {
    public partial class ResidualFormLive : Form {
        //The Residuallog to be plotted
        ResidualLog log;
        //The time after which the data is replotted
        int millisecondsPolling;
        //Local or Global plot 
        bool local = false;
        string titleOfPlot = "\"Residual Plot - Global\"";

        /// <summary>
        /// Plots the Residuals live.
        /// </summary>
        /// <param name="log">Residual Log to be Plotted</param>
        /// <param name="millisecondsPolling">Refreshing rate</param>
        public ResidualFormLive(ResidualLog log, int millisecondsPolling) {
            this.log = log;
            this.millisecondsPolling = millisecondsPolling;

            //Initialize Form and plot the graph.
            InitializeComponent();
            refreshPlot();
        }

        private void m_Stop_Click(object sender, EventArgs e) {
            timer1.Enabled = !timer1.Enabled;
        }

        private void timer1_Tick(object sender, EventArgs e) {
            refreshPlot();
        } 
        //Change the view from local to global and viceversa.
        private void m_toggleViewspace_Click(object sender, EventArgs e) {
            local = !local;
            if (local) {
                titleOfPlot = "\"Residual Plot - Local\"";
            } else {
                titleOfPlot = titleOfPlot = "\"Residual Plot - Global\"";
            }
            refreshPlot();
        }

        /// <summary>
        /// Updates the logfile, which contains the Residualdata and plots the graph Image.
        /// </summary>
        private void refreshPlot() {
            if (local) {
                log.ReadResiduals();
            } else {
                log.UpdateResiduals();
            }

            Image graph = drawGnuplotGraph();
            pictureBox1.Image = graph;
            pictureBox1.Refresh();
        }

        /// <summary>
        /// Plots a line plot from a given residual log.
        /// The y-axis will be scaled logarithmically.
        /// </summary>
        private Image drawGnuplotGraph() {
            Image graph;

            using (Gnuplot gp = new Gnuplot()) {
                gp.SetXLabel("Iteration");
                gp.SetYLabel(log.Norm + " Norm");
                //gp.Cmd("set terminal wxt noraise title 'SessionGuid: " + log.SessionGuid.ToString() + "'");
                gp.Cmd("set title " + titleOfPlot);
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

                graph = gp.PlotGIF();
            }
            
            return(graph);
        }
    }
}
