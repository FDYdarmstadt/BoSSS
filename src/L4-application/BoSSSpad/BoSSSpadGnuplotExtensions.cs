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

using ilPSP;
using BoSSS.Platform;
using BoSSS.Solution.Gnuplot;
using System.Drawing;
using System.IO;
using System;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Default output of <see cref="BoSSSpadGnuplotExtensions.PlotNow(Gnuplot)"/> (and variants)
    /// </summary>
    public enum PlotNowMode {
        
        /// <summary>
        /// Deprecated
        /// </summary>
        GIF = 0,

        /// <summary>
        /// <see cref="GnuplotExtensions.PlotCairolatex(Gnuplot, string, double, double)"/>,
        /// used for LaTeX output
        /// </summary>
        CairoLatex = 1,

        /// <summary>
        /// <see cref="GnuplotExtensions.PlotSVG"/>,
        /// main option for Jupyter notebooks
        /// </summary>
        SVG = 2
    }


    /// <summary>
    /// Spice up Gnuplot.
    /// </summary>
    public static class BoSSSpadGnuplotExtensions {




        /// <summary>
        /// Hack for <see cref="PlotNow(Gnuplot)"/>
        /// </summary>
        public static PlotNowMode PlotMode = PlotNowMode.GIF;

        /// <summary>
        /// Gnuplot plotting, automatic choice of gnuplot driver depending on
        /// the current value of <see cref="PlotMode"/>.
        /// </summary>
        public static object PlotNow(this Gnuplot gp) {
            
            if (PlotMode == PlotNowMode.CairoLatex) {
                return gp.PlotCairolatex();
            } else if(PlotMode == PlotNowMode.GIF) {
                throw new NotImplementedException("Is this still used by anyone?");
                //return gp.PlotGIF();
            } else if(PlotMode == PlotNowMode.SVG) {
                return gp.PlotSVG();
            } else {
                throw new NotImplementedException("Unknown option: " + PlotMode);
            }
        }

        /// <summary>
        /// Gnuplot plotting (single plot), automatic choice of gnuplot driver depending on
        /// the current value of <see cref="PlotMode"/>.
        /// </summary>
        public static object PlotNow(this Plot2Ddata _2DData) {
            using (Gnuplot gp = _2DData.ToGnuplot()) {

                return PlotNow(gp);
            }
        }

        /// <summary>
        /// Gnuplot plotting (multi-plot), automatic choice of Gnuplot driver depending on
        /// the current value of <see cref="PlotMode"/>.
        /// </summary>
        public static object PlotNow(this Plot2Ddata[,] _2DData) {
            using (Gnuplot gp = _2DData.ToGnuplot()) {

                return PlotNow(gp);
            }
        }

        /// <summary>
        /// Gnuplot plotting (multiplot), automatic choice of gnuplot driver depending on
        /// the current value of <see cref="PlotMode"/>.
        /// </summary>
        public static object PlotNow(this CairolatexContainer cl) {

            if (PlotMode == PlotNowMode.CairoLatex) {
                return cl;
            } else {
                return cl.PreviewSVG();
            }
        }

        /*
        /// <summary>
        /// Single plot window:
        /// Converts <see cref="Plot2Ddata"/> into an alive Gnuplot object.
        /// </summary>
        public static Gnuplot ToGnuplot(this Plot2Ddata _2DData, GnuplotPageLayout layout = null) {
            if (layout != null)
                throw new NotImplementedException("todo");

            Gnuplot gp = new Gnuplot();

            _2DData.ToGnuplot(gp);
            return gp;
        }

        /// <summary>
        /// Multiple plot windows:
        /// Converts <see cref="Plot2Ddata"/> into an alive Gnuplot object.
        /// </summary>
        public static Gnuplot ToGnuplot(this Plot2Ddata[,] _2DData, GnuplotPageLayout layout = null) {
            if (layout != null)
                throw new NotImplementedException("todo");
            if (_2DData.GetLowerBound(0) != 0)
                throw new ArgumentException();
            if (_2DData.GetLowerBound(1) != 0)
                throw new ArgumentException();
            if (_2DData.GetLength(0) <= 0)
                throw new ArgumentException();
            if (_2DData.GetLength(1) <= 0)
                throw new ArgumentException();

            Gnuplot gp = new Gnuplot();

            gp.SetMultiplot(_2DData.GetLength(0), _2DData.GetLength(1));

            for (int iRow = 0; iRow < _2DData.GetLength(0); iRow++) {
                for (int iCol = 0; iCol < _2DData.GetLength(1); iCol++) {
                    if (_2DData[iRow, iCol] != null) {
                        gp.SetSubPlot(iRow, iCol);
                        _2DData[iRow, iCol].ToGnuplot(gp);
                    }
                }
            }
            return gp;
        }
        */

    }
}
