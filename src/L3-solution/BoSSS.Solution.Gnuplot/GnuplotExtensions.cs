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
using System.Collections.Generic;
//using Microsoft.AspNetCore.Html;

namespace BoSSS.Solution.Gnuplot {


    /// <summary>
    /// Spice up <see cref="Plot2Ddata"/>
    /// </summary>
    public static class GnuplotExtensions {
        
        
        /// <summary>
        /// Executes a pre-configured Gnuplot object and saves the 
        /// plot to a PNG file (`set terminal png`).
        /// </summary>
        /// <param name="gp">
        /// pre-configured Gnuplot; must already contain the data to plot.
        /// </param>
        /// <param name="xRes">Horizontal resolution in pixels.</param>
        /// <param name="yRes">Vertical resolution in pixels.</param>
        /// <param name="OutfileName">Path/filename for output PNG.</param>
        static public void SaveToGIF(this Gnuplot gp, string OutfileName, int xRes = 800, int yRes = 600) {

            if(xRes <= 0)
                throw new ArgumentOutOfRangeException();
            if(yRes <= 0)
                throw new ArgumentOutOfRangeException();

            // set terminal
            gp.Terminal = string.Format("pngcairo size {0},{1}", xRes, yRes);

            // set output file
            gp.OutputFile = OutfileName;

            // call gnuplot
            int exCode = gp.RunAndExit(); // run & close gnuplot
            if (exCode != 0) {
                throw new IOException("Gnuplot-internal error: exit code " + exCode);
            }

            // return image
            var fi = (new FileInfo(OutfileName));
            if (fi.Exists && fi.Length > 0) {
                
            } else {
                throw new IOException($"Gnuplot output file empty or non-existent.");
            }
        }
        
        /// <summary>
        /// Plot to a scalable vector graphics file (`set terminal svg`).
        /// </summary>
        /// <param name="gp"></param>
        /// <param name="xRes">Horizontal resolution in pixels.</param>
        /// <param name="yRes">Vertical resolution in pixels.</param>
        /// <param name="OutfileName">Path/filename for output SVG.</param>
        static public void SaveToSVG(this Gnuplot gp, string OutfileName, int xRes = 800, int yRes = 600) {
            if(xRes <= 0)
                throw new ArgumentOutOfRangeException();
            if(yRes <= 0)
                throw new ArgumentOutOfRangeException();


            // set terminal
            gp.Terminal = string.Format("svg enhanced background rgb 'white' size {0},{1}", xRes, yRes);

            // set output file
            gp.OutputFile = OutfileName;

            // call gnuplot
            int exCode = gp.RunAndExit(); // run & close gnuplot
            if (exCode != 0) {
                throw new IOException("Gnuplot-internal error: exit code " + exCode);
            }

            // return image
            var fi = (new FileInfo(OutfileName));
            if (fi.Exists && fi.Length > 0) {
            } else {
                throw new IOException($"Gnuplot output file empty or non-existent.");
            }
        }


        /*
        /// <summary>
        /// Plotting using Gnuplot with Cairolatex output.
        /// </summary>
        /// <param name="xSize">Horizontal size in centimeters, ignored if <paramref name="xSize"/> or <paramref name="ySize"/> is negative.</param>
        /// <param name="ySize">Vertical size in centimeters.</param>
        /// <param name="Options">
        /// Options for gnuplot cairolatex terminal
        /// </param>
        /// <param name="plot"></param>
        /// <returns>
        /// A memory-image of gnuplot cairolatex output.
        /// </returns>
        static public CairolatexContainer PlotCairolatex(this Plot2Ddata plot,
            //string Options = " pdf input noheader blacktext nobackground noenhanced fontscale 0.6 ", 
            string Options = " pdf  ",
            double xSize = 14, double ySize = 10.5) {


            using(var gp = plot.ToGnuplot()) {
                return gp.PlotCairolatex(Options, xSize, ySize);
            }

        }
        

        /// <summary>
        /// Plotting using Gnuplot with Cairolatex output.
        /// </summary>
        /// <param name="xSize">Horizontal size in centimeters, ignored if <paramref name="xSize"/> or <paramref name="ySize"/> is negative.</param>
        /// <param name="ySize">Vertical size in centimeters.</param>
        /// <param name="Options">
        /// Options for gnuplot cairolatex terminal
        /// </param>
        /// <param name="plots"></param>
        /// <returns>
        /// A memory-image of gnuplot cairolatex output.
        /// </returns>
        static public CairolatexContainer PlotCairolatex(this Plot2Ddata[,] plots,
            //string Options = " pdf input noheader blacktext nobackground noenhanced fontscale 0.6 ", 
            string Options = " pdf  ",
            double xSize = 14, double ySize = 10.5) {

           

            using(var gp = plots.ToGnuplot()) {
                return gp.PlotCairolatex(Options, xSize, ySize);
            }

        }
        
     
        /// <summary>
        /// Plotting using Gnuplot with Cairolatex output.
        /// </summary>
        /// <param name="xSize">Horizontal size in centimeters, ignored if <paramref name="xSize"/> or <paramref name="ySize"/> is negative.</param>
        /// <param name="ySize">Vertical size in centimeters.</param>
        /// <param name="Options">
        /// Options for gnuplot cairolatex terminal
        /// </param>
        /// <param name="gp"></param>
        /// <returns>
        /// A memory-image of gnuplot cairolatex output.
        /// </returns>
        static public CairolatexContainer PlotCairolatex(this Gnuplot gp,
            //string Options = " pdf input noheader blacktext nobackground noenhanced fontscale 0.6 ", 
            string Options = " pdf  ",
            double xSize = 14, double ySize = 10.5) {

            // return object
            var clc = new CairolatexContainer();

            // set terminal
            if (xSize >= 0 && ySize >= 0)
                gp.Terminal = string.Format("cairolatex {0} size {1}cm,{2}cm", Options != null ? Options : " ", xSize.ToStringDot(), ySize.ToStringDot());
            else
                gp.Terminal = string.Format("cairolatex {0} size {1}cm,{2}cm", Options != null ? Options : " ");

            // set output file
            string baseName = MyGetTempFileName();
            baseName = Path.Combine(Path.GetDirectoryName(baseName), Path.GetFileNameWithoutExtension(baseName));
            string TexOutfileName = baseName + ".tex";
            gp.OutputFile = TexOutfileName;

            string GraphisOut = Path.Combine(Path.GetDirectoryName(TexOutfileName), baseName);

            // gnuplot script
            {
                string[] tmpFiles = gp.TempFilesPath.ToArray();

                clc.DataFiles = new string[tmpFiles.Length];
                clc.DataFileNames = new string[tmpFiles.Length];

                string AllCommands = gp.GetAllCommandsString();
                
                for (int i = 0; i < tmpFiles.Length; i++) {
                    string tmpFile = tmpFiles[i];
                    
                    string scriptName = tmpFile.Replace(Path.DirectorySeparatorChar, '/');
                    string newName = baseName + "_data_" + i + ".csv";

                    AllCommands = AllCommands.Replace(scriptName, newName);

                    clc.DataFiles[i] = File.ReadAllText(tmpFile);
                    clc.DataFileNames[i] = newName;
                }

                AllCommands = AllCommands.Replace(TexOutfileName, Path.GetFileName(TexOutfileName));

                AllCommands = AllCommands + System.Environment.NewLine + "exit" + System.Environment.NewLine;

                clc.GnuplotScript = AllCommands;
            }

            // call gnuplot
            {
                int exCode = gp.RunAndExit(); // run & close gnuplot
                if (exCode != 0) {
                    Console.WriteLine("Gnuplot-internal error: exit code " + exCode);
                    return null;
                }
            }

            // return graphics
            {
                string TexContent = File.ReadAllText(TexOutfileName);

                string GraphisOut_ext = null;
                if (File.Exists(GraphisOut + ".eps"))
                    GraphisOut_ext = ".eps";
                else if (File.Exists(GraphisOut + ".pdf"))
                    GraphisOut_ext = ".pdf";
                else
                    throw new FileNotFoundException(string.Format("Unable to find either eps or pdf file: {0}.pdf or {0}.eps.", GraphisOut));

                byte[] PdfOrEps = File.ReadAllBytes(GraphisOut + GraphisOut_ext);
                TexContent = TexContent.Replace(GraphisOut, Path.GetFileNameWithoutExtension(baseName)); // hack replacement of absolute path.
                clc.LatexCode = TexContent;
                clc.GraphicsData = PdfOrEps;
                clc.GraphicsFilename = Path.GetFileName(GraphisOut) + GraphisOut_ext;
            }

            return clc;  
        }
        */

        /// <summary>
        /// Add a Log slope to the <see cref="Gnuplot"/>
        /// </summary>
        /// <param name="gp"></param>
        /// <param name="_2DData"></param>
        /// <param name="mode">'a' = average, '+' = max, '-' = min</param>
        /// <param name="round"></param>
        /// <param name="format"></param>
        public static void PlotLogSlope(this Gnuplot gp, Plot2Ddata _2DData, char mode = 'a', PlotFormat format = null, double round = 0.5) {

            var slopes = _2DData.Regression();
            double slope = 0.0;
            if(mode == '+') {
                slope = slopes.Max(kv => kv.Value);
            }else if (mode == '-') {
                slope = slopes.Min(kv => kv.Value);
            } else {
                slope = slopes.Average(kv => kv.Value);
            }

            double yMin_mag = double.MaxValue;
            double yMax_mag = double.MinValue;
            double xMin_mag = double.MaxValue;
            double xMax_mag = double.MinValue;

            foreach (var grp in _2DData.dataGroups) {
                yMin_mag = Math.Min(yMin_mag, Math.Log(grp.Values.Min(), 10));
                yMax_mag = Math.Max(yMax_mag, Math.Log(grp.Values.Max(), 10));
                xMin_mag = Math.Min(xMin_mag, Math.Log(grp.Abscissas.Min(), 10));
                xMax_mag = Math.Max(xMax_mag, Math.Log(grp.Abscissas.Max(), 10));
            }

            double[] point = new double[2];
            point[0] = Math.Pow(10, xMin_mag + (1.0 / 3.0) * (xMax_mag - xMin_mag));
            point[1] = slope > 0 ? Math.Pow(10, yMin_mag) : Math.Pow(10, yMax_mag);
            double size = Math.Pow(10, xMin_mag + (1.0 / 3.0 + 1.0 / 10.0) * (xMax_mag - xMin_mag)) - point[0];

            gp.PlotLogSlope(slope, point, null, format, true, false, false, slope != 0.0, size, 'a', round);
        }

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

        /// <summary>
        /// Multiple plot windows:
        /// Converts multiple <see cref="Plot2Ddata"/> into an alive Gnuplot object.
        /// </summary>
        /// <param name="_2DData">a list of plots</param>
        /// <param name="layout">
        /// describes the 2d-multiplot arrangement
        /// - 1st index: multiplot row index
        /// - 2nd index: multiplot column index
        /// - content: index into <paramref name="_2DData"/>
        /// </param>
        public static Gnuplot ToGnuplot(this IEnumerable<Plot2Ddata> _2DData, int[,] layout = null) {

            if(layout == null) {
                layout = new int[_2DData.Count(), 1];
                for(int iRow = 0; iRow < layout.GetLength(0); iRow++) {
                    layout[iRow, 0] = iRow;
                }
            }
   

            Gnuplot gp = new Gnuplot();

            gp.SetMultiplot(layout.GetLength(0), layout.GetLength(1));
           
            for(int iRow = 0; iRow < layout.GetLength(0); iRow++) {
                for(int iCol = 0; iCol < layout.GetLength(1); iCol++) {
                    int elem = layout[iRow, iCol];
                    if(elem >= 0 && elem < _2DData.Count()) {
                        gp.SetSubPlot(iRow, iCol);
                        _2DData.ElementAt(elem).ToGnuplot(gp);
                    }
                }
            }
            return gp;
        }


        /// <summary>
        /// <see cref="Plot2Ddata"/> into an alive Gnuplot object and executes Gnuplot interactively
        /// </summary>
        /// <param name="_2DData"></param>
        /// <param name="layout"></param>
        public static void PlotInteractive(this Plot2Ddata[,] _2DData, GnuplotPageLayout layout = null) {
            using(var gp = ToGnuplot(_2DData, layout)) {
                Console.WriteLine("Executing Gnulpot...");
                gp.Execute();

                Console.WriteLine("Press any key to continue...");
                Console.ReadKey(true);
                Console.WriteLine("killing gnuplot...");

            }
        }

        /// <summary>
        /// <see cref="Plot2Ddata"/> into an alive Gnuplot object and executes Gnuplot interactively
        /// </summary>
        public static void PlotInteractive(this Plot2Ddata _2DData, GnuplotPageLayout layout = null) {
            using(var gp = ToGnuplot(_2DData, layout)) {
                Console.WriteLine("Executing Gnulpot...");
                gp.Execute();

                Console.WriteLine("Press any key to continue...");
                Console.ReadKey(true);
                Console.WriteLine("killing gnuplot...");

            }
        }

        /// <summary>
        /// Writes a 2D multi-plot to Gnuplot and saves the 
        /// plot to a PNG file ('set terminal png').
        /// </summary>
        static public void SaveToGIF(this Plot2Ddata[,] _2DData, string OutfileName, int xRes = 800, int yRes = 600) {
            using (var gp = _2DData.ToGnuplot()) {
                gp.SaveToGIF(OutfileName, xRes, yRes);
            }
        }

        /// <summary>
        /// Writes a 2D plot to Gnuplot and saves the 
        /// plot to a PNG file ('set terminal png').
        /// </summary>
        /// <param name="_2DData">data to plot</param>
        /// <param name="xRes">Horizontal resolution in pixels.</param>
        /// <param name="yRes">Vertical resolution in pixels.</param>
        /// <param name="OutfileName">Path/filename for output PNG.</param>
        static public void SaveToGIF(this Plot2Ddata _2DData, string OutfileName, int xRes = 800, int yRes = 600) {
            using (var gp = _2DData.ToGnuplot()) {
                gp.SaveToGIF(OutfileName, xRes, yRes);
            }
        }

        /// <summary>
        /// Writes a 2D multi-plot to Gnuplot and saves the 
        /// plot to a PNG file ('set terminal svg').
        /// </summary>
        static public void SaveToSVG(this Plot2Ddata[,] _2DData, string OutfileName, GnuplotPageLayout layout = null) {
            using (var gp = _2DData.ToGnuplot(layout)) {
                gp.SaveToSVG(OutfileName);
            }
        }

        /// <summary>
        /// Writes a 2D plot to Gnuplot and saves the 
        /// plot to a PNG file ('set terminal svg').
        /// </summary>
        /// <param name="_2DData">data to plot</param>
        /// <param name="xRes">Horizontal resolution in pixels.</param>
        /// <param name="yRes">Vertical resolution in pixels.</param>
        /// <param name="OutfileName">Path/filename for output PNG.</param>
        static public void SaveToSVG(this Plot2Ddata _2DData, string OutfileName, int xRes = 800, int yRes = 600) {
            using (var gp = _2DData.ToGnuplot()) {
                gp.SaveToSVG(OutfileName, xRes, yRes);
            }
        }
    }
}
