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
//using Microsoft.AspNetCore.Html;

namespace BoSSS.Application.BoSSSpad {


    /// <summary>
    /// Spice up <see cref="Plot2Ddata"/>
    /// </summary>
    public static class GnuplotExtensions {
        
        /*
        /// <summary>
        /// Plot to a gif file (`set terminal png`).
        /// </summary>
        /// <param name="gp"></param>
        /// <param name="xRes">Horizontal resolution in pixels.</param>
        /// <param name="yRes">Vertical resolution in pixels.</param>
        static public Image PlotGIF(this Gnuplot gp, int xRes = 800, int yRes = 600) {

          
            // set output file
            string OutfileName = Path.GetTempFileName();
            
            try {
                gp.SaveToGIF(OutfileName, xRes, yRes);
            } catch (IOException ioe) {
                Console.Error.WriteLine($"Gnuplot wrapper has thrown {ioe.GetType().Name}: {ioe.Message}");
            }


            // return image
            var fi = (new FileInfo(OutfileName));
            if (fi.Exists && fi.Length > 0) {
                byte[] IOMmem = File.ReadAllBytes(OutfileName);
                File.Delete(OutfileName);
                return Image.FromStream(new MemoryStream(IOMmem));
                //return Image.FromFile(OutfileName); // it seems, the image object does not work anymore when the file is deleted
            } else {
                Console.Error.WriteLine("Gnuplot output file empty or non-existent.");
                return null;
            }
        }//*/
        
        /// <summary>
        /// Plot to a scalable vector graphics file (`set terminal svg`).
        /// </summary>
        /// <param name="gp"></param>
        /// <param name="xRes">Horizontal resolution in pixels.</param>
        /// <param name="yRes">Vertical resolution in pixels.</param>
        static public Microsoft.AspNetCore.Html.HtmlString PlotSVG(this Gnuplot gp, int xRes = 800, int yRes = 600) {
            string OutfileName = Path.GetTempFileName();

            Console.WriteLine("Note: In a Jupyter Worksheet, you must NOT have a trailing semicolon in order to see the plot on screen; otherwise, the output migth be surpressed.!");

            try {
                gp.SaveToSVG(OutfileName, xRes, yRes);
            } catch (IOException ioe) {
                Console.Error.WriteLine($"Gnuplot wrapper has thrown {ioe.GetType().Name}: {ioe.Message}");
            }

            // return image
            var fi = (new FileInfo(OutfileName));
            if (fi.Exists && fi.Length > 0) {
                string SVGtext = File.ReadAllText(OutfileName);
                File.Delete(OutfileName);
                return new Microsoft.AspNetCore.Html.HtmlString(SVGtext);
                //return Image.FromFile(OutfileName); // it seems, the image object does not work anymore when the file is deleted
            } else {
                Console.Error.WriteLine("Gnuplot output file empty or non-existent.");
                return null;
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
                gp.Terminal = string.Format("cairolatex {0}", Options != null ? Options : " ");

            // set output file
            string baseName = Path.GetTempFileName();
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

        /* moved to L3;

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

        */
    }
}
