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
    /// Spice up Gnuplot.
    /// </summary>
    public static class GnuplotExtensions {

        /// <summary>
        /// Plot to a gif file ('set terminal gif').
        /// </summary>
        /// <param name="gp"></param>
        /// <param name="xRes">Horizontal resolution in pixels.</param>
        /// <param name="yRes">Vertical resolution in pixels.</param>
        static public Image PlotGIF(this Gnuplot gp, int xRes = 800, int yRes = 600) {

            // set terminal
            gp.Terminal = string.Format("pngcairo size {0},{1}", xRes, yRes);

            // set output file
            //string OutfileName = null;
            //while (OutfileName == null || File.Exists(OutfileName)) {
            //    OutfileName = Path.GetTempFileName() + ".gif";
            //}
            string OutfileName = Path.GetTempFileName();
            gp.OutputFile = OutfileName;

            // call gnuplot
            int exCode = gp.RunAndExit(); // run & close gnuplot
            if (exCode != 0) {
                Console.WriteLine("Gnuplot-internal error: exit code " + exCode);
                return null;
            }


            // return image
            var fi = (new FileInfo(OutfileName));
            if (fi.Exists && fi.Length > 0) {
                return Image.FromFile(OutfileName);
            } else {
                Console.WriteLine("Gnuplot output file empty or non-existent.");
                return null;
            }
        }

        /// <summary>
        /// Memory-image of gnuplot cairolatex output.
        /// </summary>
        [Serializable]
        public class CairolatexContainer {

            /// <summary>
            /// Content of data files used in <see cref="GnuplotScript"/>.
            /// </summary>
            public string[] DataFiles;

            /// <summary>
            /// Names of <see cref="DataFiles"/>.
            /// </summary>
            public string[] DataFileNames;
            
            /// <summary>
            /// The gnuplot source code.
            /// </summary>
            public string GnuplotScript;
            
            /// <summary>
            /// LaTeX-code which is read via the '\input'
            /// </summary>
            public string LatexCode;

            /// <summary>
            /// The content of the pdf or eps file.
            /// </summary>
            public byte[] GraphicsData;

            /// <summary>
            /// The filename of the pdf or eps file.
            /// </summary>
            public string GraphicsFilename;

            /// <summary>
            /// Saves the LaTeX-code and graphics data (<see cref="GraphicsData"/>) to a specific
            /// file/directory.
            /// </summary>
            /// <param name="TexInputFilePath">
            /// Name and path of LaTeX-file, which can be imported into a main file via the 
            /// `\input{...}`-command in LaTeX, e.g. `C:\Users\someone\grafik.tex`.
            /// </param>
            public void SaveTo(string TexInputFilePath) {
                if (!TexInputFilePath.EndsWith(".tex"))
                    throw new ArgumentException("Expecting name to end with '.tex'.");
                
                string OutDir = Path.GetDirectoryName(TexInputFilePath);
                string name = Path.GetFileNameWithoutExtension(TexInputFilePath);
                string tempName = Path.GetFileNameWithoutExtension(this.GraphicsFilename);
                string graphicsExt = Path.GetExtension(this.GraphicsFilename);

                string gpCode = this.GnuplotScript;
                for (int i = 0; i < this.DataFiles.Length; i++) {
                    string csvTempName = this.DataFileNames[i];
                    Debug.Assert(csvTempName.Contains(tempName));

                    string csvName = name + "_data_" + i + ".csv";
                    gpCode = gpCode.Replace(csvTempName, csvName);

                    File.WriteAllText(Path.Combine(OutDir, csvName), this.DataFiles[i]);
                }
                gpCode = gpCode.Replace(tempName + ".tex", name + ".tex");

                File.WriteAllText(Path.Combine(OutDir, name + ".gp"), gpCode);

                File.WriteAllText(TexInputFilePath, this.LatexCode.Replace(tempName, name));
                File.WriteAllBytes(Path.Combine(OutDir, name + graphicsExt), this.GraphicsData);

            }

            /// <summary>
            /// `True`, if the graphics data <see cref="GraphicsData"/> represents a pdf-file, `false` if it represents an eps-file.
            /// </summary>
            public bool PdfLatex {
                get {
                    return this.GraphicsFilename.ToLowerInvariant().EndsWith(".pdf");
                }
            }


            /// <summary>
            /// 
            /// </summary>
            /// <param name="MinimalExampleFilePath">
            /// Path and filename of the LaTeX-main file, e.g. `c:\temp\miniex.tex`.
            /// </param>
            /// <param name="GraphicsFile">
            /// Name of the Cairolatex file, e.g. `figure.tex`. 
            /// </param>
            /// <param name="PerformLatexCompilation">
            /// If true, the `latex` or `pdflatex` compiler will be called.
            /// </param>
            public void WriteMinimalCompileableExample(string MinimalExampleFilePath, string GraphicsFile = null, bool PerformLatexCompilation = true) {
                string outDir = Path.GetDirectoryName(MinimalExampleFilePath);
                string fileName = Path.GetFileNameWithoutExtension(MinimalExampleFilePath);
                if(GraphicsFile.IsEmptyOrWhite())
                    GraphicsFile = "cairolatex.tex";

                // save graphic
                //////////////////////////////////
                this.SaveTo(Path.Combine(outDir, GraphicsFile));
                
                // write minimal main LaTeX file
                //////////////////////////////////
                using (var stw = new StreamWriter(MinimalExampleFilePath)) {


                    stw.WriteLine(@"\documentclass{article}");

                    if (!this.PdfLatex)
                        stw.WriteLine(@"\usepackage[dvips]{graphicx}");
                    else
                        stw.WriteLine(@"\usepackage{graphicx}");
                    if (!this.PdfLatex) {
                        stw.WriteLine(@"\usepackage{psfrag}");
                        stw.WriteLine(@"\usepackage{epsfig}");
                    }
                    stw.WriteLine(@"\usepackage{amsmath}");
                    stw.WriteLine(@"\usepackage{amssymb}");
                    //stw.WriteLine(@"\usepackage{cmbright}");
                                        
                    stw.WriteLine(@"\begin{document}");
                    stw.WriteLine(@"\begin{figure}");
                    stw.WriteLine(@"\begin{center}");
                    stw.WriteLine("\\input{{{0}}}", GraphicsFile);
                    stw.WriteLine(@"\end{center}");
                    stw.WriteLine(@"\rule{0cm}{1.1cm}");
                    stw.WriteLine(@"\caption{");
                    stw.WriteLine(@"ToDo");
                    stw.WriteLine(@"}");
                    stw.WriteLine(@"\end{figure}");


                    stw.WriteLine(@"\end{document}");


                    stw.Flush();
                }

                string LatexExe = this.PdfLatex ? "pdflatex" : "latex", DvipsExe = "dvips";
                if (Path.DirectorySeparatorChar == '\\') {
                    LatexExe += ".exe";
                }


                // compile LaTeX
                //////////////////////
                {
                    ProcessStartInfo psi = new ProcessStartInfo();
                    psi.WorkingDirectory = outDir;
                    psi.Arguments = fileName + ".tex";
                    psi.FileName = LatexExe;


                    var p = new Process();
                    p.StartInfo = psi;
                    p.Start();
                    p.WaitForExit();

                    if (p.ExitCode != 0) {
                        Console.WriteLine("Unable to compile, exiting.");
                        Console.WriteLine("'" + LatexExe + "' exited with: " + p.ExitCode);
                        return;
                    }
                }

                // use dvips
                ////////////////////////
                if(!this.PdfLatex) {
                    ProcessStartInfo psi = new ProcessStartInfo();
                    psi.WorkingDirectory = outDir;
                    psi.Arguments = fileName + ".dvi";
                    psi.FileName = DvipsExe;


                    var p = new Process();
                    p.StartInfo = psi;
                    p.Start();
                    p.WaitForExit();

                    if (p.ExitCode != 0) {
                        Console.WriteLine("Unable to compile, exiting.");
                        Console.WriteLine("'" + DvipsExe + "' exited with: " + p.ExitCode);
                        return;
                    }
                }
            }
        }

        /// <summary>
        /// Plot to a gif file ('set terminal gif').
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

        /// <summary>
        /// Hack for <see cref="PlotNow(Gnuplot)"/>
        /// </summary>
        public static bool UseCairoLatex = false;

        /// <summary>
        /// Gnuplot plotting, automatic choice of gnuplot driver depending on
        /// the current value of <see cref="UseCairoLatex"/>.
        /// </summary>
        /// <param name="gp"></param>
        /// <returns></returns>
        public static object PlotNow(this Gnuplot gp) {

            if (UseCairoLatex) {
                return gp.PlotCairolatex();
            } else {
                return gp.PlotGIF();
            }
        }
    }
}
