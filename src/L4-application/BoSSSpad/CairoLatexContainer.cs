using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;

namespace BoSSS.Application.BoSSSpad {
   
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
            if(!TexInputFilePath.EndsWith(".tex"))
                throw new ArgumentException("Expecting name to end with '.tex'.");

            string OutDir = Path.GetDirectoryName(TexInputFilePath);
            string name = Path.GetFileNameWithoutExtension(TexInputFilePath);
            string tempName = Path.GetFileNameWithoutExtension(this.GraphicsFilename);
            string graphicsExt = Path.GetExtension(this.GraphicsFilename);

            string gpCode = this.GnuplotScript;
            for(int i = 0; i < this.DataFiles.Length; i++) {
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
            GraphicsFile = Path.GetFileName(GraphicsFile);
            if(GraphicsFile.IsEmptyOrWhite())
                GraphicsFile = "cairolatex.tex";

            // save graphic
            //////////////////////////////////
            this.SaveTo(Path.Combine(outDir, GraphicsFile));

            // write minimal main LaTeX file
            //////////////////////////////////
            using(var stw = new StreamWriter(MinimalExampleFilePath)) {


                stw.WriteLine(@"\documentclass{article}");
                //stw.WriteLine(@"\documentclass[preview]{standalone}");

                if(!this.PdfLatex)
                    stw.WriteLine(@"\usepackage[dvips]{graphicx}");
                else
                    stw.WriteLine(@"\usepackage{graphicx}");
                if(!this.PdfLatex) {
                    stw.WriteLine(@"\usepackage{psfrag}");
                    stw.WriteLine(@"\usepackage{epsfig}");
                }
                stw.WriteLine(@"\usepackage{amsmath}");
                stw.WriteLine(@"\usepackage{amssymb}");
                //stw.WriteLine(@"\usepackage{cmbright}");

                stw.WriteLine(@"\begin{document}");
                stw.WriteLine(@"\pagenumbering{gobble}"); // surpress page number
                stw.WriteLine(@"\begin{figure}");
                stw.WriteLine(@"\begin{center}");
                stw.WriteLine("\\input{{{0}}}", GraphicsFile);
                stw.WriteLine(@"\end{center}");
                stw.WriteLine(@"\rule{0cm}{1.1cm}");
                stw.WriteLine(@"\caption{");
                stw.WriteLine(@"\LaTeX figure caption goes here.");
                stw.WriteLine(@"}");
                stw.WriteLine(@"\end{figure}");


                stw.WriteLine(@"\end{document}");


                stw.Flush();
            }

            if(PerformLatexCompilation) {
                string LatexExe = this.PdfLatex ? "pdflatex" : "latex", DvipsExe = "dvips";
                if(Path.DirectorySeparatorChar == '\\') {
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

                    if(p.ExitCode != 0) {
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

                    if(p.ExitCode != 0) {
                        Console.WriteLine("Unable to compile, exiting.");
                        Console.WriteLine("'" + DvipsExe + "' exited with: " + p.ExitCode);
                        return;
                    }
                }
            }
        }

        /*
        /// <summary>
        /// Compiles the Latex document (to pdf or eps, depending on <see cref="CairolatexContainer.PdfLatex"/>),
        /// and converts the output to a png image, which can be previewed in (old) BoSSSpad.
        /// </summary>
        public Image PreviewPNG(bool trimPage = true, int dpi = 200) {
            if(dpi < 10 || dpi > 10000)
                throw new ArgumentOutOfRangeException();

            // get a temporary directory
            // =========================
            DirectoryInfo WorkingDirectory;
            {
                var rnd = new Random();
                bool Exists = false;
                do {
                    var tempPath = Path.GetTempPath();
                    var tempDir = rnd.Next().ToString();
                    WorkingDirectory = new DirectoryInfo(Path.Combine(tempPath, tempDir));
                    Exists = WorkingDirectory.Exists;
                    if(!Exists) {
                        WorkingDirectory.Create();
                    }
                } while(Exists == true);
            }


            // write data & compile
            // ====================
            string mainTexFile = Path.Combine(WorkingDirectory.FullName, "main.tex");
            WriteMinimalCompileableExample(mainTexFile, PerformLatexCompilation: true);



            // try to read an image
            // ====================
            if(this.PdfLatex) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Pdf Output: using 'imagemagick' for conversion to a png
                // magick.exe -verbose -density 300  main.pdf -trim PNG:main.png
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                string mainPdfFile = Path.ChangeExtension(mainTexFile, ".pdf");
                if(!File.Exists(mainPdfFile))
                    throw new IOException("Unable to find PDF output.");

                string mainPngFile = Path.ChangeExtension(mainTexFile, ".png");


                string ImageMagikTool = "magick";
                if (Path.DirectorySeparatorChar == '\\') {
                    ImageMagikTool += ".exe";
                }

                ProcessStartInfo psi = new ProcessStartInfo();
                psi.WorkingDirectory = WorkingDirectory.FullName;
                psi.Arguments = string.Format($" -verbose -density {dpi} {Path.GetFileName(mainPdfFile)} {(trimPage ? "-trim" : "")} PNG:{Path.GetFileName(mainPngFile)}");
                psi.FileName = ImageMagikTool;

                var p = new Process();
                try {
                    p.StartInfo = psi;
                    p.Start();
                    p.WaitForExit();

                    if(p.ExitCode != 0 || !File.Exists(mainPngFile))
                        throw new Exception("'" + ImageMagikTool + "' exited with: " + p.ExitCode + ",  directory is " + WorkingDirectory.FullName);
                } catch(Exception e) {

                    Console.Error.WriteLine("Unable to convert to png: " + e.Message + "  (" + e.GetType().Name + ")");
                    return null;
                }


                byte[] IOMmem = File.ReadAllBytes(mainPngFile);
                //File.Delete(mainPngFile);
                return Image.FromStream(new MemoryStream(IOMmem));

            } else {
                throw new NotImplementedException("Todo: implement conversion using 'dvipng' tool.");
            }
        }
        */

        /// <summary>
        /// Compiles the Latex document (to pdf or eps, depending on <see cref="CairolatexContainer.PdfLatex"/>),
        /// and converts the output to a SVG image, which can be previewed in Jupyter.
        /// </summary>
        /// <remarks>
        /// Requires the following programs to be correctly installed and in the PATH environment variable:
        /// - pdflatex
        /// - pdfcrop (if <paramref name="trimPage"/> is true)
        /// - pdf2svg https://github.com/dawbarton/pdf2svg (binaries for windows see: https://github.com/jalios/pdf2svg-windows)
        /// </remarks>
        public Microsoft.AspNetCore.Html.HtmlString PreviewSVG(bool trimPage = true) {
            
            // get a temporary directory
            // =========================
            DirectoryInfo WorkingDirectory;
            {
                var rnd = new Random();
                bool Exists = false;
                do {
                    var tempPath = Path.GetTempPath();
                    var tempDir = rnd.Next().ToString();
                    WorkingDirectory = new DirectoryInfo(Path.Combine(tempPath, tempDir));
                    Exists = WorkingDirectory.Exists;
                    if (!Exists) {
                        WorkingDirectory.Create();
                    }
                } while (Exists == true);
            }


            // write data & compile
            // ====================
            string mainTexFile = Path.Combine(WorkingDirectory.FullName, "main.tex");
            WriteMinimalCompileableExample(mainTexFile, PerformLatexCompilation: true);



            // try to read an image
            // ====================
            if (this.PdfLatex) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Pdf Output: using 'pdfcrop' to crop
                // using small too pdf2svg to get an svg
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                string mainPdfFile = Path.ChangeExtension(mainTexFile, ".pdf");
                if (!File.Exists(mainPdfFile))
                    throw new IOException("Unable to find PDF output.");

                string cropPdfFile = Path.Combine(Path.GetDirectoryName(mainTexFile), Path.GetFileNameWithoutExtension(mainPdfFile) + "-crop.pdf");

                string mainSvgFile = Path.ChangeExtension(mainTexFile, ".svg");

                try {

                    if (trimPage)
                        CallConversion(WorkingDirectory, "pdfcrop", $" {mainPdfFile} {cropPdfFile} ", cropPdfFile);
                    else
                        cropPdfFile = mainPdfFile;

                    CallConversion(WorkingDirectory, "pdf2svg", $" {Path.GetFileName(cropPdfFile)} {Path.GetFileName(mainSvgFile)} 1", mainSvgFile);
                //if (Path.DirectorySeparatorChar == '\\') {
                //    ConversionTool += ".exe";
                //}

                //ProcessStartInfo psi = new ProcessStartInfo();
                //psi.WorkingDirectory = WorkingDirectory.FullName;
                //psi.Arguments = string.Format($" {Path.GetFileName(mainPdfFile)} {Path.GetFileName(mainSvgFile)} 1");
                //psi.FileName = ConversionTool;

                //var p = new Process();
               
                //    p.StartInfo = psi;
                //    p.Start();
                //    p.WaitForExit();

                //    if (p.ExitCode != 0 || !File.Exists(mainSvgFile))
                //        throw new Exception("'" + ConversionTool + "' exited with: " + p.ExitCode + ",  directory is " + WorkingDirectory.FullName);


                } catch (Exception e) {

                    Console.Error.WriteLine("Unable to convert to svg: " + e.Message + "  (" + e.GetType().Name + ")");
                    return null;
                }

                // return image
                var fi = (new FileInfo(mainSvgFile));
                if (fi.Exists && fi.Length > 0) {
                    string SVGtext = File.ReadAllText(mainSvgFile);
                    File.Delete(mainSvgFile);
                    return new Microsoft.AspNetCore.Html.HtmlString(SVGtext);
                    //return Image.FromFile(OutfileName); // it seems, the image object does not work anymore when the file is deleted
                } else {
                    Console.Error.WriteLine("Gnuplot output file empty or non-existent.");
                    return null;
                }

            } else {
                throw new NotImplementedException("Todo: implement conversion using 'dvi -> ps -> pdf -> svg' tools.");
            }
        }

        static void CallConversion(DirectoryInfo WorkingDirectory, string ConversionTool, string arguments, string resultFile) {
            if (Path.DirectorySeparatorChar == '\\') {
                ConversionTool += ".exe";
            }

            ProcessStartInfo psi = new ProcessStartInfo();
            psi.WorkingDirectory = WorkingDirectory.FullName;
            psi.Arguments = arguments;
            psi.FileName = ConversionTool;



            var p = new Process();

            p.StartInfo = psi;
            p.Start();
            p.WaitForExit();

            if (p.ExitCode != 0 || !File.Exists(resultFile))
                throw new Exception("'" + ConversionTool + "' exited with: " + p.ExitCode + ",  directory is " + WorkingDirectory.FullName);



        }

    }


}
