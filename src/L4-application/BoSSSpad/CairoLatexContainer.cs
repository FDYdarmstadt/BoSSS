using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

        
        /// <summary>
        /// Compiles the Latex document (to pdf or eps, depending on <see cref="CairolatexContainer.PdfLatex"/>),
        /// and converts the output to a png image, which can be previewed in BoSSSpad.
        /// </summary>
        public Image Preview(bool trimPage = true, int dpi = 200) {
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


                string ImageMagikTool = "magick.exe";

                ProcessStartInfo psi = new ProcessStartInfo();
                psi.WorkingDirectory = WorkingDirectory.FullName;
                psi.Arguments = string.Format(" -verbose -density {3} {0} {2} PNG:{1}", Path.GetFileName(mainPdfFile), Path.GetFileName(mainPngFile), trimPage ? "-trim" : "", dpi);
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
        
    }


}
