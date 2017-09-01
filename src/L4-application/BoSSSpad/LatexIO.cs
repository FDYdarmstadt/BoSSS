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

using BoSSS.Foundation.IO;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// IO routines for storing <see cref="Document"/>-object in LaTeX-files.
    /// </summary>
    static class LatexIO {

        // "tokens"
        static readonly string TeXlink_BoSSScmd = "\\BoSSScmd";
        static readonly string TeXlink_BoSSScmdSilent = "\\BoSSScmdSilent";

        /// <summary>
        /// Parses a TeX-File and extracts LaTeX-code and C#-commands.
        /// </summary>
        /// <param name="TeXfilePath"></param>
        /// <param name="TexLines">
        /// The LaTeX code in the file.
        /// </param>
        /// <param name="doc">
        /// The C#-commands extracted from the file.
        /// </param>
        public static void SplitTexFile(string TeXfilePath, out List<string> TexLines, out Document doc)
        {


            // Load all TeX-file Lines
            // =======================
            string Tex = File.ReadAllText(TeXfilePath);


            // Remove old BoSSSexe and BoSSSexeSilent
            // =======================
            while (Tex.Contains("\\BoSSSexeSilent")) {
                var position2 = Tex.IndexOf("\\BoSSSexeSilent");
                Tex = Tex.CloneAs().Remove(position2, 15);
            }

            while (Tex.Contains("\\BoSSSexe")) {
                var position = Tex.IndexOf("\\BoSSSexe");
                Tex = Tex.CloneAs().Remove(position, 9);
            }



            // Parsing...
            // ==========
            TexLines = new List<string>();
            List<string> csCommands = new List<string>();
            string NextPart = Tex;
            while (NextPart.Length > 0) {
                string orgLine = NextPart;
                int SplitIdx1 = orgLine.IndexOf(TeXlink_BoSSScmd);
                int SplitIdx2 = orgLine.IndexOf(TeXlink_BoSSScmdSilent);
                //int SplitIdx = Math.Min(SplitIdx1, SplitIdx2
                //int SkipLength;
                int SplitIdx;
                string NextToken;
                if (SplitIdx1 < 0 && SplitIdx2 < 0) {
                    SplitIdx = -1;
                    //SkipLength = int.MinValue;
                    NextToken = null;
                }
                else if (SplitIdx1 < 0 && SplitIdx2 >= 0) {
                    SplitIdx = SplitIdx2;
                    //SkipLength = TeXlink_BoSSScmdSilent.Length;
                    NextToken = TeXlink_BoSSScmdSilent;
                }
                else if (SplitIdx2 < 0 && SplitIdx1 >= 0) {
                    SplitIdx = SplitIdx1;
                    //SkipLength = TeXlink_BoSSScmd.Length;
                    NextToken = TeXlink_BoSSScmd;
                }
                else {
                    if (SplitIdx1 < SplitIdx2) {
                        SplitIdx = SplitIdx1;
                        //SkipLength = TeXlink_BoSSScmd.Length;
                        NextToken = TeXlink_BoSSScmd;
                    }
                    else {
                        SplitIdx = SplitIdx2;
                        //SkipLength = TeXlink_BoSSScmdSilent.Length;
                        NextToken = TeXlink_BoSSScmdSilent;
                    }
                }

                if (SplitIdx >= 0) {
                    string Before = orgLine.Substring(0, SplitIdx);
                    string After = orgLine.Substring(SplitIdx);
                    if (!Before.IsEmptyOrWhite()) {
                        TexLines.Add(Before);
                        NextPart = After;
                        continue;
                    }

                    // trim the BoSSScmd token
                    Debug.Assert(After.StartsWith(NextToken));
                    After = After.Substring(NextToken.Length);

                    // searching for opening brace...
                    // now, in 'After' i am expecting either whitespaces or '{'...
                    // skipping whitespaces, searching for opening brace ...
                    After = After.TrimStart(' ', '\t', '\n', '\r');
                    if (!After.StartsWith("{")) {
                        int errLine, errCol;
                        GetParsingErrorPosition(Tex, After, out errLine, out errCol);
                        throw new ArgumentException("line " + errLine + ", column " + errCol + ": Parsing error: after \\BoSSScmd, a '{' - brace is expected.");
                    }

                    // searching for closing brace...
                    int idxClose = -1;
                    bool inComment = false;
                    for (int i = 1; i < After.Length; i++) {
                        if ((i >= 3) && After[i] == '/' && After[i - 1] == '/' && After[i - 2] == '/') {
                            // found a comment
                            inComment = true;
                        }

                        if (inComment == true && After[i] == '\n')
                            // end of comment
                            inComment = false;

                        if (!inComment) {
                            if (After[i] == '{' && After[i - 1] != '\\') {
                                int errLine, errCol;
                                GetParsingErrorPosition(Tex, After.Substring(i), out errLine, out errCol);
                                throw new ArgumentException("line " + errLine + ", column " + errCol + ": Parsing error: found un-escaped '{'-brace in \\BoSSScmd. Missing closing brace '{' or backslash '\\{'?");
                            }

                            if (After[i] == '}' && After[i - 1] != '\\') {
                                // found closing brace of \BoSSScmd - command
                                idxClose = i + 1;
                                break;
                            }
                        }
                    }

                    if (idxClose < 0) {
                        int errLine, errCol;
                        GetParsingErrorPosition(Tex, After, out errLine, out errCol);
                        throw new ArgumentException("line " + errLine + ", column " + errCol + ": Parsing error: missing closing brace '}' for \\BoSSScmd.");
                    }

                    string csCommand;
                    if (idxClose > 2) {
                        csCommand = After.Substring(1, idxClose - 2);

                        // remove leading and trailing emply lines
                        csCommand = csCommand.TrimEnd(' ', '\t', '\n', '\r');
                        List<string> csCommandLines = csCommand.Split(new char[] { '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries).ToList();
                        if (csCommandLines.Count > 0) {
                            while (csCommandLines.Count > 0 && csCommandLines[0].IsEmptyOrWhite()) {
                                csCommandLines.RemoveAt(0);
                            }
                            using (var stw = new StringWriter()) {
                                for (int z1 = 0; z1 < csCommandLines.Count - 1; z1++)
                                    stw.WriteLine(csCommandLines[z1]);
                                stw.Write(csCommandLines.Last());
                                csCommand = stw.ToString();
                            }
                        }
                        else {
                            csCommand = "";
                        }
                    }
                    else {
                        csCommand = "";
                    }
                    csCommands.Add(csCommand);

                    // proceed with the rest after the closing brace...
                    NextPart = After.Substring(idxClose);
                    if (NextPart.Length >= 0) {
                        NextPart = NextPart.TrimStart(' ', '\t', '\n', '\r');
                    }
                }
                else {
                    TexLines.Add(NextPart);
                    NextPart = "";
                }
            }

            if (TexLines.Count == 0)
                TexLines.Add("");

            // create doc
            // ==========
            doc = new Document();
            doc.CommandAndResult = csCommands.Select(csCommand => new Document.Tuple() { Command = csCommand, InterpreterTextOutput = "" }).ToList();

            // convert all TeX-escaped C# to real C#:
            // ======================================
            foreach (var t in doc.CommandAndResult) {
                string Command_TexEsc = t.Command;
                t.Command = Tex2Bws(Command_TexEsc);
            }
        }


        /// <summary>
        /// Integrates the content of a document into a LaTeX-file, using commands provided by package 'wsDBEtex.sty'
        /// </summary>
        public static void UpdateTexFile(string TeXfilePath, Document doc)
        {

            // open TeX-file
            // =============

            // first, create a backup of the TeX-file
            File.Copy(TeXfilePath, TeXfilePath + "~", true); // emacs convention.

            // open file
            Document dummy;
            List<string> TexLines;
            SplitTexFile(TeXfilePath, out TexLines, out dummy);
            dummy = null;

            // insert new commands into TeX
            // ============================
            using (var wrt = new StreamWriter(new FileStream(TeXfilePath, FileMode.Create, FileAccess.Write))) {
                int NoOfCmds = doc.CommandAndResult.Count;
                while (NoOfCmds >= 0 && doc.CommandAndResult[NoOfCmds - 1].Command.IsEmptyOrWhite())
                    NoOfCmds--;

                // all boxes before restart should get an BoSSSexeSilent
                bool restartOcurr = false;
                for (int i = 0; i < Math.Max(NoOfCmds, TexLines.Count); i++) {
                    bool haveTex = i < TexLines.Count;
                    bool haveCmd = i < NoOfCmds;

                    if (haveTex) {

                        wrt.Write(TexLines[i]);
                    }
                    if (haveCmd) {
                        string TexEscapedCommand = Bws2Tex(doc.CommandAndResult[i].Command);

                        if (TexEscapedCommand.Contains("BoSSScmdSilent")) {
                            wrt.WriteLine(TeXlink_BoSSScmdSilent + "{");
                            wrt.WriteLine(TexEscapedCommand);
                            wrt.WriteLine(" }");
                        }
                        else {
                            wrt.WriteLine(TeXlink_BoSSScmd + "{");
                            wrt.WriteLine(TexEscapedCommand);
                            wrt.WriteLine(" }");
                        }

                        // BoSSSexeSilent for restart and for all boxes before box 'restart'
                        if ((TexEscapedCommand == "restart") || (TexEscapedCommand =="restart;")) { 
                            wrt.WriteLine("\\BoSSSexeSilent");
                            restartOcurr = true;
                        }else if (!restartOcurr) {
                            wrt.WriteLine("\\BoSSSexeSilent");
                        }
                        else {
                            if (!haveTex){
                                if (TexEscapedCommand.Contains("BoSSSexeSilent")) {
                                    wrt.WriteLine("\\BoSSSexeSilent");
                                }
                                else {
                                    wrt.WriteLine("\\BoSSSexe");
                                }
                            }
                            else { wrt.WriteLine("\\BoSSSexe"); }
                        }
                    }
                }
            }

            // save texbatch output
            // ====================
            {
                string otdir = Path.GetDirectoryName(TeXfilePath);
                string docnm = Path.GetFileNameWithoutExtension(TeXfilePath) + ".texbatch";
                Save_Texbatch(otdir, docnm, doc);
            }
        }

        static void GetParsingErrorPosition(string org, string rest, out int iLine, out int iCol)
        {
            int idx = org.IndexOf(rest);
            if (idx <= 0) {
                iLine = -1;
                iCol = -1;
                return;
            }

            iLine = 1;
            iCol = 0;
            for (int i = 0; i < idx; i++) {
                if (org[i] != '\r')
                    iCol++;
                if (org[i] == '\n') {
                    iCol = 0;
                    iLine++;
                }
            }
        }


        static bool IsEmptyLinePlaceHolder(string s)
        {
            int L = s.Length;
            for (int l = 0; l < L; l++) {
                char sl = s[l];
                if ((!char.IsWhiteSpace(sl)) && (sl != '%'))
                    return false;
            }

            return true;
        }


        /// <summary>
        /// Converts 'LaTeX-escaped C#' charaters to real C#
        /// </summary>
        public static string Tex2Bws(string CompleteCommand)
        {
            // split Command into multiple lines.
            string[] CommandLines = CompleteCommand.Split(new string[] { System.Environment.NewLine }, StringSplitOptions.None);

            using (var stw = new StringWriter()) {
                for (int k = 0; k < CommandLines.Length; k++) {
                    string Command = CommandLines[k];
                    string lineWoW = Command.TrimStart(' ', '\t');
                    if (lineWoW.StartsWith("///")) {
                        // LaTeX - verbatim: no replacement


                    }
                    else if (IsEmptyLinePlaceHolder(Command)) {
                        Command = "";
                    }
                    else {

                        // a bunch of replacements to embedd C# into LaTeX:
                        Command = Command.Replace("\\rule {0.5cm}{0.0cm}", "   ");
                        Command = Command.Replace("\\btab", "   ");
                        Command = Command.Replace("\\newline ", "");
                        Command = Command.Replace("\\_", "_");
                        Command = Command.Replace("\\%", "%");
                        Command = Command.Replace("\\par", System.Environment.NewLine);
                        Command = Command.Replace("\\{", "{");
                        Command = Command.Replace("\\}", "}");
                        Command = Command.Replace("\\textbackslash ", "\\");
                    }

                    if (k < CommandLines.Length - 1)
                        stw.WriteLine(Command);
                    else
                        stw.Write(Command);

                }

                return stw.ToString();
            }
        }

        /// <summary>
        /// Converts C# to 'LaTeX-escaped C#'
        /// </summary>
        public static string Bws2Tex(string Command)
        {

            if (Command == null)
                Command = "";

            // split Command into multiple lines.
            string[] CommandLines = Command.Split(new string[] { System.Environment.NewLine }, StringSplitOptions.None);

            // remove trailing empty lines 
            while (CommandLines.Length > 0) {
                if (CommandLines[CommandLines.Length - 1].IsEmptyOrWhite()) {
                    CommandLines = ArrayTools.GetSubVector(CommandLines, 0, CommandLines.Length - 1);
                }
                else {
                    break;
                }
            }

            // write command
            using (var stw = new StringWriter()) {
                for (int k = 0; k < CommandLines.Length; k++) {
                    string ln = CommandLines[k];
                    string lineWoW = ln.TrimStart(' ', '\t');
                    if (lineWoW.StartsWith("///")) {
                        // ++++++++++++++++
                        // LaTeX - verbatim
                        // ++++++++++++++++

                        // write line
                        if (k < CommandLines.Length - 1)
                            stw.WriteLine(ln);
                        else
                            stw.Write(ln);
                    }
                    else {

                        // global operations
                        // =================
                        ln = ln.Replace("\\", "\\textbackslash ");
                        ln = ln.Replace("{", "\\{");
                        ln = ln.Replace("}", "\\}");
                        ln = ln.Replace("\t", "\btab ");
                        ln = ln.Replace("\r", "");
                        //ln = ln.Replace("\n", "\\newline ");
                        ln = ln.Replace("_", "\\_");
                        ln = ln.Replace("%", "\\%");

                        // per-line operations
                        // ===================


                        // leading tabs/white-spaces: replace with "\btab"
                        //for (int i = 0; i < CommandLines.Length; i++) {
                        if (ln.IsEmptyOrWhite()) {
                            ln = " ";
                        }

                        int NoOfIntendTabs = 0;
                        while (ln.StartsWith("    ")) {
                            ln = ln.Substring(4);
                            NoOfIntendTabs++;
                        }

                        for (int j = 0; j < NoOfIntendTabs; j++) {
                            ln = "\\btab " + ln;
                        }
                        //}

                        // write line
                        if (k < CommandLines.Length - 1)
                            stw.WriteLine(ln + "\\newline ");
                        else
                            stw.Write(ln);
                    }

                }

                // return
                // ======
                string ret = stw.ToString();
                if (ret.IsEmptyOrWhite())
                    // completly empty commands
                    return " % ";
                else
                    return ret;
            }
        }

        static string lstset_commands =
            @"basicstyle=\footnotesize\ttfamily\color{darkblue}, " +
            @"language=C, " +
            @"keywordstyle=\color{darkgreen}, " +
            @"stringstyle=\color{darkgrey}, " +
            @"commentstyle=\color{lightergray}";


        /// <summary>
        /// Saves document in a '.texbatch'--directory so that it can be imported into LaTeX.
        /// </summary>
        public static void Save_Texbatch(string _OutDir, string DocName, Document doc)
        {

            if (!DocName.EndsWith(".texbatch"))
                throw new ArgumentException();
            string OutDir = Path.Combine(_OutDir, DocName);

            // ==================================
            // create output dir, if not existent
            // ==================================
            if (!Directory.Exists(OutDir)) {
                Directory.CreateDirectory(OutDir);
            }
            if (!Directory.Exists(OutDir)) {
                throw new IOException("Unable to create output directory: '" + OutDir + "'");
            }

            DirectoryInfo di = new DirectoryInfo(OutDir);


            // ================
            // clear output dir
            // ================
            string[] files = Directory.GetFiles(OutDir);
            foreach (var f in files)
                File.Delete(f);

            // ===================================
            // save every entry in a separate file
            // ===================================
            int linecounter = 1;
            for (int iEntry = 0; iEntry < doc.CommandAndResult.Count; iEntry++) {
                var Entry = doc.CommandAndResult[iEntry];

                // --------------------------------
                // write the C#-statement / in-file
                // --------------------------------


                string inName = Path.Combine(OutDir, "in" + iEntry + ".tex");
                //File.WriteAllText(inName, Entry.Command);
                using (var inFile = new StreamWriter(inName)) {
                    using (var cmdReader = new StringReader(Entry.Command != null ? Entry.Command : "")) {
                        int iEntrySub = 0;
                        StreamWriter inFileSub = null;
                        for (string line = cmdReader.ReadLine(); line != null; line = cmdReader.ReadLine()) {
                            string lineWoW = line.TrimStart(' ', '\t');
                            if (lineWoW.StartsWith("///")) {
                                if (inFileSub != null) {
                                    inFileSub.Flush();
                                    inFileSub.Close();
                                    inFileSub.Dispose();
                                    inFileSub = null;
                                }

                                lineWoW = lineWoW.Substring(3);
                                inFile.WriteLine(lineWoW);
                            }
                            else {
                                if (inFileSub == null) {
                                    string inNameSub = "in" + iEntry + "-" + iEntrySub + ".txt";
                                    iEntrySub++;
                                    inFile.WriteLine(@"\lstinputlisting[" + lstset_commands + "]{" + DocName + "/" + inNameSub + "}");
                                    inFileSub = new StreamWriter(Path.Combine(OutDir, inNameSub));
                                }
                                inFileSub.Write("{0}: ", linecounter);
                                inFileSub.WriteLine(line);
                                linecounter++;
                            }
                        }

                        if (inFileSub != null) {
                            inFileSub.Flush();
                            inFileSub.Close();
                            inFileSub.Dispose();
                            inFileSub = null;
                        }
                    }
                }

                // ------------------------------------
                // write result / out-file
                // ------------------------------------

                if (Entry.Result is GnuplotExtensions.CairolatexContainer) {
                    // graphical output
                    var CairoHack = (GnuplotExtensions.CairolatexContainer)Entry.Result;
                    //File.WriteAllText(Path.Combine(OutDir, "out" + iEntry + ".tex"), CairoHack.LatexCode);
                    File.WriteAllText(Path.Combine(OutDir, "empty" + iEntry + ".txt"), "");
                    //File.WriteAllBytes(Path.Combine(OutDir, CairoHack.GraphicsFilename), CairoHack.GraphicsData);
                    CairoHack.SaveTo(Path.Combine(OutDir, "out" + iEntry + ".tex"));
                }
                else {
                    if (Entry.InterpreterTextOutput.IsEmptyOrWhite()) {
                        File.WriteAllText(Path.Combine(OutDir, "empty" + iEntry + ".txt"), "");
                    }
                    else {
                        // text output
                        //File.WriteAllText(Path.Combine(OutDir, "out" + iEntry + ".txt"), Entry.InterpreterTextOutput);
                        using (var outReader = new StringReader(Entry.InterpreterTextOutput)) {
                            using (var outFile = new StreamWriter(Path.Combine(OutDir, "out" + iEntry + ".txt"))) {
                                for (string line = outReader.ReadLine(); line != null; line = outReader.ReadLine()) {
                                    outFile.WriteLine("> {0}", line);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
