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

using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// storing 
    /// </summary>
    [Serializable]
    class Document {

        public Document() {
            CommandAndResult = new List<Tuple>();
        }

        /// <summary>
        /// Rem.: the items of this tuple are mutable 
        /// </summary>
        [Serializable]
        public class Tuple {

            /// <summary>
            /// C# command entered by the user
            /// </summary>
            public string Command;

            /// <summary>
            /// The interpreter text output.
            /// </summary>
            public string InterpreterTextOutput;

            /// <summary>
            /// The actual result of the evaluation.
            /// </summary>
            [NonSerialized]
            public object Result;

            /// <summary>
            /// 
            /// </summary>
            [NonSerialized]
            public Assembly AssemblyProduced;

            /// <summary>
            /// Evaluates this command and updated <see cref="InterpreterTextOutput"/>.
            /// </summary>
            public bool Evaluate() {
                StringWriter stw = new StringWriter();
                bool supOut = ilPSP.Environment.StdOut.surpressStream0;
                //ilPSP.Environment.StdOut.surpressStream0 = true;
                bool superr = ilPSP.Environment.StdErr.surpressStream0;
                //ilPSP.Environment.StdErr.surpressStream0 = true;

                ilPSP.Environment.StdOut.WriterS.Add(stw);
                ilPSP.Environment.StdErr.WriterS.Add(stw);

                if (this.Command != null && this.Command.Length > 0) {
                    this.Result = ReadEvalPrintLoop.EvalPrint(this.Command, out AssemblyProduced);
                } else {
                    this.Result = null;
                }

                Console.Out.Flush();
                Console.Error.Flush();

                this.InterpreterTextOutput = stw.ToString();
                ilPSP.Environment.StdOut.WriterS.Remove(stw);
                ilPSP.Environment.StdErr.WriterS.Remove(stw);

                //ilPSP.Environment.StdOut.surpressStream0 = supOut;
                //ilPSP.Environment.StdErr.surpressStream0 = superr;

                return (InteractiveShell.LastError == null 
                    && ReadEvalPrintLoop.cmpCont != null // wenn das nix ist eh irgendwas oberfaul
                    && ReadEvalPrintLoop.cmpCont.Report.Errors == 0);

            }
        }

        /// <summary>
        /// A list of command/result pairs that make up the document. <br/>
        /// 1st item: C# command <br/>
        /// 2nd item: interpreter answer
        /// </summary>
        public List<Tuple> CommandAndResult;

        /// <summary>
        /// Marks the beginning of a result in a worksheet-textfile.
        /// </summary>
        readonly static public string ResultStartMarker = "**************";


        /// <summary>
        /// Marks the ending of a result in a worksheet-textfile.
        /// </summary>
        readonly static public string ResultEndMarker = "==============";


        /// <summary>
        /// Saves this document to a file.
        /// </summary>
        public void Serialize(string DocumentPath) {

            using (StreamWriter fs = File.CreateText(DocumentPath)) {
                //fs.Write(JsonConvert.SerializeObject(this));

                foreach (Tuple t in CommandAndResult) {
                    if (t.Command != null && t.Command.Length > 0)
                        fs.WriteLine(t.Command);
                    fs.WriteLine(ResultStartMarker);
                    if (t.InterpreterTextOutput != null && t.InterpreterTextOutput.Length > 0)
                        fs.WriteLine(t.InterpreterTextOutput);
                    fs.WriteLine(ResultEndMarker);
                }
            }
        }

        /// <summary>
        /// Loads a document from a file.
        /// </summary>
        public static Document Deserialize(string DocumentPath) {


            using (StreamReader fs = File.OpenText(DocumentPath)) {

                List<Tuple> CRlist = new List<Tuple>();
                Tuple current = new Tuple();
                StringWriter currentWriter = new StringWriter();

                bool ReadingResult = false;
                for (string line = fs.ReadLine(); line != null; line = fs.ReadLine()) {
                    if (line.StartsWith(ResultStartMarker)) {
                        Debug.Assert(current.Command == null);
                        currentWriter.Flush();
                        current.Command = currentWriter.ToString();
                        currentWriter.Dispose();
                        currentWriter = new StringWriter();

                        ReadingResult = true;
                    } else if (line.StartsWith(ResultEndMarker)) {
                        Debug.Assert(current.InterpreterTextOutput == null);
                        currentWriter.Flush();
                        current.InterpreterTextOutput = currentWriter.ToString();
                        currentWriter.Dispose();
                        currentWriter = new StringWriter();

                        ReadingResult = false;
                        CRlist.Add(current);
                        current = new Tuple();
                    } else {
                        if (currentWriter.ToString().Length > 0)
                            currentWriter.WriteLine();
                        currentWriter.Write(line);

                        //if (!ReadingResult) {
                        //    if (current.Command == null)
                        //        current.Command = "";
                        //    current.Command = current.Command + "\n" + line;
                        //} else {
                        //    if (current.Result == null)
                        //        current.Result = "";
                        //    current.Result = current.Result + "\n" + line;
                        //}
                    }
                }

                if (CRlist.Count <= 0)
                    // document must have at least one entry
                    CRlist.Add(current);
                currentWriter.Flush();
                if (currentWriter.ToString().Length >= 0) {
                    // unfinished tuple at end of file.
                    CRlist.Add(current);

                    if (!ReadingResult) {
                        Debug.Assert(current.Command == null);
                        current.Command = currentWriter.ToString();
                    } else {
                        current.InterpreterTextOutput = currentWriter.ToString();
                    }
                }
                currentWriter.Dispose();

                Document R = new Document();
                R.CommandAndResult = CRlist;
                return R;
            }
        }
    }

}

