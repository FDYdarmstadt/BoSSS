using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Diagnostics;
using System.Globalization;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.BoSSSpad{

    /// <summary>
    /// Entrypoint used by <see cref="ElectronWorksheet"/> project. 
    /// Realizes communication between electron BoSSSpad and C# BoSSSpad.
    /// </summary>
    public sealed class ElectronWorksheet : ResolvableAssembly {

        /// <summary>
        /// Will only work for one instance
        /// </summary>
        /// <param name="BoSSSpath">
        /// Path to the ElectronWorksheet.dll, ElectronBoSSSpad.exe and affiliated DLLs
        /// </param>
        public ElectronWorksheet(string BoSSSpath) : base(BoSSSpath){
            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            //Setup Environment
            ilPSP.Environment.Bootstrap(
                new string[0],
                Utils.GetBoSSSInstallDir(),
                out bool mpiInitialized
            );
        }

        public Tuple<string, string> RunCommand(string command) {
            
            Document.Tuple singleCommandAndResult = new Document.Tuple {
                Command = command
            };
            singleCommandAndResult.Evaluate();
            String base64Result = null;
            
            if (singleCommandAndResult.Result != null 
                && singleCommandAndResult.Result as System.Drawing.Image != null)
            {
                Byte[] result = null;
                System.Drawing.Image img = (System.Drawing.Image)singleCommandAndResult.Result;
                using (System.IO.MemoryStream ms = new System.IO.MemoryStream())
                {
                    img.Save(ms, System.Drawing.Imaging.ImageFormat.Png);
                    result = ms.ToArray();
                    base64Result = Convert.ToBase64String(result);
                };
            }
            
            return new Tuple<string, string>(
                singleCommandAndResult.InterpreterTextOutput,
                base64Result);
        }

        public void Save(string path, string[] commands, string[] results) {
            //build document 
            Document document = new Document();
            for (int i = 0; i < commands.Length; ++i)
            {
                Document.Tuple commandBox = new Document.Tuple()
                {
                    Command = commands[i],
                    InterpreterTextOutput = results[i]
                };
                document.CommandAndResult.Add(commandBox);
            }

            //Save document
            document.Serialize(path);
        }

        public Tuple<string[], string[]> Load(string path){

            Document document = Document.Deserialize(path);
            int numberOfBoxes = document.CommandAndResult.Count;
            string[] commands = new string[numberOfBoxes];
            string[] results = new string[numberOfBoxes];
            for(int i = 0; i < numberOfBoxes; ++i){
                commands[i] = document.CommandAndResult[i].Command;
                results[i] = document.CommandAndResult[i].InterpreterTextOutput;
            }

            return new Tuple<string[], string[]>(commands, results);
        }

        public string[] GetAutoCompleteSuggestions(string textToBeCompleted){
            
            string[] completions = null;
            string originalPrefix = null;
            int timeout = 1000;

            if (ReadEvalPrintLoop.eval != null){
                bool completed = ReadEvalPrintLoop.eval.TryGetCompletions(
                    textToBeCompleted, out completions, out originalPrefix, timeout);   
            }

            if (completions != null){
                for (int i = 0; i < completions.Length; ++i){
                    completions[i] = originalPrefix + completions[i];
                }
            }

            return completions;
        }
    }
}

