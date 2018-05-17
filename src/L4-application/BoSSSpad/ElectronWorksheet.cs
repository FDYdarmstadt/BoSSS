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
    /// Singleton class; 
    /// </summary>
    public sealed class ElectronWorksheet {
        Document document;
        private static readonly ElectronWorksheet instance = new ElectronWorksheet();

        private ElectronWorksheet() {

            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            // launch the app
            // ==============
            ilPSP.Environment.Bootstrap(
                new string[0],
                Utils.GetBoSSSInstallDir(),
                out bool mpiInitialized
            );
        }

        public static ElectronWorksheet Instance {
            get {
                return instance;
            }
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
                using (System.IO.MemoryStream ms = new System.IO.MemoryStream()) {
                    img.Save(ms, System.Drawing.Imaging.ImageFormat.Png);
                    result = ms.ToArray();
                    base64Result = Convert.ToBase64String(result);
                };
            }
            
            return new Tuple<string, string>(
                singleCommandAndResult.InterpreterTextOutput,
                base64Result
                );
        }

        public void Save(string path, string[] commands, string[] results) {
            //build document 
            document = new Document();
            for (int i = 0; i < commands.Length; ++i) {
                Document.Tuple commandBox = new Document.Tuple() {
                    Command = commands[i],
                    InterpreterTextOutput = results[i]
                };
                document.CommandAndResult.Add(commandBox);
            }

            //Save document
            document.Serialize(path);
        }

        public Tuple<string[], string[]> Load(string path){

            document = Document.Deserialize(path);
            int numberOfBoxes = document.CommandAndResult.Count;
            string[] commands = new string[numberOfBoxes];
            string[] results = new string[numberOfBoxes];
            for(int i = 0; i < numberOfBoxes; ++i){
                commands[i] = document.CommandAndResult[i].Command;
                results[i] = document.CommandAndResult[i].InterpreterTextOutput;
            }

            return new Tuple<string[], string[]>(commands, results);
        }
    }
}
