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
<<<<<<< HEAD
    /// Entrypoint used by <see cref="ElectronWorksheet"/> project to communicate between electron BoSSSpad and C# BoSSSpad   
    /// </summary>
    public sealed class ElectronWorksheet {

        /// <summary>
        /// Will only work for one instance
        /// </summary>
        /// <param name="BoSSSpath">
        /// Path to the ElectronWorksheet.dll, ElectronBoSSSpad.exe and affiliated DLLs
        /// </param>
        public ElectronWorksheet(string BoSSSpath) {
            path = BoSSSpath;
=======
    /// Singleton class; 
    /// </summary>
    public sealed class ElectronWorksheet {

        private static readonly ElectronWorksheet instance = new ElectronWorksheet();

        private ElectronWorksheet() {

>>>>>>> 6fa3faeb05dbc80532ef484f3623d6be69e8da96
            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            // launch the app
            // ==============
            ilPSP.Environment.Bootstrap(
                new string[0],
                Utils.GetBoSSSInstallDir(),
                out bool mpiInitialized
            );
<<<<<<< HEAD
            //Find dlls in own folder if called from ElectronBoSSSpad
            AppDomain.CurrentDomain.AssemblyResolve += CurrentDomain_AssemblyResolve;
        }

        string path = null;

        /// <summary>
        /// Resolve assembly not found exceptions. 
        /// In this case it happens when electron looks for dlls of BoSSSpad.exe in the folder of electron.exe
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="args"></param>
        /// <returns></returns>
        private System.Reflection.Assembly CurrentDomain_AssemblyResolve(object sender, ResolveEventArgs args)
        {
            // Ignore missing resources
            if (args.Name.Contains(".resources"))
                return null;

            // check for assemblies already loaded
            
            System.Reflection.Assembly assembly = AppDomain.CurrentDomain.GetAssemblies().
                FirstOrDefault(a => a.FullName == args.Name);
            if (assembly != null)
                return assembly;
            
            // Try to load by filename - split out the filename of the full assembly name
            // and append the base path of the original assembly (ie. look in the same dir)
            string filename = args.Name.Split(',')[0] + ".dll".ToLower();

            string asmFile = System.IO.Path.Combine(path, filename);

            try
            {
                return System.Reflection.Assembly.LoadFrom(asmFile);
            }
            catch (Exception ex)
            {
                return null;
=======

        }

        public static ElectronWorksheet Instance {
            get {
                return instance;
>>>>>>> 6fa3faeb05dbc80532ef484f3623d6be69e8da96
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
