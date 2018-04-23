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

    //Singleton class
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
            this.document = new Document();
        }

        public static ElectronWorksheet Instance{
            get{
                return instance; 
            }
        }

        public string RunCommand(string command){

            Document.Tuple singleCommandAndResult = new Document.Tuple{
                Command = command
            };

            singleCommandAndResult.Evaluate();
            return singleCommandAndResult.InterpreterTextOutput;
        }

        //save

        //load
    }
}
