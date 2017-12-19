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
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Wrapper class for residual data and meta-data.
    /// </summary>
    public class ResidualLog : IDisposable {

        /// <summary>
        /// Residual values
        /// Key: variable name (the very first key is always the iteration number).
        /// Value: the residual values.
        /// </summary>                
        public Dictionary<string, IList<double>> Values {
            get;
            private set;
        }

        /// <summary>
        /// The identifier of the applied norm (e.g., "L2")
        /// </summary>
        public string Norm {
            get;
            set;
        }

        private StreamReader reader;
        private IList<int> columnIndices;
        private int stride;
        private string[] variables;
        public ISessionInfo session;
        private string CurrentLine = null;

        /// <summary>
        /// The id of the session the residuals belong to
        /// </summary>
        public Guid SessionGuid {
            get;
            private set;
        }

        /// <summary>
        /// Calls <see cref="Dispose"/>
        /// </summary>
        ~ResidualLog() {
            this.Dispose();
        }

        /// <summary>
        /// Reads the residual log and extracts the data (e.g. for plotting).
        /// </summary>
        /// <param name="session">The session in question.</param>
        /// <param name="norm">
        /// Identifier of the norm to be used. Use e.g. "L2_abs" if you want
        /// to read the residuals from file "residual-L2_abs.txt".
        /// </param>
        /// <param name="variables">
        /// The names of the variables from which the residuals are to be
        /// extracted. If empty, all variables will be used.
        /// This field will be overwritten with the actual names found in the
        /// file header.
        /// </param>
        /// <param name="stride">Determines in between how many iterations the parser skips.</param>
        public ResidualLog(ISessionInfo session, string norm, string[] variables, int stride) {
            this.stride = stride;
            this.SessionGuid = session.ID;
            this.variables = variables;
            this.session = session;
            // Finally, finish the object by setting all the properties.            
            this.Norm = norm;

            // Read the Data into the Logger
            try {
                ReadResiduals();
            }
            catch (Exception e) {
                Console.WriteLine(e);
                return;
            }
        }

        private void InitializeReaderAndValues() {
            // Get rid of bogus variables
            variables = variables.Where(v => !String.IsNullOrWhiteSpace(v)).ToArray();

            // Filename by convention
            string logPath = Path.Combine(
                DatabaseDriver.GetSessionDirectory(session),
                "residual-" + Norm + ".txt");

            reader = new StreamReader(new FileStream(
                    logPath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

            // Read the zeroth line ==> header
            string headerline = reader.ReadLine();
            string[] headerParts = headerline.Split(new char[] { '\t', ',' },
                StringSplitOptions.RemoveEmptyEntries);

            for (int i = 0; i < headerParts.Length; i++)
                headerParts[i] = headerParts[i].TrimStart('\t', ' ').TrimEnd('\t', ' '); // remove leading and trailing whitespaces

            // This specifies the column numbers we need to parse in each line.
            columnIndices = new List<int> { 0 }; // always extract the iteration number

            // Determine columnIndices:
            if (variables.Length == 0) {
                // No variables specified => use all columns
                columnIndices = Enumerable.Range(0, headerParts.Length).ToArray();
            } else {
                //columnIndices = new int[variables.Length + 1];  //+1 because of iteration numbers

                // Compare headerParts to variables and receive the column
                // indices we need to parse.
                for (int i = 0; i < variables.Length; i++) { // Loop over variables
                    string variable = variables[i];

                    // FindAll is useful if someone just wants to see all
                    // the residuals that contain "u", e.g. "u_1", "u_2", ...
                    int[] foundIndices = Array.FindAll(headerParts, col
                        => col.Contains(variable))
                        .Select(entry => Array.IndexOf(headerParts, entry))
                        .ToArray();
                    foreach (int idx in foundIndices) {
                        // index cannot be first column or nonexisting
                        if (idx > 0) {
                            columnIndices.Add(idx);
                        } else {
                            throw new ArgumentException("Invalid variable identifiers detected.");
                        }
                    }
                }
            }

            // Get rid of double entries
            columnIndices = columnIndices.Distinct().ToList();

            // Overwrite this so we see the proper names of the variables used,
            // excluding the iteration number
            variables = columnIndices.Where(idx => idx > 0)
                .Select(idx => headerParts[idx]).ToArray();

            // Initialize data fields, now that we know the number of variables
            Values = new Dictionary<string, IList<double>>();
            Values.Add(headerParts[0], new List<double>());
            for (int i = 1; i < columnIndices.Count; i++) {
                Values.Add(variables[i - 1], new List<double>());
            }
        }

        /// <summary>
        /// Read residuals from text file. Deletes previous Residuals.
        /// </summary>        
        public void ReadResiduals() {
            if (reader == null) {
                try {
                    InitializeReaderAndValues();
                } catch (FileNotFoundException e) {
                    throw new FileNotFoundException("The logfile has not been created yet.", e);
                }    
            }

            // delete previous residuals
            if (Values.ElementAt(0).Value.Count() > 1) {
                for (int col = 0; col < columnIndices.Count; col++) {
                    IList<double> CurrentCol = Values.ElementAt(col).Value;
                    double recycled = CurrentCol.Last();
                    CurrentCol.Clear();
                    CurrentCol.Add(recycled);
                }
            }

            // Read new residuals
            UpdateResiduals();
        }

        /// <summary>
        /// Adds freshly computed Residuals without deleting previous Residuals.
        /// </summary>
        public void UpdateResiduals() {
            if (reader == null) {
                try {
                    InitializeReaderAndValues();
                } catch (FileNotFoundException e) {
                    throw new FileNotFoundException("The logfile has not been created yet.",e);
                }
            }
            // read new residuals
            IList<string> NewLines = new List<string>();
            int NextChar;

            do {
                NextChar = reader.Read();
                if (NextChar >= 0) {
                    CurrentLine += (char)NextChar;
                    if (CurrentLine.Contains("\n")) {
                        NewLines.Add(CurrentLine);
                        CurrentLine = null;
                    }
                }
            } while (NextChar >= 0);

            // add new residuals to values list
            int lineCount = 0;

            foreach (string line in NewLines) {
                IList<string> lineParts = line.Split('\t', ',');
                if ((lineCount % stride == 0) || (lineCount == (NewLines.Count - 1))) {
                    for (int col = 0; col < columnIndices.Count; col++) {
                        // avoid confusion about which decimal separator to use.
                        Values.ElementAt(col).Value.Add(Double.Parse(
                            lineParts[columnIndices[col]], CultureInfo.InvariantCulture));
                    }
                }
                lineCount++;
            }
        }

        /// <summary>
        /// Closes the connection to the current file
        /// </summary>
        public void Dispose() {
            if (reader != null) {
                reader.Dispose();
            }
        }

        /// <summary>
        /// Plots the residuals into a graphic window. The Plot runs in a new Thread.
        /// </summary>
        public void Plot() {
            //Create a Plotter and start it in a new Thread.
            
            //try {
            //    ReadResiduals();
            //} catch (Exception e) {
            //    Console.WriteLine(e);
            //    return;
            //}

            Func<Form> Plotter = () => new ResidualForm(this);
            Thread threadPlotter = new Thread(() => AutonomuousPlotter.DisplayWindow(Plotter));
            threadPlotter.Start();
        }

        /// <summary>
        /// Live plotting of residuals during simulation. The LivePlotting tool runs in a new Thread.
        /// </summary>
        /// <param name="millisecondsPolling"></param>
        public void PlotLive(int millisecondsPolling = 500) {
            //Create a Plotter and run it in a new Thread.
            try {
                ReadResiduals();
            }catch(Exception e) {
                Console.WriteLine(e);
                return;
            }
 
            Func<Form> LivePlotter = () => new ResidualFormLive(this, millisecondsPolling);
            Thread threadLivePlotter = new Thread(() => AutonomuousPlotter.DisplayWindow(LivePlotter));
            threadLivePlotter.Start();
        }

        /// <summary>
        /// Displays the last line(s) of the residual log.
        /// </summary>
        /// <param name="lines">The number of lines to be printed, counting from end of file.</param>
        public void Tail(int lines = 1) {
            StringBuilder sb;
            try {
                sb = InitTail();
            }catch(Exception e) {
                Console.WriteLine(e);
                return;
            }

            for (int line = Values.ElementAt(0).Value.Count - lines; line < Values.ElementAt(0).Value.Count; line++) {
                // iteration number
                sb.Append(Values.ElementAt(0).Value[line] + "\t");
                for (int col = 1; col < Values.Count; col++) {
                    // residual value
                    sb.Append(Values.ElementAt(col).Value[line].ToString(
                        "E", NumberFormatInfo.InvariantInfo) + "\t");
                }
                sb.Append("\n");
            }

            Console.WriteLine(sb.ToString());
        }

        /// <summary>
        /// Show residuals live during simulation.
        /// </summary>
        public void TailLive() {
            StringBuilder sb;
            try {
                sb = InitTail();
            } catch (Exception e) {
                Console.WriteLine(e);
                return;
            }

            Console.WriteLine("Press any key to quit: ");
            while (Console.KeyAvailable == false) {
                Thread.Sleep(500);
                ReadResiduals();
                if (Values.ElementAt(0).Value.Count > 1) {
                    for (int line = 1; line < Values.ElementAt(0).Value.Count; line++) {
                        // iteration number
                        sb.Append(Values.ElementAt(0).Value[line] + "\t");
                        for (int col = 1; col < Values.Count; col++) {
                            // residual value
                            sb.Append(Values.ElementAt(col).Value[line].ToString(
                                "E", NumberFormatInfo.InvariantInfo) + "\t");
                        }
                        sb.Append("\n");
                    }
                    Console.Write(sb.ToString());
                    sb.Clear();
                }
            }
        }

        /// <summary>
        /// Init header for tail and read residuals once.
        /// </summary>
        /// <returns></returns>
        private StringBuilder InitTail() {
            StringBuilder sb = new StringBuilder();
            ReadResiduals();

            foreach (string variable in Values.Keys) {
                sb.Append("\t" + variable);
            }
            sb.Append("\n");

            return sb;
        }
    }
}
