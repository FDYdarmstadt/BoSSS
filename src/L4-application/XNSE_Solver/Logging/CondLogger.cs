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

using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Logging {

    /// <summary>
    /// in-situ post-processing for <see cref="Solution.IApplication.OperatorAnalysis"/>
    /// </summary>
    [Serializable]
    public class CondLogger : CondLogger<XNSE_Control> {
        /// <summary>
        /// Logger for the condition number.
        /// </summary>
        /// <param name="config">Config object for OperatorAnalysis.  <see cref="BoSSS.Solution.OperatorAnalysisConfig"/></param>
        /// <param name="loggingStyle">Logging style, see <see cref="BoSSS.Application.XNSE_Solver.Logging.CondLogger{T}.logStyle"/>.</param>
        public CondLogger(OperatorAnalysisConfig config, int loggingStyle = 0) {
            SolverStage = 2;
            OperatorAnalysisConfiguration = config;
            logStyle = loggingStyle;
        }
    }

    /// <summary>
    /// For single droplet/bubble computations with origin in (0,0,0):
    /// extraction of the major axis lengths (theta = 0 and theta = 90, where theta describes the inclination angle) and volume  
    /// </summary>
    [Serializable]
    public class CondLogger<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => "CondNumbers";

        /// <summary>
        /// style of logfile, 0: writes variables in columns, 1: writes variables in rows.
        /// </summary>
        public int logStyle = 0;

        public OperatorAnalysisConfig OperatorAnalysisConfiguration { get => operatorAnalysisConfiguration; set => operatorAnalysisConfiguration = value; }

        /// <summary>
        /// CSV header
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            if (logStyle == 1) { 
            string header = String.Format("{0}\t{1}\t{2}", "#timestep", "varName", "value");
            Log.WriteLine(header);
            }
        }


        private OperatorAnalysisConfig operatorAnalysisConfiguration = new OperatorAnalysisConfig();

        private void WriteOutStuff(int TimestepNo, double PhysTime, IDictionary<string, double> dict) {
            if (logStyle == 0) {
                if (TimestepNo == 1) {
                    string header = "#timestep\tphysicalTime";
                    foreach (var kv in dict) {
                        header += $"\t{kv.Key}";
                    }
                    AppendLineToLog(header);
                }

                if (!String.IsNullOrEmpty(LogFileName)) {
                    string line = $"{TimestepNo}\t{PhysTime}";
                    foreach (var kv in dict) {
                        line += $"\t{kv.Value}";
                    }
                    AppendLineToLog(line);

                }
            } else if (logStyle == 1) {
                if (!String.IsNullOrEmpty(LogFileName)) {
                    foreach (var kv in dict) {
                        base.AppendLineToLog($"{TimestepNo}\t{kv.Key}\t{kv.Value}");           
                    }
                }
            } else { 
                throw new NotImplementedException("No such logging style in CondLogger");
            }

        }

        /// <summary>
        /// Driver routine for the application to call the post-processing
        /// </summary>
        override public void DriverTimestepPostProcessing(int iTimestep, double PhysTime) {
            if (iTimestep % LogPeriod == 0) {
                PerformTimestepPostProcessing(iTimestep, PhysTime);

                if (Log != null) {
                    Log.Flush();
                }
            }
        }

        internal void ComputeFullAnalysis(int TimestepNo, double PhysTime) {
            var dict = (this.SolverMainOverride as XNSE<XNSE_Control>).OperatorAnalysis(OperatorAnalysisConfiguration);
            WriteOutStuff(TimestepNo, PhysTime, dict);
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double physTime) {
            using (new FuncTrace()) {

                ComputeFullAnalysis(TimestepNo, physTime);

            }
        }
    }

}
