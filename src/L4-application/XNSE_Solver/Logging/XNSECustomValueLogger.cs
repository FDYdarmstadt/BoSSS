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


using System;
using System.Collections.Generic;
using System.IO;


namespace BoSSS.Application.XNSE_Solver.Logging {

    /// <summary>
    /// Logger for custom logging values, <see cref="XNSECustomValueLogger"/>
    /// </summary>
    [Serializable]
    public class XNSECustomValueLogger : XNSECustomValueLogger<XNSE_Control> {
    }

    /// <summary>
    /// customizable module for adding various logging value over time
    /// </summary>
    [Serializable]
    public class XNSECustomValueLogger<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {


        protected override string LogFileName => "XNSE_LogValues";

        /// <summary>
        /// empty ctr
        /// </summary>
        public XNSECustomValueLogger() {
            this.LoggingValues = new List<XNSELogValue>();
        }


        protected List<XNSELogValue> LoggingValues;


        /// <summary>
        /// ctr defining custom LogValues
        /// </summary>
        //public XNSECustomValueLogger(List<XNSELogValue> LogValues) {
        //    this.LoggingValues = LogValues;
        //}


        public void AddLogValue(XNSELogValue newLogValue) {
            this.LoggingValues.Add(newLogValue);
        }


        protected override void WriteHeader(TextWriter textWriter) {
            string header = "#timesteps";
            foreach(var logVal in LoggingValues) {
                foreach(var valName in logVal.LogValues)
                    header += "\t" + valName;
            }
            Log.WriteLine(header);
            Log.Flush();
        }


        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            string line = "";
            foreach(var logVal in LoggingValues) {
                double[] retValues = logVal.ComputeLogValues(this);
                foreach(double rV in retValues) { 
                    line += $"\t{rV}";
                }
            }
            Log.Write(line);
            Log.Flush();
        }
    }

}
