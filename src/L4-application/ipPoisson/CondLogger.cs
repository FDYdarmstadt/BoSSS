using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution;

namespace BoSSS.Application.SipPoisson {
    
    [Serializable]
    public class CondLogger : CondLogger<SipControl> { }

    /// <summary>
    /// Logs condition number of every timestep.
    /// Caution! Activating this will have a severe impact on performance!
    /// </summary>
    public class CondLogger<T> : InSituPostProcessingModule where T : SipControl, new() {

        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t", "#timestep", "type","value");
            Log.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// les gird
        /// </summary>
        protected BoSSS.Foundation.Grid.Classic.GridData GridData {
            get {
                return (Foundation.Grid.Classic.GridData)(this.SolverMainOverride.GridData);
            }
        }


        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            ComputeFullAnalysis(iTimestep);
        }

        protected override string LogFileName => "OpAnalysis";

        internal void ComputeFullAnalysis(int TimestepNo) {
            var dict = (this.SolverMainOverride as SipPoissonMain).OperatorAnalysis(new OperatorAnalysisConfig());
            WriteOutStuff(TimestepNo, dict);
        }

        private void WriteOutStuff(int TimestepNo, IDictionary<string,double> dict) {
            if (TimestepNo == 0) {
                foreach (var kv in dict) {
                    try {
                        base.QueryHandler.ValueQuery(kv.Key, kv.Value, false);
                    } catch (Exception ex) {
                        Console.WriteLine(ex.Message);
                    }
                }
            }
            if (!String.IsNullOrEmpty(LogFileName)) {
                foreach (var kv in dict) {
                    base.AppendToLog($"{TimestepNo}\t{kv.Key}\t{kv.Value}\t");
                }
            }
        }

        protected override char ColumnSeperator => '\n';

        /// <summary>
        /// reference to solver application class
        /// </summary>
        protected SipPoissonMain SolverMainOverride {
            get {
                return (SipPoissonMain)base.SolverMain;
            }
        }

        /// <summary>
        /// control object
        /// </summary>
        new protected T Control {
            get {
                return (T)(base.Control);
            }
        }
    }
}
