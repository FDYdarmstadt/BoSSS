using ilPSP;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Logging for Channel flow
    /// </summary>
    [Serializable]
    public class ChannelFlowLogging : ChannelFlowLogging<XNSE_Control> { }

    /// <summary>
    /// Logging for Channel flow
    /// </summary>
    [Serializable]
    public class ChannelFlowLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {
        
        /// <summary>
        /// 
        /// </summary>
        protected override string LogFileName => "ChannelFlow";

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#timestep", "time", "U_max", "derivedKinEnergy_max", "KinEnergy_max");
            Log.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {
            double C_halfw = 0.5;

            
            double Umax = base.CurrentVel[0].ProbeAt(new double[] { 1.01, C_halfw });
            double derivedKinEmax = 0.5 * Umax.Pow2();
            double KinEmax = this.KineticEnergy.ProbeAt(new double[] { 1.01, C_halfw });

            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, Umax, derivedKinEmax, KinEmax);

            Log.WriteLine(line);
            Log.Flush();
        }
    }
}
