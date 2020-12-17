using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {


    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    class InterfacePointsLogger : XNSEinSituPostProcessingModule {
        
        /// <summary>
        /// 
        /// </summary>
        protected override string LogFileName => "InterfaceP";

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {

            string header = String.Format("{0}\t{1}\t{2}", "#timestep", "#time", "interfacePoints");
            Log.WriteLine(header);
            Log.Flush();

        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {
            using(new FuncTrace()) {
                double[] interfaceP;
                if(Fourier_LevSet != null) {
                    interfaceP = Fourier_LevSet.current_interfaceP.To1DArray();
                } else {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                    interfaceP = interP.ResizeShallow(interP.Length).To1DArray();
                }
                
                {
                    string logline = String.Format("{0}\t{1}", TimestepNo, phystime);
                    //for (int ip = 0; ip < interfaceP.Length; ip++) {
                    //    logline = logline + "\t" + interfaceP[ip].ToString();
                    //}
                    logline = logline + "\t" + String.Join("\t", interfaceP.Select(ip => ip.ToString()).ToArray());
                    Log.WriteLine(logline);
                    Log.Flush();
                }

                {
                    base.QueryResultTable.LogValue("time", phystime);
                    base.QueryResultTable.LogValue("interfaceP", interfaceP);
                }
            }
        }
    }
}
