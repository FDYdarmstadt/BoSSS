using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {


    /// <summary>
    /// Logging of evaporation, <see cref="HeatedWall"/>
    /// </summary>
    [Serializable]
    public class EvaporationLogging : XNSEinSituPostProcessingModule {
        
        /// <summary>
        /// Probably some specialization for the Fourier level set.
        /// </summary>
        public enum Mode {

            /// <summary>
            /// Line interface ?
            /// </summary>
            LineInterface = 1,

            /// <summary>
            /// circular interface ?
            /// </summary>
            CircleInterface = 2
        }


        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public Mode mode = Mode.LineInterface;


        /// <summary>
        /// Evaporation
        /// </summary>
        protected override string LogFileName => "Evaporation";

        /// <summary>
        /// CSV first line
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "'#timestep", "time", "interfacePosition", "meanInterfaceVelocity", "meanMassFlux"); //, "Temperatureprofile");
            Log.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {
            using(new FuncTrace()) {

                MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                double nNodes = InterfacePoints.GetLength(0);

                double posI = 0.0;
                if(this.mode == Mode.LineInterface)
                    posI = InterfacePoints.ExtractSubArrayShallow(-1, 1).To1DArray().Sum() / nNodes;

                double hVap = this.Control.ThermalParameters.hVap;
                double MassFlux = EvapVelocMean * hVap;

                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, posI, EvapVelocMean, MassFlux);

                //// temperature profile
                //int N = 100;
                //double[] tempP = new double[N + 1];

                //double[] probe = new double[N + 1];
                ////if (this.Control.LogValues == XNSE_Control.LoggingValues.EvaporationL) {
                //double L = this.Control.AdditionalParameters[0];
                //double x_probe = this.Control.AdditionalParameters[1];
                //for (int i = 0; i <= N; i++) {
                //    probe = new double[] { x_probe, i * (L / (double)N) };
                //    try {
                //        tempP[i] = this.Temperature.ProbeAt(probe);
                //    } catch {
                //        tempP[i] = 0.0;
                //    }
                //}
                ////}

                //line = line + "\t" + String.Join("\t", tempP.Select(ip => ip.ToString()).ToArray());

                Log.WriteLine(line);
                Log.Flush();


            }
        }
    }
}
