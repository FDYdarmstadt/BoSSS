using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// in-situ post-processing for <see cref="CapillaryRise"/>
    /// </summary>
    [Serializable]
    public class CapillaryHeightLogging : CapillaryHeightLogging<XNSE_Control> { 
    }

    /// <summary>
    /// in-situ post-processing for <see cref="CapillaryRise"/>
    /// </summary>
    [Serializable]
    public class CapillaryHeightLogging<T> : XNSEinSituPostProcessingModule<T> where T:XNSE_Control, new() {
        
        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => "CapillaryHeight";

        /// <summary>
        /// CSV header
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {

            string header = String.Format("{0}\t{1}\t{2}\t{3}", "#timestep", "time", "capillary-height", "at-PositionX");
            Log.WriteLine(header);
                             
        }


        /// <summary>
        /// minimal y-coordinate of the interface
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {

            using(new FuncTrace()) {

                MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);

                double h_min = double.MaxValue, x_pos = 0.0;
                for(int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                    if(InterfacePoints[i, 1] < h_min) {
                        h_min = InterfacePoints[i, 1];
                        x_pos = InterfacePoints[i, 0];
                    }
                }

                string line = String.Format("{0}\t{1}\t{2}\t{3}", TimestepNo, phystime, h_min, x_pos);
                Log.WriteLine(line);
                Log.Flush();
            }
        }


    }
}
