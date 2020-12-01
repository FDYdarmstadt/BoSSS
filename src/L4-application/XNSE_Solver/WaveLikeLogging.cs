using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// 
    /// </summary>
    public class WaveLikeLogging : XNSEinSituPostProcessingModule {
        protected override string LogFileName => throw new NotImplementedException();

        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {
            using(new FuncTrace()) {

                Complex DFT_k;
                int numP;
                if(Fourier_LevSet != null) {
                    //amplitude = 2.0 * (Fourier_LevSet.DFT_coeff[1].Magnitude / Fourier_LevSet.current_samplP.Length);
                    DFT_k = Fourier_LevSet.DFT_coeff[(int)this.Control.AdditionalParameters[0]];
                    numP = Fourier_LevSet.current_samplP.Length;
                } else {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                    Complex[] DFT_coeff = DiscreteFourierTransformation.TransformForward_nonequidistant(interP, this.Control.AdditionalParameters[1]);
                    DFT_k = DFT_coeff[(int)this.Control.AdditionalParameters[0]];
                    numP = interP.Lengths[0];
                    //amplitude = -2.0 * DFT_coeff[1].Imaginary / (double)interP.Lengths[0];
                    //amplitude = DiscreteFourierTransformation.SingleSidedPowerSpectrum(DFT_coeff)[1];
                }
                string logline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, 2.0 * DFT_k.Magnitude / numP, 2.0 * DFT_k.Real / numP, -2.0 * DFT_k.Imaginary / numP);
                Log.WriteLine(logline);
                Log.Flush();

                return;
            }
        }
    }
}
