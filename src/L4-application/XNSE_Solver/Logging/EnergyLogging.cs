using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Log total Energy
    /// </summary>
    [Serializable]
    public class EnergyLogging : EnergyLogging<XNSE_Control> { }

    /// <summary>
    /// Benchmark quantities for droplet-testcases, <see cref="Droplet"/>
    /// </summary>
    [Serializable]
    public class EnergyLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// SemiAxis
        /// </summary>
        public const string LogfileName = "Energy";

        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => LogfileName;

        public EnergyLogging() {
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#timestep", "time", "total energy", "kinetic energy", "surface energy");
            textWriter.WriteLine(header);

        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {

            var result = ComputeEnergyProperties();

            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, result[0], result[1], result[2]);
            Log.WriteLine(line);
            Log.Flush();

            return;
        }

        double[] ComputeEnergyProperties() {

            double Etot, Ekin, Eint;

            try {
                Ekin = EnergyUtils.GetKineticEnergy(LsTrk, base.CurrentVel, new[] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B }, this.m_HMForder);
                Eint = EnergyUtils.GetSurfaceEnergy(LsTrk, Control.PhysicalParameters.Sigma, this.m_HMForder);
                Etot = Ekin + Eint;
            } catch {
                Ekin = -1;
                Etot = -1;
                Eint = -1;
            }

            return new double[] { Etot, Ekin, Eint };

        }
    }
}
