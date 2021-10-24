using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.ZLSinSituPostProcessing
{
    /// <summary>
    /// Post-processing specific to <see cref="ConvergenceTest"/>
    /// </summary>
    [Serializable]
    public class BeamPositionLogging : ZLSinSituPostProcessingModule
    {

        public BeamPositionLogging()
        {

        }
        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "BenchmarkQuantities_Beam_ZLS";

        /// <summary>
        /// hard-coded name for the Beam_ZLS
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// Header for the Beam_ZLS log
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter)
        {
            string header = String.Format("{0}\t{1}\t{2}\t{3}", "#timestep", "Time", "X", "Y");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// compute and log
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime)
        {
            using (new FuncTrace())
            {
                double[] interfaceP;
                if (Fourier_LevSet != null)
                {
                    //interfaceP = Fourier_LevSet.current_interfaceP.To1DArray();
                    throw new DirectoryNotFoundException();
                }
                else
                {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet1);
                    int j = 0;
                    for (int i = 0; i < interP.Lengths[0]; i++)
                    {
                        if (interP[i, 1] > interP[j, 1])
                        {
                            j = i;
                        }
                    }
                    interfaceP = new double[2] { interP[j, 0], interP[j, 1] };
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
