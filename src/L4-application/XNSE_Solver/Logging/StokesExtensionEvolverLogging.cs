using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Post-processing specific to <see cref="StokesExtensionEvolver"/>
    /// </summary>
    [Serializable]
    public class StokesExtensionEvolverLogging : StokesExtensionEvolverLogging<XNSE_Control> {
    }


    /// <summary>
    /// Post-processing specific to <see cref="StokesExtensionEvolver"/>
    /// </summary>
    [Serializable]
    public class StokesExtensionEvolverLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "StokesExtensionEvolver_Logging";

        /// <summary>
        /// Filename for logging ReInit values within the Stokes extension solver
        /// </summary>
        protected override string LogFileName => LogfileName;

        string LogFormat = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}";

        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format(LogFormat, "#timestep", "time", "L2-Norm_NB (rel)", "L2-Norm_CC (rel)", "L2minVal(cell,rel)", "L2maxVal(cell,rel)",  "MeanTotalValue_NB", "MeanTotalValue_CC", "InterfaceDist(rel)");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;

                AppendToLog(iTimestep);
                AppendToLog(PhysTime);
                LevelSetUpdater LSU = (LevelSetUpdater)SolverMainOverride.GetLevelSetUpdater();
                if(LSU.Updaters.TryGetValue("Phi", out LevelSetUpdater.SingleLevelSetUpdater evolver)) {
                    if(evolver.GetLsMover().GetType() == typeof(StokesExtensionEvolver)) {
                        double[] LogValues = ((StokesExtensionEvolver)evolver.GetLsMover()).GetLogValues();
                        if(LogValues != null) {
                            AppendToLog(LogValues);
                        }
                    } else {
                        tr.Warning("StokesExtensionSolverLogging is active but no StokesExtensionEvolver is used!");
                    }
                }
            }
        }
    }


}
