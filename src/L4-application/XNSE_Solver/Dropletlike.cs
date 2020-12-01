using BoSSS.Solution;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Benchmark quantities for droplet-testcases, <see cref="Droplet"/>
    /// </summary>
    [Serializable]
    class Dropletlike : InSituPostProcessingModule<XNSE_SolverMain, XNSE_Control> {
        protected override string LogFileName => "SemiAxis";


        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", "#timestep", "time", "semi axis x", "semi axis y", "area", "perimeter");
            textWriter.WriteLine(header);

        }

        public override void PerformPostProcessing(int iTimestep, double PhysTime, double dt) {
            throw new NotImplementedException();
        }
    }
}
