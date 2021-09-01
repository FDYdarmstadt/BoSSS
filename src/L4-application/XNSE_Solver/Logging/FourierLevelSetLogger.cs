using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Logging of stuff from Fourier Level-Set
    /// </summary>
    [Serializable]
    public class FourierLevelSetLogger : FourierLevelSetLogger<XNSE_Control> { }

    /// <summary>
    /// Logging of stuff from Fourier Level-Set
    /// </summary>
    [Serializable]
    public class FourierLevelSetLogger<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {
        protected override string LogFileName => throw new NotImplementedException();

        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            throw new NotImplementedException();
        }
    }
}
