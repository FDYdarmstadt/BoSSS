using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Logging of stuff from Fourier Level-Set
    /// </summary>
    [Serializable]
    public class FourierLevelSetLogger : XNSEinSituPostProcessingModule {
        protected override string LogFileName => throw new NotImplementedException();

        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            throw new NotImplementedException();
        }
    }
}
