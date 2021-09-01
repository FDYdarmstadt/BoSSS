using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases;

/// <summary>
/// Collection of overrides for the generic type in <see cref="BoSSS.Application.XNSE_Solver.XNSEinSituPostProcessingModule{T}"/>
/// </summary>
namespace BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// in-situ post-processing for <see cref="CapillaryRise"/>
    /// </summary>
    [Serializable]
    public class CapillaryHeightLogging : CapillaryHeightLogging<XNSFE_Control> {
    }

    /// <summary>
    /// Logging for Channel flow
    /// </summary>
    [Serializable]
    public class ChannelFlowLogging : ChannelFlowLogging<XNSFE_Control> { 
    }

    /// <summary>
    /// Benchmark quantities for droplet-testcases, <see cref="Droplet"/>
    /// </summary>
    [Serializable]
    public class Dropletlike : Dropletlike<XNSFE_Control> { 
    }

    /// <summary>
    /// Logging of stuff from Fourier Level-Set
    /// </summary>
    [Serializable]
    public class FourierLevelSetLogger : FourierLevelSetLogger<XNSFE_Control> { 
    }

    /// <summary>
    /// Errors against exact solution (<see cref="XNSE_Control.ExactSolutionVelocity"/>, <see cref="XNSE_Control.ExactSolutionPressure"/>)
    /// </summary>
    [Serializable]
    public class L2ErrorLogger : L2ErrorLogger<XNSFE_Control> { 
        // implement extension for temperature!!!
    }

    /// <summary>
    /// Used by <see cref="HeatedWall"/>, <see cref="TwoPhaseCouetteFlow"/>, etc.
    /// </summary>
    [Serializable]
    public class MovingContactLineLogging : MovingContactLineLogging<XNSFE_Control> { 
    }

    /// <summary>
    /// Post-processing specific to <see cref="RisingBubble"/>
    /// </summary>
    [Serializable]
    public class RisingBubble2DBenchmarkQuantities : RisingBubble2DBenchmarkQuantities<XNSFE_Control> { 
    }

    /// <summary>
    /// 
    /// </summary>
    public class WaveLikeLogging : WaveLikeLogging<XNSFE_Control> { 
    }

}
