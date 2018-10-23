using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// Codes for different boundary condition types
    /// </summary>
    public enum CompressibleBcType {
        /// <summary>
        /// %
        /// </summary>
        adiabaticSlipWall,
        
        /// <summary>
        /// %
        /// </summary>
        symmetryPlane,

        /// <summary>
        /// %
        /// </summary>
        adiabaticWall,

        /// <summary>
        /// %
        /// </summary>
        isothermalWall, 
        
        /// <summary>
        /// %
        /// </summary>
        subsonicInlet,
        
        /// <summary>
        /// %
        /// </summary>
        subsonicPressureInlet,
        
        /// <summary>
        /// %
        /// </summary>
        subsonicOutlet,
        
        /// <summary>
        /// %
        /// </summary>
        supersonicInlet,
        
        /// <summary>
        /// %
        /// </summary>
        supersonicOutlet
    }
}
