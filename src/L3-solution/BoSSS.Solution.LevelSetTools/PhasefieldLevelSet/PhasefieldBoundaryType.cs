using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Solution.Control.AppControl;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// All supported boundary condition types
    /// </summary>
    public enum BoundaryType
    {

        /// <summary>
        /// on wall n*grad(phi) = 0, n*grad(mu) = 0; Neumann B.C.
        /// </summary>
        Wall = 0,

        /// <summary>
        /// at inlet the concentration is given, Dirichlet B.C.
        /// </summary>
        Inlet = 1,

        /// <summary>
        /// at outlet the inner values are used
        /// </summary>
        Outlet = 2,

        /// <summary>
        /// same as outlet
        /// </summary>
        Outflow = 3,

        /// <summary>
        /// Upwinding, could be inflow or outflow
        Pressure_Dirichlet = 4,

        /// <summary>
        /// like normal wall, all slip b.c. are treated in the same way
        Slip = 5,

        /// <summary>
        /// symmetry boundary condition with contact angle (90° contact angle) and free slip
        /// </summary>
        SlipSymmetry = 7,

    }

    partial class Phasefield
    {
        /// <summary>
        /// Translate the Boundary Values of the Base Control to the appropriate Cahn-Hilliard Boundary Values
        /// </summary>
        /// <param name="_BndIn"></param>
        /// <returns></returns>
        private IDictionary<string, BoundaryValueCollection> BoundaryTranslator(in IDictionary<string, BoundaryValueCollection> _BndIn)
        {
            // automatische generierung passender Cahn Hilliard Bedingungen, wie funktioniert das mit den EdgeTagNames und Boundaryfunctions etc.?
            IDictionary<string, BoundaryValueCollection> BndOut = new Dictionary<string, BoundaryValueCollection>();
            BndOut = _BndIn;
            return BndOut;
        }
    }

}
