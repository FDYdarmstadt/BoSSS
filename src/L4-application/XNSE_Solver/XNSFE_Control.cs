using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// control object for <see cref="XNSFE"/>
    /// </summary>
    public class XNSFE_Control : XNSE_Control {

        /// <summary>
        /// type for solver factory
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XNSFE);
        }

        /// <summary>
        /// How is the massflux coupled to the flow equations?
        /// </summary>
        public enum Coupling {
            weak = 0, // massflux as Parameter updated once per timestep
            iterative = 1, // massflux as Parameter updated iteratively per timestep
            strong = 2 // massflux evaluated from domain variable T, this makes even Stokes Problem nonlinear due to recoil pressure
        }

        /// <summary>
        /// Massflux Coupling option, standard is massflux as Parameter updated once per timestep
        /// </summary>
        public Coupling MassfluxCoupling = Coupling.weak;
    }
}
