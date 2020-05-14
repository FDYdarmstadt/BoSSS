using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.SipPoisson {


    /// <summary>
    /// All supported boundary condition types
    /// </summary>
    public enum BoundaryType {

        /// <summary>
        /// Dirichlet type boundary condition
        /// </summary>
        Dirichlet = 1,

        /// <summary>
        /// Neumann type boundary condition
        /// </summary>
        Neumann = 2
    }
}
