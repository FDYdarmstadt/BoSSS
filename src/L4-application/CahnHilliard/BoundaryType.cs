using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.CahnHilliard {

    

    /// <summary>
    /// All supported boundary condition types
    /// </summary>
    public enum BoundaryType {

        /// <summary>
        /// Neumann type boundary condition for Walls
        /// </summary>
        Wall = 1,

        ///// <summary>
        ///// Flow Boundary Condition
        ///// </summary>
        Flow = 2

        ///// <summary>
        ///// Periodic Boundary Condition
        ///// </summary>
        //Periodic = 3
    }
    

}
