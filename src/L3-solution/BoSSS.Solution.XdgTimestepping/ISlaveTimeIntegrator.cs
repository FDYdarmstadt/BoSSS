using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// Interface for a timestepper which 'follows' some Master-Timestepper.
    /// 
    /// </summary>
    public interface ISlaveTimeIntegrator {


        /// <summary>
        /// Routine to update the topmost time-slab.
        /// </summary>
        /// <param name="CurrentState">
        /// most recent solution state of the master timestepper.
        /// </param>
        /// <param name="time">
        /// Actual simulation time for the known value;
        /// </param>
        /// <param name="dt">
        /// Timestep size.
        /// </param>
        /// <param name="UnderRelax">
        /// </param>
        /// <returns>
        /// Some kind of residual in order to check convergence in a fully coupled simulation
        /// </returns>
        double Update(DGField[] CurrentState, double time, double dt, double UnderRelax, bool incremental);

        /// <summary>
        /// Accepts the current state and pushes a clone onto the stack
        /// </summary>
        void Push();


        /// <summary>
        /// discards the most recent solution and resets to the previous solution
        /// </summary>
        void Pop();




    }
}
