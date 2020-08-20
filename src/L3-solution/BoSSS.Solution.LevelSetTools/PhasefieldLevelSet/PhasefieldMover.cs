using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// Movement of the Phasefield
    /// </summary>
    partial class Phasefield : Application
    {
        /// <summary>
        /// Move Phasefield
        /// </summary>
        /// <param name="_TimestepNo"></param>
        /// <param name="_dt"></param>
        /// <param name="_phystime"></param>
        public void MovePhasefield(int _TimestepNo,double _dt, double _phystime)
        {
            RunSolverOneStep(_TimestepNo, _dt, _phystime);
        }

        /// <summary>
        /// Initial step to initialize level set near its equilibrium
        /// TODO make this to calculate initial steady state without convection
        /// </summary>
        /// <param name="_TimestepNo"></param>
        /// <param name="_dt"></param>
        /// <param name="_phystime"></param>
        public void RelaxationStep(int _TimestepNo = -1, double _dt = 0.01, double _phystime = 0.0)
        {
            this.Velocity.Clear();
            RunSolverOneStep(_TimestepNo, _dt, _phystime);
        }

        protected override double RunSolverOneStep(int _TimestepNo, double _dt, double _phystime)
        {
            using (new FuncTrace())
            {
                if (_TimestepNo == -1.0)
                {
                    Console.WriteLine("Initializing Phasefield");
                }
                else
                {
                    Console.WriteLine("Moving Phasefield in timestep #{0}, t = {1}, dt = {2} ...", _TimestepNo, _phystime, _dt);
                }

                // Perform timestep
                // ================

                this.m_Timestepper.m_ResLogger = new ResidualLogger(this.MPIRank, null, new Guid());
                this.m_Timestepper.m_ResidualNames = new string[] { "Res_phi", "Res_mu" };

                //PlotCurrentState(_phystime, new Foundation.IO.TimestepNumber(new int[] { _TimestepNo , 0}), 2);
                this.m_Timestepper.Solve(_phystime, _dt);
                PlotCurrentState(_phystime, new Foundation.IO.TimestepNumber(new int[] { _TimestepNo }), 2);

                // update DG LevelSet
                DGLevSet.Clear();
                DGLevSet.Acc(1.0, phi);

                // return
                // ======

                Console.WriteLine("done moving Phasefield in timestep #{0}.", _TimestepNo);

                return _dt;
            }
        }



    }
}
