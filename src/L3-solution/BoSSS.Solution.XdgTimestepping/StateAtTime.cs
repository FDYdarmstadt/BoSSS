using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.XdgTimestepping {



    /// <summary>
    /// Full state of solver at a specific point of time
    /// </summary>
    public class StateAtTime {

        public double[] SolutionState;
        public LevelSetTracker.TrackerBackup trkBkup;

        public double time {
            get {
                return trkBkup.time;
            }
        }


        /// <summary>
        /// Creates a full backup from the present state
        /// </summary>
        public static StateAtTime Obtain(XdgTimestepping mainObj, double time) {
            int HistIdx = mainObj.LsTrk.PopulatedHistoryIndices.ElementAtMin(idx => Math.Abs(time - mainObj.LsTrk.RegionsHistory[idx].Time));

            return new StateAtTime() {
                SolutionState = (new CoordinateVector(mainObj.CurrentState)).ToArray(),
                trkBkup = mainObj.LsTrk.BackupTimeLevel(HistIdx).CloneAs()
            };
        }

        /// <summary>
        /// Resets the time-integrator (<see cref="XdgTimestepping.TimesteppingBase"/>) to 
        /// start with this temporal state.
        /// </summary>
        /// <param name="mainObj"></param>
        public void Apply(XdgTimestepping mainObj) {
            mainObj.ResetTimestepper();
            mainObj.LsTrk.PopStacks();
            mainObj.LsTrk.ReplaceCurrentTimeLevel(trkBkup.CloneAs());
            //mainObj.LsTrk.UpdateTracker();
            mainObj.LsTrk.ObserverHack();

            (new CoordinateVector(mainObj.CurrentState)).SetV(SolutionState);

        }
    }
}
