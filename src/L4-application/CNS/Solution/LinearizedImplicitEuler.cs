using System;
using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Solution.Timestepping;
using ilPSP.LinSolvers;

namespace CNS.Solution {

    public class LinearizedImplicitEuler : ITimeStepper {

        private ImplicitEuler implicitEuler;

        private int updateInterval;

        private int updateCounter;

        public LinearizedImplicitEuler(SpatialOperator op, ISparseSolverExt solver, CNSFieldSet workingSet, int updateInterval) {
            if (updateInterval < 1) {
                throw new ArgumentException("Update interval must be positive", "updateInterval");
            }

            this.updateCounter = updateInterval;
        }

        #region ITimeStepper Members

        public double Time {
            get {
                return implicitEuler.Time;
            }
        }

        public void ResetTime(double NewTime) {
            implicitEuler.ResetTime(NewTime);
        }

        public void Perform(double dt) {
            if (updateCounter == 0) {
                implicitEuler = new ImplicitEuler();
                updateCounter += updateInterval;
            }

            implicitEuler.Perform(dt);
            updateCounter--;
        }

        public CoordinateMapping Mapping {
            get {
                return implicitEuler.Mapping;
            }
        }

        public CoordinateVector DGCoordinates {
            get {
                return implicitEuler.DGCoordinates;
            }
        }

        #endregion
    }
}
