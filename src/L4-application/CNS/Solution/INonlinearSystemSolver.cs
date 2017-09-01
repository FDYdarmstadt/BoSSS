using System;
using BoSSS.Foundation;

namespace CNS.Solution {

    interface INonlinearSystemSolver : IDisposable {

        int CurrentIteration {
            get;
        }

        CoordinateVector DGCoordinates {
            get;
        }

        CoordinateMapping CoordinateMapping {
            get;
        }

        SpatialOperator Operator {
            get;
        }

        void PerformIteration();

        void ReInitialize();
    }
}
