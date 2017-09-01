/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using System.Diagnostics;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using ilPSP.Tracing;
using BoSSS.Platform;
using MPI.Wrappers;

namespace NSE_SIMPLE {

    /// <summary>
    /// Utility functions for SIMPLESolvers, which might be used for incompressible flows
    /// as well as low Mach number flows.
    /// </summary>
    public static class SolverUtils {

        /// <summary>
        /// Mass defect related to divergence in continuity equation.
        /// </summary>
        /// <param name="Divergence"></param>
        /// <param name="Velocity"></param>
        /// <param name="MassDefect"></param>
        public static void CalculateMassDefect_Divergence(SIMPLEOperator[] Divergence, VectorField<SinglePhaseField> Velocity, SinglePhaseField MassDefect) {
            Debug.Assert(Divergence.Length == Velocity.Dim);

            for (int comp = 0; comp < Divergence.Length; comp++) {
                Divergence[comp].OperatorMatrix.SpMVpara(1.0, Velocity[comp].CoordinateVector, 1.0, MassDefect.CoordinateVector);
                MassDefect.CoordinateVector.Acc(1.0, Divergence[comp].OperatorAffine);
            }
        }

        
    }
}
