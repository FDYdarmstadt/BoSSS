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
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Calculation of NusseltNumber.
    /// In fact, \int lambda * dT/dn ds is calculated.
    /// </summary>
    class NusseltNumber {

        GridData GridDat;

        SinglePhaseField Temperature;
        SinglePhaseField dTdx;
        SinglePhaseField dTdy;

        EdgeIntegral[] NusseltIntegrals;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="GridDat"></param>
        /// <param name="Temperature"></param>
        /// <param name="EoS"></param>
        /// <param name="edgeTagNames"></param>
        public NusseltNumber(GridData GridDat, SinglePhaseField Temperature, MaterialLaw EoS, string[] edgeTagNames) {
            this.GridDat = GridDat;
            this.Temperature = Temperature;

            //Basis BasisDerivative = new Basis(GridDat, Temperature.Basis.Degree - 1);
            Basis BasisDerivative = new Basis(GridDat, Temperature.Basis.Degree);
            dTdx = new SinglePhaseField(BasisDerivative);
            dTdy = new SinglePhaseField(BasisDerivative);

            NusseltIntegrals = new EdgeIntegral[edgeTagNames.Length];
            Nusselt = new double[edgeTagNames.Length];

            for (int bc = 0; bc < edgeTagNames.Length; bc++) {
                NusseltIntegrals[bc] = new EdgeIntegral(GridDat,
                    edgeTagNames[bc],
                    new NusseltFlux2D(EoS),
                    new CoordinateMapping(dTdx, dTdy, Temperature),
                    20);
            }
        }

        /// <summary>
        /// The calculated Nusselt number.
        /// In fact, \int lambda * dT/dn ds is calculated.
        /// </summary>
        public double[] Nusselt {
            get;
            private set;
        }

        /// <summary>
        /// Calculate current Nusselt number.
        /// </summary>
        public void CalculateNusseltNumber() {
            dTdx.Clear();
            //dTdx.Derivative(1.0, Temperature, 0, GridDat.BoundaryCells.VolumeMask);
            dTdx.DerivativeByFlux(1.0, Temperature, 0);
            dTdy.Clear();
            //dTdy.Derivative(1.0, Temperature, 1, GridDat.BoundaryCells.VolumeMask);
            dTdy.DerivativeByFlux(1.0, Temperature, 1);

            for (int bc = 0; bc < Nusselt.Length; bc++) {
                double LocalNusselt = NusseltIntegrals[bc].Evaluate();
                double GlobalNusselt = 0.0;

                unsafe {
                    MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&LocalNusselt), (IntPtr)(&GlobalNusselt), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.DOUBLE, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                }

                Nusselt[bc] = GlobalNusselt;
            }
        }
    }

    /// <summary>
    /// Nusselt 2D.
    /// </summary>
    class NusseltFlux2D : EdgeIntegral.EdgeFlux {

        MaterialLaw EoS;

        public NusseltFlux2D(MaterialLaw EoS) {
            this.EoS = EoS;
        }

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            double dTdx = Uin[0];
            double dTdy = Uin[1];
            double Temperature = Uin[2];
            double lambda = EoS.GetViscosity(Temperature);

            return lambda * (dTdx * normal[0] + dTdy * normal[1]);
        }

        public override IList<string> ArgumentOrdering {
            get { return new string[] { "dTdx", "dTdy", VariableNames.Temperature }; }
        }
    }
}
