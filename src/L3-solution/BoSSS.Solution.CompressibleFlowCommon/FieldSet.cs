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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.CompressibleFlowCommon {
    public class FieldSet {

        /// <summary>
        /// The omnipresent context;
        /// </summary>
        protected IGridData gridData;

        /// <summary>
        /// Fields representing the derivatives of the primal variables
        /// </summary>
        /// <summary>
        /// The density $\rho$
        /// </summary>
        public DGField Density;

        /// <summary>
        /// The momentum field $\rho \vec{u}$
        /// </summary>
        public VectorField<DGField> Momentum;

        /// <summary>
        /// The total energy per volume $\rho E$
        /// </summary>
        public DGField Energy;

        /// <summary>
        /// Vector representation of <see cref="Density"/>,
        /// <see cref="Momentum"/> and <see cref="Energy"/>.
        /// </summary>
        public DGField[] ConservativeVariables {
            get {
                DGField[] fields = new DGField[CompressibleEnvironment.NumberOfDimensions + 2];

                fields[CompressibleEnvironment.PrimalArgumentToIndexMap[Variables.Density]] = Density;
                for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                    fields[CompressibleEnvironment.PrimalArgumentToIndexMap[Variables.Momentum[d]]] = Momentum[d];
                }
                fields[CompressibleEnvironment.PrimalArgumentToIndexMap[Variables.Energy]] = Energy;

                return fields;
            }
        }

        public FieldSet(IGridData gridData) {
            this.gridData = gridData;
        }

        public FieldSet(IGridData gridData, CompressibleControl config) {
            this.gridData = gridData;

            int numberOfDimensions = CompressibleEnvironment.NumberOfDimensions;

            SinglePhaseField[] momentumFields = new SinglePhaseField[numberOfDimensions];
            Basis momentumBasis = new Basis(gridData, config.MomentumDegree);

            // Mandatory fields
            Density = new SinglePhaseField(
                new Basis(gridData, config.DensityDegree),
                Variables.Density);

            for (int d = 0; d < numberOfDimensions; d++) {
                string variableName = Variables.Momentum[d];
                momentumFields[d] = new SinglePhaseField(momentumBasis, variableName);
            }
            Momentum = new VectorField<DGField>(momentumFields);

            Energy = new SinglePhaseField(
                new Basis(gridData, config.EnergyDegree), Variables.Energy);
        }

        //public abstract void UpdateDerivedVariables() {

        //}
    }
}
