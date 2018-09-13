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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;

namespace CNS.IBM {

    /// <summary>
    /// Extensions of <see cref="CNSFieldSet"/> for IBM simulations
    /// </summary>
    public class IBMFieldSet : CNSFieldSet {

        private new IBMControl config;

        /// <summary>
        /// Optional level set field defining an immersed interface
        /// </summary>
        public readonly LevelSet LevelSet;

        /// <summary>
        /// Gradient of <see cref="LevelSet"/>
        /// </summary>
        public readonly DGField[] LevelSetGradient;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="config"></param>
        public IBMFieldSet(IGridData gridData, IBMControl config)
            : base(gridData, config) {
            this.config = config;

            if (config.RestartInfo == null && !config.FieldOptions.ContainsKey(IBMVariables.LevelSet)) {
                throw new Exception(
                    "Field 'levelSet' is required for IBM applications");
            }

            LevelSet = new LevelSet(
                new Basis(gridData, config.FieldOptions[IBMVariables.LevelSet].Degree),
                IBMVariables.LevelSet);

            LevelSetGradient = new DGField[CNSEnvironment.NumberOfDimensions];
            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                LevelSetGradient[d] = DerivedFields[IBMVariables.LevelSetGradient[d]];
            }
        }

        private IBMFieldSet(IGridData gridData, CNSControl config, IBMFieldSet template)
            : base(gridData, config, template) {
            LevelSet = template.LevelSet.CloneAs();
        }

        /// <summary>
        /// Returns the fields contained in <see cref="LevelSetGradient"/>
        /// since they are required for the calculation of normal vectors in
        /// interface cells
        /// </summary>
        public override DGField[] ParameterFields {
            get {
                return LevelSetGradient.Concat(base.ParameterFields).ToArray();
            }
        }

        /// <summary>
        /// Union of <see cref="CNSFieldSet.ConservativeVariables"/>,
        /// <see cref="CNSFieldSet.DerivedFields"/> and <see cref="LevelSet"/>
        /// </summary>
        public override IEnumerable<DGField> AllFields {
            get {
                return base.AllFields.Concat(LevelSet);
            }
        }

        /// <summary>
        /// Projects the initial level set after calling
        /// <see cref="CNSFieldSet.ProjectInitialValues"/>
        /// </summary>
        /// <param name="speciesMap"></param>
        /// <param name="initialValues"></param>
        public override void ProjectInitialValues(ISpeciesMap speciesMap, IDictionary<string, Func<double[], double>> initialValues) {
            LevelSet.ProjectField(NonVectorizedScalarFunction.Vectorize(
                X => config.LevelSetFunction(X, 0.0)));
            speciesMap.As<ImmersedSpeciesMap>().Tracker.UpdateTracker();

            base.ProjectInitialValues(speciesMap, initialValues);
        }

        /// <summary>
        /// Deep clone
        /// </summary>
        /// <returns></returns>
        public override object Clone() {
            return new IBMFieldSet(gridData, config, this);
        }
    }
}
