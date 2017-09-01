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

using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.LinSolvers;
using NSE_SIMPLE.BaseVariableDensity;
using System;

namespace NSE_SIMPLE.Multiphase {

    /// <summary>
    /// 
    /// </summary>
    public class MultiphaseSIMPLEControl : VariableDensitySIMPLEControl {

        /// <summary>
        /// 
        /// </summary>
        public new MaterialLawMultiphase EoS {
            get {
                MaterialLawMultiphase eos = base.EoS as MaterialLawMultiphase;
                return eos;
            }
            set {
                base.EoS = value;
            }
        }

        /// <summary>
        /// Under-relaxation factor level-set.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        [InclusiveUpperBound(1.0)]
        public double RelaxationFactorLevelSet;

        /// <summary>
        /// Modus for under-relaxation of level-set.
        /// </summary>
        public RelaxationTypes LevelSetRelaxationType;

        /// <summary>
        /// Convergence criterion level-set.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double L2NormLevelSetResidual;

        /// <summary>
        /// Linear solver configuration level-set.
        /// </summary>
        [NotNull]
        public Func<ISparseSolver> LevelSetSolverFactory;

        /// <summary>
        /// Analytic solution level-set.
        /// </summary>
        public Func<double[], double> AnalyticLevelSet = null;

        /// <summary>
        /// 
        /// </summary>
        public override void Verify() {
            base.Verify();

            if (base.EoS is MaterialLawMultiphase == false) {
                throw new Exception();
            }
        }
    }
}
