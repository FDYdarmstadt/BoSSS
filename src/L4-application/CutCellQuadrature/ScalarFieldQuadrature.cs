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
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using System;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CutCellQuadrature {

    class ScalarFieldQuadrature : CellQuadrature {

        public double Result;

        private DGField field;

        public ScalarFieldQuadrature(IGridData gridData, DGField field, CellQuadratureScheme quadInstr, int quadOrder)
            : base(new int[] { 1 }, gridData, quadInstr.Compile(gridData, quadOrder)) {
            this.field = field;
        }

     
        protected override void Evaluate(int i0, int Length, QuadRule qr, MultidimensionalArray EvalResult) {
            if (field is XDGField) {
                SpeciesId speciesB = new SpeciesId() {
                    cntnt = 11112
                };

                MultidimensionalArray result = EvalResult.ExtractSubArrayShallow(-1, -1, 0);
                ((XDGField)field).GetSpeciesShadowField(speciesB).Evaluate(i0, Length, qr.Nodes, result, 0, 0.0);
            } else {
                field.Evaluate(i0, Length, qr.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
            }
        }

        protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
            for (int i = 0; i < Length; i++) {
                Result += ResultsOfIntegration[i, 0];
            }
        }
    }
}
