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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
using System.Linq;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// refinement of cells based on a gradient of a variable
    /// </summary>
    [Serializable]
    public class AMR_BasedOnFieldGradient : AMRLevelIndicatorWithLevelset {
        public string FieldName;

        public double Treshhold = 0.6;

        public AMR_BasedOnFieldGradient() {
        }

        public AMR_BasedOnFieldGradient(int _maxRefinementLevelval, double _refinementTreshhold, string _FieldName) {
            maxRefinementLevel = _maxRefinementLevelval;
            Treshhold = _refinementTreshhold;
            FieldName = _FieldName;
        }

        protected DGField FieldForRefinement {
            get {
                return ((XNSEC)this.SolverMain).CurrentStateVector.Fields.Where(f => f.Identification == FieldName).Single();
            }
        }

        public override int[] DesiredCellChanges() {
            //==========================================
            //Refinement based on gradient of a variable
            //==========================================
            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];
            int D = GridData.SpatialDimension;

            var field = this.FieldForRefinement;
            VectorField<XDGField> GradientOfField = new VectorField<XDGField>(D.ForLoop(d => new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Temperature].Degree), "Grad[" + d + "]")));
            GradientOfField.Gradient(1.0, field);

            var MagnitudeGradient = new SinglePhaseField(new Basis(base.GridData, 0), "Magnitude_Grad");
            MagnitudeGradient.Clear();
            MagnitudeGradient.ProjectFunction(1.0,
                (ilPSP.Vector X, double[] U, int jCell) => Math.Sqrt(U[0].Pow2() + U[1].Pow2()),
                new Foundation.Quadrature.CellQuadratureScheme(),
                GradientOfField.ToArray());

            double globalMinVal; double globalMaxVal;
            MagnitudeGradient.GetExtremalValues(out globalMinVal, out globalMaxVal);

            int cellsToRefine = 0;
            int cellsToCoarse = 0;
            Cell[] cells = GridData.Grid.Cells;
            for(int j = 0; j < J; j++) {
                double localMinVal, localMaxVal;
                MagnitudeGradient.GetExtremalValuesInCell(out localMinVal, out localMaxVal, j);

                int currentLevel = cells[j].RefinementLevel;
                if(localMaxVal / globalMaxVal > Treshhold && localMaxVal > 1e-10 && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                    cellsToRefine++;
                } 
                //else if(localMaxVal / globalMaxVal <= Treshhold  && currentLevel > 0) {
                //    levels[j] = -1;
                //    cellsToCoarse++;
                //}
            }
            Console.WriteLine("Refining {0} cells based on the gradient of {1}", cellsToRefine, FieldName);
            //Console.WriteLine("Coarsening {0} cells based on the gradient of {1}", cellsToCoarse, FieldName);

            return levels;
        }
    }
}