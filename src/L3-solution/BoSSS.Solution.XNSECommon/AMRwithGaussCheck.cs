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
using System.Text;
using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using BoSSS.Solution.LevelSetTools;
using System.Linq;
using ilPSP.Tracing;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// checks the Gauss theorem on cut cells and refines if a given error threshold is exceeded
    /// </summary>
    [Serializable]
    public class AMRwithGaussCheck : LevelSetTools.SolverWithLevelSetUpdater.AMRLevelIndicatorWithLevelset {

        public int levelSet = 0; // level set this Indicator should be active on

        public double errThreshhold = 1.0e-10;

        public override int[] DesiredCellChanges() {
            using (var tr = new FuncTrace("AMRcheckGauss")) {
                tr.InfoToConsole = true;

                int J = GridData.CellPartitioning.LocalLength;
                int[] levels = new int[J];

                var testField = GetTestField(GridData);

                CellMask CCmask = this.LsTrk.Regions.GetCutCellMask4LevSet(levelSet);

                int cellsToRefine = 0;
                int cellsToCoarse = 0;
                Cell[] cells = GridData.Grid.Cells;
                for (int j = 0; j < J; j++) {
                    int currentLevel = cells[j].RefinementLevel;
                    if (!CCmask.Contains(j))
                        continue;

                    var result = XNSEUtils.CheckGaussInCutCell(j, LsTrk, LsTrk.GetSpeciesId("A"), testField, LsTrk.GetCachedOrders().Max());
                    bool gaussViolated = result.error.Abs() > errThreshhold;
                    if (gaussViolated) {
                        tr.Info($"Gauss theorem error in cell {j} above threshhold {errThreshhold}: {result.error} ");
                    }
                    if (gaussViolated && currentLevel < maxRefinementLevel) {
                        levels[j] = 1;
                        cellsToRefine++;
                    } else if (!gaussViolated && currentLevel > 0) {
                        levels[j] = -1;
                        cellsToCoarse++;
                    }
                }

            return levels;
            }
            
        }


        public VectorField<SinglePhaseField> GetTestField(GridData grdDat) {

            Basis ScalarBasis = new Basis(grdDat, 0);

            SinglePhaseField testFieldX = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fx = (X => 1.0);
            testFieldX.ProjectField(fx);
            SinglePhaseField testFieldY = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fy = (X => 1.0);
            testFieldY.ProjectField(fy);
            SinglePhaseField testFieldZ = new SinglePhaseField(ScalarBasis);
            Func<double[], double> fz = (X => 1.0);
            testFieldZ.ProjectField(fz);

            if (grdDat.SpatialDimension == 2) {
                return new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY });
            } else {
                return new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY, testFieldZ });
            }

        }


        public override bool Equals(object obj) {
            if (!base.Equals(obj))
                return false;
            var other = obj as AMRwithGaussCheck;
            if (other == null)
                return false;
            if (other.levelSet != this.levelSet)
                return false;
            if (other.errThreshhold != this.errThreshhold)
                return false;
            return true;
        }


    }


}
