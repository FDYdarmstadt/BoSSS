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
using BoSSS.Solution;

namespace BoSSS.Application.XNSE_Solver.LoadBalancing {

    /*
    public enum ClassifierType {
        VoidCutNormal = 1,
        CutCells = 2,
        Species = 3
    }
    
    public static class CellClassifier{ 
        public static (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IApplication<XNSE_Control> program, ClassifierType CType) {
            if (program.LsTrk == null)
                throw new ArgumentNullException("LsTrk not initialized! Not good, I need it!");

            switch (CType) {
                case ClassifierType.Species:
                    return SpeciesClassification(program);
                case ClassifierType.CutCells:
                    return CutCellClassification(program);
                case ClassifierType.VoidCutNormal:
                    return VoidCutNormalClassification(program);
                default:
                    throw new NotSupportedException("Type is not supported");
            }
        }

        private static (int noOfClasses, int[] cellToPerformanceClassMap) SpeciesClassification(IApplication<XNSE_Control> program) {
            var LsTrk = program.LsTrk;
            if (LsTrk == null)
                throw new NotSupportedException("Needs Information of Levelset tracker");

            int noOfClasses = LsTrk.TotalNoOfSpecies;
            int J = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellToPerformanceClassMap = new int[J];
            for (int iCell = 0; iCell < J; iCell++) {
                cellToPerformanceClassMap[iCell] = LsTrk.Regions.GetNoOfSpecies(iCell);
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }


        private static (int noOfClasses, int[] cellToPerformanceClassMap) VoidCutNormalClassification(IApplication<XNSE_Control> program) {
            var LsTrk = program.LsTrk;
            if (LsTrk == null)
                throw new NotSupportedException("Needs Information of Levelset tracker");

            int J = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellToPerformanceClassMap = new int[J];

            (int, int[]) tyield = CutCellClassification(program);
            int[] cutcellcostcluster = tyield.Item2;
            int noOfClasses = tyield.Item1 + 1;

            for (int iCell = 0; iCell < J; iCell++) {
                bool AtLeastOneSpecies = LsTrk.Regions.GetNoOfSpecies(iCell) >= 1;
                cellToPerformanceClassMap[iCell] = AtLeastOneSpecies ? 1 : 0;
                cellToPerformanceClassMap[iCell] += cutcellcostcluster[iCell];
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }

        /// <summary>
        /// 0: non cutcell
        /// 1: cutcell
        /// </summary>
        /// <param name="program"></param>
        /// <returns></returns>
        private static (int noOfClasses, int[] cellToPerformanceClassMap) CutCellClassification(IApplication<XNSE_Control> program) {
            var LsTrk = program.LsTrk;
            if (LsTrk == null)
                throw new NotSupportedException("Needs a Information of Levelset tracker");

            int noOfClasses = 2; // we distinguish only between cutcells and non-cutcells
            int J = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellToPerformanceClassMap = new int[J];
            foreach (int j in LsTrk.Regions.GetCutCellMask().ItemEnum) {
                cellToPerformanceClassMap[j] = 1;
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
    */
}
