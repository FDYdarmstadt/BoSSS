﻿/* =======================================================================
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
using ilPSP.Utils;
using System;

namespace BoSSS.Solution.LoadBalancing {

    public static class CellCostEstimatorLibrary {

        public static ICellCostEstimator[] AllCellsAreEqual {
            get {
                int[] performanceClassToCostMap = new int[1];
                performanceClassToCostMap.SetAll(1);
                return new ICellCostEstimator[] {
                    new StaticCellCostEstimator(performanceClassToCostMap)
                };
            }
        }

        public static ICellCostEstimator[] MeasureCostOfExplicitOperatorEvaluation {
            get {
                return new ICellCostEstimator[] { 
                    new RuntimeCellCostEstimator(new string[][] {
                        new[] { "*RunSolverOneStep*", "Volume_Integration_NonLin" },
                        new[] { "*RunSolverOneStep*", "Edge_Integration_NonLin" },
                    })
                };
            }
        }

        public static ICellCostEstimator[] OperatorAssemblyAndCutCellQuadrules {
            get {
                return new ICellCostEstimator[] {
                    new RuntimeCellCostEstimator(new string[][] {
                        new[] { "*RunSolverOneStep*", "*LevelSetComboRuleFactory2.GetQuadRuleSet_Internal*" },
                        new[] { "*RunSolverOneStep*", "*Edge_Integration*" },
                        new[] { "*RunSolverOneStep*", "*Volume_Integration*" },
                        new[] { "*RunSolverOneStep*", "FIND_CUT_CELLS"},
                        new[] { "*RunSolverOneStep*", "*ComputeMassMatrixBlocks*"}
                    })
                };
            }
        }
    }
}
