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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using CNS.IBM;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CNS.ShockCapturing {

    /// <summary>
    /// See Zhang, JCP, 2017
    /// </summary>
    class ZhangPosivityPreservingLimiter : ILimiter {

        private double epsilon;

        private ImmersedSpeciesMap speciesMap;

        private IEnumerable<IChunkRulePair<QuadRule>> quadRuleSet;

        public IShockSensor Sensor {
            get {
                return null;
            }
        }

        public ZhangPosivityPreservingLimiter(double epsilon, ImmersedSpeciesMap speciesMap, IEnumerable<IChunkRulePair<QuadRule>> quadRuleSet) {
            this.epsilon = epsilon;
            this.speciesMap = speciesMap;
            this.quadRuleSet = quadRuleSet;
        }

        public void LimitFieldValues(IEnumerable<DGField> _fieldSet) {
            IProgram<CNSControl> program;
            CNSFieldSet fieldSet = program.WorkingSet;

            foreach (var chunkRulePair in quadRuleSet) {
                if (chunkRulePair.Chunk.Len > 1) {
                    throw new System.Exception();
                }

                MultidimensionalArray densityValues = MultidimensionalArray.Create(chunkRulePair.Chunk.Len, chunkRulePair.Rule.NoOfNodes);
                fieldSet.Density.Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, densityValues);
                MultidimensionalArray m0Values = MultidimensionalArray.Create(chunkRulePair.Chunk.Len, chunkRulePair.Rule.NoOfNodes);
                fieldSet.Momentum[0].Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, m0Values);
                MultidimensionalArray m1Values = MultidimensionalArray.Create(chunkRulePair.Chunk.Len, chunkRulePair.Rule.NoOfNodes);
                fieldSet.Momentum[1].Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, m1Values);
                MultidimensionalArray energyValues = MultidimensionalArray.Create(chunkRulePair.Chunk.Len, chunkRulePair.Rule.NoOfNodes);
                fieldSet.Energy.Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, energyValues);

                for (int i = 0; i < chunkRulePair.Chunk.Len; i++) {
                    int cell = i + chunkRulePair.Chunk.i0;

                    CellMask singleCellMask = new CellMask(speciesMap.Tracker.GridDat, chunkRulePair.Chunk);
                    CellQuadratureScheme singleCellScheme = new CellQuadratureScheme(
                        new FixedRuleFactory<QuadRule>(chunkRulePair.Rule), singleCellMask);
                    var singleCellRule = singleCellScheme.Compile(speciesMap.Tracker.GridDat, 99);

                    SpeciesId species = speciesMap.Tracker.GetSpeciesId(speciesMap.Control.FluidSpeciesName);
                    double volume = speciesMap.QuadSchemeHelper.NonAgglomeratedMetrics.CutCellVolumes[species][cell];

                    double integralDensity = fieldSet.Density.LxError((ScalarFunction)null, (X, a, b) => a, singleCellRule);
                    double meanDensity = integralDensity / volume;

                    if (meanDensity < epsilon) {
                        throw new System.Exception();
                    }

                    double minDensity = double.MaxValue;
                    for (int j = 0; j < chunkRulePair.Rule.NoOfNodes; j++) {
                        minDensity = Math.Min(minDensity, densityValues[i, j]);
                    }

                    double thetaDensity = Math.Min(1.0, (meanDensity - epsilon) / (meanDensity - minDensity));
                    if (thetaDensity < 1.0) {
                        Console.WriteLine("Scaled density in cell {0}!", cell);
                    }
                    
                    // Scale for positive density (Beware: Basis is not orthonormal on sub-volume!)
                    for (int j = 0; j < fieldSet.Density.Basis.Length; j++) {
                        fieldSet.Density.Coordinates[cell, j] *= thetaDensity;
                    }
                    fieldSet.Density.AccConstant(meanDensity * (1.0 - thetaDensity), singleCellMask);

                    // Re-evaluate since inner energy has changed
                    densityValues.Clear();
                    fieldSet.Density.Evaluate(cell, 1, chunkRulePair.Rule.Nodes, densityValues);

#if DEBUG
                    // Probe 1
                    for (int j = 0; j < chunkRulePair.Rule.NoOfNodes; j++) {
                        if (densityValues[i, j] - epsilon < -1e-15) {
                            throw new System.Exception();
                        }
                    }
#endif

                    double integralMomentumX = fieldSet.Momentum[0].LxError((ScalarFunction)null, (X, a, b) => a, singleCellRule);
                    double meanMomentumX = integralMomentumX / volume;
                    double integralMomentumY = fieldSet.Momentum[1].LxError((ScalarFunction)null, (X, a, b) => a, singleCellRule);
                    double meanMomentumY = integralMomentumY / volume;
                    double integralEnergy = fieldSet.Energy.LxError((ScalarFunction)null, (X, a, b) => a, singleCellRule);
                    double meanEnergy = integralEnergy / volume;

                    double meanInnerEnergy = meanEnergy - 0.5 * (meanMomentumX * meanMomentumX + meanMomentumY * meanMomentumY) / meanDensity;
                    if (meanInnerEnergy < epsilon) {
                        throw new System.Exception();
                    }

                    double minInnerEnergy = double.MaxValue;
                    for (int j = 0; j < chunkRulePair.Rule.NoOfNodes; j++) {
                        double innerEnergy = energyValues[i, j] - 0.5 * (m0Values[i, j] * m0Values[i, j] + m1Values[i, j] * m1Values[i, j]) / densityValues[i, j];
                        minInnerEnergy = Math.Min(minInnerEnergy, innerEnergy);
                    }

                    // Scale for positive inner energy (Beware: Basis is not orthonormal on sub-volume!)
                    double thetaInnerEnergy = Math.Min(1.0, (meanInnerEnergy - epsilon) / (meanInnerEnergy - minInnerEnergy));
                    if (thetaInnerEnergy < 1.0) {
                        Console.WriteLine("Scaled inner energy in cell {0}!", cell);
                    }

                    
                    for (int j = 0; j < fieldSet.Density.Basis.Length; j++) {
                        fieldSet.Density.Coordinates[chunkRulePair.Chunk.i0 + i, j] *= thetaInnerEnergy;
                        fieldSet.Momentum[0].Coordinates[chunkRulePair.Chunk.i0 + i, j] *= thetaInnerEnergy;
                        fieldSet.Momentum[1].Coordinates[chunkRulePair.Chunk.i0 + i, j] *= thetaInnerEnergy;
                        fieldSet.Energy.Coordinates[chunkRulePair.Chunk.i0 + i, j] *= thetaInnerEnergy;
                    }
                    fieldSet.Density.AccConstant(meanDensity * (1.0 - thetaInnerEnergy), singleCellMask);
                    fieldSet.Momentum[0].AccConstant(meanMomentumX * (1.0 - thetaInnerEnergy), singleCellMask);
                    fieldSet.Momentum[1].AccConstant(meanMomentumY * (1.0 - thetaInnerEnergy), singleCellMask);
                    fieldSet.Energy.AccConstant(meanEnergy * (1.0 - thetaInnerEnergy), singleCellMask);


#if DEBUG
                    // Probe 2
                    densityValues.Clear();
                    fieldSet.Density.Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, densityValues);
                    m0Values.Clear();
                    fieldSet.Momentum[0].Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, m0Values);
                    m1Values.Clear();
                    fieldSet.Momentum[1].Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, m1Values);
                    energyValues.Clear();
                    fieldSet.Energy.Evaluate(chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len, chunkRulePair.Rule.Nodes, energyValues);

                    for (int j = 0; j < chunkRulePair.Rule.NoOfNodes; j++) {
                        StateVector state = new StateVector(
                            speciesMap.GetMaterial(1.0),
                            densityValues[i, j],
                            new Vector(m0Values[i, j], m1Values[i, j], 0.0),
                            energyValues[i, j]);
                        if (!state.IsValid) {
                            throw new System.Exception();
                        }
                    }
#endif
                }
            }
        }
    }
}
