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
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;

namespace CNS.IBM {

    /// <summary>
    /// Wrapper around <see cref="MassMatrixFactory"/> that ensures that the
    /// mass matrix is also invertible in void cells
    /// </summary>
    public class IBMMassMatrixFactory : IObserver<LevelSetTracker.LevelSetRegionsInfo> {

        /// <summary>
        /// 
        /// </summary>
        public readonly CoordinateMapping Mapping;

        /// <summary>
        /// 
        /// </summary>
        public readonly ImmersedSpeciesMap SpeciesMap;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="speciesMap"></param>
        /// <param name="mapping"></param>
        public IBMMassMatrixFactory(ImmersedSpeciesMap speciesMap, CoordinateMapping mapping) {
            this.Mapping = mapping;
            this.SpeciesMap = speciesMap;
            speciesMap.Tracker.Subscribe(this);
        }

        private MassMatrixFactory baseFactory;

        /// <summary>
        /// 
        /// </summary>
        public MassMatrixFactory BaseFactory {
            get {
                if (baseFactory == null) {
                    Basis maxBasis = Mapping.BasisS.ElementAtMax(b => b.Degree);
                    baseFactory = new MassMatrixFactory(
                        maxBasis,
                        SpeciesMap.QuadSchemeHelper.CellAgglomeration);
                }
                return baseFactory;
            }
        }

        /// <summary>
        /// Backing field for <see cref="MassMatrix"/>
        /// </summary>
        private BlockMsrMatrix massMatrix;

        /// <summary>
        /// 
        /// </summary>
        public BlockMsrMatrix MassMatrix {
            get {
                if (massMatrix == null) {
                    massMatrix = BaseFactory.GetMassMatrix(Mapping, false);
                    
                    // Make void part 0 instead of -1
                    CellMask fluidCells = SpeciesMap.SubGrid.VolumeMask;
                    foreach (Chunk chunk in fluidCells.Complement()) {
                        foreach (int cell in chunk.Elements) {
                            for (int fieldIndex = 0; fieldIndex < Mapping.BasisS.Count; fieldIndex++) {
                                for (int i = 0; i < Mapping.BasisS[fieldIndex].MaximalLength; i++) {
                                    //int localIndex = Mapping.LocalUniqueCoordinateIndex(fieldIndex, 0, i);
                                    //massMatrix[cell, localIndex, localIndex] = 0.0;
                                    int globalIndex = Mapping.GlobalUniqueCoordinateIndex(fieldIndex, cell, i);
                                    massMatrix[globalIndex, globalIndex] = 0.0;
                                }
                            }
                        }
                    }
                }

                return massMatrix;
            }
        }

        /// <summary>
        /// Backing field for <see cref="InverseMassMatrix"/>
        /// </summary>
        private BlockMsrMatrix inverseMassMatrix;

        /// <summary>
        /// 
        /// </summary>
        public BlockMsrMatrix InverseMassMatrix {
            get {
                if (inverseMassMatrix == null) {
                    inverseMassMatrix = MassMatrix.InvertBlocks(ignoreEmptyBlocks: true, Subblocks:true, SymmetricalInversion:true, OnlyDiagonal:true);
                }

                return inverseMassMatrix;
            }
        }

        BlockMsrMatrix nonAgglomeratedMassMatrix;
        
        public BlockMsrMatrix NonAgglomeratedMassMatrix {
            get {
                if (nonAgglomeratedMassMatrix == null) {
                    SpeciesId speciesId = SpeciesMap.Tracker.GetSpeciesId(SpeciesMap.Control.FluidSpeciesName);
                    nonAgglomeratedMassMatrix = BaseFactory.GetMassMatrix(
                        Mapping,
                        new Dictionary<SpeciesId, IEnumerable<double>>() {
                            { speciesId, Enumerable.Repeat(1.0, Mapping.NoOfVariables) } },
                        inverse: false,
                        VariableAgglomerationSwitch: new bool[Mapping.Fields.Count]);

                    // Make void part 0 instead of -1
                    CellMask fluidCells = SpeciesMap.SubGrid.VolumeMask;
                    foreach (Chunk chunk in fluidCells.Complement()) {
                        foreach (int cell in chunk.Elements) {
                            for (int fieldIndex = 0; fieldIndex < Mapping.BasisS.Count; fieldIndex++) {
                                for (int i = 0; i < Mapping.BasisS[fieldIndex].MaximalLength; i++) {
                                    //int localIndex = Mapping.LocalUniqueCoordinateIndex(fieldIndex, 0, i);
                                    //nonAgglomeratedMassMatrix[cell, localIndex, localIndex] = 0.0;
                                    int globalIndex = Mapping.GlobalUniqueCoordinateIndex(fieldIndex, cell, i);
                                    nonAgglomeratedMassMatrix[globalIndex, globalIndex] = 0.0;
                                }
                            }
                        }
                    }
                }

                return nonAgglomeratedMassMatrix;
            }
        }

        /// <summary>
        /// Discards all saved mass matrices
        /// </summary>
        /// <param name="value"></param>
        public void OnNext(LevelSetTracker.LevelSetRegionsInfo value) {
            this.baseFactory = null;
            this.nonAgglomeratedMassMatrix = null;
            this.massMatrix = null;
            this.inverseMassMatrix = null;
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="error"></param>
        public void OnError(System.Exception error) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void OnCompleted() {
        }
    }
}
