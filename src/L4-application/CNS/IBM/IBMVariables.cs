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
using BoSSS.Platform;
using System;
using static CNS.Variable;

namespace CNS.IBM {

    /// <summary>
    /// Variables that are specific for IBM runs
    /// </summary>
    public class IBMVariables {

        /// <summary>
        /// The level set function defining the immersed boundary
        /// </summary>
        public static readonly Variable LevelSet = new Variable("levelSet", VariableTypes.Other);

        /// <summary>
        /// Gradient of <see cref="LevelSet"/> which defines the normal field
        /// in boundary cells
        /// </summary>
        public static readonly Vector<DerivedVariable> LevelSetGradient = new Vector<DerivedVariable>(
            d => new DerivedVariable(
                "levelSetGradient" + d,
                VariableTypes.Other,
                delegate (DGField gradientField, CellMask cellMask, IProgram<CNSControl> program) {
                    IBMControl control = program.Control as IBMControl;
                    if (control == null) {
                        throw new Exception(
                            "Level set gradient can only be computed in immersed boundary runs");
                    }

                    //// DON'T do that (at the moment). This will fail because
                    //// the update needs to be done at least once in the
                    //// start-up phase. We could ensure that, but the 
                    //// performance gain is marginal, so why bother?
                    //if (control.DomainType != DomainTypes.MovingImmersedBoundary) {
                    //    // No reason to compute anything
                    //    return;
                    //}

                    IBMFieldSet fieldSet = program.WorkingSet as IBMFieldSet;
                    CellMask cutCells = program.SpeciesMap.As<ImmersedSpeciesMap>().Tracker.Regions.GetCutCellMask();

                    gradientField.Clear();
                    gradientField.Derivative(
                        -1.0 * (double)control.FluidSpeciesSign,
                        fieldSet.LevelSet,
                        d,
                        cutCells);
                })
            );

        public static readonly DerivedVariable FluidCells = new DerivedVariable(
            "fluidCells",
            VariableTypes.Other,
            delegate (DGField fluidField, CellMask cellMask, IProgram<CNSControl> program) {
                fluidField.Clear();

                IBMControl control = program.Control as IBMControl;
                if (control == null) {
                    throw new Exception(
                        "Fluid cells can only be computed in immersed boundary runs");
                }

                foreach (Chunk chunk in cellMask) {
                    foreach (int cell in chunk.Elements) {
                        fluidField.SetMeanValue(cell, 1);
                    }
                }
            });

        public static readonly DerivedVariable CutCells = new DerivedVariable(
            "cutCells",
            VariableTypes.Other,
            delegate (DGField cutCellField, CellMask cellMask, IProgram<CNSControl> program) {
                cutCellField.Clear();

                IBMControl control = program.Control as IBMControl;
                if (control == null) {
                    throw new Exception(
                        "cutCells can only be computed in immersed boundary runs");
                }

                ImmersedSpeciesMap ibmSpeciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                CellMask cutCellMask = ibmSpeciesMap.Tracker.Regions.GetCutCellMask();

                foreach (Chunk chunk in cutCellMask) {
                    foreach (int cell in chunk.Elements) {
                        cutCellField.SetMeanValue(cell, 1);
                    }
                }
            });

        public static readonly DerivedVariable CutCellsWithoutSourceCells = new DerivedVariable(
            "cutCellsWithoutSourceCells",
            VariableTypes.Other,
            delegate (DGField cutCellWithoutSourceCellField, CellMask cellMask, IProgram<CNSControl> program) {
                cutCellWithoutSourceCellField.Clear();

                IBMControl control = program.Control as IBMControl;
                if (control == null) {
                    throw new Exception(
                        "cutCellsWithoutSourceCells can only be computed in immersed boundary runs");
                }

                ImmersedSpeciesMap ibmSpeciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                CellMask cutCellMask = ibmSpeciesMap.Tracker.Regions.GetCutCellMask().Except(ibmSpeciesMap.Agglomerator.AggInfo.SourceCells);

                foreach (Chunk chunk in cutCellMask) {
                    foreach (int cell in chunk.Elements) {
                        cutCellWithoutSourceCellField.SetMeanValue(cell, 1);
                    }
                }
            });

        public static readonly DerivedVariable FluidCellsWithoutSourceCells = new DerivedVariable(
            "fluidCellsWithoutSourceCells",
            VariableTypes.Other,
            delegate (DGField fluidCellWithoutSourceCellField, CellMask cellMask, IProgram<CNSControl> program) {
                fluidCellWithoutSourceCellField.Clear();

                IBMControl control = program.Control as IBMControl;
                if (control == null) {
                    throw new Exception(
                        "cutCellsWithoutSourceCells can only be computed in immersed boundary runs");
                }

                ImmersedSpeciesMap ibmSpeciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                CellMask cutCellMask = cellMask.Except(ibmSpeciesMap.Agglomerator.AggInfo.SourceCells);

                foreach (Chunk chunk in cutCellMask) {
                    foreach (int cell in chunk.Elements) {
                        fluidCellWithoutSourceCellField.SetMeanValue(cell, 1);
                    }
                }
    });

        public static readonly DerivedVariable SourceCells = new DerivedVariable(
            "sourceCells",
            VariableTypes.Other,
            delegate (DGField sourceCellField, CellMask cellMask, IProgram<CNSControl> program) {
                sourceCellField.Clear();

                IBMControl control = program.Control as IBMControl;
                if (control == null) {
                    throw new Exception(
                        "sourceCells can only be computed in immersed boundary runs");
                }

                ImmersedSpeciesMap ibmSpeciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                CellMask sourceCellMask = ibmSpeciesMap.Agglomerator.AggInfo.SourceCells;

                foreach (Chunk chunk in sourceCellMask) {
                    foreach (int cell in chunk.Elements) {
                        sourceCellField.SetMeanValue(cell, 1);
                    }
                }
            });

        //public static readonly DerivedVariable VoidCells = new DerivedVariable(
        //    "voidCells",
        //    VariableTypes.Other,
        //    delegate (DGField voidCells, CellMask cellMask, IProgram<CNSControl> program) {
        //        voidCells.Clear();

        //        IBMControl control = program.Control as IBMControl;
        //        if (control == null) {
        //            throw new Exception(
        //                "voidCells can only be computed in immersed boundary runs");
        //        }

        //        ImmersedSpeciesMap ibmSpeciesMap = program.SpeciesMap as ImmersedSpeciesMap;
        //        CellMask voidCellMask = .Except(cellMask); // TODO

        //        foreach (Chunk chunk in voidCellMask) {
        //            foreach (int cell in chunk.Elements) {
        //                voidCells.SetMeanValue(cell, 1);
        //            }
        //        }
        //    });
    }
}
