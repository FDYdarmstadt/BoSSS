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
using BoSSS.Solution.CompressibleFlowCommon;
using ilPSP;
using System;
using System.Collections.Generic;

namespace CNS.Source {

    /// <summary>
    /// Implements a non-reflecting boundary condition by damping all
    /// deviations from a given reference state using a quadratic forcing
    /// function with user-definable strength. For details, see
    /// Mani, "Analysis and optimization of numerical sponge layers as a
    /// nonreflective boundary treatment", 2012
    /// </summary>
    public class SpongeLayerSource : INonlinearSource {

        private int variableIndex;

        private double referenceValue;

        private int dampingDirection;

        private double dampingStart;

        private double dampingEnd;

        private double dampingStrength;

        /// <summary>
        /// <see cref="ArgumentOrdering"/>
        /// </summary>
        private string[] argumentOrdering;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="variableIndex"></param>
        /// <param name="referenceValue"></param>
        /// <param name="dampingDirection"></param>
        /// <param name="dampingStart"></param>
        /// <param name="dampingEnd"></param>
        /// <param name="dampingStrength"></param>
        public SpongeLayerSource(int variableIndex, double referenceValue, int dampingDirection, double dampingStart, double dampingEnd, double dampingStrength) {
            this.variableIndex = variableIndex;
            this.referenceValue = referenceValue;
            this.dampingDirection = dampingDirection;

            this.dampingStart = dampingStart;
            this.dampingEnd = dampingEnd;
            this.dampingStrength = dampingStrength;

            this.argumentOrdering = CompressibleEnvironment.PrimalArgumentOrdering;
        }

        /// <summary>
        /// Applies a source term given by
        /// \f[ w \cdot s^2 \cdot (c - c_\text{ref})\f], where \f[ w \f] is a
        /// user-defined damping strength (see <see cref="dampingStrength"/>),
        /// \f[ s \f] is the local coordinate within the damping zone (0 at the
        /// start, 1 at the end), \f[ c \f] is the field variable to damped and
        /// \f[ c_\text{ref} \f] is the corresponding reference state.
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="U"></param>
        /// <param name="IndexOffset"></param>
        /// <param name="FirstCellInd"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
        public void Source(double time, MultidimensionalArray x, MultidimensionalArray[] U, int IndexOffset, int FirstCellInd, int Lenght, MultidimensionalArray Output) {
            int D = CompressibleEnvironment.NumberOfDimensions;
            int NoOfNodes = x.GetLength(1);
            double dampingWidth = dampingEnd - dampingStart;

            for (int i = 0; i < Lenght; i++) {
                for (int n = 0; n < NoOfNodes; n++) {
                    double actualValue = U[variableIndex][i + IndexOffset, n];
                    double xi = x[i + IndexOffset, n, dampingDirection];

                    double s = Math.Max((xi - dampingStart) / dampingWidth, 0.0);
                    Output[i + IndexOffset, n] += dampingStrength * s * s * (actualValue - referenceValue);
                }
            }
        }

        /// <summary>
        /// <see cref="CompressibleEnvironment.PrimalArgumentOrdering"/>
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return argumentOrdering;
            }
        }

        /// <summary>
        /// No parameters
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }
    }
}
