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


namespace CutCellQuadrature.TestCases {

    /// <summary>
    /// Generic base class for a shift of a given level set function.
    /// </summary>
    public abstract class Shift {

        /// <summary>
        /// Introduces a scaling of the magnitude of the shift. Useful to adapt
        /// the size of the shift to the grid resolution/mesh size
        /// </summary>
        /// <param name="scaling">
        /// A scaling factor
        /// </param>
        public abstract void Scale(double scaling);
    }

    /// <summary>
    /// Two-dimensional shift
    /// </summary>
    class Shift2D : Shift {

        /// <summary>
        /// Offset along the x-axis of the grid
        /// </summary>
        public double OffsetX {
            get;
            set;
        }

        /// <summary>
        /// Offset along the y-axis of the grid
        /// </summary>
        public double OffsetY {
            get;
            set;
        }

        /// <summary>
        /// <see cref="Shift.Scale"/>
        /// </summary>
        /// <param name="scaling">
        /// <see cref="Shift.Scale"/>
        /// </param>
        public override void Scale(double scaling) {
            OffsetX *= scaling;
            OffsetY *= scaling;
        }
    }

    /// <summary>
    /// Three-dimensional shift
    /// </summary>
    class Shift3D : Shift {

        /// <summary>
        /// Offset along the x-axis of the grid
        /// </summary>
        public double OffsetX {
            get;
            set;
        }

        /// <summary>
        /// Offset along the y-axis of the grid
        /// </summary>
        public double OffsetY {
            get;
            set;
        }

        /// <summary>
        /// Offset along the z-axis of the grid
        /// </summary>
        public double OffsetZ {
            get;
            set;
        }

        /// <summary>
        /// <see cref="Shift.Scale"/>
        /// </summary>
        /// <param name="scaling">
        /// <see cref="Shift.Scale"/>
        /// </param>
        public override void Scale(double scaling) {
            OffsetX *= scaling;
            OffsetY *= scaling;
            OffsetZ *= scaling;
        }
    }
}