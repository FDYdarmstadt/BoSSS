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


namespace BoSSS.Foundation.Grid.RefElements {

    /// <summary>
    /// Encoding of elementary cell types.
    /// </summary>
    public enum CellType {

        /// <summary>
        /// Linear quad (3 points)
        /// </summary>
        Square_Linear,

        /// <summary>
        /// linear triangle
        /// </summary>
        Triangle_3,

        /// <summary>
        /// linear cube (4 points)
        /// </summary>
        Cube_Linear,

        /// <summary>
        /// linear (4-point) tetra
        /// </summary>
        Tetra_Linear,

        /// <summary>
        /// Quadratic tetra
        /// </summary>
        Tetra_10,

        /// <summary>
        /// Cubic tetra
        /// </summary>
        Tetra_20,

        /// <summary>
        /// Quartic tetra
        /// </summary>
        Tetra_35,

        /// <summary>
        /// Qunitic tetra
        /// </summary>
        Tetra_56,

        /// <summary>
        /// Bi-linear square
        /// </summary>
        Square_4,

        /// <summary>
        /// Quadratic square
        /// </summary>
        Square_8,

        /// <summary>
        /// Bi-quadratic square
        /// </summary>
        Square_9,

        /// <summary>
        /// Cubic square
        /// </summary>
        Square_12,

        /// <summary>
        /// Bi-cubic square
        /// </summary>
        Square_16,

        /// <summary>
        /// Bi-quartic square
        /// </summary>
        Square_25,

        /// <summary>
        /// Bi-quintic square
        /// </summary>
        Square_36,

        /// <summary>
        /// Bi-sextic square
        /// </summary>
        Square_49,

        /// <summary>
        /// Bi-septic square
        /// </summary>
        Square_64,

        /// <summary>
        /// Bi-octic square
        /// </summary>
        Square_81,

        /// <summary>
        /// Bi-nonic square (if such a thing exists...)
        /// </summary>
        Square_100,

        /// <summary>
        /// Cubic triangle
        /// </summary>
        Triangle_9,

        /// <summary>
        /// Quadratic triangle
        /// </summary>
        Triangle_6,

        /// <summary>
        /// ?
        /// </summary>
        Triangle_21,

        /// <summary>
        /// ?
        /// </summary>
        Triangle_15,

        /// <summary>
        /// What the hell?
        /// </summary>
        Triangle_152,

        /// <summary>
        /// ?
        /// </summary>
        Triangle_12,

        /// <summary>
        /// ?
        /// </summary>
        Triangle_10,

        /// <summary>
        /// Tri-linear cube
        /// </summary>
        Cube_8,

        /// <summary>
        /// Cubic cube (lol)
        /// </summary>
        Cube_20,

        /// <summary>
        /// Tri-cubic cube (lol again)
        /// </summary>
        Cube_27,

        /// <summary>
        /// Tri-quartic cube
        /// </summary>
        Cube_64,

        /// <summary>
        /// Tri-quintic cube
        /// </summary>
        Cube_125,

        /// <summary>
        /// Tri-sextic cube
        /// </summary>
        Cube_216,

        ///// <summary> </summary>
        //Cube_343,

        ///// <summary> </summary>
        //Cube_512,

        ///// <summary> </summary>
        //Cube_729,

        ///// <summary> </summary>
        //Cube_1000,

        /// <summary>
        /// Linear line element
        /// </summary>
        Line_2,

        /// <summary>
        /// Quadratic line element
        /// </summary>
        Line_3,

        /// <summary>
        /// Cubic line element
        /// </summary>
        Line_4,

        /// <summary>
        /// Quartic line element
        /// </summary>
        Line_5,

        /// <summary>
        /// Quintic line element
        /// </summary>
        Line_6,

        /// <summary>
        /// Zero-dimensional element
        /// </summary>
        Point
    }

    /// <summary>
    /// Extension methods
    /// </summary>
    static public class CellTypeExtensions {

        /// <summary>
        /// True, if <paramref name="type"/> represents an element with a
        /// linear transformation to physical space.
        /// </summary>
        static public bool IsLinear(this CellType type) {
            return (type == CellType.Square_Linear
                || type == CellType.Tetra_Linear
                || type == CellType.Cube_Linear
                || type == CellType.Triangle_3
                || type == CellType.Line_2
                || type == CellType.Point);
        }
    }
}
