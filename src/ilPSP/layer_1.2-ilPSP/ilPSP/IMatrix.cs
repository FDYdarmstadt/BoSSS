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

using System.Collections.Generic;

namespace ilPSP.Utils {

    /// <summary>
    /// minimal interface that a non-distributed (non-parallel) matrix must supply.
    /// </summary>
    public interface IMatrix {

        /// <summary>
        /// Number of Columns
        /// </summary>
        int NoOfCols {
            get;
        }

        /// <summary>
        /// Number of Rows
        /// </summary>
        int NoOfRows {
            get;
        }


        /// <summary>
        /// set/get an entry
        /// </summary>
        /// <param name="i">row index</param>
        /// <param name="j">column index</param>
        /// <returns></returns>
        double this[int i, int j] {
            get;
            set;
        }
        
        /// <summary>
        /// sets all entries to zero 
        /// </summary>
        void Clear();

        /// <summary>
        /// multiplies all entries of the matrix by <paramref name="a"/>
        /// </summary>
        /// <param name="a"></param>
        void Scale(double a);

        /// <summary>
        /// accumulates <paramref name="a"/>*<paramref name="M"/> to this matrix
        /// </summary>
        void Acc(double a, IMatrix M);

        /// <summary>
        /// Copies values from this array into a 2-dimensional .NET-array;
        /// </summary>
        /// <param name="dest">
        /// target for the copy operation;
        /// the [0,0]-entry of this matrix will end up in the [<paramref name="i0"/>,<paramref name="i1"/>]-entry
        /// of <paramref name="dest"/>;
        /// </param>
        /// <param name="i0">offset into this array, 1st index;</param>
        /// <param name="i1">offset into this array, 2nd index;</param>
        void CopyTo(double[,] dest, int i0, int i1);

        /// <summary>
        /// Copies to a one-dimensional array/list
        /// </summary>
        /// <param name="array"></param>
        /// <param name="RowWise">
        /// if true, elements are taken row-by-row (column index rotating fastest) and copied to <paramref name="array"/>;
        /// if false, elements are taken column-by-column;
        /// </param>
        /// <param name="arrayoffset"></param>
        void CopyTo<T>(T array, bool RowWise, int arrayoffset)
            where T : IList<double>;
    }

}
