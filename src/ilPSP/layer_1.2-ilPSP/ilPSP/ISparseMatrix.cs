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

using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Text;

namespace ilPSP.LinSolvers {
    
    /// <summary>
    /// defines the minimum requirements for a
    /// sparse matrix data structure that is used in the sparse solvers
    /// </summary>
    public interface ISparseMatrix {

        /// <summary>
        /// returns the diagonal element in the <paramref name="row"/>-th row.
        /// </summary>
        /// <param name="row">
        /// global row/column index; 
        /// must be in the index range which is defined by the row partition <see cref="RowPartitioning"/>,
        /// i.e. it must be greater or equal to <see cref="Partitioning.i0"/> 
        /// and smaller than (<see cref="Partitioning.LocalLength"/>+<see cref="Partitioning.i0"/>)!
        /// </param>
        /// <returns>value of diagonal element</returns>
        /// <remarks>
        /// this method may have better performance than e.g. <see cref="IMutableMatrix.GetValues"/>, ...
        /// </remarks>
        double GetDiagonalElement(int row);

        /// <summary>
        /// sets the diagonal element in the <paramref name="row"/>-th row
        /// to value <paramref name="val"/>
        /// </summary>
        /// <param name="row">global row/column index</param>
        /// <param name="val">new value of diagonal element</param>
        /// <remarks>
        /// this method may have better performance than e.g. <see cref="IMutableMatrix.SetValues"/>, ...
        /// </remarks>
        void SetDiagonalElement(int row, double val);

        /// <summary>
        /// row partition; this defines how the lines of the matrix are distributed over 
        /// the MPI processes;
        /// </summary>
        IPartitioning RowPartitioning { get; }


        /// <summary>
        /// column partition; this may not directly affect the storage scheme, 
        /// but the behavior of several operations, like transpose (where it defines the row partition of the result)
        /// or the Matrix/Vector product (where it defines the partition of the vector).
        /// </summary>
        IPartitioning ColPartition { get; }


        /// <summary>
        /// Number of rows over all MPI processes.
        /// </summary>
        int NoOfRows { get; }

        /// <summary>
        /// Number of columns over all MPI processes.
        /// </summary>
        int NoOfCols { get; }


        /// <summary>
        /// 'Sp'arse 'M'atrix/'V'ector 'M'ultiplication;<br/>
        /// Performs the calculation
        /// <paramref name="acc"/> = <paramref name="acc"/>*<paramref name="beta"/> + this*<paramref name="a"/>*<paramref name="alpha"/>;
        /// </summary>
        /// <typeparam name="VectorType1"></typeparam>
        /// <typeparam name="VectorType2"></typeparam>
        /// <param name="alpha"></param>
        /// <param name="a">
        /// length must be equal to the local length 
        /// (see <see cref="Partitioning.LocalLength"/>) of the column partition,
        /// <see cref="ColPartition"/>;
        /// </param>
        /// <param name="beta"></param>
        /// <param name="acc">
        /// length of accumulator must be at least equal to the local length (see <see cref="Partitioning.LocalLength"/>)
        /// of the row partition <see cref="RowPartitioning"/>
        /// </param>
        /// <remarks>
        /// The implementer is responsible for MPI-communication of Vector <paramref name="a"/>.
        /// </remarks>
        void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
            where VectorType1 : IList<double>
            where VectorType2 : IList<double>;

        /// <summary>
        /// MPI Communicator on which this object lives on
        /// </summary>
        MPI_Comm MPI_Comm {
            get;
        }
    }

    /// <summary>
    /// in contrast to the base interface, this matrix is also able to
    /// alter values;
    /// </summary>
    public interface IMutableMatrix : ISparseMatrix {

        /// <summary>
        /// gets multiple values from row <paramref name="RowIndex"/> at once
        /// </summary>
        /// <param name="RowIndex"></param>
        /// <param name="ColumnIndices"></param>
        /// <returns></returns>
        double[] GetValues(int RowIndex, int[] ColumnIndices);

        /// <summary>
        /// sets multiple values in one row at once
        /// </summary>
        /// <param name="RowIndex">
        /// index of the row that should be altered
        /// </param>
        /// <param name="ColumnIndices">
        /// the indices of the rows that should be set with new values;
        /// must have the same length as <paramref name="newValues"/>
        /// </param>
        /// <param name="newValues">
        /// new values for row number <paramref name="RowIndex"/>;
        /// must have the same length as <paramref name="ColumnIndices"/>
        /// </param>
        void SetValues(int RowIndex, int[] ColumnIndices, double[] newValues);

        /// <summary>
        /// gets/sets a specific entry of the matrix;
        /// </summary>
        /// <param name="i">global row index</param>
        /// <param name="j">global column index</param>
        /// <returns></returns>
        /// <remarks>
        /// for setting/getting the diagonal element, the methods <see cref="ISparseMatrix.GetDiagonalElement"/> and <see cref="ISparseMatrix.SetDiagonalElement"/>
        /// may have better performance
        /// </remarks>
        double this[int i, int j] {
            get;
            set;
        }
        
        /// <summary>
        /// True, if the occupation pattern may change due to a call to e.g. 
        /// <see cref="GetValues"/> or <see cref="SetValues"/>, i.e.
        /// that zero entries can be set to nonzero; However, it should be always possible
        /// to set occupied entries to 0.0, but this may not release the occupied memory.
        /// </summary>
        bool OccupationMutable {
            get;
        }

        /// <summary>
        /// Accumulates a block of entries to this matrix.
        /// </summary>
        /// <param name="i0">Row index offset.</param>
        /// <param name="j0">Column index offset.</param>
        /// <param name="alpha">Scaling factor for the accumulation.</param>
        /// <param name="Block">Block to add.</param>
        void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block);

        /// <summary>
        /// Accumulates a block of entries to this matrix.
        /// </summary>
        /// <param name="i0">Row index offset.</param>
        /// <param name="j0">Column index offset.</param>
        /// <param name="alpha">Scaling factor for the accumulation.</param>
        /// <param name="Block">Block to add.</param>
        /// <param name="beta">Scaling applied to this matrix before accumulation</param>
        void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block, double beta);
    }

    /// <summary>
    /// In contrast to the base interface, this matrix is also able to
    /// provide information about which entries are occupied. 
    /// One of the main advantages is that, with this info, sparse matrix-matrix - multiplication and other matrix-operations become possible,
    /// see <see cref="MsrMatrix.Multiply"/>.
    /// </summary>
    public interface IMutableMatrixEx : IMutableMatrix, ICloneable {

        /// <summary>
        /// Returns a collection of all occupied columns in a the row <paramref name="RowIndex"/>;
        /// </summary>
        /// <param name="RowIndex">Row index.</param>
        /// <param name="ColumnIndices">
        /// Ouput, the column indices of non-zero elements.
        /// The caller may provide a pre-allocated array in order to prevent the allocation of multiple integer arrays 
        /// on multiple calls to this method. 
        /// If the input array is to small, or null, a sufficiently large array is allocated.
        /// </param>
        /// <returns>
        /// Number of entries used in <paramref name="ColumnIndices"/>.
        /// </returns>
        int GetOccupiedColumnIndices(int RowIndex, ref int[] ColumnIndices);

        /// <summary>
        /// Returns a non-shallow copy of the row <paramref name="RowIndex"/>.
        /// </summary>
        /// <param name="RowIndex">Row index.</param>
        /// <param name="ColumnIndices">
        /// Ouput, the column indices of non-zero elements.
        /// The caller may provide a pre-allocated array in order to prevent the allocation of multiple integer arrays 
        /// on multiple calls to this method. 
        /// If the input array is to small, or null, a sufficiently large array is allocated.
        /// </param>
        /// <param name="Values">
        /// Ouput, the values for the respective columns in <paramref name="ColumnIndices"/>,
        /// allocation works in the same way as for <see cref="ColumnIndices"/>.
        /// </param>
        /// <returns>
        /// Number of entries used in <paramref name="ColumnIndices"/> and <paramref name="Values"/>.
        /// </returns>
        int GetRow(int RowIndex, ref int[] ColumnIndices, ref double[] Values);

        /// <summary>
        /// Sets all entries to 0.0; 
        /// </summary>
        void Clear();
    }

}
