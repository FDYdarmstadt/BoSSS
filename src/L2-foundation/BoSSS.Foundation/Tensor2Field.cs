/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System;

namespace BoSSS.Foundation {

    /// <summary>
    /// A 2<sup>nd</sup> stage tensor field (matrix field) as composition
    /// of scalar fields (<see cref="Field"/>);
    /// </summary>
    public class Tensor2Field<T> where T : Field {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ctx"></param>
        /// <param name="D"></param>
        /// <param name="b"></param>
        /// <param name="id"></param>
        public Tensor2Field(Context ctx, int D, Basis b, string id) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ctx"></param>
        /// <param name="D"></param>
        /// <param name="b"></param>
        public Tensor2Field(Context ctx, int D, Basis b)
            : this(ctx, D, b, "") { }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="components"></param>
        public Tensor2Field(T[,] components) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// dimension (number of vector components) of the vector field.
        /// </summary>
        public int Dim {
            get {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// access the compontent in the <paramref name="d"/>-th row and <paramref name="l"/>-th column
        /// </summary>
        /// <param name="d">row index</param>
        /// <param name="l">column index</param>
        /// <returns></returns>
        public T this[int d, int l] {
            get {
                throw new NotImplementedException();
            }
        }
        
        /// <summary>
        /// returns the <paramref name="d"/>-th column as some <see cref="VectorField{T}"/>-object
        /// </summary>
        /// <param name="d">row index</param>
        /// <returns></returns>
        public VectorField<T> GetRow(int d) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// returns the <paramref name="d"/>-th column as an array of <see cref="Field"/>-objects
        /// </summary>
        /// <param name="d">row index</param>
        /// <returns></returns>
        public T[] GetRowAsArray(int d) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// returns the <paramref name="l"/>-th row as some <see cref="VectorField{T}"/>-object
        /// </summary>
        /// <param name="l">column index</param>
        /// <returns></returns>
        public VectorField<T> GetColumn(int l) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// returns the <paramref name="l"/>-th as an array of <see cref="Field"/>-objects
        /// </summary>
        /// <param name="l">column index</param>
        /// <returns></returns>
        public T[] GetColumnAsArray(int l) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// returns the components of this stage-2 
        /// tensor field as some array of <see cref="Field"/>-objects
        /// </summary>
        /// <returns></returns>
        public T[,] AsArray() {
            throw new NotImplementedException();
        }


    }
}
