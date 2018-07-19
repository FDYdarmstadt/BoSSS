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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Tracing;
using ilPSP;

namespace BoSSS.Foundation {

    /// <summary>
    /// used by <see cref="DGField.ProjectFunction"/>
    /// </summary>
    /// <param name="U">
    /// value of some DG fields
    /// </param>
    /// <param name="X">
    /// spatial coordinate
    /// </param>
    /// <param name="jCell">
    /// cell index
    /// </param>
    /// <returns></returns>
    public delegate double Func(double[] X, double[] U, int jCell);

    partial class DGField {

        /// <summary>
        /// multiplies all DG coordinates of this field with a factor <paramref name="a"/>;
        /// </summary>
        /// <param name="a"></param>
        virtual public void Scale(double a) {
            Coordinates.Scale(a);
        }

        /// <summary>
        /// multiplies all DG coordinates of this field,
        /// within the cells in cell mask <paramref name="cm"/>, 
        /// by a factor <paramref name="a"/>;
        /// </summary>
        virtual public void Scale(double a, CellMask cm) {
            if (cm == null) {
                Scale(a);
            } else {

                BoSSS.Foundation.Basis b = this.Basis;
                var coord = this.Coordinates;

                foreach (Chunk c in cm) {
                    int JE = c.Len + c.i0;
                    for (int j = c.i0; j < JE; j++) {
                        int N = b.GetLength(j);
                        for (int n = 0; n < N; n++)
                            coord[j, n] *= a;
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="Acc(double,DGField,CellMask)"/>;
        /// </summary>
        /// <remarks>
        /// the accumulation is performed in all cells;
        /// </remarks>
        virtual public void Acc(double mult, DGField a) {
            Acc(mult, a, null);
        }

        /// <summary>
        /// accumulates an other field, <paramref name="a"/>*<paramref name="mult"/> to this field;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="mult"></param>
        /// <param name="cm">
        /// optional cell mask, null indicates the whole domain
        /// </param>
        /// <remarks>
        /// the basis of <paramref name="a"/> and of this field must be equal!
        /// </remarks>
        virtual public void Acc(double mult, DGField a, CellMask cm) {
            using (new FuncTrace()) {
                if (!a.Basis.Equals(this.Basis))
                    throw new ArgumentException("Basis of 'a' must be equal to basis of this field", "a");

                AccLaidBack(mult, a, cm);
            }
        }
        
        /// <summary>
        /// see <see cref="AccLaidBack(double,DGField,CellMask)"/>;
        /// </summary>
        public virtual void AccLaidBack(double mult, DGField a) {
            AccLaidBack(mult, a, null);
        }

        /// <summary>
        /// 'laxly' accumulation of another DG field to this field.
        /// </summary>
        /// <param name="mult"></param>
        /// <param name="a"></param>
        /// <param name="cm">
        /// optional cell mask, null indicates the whole domain
        /// </param>
        /// <remarks>
        /// in comparison to <see cref="Acc(double,DGField,CellMask)"/>, the DG
        /// basis of <paramref name="a"/> and this basis are not required to match exactly.
        /// </remarks>
        public virtual void AccLaidBack(double mult, DGField a, CellMask cm) {
            using (new FuncTrace()) {
                int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                IsCompatible(a);

                if (cm == null) {
                    for (int j = 0; j < J; j++) {
                        int N = Math.Min(this.Basis.GetLength(j), a.Basis.GetLength(j));
                        for (int n = 0; n < N; n++)
                            this.Coordinates[j, n] += a.Coordinates[j, n] * mult;
                    }
                } else {
                    foreach (Chunk c in cm) {
                        int IE = c.i0 + c.Len;
                        for (int j = c.i0; j < IE; j++) {
                            int N = Math.Min(this.Basis.GetLength(j), a.Basis.GetLength(j));
                            for (int n = 0; n < N; n++)
                                this.Coordinates[j, n] += a.Coordinates[j, n] * mult;
                        }
                    }
                }
            }
        }

        private void IsCompatible(DGField other) {
            if (!other.Basis.IsSubBasis(this.Basis) && !this.Basis.IsSubBasis(other.Basis))
                throw new ApplicationException("cant copy field values, because basis of other field isn't contained in basis of this field.");
            if (!object.ReferenceEquals(this.Basis.GridDat, other.Basis.GridDat))
                throw new ArgumentException("other field is different on different grid.", "other");

            int J = other.Coordinates.NoOfRows;

            if (J != this.Coordinates.NoOfRows)
                throw new ApplicationException("unknown application error.");
        }


        /// <summary>
        /// adds the constant <paramref name="a"/> to this field.
        /// </summary>
        /// <param name="a"></param>
        virtual public void AccConstant(double a) {
            AccConstant(a, null);
        }

        /// <summary>
        /// adds the constant <paramref name="a"/> to this field, but only within the cells in <paramref name="cm"/>.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="cm">
        /// optional cell-mask; if null, the operation is carried out in all cells
        /// </param>
        virtual public void AccConstant(double a, CellMask cm) {
            using (new FuncTrace()) {
                if (cm == null) {
                    int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                    for (int j = 0; j < J; j++)
                        SetMeanValue(j, GetMeanValue(j) + a);
                } else {
                    foreach (var c in cm) {
                        for (int j = 0; j < c.Len; j++)
                            SetMeanValue(j + c.i0, GetMeanValue(j + c.i0) + a);
                    }
                }
            }
        }

        /// <summary>
        /// adds two fields;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <remarks>
        /// Using this method will allocate a new Field;
        /// to avoid memory allocation, use <see cref="Acc(double,DGField)"/>;
        /// </remarks>
        public static DGField operator +(DGField a, DGField b) {
            DGField sum;

            if (a.Basis.IsSubBasis(b.Basis)) {
                sum = (DGField)b.Clone();
                sum.AccLaidBack(1.0, a);
            } else if (b.Basis.IsSubBasis(a.Basis)) {
                sum = (DGField)a.Clone();
                sum.AccLaidBack(1.0, b);
            } else {
                throw new ApplicationException("can't add the two fields, because their basis are incompatible");
            }

            if (a.Identification == null || a.Identification.Length <= 0
                || b.Identification == null || b.Identification.Length <= 0)
                sum.m_Identification = null;
            else
                sum.m_Identification = "(" + a.Identification + " + " + b.Identification + ")";

            return sum;
        }

        /// <summary>
        /// subtracts two field;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <remarks>
        /// Using this method will allocate a new Field;
        /// to avoid memory allocation, use <see cref="Acc(double,DGField)"/>;
        /// </remarks>
        public static DGField operator -(DGField a, DGField b) {
            DGField diff;

            if (a.Basis.IsSubBasis(b.Basis)) {
                diff = (DGField)b.Clone();
                diff.Scale(-1.0);
                diff.AccLaidBack(1.0, a);
            } else if (b.Basis.IsSubBasis(a.Basis)) {
                diff = (DGField)a.Clone();
                diff.AccLaidBack(-1.0, b);
            } else {
                throw new ApplicationException("can't subtract the two fields, because their basis are incompatible");
            }

            if (a.Identification == null || a.Identification.Length <= 0
                || b.Identification == null || b.Identification.Length <= 0)
                diff.m_Identification = null;
            else
                diff.m_Identification = "(" + a.Identification + " - " + b.Identification + ")";

            return diff;
        }

        /// <summary>
        /// multiplies a dg field with a constant
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        /// <remarks>
        /// Using this method will allocate a new Field;
        /// to avoid memory allocation, use <see cref="Scale(double)"/>;
        /// </remarks>
        public static DGField operator *(double alpha, DGField f) {
            DGField r = (DGField)f.Clone();
            r.Scale(alpha);
            return r;
        }

        /// <summary>
        /// multiplies a dg field with a constant
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        /// <remarks>
        /// Using this method will allocate a new Field;
        /// to avoid memory allocation, use <see cref="Scale(double)"/>;
        /// </remarks>
        public static DGField operator *(DGField f, double alpha) {
            return (alpha * f);
        }

        /// <summary>
        /// multiplies a dg field by 1 over <paramref name="alpha"/>
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        /// <remarks>
        /// Using this method will allocate a new Field;
        /// to avoid memory allocation, use <see cref="Scale(double)"/>;
        /// </remarks>
        public static DGField operator /(DGField f, double alpha) {
            return ((1.0 / alpha) * f);
        }

        /// <summary>
        /// see <see cref="ProjectAbs(double,CellMask, DGField[])"/>;
        /// </summary>
        virtual public void ProjectAbs<T>(double alpha, VectorField<T> vec) where T : DGField {
            ProjectAbs(alpha, null, vec.ToArray());
        }

        /// <summary>
        /// accumulates the projection of some vector field to this field, i.e.
        /// \f[ 
        ///   this = this + \alpha \cdot \| \vec{vec} \|.
        /// \f]
        /// </summary>
        /// <param name="alpha">factor \f$ \alpha \f$ </param>
        /// <param name="vec">vector field \f$ \vec{vec} \f$ </param>
        /// <param name="em">
        /// An optional restriction to the domain in which the projection is computed (it may, e.g.
        /// be only required in boundary cells, so a computation over the whole domain 
        /// would be a waste of computation power. If null, the computation is carried out in the whole domain;
        /// </param>
        virtual public void ProjectAbs(double alpha, CellMask em, params DGField[] vec) {
            int K = vec.Length;
            string[] args = new string[K];
            for (int k = 0; k < K; k++)
                args[k] = "_" + k;

            SpatialOperator powOp = new SpatialOperator(args, new string[] { "res" }, QuadOrderFunc.SumOfMaxDegrees());
            powOp.EquationComponents["res"].Add(new AbsSource(args));
            powOp.Commit();

            CoordinateVector coDom = new CoordinateVector(this);

            var ev = powOp.GetEvaluatorEx(
                new CoordinateMapping(vec),
                null,
                coDom.Mapping,
                edgeQrCtx: new EdgeQuadratureScheme(true, EdgeMask.GetEmptyMask(this.Basis.GridDat)),
                volQrCtx: new CellQuadratureScheme(true, em));

            ev.Evaluate<CoordinateVector>(alpha, 1.0, coDom); // only sources, no edge integrals required
        }

        /// <summary>
        /// a nonlinear source, which takes the absolute value of some DG fields
        /// used by <see cref="ProjectAbs(double,CellMask, DGField[])"/>;
        /// </summary>
        protected class AbsSource : INonlinearSource {

            /// <summary>
            /// ctor
            /// </summary>
            public AbsSource(string[] _args) {
                args = _args;
            }

            #region INonlinearSource Members

            void INonlinearSource.Source(double time, MultidimensionalArray x, MultidimensionalArray[] U, int IndexOffset, int j0, int Lenght, MultidimensionalArray Output) {
                int NoOfNodes = x.GetLength(1);
                int K = U.Length;

                for (int e = 0; e < Lenght; e++) {
                    for (int n = 0; n < NoOfNodes; n++) {

                        double acc = 0;
                        for (int k = 0; k < K; k++) {
                            double f = U[k][e + IndexOffset, n];
                            acc += f * f;
                        }
                        Output[e + IndexOffset, n] += Math.Sqrt(acc);
                    }
                }

            }

            #endregion

            #region IEquationComponent Members

            string[] args;

            IList<string> IEquationComponent.ArgumentOrdering {
                get {
                    return args;
                }
            }

            IList<string> IEquationComponent.ParameterOrdering {
                get {
                    return null;
                }
            }

            #endregion
        }

        /// <summary>
        /// Accumulates the DG-projection (with respect to the DG-basis
        /// of this field, <see cref="Basis"/>) of <br/>
        /// <paramref name="alpha"/>*<paramref name="f"/>^<paramref name="pow"/> to this field;
        /// </summary>
        /// <param name="alpha">scaling for accumulation</param>
        /// <param name="f">operand</param>
        /// <param name="pow">exponent</param>
        virtual public void ProjectPow(double alpha, DGField f, double pow) {
            ProjectPow(alpha, f, pow, null);
        }

        /// <summary>
        /// Accumulates the DG-projection (with respect to the DG-basis
        /// of this field, <see cref="Basis"/>) of 
        /// <paramref name="alpha"/>*<paramref name="f"/>^<paramref name="pow"/> to this field;
        /// </summary>
        /// <param name="alpha">scaling for accumulation</param>
        /// <param name="f">operand</param>
        /// <param name="pow">exponent</param>
        /// <param name="em">
        /// An optional restriction to the domain in which the projection is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computation
        /// power. A proper execution mask would be see e.g. 
        /// <see cref="GridData.BoundaryCells"/>.)
        /// if null, the computation is carried out in the whole domain;
        /// </param>
        virtual public void ProjectPow(double alpha, DGField f, double pow, CellMask em) {
            if (!object.ReferenceEquals(f.Basis.GridDat, this.Basis.GridDat))
                throw new ArgumentException("field is associated to another context.", "a");

            SpatialOperator powOp = new SpatialOperator(new string[] { "f" },
                                                        new string[] { "res" },
                                                        QuadOrderFunc.SumOfMaxDegrees());
            powOp.EquationComponents["res"].Add(new PowSource(pow));
            powOp.Commit();

            CoordinateVector coDom = this.CoordinateVector;

            var ev = powOp.GetEvaluatorEx(
                new CoordinateMapping(f), null, coDom.Mapping,
                edgeQrCtx: new EdgeQuadratureScheme(true, EdgeMask.GetEmptyMask(this.Basis.GridDat)),
                volQrCtx: new CellQuadratureScheme(true, em));

            ev.Evaluate<CoordinateVector>(alpha, 1.0, coDom); // only sources, no edge integrals required
        }

        /// <summary>
        /// a nonlinear source, which takes the power of some DG field, 
        /// used by <see cref="ProjectPow(double, DGField,double,CellMask)"/>;
        /// </summary>
        protected class PowSource : INonlinearSource {

            /// <summary>
            /// no text.
            /// </summary>
            /// <param name="pow">exponent</param>
            public PowSource(double pow) {
                m_Pow = pow;
            }

            double m_Pow = double.NaN;

            #region INonlinearSource Members

            /// <summary>
            /// multiplication of <paramref name="U"/>[0] and <paramref name="U"/>[1].
            /// </summary>
            public void Source(double time, MultidimensionalArray x, MultidimensionalArray[] U, int IndexOffset, int j0, int Lenght, MultidimensionalArray Output) {
                int NoOfNodes = x.GetLength(1);

                for (int e = 0; e < Lenght; e++) {
                    for (int n = 0; n < NoOfNodes; n++) {

                        double f = U[0][e + IndexOffset, n];
                        Output[e + IndexOffset, n] += Math.Pow(f, m_Pow);
                    }
                }
            }

            #endregion

            #region IEquationComponent Members

            /// <summary>
            /// always equal to { "a", "b" }, i.e. the variables that are multiplied.
            /// </summary>
            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "f" };
                }
            }

            #endregion

            /// <summary>
            /// no parameters, always null.
            /// </summary>
            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }

        /// <summary>
        /// used by <see cref="ProjectFunction"/>
        /// </summary>
        class ProjectFunctionSource : INonlinearSource {

            /// <summary>
            /// Domain variables
            /// </summary>
            private string[] Dom;

            /// <summary>
            /// The function to be applied.
            /// </summary>
            private Func m_f;

            /// <summary>
            /// Constructs a new source
            /// </summary>
            /// <param name="_Dom">
            /// Domain variables.
            /// </param>
            /// <param name="f">
            /// The function to be applied.
            /// </param>
            public ProjectFunctionSource(string[] _Dom, Func f) {
                Dom = _Dom;
                m_f = f;
            }

            /// <summary>
            /// The domain variables, see
            /// <see cref="ProjectFunctionSource.ProjectFunctionSource"/>
            /// </summary>
            public IList<string> ArgumentOrdering {
                get {
                    return Dom;
                }
            }

            /// <summary>
            /// Empty
            /// </summary>
            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }

            /// <summary>
            /// Applies the given in each evaluation point.
            /// </summary>
            /// <param name="time"></param>
            /// <param name="_x"></param>
            /// <param name="_U"></param>
            /// <param name="IndexOffset"></param>
            /// <param name="j0"></param>
            /// <param name="Lenght"></param>
            /// <param name="Output"></param>
            public void Source(double time,
                MultidimensionalArray _x,
                MultidimensionalArray[] _U,
                int IndexOffset, int j0, int Lenght,
                MultidimensionalArray Output) {

                int D = _x.GetLength(2);
                int J = _x.GetLength(0);
                int N = _x.GetLength(1);
                int Lambda = Dom.Length;

                Func f = m_f;
                double[] U = new double[Lambda];
                double[] X = new double[D];

                for (int j = 0; j < J; j++) {
                    for (int n = 0; n < N; n++) {
                        for (int delta = 0; delta < Lambda; delta++) {
                            U[delta] = _U[delta][j, n];
                        }

                        for (int d = 0; d < D; d++) {
                            X[d] = _x[j, n, d];
                        }

                        Output[j, n] = f(X, U, j + j0);
                    }
                }
            }
        }

        /// <summary>
        /// accumulates
        /// <paramref name="alpha"/>*<paramref name="f"/>(<paramref name="U"/>)
        /// to this field;
        /// </summary>
        /// <param name="alpha">scaling</param>
        /// <param name="f">some function</param>
        /// <param name="cqs">
        /// cell quadrature scheme: domain and quadrature rule
        /// </param>
        /// <param name="U">
        /// arguments for <paramref name="f"/>
        /// </param>
        public void ProjectFunction(double alpha, Func f, CellQuadratureScheme cqs, params DGField[] U) {

            string[] Dom = new string[U.Length];
            for (int i = 0; i < Dom.Length; i++)
                Dom[i] = "_" + i;

            string[] Cod = new string[] { "res" };

            SpatialOperator src = new SpatialOperator(Dom, Cod, QuadOrderFunc.NonLinear(3));
            src.EquationComponents[Cod[0]].Add(new ProjectFunctionSource(Dom, f));
            src.Commit();

            var ev = src.GetEvaluatorEx(
                new CoordinateMapping(U), null, this.Mapping,
                edgeQrCtx: new EdgeQuadratureScheme(false, EdgeMask.GetEmptyMask(this.Basis.GridDat)),
                volQrCtx: cqs);

            ev.Evaluate(alpha, 1.0, this.CoordinateVector);
        }

        /// <summary>
        /// projects the product of DG fields <paramref name="a"/>
        /// and <paramref name="b"/> onto this field;
        /// </summary>
        /// <param name="a">1st multiplicand</param>
        /// <param name="b">2nd multiplicand</param>
        /// <param name="alpha">
        /// scaling for <paramref name="a"/>*<paramref name="b"/>
        /// </param>
        virtual public void ProjectProduct(double alpha, DGField a, DGField b) {
            ProjectProduct(alpha, a, b, null);
        }

        /// <summary>
        /// Accumulates the DG-projection (with respect to the DG-basis
        /// of this field, <see cref="Basis"/>) of <br/>
        /// <paramref name="alpha"/>*<paramref name="a"/>*<paramref name="b"/> to this field;
        /// </summary>
        /// <param name="a">1st multiplicand</param>
        /// <param name="b">2nd multiplicand</param>
        /// <param name="alpha">scaling for <paramref name="a"/>*<paramref name="b"/></param>
        /// <param name="em">
        /// An optional restriction to the domain in which the projection is computed (it may, e.g.
        /// be only required in boundary cells, so a computation over the whole domain 
        /// would be a waste of computation power. A proper execution mask would be see e.g. 
        /// <see cref="GridData.BoundaryCells"/>.)<br/>
        /// if null, the computation is carried out in the whole domain;
        /// </param>
        public void ProjectProduct(double alpha, DGField a, DGField b, CellMask em) {
            ProjectProduct(alpha, a, b, em, true);
        }

        /// <summary>
        /// Calculates the DG-projection (with respect to the DG-basis
        /// of this field, <see cref="Basis"/>) of 
        /// <paramref name="alpha"/>*<paramref name="a"/>*<paramref name="b"/>
        /// and, depending on the value of <paramref name="accumulateResult"/>,
        /// either adds or assigns it to this field.
        /// </summary>
        /// <param name="a">1st multiplicand</param>
        /// <param name="b">2nd multiplicand</param>
        /// <param name="alpha">scaling for <paramref name="a"/>*<paramref name="b"/></param>
        /// <param name="em">
        /// An optional restriction to the domain in which the projection is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask would be see e.g. 
        /// <see cref="GridData.BoundaryCells"/>.)<br/>
        /// if null, the computation is carried out in the whole domain;
        /// </param>
        /// <param name="accumulateResult">
        /// Tells this method whether to accumulate (true) or not (false)
        /// </param>
        virtual public void ProjectProduct(double alpha, DGField a, DGField b, CellMask em, bool accumulateResult) {
            if (!object.ReferenceEquals(a.Basis.GridDat, this.Basis.GridDat))
                throw new ArgumentException("field is associated to another grid.", "a");
            if (!object.ReferenceEquals(b.Basis.GridDat, this.Basis.GridDat))
                throw new ArgumentException("field is associated to another grid.", "b");

            if (!accumulateResult) {
                if (a == this) {
                    a = (DGField)a.Clone();

                    if (b == this) {
                        b = a;
                    }
                } else if (b == this) {
                    b = (DGField)b.Clone();
                }

                this.Clear();
            }

            SpatialOperator multOp = new SpatialOperator(new string[] { "a", "b" },
                                                          new string[] { "res" },
                                                          QuadOrderFunc.NonLinear(2));
            multOp.EquationComponents["res"].Add(new MultiplySource());
            multOp.Commit();

            var ev = multOp.GetEvaluatorEx(
                new CoordinateMapping(a, b), null, this.Mapping,
                edgeQrCtx: new EdgeQuadratureScheme(true, EdgeMask.GetEmptyMask(this.Basis.GridDat)),
                volQrCtx: new CellQuadratureScheme(true, em));

            ev.Evaluate<CoordinateVector>(alpha, 1.0, this.CoordinateVector);
        }

        /// <summary>
        /// a nonlinear source, which multiplies two DG fields, used by
        /// <see cref="ProjectProduct(double, DGField, DGField, CellMask)"/>;
        /// </summary>
        protected class MultiplySource : INonlinearSource {

            #region INonlinearSource Members

            /// <summary>
            /// multiplication of <paramref name="U"/>[0] and <paramref name="U"/>[1].
            /// </summary>
            public void Source(double time, MultidimensionalArray x, MultidimensionalArray[] U, int IndexOffset, int j0, int Lenght, MultidimensionalArray Output) {
                int NoOfNodes = x.GetLength(1);

                for (int e = 0; e < Lenght; e++) {
                    for (int n = 0; n < NoOfNodes; n++) {

                        double a = U[0][e + IndexOffset, n];
                        double b = U[1][e + IndexOffset, n];

                        Output[e + IndexOffset, n] += a * b;
                    }
                }
            }

            #endregion

            #region IEquationComponent Members

            /// <summary>
            /// always equal to { "a", "b" }, i.e. the variables that are multiplied.
            /// </summary>
            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "a", "b" };
                }
            }

            #endregion

            /// <summary>
            /// no parameters, always null.
            /// </summary>
            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }

        /// <summary>
        /// projects the quotient of DG fields <paramref name="a"/>
        /// and <paramref name="b"/> times <paramref name="alpha"/> onto this field;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="alpha">scaling</param>
        public void ProjectQuotient(double alpha, DGField a, DGField b) {
            ProjectQuotient(alpha, a, b, null, true);
        }

        /// <summary>
        /// Calculates the DG-projection (with respect to the DG-basis
        /// of this field, <see cref="Basis"/>) of 
        /// <paramref name="alpha"/>*<paramref name="a"/>/<paramref name="b"/>
        /// and, depending on the value of <paramref name="accumulateResult"/>,
        /// either adds or assigns it to this field.
        /// </summary>
        /// <param name="a">1st multiplicand</param>
        /// <param name="b">2nd multiplicand</param>
        /// <param name="alpha">scaling for <paramref name="a"/>*<paramref name="b"/></param>
        /// <param name="accumulateResult">
        /// Tells this method whether to accumulate (true) or not (false)
        /// </param>
        /// <param name="cm">
        /// optional restriction to computational domain
        /// </param>
        virtual public void ProjectQuotient(
            double alpha, DGField a, DGField b, CellMask cm, bool accumulateResult) {

            if (!object.ReferenceEquals(a.Basis.GridDat, this.Basis.GridDat))
                throw new ArgumentException("field is associated to another grid.", "a");
            if (!object.ReferenceEquals(b.Basis.GridDat, this.Basis.GridDat))
                throw new ArgumentException("field is associated to another grid.", "b");

            if (!accumulateResult) {
                if (a == this) {
                    a = (DGField)a.Clone();

                    if (b == this) {
                        b = a;
                    }
                } else if (b == this) {
                    b = (DGField)b.Clone();
                }

                this.Clear();
            }

            SpatialOperator fracOp = new SpatialOperator(new string[] { "a", "b" },
                                                          new string[] { "res" },
                                                          QuadOrderFunc.Linear());
            fracOp.EquationComponents["res"].Add(new QuotientSource());
            fracOp.Commit();

            CoordinateVector coDom = this.CoordinateVector;

            var ev = fracOp.GetEvaluatorEx(
                new CoordinateMapping(a, b), null, coDom.Mapping,
                edgeQrCtx: new EdgeQuadratureScheme(true, EdgeMask.GetEmptyMask(this.Basis.GridDat)),
                volQrCtx: new CellQuadratureScheme(true, cm));

            ev.Evaluate<CoordinateVector>(alpha, 1.0, coDom);
        }

        /// <summary>
        /// a nonlinear source, which multiplies tow DG fields, used by
        /// <see cref="ProjectProduct(double, DGField, DGField, CellMask)"/>;
        /// </summary>
        protected class QuotientSource : INonlinearSource {

            #region INonlinearSource Members

            /// <summary>
            /// <paramref name="U"/>[0] / <paramref name="U"/>[1]
            /// </summary>
            public void Source(
                double time,
                MultidimensionalArray x,
                MultidimensionalArray[] U,
                int IndexOffset,
                int j0,
                int Lenght,
                MultidimensionalArray Output) {

                int NoOfNodes = x.GetLength(1);

                for (int e = 0; e < Lenght; e++) {
                    for (int n = 0; n < NoOfNodes; n++) {

                        double a = U[0][e + IndexOffset, n];
                        double b = U[1][e + IndexOffset, n];

                        Output[e + IndexOffset, n] += a / b;
                    }
                }
            }

            #endregion

            #region IEquationComponent Members

            /// <summary>
            /// always equal to { "a", "b" }, i.e. the variables form which the
            /// fraction a/b is computed.
            /// </summary>
            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "a", "b" };
                }
            }

            #endregion

            /// <summary>
            /// no parameters, always null.
            /// </summary>
            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }
    }
}
