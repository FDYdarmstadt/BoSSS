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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation {

    public partial class DGField {

        /// <summary>
        /// computes the inner product of fields <paramref name="a"/> and <paramref name="b"/>
        /// </summary>
        static public double InnerProduct(DGField a, DGField b, CellQuadratureScheme quadScheme = null) {
            using (new FuncTrace()) {
                if (!Object.ReferenceEquals(a.GridDat, b.GridDat))
                    throw new ApplicationException("fields must be defined on the same grid.");
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                InnerProductQuadrature ipq = new InnerProductQuadrature(a, b, quadScheme, a.Basis.Degree + b.Basis.Degree);
                ipq.Execute();
                unsafe
                {
                    double innerProdTot = double.NaN;
                    double InnerProdLoc = ipq.m_InnerProd;
                    csMPI.Raw.Allreduce((IntPtr)(&InnerProdLoc), (IntPtr)(&innerProdTot), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                    return innerProdTot;
                }
            }
        }

        /// <summary>
        /// quadrature that calculates inner product of two fields, see <see cref="InnerProduct(DGField, DGField,CellQuadratureScheme)"/>
        /// </summary>
        class InnerProductQuadrature : BoSSS.Foundation.Quadrature.CellQuadrature {

            /// <summary>
            /// ctor
            /// </summary>
            /// <param name="a">1st field</param>
            /// <param name="b">2nd field</param>
            /// <param name="quaddegree">
            /// for the quadrature rule to chose, the 
            /// the degree of a polynomial that can be exactly integrated;
            /// </param>
            /// <param name="quadScheme"></param>
            public InnerProductQuadrature(DGField a, DGField b, CellQuadratureScheme quadScheme, int quaddegree)
                : base(new int[] { 1 }, a.GridDat, quadScheme.SaveCompile(a.GridDat, quaddegree)) {
                if (!object.ReferenceEquals(a.GridDat, b.GridDat))
                    throw new ArgumentException("Fields cannot be assigned to different grids.");

                m_FieldA = a;
                m_FieldB = b;
            }

            DGField m_FieldA;
            DGField m_FieldB;



            protected override void AllocateBuffers(int NoOfItems, NodeSet ruleNodes) {
                base.AllocateBuffers(NoOfItems, ruleNodes);
                int NoOfNodes = ruleNodes.GetLength(0);
                fieldvalsB.Allocate(NoOfItems, NoOfNodes);
            }


            MultidimensionalArray fieldvalsB = new MultidimensionalArray(2);

            protected override void Evaluate(int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) {
                NodeSet NodesUntransformed = rule.Nodes;
                int M = NodesUntransformed.GetLength(0);
                //int D = NodesUntransformed.GetLength(1);

                MultidimensionalArray fieldvalsA = EvalResult.ResizeShallow(new int[] { Length, NodesUntransformed.GetLength(0) });
                fieldvalsA.Clear();
                m_FieldA.Evaluate(i0, Length, NodesUntransformed, fieldvalsA, 1.0);
                fieldvalsB.Clear();
                m_FieldB.Evaluate(i0, Length, NodesUntransformed, fieldvalsB, 1.0);

                for (int j = 0; j < Length; j++)
                    for (int m = 0; m < M; m++) {
                        EvalResult[j, m, 0] = fieldvalsA[j, m] * fieldvalsB[j, m];
                    }
            }

            public double m_InnerProd = 0.0;


            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int j = 0; j < Length; j++) {
                    m_InnerProd += ResultsOfIntegration[j, 0];
                }
            }

        }

        /// <summary>
        /// see <see cref="GetExtremalValues(out double,out double,out int,out int,CellMask)"/>;
        /// </summary>
        public void GetExtremalValues(out double min, out double max) {
            int j, k;
            GetExtremalValues(out min, out max, out j, out k, null);
        }

        /// <summary>
        /// used by <see cref="GetExtremalValues(out double,out double,out int,out int,CellMask)"/> and <see cref="GetExtremalValuesInCell"/>
        /// to find the extremal values in 
        /// </summary>
        [NonSerialized]
        NodeSet[] m_ExtremalProbeNS;

        private void InitExtremalProbeNS() {
            if (m_ExtremalProbeNS == null) {
                // nodes for the evaluation of the velocity vector
                // (all vertices an 0)
                int D = Basis.GridDat.SpatialDimension;

                var KrefS = Basis.GridDat.iGeomCells.RefElements;
                m_ExtremalProbeNS = new NodeSet[KrefS.Length];
                for (int i = 0; i < KrefS.Length; i++) {
                    m_ExtremalProbeNS[i] = new NodeSet(KrefS[i], KrefS[i].Vertices);
                }
            }
        }


        /// <summary>
        /// Finds minimum and maximum value in cell <paramref name="jL"/>.
        /// </summary>
        /// <param name="min">On exit, the minimum value.</param>
        /// <param name="max">On exit, the maximum value.</param>
        /// <param name="jL">Local cell index.</param>
        public void GetExtremalValuesInCell(out double min, out double max, int jL) {
            // create node set
            InitExtremalProbeNS();

            double LocMax = -double.MaxValue, LocMin = double.MaxValue;

            foreach (int j in this.GridDat.GetGeometricCellIndices(jL)) {
                // les storage
                int iKref = Basis.GridDat.iGeomCells.GetRefElementIndex(j);
                int N = m_ExtremalProbeNS[iKref].NoOfNodes;
                MultidimensionalArray fieldValues = MultidimensionalArray.Create(1, N);

                // evaluate DG field
                this.Evaluate(j, 1, m_ExtremalProbeNS[iKref], fieldValues, 0.0);

                // loop over nodes ...
                for (int n = 0; n < N; n++) {
                    double vel = 0;
                    vel = fieldValues[0, n];
                    if (vel > LocMax) {
                        LocMax = vel;
                    }
                    if (vel < LocMin) {
                        LocMin = vel;

                    }
                }
            }

            // return
            min = LocMin;
            max = LocMax;
        }

        /// <summary>
        /// Returns the minimum and the maximum value of the DG field
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        /// <param name="max">
        /// on all invoking MPI processes, the maximum value (over all processors) of
        /// this field;
        /// </param>
        /// <param name="min">
        /// on all invoking MPI processes, the minimum value (over all processors) of
        /// this field;
        /// </param>
        /// <remarks>
        /// to find the maximum value, this field is evaluated on each vertex and the center of the simplex.
        /// </remarks>
        /// <param name="cm">
        /// optional domain restriction
        /// </param>
        /// <param name="jMaxGlob">
        /// global cell index in which the maximum value <paramref name="max"/> was found
        /// </param>
        /// <param name="jMinGlob">
        /// global cell index in which the maximum value <paramref name="min"/> was found
        /// </param>
        public void GetExtremalValues(out double min, out double max, out int jMinGlob, out int jMaxGlob, CellMask cm = null) {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            using (new FuncTrace()) {
                int J = Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                // find maximum on this processor
                // ------------------------------

                // create node set
                InitExtremalProbeNS();

                // vectorisation of this method
                int VectorSize = -1;
                int N0 = m_ExtremalProbeNS[0].GetLength(0);   // number of nodes
                MultidimensionalArray fieldValues = new MultidimensionalArray(2);


                IEnumerable<Chunk> all_chunks;
                if (cm == null) {
                    var _ch = new Chunk[1];
                    _ch[0].i0 = 0;
                    _ch[0].Len = J;
                    all_chunks = _ch;
                } else {
                    all_chunks = cm;
                }

                double LocMax = -double.MaxValue, LocMin = double.MaxValue;
                int jMinLoc = int.MaxValue, jMaxLoc = int.MaxValue;
                foreach (Chunk chk in all_chunks) {
                    VectorSize = this.GridDat.iGeomCells.GetNoOfSimilarConsecutiveCells(CellInfo.RefElementIndex_Mask, chk.i0, Math.Min(100, chk.Len)); // try to be a little vectorized;

                    int _J = chk.Len + chk.i0;
                    for (int j = chk.i0; j < _J; j += VectorSize) {
                        if (j + VectorSize > _J)
                            VectorSize = _J - j;
                        int iKref = Basis.GridDat.iGeomCells.GetRefElementIndex(j);
                        int N = m_ExtremalProbeNS[iKref].GetLength(0);

                        if (fieldValues.GetLength(0) != VectorSize || fieldValues.GetLength(1) != N)
                            fieldValues = MultidimensionalArray.Create(VectorSize, N);

                        this.Evaluate(j, VectorSize, m_ExtremalProbeNS[iKref], fieldValues, 0.0);

                        // loop over cells ...
                        for (int jj = j; jj < j + VectorSize; jj++) {

                            // loop over nodes ...
                            for (int n = 0; n < N; n++) {
                                double vel = 0;
                                vel = fieldValues[jj - j, n];
                                if (vel > LocMax) {
                                    LocMax = vel;
                                    jMaxLoc = jj;
                                }
                                if (vel < LocMin) {
                                    LocMin = vel;
                                    jMinLoc = jj;
                                }
                            }
                        }
                    }
                }


                // find the maximum over all processes via MPI and return
                // ------------------------------------------------------
                if (Basis.GridDat.CellPartitioning.MpiSize > 1) {
                    double[] total = new double[2];
                    double[] local = new double[] { LocMax, -LocMin };
                    unsafe
                    {
                        fixed (double* ploc = local, ptot = total) {
                            csMPI.Raw.Allreduce((IntPtr)ploc, (IntPtr)ptot, 2, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                        }
                    }

                    max = total[0];
                    min = -total[1];

                    int i0 = Basis.GridDat.CellPartitioning.i0;
                    int[] jGlob = new int[2];
                    if (max == LocMax) {
                        // (at least one) global maximum on this processor
                        jGlob[0] = jMaxLoc + i0;
                    } else {
                        // maximum on some other processor
                        jGlob[0] = int.MaxValue;
                    }

                    if (min == LocMin) {
                        // (at least one) global minimum on this processor
                        jGlob[1] = jMinLoc + i0;
                    } else {
                        // minimum on some other processor
                        jGlob[1] = int.MaxValue;
                    }


                    // in case of multiple global minimums/maximums, e.g. for a constant field, we return the lowest (jMaxGlob,jMinGlob) 
                    int[] jGlobM = new int[2];
                    unsafe
                    {
                        fixed (int* ploc = jGlob, ptot = jGlobM) {
                            csMPI.Raw.Allreduce((IntPtr)ploc, (IntPtr)ptot, 2, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                        }
                    }
                    jMaxGlob = jGlob[0];
                    jMinGlob = jGlob[1];

                } else {
                    min = LocMin;
                    max = LocMax;

                    jMaxGlob = jMaxLoc;
                    jMinGlob = jMinLoc;
                }
            }
        }

        /// <summary>
        /// Returns the minimum and the maximum value of <paramref name="mod"/>(this DG field), for each cell;
        /// </summary>
        /// <param name="max">
        /// Vector, with the same number of entries as the cell mask <paramref name="cm"/>;
        /// On exit, the approximate local maximum within the cell.
        /// </param>
        /// <param name="min">
        /// Vector, with the same number of entries as the cell mask <paramref name="cm"/>;
        /// On exit, the approximate local minimum within the cell.
        /// </param>
        /// <param name="cm">
        /// optional domain restriction
        /// </param>
        /// <param name="mod">
        /// optimal modifier, which e.g. may select the absolute value; if null, the identity function;
        /// </param>
        public void GetCellwiseExtremalValues<T1, T2>(T1 min, T2 max, CellMask cm = null, Func<double, double> mod = null)
            where T1 : IList<double>
            where T2 : IList<double> //
        {
            using (new FuncTrace()) {
                int J = Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                int Jcm = cm == null ? J : cm.NoOfItemsLocally;
                if (min.Count != Jcm) {
                    throw new ArgumentException("wrong length of output vector", "min");
                }
                if (max.Count != Jcm) {
                    throw new ArgumentException("wrong length of output vector", "max");
                }
                min.SetAll(double.MaxValue);
                max.SetAll(double.MinValue);

                if (mod == null)
                    mod = (x => x);

                //int N = m_context.Grid.GridSimplex.NoOfVertices + 1;
                // find maximum on this processor
                // ------------------------------
                // create node set
                if (m_ExtremalProbeNS == null) {
                    // nodes for the evaluation of the velocity vector
                    // (all vertices an 0)
                    int D = Basis.GridDat.SpatialDimension;

                    var KrefS = Basis.GridDat.iGeomCells.RefElements;
                    m_ExtremalProbeNS = new NodeSet[KrefS.Length];
                    for (int i = 0; i < KrefS.Length; i++) {
                        m_ExtremalProbeNS[i] = new NodeSet(KrefS[i], KrefS[i].Vertices);
                    }
                }

                // vectorisation of this method
                int VectorSize = -1;
                int N0 = m_ExtremalProbeNS[0].NoOfNodes;   // number of nodes
                MultidimensionalArray fieldValues = new MultidimensionalArray(2);


                IEnumerable<Chunk> all_chunks;
                if (cm == null) {
                    var _ch = new Chunk[1];
                    _ch[0].i0 = 0;
                    _ch[0].Len = J;
                    all_chunks = _ch;
                } else {
                    all_chunks = cm;
                }

                int jSub = 0;
                foreach (Chunk chk in all_chunks) {
                    VectorSize = this.GridDat.iGeomCells.GetNoOfSimilarConsecutiveCells(CellInfo.RefElementIndex_Mask, chk.i0, Math.Min(100, chk.Len)); // try to be a little vectorized;

                    int _J = chk.Len + chk.i0;
                    for (int j0 = chk.i0; j0 < _J; j0 += VectorSize) {
                        if (j0 + VectorSize > _J)
                            VectorSize = _J - j0;
                        int iKref = Basis.GridDat.iGeomCells.GetRefElementIndex(j0);
                        int N = m_ExtremalProbeNS[iKref].NoOfNodes;

                        if (fieldValues.GetLength(0) != VectorSize || fieldValues.GetLength(1) != N)
                            fieldValues = MultidimensionalArray.Create(VectorSize, N);

                        this.Evaluate(j0, VectorSize, m_ExtremalProbeNS[iKref], fieldValues, 0.0);

                        // loop over cells ...
                        for (int jCell = j0; jCell < j0 + VectorSize; jCell++) {

                            // loop over nodes ...
                            for (int n = 0; n < N; n++) {
                                double vel = 0;
                                vel = mod(fieldValues[jCell - j0, n]);
                                min[jSub] = Math.Min(min[jSub], vel);
                                max[jSub] = Math.Max(max[jSub], vel);
                            }

                            jSub++;
                        }
                    }
                }

            }
        }

        /// <summary>
        /// Computes L2 distance between this field and
        /// <paramref name="other"/>
        /// </summary>
        /// <param name="other"></param>
        /// <param name="cm">
        /// Optional restriction of domain
        /// </param>
        public double L2Error(DGField other, CellMask cm = null) {
            int order = Math.Max(this.Basis.Degree, other.Basis.Degree) * 2;
            CellQuadratureScheme cqs = new CellQuadratureScheme(domain: cm);
            return LxError(other.Evaluate, null, cqs.Compile(this.GridDat, order));
        }

        /// <summary>
        /// Computes L2 measure with default quadrature
        /// </summary>
        public double L2Error(ScalarFunction function, int order, CellQuadratureScheme scheme = null) {
            return LxError(function, null, scheme.SaveCompile(this.GridDat, order)).Sqrt();
        }

        /// <summary>
        /// Computes L2 measure with default quadrature
        /// </summary>
        /// <param name="function"></param>
        /// <param name="scheme"></param>
        /// <returns></returns>
        public double L2Error(ScalarFunction function, CellQuadratureScheme scheme = null) {
            return LxError(function, null, scheme.SaveCompile(this.GridDat, this.Basis.Degree * 2 + 1)).Sqrt();
        }

        /// <summary>
        /// This call computes an integral measure which may depend on 
        /// this <see cref="DGField"/> an the given <paramref name="function"/>;
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        /// <param name="function">
        /// Optional: Reference function for error computation. If null, zero
        /// is taken as the reference function
        /// </param>
        /// <param name="rule">
        /// composite quadrature rule.
        /// </param>
        /// <param name="Map">
        /// Arbitrary mapping applied to the values of this field and
        /// <paramref name="function"/> at some point, which is finally integrated.
        /// E.g., the mapping for an L2-error would be \f$ (a,b) => (a - b)^2 \f$, 
        /// where \f$ a \f$ is the value of this field at some point \f$ \vec{x} \f$ and
        /// \f$ b \f$ is the value of <paramref name="function"/> at \f$ \vec{x} \f$.
        /// </param>
        /// <returns>
        /// on all invoking MPI processes, the L2 norm of
        /// this field minus <paramref name="function"/>, approximated by the
        /// quadrature rule <paramref name="rule"/>
        /// </returns>
        public double LxError(ScalarFunction function, Func<double, double, double> Map, ICompositeQuadRule<QuadRule> rule) {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            using (new FuncTrace()) {
                LxNormQuadrature l2nq = new LxNormQuadrature(this, function, Map, rule);
                l2nq.Execute();

                double nrmtot = l2nq.LxNorm.MPISum();
                return nrmtot;
            }
        }

        /// <summary>
        /// This call computes an integral measure which may depend on 
        /// this <see cref="DGField"/> an the given <paramref name="function"/>;
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        /// <param name="function"></param>
        /// <param name="rule">
        /// composite quadrature rule.
        /// </param>
        /// <param name="Map">
        /// Arbiter mapping applied to the values of this field and
        /// <paramref name="function"/> at some point, which is finally integrated.
        /// E.g., the mapping for an L2-error would be \f$ (a,b) => (a - b)^2 \f$, 
        /// where \f$ a \f$ is the value of this field at some point \f$ \vec{x} \f$ and
        /// \f$ b \f$ is the value of <paramref name="function"/> at \f$ \vec{x} \f$.
        /// </param>
        /// <returns>
        /// on all invoking MPI processes, the L2 norm of
        /// this field minus <paramref name="function"/>
        /// </returns>
        public double LxError(ScalarFunctionEx function, Func<double, double, double> Map, ICompositeQuadRule<QuadRule> rule) {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            using (new FuncTrace()) {
                LxNormQuadrature l2nq = new LxNormQuadrature(this, function, Map, rule);
                l2nq.Execute();

                double nrmtot = l2nq.LxNorm.MPISum();
                return nrmtot;
            }
        }

        /// <summary>
        /// Variant of
        /// <see cref="LxError(ScalarFunction, Func{double, double, double}, ICompositeQuadRule{QuadRule})"/>
        /// that computes the cell-local Lx norm (and does thus not include MPI
        /// communication)
        /// </summary>
        /// <param name="function">
        /// Optional: Reference function for error computation. If null, zero
        /// is taken as the reference function
        /// </param>
        /// <param name="quadRule">
        /// composite quadrature rule.
        /// </param>
        /// <param name="Map">
        /// Arbitrary mapping applied to the values of this field and
        /// <paramref name="function"/> at some point, which is finally
        /// integrated. E.g., the mapping for an L2-error would be
        /// \f$ (a,b) => (a - b)^2 \f$, where \f$ a \f$ is the value of this
        /// field at some point \f$ \vec{x} \f$ and \f$ b \f$ is the value of
        /// <paramref name="function"/> at \f$ \vec{x} \f$.
        /// </param>
        /// <returns>
        /// The cell-local Lx norm of this field minus
        /// <paramref name="function"/>, approximated by the quadrature rule
        /// <paramref name="quadRule"/>
        /// </returns>
        public double[] LocalLxError(ScalarFunction function, Func<double, double, double> Map, ICompositeQuadRule<QuadRule> quadRule) {
            using (new FuncTrace()) {
                LocalLxNormQuadrature quacdrature = new LocalLxNormQuadrature(this, function, Map, quadRule);
                quacdrature.Execute();

                return quacdrature.LocalLxNorms;
            }
        }

        /// <summary>
        /// L2 - norm of this field;
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        virtual public double L2Norm() {
            return L2Norm(default(CellMask));
        }

        /// <summary>
        /// L2 - norm of this field;
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        /// <returns>
        /// on all invoking MPI processes, the L2 norm of
        /// this field;
        /// </returns>
        /// <param name="CM">
        /// optional cell mask
        /// </param>
        virtual public double L2Norm(CellMask CM) {
            CellQuadratureScheme scheme = null;
            if (CM != null)
                scheme = new CellQuadratureScheme(true, CM);
            double r = LxError((ScalarFunction)null, null, scheme.SaveCompile(this.GridDat, 2 * m_Basis.Degree)).Sqrt();
            return r;
        }

        /// <summary>
        /// L1 - norm of this field;
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        virtual public double L1Norm() {
            return L1Norm(default(CellMask));
        }

        /// <summary>
        /// L1 - norm of this field;
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        /// <returns>
        /// on all invoking MPI processes, the L1 norm of
        /// this field;
        /// </returns>
        /// <param name="CM">
        /// optional cell mask
        /// </param>
        virtual public double L1Norm(CellMask CM) {
            CellQuadratureScheme scheme = null;
            if (CM != null)
                scheme = new CellQuadratureScheme(true, CM);
            double r = LxError((ScalarFunction)null, (a, b) => a.Abs(), scheme.SaveCompile(this.GridDat, 2 * m_Basis.Degree));
            return r;
        }

        /// <summary>
        /// L-infinity - norm of this field;<br/>
        /// This is a collective call, it must be invoked by all 
        /// MPI processes within the communicator; internally, it invokes MPI_Allreduce;
        /// </summary>
        /// <returns>
        /// on all invoking MPI processes, the L-inv norm of
        /// this field (over all processors);
        /// </returns>
        /// <remarks>
        /// to find the maximum value, this field is evaluated on each vertex and the center of the simplex.
        /// </remarks>
        public double LinfNorm(CellMask cm = null) {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            using (new FuncTrace()) {
                double Min, Max;
                int j, k;
                GetExtremalValues(out Min, out Max, out j, out k, cm);

                return Math.Max(Math.Abs(Max), Math.Abs(Min));
            }
        }

        /// <summary>
        /// quadrature that calculates some norm such as the L2 norm
        /// </summary>
        class LxNormQuadrature : CellQuadrature {

            /// <summary>
            /// ctor.
            /// </summary>
            public LxNormQuadrature(DGField owner, ScalarFunction func, Func<double, double, double> Map, ICompositeQuadRule<QuadRule> rule)
                : base(new int[] { 1 }, owner.Basis.GridDat, rule) //
            {
                m_func = func;
                m_Owner = owner;
                m_Map = Map;
            }

            /// <summary>
            /// ctor.
            /// </summary>
            public LxNormQuadrature(DGField owner, ScalarFunctionEx func, Func<double, double, double> Map, ICompositeQuadRule<QuadRule> rule)
                : base(new int[] { 1 }, owner.Basis.GridDat, rule) //
            {
                m_funcEx = func;
                m_Owner = owner;
                m_Map = Map;
            }

            DGField m_Owner;

            ScalarFunction m_func;
            ScalarFunctionEx m_funcEx;

            Func<double, double, double> m_Map;

            /// <summary>
            /// 1st index: cell index (minus some offset);
            /// 2nd index: node index;
            /// 3rd index; spatial coordinate;
            /// </summary>
            MultidimensionalArray m_NodesTransformed = new MultidimensionalArray(3);

            protected override void AllocateBuffers(int NoOfItems, NodeSet rule) {
                base.AllocateBuffers(NoOfItems, rule);

                int NoOfNodes = rule.GetLength(0);
                if (m_func != null) {
                    m_NodesTransformed.Allocate(new int[] { NoOfItems, NoOfNodes, GridDat.SpatialDimension });
                }
            }

            /// <summary>
            /// Integrand evaluation.
            /// </summary>
            protected override void Evaluate(int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) {
                NodeSet NodesUntransformed = rule.Nodes;
                int M = NodesUntransformed.GetLength(0);
                int D = NodesUntransformed.GetLength(1);

                // evaluate scalar function ans store result in 'EvalResult'
                // =========================================================
                Debug.Assert(!((m_func != null) && (m_funcEx != null)));
                if (m_func != null) {
                    GridDat.TransformLocal2Global(NodesUntransformed, i0, Length, m_NodesTransformed, 0);

                    MultidimensionalArray inp = m_NodesTransformed.ResizeShallow(new int[] { Length * M, D });
                    MultidimensionalArray outp = EvalResult.ResizeShallow(new int[] { Length * M });
                    m_func(inp, outp);
                }
                if (m_funcEx != null) {
                    m_funcEx(i0, Length, NodesUntransformed, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                }


                if (m_Map == null) {
                    // L2 error branch
                    // +++++++++++++++

                    MultidimensionalArray fieldvals = EvalResult.ResizeShallow(new int[] { Length, NodesUntransformed.GetLength(0) });
                    m_Owner.Evaluate(i0, Length, NodesUntransformed, fieldvals, -1.0);

                    for (int j = 0; j < Length; j++) {
                        for (int m = 0; m < M; m++) {
                            double e;
                            e = EvalResult[j, m, 0];
                            EvalResult[j, m, 0] = e * e;
                        }
                    }
                } else {
                    // arbitrary mapping branch
                    // ++++++++++++++++++++++++

                    MultidimensionalArray fieldvals = MultidimensionalArray.Create(new int[] { Length, NodesUntransformed.GetLength(0) });
                    m_Owner.Evaluate(i0, Length, NodesUntransformed, fieldvals, -1.0);


                    for (int j = 0; j < Length; j++) {
                        for (int m = 0; m < M; m++) {
                            double e;
                            e = EvalResult[j, m, 0];
                            EvalResult[j, m, 0] = this.m_Map(fieldvals[j, m], e);
                        }
                    }
                }
            }

            double m_L2pow2;

            /// <summary>
            /// returns the L2 Norm to the power of 2
            /// </summary>
            public double LxNorm {
                get {
                    return m_L2pow2;
                }
            }

            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int j = 0; j < Length; j++) {
                    m_L2pow2 += ResultsOfIntegration[j, 0];
                }
            }
        }

        /// <summary>
        /// Quadrature for
        /// <see cref="LocalLxError(ScalarFunction, Func{double, double, double}, ICompositeQuadRule{QuadRule})"/>
        /// </summary>
        class LocalLxNormQuadrature : LxNormQuadrature {

            public LocalLxNormQuadrature(DGField owner, ScalarFunction func, Func<double, double, double> Map, ICompositeQuadRule<QuadRule> rule)
                : base(owner, func, Map, rule) {
                m_localLxNorms = new double[rule.NumberOfItems];
            }

            public override void Execute() {
                m_localLxNorms.Clear();
                subgridCellIndex = 0;
                base.Execute();
            }

            double[] m_localLxNorms;
            
            public double[] LocalLxNorms {
                get {
                    return m_localLxNorms;
                }
            }

            private int subgridCellIndex = 0;

            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int j = 0; j < Length; j++) {
                    m_localLxNorms[subgridCellIndex] += ResultsOfIntegration[j, 0];
                    subgridCellIndex++;
                }
            }
        }


        /// <summary>
        /// integrates this field over the domain specified in <paramref name="volumemask"/>
        /// </summary>
        /// <param name="volumemask">
        /// an optional volume mask; if null, the whole grid is taken;
        /// </param>
        public double IntegralOver(CellMask volumemask) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                if (volumemask == null) {
                    volumemask = CellMask.GetFullMask(Basis.GridDat);
                }

                double acc = 0;
                foreach (var chunk in volumemask) {
                    int iE = chunk.i0 + chunk.Len;
                    for (int j = chunk.i0; j < iE; j++) {
                        double mv = GetMeanValue(j);
                        double vol = Basis.GridDat.iLogicalCells.GetCellVolume(j);
                        acc += mv * vol;
                    }
                }

                double accglob = double.NaN;
                unsafe
                {
                    csMPI.Raw.Allreduce((IntPtr)(&acc), (IntPtr)(&accglob), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                }

                return accglob;
            }
        }

        /// <summary>
        /// used by method <see cref="IntegralOverEx"/>
        /// </summary>
        class IntegralOverExQuadrature : BoSSS.Foundation.Quadrature.CellQuadrature {

            internal IntegralOverExQuadrature(IGridData dat, DGField[] fields, ICompositeQuadRule<QuadRule> qr, Func f)
                : base(new int[] { 1 }, dat, qr) {
                m_fields = fields;
                m_f = f;
            }

            DGField[] m_fields;



            MultidimensionalArray[] evalRes;

            MultidimensionalArray X;

            protected override void AllocateBuffers(int NoOfItems, NodeSet rule) {
                base.AllocateBuffers(NoOfItems, rule);
                int NoOfNodes = rule.GetLength(0);
                if (evalRes == null)
                    evalRes = new MultidimensionalArray[m_fields.Length];
                for (int f = 0; f < m_fields.Length; f++)
                    evalRes[f] = MultidimensionalArray.Create(NoOfItems, NoOfNodes);
                int D = GridDat.SpatialDimension;
                X = MultidimensionalArray.Create(NoOfItems, NoOfNodes, D);
            }

            Func m_f;

            protected override void Evaluate(int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) {
                NodeSet Nodes = rule.Nodes;
                int F = m_fields.Length;
                int D = GridDat.SpatialDimension;

                for (int f = 0; f < F; f++) {
                    m_fields[f].Evaluate(i0, Length, Nodes, evalRes[f]);
                }

                GridDat.TransformLocal2Global(Nodes, i0, Length, X, 0);
                int NoOfNodes = Nodes.NoOfNodes;

                double[] _x = new double[D];
                double[] args = new double[F];
                for (int i = 0; i < Length; i++) {
                    for (int n = 0; n < NoOfNodes; n++) {
                        for (int f = 0; f < F; f++)
                            args[f] = evalRes[f][i, n];
                        for (int d = 0; d < D; d++) {
                            _x[d] = X[i, n, d];
                        }

                        EvalResult[i, n, 0] = m_f(_x, args, i + i0);
                    }
                }
            }

            public double result = 0;

            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++)
                    result += ResultsOfIntegration[i, 0];
            }
        }

        /// <summary>
        /// integrates some arbitrary function <paramref name="f"/> (which
        /// depends on DG fields <paramref name="Fields"/>) over some domain
        /// implied by <paramref name="scheme"/>
        /// </summary>
        /// <param name="scheme">
        /// specification of quadrature rule and integration domain
        /// </param>
        /// <param name="f">
        /// integrand
        /// </param>
        /// <param name="Fields">
        /// input arguments for function <paramref name="f"/>
        /// </param>
        /// <returns>
        /// the integral of <paramref name="f"/> over the domain implied by
        /// <paramref name="scheme"/>
        /// </returns>
        /// <param name="OverIntegrationMultiplier">
        /// Multiplier for Quadrature Order
        /// </param>
        public static double IntegralOverEx(CellQuadratureScheme scheme, Func f, int OverIntegrationMultiplier = 2, params DGField[] Fields) {
            using (new FuncTrace()) {
                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                var g = Fields[0].Basis.GridDat;

                int order = Fields.Max(x => x.Basis.Degree) * OverIntegrationMultiplier;

                for (int i = 1; i < Fields.Length; i++)
                    if (!Object.ReferenceEquals(g, Fields[i].Basis.GridDat))
                        throw new ArgumentException("all fields must be defined on the same grid");



                IntegralOverExQuadrature q = new IntegralOverExQuadrature(g, Fields, scheme.SaveCompile(g, order), f);
                q.Execute();

                unsafe
                {
                    double locRes = q.result, glRes = 0;
                    csMPI.Raw.Allreduce((IntPtr)(&locRes), (IntPtr)(&glRes), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);

                    return glRes;
                }

            }
        }
    }
}
