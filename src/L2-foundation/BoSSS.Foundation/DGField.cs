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
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation {

    /// <summary>
    /// This class represents the base class of all DG fields.
    /// </summary>
    abstract public partial class DGField : ICloneable {

        /// <summary>
        /// constructs a new field. 
        /// </summary>
        /// <param name="__Basis">The basis that is used for this field</param>
        /// <param name="__Identification">
        /// identification string for this field;
        /// This can be null or empty, 
        /// however, if IO should be preformed for this object, the identification must be unique 
        /// within the given <see cref="GridData"/>-instance
        /// </param>
        protected DGField(Basis __Basis, String __Identification) {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            // set members
            m_Basis = __Basis;

            if (__Identification != null && __Identification.Length < 1) {
                __Identification = null;
            }
            m_Identification = __Identification;
        }

        /// <summary>
        /// Degrees-of-freedom for storing this field.
        /// </summary>
        public int DOFLocal {
            get {
                int max = m_Basis.MaximalLength;
                int min = m_Basis.MinimalLength;
                int J = Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                if (min == max) {
                    return J * max;
                } else {
                    int sum = 0;
                    for (int j = 0; j < J; j++) {
                        sum += m_Basis.GetLength(j);
                    }
                    return sum;
                }
            }
        }

        /// <summary>
        /// projects a scalar function onto the DG-space;
        /// It's used by
        /// <see cref="DGField.ProjectField(double,ScalarFunction,CellQuadratureScheme)"/>.
        /// </summary>
        protected class ProjectionQuadrature : BoSSS.Foundation.Quadrature.CellQuadrature {

            /// <summary>
            /// 
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="func"></param>
            /// <param name="alpha"></param>
            /// <param name="qr">quad. rule to be used</param>
            public ProjectionQuadrature(DGField owner, double alpha, ScalarFunction func, ICompositeQuadRule<QuadRule> qr)
                : base(new int[] { owner.m_Basis.Length }, owner.Basis.GridDat, qr) {
                m_func = func;
                m_Owner = owner;
                m_alpha = alpha;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="func"></param>
            /// <param name="alpha"></param>
            /// <param name="qr">quad. rule to be used</param>
            public ProjectionQuadrature(DGField owner, double alpha, ScalarFunctionEx func, ICompositeQuadRule<QuadRule> qr)
                : base(new int[] { owner.m_Basis.Length }, owner.Basis.GridDat, qr) {
                m_func2 = func;
                m_Owner = owner;
                m_alpha = alpha;
            }

            /// <summary>
            /// 
            /// </summary>
            protected DGField m_Owner;

            /// <summary>
            /// the function to project onto <see cref="m_Owner"/>
            /// </summary>
            protected ScalarFunction m_func;

            /// <summary>
            /// the function to project onto <see cref="m_Owner"/>
            /// </summary>
            protected ScalarFunctionEx m_func2;

            /// <summary>
            /// a constant to multiply <see cref="m_func"/> or <see cref="m_func2"/> with;
            /// </summary>
            protected double m_alpha;




            /// <summary>
            /// 1st index: cell index (minus some offset);
            /// 2nd index: node index;
            /// 3rd index; spatial coordinate;
            /// </summary>
            protected MultidimensionalArray m_NodesTransformed = new MultidimensionalArray(3);

            /// <summary>
            /// results of function evaluation
            /// 1st index: cell index (minus some offset);
            /// 2nd index: node index;
            /// </summary>
            protected MultidimensionalArray m_FunctionValues = new MultidimensionalArray(2);

            /// <summary>
            /// Allocates memory for the global coordinates and the function values
            /// </summary>
            /// <param name="NoOfItems"></param>
            /// <param name="rule"></param>
            protected override void AllocateBuffers(int NoOfItems, NodeSet rule) {
                base.AllocateBuffers(NoOfItems, rule);
                int NoOfNodes = rule.GetLength(0);
                m_NodesTransformed.Allocate(new int[] { NoOfItems, NoOfNodes, GridDat.SpatialDimension });
                m_FunctionValues.Allocate(new int[] { NoOfItems, NoOfNodes });
            }

            /// <summary>
            /// Integrand evaluation.
            /// </summary>
            protected override void Evaluate(int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) {
                NodeSet NodesUntransformed = rule.Nodes;
                if (Length != EvalResult.GetLength(0))
                    throw new ApplicationException();
                if (Length != m_NodesTransformed.GetLength(0))
                    throw new ApplicationException();


                int N = NodesUntransformed.GetLength(0); // number of nodes per cell
                int D = NodesUntransformed.GetLength(1); // spatial dimension
                int M = EvalResult.GetLength(2); // number of integrals per Cell

                if (m_func != null) {
                    // transform nodes
                    GridDat.TransformLocal2Global(NodesUntransformed, i0, Length, m_NodesTransformed, 0);

                    // evaluate function
                    m_func(m_NodesTransformed.ResizeShallow(new int[] { N * Length, D }),
                        m_FunctionValues.ResizeShallow(new int[] { N * Length }));
                } else {
                    // evaluate function
                    m_func2(i0, Length, NodesUntransformed, m_FunctionValues);
                }

                // get basis values
                MultidimensionalArray BasisValues = m_Owner.Basis.CellEval(NodesUntransformed, i0, Length);

                // (ScalarField) * (m-th Basis function)                
                for (int j = 0; j < Length; j++) {   // loop over cells
                    for (int n = 0; n < N; n++) {    // loop over quadrature nodes
                        for (int m = 0; m < M; m++) { // loop over coordinates
                            EvalResult[j, n, m] = m_FunctionValues[j, n] * BasisValues[j, n, m];
                        }
                    }
                }

            }

            /// <summary>
            /// performs the accumulation (multiplication of
            /// integration result by <see cref="m_alpha"/> and addition to)
            /// the DG coordinates of <see cref="m_Owner"/>;
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="ResultsOfIntegration"></param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                int N = m_Owner.Basis.MaximalLength; // number of integrals per Cell
                if (N != ResultsOfIntegration.GetLength(1))
                    throw new ApplicationException("internal error.");

                for (int j = 0; j < Length; j++) {
                    for (int n = 0; n < N; n++) {
                        m_Owner.Coordinates[j + i0, n] = m_Owner.Coordinates[j + i0, n] + m_alpha * ResultsOfIntegration[j, n];
                    }
                }
            }
        }

        /// <summary>
        /// Initializes the field to be the projection of a 
        /// given <see cref="ScalarFunction"/> onto the Discontinuous Galerkin Space.<br/>
        /// The original content of this field is cleared;
        /// </summary>
        /// <param name="func"></param>
        public void ProjectField(ScalarFunction func) {
            this.Clear();
            ProjectField(1.0, func, default(CellQuadratureScheme));
        }

        /// <summary>
        /// Initializes the field to be the projection of a 
        /// given <see cref="ScalarFunction"/> onto the Discontinuous Galerkin Space.<br/>
        /// The original content of this field is cleared;
        /// </summary>
        /// <param name="func"></param>
        public void ProjectField(ScalarFunctionEx func) {
            this.Clear();
            ProjectField(1.0, func, default(CellQuadratureScheme));
        }

        /// <summary>
        /// creates a <see cref="ProjectionQuadrature"/> which is used by
        /// <see cref="ProjectField(double,ScalarFunction,CellQuadratureScheme)"/>;
        /// </summary>
        virtual protected ProjectionQuadrature GetProjectionQuadrature(double alpha, ScalarFunction func, ICompositeQuadRule<QuadRule> qr) {
            return new ProjectionQuadrature(this, alpha, func, qr);
        }

        /// <summary>
        /// creates a <see cref="ProjectionQuadrature"/> which is used by 
        /// <see cref="ProjectField(double,ScalarFunction,CellQuadratureScheme)"/>;
        /// </summary>
        virtual protected ProjectionQuadrature GetProjectionQuadrature(double alpha, ScalarFunctionEx func, ICompositeQuadRule<QuadRule> qr) {
            return new ProjectionQuadrature(this, alpha, func, qr);
        }

        /// <summary>
        /// Accumulates the DG projection (with respect to <see cref="Basis"/>)
        /// times <paramref name="alpha"/> to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>*<paramref name="func"/>;
        /// </summary>
        /// <param name="func"></param>
        /// <param name="alpha">scaling of <paramref name="func"/></param>
        /// <param name="scheme"></param>
        virtual public void ProjectField(double alpha, ScalarFunction func, CellQuadratureScheme scheme = null) {
            using (new FuncTrace()) {
                int order = this.Basis.Degree * 2 + 2;
                ProjectionQuadrature pq = GetProjectionQuadrature(alpha, func, scheme.SaveCompile(this.Basis.GridDat, order));
                pq.Execute();
            }
        }

        /// <summary>
        /// Accumulates the DG projection (with respect to <see cref="Basis"/>)
        /// times <paramref name="alpha"/> to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>*<paramref name="func"/>;
        /// </summary>
        /// <param name="func"></param>
        /// <param name="alpha">scaling of <paramref name="func"/></param>
        /// <param name="rule">
        /// quadrature rule
        /// </param>
        virtual public void ProjectField(double alpha, ScalarFunction func, ICompositeQuadRule<QuadRule> rule) {
            using (new FuncTrace()) {
                ProjectionQuadrature pq = GetProjectionQuadrature(alpha, func, rule);
                pq.Execute();
            }
        }

        /// <summary>
        /// Accumulates the DG projection (with respect to <see cref="Basis"/>)
        /// times <paramref name="alpha"/> to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>*<paramref name="func"/>;
        /// </summary>
        /// <param name="func"></param>
        /// <param name="alpha">scaling of <paramref name="func"/></param>
        /// <param name="rule">
        /// quadrature rule
        /// </param>
        public void ProjectField(double alpha, ScalarFunctionEx func, ICompositeQuadRule<QuadRule> rule) {
            using (new FuncTrace()) {
                ProjectionQuadrature pq = GetProjectionQuadrature(alpha, func, rule);
                pq.Execute();
            }
        }

        /// <summary>
        /// Accumulates the DG projection (with respect to <see cref="Basis"/>)
        /// times <paramref name="alpha"/> to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>*<paramref name="func"/>;
        /// </summary>
        /// <param name="func"></param>
        /// <param name="alpha">scaling of <paramref name="func"/></param>
        /// <param name="scheme"></param>
        public void ProjectField(double alpha, ScalarFunctionEx func, CellQuadratureScheme scheme = null) {
            using (new FuncTrace()) {
                int order = this.Basis.Degree * 2 + 2;
                ProjectionQuadrature pq = GetProjectionQuadrature(alpha, func, scheme.SaveCompile(this.Basis.GridDat, order));
                pq.Execute();
            }
        }

        /// <summary>
        /// Performs a nodal projection, i.e. accumulates a DG field, defined by
        /// \f[ 
        ///     \textrm{in each cell } K_j: \ u(\vec{\xi}_i) = \alpha f(\vec{\xi}_i) \quad \forall \vec{\xi}_i
        /// \f]
        /// to this field.
        /// If the number of nodes is not equal to the degrees-of-freedom in a specific cell,
        /// a least-square projection is performed.
        /// </summary>
        /// <param name="alpha">scaling of <paramref name="func"/></param>
        /// <param name="func">function to be projected</param>
        /// <param name="NodeSet">
        /// cell-local coordinates \f$ (\vec{xi}_1, \ldots , \vec{xi}_M) \f$
        /// </param>
        public void ProjectNodal(double alpha, ScalarFunction func, NodeSet NodeSet) {
            using (new FuncTrace()) {
                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                var B = this.Basis;
                int D = this.GridDat.SpatialDimension;


                int K = NodeSet.GetLength(0);

                var evalB = MultidimensionalArray.Create(1, K);
                var NodesGlob = MultidimensionalArray.Create(K, D);

                for (int j = 0; j < J; j++) { // loop over cells...
                    int N = B.GetLength(j);
                    evalB.Clear();
                    this.GridDat.TransformLocal2Global(NodeSet, NodesGlob, j);
                    func(NodesGlob, evalB.ExtractSubArrayShallow(0, -1));

                    var CellVal = B.CellEval(NodeSet, j, 1).ExtractSubArrayShallow(0, -1, -1);
                    MultidimensionalArray system = CellVal;
                    Debug.Assert(system.NoOfRows == K);
                    Debug.Assert(system.NoOfCols == N);

                    double[] RHS = evalB.GetRow(0);
                    Debug.Assert(RHS.Length == system.NoOfRows);

                    double[] X = new double[N];

                    if (system.NoOfCols == system.NoOfRows) {
                        system.Solve(X, RHS);
                    } else {
                        system.LeastSquareSolve(X, RHS);
                    }

                    for (int n = 0; n < N; n++)
                        this.Coordinates[j, n] += alpha * X[n];
                }

            }
        }

        /// <summary>
        /// Performs a nodal projection, i.e. accumulates a DG field, defined by
        /// \f[ 
        ///     \textrm{in each cell } K_j: \ u(\vec{\xi}_i) = \alpha f(\vec{\xi}_i) \quad \forall \vec{\xi}_i
        /// \f]
        /// to this field. If the number of nodes is not equal to the
        /// degrees-of-freedom in a specific cell, a least-square projection is
        /// performed.
        /// </summary>
        /// <param name="alpha">scaling of <paramref name="func"/></param>
        /// <param name="func">function to be projected</param>
        /// <param name="NodeSet">
        /// cell-local coordinates \f$ (\vec{xi}_1, \ldots , \vec{xi}_M) \f$
        /// </param>
        public void ProjectNodal(double alpha, ScalarFunctionEx func, NodeSet NodeSet) {
            using (new FuncTrace()) {
                int noOfCells = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                int D = this.GridDat.SpatialDimension;
                int noOfNodesPerCell = NodeSet.GetLength(0);

                MultidimensionalArray evalResult = MultidimensionalArray.Create(noOfCells, noOfNodesPerCell);
                func(0, noOfCells, NodeSet, evalResult);

                for (int j = 0; j < noOfCells; j++) {
                    int N = Basis.GetLength(j);

                    MultidimensionalArray matrix =
                        Basis.CellEval(NodeSet, j, 1).ExtractSubArrayShallow(0, -1, -1);
                    Debug.Assert(matrix.NoOfRows == noOfNodesPerCell);
                    Debug.Assert(matrix.NoOfCols == N);

                    double[] RHS = evalResult.GetRow(j);
                    Debug.Assert(RHS.Length == matrix.NoOfRows);

                    double[] X = new double[N];

                    if (matrix.NoOfCols == matrix.NoOfRows) {
                        matrix.Solve(X, RHS);
                    } else {
                        matrix.LeastSquareSolve(X, RHS);
                    }

                    for (int n = 0; n < N; n++) {
                        this.Coordinates[j, n] += alpha * X[n];
                    }
                }
            }
        }

        /// <summary>
        /// Similar to <see cref="ProjectNodal(double, ScalarFunctionEx, NodeSet)"/>,
        /// but uses the DG coordinates associated with multi-linear basis functions
        /// (i.e., constant and linear term in 1D; constant, linear and
        /// bilinear terms in 2D; constant, linear, bilinear and trilinear
        /// terms in 3D)
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="func"></param>
        /// <param name="NodeSet"></param>
        public void ProjectNodalMultilinear(double alpha, ScalarFunctionEx func, NodeSet NodeSet) {
            using (new FuncTrace()) {
                int D = this.GridDat.SpatialDimension;
                if (this.Basis.Degree < D) {
                    throw new Exception();
                }

                int noOfCells = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                int noOfNodesPerCell = NodeSet.GetLength(0);

                MultidimensionalArray evalResult = MultidimensionalArray.Create(noOfCells, noOfNodesPerCell);
                func(0, noOfCells, NodeSet, evalResult);

                var linearIndices = new List<int>();
                for (int n = 0; n < Basis.Length; n++) {
                    bool isMultilinear = true;
                    for (int c = 0; c < Basis.Polynomials[0][n].Exponents.GetLength(0); c++) {
                        for (int d = 0; d < D; d++) {
                            if (Basis.Polynomials[0][n].Exponents[c, d] > 1) {
                                isMultilinear = false;
                                break;
                            }
                        }

                        if (!isMultilinear) {
                            break;
                        }
                    }

                    if (isMultilinear) {
                        linearIndices.Add(n);
                    }
                }

                for (int i = 0; i < noOfCells; i++) {
                    MultidimensionalArray basisValues = Basis.CellEval(NodeSet, i, 1);
                    MultidimensionalArray matrix = MultidimensionalArray.Create(
                        NodeSet.NoOfNodes, linearIndices.Count);
                    for (int j = 0; j < NodeSet.NoOfNodes; j++) {
                        for (int k = 0; k < linearIndices.Count; k++) {
                            matrix[j, k] = basisValues[0, j, linearIndices[k]];
                        }
                    }

                    double[] RHS = evalResult.GetRow(i);

                    Debug.Assert(matrix.NoOfCols == matrix.NoOfRows);
                    Debug.Assert(RHS.Length == matrix.NoOfRows);

                    double[] X = new double[linearIndices.Count];
                    matrix.Solve(X, RHS);

                    for (int n = 0; n < NodeSet.NoOfNodes; n++) {
                        this.Coordinates[i, linearIndices[n]] += alpha * X[n];
                    }
                }
            }
        }

        /// <summary>
        /// overwrites the mean value in one cell
        /// </summary>
        /// <param name="j">a local cell index</param>
        /// <param name="v">the new mean value in this cell</param>
        virtual public void SetMeanValue(int j, double v) {
            if (j < 0 || j >= Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells)
                throw new ArgumentException("j must be in the range of locally updated cells.", "j");

            int iKref = this.Basis.GridDat.iGeomCells.GetRefElementIndex(j);
            double bv = m_Basis.Polynomials[iKref][0].Coeff[0];
            Debug.Assert(m_Basis.Polynomials[iKref][0].AbsoluteDegree == 0, "Polynomial degree of 0-th polynomial is expected to be 0.");
            Debug.Assert(m_Basis.Polynomials[iKref][0].Coeff.Length == 1, "Polynomial degree of 0-th polynomial is expected to be 0.");

            double sc;
            if (this.Basis.GridDat.iGeomCells.IsCellAffineLinear(j)) {
                sc = Basis.Data.Scaling[j];
            } else {
                sc = Basis.Data.OrthonormalizationTrafo.GetValue_Cell(j, 1, 0)[0, 0, 0];
            }

            this.Coordinates[j, 0] = v / (bv * sc);
        }

        /// <summary>
        /// evaluates the mean value of the DG field in one cell;
        /// </summary>
        /// <param name="j">a local cell index</param>
        /// <returns>%</returns>
        virtual public double GetMeanValue(int j) {
            if (j < 0 || j > Basis.GridDat.iLogicalCells.NoOfCells)
                throw new ArgumentException("cell index out of range.", "j");

            int iKref = this.Basis.GridDat.iGeomCells.GetRefElementIndex(j);
            double bv = m_Basis.Polynomials[iKref][0].Coeff[0];
            Debug.Assert(m_Basis.Polynomials[iKref][0].AbsoluteDegree == 0, "Polynomial degree of 0-th polynomial is expected to be 0.");
            Debug.Assert(m_Basis.Polynomials[iKref][0].Coeff.Length == 1, "Polynomial degree of 0-th polynomial is expected to be 0.");

            double sc;
            if (this.Basis.GridDat.iGeomCells.IsCellAffineLinear(j)) {
                sc = Basis.Data.Scaling[j];
            } else {
                sc = Basis.Data.OrthonormalizationTrafo.GetValue_Cell(j, 1, 0)[0, 0, 0];
            }

            return (bv * sc * this.Coordinates[j, 0]);
        }

        /// <summary>
        /// computes the mean value of this DG field.
        /// </summary>
        /// <param name="cm">optional restriction to computational domain</param>
        /// <param name="mean">
        /// If false, the return value equals \f$ \int_{\Omega} u \dV \f$, 
        /// otherwise it equals \f$ \frac{  \int_{\Omega} u \dV }{  \int_{\Omega} 1 \dV  } \f$.
        /// </param>
        virtual public double GetMeanValueTotal(CellMask cm, bool mean = true) {
            using (new FuncTrace()) {
                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                IEnumerable<Chunk> chunks;
                if (cm == null) {
                    Chunk[] _cm = new Chunk[1];
                    _cm[0].i0 = 0;
                    _cm[0].Len = Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                    chunks = _cm;
                } else {
                    chunks = cm;
                }

                double totalVolume = 0;
                double Integral = 0;
                foreach (var c in chunks) { // loop over chunks (of logical cells)...
                    int JE = c.i0 + c.Len;
                    for (int jL = c.i0; jL < JE; jL++) { // loop over logical cells in chunk....
                        double vol = Basis.GridDat.iLogicalCells.GetCellVolume(jL);
                        totalVolume += vol;
                        Integral += GetMeanValue(jL) * vol;
                    }
                }

                unsafe
                {
                    double[] sendBuf = new double[] { Integral, totalVolume };
                    double[] rvcBuf = new double[2];

                    fixed (double* pSend = sendBuf, pRcv = rvcBuf) {
                        csMPI.Raw.Allreduce((IntPtr)pSend, (IntPtr)pRcv, 2, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                    }

                    if (mean == false)
                        rvcBuf[1] = 1;

                    return rvcBuf[0] / rvcBuf[1];
                }
            }
        }

        /// <summary>
        /// computes the integral value of this DG field over the domain.
        /// </summary>
        /// <param name="cm">optional restriction to computational domain</param>
        virtual public double GetIntegral(CellMask cm) {
            using (new FuncTrace()) {
                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                IEnumerable<Chunk> chunks;
                if (cm == null) {
                    Chunk[] _cm = new Chunk[1];
                    _cm[0].i0 = 0;
                    _cm[0].Len = Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                    chunks = _cm;
                } else {
                    chunks = cm;
                }

                double Integral = 0;
                foreach (var c in chunks) {
                    int JE = c.i0 + c.Len;
                    for (int j = c.i0; j < JE; j++) {
                        double vol = Basis.GridDat.iLogicalCells.GetCellVolume(j);
                        Integral += GetMeanValue(j) * vol;
                    }
                }

                unsafe
                {
                    double[] sendBuf = new double[] { Integral };
                    double[] rvcBuf = new double[1];

                    fixed (double* pSend = sendBuf, pRcv = rvcBuf) {
                        csMPI.Raw.Allreduce((IntPtr)pSend, (IntPtr)pRcv, 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                    }

                    return rvcBuf[0];
                }
            }

        }

        /// <summary>
        /// Currently not implemented.
        /// </summary>
        /// <remarks>
        /// This method is not very performance-efficient and shouldn't be
        /// used by performance-critical sections;
        /// </remarks>
        public double GetMeanDerivative(int j, int d) {
            throw new NotImplementedException("really necessary?");

            /*
            int D = m_context.Grid.SpatialDimension;
            if (d < 0 || d >= D)
                throw new ArgumentException("spatial dimension out of range.", "d");

            FullMatrix[] derivMtx = this.m_Basis.DerivativeMatrices;

            int J = m_context.GridDat.Cells.NoOfCells;
            int M = this.m_Basis.GetLength(j);
            double[,] invTrafo = m_context.GridDat.InverseTransformation;

            FullMatrix Qd = new FullMatrix(derivMtx[0].NoOfRows, derivMtx[0].NoOfCols);

            Qd.Clear();
            for (int dd = 0; dd < D; dd++) {
                Qd.Acc(derivMtx[dd], invTrafo[j, d * D + dd]);
            }

            double coord = 0;
            for (int m = 0; m < M; m++) {
                coord += Qd[0, m] * this.Coordinates[j, m];
            }


            double bv = m_Basis.Polynomials[0].Coeff[0];
            double sc = m_context.GridDat.OneOverSqrt_AbsDetTransformation[j];
            return (bv * sc * coord);
            */
        }

        /// <summary>
        /// Currently not implemented.
        /// </summary>
        /// <param name="j"></param>
        /// <param name="gradVec"></param>
        void SetGradient(int j, double[] gradVec) {
            int D = Basis.GridDat.SpatialDimension;
            if (gradVec.Length != D)
                throw new ArgumentException("length must match spatial dimension.", "gradVec");

            double meanVal = GetMeanValue(j);

            throw new NotImplementedException("will come soon.");
        }

        /// <summary>
        /// See <see cref="Evaluate(int,int,NodeSet,MultidimensionalArray,int,double)"/>.
        /// </summary>
        virtual public void Evaluate(int j0, int Len, NodeSet N, MultidimensionalArray result, double ResultPreScale) {
            if (Len > result.GetLength(0))
                throw new ArgumentOutOfRangeException("mismatch between Len and 0-th length of result");
            Evaluate(j0, Len, N, result, 0, ResultPreScale);
        }

        /// <summary>
        /// See <see cref="Evaluate(int,int,NodeSet,MultidimensionalArray,int,double)"/>.
        /// </summary>
        virtual public void Evaluate(int j0, int Len, NodeSet N, MultidimensionalArray result) {
            Evaluate(j0, Len, N, result, 0.0);
        }

        /// <summary>
        /// Evaluates the field;
        /// </summary>
        /// <param name="ResultCellindexOffset">
        /// An offset for the first index of <paramref name="result"/>,
        /// i.e. the first result will be written to
        /// <paramref name="result"/>[<paramref name="ResultCellindexOffset"/>,*].
        /// </param>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="result">
        /// The output: 
        /// On exit, the value of the DG field at the given nodes are
        /// <em>accumulated</em> (!) there. Before the values are added, 
        /// the original content is scaled by <paramref name="ResultPreScale"/>.<br/>
        /// The array is 2-dimensional:
        /// <list type="bullet">
        ///   <item>
        ///   1st index: cell index <i>j</i> - <paramref name="j0"/>;
        ///   </item>
        ///   <item>
        ///   2nd index: node index <i>k</i>, corresponds with 1st index of
        ///   the node set <paramref name="NodeSet"/>;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="ResultPreScale">
        /// Scaling that is applied to <paramref name="result"/> before
        /// the field evaluation is added
        /// </param>
        /// <param name="NodeSet">
        /// nodes to evaluate at
        /// </param>
        /// <remarks>
        /// This method is vectorized: Here, it means that the Points at which
        /// the DG field should be evaluated, are given for one cell in
        /// reference coordinates, but the evaluation is performed for
        /// <paramref name="Len"/> cells at once.
        /// </remarks>
        public abstract void Evaluate(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale);

        /// <summary>
        /// Evaluates the field along edges.
        /// </summary>
        /// <param name="ResultIndexOffset">
        /// An offset for the first index of <paramref name="ValueIN"/> resp.
        /// <paramref name="ValueOUT"/>, i.e. the first result will be written
        /// to
        /// <paramref name="ValueIN"/>[<paramref name="ResultIndexOffset"/>,*].
        /// </param>
        /// <param name="e0">Index of the first edge to evaluate.</param>
        /// <param name="Len">Number of edges to evaluate</param>
        /// <param name="NS">
        /// nodes to evaluate at
        /// </param>
        /// <param name="ValueIN">
        /// If not null, contains the following output: 
        /// On exit, the value of the DG field at the given nodes are
        /// <em>accumulated</em> (!) there. Before the values are added, 
        /// the original content is scaled by <paramref name="ResultPreScale"/>.<br/>
        /// The array is 2-dimensional:
        /// <list type="bullet">
        ///   <item>
        ///   1st index: edge index <i>j</i> - <paramref name="e0"/>;
        ///   </item>
        ///   <item>
        ///   2nd index: node index <i>k</i>, corresponds with 1st index of
        ///   the node set <paramref name="NS"/>;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="ValueOUT">
        /// Same as <paramref name="ValueIN"/>.
        /// </param>
        /// <param name="MeanValueIN">
        /// If not null, contains the following output: 
        /// On exit, the mean values of the DG field at the given edges are
        /// <em>accumulated</em> (!) there. Before the values are added, 
        /// the original content is scaled by <paramref name="ResultPreScale"/>.<br/>
        /// The array is 2-dimensional:
        /// <list type="bullet">
        ///   <item>
        ///   1st index: edge index <i>j</i> - <paramref name="e0"/>;
        ///   </item>
        ///   <item>
        ///   2nd index: node index <i>k</i>, corresponds with 1st index of
        ///   the node set <paramref name="NS"/>;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="MeanValueOT">
        /// Same as <paramref name="MeanValueIN"/>.
        /// </param>
        /// <param name="GradientIN">
        /// If not null, contains the following output: 
        /// On exit, the value of the gradient of DG field at the given nodes
        /// are <em>accumulated</em> (!) there. Before the values are added, 
        /// the original content is scaled by <paramref name="ResultPreScale"/>.<br/>
        /// The array is 3-dimensional:
        /// <list type="bullet">
        ///   <item>
        ///   1st index: edge index <i>j</i> - <paramref name="e0"/>;
        ///   </item>
        ///   <item>
        ///   2nd index: node index <i>k</i>, corresponds with 1st index of
        ///   the node set <paramref name="NS"/>;
        ///   </item>
        ///   <item>
        ///   2rd index: spatial dimension;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="GradientOT">
        /// Same as <paramref name="GradientIN"/>.
        /// </param>
        /// <param name="ResultPreScale">
        /// Scaling that is applied to <paramref name="ValueIN"/> and <paramref name="ValueOUT"/> before
        /// the field evaluation is added
        /// </param>
        public abstract void EvaluateEdge(int e0, int Len, NodeSet NS,
            MultidimensionalArray ValueIN, MultidimensionalArray ValueOUT,
            MultidimensionalArray MeanValueIN, MultidimensionalArray MeanValueOT,
            MultidimensionalArray GradientIN, MultidimensionalArray GradientOT,
            int ResultIndexOffset, double ResultPreScale);


        /// <summary>
        /// Evaluates the gradient of the field;
        /// </summary>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="NodeSet">
        /// as usual, the nodes to evaluate at;
        /// </param>
        /// <param name="result">
        /// on exit, result of the evaluations are accumulated there;
        /// the original content is scaled by <paramref name="ResultPreScale"/>;
        /// 1st index: cell index minus <paramref name="j0"/>;
        /// 2nd index: node index;
        /// 3rd index: spatial coordinate;
        /// </param>
        /// <param name="ResultCellindexOffset">
        /// an offset for the first index of <paramref name="result"/>;
        /// </param>
        /// <param name="ResultPreScale">
        /// see <paramref name="result"/>
        /// </param>
        public abstract void EvaluateGradient(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset = 0, double ResultPreScale = 0.0);

        /// <summary>
        /// evaluates the mean value over a cell;
        /// of course, the mean value doesn't depend on node set or anything like that,
        /// so no information about that has to be provided.
        /// </summary>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="result">
        /// on exit, result of the evaluations are accumulated there;
        /// the original content is scaled by <paramref name="ResultPreScale"/>;
        /// 1st index: cell index minus <paramref name="j0"/>;
        /// </param>
        /// <param name="ResultCellindexOffset">
        /// an offset for the first index of <paramref name="result"/>;
        /// </param>
        /// <param name="ResultPreScale">
        /// see <paramref name="result"/>
        /// </param>
        public abstract void EvaluateMean(int j0, int Len, MultidimensionalArray result, int ResultCellindexOffset = 0, double ResultPreScale = 0.0);

        /// <summary>
        /// evaluates this field at a point given in global coordinates; this
        /// function is very expensive, read remarks!!! This method is
        /// MPI-collective.
        /// </summary>
        /// <param name="_inp">
        /// global (!!!) coordinates
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// This function is very expensive, and it should be only used for
        /// probing (in post-processing or for debugging) at a handful of
        /// points.
        /// </remarks>
        public double ProbeAt(params double[] _inp) {

            // locate cell
            long globalID, globalIndex;
            bool inside, dummy;
            Basis.GridDat.LocatePoint(_inp, out globalID, out globalIndex, out inside, out dummy);
            if (!inside)
                throw new ArgumentException("specified point is outside of grid");

            // process which own the cell
            int ownerProc = Basis.GridDat.CellPartitioning.FindProcess((int)globalIndex);

            double value = double.NaN;
            int D = Basis.GridDat.SpatialDimension;
            if (ownerProc == Basis.GridDat.CellPartitioning.MpiRank) {
                int j = (int)(globalIndex - Basis.GridDat.CellPartitioning.i0);

                MultidimensionalArray inp = MultidimensionalArray.Create(1, D);
                for (int i = 0; i < D; i++) {
                    inp[0, i] = _inp[i];
                }
                MultidimensionalArray outp = MultidimensionalArray.Create(1, 1, D);
                Basis.GridDat.TransformGlobal2Local(inp, outp, j, 1, 0);

                MultidimensionalArray locVtx = outp.ExtractSubArrayShallow(0, -1, -1);
                RefElement Kref = this.GridDat.iGeomCells.GetRefElement(j);
                NodeSet cont = new NodeSet(Kref, locVtx);

                MultidimensionalArray res = MultidimensionalArray.Create(1, 1);
                this.Evaluate(j, 1, cont, res);
                value = res[0, 0];
            }

            // communicate and return;
            value = value.MPIBroadcast(ownerProc);
            return value;
        }

        /// <summary>
        /// adds a vector to the coordinates vector of this field;
        /// Extended Options version;
        /// </summary>
        /// <param name="x">
        /// vector which is added to the coordinates of this field.
        /// </param>
        /// <param name="xSkip">
        /// Skip (offset) into vector <paramref name="x"/>;
        /// </param>
        /// <param name="xStride">
        /// Stride of vector <paramref name="x"/>;
        /// </param>
        /// <param name="OnlyLocalUpdated">
        /// if true, only coordinates of locally updated cells are changed; if
        /// false, also external cells are affected.
        /// </param>
        /// <param name="alpha">
        /// scaling of <paramref name="x"/>
        /// </param>
        internal void _Acc<T>(double alpha, T x, int xSkip, int xStride, bool OnlyLocalUpdated) where T : IList<double> {
            int J;
            if (OnlyLocalUpdated)
                J = Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            else
                J = Basis.GridDat.iLogicalCells.NoOfCells;
            IMatrix Coords = Coordinates;

            if (m_Basis.MaximalLength == m_Basis.MinimalLength) {
                int N = Coordinates.NoOfCols;
                xStride -= N;
                int iSrc = xSkip;
                for (int j = 0; j < J; j++) {
                    for (int n = 0; n < N; n++) {
                        Coords[j, n] += alpha * x[iSrc];
                        iSrc++;
                    }
                    iSrc += xStride;
                }
            } else {
                int N = Coordinates.NoOfCols;
                xStride -= N;
                int iSrc = xSkip;
                for (int j = 0; j < J; j++) {
                    int Nj = m_Basis.GetLength(j);
                    for (int n = 0; n < Nj; n++) {
                        Coords[j, n] += alpha * x[iSrc];
                        iSrc++;
                    }
                    iSrc += xStride;
                }
            }
        }

        /// <summary>
        /// initializes this field to be a (non-shallow) copy the field
        /// <paramref name="other"/>; The basis of the <paramref name="other"/>
        /// field must be contained in the basis of this field (see
        /// <see cref="Basis"/>);
        /// </summary>
        /// <param name="other"></param>
        virtual public void CopyFrom(DGField other) {
            if (!other.Basis.Equals(this.Basis)) {
                throw new ApplicationException(
                    "unable to copy, because the DG polynomial basis of other field is different.");
            }

            int J = this.Coordinates.NoOfRows;
            int N = this.Coordinates.NoOfCols;
            for (int j = 0; j < J; j++)
                for (int n = 0; n < N; n++)
                    this.Coordinates[j, n] = other.Coordinates[j, n];
        }

        /// <summary>
        /// sets all coordinates of this field to 0;
        /// </summary>
        public void Clear() {
            Coordinates.Clear();
        }

        /// <summary>
        /// sets all coordinates of cells in <paramref name="_cellMask"/> of this field to 0;
        /// </summary>
        /// <param name="_cellMask">
        /// optional specification of cell mask; null corresponds to all cells.
        /// </param>
        public void Clear(CellMask _cellMask = null) {

            CellMask cellMask = _cellMask;
            if (cellMask == null) {
                cellMask = CellMask.GetFullMask(this.Basis.GridDat);
            }

            foreach (var chunk in cellMask) {
                for (int j = 0; j < chunk.Len; j++) {

                    int jCell = j + chunk.i0;
                    int N = Basis.GetLength(jCell);
                    for (int n = 0; n < N; n++) {
                        Coordinates[jCell, n] = 0.0;
                    }
                }
            }
        }

        /// <summary>
        /// <see cref="Basis"/>
        /// </summary>
        protected Basis m_Basis;

        /// <summary>
        /// identification string
        /// </summary>
        protected String m_Identification;

        [NonSerialized]
        CoordinateMapping m_Mapping;

        /// <summary>
        /// creates a new <see cref="CoordinateMapping"/> which contains only this field
        /// </summary>
        public CoordinateMapping Mapping {
            get {
                if (m_Mapping == null) {
                    m_Mapping = new CoordinateMapping(this);
                }
                return m_Mapping;
            }
        }

        [NonSerialized]
        CoordinateVector m_CoordinatesAsVector;

        /// <summary>
        /// creates a new <see cref="CoordinateVector"/> which contains the DG
        /// coordinates of this field as one column vector
        /// </summary>
        public CoordinateVector CoordinateVector {
            get {
                if (m_CoordinatesAsVector == null) {
                    m_CoordinatesAsVector = new CoordinateVector(Mapping);
                }
                return m_CoordinatesAsVector;
            }
        }

        /// <summary>
        /// DG coordinates of this field;
        /// </summary>
        /// <remarks>
        /// Indices:
        /// <list type="bullet">
        ///   <item>1st index/row index: local cell index;</item>
        ///   <item>
        ///         2nd index/column index: basis function index, corresponds with the order 
        ///         of polynomials
        ///         in the associated <see cref="Basis"/>, <see cref="Basis"/>.
        ///   </item>
        /// </list>
        /// 
        /// </remarks>
        abstract public IMatrix Coordinates {
            get;
        }

        /// <summary>
        /// User given string identification for this field;
        /// </summary>
        public string Identification {
            get {
                return m_Identification;
            }
            set {
                m_Identification = value;
            }
        }

        /// <summary>
        /// the Basis assigned to this field
        /// </summary>
        virtual public Basis Basis {
            get {
                return m_Basis;
            }
        }

        /// <summary>
        /// the gird on which this DG field is defined
        /// </summary>
        public IGridData GridDat {
            get {
                return Basis.GridDat;
            }
        }


        /// <summary>
        /// returns the identification string of this field
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return Identification;
        }

        /// <summary>
        /// checks whether NAS's or INF's are present in this field
        /// </summary>
        /// <param name="CheckForInf"></param>
        /// <param name="CheckForNan"></param>
        /// <param name="ExceptionIfFound">
        /// if this is true, an ArithmeticException is thrown in case of a
        /// positive result
        /// </param>
        /// <returns>
        /// The GlobalId of the first cell (depending on the actual
        /// GlobalID-Permutation <see cref="Element.GlobalID"/>)
        /// within the actual MPI-process in which an 'illegal value'
        /// (depending on <paramref name="CheckForInf"/> and
        /// <paramref name="CheckForNan"/>) is found; otherwise, a negative
        /// number;
        /// </returns>
        public long CheckForNanOrInf(
            bool CheckForInf = true,
            bool CheckForNan = true,
            bool ExceptionIfFound = true) {

            IMatrix coord = Coordinates;
            int J = Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int N = Coordinates.NoOfCols;

            long ret = -1;
            for (int j = 0; j < J; j++) {
                for (int n = 0; n < N; n++) {

                    double u_jn = coord[j, n];

                    if (CheckForNan)
                        if (double.IsNaN(u_jn)) {
                            ret = Basis.GridDat.iLogicalCells.GetGlobalID(j);
                            if (ExceptionIfFound)
                                throw new ArithmeticException("NaN found in field '" + Identification
                                    + "' in Cell (GlobalID = " + ret + ", local index = " + j
                                    + "), DG coordinate index: " + n + ";");
                        }

                    if (CheckForInf)
                        if (double.IsInfinity(u_jn)) {
                            ret = Basis.GridDat.iLogicalCells.GetGlobalID(j);
                            if (ExceptionIfFound)
                                throw new ArithmeticException("Inf found in field '" + Identification
                                    + "' in Cell (GlobalID = " + ret + ", local index = " + j
                                    + "), DG coordinate index: " + n + ";");
                        }

                    if (ret >= 0)
                        return ret;
                }
            }

            return -1;
        }

        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// See <see cref="ICloneable.Clone"/>
        /// </summary>
        abstract public object Clone();

        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// See <see cref="ICloneable.Clone"/>
        /// </summary>
        /// <returns>A new object that is a copy of this instance.</returns>
        public DGField CloneAs() {
            return (DGField)Clone();
        }

        /// <summary>
        /// For MPI communication; <br/>
        /// Determines the size of the receive buffer that some <see cref="Comm.Transceiver"/>
        /// must allocate;
        /// </summary>
        /// <param name="proc"></param>
        /// <returns>
        /// number of double's in the send buffer for processor <paramref name="proc"/>;
        /// </returns>
        /// <remarks>
        /// Of course, the size of the send buffer on this process
        /// must be equal to the size of the receive buffer (<see cref="GetMPIRecvBufferSize"/>)
        /// on process <paramref name="proc"/>. However, this is not tested, and a 
        /// violation of this rule will possibly result in an MPI error.
        /// </remarks>
        public abstract int GetMPISendBufferSize(int proc);

        /// <summary>
        /// For MPI communication; <br/>
        /// Determines the size of the receive buffer that some <see cref="Comm.Transceiver"/>
        /// must allocate;
        /// </summary>
        /// <param name="proc">
        /// MPI processor rank.
        /// </param>
        /// <returns>
        /// number of double's in the receive buffer for processor <paramref name="proc"/>;
        /// </returns>
        /// <remarks>
        /// Of course, the size of the receive buffer on this process
        /// must be equal to the size of the send buffer (<see cref="GetMPISendBufferSize"/>)
        /// on process <paramref name="proc"/>. However, this is not tested, and a 
        /// violation of this rule will possibly result in an MPI error.
        /// </remarks>
        public abstract int GetMPIRecvBufferSize(int proc);

        /// <summary>
        /// Copies the DG coordinates which must be send to process
        /// <paramref name="proc"/> to the send buffer
        /// (<paramref name="Buffer"/>); This coordinates are those of the
        /// cells (with local index)
        /// <see cref="GridData.Parallelization.m_SendCommLists"/>[<paramref name="proc"/>];
        /// </summary>
        /// <param name="Buffer">
        /// send buffer; the first item to send must be inserted at
        /// index <paramref name="i0"/>;
        /// </param>
        /// <param name="i0">
        /// index offset into to the send buffer;
        /// </param>
        /// <param name="proc">
        /// MPI processor rank of the process to send data to.
        /// </param>
        /// <returns>
        /// The number of double's that were copied to <see cref="Buffer"/>,
        /// equal to the return value of
        /// <see cref="GetMPISendBufferSize"/>(<paramref name="proc"/>,...).
        /// </returns>
        public abstract int FillMPISendBuffer(int proc, double[] Buffer, int i0);

        /// <summary>
        /// Copies the DG coordinates which were received from processor <paramref name="proc"/>
        /// to the internal storage of this DG field.
        /// </summary>
        /// <param name="proc">
        ///  MPI processor rank of the process from which data was received;
        /// </param>
        /// <param name="Buffer">
        /// receive buffer; the first received item, which belongs to this DG field is 
        /// located at index <paramref name="i0"/>;
        /// </param>
        /// <param name="i0">
        /// index offset into the <paramref name="Buffer"/>
        /// </param>
        /// <returns>
        /// The number of double's that were copied from <see cref="Buffer"/>,
        /// equal to the return value of <see cref="GetMPIRecvBufferSize"/>(<paramref name="proc"/>,...).
        /// </returns>
        public abstract int CopyFromMPIrecvBuffer(int proc, double[] Buffer, int i0);

        /// <summary>
        /// MPI update of this field
        /// </summary>
        public void MPIExchange() {
            var trx = new Comm.Transceiver(this);
            trx.TransceiveStartImReturn();
            trx.TransceiveFinish();
        }
    }
}
