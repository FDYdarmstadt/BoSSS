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
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Tracing;
using MPI.Wrappers;
using ilPSP;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// a level set that is implemented as a DG field, i.e. a <see cref="XDGField"/>
    /// </summary>
    public partial class LevelSet : SinglePhaseField, ILevelSet {

        /// <summary>
        /// Constructor
        /// </summary>
        public LevelSet(Basis b, string id)
            : base(b, id) {
            if (b.Degree <= 0) {
                throw new ArgumentException(
                    "Level Set Polynomial degree of 0 is not suitable:"
                    + " the level set tracker detects cut cells if the sign of"
                    + " the level set changes within one cell -- impossible if"
                    + " polynomial degree is 0, i.e. the level set is"
                    + " constant within one cell", "b");
            }
        }

        #region ILevelSet Members

        /// <summary>
        /// as defined by <see cref="ILevelSet.EvaluateGradient"/>
        /// </summary>
        public virtual void EvaluateGradient(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result) {
            base.EvaluateGradient(j0, Len, NodeSet, result, 0, 0.0);
        }

        

        #endregion


        /// <summary>
        /// Evaluates the broken Laplacian (by cell-local analytic derivation of
        /// the basis polynomials) of this field; this method my move to the
        /// <see cref="DGField"/>-class in future;
        /// </summary>
        /// <param name="j0"></param>
        /// <param name="Len"></param>
        /// <param name="NodeSet"></param>
        /// <param name="result">
        /// <list type="bullet">
        ///   <item>1st index: cell index <em>j</em></item>
        ///   <item>2nd index: node index <em>m</em> into nodeset #<paramref name="NodeSetIndex"/></item>
        /// </list>
        /// So, the entry [j,m] = \f$ \sum_{d=1}^{D} \frac{\partial}{\partial x_d} \varphi (\vec{\xi}_m)\f$ 
        /// where \f$ \vec{xi}_m\f$  is the <em>m</em>-th vector in the nodeset #<paramref name="NodeSetIndex"/>,
        /// in the <em>j</em>-th cell.
        /// </param>
        /// <remarks>
        /// Because of 2 derivatives taken, this field needs to be at least of DG degree 2 to get a non-zero result
        /// from this method.
        /// </remarks>
        public void EvaluateLaplacian(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result) {
            int D = this.GridDat.SpatialDimension;
            int K = NodeSet.NoOfNodes; // No of nodes

            if (result.Dimension != 2)
                throw new ArgumentException();
            if (result.GetLength(0) != Len)
                throw new ArgumentException();
            if (result.GetLength(1) != K)
                throw new ArgumentException();


            MultidimensionalArray Hess = MultidimensionalArray.Create(Len, K, D, D);
            this.EvaluateHessian(j0, Len, NodeSet, Hess);


            for (int i = 0; i < Len; i++) {
                for (int k = 0; k < K; k++) {
                    double acc = 0;
                    for (int d = 0; d < D; d++) {
                        acc += Hess[i, k, d, d];
                    }
                    result[i, k] = acc;
                }
            }

        }


        /// <summary>
        /// Assigns the normalized gradient of the level set to the Output
        /// vector
        /// </summary>
        /// <param name="Output">Normal vector</param>
        /// <param name="optionalSubGrid">
        /// Restriction of the computations to a an optional subgrid
        /// </param>
        /// <param name="bndMode"></param>
        public void ComputeNormalByFlux(VectorField<SinglePhaseField> Output, SubGrid optionalSubGrid = null, SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {
            if (this.m_Basis.Degree < 1)
                throw new ArgumentException("For correct computation of these level set quantities, the level set has to be at least of degree 1!");
            SinglePhaseField absval = new SinglePhaseField(Output[0].Basis);
            //Output.Clear();
            for (int i = 0; i < Output.Dim; i++) {
                Output[i].DerivativeByFlux(1.0, this, i, optionalSubGrid, bndMode);
            }

            absval.ProjectAbs(1.0, Output);
            for (int i = 0; i < Output.Dim; i++) {
                Output[i].ProjectQuotient(1.0, Output[i], absval, null, false);
            }
        }

        /// <summary>
        /// Determines the total curvature which is twice the mean curvature
        /// Please pay attention to the sign of this expression!
        /// </summary>
        /// <param name="Output">The total curvature</param>
        /// <param name="optionalSubGrid">
        /// Subgrid which can be defined, for example for carrying out
        /// computations on a narrow band
        /// </param>
        /// <param name="bndMode">
        /// Definition of the behavior at subgrid boundaries
        /// </param>
        public void ComputeTotalCurvatureByFlux(SinglePhaseField Output,
            SubGrid optionalSubGrid = null,
            SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {

            if (this.m_Basis.Degree <= 1)
                throw new ArgumentException("For correct computation of these level set quantities, the level set has to be at least of degree 2!");

            Basis basisForNormalVec = new Basis(this.GridDat, this.m_Basis.Degree - 1);
            Basis basisForCurvature = new Basis(this.GridDat, this.m_Basis.Degree - 2);

            Func<Basis,string,SinglePhaseField> fac = (Basis b, string id) => new SinglePhaseField(b, id);
            VectorField<SinglePhaseField> normalVector = new VectorField<SinglePhaseField>(this.GridDat.SpatialDimension, basisForNormalVec, fac);
            ComputeNormalByFlux(normalVector, optionalSubGrid, bndMode);
            VectorField<SinglePhaseField> secondDerivatives = new VectorField<SinglePhaseField>(this.GridDat.SpatialDimension, basisForCurvature, fac);

            Output.Clear();
            for (int i = 0; i < normalVector.Dim; i++) {
                secondDerivatives[i].DerivativeByFlux(1.0, normalVector[i], i, optionalSubGrid, bndMode);

                Output.Acc(-1.0, secondDerivatives[i], optionalSubGrid.VolumeMask);
            }

        }

        /// <summary>
        /// Vector of the d partial derivatives of the i-th component of the
        /// normal vector 
        /// </summary>
        /// <param name="Output">Vector of the partial derivatives</param>
        /// <param name="componentOfNormalVec">
        /// specifies the component i of the normal vector
        /// </param>
        /// <param name="optionalSubGrid"></param>
        /// <param name="bndMode"></param>
        public void ComputeDerivativesOfTheNormalByFlux(
            VectorField<SinglePhaseField> Output,
            int componentOfNormalVec,
            SubGrid optionalSubGrid = null,
            SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {

            if (this.m_Basis.Degree <= 1)
                throw new ArgumentException("For correct computation of these level set quantities, the level set has to be at least of degree 2!");
            Basis basisForNormalVec = new Basis(this.GridDat, this.m_Basis.Degree - 1);

            Func<Basis,string,SinglePhaseField> fac = (Basis b, string id) => new SinglePhaseField(b, id);
            VectorField<SinglePhaseField> normalVector = new VectorField<SinglePhaseField>(this.GridDat.SpatialDimension, basisForNormalVec, fac);
            ComputeNormalByFlux(normalVector, optionalSubGrid, bndMode);

            for (int i = 0; i < Output.Dim; i++) {
                Output[i].DerivativeByFlux(1.0, normalVector[componentOfNormalVec], i, optionalSubGrid, bndMode);
            }
        }

        /// <summary>
        /// Assigns the normalized gradient of the level set to the Output
        /// vector
        /// </summary>
        /// <param name="Output">Normal Vector</param>
        /// <param name="optionalCellMask">
        /// Cell mask used when computing the derivatives
        /// </param>
        public void ComputeNormal(VectorField<SinglePhaseField> Output, CellMask optionalCellMask) {
            if (this.m_Basis.Degree < 1)
                throw new ArgumentException("For correct computation of these level set quantities, the level set has to be at least of degree 1!");
            SinglePhaseField absval = new SinglePhaseField(Output[0].Basis);
            Output.Clear();
            for (int i = 0; i < Output.Dim; i++) {
                Output[i].Derivative(1.0, this, i, optionalCellMask);
            }

            absval.ProjectAbs(1.0, Output);

            for (int i = 0; i < Output.Dim; i++) {
                Output[i].ProjectQuotient(1.0, Output[i], absval, null, false);
            }
        }

        /// <summary>
        /// Computation of mean curvature according to Bonnet's formula.
        /// </summary>
        /// <param name="scale"></param>
        /// <param name="Output">
        /// output.
        /// </param>
        /// <param name="quadScheme"></param>
        /// <param name="UseCenDiffUpTo">
        /// Either 0, 1, or 2:
        /// If 0, all derivatives are computed locally (broken derivative); 
        /// if 1, the first order derivatives are computed by central
        /// differences, while the second order ones are computed locally,
        /// based on the first order ones;
        /// if 2, all derivatives are computed by central differences.
        /// </param>
        /// <param name="_1stDerivDegree">
        /// Relative DG polynomial degree for the 1st order derivatives, i.e.
        /// degree is <paramref name="_1stDerivDegree"/>+<em>p</em>, where
        /// <em>p</em> is the degree of this field.
        /// Only active if <paramref name="UseCenDiffUpTo"/> is greater than 0.
        /// </param>
        /// <param name="_2ndDerivDegree">
        /// Relative DG polynomial degree for the 2nd order derivatives, i.e.
        /// degree is <paramref name="_2ndDerivDegree"/>+<em>p</em>, where
        /// <em>p</em> is the degree of this field. Only active if
        /// <paramref name="UseCenDiffUpTo"/> is greater than 1.
        /// </param>
        /// <remarks>
        /// using central differences causes memory allocation: <em>D</em>
        /// fields for <paramref name="UseCenDiffUpTo"/>=1, and
        /// <em>D</em>*(<em>D</em>+1) for <paramref name="UseCenDiffUpTo"/>=2,
        /// where <em>D</em> notates the spatial dimension.
        /// </remarks>
        public void ProjectTotalcurvature2(
            double scale,
            SinglePhaseField Output,
            int UseCenDiffUpTo,
            int _1stDerivDegree = 0,
            int _2ndDerivDegree = 0,
            CellQuadratureScheme quadScheme = null) {

            using (new FuncTrace()) {
                if (UseCenDiffUpTo < 0 || UseCenDiffUpTo > 2)
                    throw new ArgumentOutOfRangeException();

                //int M = Output.Basis.Length;
                //int N = this.Basis.Length;
                int D = this.GridDat.SpatialDimension;
                //var NSC = m_context.NSC;

                SubGrid sgrd = null;
                SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.InnerEdge;
                if (UseCenDiffUpTo >= 1 && quadScheme != null && quadScheme.Domain != null) {
                    sgrd = new SubGrid(quadScheme.Domain);
                    bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary;
                }

                // compute 1st order derivatives by central differences, if desired
                // ================================================================
                Basis B2 = new Basis(this.GridDat, this.Basis.Degree + _1stDerivDegree);

                SinglePhaseField[] GradientVector = null;
                if (UseCenDiffUpTo >= 1) {
                    GradientVector = new SinglePhaseField[D];
                    for (int d = 0; d < D; d++) {
                        GradientVector[d] = new SinglePhaseField(B2);
                        GradientVector[d].DerivativeByFlux(1.0, this, d, sgrd, bndMode);
                    }
                }

                // compute 2nd order derivatives by central differences, if desired
                // ===============================================================
                Basis B3 = new Basis(this.GridDat, this.Basis.Degree + _2ndDerivDegree);

                SinglePhaseField[,] HessianTensor = null;
                if (UseCenDiffUpTo >= 2) {
                    HessianTensor = new SinglePhaseField[D, D];
                    for (int d1 = 0; d1 < D; d1++) {
                        for (int d2 = 0; d2 < D; d2++) {
                            HessianTensor[d1, d2] = new SinglePhaseField(B3);
                            HessianTensor[d1, d2].DerivativeByFlux(1.0, GradientVector[d1], d2, sgrd, bndMode);
                        }
                    }
                }

                // compute and project 
                // ===================

                // buffers:
                MultidimensionalArray Phi = new MultidimensionalArray(2);
                MultidimensionalArray GradPhi = new MultidimensionalArray(3);
                MultidimensionalArray HessPhi = new MultidimensionalArray(4);

                MultidimensionalArray ooNormGrad = new MultidimensionalArray(2);
                MultidimensionalArray Laplace = new MultidimensionalArray(2);
                MultidimensionalArray Q = new MultidimensionalArray(3);

                // evaluate/project:
                //double Erracc = 0;
                Output.ProjectField(scale,
                    (ScalarFunctionEx)delegate(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes


                    // alloc buffers
                    // -------------

                    if (Phi.GetLength(0) != Len || Phi.GetLength(1) != K) {
                        Phi.Allocate(Len, K);
                        GradPhi.Allocate(Len, K, D);
                        HessPhi.Allocate(Len, K, D, D);
                        ooNormGrad.Allocate(Len, K);
                        Laplace.Allocate(Len, K);
                        Q.Allocate(Len, K, D);
                    } else {
                        Phi.Clear();
                        GradPhi.Clear();
                        HessPhi.Clear();
                        ooNormGrad.Clear();
                        Laplace.Clear();
                        Q.Clear();
                    }

                    // evaluate Gradient and Hessian
                    // -----------------------------

                    if (UseCenDiffUpTo >= 1) {
                        for (int d = 0; d < D; d++)
                            GradientVector[d].Evaluate(j0, Len, NodeSet, GradPhi.ExtractSubArrayShallow(-1, -1, d));
                    } else {
                        this.EvaluateGradient(j0, Len, NodeSet, GradPhi);
                    }

                    if (UseCenDiffUpTo == 2) {
                        for (int d1 = 0; d1 < D; d1++)
                            for (int d2 = 0; d2 < D; d2++)
                                HessianTensor[d1, d2].Evaluate(j0, Len, NodeSet, HessPhi.ExtractSubArrayShallow(-1, -1, d1, d2));
                    } else if (UseCenDiffUpTo == 1) {
                        for (int d = 0; d < D; d++) {
                            var GradientVector_d = GradientVector[d];
                            GradientVector_d.EvaluateGradient(j0, Len, NodeSet, HessPhi.ExtractSubArrayShallow(-1, -1, d, -1), 0, 0.0);
                        }
                    } else if (UseCenDiffUpTo == 0) {
                        this.EvaluateHessian(j0, Len, NodeSet, HessPhi);
                    } else
                        Debug.Assert(false);

                    // compute the monstrous formula
                    // -----------------------------

                    // norm of Gradient:
                    for (int d = 0; d < D; d++) {
                        var GradPhi_d = GradPhi.ExtractSubArrayShallow(-1, -1, d);
                        ooNormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
                    }
                    ooNormGrad.ApplyAll(x => 1.0 / Math.Sqrt(x));

                    // laplacian of phi:
                    for (int d = 0; d < D; d++) {
                        var HessPhi_d_d = HessPhi.ExtractSubArrayShallow(-1, -1, d, d);
                        Laplace.Acc(1.0, HessPhi_d_d);
                    }

                    // result = Laplacian(phi)/|Grad phi|
                    result.Multiply(1.0, Laplace, ooNormGrad, 0.0, "ik", "ik", "ik");


                    // result = Grad(1/|Grad(phi)|)
                    for (int d1 = 0; d1 < D; d1++) {
                        var Qd = Q.ExtractSubArrayShallow(-1, -1, d1);

                        for (int d2 = 0; d2 < D; d2++) {
                            var Grad_d2 = GradPhi.ExtractSubArrayShallow(-1, -1, d2);
                            var Hess_d2_d1 = HessPhi.ExtractSubArrayShallow(-1, -1, d2, d1);

                            Qd.Multiply(-1.0, Grad_d2, Hess_d2_d1, 1.0, "ik", "ik", "ik");
                        }
                    }

                    ooNormGrad.ApplyAll(x => x * x * x);

                    result.Multiply(1.0, GradPhi, Q, ooNormGrad, 1.0, "ik", "ikd", "ikd", "ik");

                    //for (int i = 0; i < Len; i++) {
                    //    for (int k = 0; k < K; k++) {
                    //        double acc = 0;
                    //        for (int d = 0; d < D; d++) {
                    //            acc += GradPhi[i,k,d]*Q[i,k,d]*ooNormGrad[i,k];
                    //        }

                    //        Erracc += (acc - result[i,k]).Abs();
                    //    }
                    //}

                },
                    quadScheme.SaveCompile(this.GridDat, (Output.Basis.Degree + this.m_Basis.Degree * (this.m_Basis.Degree - 1) * D) * 2)
                    );
            }
        }

        /// <summary>
        /// curvature computation according to Bonnet's formula
        /// </summary>
        public void EvaluateTotalCurvature(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result) {
            // (siehe FK, persoenliche Notizen, 08mar13)


            // checks
            // ------

            int K = NodeSet.NoOfNodes;
            if (result.Dimension != 2)
                throw new ArgumentException();
            if (result.GetLength(0) != Len)
                throw new ArgumentException();
            if (result.GetLength(1) != K)
                throw new ArgumentException();
            int D = NodeSet.SpatialDimension;
            Debug.Assert(D == this.GridDat.SpatialDimension);


            // buffers:
            // --------
            //MultidimensionalArray Phi = new MultidimensionalArray(2);
            MultidimensionalArray GradPhi = new MultidimensionalArray(3);
            MultidimensionalArray HessPhi = new MultidimensionalArray(4);

            MultidimensionalArray ooNormGrad = new MultidimensionalArray(2);
            MultidimensionalArray Laplace = new MultidimensionalArray(2);
            MultidimensionalArray Q = new MultidimensionalArray(3);

            //Phi.Allocate(Len, K);
            GradPhi.Allocate(Len, K, D);
            HessPhi.Allocate(Len, K, D, D);
            ooNormGrad.Allocate(Len, K);
            Laplace.Allocate(Len, K);
            Q.Allocate(Len, K, D);

            // derivatives
            // -----------

            // evaluate gradient
            this.EvaluateGradient(j0, Len, NodeSet, GradPhi);

            // evaluate Hessian
            this.EvaluateHessian(j0, Len, NodeSet, HessPhi);


            // compute the monstrous formula
            // -----------------------------

            // norm of Gradient:
            for (int d = 0; d < D; d++) {
                var GradPhi_d = GradPhi.ExtractSubArrayShallow(-1, -1, d);
                ooNormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
            }
            ooNormGrad.ApplyAll(x => 1.0 / Math.Sqrt(x));

            // laplacian of phi:
            for (int d = 0; d < D; d++) {
                var HessPhi_d_d = HessPhi.ExtractSubArrayShallow(-1, -1, d, d);
                Laplace.Acc(1.0, HessPhi_d_d);
            }

            // result = Laplacian(phi)/|Grad phi|
            result.Multiply(1.0, Laplace, ooNormGrad, 0.0, "ik", "ik", "ik");


            // result += Grad(1/|Grad(phi)|)
            for (int d1 = 0; d1 < D; d1++) {
                var Qd = Q.ExtractSubArrayShallow(-1, -1, d1);

                for (int d2 = 0; d2 < D; d2++) {
                    var Grad_d2 = GradPhi.ExtractSubArrayShallow(-1, -1, d2);
                    var Hess_d2_d1 = HessPhi.ExtractSubArrayShallow(-1, -1, d2, d1);

                    Qd.Multiply(-1.0, Grad_d2, Hess_d2_d1, 1.0, "ik", "ik", "ik");
                }
            }

            ooNormGrad.ApplyAll(x => x * x * x);

            result.Multiply(1.0, GradPhi, Q, ooNormGrad, 1.0, "ik", "ikd", "ikd", "ik");

        }


        /// <summary>
        /// Determines the total curvature which is twice the mean curvature
        /// Please pay attention to the sign of this expression!
        /// </summary>
        /// <param name="Output">Field representing the total curvature</param>
        /// <param name="optionalCellMask">
        /// Cell Mask for the derivatives and their accumulation regarding the
        /// curvature
        /// </param>
        public void ComputeTotalCurvature(SinglePhaseField Output, CellMask optionalCellMask) {
            if (this.m_Basis.Degree <= 1)
                throw new ArgumentException("For correct computation of these level set quantities, the level set has to be at least of degree 2!");

            Basis basisForNormalVec = new Basis(this.GridDat, this.m_Basis.Degree - 1);
            Basis basisForCurvature = new Basis(this.GridDat, this.m_Basis.Degree - 2);

            Func<Basis,string,SinglePhaseField> fac = (Basis b, string id) => new SinglePhaseField(b, id);
            VectorField<SinglePhaseField> normalVector = new VectorField<SinglePhaseField>(this.GridDat.SpatialDimension, basisForNormalVec, fac);
            ComputeNormal(normalVector, optionalCellMask);
            VectorField<SinglePhaseField> secondDerivatives = new VectorField<SinglePhaseField>(this.GridDat.SpatialDimension, basisForCurvature, fac);
            Output.Clear();
            for (int i = 0; i < normalVector.Dim; i++) {
                secondDerivatives[i].Derivative(1.0, normalVector[i], i, optionalCellMask);
                Output.Acc(-1.0, secondDerivatives[i], optionalCellMask);
            }
        }

        /// <summary>
        /// Vector of the d partial derivatives of the i-th component of the
        /// normal vector 
        /// </summary>
        /// <param name="Output">Vector of the partial derivatives</param>
        /// <param name="componentOfNormalVec">
        /// specifies the component i of the normal vector
        /// </param>
        /// <param name="optionalCellMask">
        /// Cell Mask for the computation
        /// </param>
        public void ComputeDerivativesOfTheNormal(
            VectorField<SinglePhaseField> Output, int componentOfNormalVec, CellMask optionalCellMask) {

            if (this.m_Basis.Degree <= 1)
                throw new ArgumentException("For correct computation of these level set quantities, the level set has to be at least of degree 2!");
            Basis basisForNormalVec = new Basis(this.GridDat, this.m_Basis.Degree - 1);

            Func<Basis,string,SinglePhaseField> fac = (Basis b, string id) => new SinglePhaseField(b, id);
            VectorField<SinglePhaseField> normalVector = new VectorField<SinglePhaseField>(this.GridDat.SpatialDimension, basisForNormalVec, fac);
            ComputeNormal(normalVector, optionalCellMask);

            for (int i = 0; i < Output.Dim; i++) {
                Output[i].Derivative(1.0, normalVector[componentOfNormalVec], i, optionalCellMask);
            }

        }

        /// <summary>
        /// Approximation of the signum function of the level set field
        /// </summary>
        /// <param name="Output"></param>
        /// <param name="gridParameter"></param>
        /// <param name="optionalCellMask"></param>
        public void ApproximateSignFunction(
            SinglePhaseField Output, double gridParameter, CellMask optionalCellMask) {

            SinglePhaseField quot = new SinglePhaseField(m_Basis, "quot");
            SinglePhaseField sqrt = new SinglePhaseField(m_Basis, "sqrt");

            sqrt.ProjectProduct(1.0, this, this, optionalCellMask);
            sqrt.AccConstant(gridParameter * gridParameter, optionalCellMask);

            quot.ProjectPow(1.0, sqrt, 0.5, optionalCellMask);
            //Output.ProjectQuotient (1.0, this, quot,false);
            Output.ProjectQuotient(1.0, this, quot, optionalCellMask, false);

        }

        /// <summary>
        /// computes the L2 norm of the Residual of the Eikonal equation in the
        /// domain <paramref name="K"/>, i.e. <br/>
        /// \f[ 
        /// \left\|
        /// \textrm{\textbf{1}}_K
        /// \cdot
        /// \left(  | \nabla \phi | - 1 \right)
        /// \right\|_2
        /// \f]
        /// </summary>
        /// <param name="K">
        /// optional restriction to a subdomain
        /// </param>
        public double SignedDistanceError(CellMask K) {
            var cqc = new CellQuadratureScheme(true, K);

            var q = new SignedDistanceErrorQuad(
                this, cqc.Compile(this.Basis.GridDat, this.Basis.Degree * 4));

            q.Execute();
            unsafe {
                double local = q.LocalErrorL2Accumulated;
                double global;
                csMPI.Raw.Allreduce((IntPtr)(&local), (IntPtr)(&global), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                return global;
            }
        }

        /// <summary>
        /// Quadrature for the computation of the deviance of the level set
        /// function from a signed distance function.
        /// </summary>
        class SignedDistanceErrorQuad : Foundation.Quadrature.CellQuadrature {

            /// <summary>
            /// ctor
            /// </summary>
            public SignedDistanceErrorQuad(LevelSet phi, ICompositeQuadRule<QuadRule> daRule)
                : base(new int[] { 1 }, phi.GridDat, daRule) {
                m_phi = phi;
                m_gradPhi = new MultidimensionalArray(3);
            }

            /// <summary>
            /// owner
            /// </summary>
            private LevelSet m_phi;

            /// <summary>
            /// Values of the gradient of <see cref="m_phi"/>
            /// </summary>
            private MultidimensionalArray m_gradPhi;

            /// <summary>
            /// Sum of L2 errors in the current process
            /// </summary>
            internal double LocalErrorL2Accumulated = 0;

           

            /// <summary>
            /// <see cref="Quadrature{U,V}.CreateNodeSetFamily"/>
            /// </summary>
            /// <param name="NoOfItems"></param>
            /// <param name="rule"></param>
            protected override void AllocateBuffers(int NoOfItems, NodeSet rule) {
                base.AllocateBuffers(NoOfItems, rule);
                int NoOfNodes = rule.GetLength(0);

                int D = m_phi.GridDat.SpatialDimension;
                m_gradPhi.Allocate(NoOfItems, NoOfNodes, D);
            }

            /// <summary>
            /// Computes the residual of the Eikonal equation.
            /// </summary>
            protected override void Evaluate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = QR.Nodes;
                var gradPhi = m_gradPhi;
                int D = m_phi.GridDat.SpatialDimension;
                int NoOfNodes = QuadNodes.NoOfNodes;

                m_phi.EvaluateGradient(i0, Length, QuadNodes, m_gradPhi);

                for (int i = 0; i < Length; i++) {
                    for (int n = 0; n < NoOfNodes; n++) {
                        // compute the residual of the Eikonal equation:
                        double acc = 0;

                        for (int d = 0; d < D; d++) {
                            double dPhi_dxd = gradPhi[i, n, d];
                            acc += dPhi_dxd * dPhi_dxd;
                        }

                        acc = Math.Sqrt(acc);
                        acc -= 1;

                        // L2-Norm of the residual:
                        EvalResult[i, n, 0] = acc * acc;
                    }
                }
            }

            /// <summary>
            /// Accumulates the result to <see cref="LocalErrorL2Accumulated"/>
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="ResultsOfIntegration"></param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {

                // sum over all cells
                double acc = 0;
                for (int i = 0; i < Length; i++) {
                    // in Cell i0 + i, the L2 norm of the residual of the Eikonal equation
                    acc += ResultsOfIntegration[i, 0];
                }

                LocalErrorL2Accumulated += acc;
            }
        }

        /// <summary>
        /// Clones this object
        /// </summary>
        new public LevelSet CloneAs() {
            return (LevelSet)Clone();
        }

        /// <summary>
        /// Clones this object
        /// </summary>
        public override object Clone() {
            var ret = new LevelSet(this.m_Basis, this.Identification);

            MultidimensionalArray dst = ((MultidimensionalArray)ret.Coordinates);
            MultidimensionalArray src = ((MultidimensionalArray)base.Coordinates);

            if (dst.IsContinious && src.IsContinious) {
                Array.Copy(src.Storage, dst.Storage, dst.Length);
            } else {
                throw new NotImplementedException("remainder");
            }

            return ret;
        }

        //TO BE CHANGED!!!!!
        public void GradientNorm(out SinglePhaseField norm, CellMask cm) {
            Basis basisForNorm = new Basis(this.GridDat, this.Basis.Degree * 2);
            norm = new SinglePhaseField(basisForNorm);

            Func<Basis,string,SinglePhaseField> fac = (Basis b, string id) => new SinglePhaseField(b, id);
            
            VectorField<SinglePhaseField> Gradient = new VectorField<SinglePhaseField>(GridDat.SpatialDimension, Basis, fac);
            for (int i = 0; i < GridDat.SpatialDimension; i++) {
                Gradient[i].Derivative(1.0, this, i, cm);
            }
            norm.ProjectAbs<SinglePhaseField>(1.0, Gradient);
        }

    }
}
