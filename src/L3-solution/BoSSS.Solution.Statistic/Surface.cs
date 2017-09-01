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
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Statistic.QuadRules {

    /// <summary>
    /// 
    /// </summary>
    public abstract class Surface : ISurfaceEvaluation {

        /// <summary>
        /// The context object
        /// </summary>
        private GridData m_Context;
        /// <summary>
        /// Basis of the field representation
        /// </summary>
        private Basis m_Basis;
        /// <summary>
        /// Basis of the gradient of the level set function
        /// </summary>
        private Basis m_gradBasis;
        /// <summary>
        /// Basis of the second derivatives of the level set function
        /// </summary>
        private Basis m_derivsBasis;
        /// <summary>
        /// Field representation of the surface
        /// </summary>
        protected SinglePhaseField m_Field;
        /// <summary>
        /// Gradient of the level set representation
        /// </summary>
        private VectorField<SinglePhaseField> m_Gradw;
        private SinglePhaseField wx, wy, wz;
        /// <summary>
        /// signum function
        /// </summary>
        private SinglePhaseField m_sign;

        private VectorField<SinglePhaseField> m_normderivs;
        /// <summary>
        ///  Curvature
        /// </summary>
        private SinglePhaseField m_Curvature;
        private SinglePhaseField dsignedwxd1, dsignedwyd2, dsignedwzd3;

        /// <summary>
        /// Constructor of an abstract surface object.
        /// </summary>
        /// <param name="cont">The context object.</param>
        /// <param name="b">Basis for the field representation of the surface</param>
        /// <param name="deltax">Parameter for approximation of the signum function</param>
        public Surface(GridData cont, Basis b, double deltax) {
            if (cont.Grid.SpatialDimension < 2)
                throw new ArgumentOutOfRangeException("Sorry. This class is created only on the basis of two- and threedimensional grids.");
            m_Context = cont;
            if (b.Degree < 2) throw new ArgumentOutOfRangeException("The basis degree required for the level set representation has to be at least two.");
            m_Basis = b;
            m_Field = new SinglePhaseField(m_Basis, "LevelSet");
            m_gradBasis = new Basis(m_Context, m_Basis.Degree - 1);
            m_derivsBasis = new Basis(m_Context, m_Basis.Degree - 2);

        }

        /// <summary>
        /// Basis of the level set. It may be independent of the basis of the solution space.
        /// </summary>

        public Basis basis {
            get {
                return m_Basis;
            }
        }

        /// <summary>
        /// Level set as single phase field
        /// </summary>

        public SinglePhaseField fieldrepresentation {
            get {
                return m_Field;
            }
        }
        /// <summary>
        /// Normalized gradient of the level set.
        /// </summary>    

        public VectorField<SinglePhaseField> normalvec {
            get {
                return m_Gradw;
            }
        }
        /// <summary>
        /// Partial derivatives of the vector normal to the level set.
        /// </summary>
        /// 
        public VectorField<SinglePhaseField> normalvecderivs {
            get {
                return m_normderivs;
            }
        }



        /// <summary>
        ///The inverse divergence of the normal vector to the surface which gives the divergence. 
        /// </summary>

        public SinglePhaseField curvature {
            get {
                return m_Curvature;
            }
        }

        /// <summary>
        /// Signum of the level set function. 
        /// </summary>

        public SinglePhaseField sign {

            get {
                return m_sign;
            }

        }

        /// <summary>
        ///  Function that describes the (initial) level set. Needs to be implemented.
        /// </summary>
        abstract public void LevSetInit(MultidimensionalArray inp, MultidimensionalArray outp);


        /// <summary>
        /// 
        /// </summary>
        /// <param name="testnodes"></param>
        /// <param name="quadwghts"></param>
        /// <param name="nw"></param>
        abstract public void CreateNodesAndWeights(out double[,] testnodes, out double[] quadwghts, int nw);
        

        /// <summary>
        ///  In this method, the vector normal to the surface as well as certain derivatives 
        /// of it that we intend to use are set up
        /// </summary>

        protected void NormalVec(double deltax) {

            wx = new SinglePhaseField(m_gradBasis, "wx");
            wy = new SinglePhaseField(m_gradBasis, "wy");

            if (m_Context.Grid.SpatialDimension == 2) {
                m_Gradw = new VectorField<SinglePhaseField>(wx, wy);
            } else if (m_Context.Grid.SpatialDimension == 3) {
                wz = new SinglePhaseField(m_gradBasis, "wz");
                m_Gradw = new VectorField<SinglePhaseField>(wx, wy, wz);
            } else {
                throw new NotSupportedException("only spatial dimension 2 and 3 are supported.");
            }
            SinglePhaseField absval = new SinglePhaseField(m_gradBasis);

            m_Gradw[0].Derivative(1.0, m_Field, 0);
            m_Gradw[1].Derivative(1.0, m_Field, 1);

            if (m_Context.Grid.SpatialDimension == 3) {
                m_Gradw[2].Derivative(1.0, m_Field, 2);
            }
            /*
             * Normalization of the gradient vector
             * Might be better implemented using a 
             * Function
             */

            absval.ProjectAbs(1.0, m_Gradw);
            m_Gradw[0].ProjectQuotient(1.0, m_Gradw[0], absval, null, false);
            m_Gradw[1].ProjectQuotient(1.0, m_Gradw[1], absval, null, false);

            if (m_Context.Grid.SpatialDimension == 3) {
                m_Gradw[2].ProjectQuotient(1.0, m_Gradw[2], absval, null, false);
            }
            SinglePhaseField n1x = new SinglePhaseField(m_derivsBasis, "n1x");
            SinglePhaseField n2x = new SinglePhaseField(m_derivsBasis, "n2x");
            SinglePhaseField n1y = new SinglePhaseField(m_derivsBasis, "n1y");
            SinglePhaseField n2y = new SinglePhaseField(m_derivsBasis, "n2y");

            if (m_Context.Grid.SpatialDimension == 3) {
                SinglePhaseField n3x = new SinglePhaseField(m_derivsBasis, "n3x");
                SinglePhaseField n3y = new SinglePhaseField(m_derivsBasis, "n3y");
                SinglePhaseField n1z = new SinglePhaseField(m_derivsBasis, "n1z");
                SinglePhaseField n2z = new SinglePhaseField(m_derivsBasis, "n2z");
                SinglePhaseField n3z = new SinglePhaseField(m_derivsBasis, "n3z");

                m_normderivs = new VectorField<SinglePhaseField>(n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z);

                /* Partial derivatives of the normal vector.
                 * These quantities are used in the product rule applied to the surface projection
                 * and to compute the mean curvature
                 */

                m_normderivs[0].Derivative(1.0, m_Gradw[0], 0);
                m_normderivs[1].Derivative(1.0, m_Gradw[1], 0);
                m_normderivs[2].Derivative(1.0, m_Gradw[2], 0);
                m_normderivs[3].Derivative(1.0, m_Gradw[0], 1);
                m_normderivs[4].Derivative(1.0, m_Gradw[1], 1);
                m_normderivs[5].Derivative(1.0, m_Gradw[2], 1);
                m_normderivs[6].Derivative(1.0, m_Gradw[0], 2);
                m_normderivs[7].Derivative(1.0, m_Gradw[1], 2);
                m_normderivs[8].Derivative(1.0, m_Gradw[2], 2);
            } else {

                m_normderivs = new VectorField<SinglePhaseField>(n1x, n2x, n1y, n2y);

                /* Partial derivatives of the normal vector.
                 * These quantities are used in the product rule applied to the surface projection
                 * and to compute the mean curvature
                 */

                m_normderivs[0].Derivative(1.0, m_Gradw[0], 0);
                m_normderivs[1].Derivative(1.0, m_Gradw[1], 0);
                m_normderivs[2].Derivative(1.0, m_Gradw[0], 1);
                m_normderivs[3].Derivative(1.0, m_Gradw[1], 1);

            }
            /*
             * Divergence of the surface normal vector yields the total curvature 
             * which is twice the mean curvature
             * ATTENTION: here with inverse sign!!!!!
             */
            m_Curvature = new SinglePhaseField(m_derivsBasis);
            m_Curvature.Acc(1.0, m_normderivs[0]);
            m_Curvature.Acc(1.0, m_normderivs[4]);
            if (m_Context.Grid.SpatialDimension == 3) {
                m_Curvature.Acc(1.0, m_normderivs[8]);
            }

            m_sign = new SinglePhaseField(m_Basis, "sign");
            SinglePhaseField quot = new SinglePhaseField(m_Basis, "quot");
            SinglePhaseField sqrt = new SinglePhaseField(m_Basis, "sqrt");

            sqrt.ProjectProduct(1.0, m_Field, m_Field);
            sqrt.AccConstant(0.01);

            quot.ProjectPow(1.0, sqrt, deltax * deltax);
            m_sign.ProjectQuotient(1.0, m_Field, quot);

        }
        /// <summary>
        /// This method is meant for the extension of quantities.
        /// In that equation \Nabla \Phi(x)*sign(\Phi(x)) is needed.
        /// </summary>

        protected void SignedNormal() {

            SinglePhaseField signedwx = new SinglePhaseField(m_gradBasis, "signedwx");
            SinglePhaseField signedwy = new SinglePhaseField(m_gradBasis, "signedwy");
            dsignedwxd1 = new SinglePhaseField(m_derivsBasis, "dsignedwxd1");
            dsignedwyd2 = new SinglePhaseField(m_derivsBasis, "dsignedwyd2");

            signedwx.ProjectProduct(1.0, m_sign, m_Gradw[0]);
            signedwy.ProjectProduct(1.0, m_sign, m_Gradw[1]);

            dsignedwxd1.Derivative(1.0, signedwx, 0);
            dsignedwyd2.Derivative(1.0, signedwy, 1);

            if (m_Context.Grid.SpatialDimension == 3) {
                SinglePhaseField signedwz = new SinglePhaseField(m_gradBasis, "signedwz");
                dsignedwzd3 = new SinglePhaseField(m_derivsBasis, "dsignedwzd3");
                signedwz.ProjectProduct(1.0, m_sign, m_Gradw[2]);
                dsignedwzd3.Derivative(1.0, signedwz, 2);
            }

        }


    }
}
