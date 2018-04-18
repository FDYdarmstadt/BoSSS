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
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;

namespace BoSSS.Foundation {

    public partial class DGField {

        /// <summary>
        /// accumulates the Laplacian of field <paramref name="f"/> times
        /// <paramref name="alpha"/> to this field.
        /// </summary>
        /// <remarks>
        /// see <see cref="Laplacian(double, DGField,CellMask)"/>;
        /// </remarks>
        public void Laplacian(double alpha, DGField f) {
            Laplacian(alpha, f, null);
        }

        /// <summary>
        /// accumulates the Laplacian of field <paramref name="f"/> times
        /// <paramref name="alpha"/> to this field.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <param name="em">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask for this case would be e.g. 
        /// <see cref="BoSSS.Foundation.Grid.GridData.BoundaryCells"/>.)<br/>
        /// if null, the computation is carried out in the whole domain
        /// </param>
        /// <remarks>
        /// This method is based on
        /// <see cref="DGField.Derivative(double,DGField,int,CellMask)"/>, i.e.
        /// it calculates derivatives by analytic cell-by-cell derivation of
        /// the DG polynomials;<br/> Note that of the Laplacian requires the
        /// allocation of a temporary DG field.
        /// </remarks>
        virtual public void Laplacian(double alpha, DGField f, CellMask em) {
            using (new FuncTrace()) {
                DGField tmp = (DGField)this.Clone();

                int D = GridDat.SpatialDimension;
                for (int d = 0; d < D; d++) {
                    tmp.Clear();
                    tmp.Derivative(1.0, f, d, em);
                    this.Derivative(alpha, tmp, d, em);
                }
            }
        }

        /// <summary>
        /// accumulates the Laplacian of field <paramref name="f"/> times
        /// <paramref name="alpha"/> to this field.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <param name="optionalSubGrid">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask would be see e.g. 
        /// <see cref="GridData.BoundaryCells"/>.)
        /// <br/>
        /// if null, the computation is carried out in the whole domain.
        /// </param>
        /// <param name="tmp">
        /// temporary storage needed for the 1st derivatives of
        /// <paramref name="f"/>; If null, a clone of <paramref name="f"/>
        /// is taken.
        /// </param>
        /// <param name="bndMode_1stDeriv"></param>
        /// <param name="bndMode_2ndDeriv"></param>
        /// <remarks>
        /// This method is based on <see cref="DGField.DerivativeByFlux"/>,
        /// i.e. it calculates derivatives by central-difference fluxes;<br/>
        /// Note that of the Laplacian requires the allocation of a temporary DG field.
        /// </remarks>
        virtual public void LaplacianByFlux(double alpha, DGField f, DGField tmp = null,
            SubGrid optionalSubGrid = null,
            SpatialOperator.SubGridBoundaryModes bndMode_1stDeriv = SpatialOperator.SubGridBoundaryModes.OpenBoundary,
            SpatialOperator.SubGridBoundaryModes bndMode_2ndDeriv = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {
            using (new FuncTrace()) {
                if (tmp == null)
                    tmp = (DGField)f.Clone();

                int D = GridDat.SpatialDimension;
                for (int d = 0; d < D; d++) {
                    tmp.Clear();
                    tmp.DerivativeByFlux(1.0, f, d, optionalSubGrid, bndMode_1stDeriv);
                    this.DerivativeByFlux(alpha, tmp, d, optionalSubGrid, bndMode_2ndDeriv);
                }
            }
        }

        /// <summary>
        /// accumulates the divergence of vector field <paramref name="vec"/>
        /// times <paramref name="alpha"/> to this field.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="vec"></param>
        /// <param name="em"></param>
        /// <remarks>
        /// This method is based on <see cref="Derivative(double,DGField,int,CellMask)"/>;
        /// </remarks>
        public void Divergence<T>(double alpha, VectorField<T> vec, CellMask em = null) where T : DGField {
            using (new FuncTrace()) {
                if (vec.Dim != GridDat.SpatialDimension)
                    throw new ArgumentException(
                        "wrong number of components in vector field.", "vec");

                int D = GridDat.SpatialDimension;
                for (int d = 0; d < D; d++)
                    this.Derivative(alpha, vec[d], d, em);
            }
        }

        /// <summary>
        /// accumulates the divergence of vector field <paramref name="vec"/>
        /// times <paramref name="alpha"/>
        /// to this field.
        /// </summary>
        /// <remarks>
        /// This method is based on <see cref="DerivativeByFlux"/>;
        /// </remarks>
        public void DivergenceByFlux<T>(double alpha, VectorField<T> vec,
            SubGrid optionalSubGrid = null, SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary) where T : DGField {
            using (new FuncTrace()) {
                if (vec.Dim != GridDat.SpatialDimension)
                    throw new ArgumentException("wrong number of components in vector field.", "vec");

                int D = GridDat.SpatialDimension;
                for (int d = 0; d < D; d++)
                    this.DerivativeByFlux(alpha, vec[d], d, optionalSubGrid, bndMode);
            }
        }

        /// <summary>
        /// accumulates the curl of 2D DG vector field <paramref name="vec"/> 
        ///  times <paramref name="alpha"/>
        /// to this vector field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* Curl(<paramref name="vec"/>)
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="vec"></param>
        /// <remarks>
        /// This method is based on
        /// <see cref="DGField.Derivative(double,DGField,int)"/>, i.e. it
        /// calculates derivatives by analytic cell-by-cell derivation of the
        /// DG polynomials;
        /// </remarks>
        /// <seealso cref="VectorField{T}.Curl3D"/>
        public void Curl2D<T>(double alpha, VectorField<T> vec) where T : DGField {
            Curl2D<T>(alpha, vec, null);
        }

        /// <summary>
        /// accumulates the curl of a 2D DG vector field <paramref name="vec"/> 
        ///  times <paramref name="alpha"/>
        /// to this vector field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* Curl(<paramref name="vec"/>)
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="vec"></param>
        /// <param name="em">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask for this case would be e.g. 
        /// <see cref="BoSSS.Foundation.Grid.GridData.BoundaryCells"/>.)<br/>
        /// if null, the computation is carried out in the whole domain
        /// </param>
        /// <remarks>
        /// This method is based on
        /// <see cref="DGField.Derivative(double,DGField,int)"/>, i.e. it
        /// calculates derivatives by analytic cell-by-cell derivation of the
        /// DG polynomials;
        /// </remarks>
        /// <seealso cref="VectorField{T}.Curl3D(double, VectorField{T}, CellMask)"/>
        public void Curl2D<T>(double alpha, VectorField<T> vec, CellMask em) where T : DGField {
            using (new FuncTrace()) {
                if (vec.Dim != 2)
                    throw new ArgumentException("vector field must be 2-dim.", "vec");

                this.Derivative(alpha, vec[1], 0, em);
                this.Derivative(-alpha, vec[0], 1, em);
            }
        }

        /// <summary>
        /// accumulates the curl of 2D DG vector field <paramref name="vec"/> 
        /// times <paramref name="alpha"/>
        /// to this vector field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* Curl(<paramref name="vec"/>)
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="vec"></param>
        /// <param name="optionalSubGrid">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask would be see e.g. 
        /// <see cref="GridData.BoundaryCells"/>.)
        /// <br/>
        /// if null, the computation is carried out in the whole domain.
        /// </param>
        /// <remarks>
        /// This method is based on <see cref="DGField.DerivativeByFlux"/>, i.e.
        /// it calculates derivatives by central-difference fluxes;
        /// </remarks>
        /// <seealso cref="VectorField{T}.Curl3DByFlux"/>
        /// <param name="bndMode"></param>
        public void Curl2DByFlux<T>(double alpha, VectorField<T> vec,
            SubGrid optionalSubGrid = null,
            SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary)
            where T : DGField {
            //diff(v(x, y), x)-(diff(u(x, y), y))
            using (new FuncTrace()) {

                if (vec.Dim != 2)
                    throw new ArgumentException("vector field must be 2-dim.", "vec");

                this.DerivativeByFlux(alpha, vec[1], 0, optionalSubGrid, bndMode);
                this.DerivativeByFlux(-alpha, vec[0], 1, optionalSubGrid, bndMode);
            }
        }

        /// <summary>
        /// accumulates the derivative of DG field <paramref name="f"/> 
        /// (along the <paramref name="d"/>-th axis) times <paramref name="alpha"/>
        /// to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* \f$ \frac{\partial}{\partial x_d} \f$ <paramref name="f"/>;
        /// </summary>
        /// <param name="f"></param>
        /// <param name="d">
        /// 0 for the x-derivative, 1 for the y-derivative, 2 for the z-derivative
        /// </param>
        /// <param name="alpha">
        /// scaling of <paramref name="f"/>;
        /// </param>
        virtual public void Derivative(double alpha, DGField f, int d) {
            Derivative(alpha, f, d, null);
        }

        /// <summary>
        /// symbolic derivation, cell-by-cell;<br/>
        /// accumulates the derivative of DG field <paramref name="f"/> 
        /// (along the <paramref name="d"/>-th axis) times <paramref name="alpha"/>
        /// to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>
        /// \f$ \cdot \frac{\partial}{\partial x_d}  \f$ <paramref name="f"/>;
        /// </summary>
        /// <param name="f"></param>
        /// <param name="d">
        /// 0 for the x-derivative, 1 for the y-derivative, 2 for the
        /// z-derivative
        /// </param>
        /// <param name="alpha">
        /// scaling of <paramref name="f"/>;
        /// </param>
        /// <param name="em">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask for this case would be e.g. 
        /// <see cref="BoSSS.Foundation.Grid.GridData.BoundaryCells"/>.)<br/>
        /// if null, the computation is carried out in the whole domain
        /// </param>
        /// <remarks>
        /// The derivative is calculated by a cell-by-cell (symbolic)
        /// derivation of the DG polynomials, therefor the (effective) DG
        /// polynomial degree is one lower than the degree of
        /// <paramref name="f"/>;<br/>
        /// In comparison to <see cref="DerivativeByFlux"/>, this method should
        /// be much faster, because no quadrature is involved;
        /// </remarks>
        abstract public void Derivative(double alpha, DGField f, int d, CellMask em);

        /// <summary>
        /// used by <see cref="DerivativeByFlux"/> to create an instance of
        /// <see cref="DerivativeFlux"/>; Override this method to use
        /// <see cref="DerivativeByFlux"/> with your own flux implementation.
        /// </summary>
        /// <param name="d">
        /// spatial direction of the derivation: 0 for derivation in
        /// x-direction, 1 for y-direction, ...
        /// </param>
        /// <param name="Identification">not used</param>
        /// <returns>
        /// </returns>
        protected virtual INonlinearFlux CreateDerivativeFlux(int d, string Identification) {
            return new DerivativeFlux(d);
        }

        /// <summary>
        /// A flux to compute du/dx, used by <see cref="DerivativeByFlux"/>;
        /// </summary>
        protected class DerivativeFlux : INonlinearFlux {

            /// <summary>
            /// 
            /// </summary>
            /// <param name="d">
            /// spatial direction of the derivation: 0 for derivation in
            /// x-direction, 1 for y-direction, ...
            /// </param>
            public DerivativeFlux(int d) {
                this.m_d = d;
            }

            /// <summary>
            /// spatial direction of the derivation: 0 for derivation in
            /// x-direction, 1 for y-direction, ...
            /// </summary>
            public int m_d;

            /// <summary>
            /// see <see cref="INonlinearFlux.BorderEdgeFlux"/>
            /// </summary>
            virtual public void BorderEdgeFlux(double time, int jEdge, MultidimensionalArray X, MultidimensionalArray normal, bool flipNormal, byte[] EdgeTags, int EdgeTagsOffset, MultidimensionalArray[] Uin, int Offset, int Lenght, MultidimensionalArray Output) {
                int NoOfNodes = Uin[0].GetLength(1);

                int _d = this.m_d;

                double normalSign = 1;// flipNormal ? -1 : 1;
                for (int e = 0; e < Lenght; e++) {
                    for (int n = 0; n < NoOfNodes; n++) {
                        //double x = X[e + Offset, n, 0];
                        //double Usoll = x * 3;
                        //double ErrIn = (Uin[0][e + Offset, n] - Usoll).Pow2();
                        
                        //if(ErrIn > 1.0e-6)
                        //    Console.WriteLine("Uin fucked up. / boundary");
                        
                        Output[e + Offset, n] += normal[e + Offset, n, _d] * Uin[0][e + Offset, n] * normalSign;
                    }
                }
            }

            /// <summary>
            /// A central-difference flux for computing the
            /// <see cref="m_d"/>-th derivative;<br/>
            /// see <see cref="INonlinearFlux.InnerEdgeFlux"/>
            /// </summary>
            virtual public void InnerEdgeFlux(double time, int jEdge, MultidimensionalArray X, MultidimensionalArray normal, MultidimensionalArray[] Uin, MultidimensionalArray[] Uot, int Offset, int Lenght, MultidimensionalArray Output) {
                int NoOfNodes = Uin[0].GetLength(1);

                int _d = this.m_d;

                for (int e = 0; e < Lenght; e++) {
                    for (int n = 0; n < NoOfNodes; n++) {
                        //double ErrJump = (Uin[0][e + Offset, n] - Uot[0][e + Offset, n]).Pow2();
                        //double x = X[e + Offset, n, 0];
                        //double y = X[e + Offset, n, 1];
                        //double Usoll = 3 * x + z;
                        //double ErrIn = (Uin[0][e + Offset, n] - Usoll).Pow2();
                        //double ErrOt = (Uot[0][e + Offset, n] - Usoll).Pow2();

                        //if(ErrJump > 1.0e-6)
                        //    Console.WriteLine("[[U]] fucked up. ");
                        //if(ErrIn > 1.0e-6)
                        //    Console.WriteLine("Uin fucked up. ");
                        //if(ErrOt > 1.0e-6)
                            //Console.WriteLine("Uot fucked up. ");
                        
                        Output[e + Offset, n] += normal[e + Offset, n, _d] * (Uin[0][e + Offset, n] + Uot[0][e + Offset, n]) * (0.5);
                    }
                }
            }

            /// <summary>
            /// see <see cref="INonlinearFlux.Flux"/>
            /// </summary>
            public void Flux(double time, MultidimensionalArray x, MultidimensionalArray[] U, int Offset, int Length, MultidimensionalArray Output) {
                int NoOfNodes = Output.GetLength(1);

                for (int j = 0; j < Length; j++) {
                    for (int n = 0; n < NoOfNodes; n++) {
                        Output[j + Offset, n, m_d] += U[0][j + Offset, n];
                    }
                }
            }

            /// <summary>
            /// see <see cref="IEquationComponent.ArgumentOrdering"/>
            /// </summary>
            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "in" };
                }
            }

            /// <summary>
            /// see <see cref="IEquationComponent.ParameterOrdering"/>
            /// </summary>
            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }

        /// <summary>
        /// accumulates the derivative of DG field <paramref name="f"/> 
        /// (along the <paramref name="d"/>-th axis) times <paramref name="alpha"/>
        /// to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* \f$ \frac{\partial}{\partial x_d} \f$ <paramref name="f"/>;
        /// </summary>
        /// <param name="f"></param>
        /// <param name="d">
        /// 0 for the x-derivative, 1 for the y-derivative, 2 for the
        /// z-derivative
        /// </param>
        /// <param name="alpha">
        /// scaling of <paramref name="f"/>;
        /// </param>
        /// <param name="optionalSubGrid">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask would be see e.g. 
        /// <see cref="GridData.BoundaryCells"/>.)
        /// <br/>
        /// if null, the computation is carried out in the whole domain.
        /// </param>
        /// <param name="bndMode">
        /// if a sub-grid is provided, this determines how the sub-grid
        /// boundary should be treated.
        /// </param>
        /// <remarks>
        /// The derivative is calculated using Riemann flux functions
        /// (central difference);<br/>
        /// In comparison to
        /// <see cref="Derivative(double, DGField, int, CellMask)"/>, this method
        /// should be slower, but produce more sane results, especially for
        /// fields of low polynomial degree (0 or 1);
        /// </remarks>
        virtual public void DerivativeByFlux(double alpha, DGField f, int d, SubGrid optionalSubGrid = null, SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {
            int D = this.Basis.GridDat.SpatialDimension;
            if (d < 0 || d >= D)
                throw new ArgumentException("spatial dimension out of range.", "d");
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            SpatialOperator d_dx = new SpatialOperator(1, 1, QuadOrderFunc.Linear(),"in", "out");
            var flux = CreateDerivativeFlux(d, f.Identification);
            d_dx.EquationComponents["out"].Add(flux);
            d_dx.Commit();

            EdgeMask emEdge = (optionalSubGrid != null) ? optionalSubGrid.AllEdgesMask : null;
            CellMask emVol = (optionalSubGrid != null) ? optionalSubGrid.VolumeMask : null;

            var ev = d_dx.GetEvaluatorEx(
                new CoordinateMapping(f), null, this.Mapping,
                edgeQrCtx: new Quadrature.EdgeQuadratureScheme(true, emEdge),
                volQrCtx: new Quadrature.CellQuadratureScheme(true, emVol),
                sgrd: optionalSubGrid, subGridBoundaryTreatment: bndMode);

            ev.Evaluate<CoordinateVector>(alpha, 1.0, this.CoordinateVector);
        }
    }
}
