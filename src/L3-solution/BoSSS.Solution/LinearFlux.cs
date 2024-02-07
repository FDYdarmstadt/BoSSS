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
using BoSSS.Foundation;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Utility class which helps the user in creating the function matrices that are needed
    /// in an implementation of Linear-Flux; The user only needs to override
    /// <see cref="BorderEdgeFlux"/>, 
    /// <see cref="InnerEdgeFlux"/> 
    /// and <see cref="Flux"/>,
    /// where he is able to implement the linear function as 
    /// an algebraic formula. All function matrices and offsets (aka. intercept)
    /// are constructed from the user-defined functions by this class.
    /// </summary>
    public abstract class LinearFlux : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {

        /// <summary>
        /// not in use, returning null
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        /// <summary>
        /// override this method to implement the Riemann flux at border edges
        /// </summary>
        protected abstract double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin);

        /// <summary>
        /// override this method to implement the Riemann flux a interior edges
        /// </summary>
        protected abstract double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout);

        /// <summary>
        /// override this method to implement the flux function.
        /// </summary>
        protected abstract void Flux(ref CommonParamsVol inp, double[] U, double[] output);


        #region IEquationComponent Member

        /// <summary>
        /// to be defined by the user/implementor;
        /// </summary>
        abstract public IList<string> ArgumentOrdering {
            get;
        }

        #endregion

        /// <summary>
        /// <see cref="IVolumeForm.VolTerms"/>
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxGradV | TermActivationFlags.GradV;
            }
        }

        //double[] buf = null; // can't be a class member because than multi-thread does not work

        /// <summary>
        /// Implementation of a Bilinear for linear volume fluxes
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            int D = GradV.Length;
            double acc = 0;
            //if (buf == null)
            //    buf = new double[D];
            var buf = new double[D];
            this.Flux(ref cpv, U, buf);
            for (int d = 0; d < D; d++)
                acc += buf[d] * GradV[d];
            return -acc;
        }

        /// <summary>
        /// <see cref="IEdgeForm.BoundaryEdgeTerms"/>
        /// </summary>
        virtual public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        /// <see cref="IEdgeForm.InnerEdgeTerms"/>
        /// </summary>
        virtual public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        /// Calls <see cref="LinearFlux.InnerEdgeFlux"/>
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            return this.InnerEdgeFlux(ref inp, _uA, _uB) * (_vA - _vB);
        }


        /// <summary>
        /// Calls <see cref="LinearFlux.BorderEdgeFlux"/>
        /// </summary>
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return this.BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
