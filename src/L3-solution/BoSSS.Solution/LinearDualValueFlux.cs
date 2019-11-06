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
using System.Text;

using BoSSS.Foundation;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Utility class which helps the user in creating the function matrices that are needed
    /// in an implementatation of LinearDualValueFlux;
    /// The user only needs to override
    /// <see cref="InnerEdgeFlux"/>
    /// and <see cref="BorderEdgeFlux_"/>)
    /// where he is able to implement the linear function as 
    /// an algebraic formula. All function matrixes and offsets (aka. intercept) are constructed from the user-defined
    /// functions by this class.
    /// </summary>
    abstract public class LinearDualValueFlux : IEdgeForm {

        /// <summary>
        /// Implementation of the LinearDualValueFlux for Inner-Edges, based on <see cref="IEdgeForm"/>
        /// </summary>
        /// <param name="inp"><see cref="IEdgeForm"/></param>
        /// <param name="_uA"><see cref="IEdgeForm"/></param>
        /// <param name="_uB"><see cref="IEdgeForm"/></param>
        /// <param name="_Grad_uA"><see cref="IEdgeForm"/></param>
        /// <param name="_Grad_uB"><see cref="IEdgeForm"/></param>
        /// <param name="_vA"><see cref="IEdgeForm"/></param>
        /// <param name="_vB"><see cref="IEdgeForm"/></param>
        /// <param name="_Grad_vA"><see cref="IEdgeForm"/></param>
        /// <param name="_Grad_vB"><see cref="IEdgeForm"/></param>
        /// <returns></returns>
        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Flx_InCell, Flx_OutCell;
            this.InnerEdgeFlux(ref inp, _uA, _uB, out Flx_InCell, out Flx_OutCell);
            return Flx_InCell * _vA + Flx_OutCell * _vB;
        }

        /// <summary>
        /// Implementation of the LinearDualValueFlux for Border-Edges, based on <see cref="IEdgeForm"/>
        /// </summary>
        /// <param name="inp"><see cref="IEdgeForm"/></param>
        /// <param name="_uA"><see cref="IEdgeForm"/></param>
        /// <param name="_Grad_uA"><see cref="IEdgeForm"/></param>
        /// <param name="_vA"><see cref="IEdgeForm"/></param>
        /// <param name="_Grad_vA"><see cref="IEdgeForm"/></param>
        /// <returns></returns>
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Flx_InCell;
            BorderEdgeFlux_(ref inp, _uA, out Flx_InCell);
            return Flx_InCell * _vA;
        }

        /// <summary>
        /// to be defined by user;
        /// </summary>
        abstract public IList<string> ArgumentOrdering { get; }

        /// <summary>
        /// null by default
        /// </summary>
        virtual public IList<string> ParameterOrdering { get { return null; } }

        
        /// <summary>
        /// override this method to implement the dual-value 'flux' at interior edges
        /// </summary>
        protected abstract void InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout, out double FluxInCell, out double FluxOutCell);

        /// <summary>
        /// override this method to implement the 'flux' at boundary edges
        /// </summary>
        protected abstract void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell);
        
        /// <summary>
        /// depends on products of trial and test function
        /// </summary>
        virtual public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        /// <summary>
        /// depends on products of trial and test function
        /// </summary>
        virtual public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; ; }
        }

        
    }
}
