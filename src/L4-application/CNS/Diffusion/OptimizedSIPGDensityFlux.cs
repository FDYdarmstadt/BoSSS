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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Diffusion;
using ilPSP;
using System;
using System.Collections.Generic;

namespace CNS.Diffusion {

    /// <summary>
    /// Dummy Class for Diffusion density flux, i.e \f$ F_v(U_1)=0\f$  
    /// </summary>
    class OptimizedSIPGDensityFlux : INonlinear2ndOrderForm {

        private CNSControl config;

        private ISpeciesMap speciesMap;

        private IBoundaryConditionMap boundaryMap;

        public bool AdiabaticWall { get; set; }

        private IGridData gridData;

        public OptimizedSIPGDensityFlux(CNSControl config, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, IGridData gridData, Func<MultidimensionalArray> cellMetric) {
            this.config = config;
            this.speciesMap = speciesMap;
            this.boundaryMap = boundaryMap;
            this.gridData = gridData;
        }

        #region IEquationComponent Members

        /// <summary>
        /// <see cref="CompressibleEnvironment.PrimalArgumentOrdering"/>
        /// </summary>
        IList<string> IEquationComponent.ArgumentOrdering {
            get {
                return CompressibleEnvironment.PrimalArgumentOrdering;
            }
        }

        /// <summary>
        /// Empty (i.e., no parameters are used)
        /// </summary>
        IList<string> IEquationComponent.ParameterOrdering {
            get {
                return null;
            }
        }
        #endregion

        #region IEdgeComponent Members
        /// <summary>
        /// <see cref="IEdgeForm"/>
        /// </summary>
        TermActivationFlags IEdgeForm.BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        /// <summary>
        /// <see cref="IEdgeForm"/>
        /// </summary>
        TermActivationFlags IEdgeForm.InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        /// <summary>
        /// not needed
        /// </summary>
        double IEdgeForm.InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            throw new NotImplementedException();
        }
        /// <summary>
        /// not needed
        /// </summary>
        double IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotImplementedException();
        }
        #endregion


        #region INonlineEdgeform_GradV Members


        void INonlinEdgeForm_GradV.InternalEdge(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout, MultidimensionalArray fIN, MultidimensionalArray fOT) {
            //Do nothing
        }

        void INonlinEdgeForm_GradV.BoundaryEdge(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin, MultidimensionalArray f) {
            //Do nothing
        }
        #endregion

        #region INonlineEdgeform_V Members
        void INonlinEdgeForm_V.InternalEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            //Do nothing
        }

        void INonlinEdgeForm_V.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {
            //Do Nothing
        }
        #endregion


        #region IVolumeForm Members
        TermActivationFlags IVolumeForm.VolTerms {
            get {
                return (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV);
            }
        }

        double IVolumeForm.VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            throw new NotImplementedException();
        }
        #endregion

        #region INonLinearVolumeForm_GradV Members

        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {
            // Do Nothing
        }
        #endregion
    }
}
