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
using CNS.Boundary;
using CNS.Exception;
using CNS.MaterialProperty;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CNS.Diffusion {

    /// <summary>
    /// Abstract class for the components of SIPG flux as described in Lecture Notes Hartmann2008 (p. 92)
    /// </summary>
    public abstract class SIPGFlux : IVolumeForm, IEdgeForm {

        /// <summary>
        /// Boundary value definition
        /// </summary>
        protected IBoundaryConditionMap boundaryMap;

        /// <summary>
        /// Configuration options
        /// </summary>
        protected CNSControl config;

        /// <summary>
        /// Mapping that determines the active species in some point.
        /// </summary>
        protected ISpeciesMap speciesMap;

        /// <summary>
        /// Common GridData
        /// </summary>
        protected GridData gridData;

        /// <summary>
        /// penalty factor used in the SIPG discretization
        /// </summary>
        protected double penaltyFactor;

        /// <summary>
        /// Dimension of the problem
        /// </summary>
        protected int dimension;

        /// <summary>
        /// Dictionary, especially needed adiabatic wall
        /// </summary>
        protected Dictionary<byte, bool> edgeTagBool = new Dictionary<byte, bool>();

        /// <summary>
        /// Frequently used constant, see AnnualReport2014_SKE
        /// </summary>
        protected double alphaPlus43;

        /// <summary>
        /// Frequently used constant, see AnnualReport2014_SKE
        /// </summary>
        protected double alphaMinus23;

        /// <summary>
        /// Frequently used constant, see AnnualReport2014_SKE
        /// </summary>
        protected double alphaPlus13;

        // [dimension, dimension, j]
        // [ k , l , j] --> indices according to Hartmann2008 or AnnualReport2014_SKE, i is omitted due to FluxComponents in BoSSS
        double[,,] GTensorIn;
        double[,,] GTensorOut;
        StateVector stateIn;
        StateVector stateOut;

        /// <summary>
        /// Flux Component
        /// </summary>
        protected abstract int Component {
            get;
        }


        private MultidimensionalArray cellMetricBack;

        private Func<MultidimensionalArray> cellMetricFunc;


        /// <summary>
        /// Metric for each cell, calculation depends on IBM or non IBM case
        /// necessary for the calculation of the penalty factor
        /// </summary>
        protected MultidimensionalArray cellMetric {
            get {
                if (cellMetricBack == null) {
                    cellMetricBack = cellMetricFunc();
                }
                return cellMetricBack;
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="config">Configuration options</param>
        /// <param name="boundaryMap">Boundary value definition</param>
        /// <param name="speciesMap">Mapping that determines the active species in some point</param>
        /// <param name="gridData"></param>
        /// <param name="cellMetric"></param>
        public SIPGFlux(CNSControl config, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, GridData gridData, Func<MultidimensionalArray> cellMetric) {
            this.config = config;
            this.boundaryMap = boundaryMap;
            this.speciesMap = speciesMap;
            this.gridData = gridData;
            this.dimension = CNSEnvironment.NumberOfDimensions;
            this.cellMetricFunc = cellMetric;

            //Fills the dictionary, to avoid later string comparison
            foreach (byte edgeTag in gridData.Edges.EdgeTags) {
                if (boundaryMap.EdgeTagNames[edgeTag].StartsWith("adiabaticWall", StringComparison.InvariantCultureIgnoreCase)) {
                    edgeTagBool[edgeTag] = true;
                } else {
                    edgeTagBool[edgeTag] = false;
                }
            }
            // calculation of the penalty factor
            double p = new int[] { config.DensityDegree, config.MomentumDegree, config.EnergyDegree }.Max();
            penaltyFactor = config.SIPGPenaltyScaling * p * p;

            //defining some constants
            double alpha = config.ViscosityRatio;
            alphaPlus43 = alpha + (4.0 / 3.0);
            alphaPlus13 = alpha + (1.0 / 3.0); 
            alphaMinus23 = alpha - (2.0 / 3.0);

            stateIn = new StateVector(new double[dimension + 2], speciesMap.GetMaterial(double.NaN));
            stateOut = new StateVector(new double[dimension + 2], speciesMap.GetMaterial(double.NaN));

            if (config.EquationOfState is IdealGas == false) {
                throw new ConfigurationException("SIPG flux currently only works for ideal gases");
            }

            GTensorIn = new double[dimension, dimension, dimension + 2];
            GTensorOut = new double[dimension, dimension, dimension + 2];
        }

        /// <summary>
        /// <see cref="CNSEnvironment.PrimalArgumentOrdering"/>
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return CNSEnvironment.PrimalArgumentOrdering.ToList();
            }
        }

        /// <summary>
        /// Empty (i.e., no parameters are used)
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return null;
            }
        }


        /// <summary>
        /// Implementation of the SIPG volume fluxes
        /// </summary>
        /// <param name="cpv">
        /// <see cref="IVolumeForm.VolumeForm(ref CommonParamsVol, double[], double[,], double, double[])"/></param>
        /// <param name="U">
        /// <see cref="IVolumeForm.VolumeForm(ref CommonParamsVol, double[], double[,], double, double[])"/></param>
        /// <param name="GradU">
        /// <see cref="IVolumeForm.VolumeForm(ref CommonParamsVol, double[], double[,], double, double[])"/></param>
        /// <param name="V">
        /// <see cref="IVolumeForm.VolumeForm(ref CommonParamsVol, double[], double[,], double, double[])"/></param>
        /// <param name="GradV">
        /// <see cref="IVolumeForm.VolumeForm(ref CommonParamsVol, double[], double[,], double, double[])"/></param>
        /// <returns>volume flux for this specific cell/node</returns>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            stateIn = new StateVector(U, speciesMap.GetMaterial(double.NaN));
            double ret = 0.0;

            UpdateTensorComponent(stateIn, false, GTensorIn, cpv.jCell);

            for (int k = 0; k < dimension; k++) {
                for (int l = 0; l < dimension; l++) {

                    for (int j = 0; j < dimension + 2; j++) {
                        ret += GTensorIn[k, l, j] * GradU[j, l] * GradV[k];
                    }
                }
            }


            return ret;
        }

        /// <summary>
        /// Weakly imposition of the boundary conditions (defined by the <see cref="CommonParamsBnd.EdgeTag"/>) by calculating the 
        /// outer value via a subclass of <see cref="BoundaryCondition"/>
        /// </summary>
        /// <param name="inp">
        /// <see cref="IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd, double[], double[,], double, double[])"/></param>
        /// <param name="_uA">
        /// <see cref="IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd, double[], double[,], double, double[])"/></param>
        /// <param name="_Grad_uA">
        /// <see cref="IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd, double[], double[,], double, double[])"/></param>
        /// <param name="_vA">
        /// <see cref="IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd, double[], double[,], double, double[])"/></param>
        /// <param name="_Grad_vA">
        /// <see cref="IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd, double[], double[,], double, double[])"/></param>
        /// <returns>Flux across the boundary edge</returns>
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            stateIn = new StateVector(_uA, speciesMap.GetMaterial(double.NaN));
            stateOut = boundaryMap.GetBoundaryState(inp.EdgeTag, 0.0, inp.X, inp.Normale, stateIn);

            bool adiabaticWall = edgeTagBool[inp.EdgeTag];
            double[] _uBC = stateOut.ToArray();

            double ret = 0.0;
            double Penalty = penalty(inp.jCellIn, -1);

            UpdateTensorComponent(stateOut, adiabaticWall, GTensorOut, inp.jCellIn);

            for (int k = 0; k < dimension; k++) {
                for (int l = 0; l < dimension; l++) {

                    for (int j = 0; j < dimension + 2; j++) {
                        // _Grad_vB=0 per definition of the Ansatz functions
                        ret -= (GTensorOut[k, l, j] * _Grad_vA[k]) * (_uA[j] - _uBC[j]) * inp.Normale[l];
                        // Term 2
                        // Grad_uB = grad_uA, implicit boundary condition for viscous stress tensor!
                        ret -= (GTensorOut[k, l, j] * _Grad_uA[j, l]) * (_vA - 0) * inp.Normale[k];
                        //ret -= (G_Out[k, l][j] * grad_uB[j]) * (_vA - 0) * inp.Normale[k];
                        //Penalty-Term
                        ret += (GTensorOut[k, l, j]) * (_uA[j] - _uBC[j]) * inp.Normale[l] * Penalty * (_vA - 0) * inp.Normale[k];
                    }
                }
            }

#if DEBUG
            if (double.IsNaN(ret)) {
                throw new System.Exception();
            }
#endif

            return ret;
        }




        public static int EVIL_HACK_CELL_INDEX = -1;



        /// <summary>
        /// Calculates the flux across an inner edge
        /// </summary>
        /// <param name="inp">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_uA">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_uB">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_Grad_uA">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_Grad_uB">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_vA">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_vB">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_Grad_vA">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <param name="_Grad_vB">
        /// <see cref="IEdgeForm.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/></param>
        /// <returns>Flux across an inner edge</returns>
        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            stateIn = new StateVector(_uA, speciesMap.GetMaterial(double.NaN));
            stateOut = new StateVector(_uB, speciesMap.GetMaterial(double.NaN));


            double ret = 0.0;
            double Penalty = penalty(inp.jCellIn, inp.jCellOut);


            UpdateTensorComponent(stateIn, false, GTensorIn, inp.jCellIn);
            UpdateTensorComponent(stateOut, false, GTensorOut, inp.jCellOut);

            for (int k = 0; k < dimension; k++) {
                for (int l = 0; l < dimension; l++) {

                    for (int j = 0; j < dimension + 2; j++) {
                        // Term 1
                        ret -= 0.5 * (GTensorIn[k, l, j] * _Grad_vA[k] + GTensorOut[k, l, j] * _Grad_vB[k]) * (_uA[j] - _uB[j]) * inp.Normale[l];
                        // Term 2
                        ret -= 0.5 * (GTensorIn[k, l, j] * _Grad_uA[j, l] + GTensorOut[k, l, j] * _Grad_uB[j, l]) * (_vA - _vB) * inp.Normale[k];
                        // Term 3: Penalty
                        ret += 0.5 * (GTensorIn[k, l, j] + GTensorOut[k, l, j]) * (_uA[j] - _uB[j]) * inp.Normale[l] * Penalty * (_vA - _vB) * inp.Normale[k];
                    }
                }
            }



#if DEBUG
            if (double.IsNaN(ret)) {
                throw new System.Exception();
            }
#endif


            return ret;
        }

        /// <summary>
        /// Updates the components of the diffusivity tensor, according to the lecture notes Hartmann2008
        /// </summary>
        /// <param name="state">The current flow state</param>
        /// <param name="adiabaticWall">Tensor needs to be adjusted, if the state is at a adiabatic wall</param>
        /// <param name="Tensor">Results are stored in this tensor</param>
        /// <param name="cellIndex"></param>
        abstract protected void UpdateTensorComponent(StateVector state, bool adiabaticWall, double[,,] Tensor, int cellIndex);

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        /// <param name="jCellIn"></param>
        /// <param name="jCellOut"></param>
        /// <returns></returns>
        protected double penalty(int jCellIn, int jCellOut) {
            if (EVIL_HACK_CELL_INDEX >= 0) {

                double cj_in_bla = cellMetric[EVIL_HACK_CELL_INDEX];
                double mu_bla = penaltyFactor * cj_in_bla;
                return mu_bla;
            }




            double cj_in = cellMetric[jCellIn];
            double mu = penaltyFactor * cj_in;
            if (jCellOut >= 0) {
                double cj_out = cellMetric[jCellOut];
                mu = Math.Max(mu, penaltyFactor * cj_out);
            }
            return mu;
        }

        /// <summary>
        /// <see cref="TermActivationFlags"/>
        /// </summary>
        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }
        /// <summary>
        /// <see cref="TermActivationFlags"/>
        /// </summary>
        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }
        /// <summary>
        /// <see cref="TermActivationFlags"/>
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV);
            }
        }


    }
}
