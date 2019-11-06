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
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.XheatCommon {

    /// <summary>
    /// 
    /// </summary>
    public class ConductivityInSpeciesBulk : swipConductivity, ISpeciesFilter {

        /// <summary>
        /// different implementations for the conductivity part (laplace operator) of the heat equation 
        /// </summary>
        public enum ConductivityMode {

            /// <summary>
            /// direct discretization of the laplace operator via symmetric interior penalty
            /// </summary>
            SIP,

            /// <summary>
            /// splitting into two first order differential equations, explicit computation of the heat flux
            /// </summary>
            LDG,

            /// <summary>
            /// splitting into two first order differential equations, explicit computation of the heat flux
            /// additional stabilization via penalty terms
            /// </summary>
            LDGstabi
        }


        public ConductivityInSpeciesBulk(double penalty, double sw, ThermalMultiphaseBoundaryCondMap bcMap, int D,
            string spcName, SpeciesId spcId, double _kA, double _kB)
            : base(penalty, D, bcMap) {

            base.m_alpha = sw;
            this.m_bcMap = bcMap;

            this.m_spcId = spcId;

            switch (spcName) {
                case "A": currentk = _kA; complementk = _kB; break;
                case "B": currentk = _kB; complementk = _kA; break;
                default: throw new ArgumentException("Unknown species.");
            }

            double muFactor = Math.Max(currentk, complementk) / currentk;
            base.m_penalty_base = penalty * muFactor;

            base.tempFunction = this.m_bcMap.bndFunction[VariableNames.Temperature + "#" + spcName];
            base.fluxFunction = D.ForLoop(d => bcMap.bndFunction[VariableNames.HeatFluxVectorComponent(d) + "#" + spcName]);

        }


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        double currentk = double.NaN;
        double complementk = double.NaN;


        ThermalMultiphaseBoundaryCondMap m_bcMap;


        protected override double Conductivity(double[] Parameters) {
            return currentk;
        }

    }


    /// <summary>
    /// 
    /// </summary>
    public class HeatFluxDivergenceInSpeciesBulk : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, ISpeciesFilter {

        int m_D;

        /// <summary>
        /// see <see cref="BoundaryCondMap{BCType}.EdgeTag2Type"/>;
        /// </summary>
        protected ThermalBcType[] EdgeTag2Type;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: spatial dimension <br/>
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[][] fluxFunction;


        public HeatFluxDivergenceInSpeciesBulk(int D, ThermalMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId) {

            this.m_D = D;

            this.m_spcId = spcId;
            //this.ksqrt = Math.Sqrt(_k);

            fluxFunction = D.ForLoop(d => bcMap.bndFunction[VariableNames.HeatFluxVectorComponent(d) + "#" + spcName]);
            EdgeTag2Type = bcMap.EdgeTag2Type;
        }


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.HeatFluxVector(m_D);
            }
        }


        public IList<string> ParameterOrdering {
            get {
                return new string[0];
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxGradV;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V; 
            }
        }


        public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0;
            for (int d = 0; d < cpv.D; d++)
                Acc += U[d] * GradV[d];

            return -Acc;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, 
            double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            double Acc = 0.0;

            for (int d = 0; d < m_D; d++) {
                Acc += 0.5 * (_uIN[d] + _uOUT[d]) * inp.Normale[d];
            }

            return Acc * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Acc = 0.0;

            ThermalBcType edgType = this.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case ThermalBcType.ConstantTemperature: {

                        for (int d = 0; d < m_D; d++) {
                            Acc += (_uA[d]) * _vA * inp.Normale[d];
                        }

                        break;
                    }
                case ThermalBcType.ZeroGradient: {

                        //for (int d = 0; d < m_D; d++) {
                        //    Acc += (_uA[d] - 0.0) * _vA * inp.Normale[d];
                        //}

                        break;
                    }
                case ThermalBcType.ConstantHeatFlux: {

                        for (int d = 0; d < inp.D; d++) {
                            double g_D = this.fluxFunction[d][inp.EdgeTag](inp.X, inp.time);

                            Acc += (g_D) * (_vA) * inp.Normale[d];
                        }

                        break;
                    }
                default:
                    throw new NotImplementedException();
            }

            return Acc;
        }


        //public void CoefficientUpdate(CoefficientSet cs, Int32[] DomainDGdeg, Int32 TestDGdeg) {
        //    throw new NotImplementedException();
        //}

    }

    /// <summary>
    /// 
    /// </summary>
    public class AuxiliaryStabilizationForm : IEdgeForm, ISpeciesFilter {


        int m_D;

        /// <summary>
        /// see <see cref="BoundaryCondMap{BCType}.EdgeTag2Type"/>;
        /// </summary>
        protected ThermalBcType[] EdgeTag2Type;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: spatial dimension <br/>
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[][] fluxFunction;


        public AuxiliaryStabilizationForm(int _D, ThermalMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId) {

            this.m_D = _D;

            EdgeTag2Type = bcMap.EdgeTag2Type;
            fluxFunction = m_D.ForLoop(d => bcMap.bndFunction[VariableNames.HeatFluxVectorComponent(d) + "#" + spcName]);

            this.m_spcId = spcId;
            //this.ksqrt = Math.Sqrt(_k);
        }


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V; ;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.HeatFluxVector(m_D);
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT,
            double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            double Acc = 0.0;

            for (int d = 0; d < m_D; d++) {
                Acc += (_uIN[d] - _uOUT[d]) * inp.Normale[d];
            }

            return -Acc * (_vIN - _vOUT);
        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {


            double Acc = 0.0;

            ThermalBcType edgType = this.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case ThermalBcType.ConstantTemperature: {

                        break;
                    }
                case ThermalBcType.ZeroGradient: {

                        for (int d = 0; d < m_D; d++) {
                            Acc += (_uA[d] - 0.0) * inp.Normale[d];
                        }

                        break;
                    }
                case ThermalBcType.ConstantHeatFlux: {

                        for (int d = 0; d < m_D; d++) {
                            double gD = fluxFunction[d][inp.EdgeTag](inp.X, inp.time);
                            Acc += (_uA[d] - gD) * inp.Normale[d];
                        }

                        break;
                    }
                default:
                    throw new NotImplementedException();
            }

            return -2.0 * Acc * _vA;
        }


    }



    /// <summary>
    /// Volume integral of identity part of auxiliary heat flux .
    /// </summary>
    public class AuxiliaryHeatFlux_Identity : IVolumeForm, IEquationComponent, ISpeciesFilter {

        private int component;

        /// <summary>
        /// Initialize identity
        /// </summary>
        public AuxiliaryHeatFlux_Identity(int component, SpeciesId spcId) {
            this.component = component;

            this.m_spcId = spcId;
        }

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        /// <summary>
        /// Choosing the required terms (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }


        /// <summary>
        /// Ordering of the dependencies
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.HeatFluxVectorComponent(component) };
            }
        }

        /// <summary>
        /// Ordering of the parameters - null at identity part
        /// </summary>
        public IList<string> ParameterOrdering { get; }

        /// <summary>
        /// Calculating the integral of the volume part
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            return U[0] * V;

        }

    }


    /// <summary>
    /// 
    /// </summary>
    public class TemperatureGradientInSpeciesBulk : LinearFlux, ISpeciesFilter {

        int m_d;    // component

        /// <summary>
        /// see <see cref="BoundaryCondMap{BCType}.EdgeTag2Type"/>;
        /// </summary>
        protected ThermalBcType[] EdgeTag2Type;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: edge tag
        /// </summary>
        protected Func<double[], double, double>[] tempFunction;


        public TemperatureGradientInSpeciesBulk(int _d, ThermalMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId, double _k) {

            this.m_d = _d;

            this.m_spcId = spcId;
            this.k = _k;

            tempFunction = bcMap.bndFunction[VariableNames.Temperature + "#" + spcName];
            EdgeTag2Type = bcMap.EdgeTag2Type;
        }

        double k;

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        public override IList<String> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }


        //public IList<String> ParameterOrdering {
        //    get {
        //        return new string[0];
        //    }
        //}


        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            int D = output.Length;
            Array.Clear(output, 0, D);
            output[m_d] = k * U[0];
        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            return 0.5 * k * (Uin[0] + Uout[0]) * inp.Normale[m_d];
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            ThermalBcType edgType = EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case ThermalBcType.ConstantTemperature:

                    double T_D = tempFunction[inp.EdgeTag](inp.X, inp.time);

                    return k * T_D * inp.Normale[m_d];

                case ThermalBcType.ZeroGradient:
                case ThermalBcType.ConstantHeatFlux:

                    return k * Uin[0] * inp.Normale[m_d];
                default:
                    throw new NotImplementedException();
            }
        }



        //public void CoefficientUpdate(CoefficientSet cs, Int32[] DomainDGdeg, Int32 TestDGdeg) {
        //    throw new NotImplementedException();
        //}

    }


    /// <summary>
    /// 
    /// </summary>
    public class TemperatureStabilizationForm : IEdgeForm, ISpeciesFilter {


        int m_d;

        /// <summary>
        /// see <see cref="BoundaryCondMap{BCType}.EdgeTag2Type"/>;
        /// </summary>
        protected ThermalBcType[] EdgeTag2Type;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: edge tag
        /// </summary>
        protected Func<double[], double, double>[] tempFunction;


        public TemperatureStabilizationForm(int _d, ThermalMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId) {

            this.m_d = _d;

            EdgeTag2Type = bcMap.EdgeTag2Type;
            tempFunction = bcMap.bndFunction[VariableNames.Temperature + "#" + spcName];

            this.m_spcId = spcId;
            //this.ksqrt = Math.Sqrt(_k);
        }


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V; ;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT,
            double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            return (_uIN[0] - _uOUT[0]) * (_vIN - _vOUT) * inp.Normale[m_d];
        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {


            ThermalBcType edgType = this.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case ThermalBcType.ConstantTemperature: {

                        double T_D = tempFunction[inp.EdgeTag](inp.X, inp.time);

                        return 2.0 * (_uA[0] - T_D) * (_vA) * inp.Normale[m_d];
                    }
                case ThermalBcType.ZeroGradient:
                case ThermalBcType.ConstantHeatFlux: {

                        return 0.0; // -(_uA[0]) * (_vA) * alpha;
                    }
                default:
                    throw new NotImplementedException();
            }
        }


    }

}
