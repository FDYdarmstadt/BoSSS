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


namespace BoSSS.Solution.XheatCommon {


    public class ConductivityInBulk : swipConductivity, IEquationComponentSpeciesNotification {


        public ConductivityInBulk(double penalty, double sw, ThermalMultiphaseBoundaryCondMap bcMap, int D, double _kA, double _kB) 
            : base(penalty, D, bcMap) {
            kA = _kA;
            kB = _kB;
            base.m_alpha = sw;
            this.m_bcMap = bcMap;
            base.tempFunction = null;
            this.m_penalty = penalty;
        }


        double kA;
        double kB;

        double currentk = double.NaN;
        double complementk = double.NaN;


        ThermalMultiphaseBoundaryCondMap m_bcMap;

        /// <summary>
        /// multiplier for the penalty computation
        /// </summary>
        double m_penalty;


        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": currentk = kA; complementk = kB; SetBndfunction("A"); break;
                case "B": currentk = kB; complementk = kA; SetBndfunction("B"); break;
                default: throw new ArgumentException("Unknown species.");
            }

            double muFactor = Math.Max(currentk, complementk) / currentk;
            base.m_penalty_base = this.m_penalty * muFactor;
        }

        void SetBndfunction(string S) {
            int D = base.m_D;
            base.tempFunction = this.m_bcMap.bndFunction[VariableNames.Temperature + "#" + S];
            base.fluxFunction = this.m_bcMap.bndFunction["HeatFlux#" + S];
        }

        protected override double Conductivity(double[] Parameters) {
            return currentk;
        }

    }


    public class swipConductivity : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IEquationComponentCoefficient {


        /// <summary>
        /// a multiplier which is applied to everything but the penalty terms
        /// </summary>
        protected double m_alpha = 1.0;

        /// <summary>
        /// see <see cref="BoundaryCondMap{BCType}.EdgeTag2Type"/>;
        /// </summary>
        protected ThermalBcType[] EdgeTag2Type;

        /// <summary>
        /// spatial dimension
        /// </summary>
        protected int m_D;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: edge tag
        /// </summary>
        protected Func<double[], double, double>[] tempFunction;

        /// <summary>
        /// flux boundary values; <br/>
        ///  - 1st index: edge tag
        /// </summary>
        protected Func<double[], double, double>[] fluxFunction;


        public swipConductivity(double _penaltyBase, int D, ThermalBoundaryCondMap bcmap) {

            this.m_penalty_base = _penaltyBase;
            this.m_D = D;

            tempFunction = bcmap.bndFunction[VariableNames.Temperature];
            fluxFunction = bcmap.bndFunction["HeatFlux"];
            EdgeTag2Type = bcmap.EdgeTag2Type;
        }


        /// <summary>
        /// in the case of constant conductivity, the value of the conductivity
        /// </summary>
        double m_constantConductivityValue = double.NaN;

        /// <summary>
        /// the thermal conductivity
        /// </summary>
        virtual protected double Conductivity(double[] Parameters) {
            return m_constantConductivityValue;            
        }


        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public virtual IList<string> ParameterOrdering {
            get {
                return new string[0];
            }
        }


        public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[0, d] * GradV[d];

            if(acc != 0.0)
                acc *= Conductivity(cpv.Parameters) * this.m_alpha;
            return -acc;
        }


        public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double kA = this.Conductivity(inp.Parameters_IN);
            double kB = this.Conductivity(inp.Parameters_OUT);


            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (kA * _Grad_uA[0, d] + kB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += 0.5 * (kA * _Grad_vA[d] + kB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            }
            Acc *= this.m_alpha;

            double muMax = (Math.Abs(kA) > Math.Abs(kB)) ? kA : kB;
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term

            return -Acc;
        }

        public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double kA = this.Conductivity(inp.Parameters_IN);
            ThermalBcType edgType = this.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case ThermalBcType.ConstantTemperature: {
                        // inhom. Dirichlet b.c.
                        // +++++++++++++++++++++

                        double g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag);

                        for(int d = 0; d < inp.D; d++) {
                            double nd = inp.Normale[d];
                            Acc += (kA * _Grad_uA[0, d]) * (_vA) * nd;
                            Acc += (kA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        Acc -= kA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                case ThermalBcType.ZeroGradient: {

                        for(int d = 0; d < inp.D; d++) {
                            double nd = inp.Normale[d];
                            //Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                            //Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        //Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                case ThermalBcType.ConstantHeatFlux: {

                        double g_D = this.g_Flux(inp.X, inp.time, inp.EdgeTag);

                        for(int d = 0; d < inp.D; d++) {
                            double nd = inp.Normale[d];
                            Acc += g_D * (_vA) * nd;
                            //Acc += (kA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        //Acc -= kA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                default:
                    throw new NotImplementedException();
            }

            return -Acc;
        }


        /// <summary>
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing,
        /// designed to match <see cref="BoSSS.Application.ipViscosity.TestSolution.U"/>
        /// </summary>
        public Func<double[], double, double> g_Diri_Override;

        /// <summary>
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing,
        /// designed to match <see cref="BoSSS.Application.ipViscosityTestSolution.dU"/>
        /// </summary>
        public Func<double[], double, double> g_Flux_Override;


        /// <summary>
        /// Dirichlet boundary value: the given temperature at the boundary.
        /// </summary>
        protected double g_Diri(double[] X, double time, int EdgeTag) {
            if(this.g_Diri_Override == null) {
                Func<double[], double, double> boundVel = this.tempFunction[EdgeTag];
                double ret = boundVel(X, time);

                return ret;
            } else {
                return g_Diri_Override(X, time);
            }
        }

        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Flux(double[] X, double time, int EdgeTag) {
            if(this.g_Flux_Override == null) {
                Func<double[], double, double> boundVel = this.fluxFunction[EdgeTag];
                double ret = boundVel(X, time);

                return ret;

            } else {
                return g_Flux_Override(X, time);
            }
        }



        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        /// <param name="cs"></param>
        /// <param name="DomainDGdeg"></param>
        /// <param name="TestDGdeg"></param>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            m_D = cs.GrdDat.SpatialDimension;
            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;

            //Lslip = (MultidimensionalArray)cs.UserDefinedValues["SlipLengths"];
        }

        /// <summary>
        /// Cell-wise length scales for the penalty computation.
        /// </summary>
        MultidimensionalArray cj;

        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double penalty(int jCellIn, int jCellOut) {

            double penaltySizeFactor_A = 1.0 / cj[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * m_penalty_base;
        }



        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }



    }

}
