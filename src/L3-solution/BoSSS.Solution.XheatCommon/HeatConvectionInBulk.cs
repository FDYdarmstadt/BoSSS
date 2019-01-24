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
using ilPSP.Utils;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;


namespace BoSSS.Solution.XheatCommon {


    public class HeatConvectionInBulk : LinearizedHeatConvection, IEquationComponentSpeciesNotification {


        public HeatConvectionInBulk(int SpatDim, ThermalMultiphaseBoundaryCondMap _bcmap, double _capA, double _capB, double _LFFA, double _LFFB, LevelSetTracker _lsTrk) 
            : base(SpatDim, _bcmap) {

            capA = _capA;
            capB = _capB;
            //varMode = _varMode;
            this.lsTrk = _lsTrk;
            this.LFFA = _LFFA;
            this.LFFB = _LFFB;
            this.m_bcmap = _bcmap;
            base.tempFunction = null;
        }

        ThermalMultiphaseBoundaryCondMap m_bcmap;
        LevelSetTracker lsTrk;

        double LFFA;
        double LFFB;

        double capA;
        double capB;
        protected double cap;

        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": this.cap = this.capA; base.LaxFriedrichsSchemeSwitch = LFFA; this.SetBndfunc("A"); break;
                case "B": this.cap = this.capB; base.LaxFriedrichsSchemeSwitch = LFFB; this.SetBndfunc("B"); break;
                default: throw new ArgumentException("Unknown species.");
            }
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(SpcId).VolumeMask.GetBitMaskWithExternal();
        }

        void SetBndfunc(string S) {
            int SpatDim = base.m_SpatialDimension;
            base.tempFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            base.tempFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Temperature + "#" + S], 0);
        }

        protected System.Collections.BitArray SubGrdMask;



        internal double IEF(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return this.InnerEdgeFlux(ref inp, Uin, Uout);
        }

        protected bool basecall = false;

        protected override double InnerEdgeFlux(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            if(basecall) {
                return base.InnerEdgeFlux(ref inp, Uin, Uout);
            } else {

                double UinBkUp = Uin[0];
                double UoutBkUp = Uout[0];
                double[] InParamsBkup = inp.Parameters_IN;
                double[] OutParamsBkup = inp.Parameters_OUT;


                // subgrid boundary handling
                // -------------------------

                if(inp.iEdge >= 0 && inp.jCellOut >= 0) {

                    bool CellIn = SubGrdMask[inp.jCellIn];
                    bool CellOut = SubGrdMask[inp.jCellOut];
                    Debug.Assert(CellIn || CellOut, "at least one cell must be in the subgrid!");

                    if(CellOut == true && CellIn == false) {
                        // IN-cell is outside of subgrid: extrapolate from OUT-cell!
                        Uin[0] = Uout[0];
                        inp.Parameters_IN = inp.Parameters_OUT.CloneAs();

                    }
                    if(CellIn == true && CellOut == false) {
                        // ... and vice-versa
                        Uout[0] = Uin[0];
                        inp.Parameters_OUT = inp.Parameters_IN.CloneAs();
                    }
                }

                // evaluate flux function
                // ----------------------

                var flx = base.InnerEdgeFlux(ref inp, Uin, Uout);
                flx *= cap;

                // cleanup mess and return
                // -----------------------

                Uout[0] = UoutBkUp;
                Uin[0] = UinBkUp;
                inp.Parameters_IN = InParamsBkup;
                inp.Parameters_OUT = OutParamsBkup;

                return flx;
            }
        }


        protected override double BorderEdgeFlux(ref BoSSS.Foundation.CommonParamsBnd inp, double[] Uin) {

            this.basecall = true;
            double flx = base.BorderEdgeFlux(ref inp, Uin);
            this.basecall = false;

            flx *= cap;

            return flx;
        }


        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            base.Flux(ref inp, U, output);
            output.ScaleV(cap);
        }



    }


    public class LinearizedHeatConvection : LinearFlux {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        ThermalBoundaryCondMap m_bcmap;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] tempFunction;


        public LinearizedHeatConvection(int SpatDim, ThermalBoundaryCondMap _bcmap) {
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;

            tempFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            tempFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Temperature], 0);

        }

        /// <summary>
        /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
        /// </summary>
        protected double LaxFriedrichsSchemeSwitch = 1.0;


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_SpatialDimension), (new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0, m_SpatialDimension));
            }
        }


        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            double rhoIn = 1.0;
            double rhoOut = 1.0;

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
            r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
            if(m_SpatialDimension == 3) {
                r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normale[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normale[2];
            }

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_SpatialDimension];
            double[] VelocityMeanOut = new double[m_SpatialDimension];
            for(int d = 0; d < m_SpatialDimension; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, true);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, true);

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = Uin[0] - Uout[0];

            r += Lambda * uJump * LaxFriedrichsSchemeSwitch;

            r *= 0.5;
            return r;
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

            ThermalBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch(edgeType) {
                case ThermalBcType.ConstantTemperature: {

                        double r = 0.0;

                        // Setup params
                        // ============
                        Foundation.CommonParams inp2;
                        inp2.GridDat = inp.GridDat;
                        inp2.Normale = inp.Normale;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;

                        // Dirichlet value for velocity
                        // ============================
                        double Uout = tempFunction[inp.EdgeTag, 0](inp.X, inp.time);


                        // Specify Parameters_OUT
                        // ======================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        // Outer values for Velocity and VelocityMean
                        for(int j = 0; j < m_SpatialDimension; j++) {

                            inp2.Parameters_OUT[j] = tempFunction[inp.EdgeTag, 0](inp.X, inp.time);

                            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
                            inp2.Parameters_OUT[m_SpatialDimension + j] = 0.0;
                        }

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================
                        r = InnerEdgeFlux(ref inp2, Uin, new double[] { Uout });

                        return r;

                    }
                case ThermalBcType.ZeroGradient:
                case ThermalBcType.ConstantHeatFlux: {

                        double r = 0.0;
                        double u1, u2, u3 = 0, u_d;

                        u_d = Uin[0];
                        u1 = inp.Parameters_IN[0];
                        u2 = inp.Parameters_IN[1];
                        if(m_SpatialDimension == 3)
                            u3 = inp.Parameters_IN[2];

                        r += u_d * (u1 * inp.Normale[0] + u2 * inp.Normale[1]);
                        if(m_SpatialDimension == 3) {
                            r += u_d * u3 * inp.Normale[2];
                        }

                        return r;
                    }
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }


        }


        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            output[0] = U[0] * inp.Parameters[0];
            output[1] = U[0] * inp.Parameters[1];
            if(m_SpatialDimension == 3) {
                output[2] = U[0] * inp.Parameters[2];
            }

        }



    }


}
