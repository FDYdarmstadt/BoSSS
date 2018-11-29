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
using System.Linq;
using System.Text;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Solution.XNSECommon.Operator.Convection {

    class ConvectionAtLevelSet_LLF : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        bool movingmesh;

        public ConvectionAtLevelSet_LLF(int _d, int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _LFFA, double _LFFB, bool _MaterialInterface, IncompressibleMultiphaseBoundaryCondMap _bcmap, bool _movingmesh) {
            m_D = _D;
            m_d = _d;
            rhoA = _rhoA;
            rhoB = _rhoB;
            m_LsTrk = LsTrk;
            //varMode = _varMode;
            MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new ConvectionInBulk_LLF(_D, _bcmap, _d, _rhoA, _rhoB, _LFFA, double.NaN, LsTrk);
            NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            PosFlux = new ConvectionInBulk_LLF(_D, _bcmap, _d, _rhoA, _rhoB, double.NaN, _LFFB, LsTrk);
            PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));

        }

        bool MaterialInterface;
        double rhoA;
        double rhoB;
        int m_D;
        int m_d;
        //EquationAndVarMode varMode;

        // Use Fluxes as in Bulk Convection
        ConvectionInBulk_LLF NegFlux;
        ConvectionInBulk_LLF PosFlux;
        
        /*
        // Flux over interface
        public override void DerivativVar_LevelSetFlux(out double FlxNeg, out double FlxPos, 
            ref CommonParamsLs  cp,
            double[] U_Neg, double[] U_Pos, double[,] GradU_Neg, double[,] GradU_Pos) {

            double[] U_NegFict, U_PosFict;
            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.ParamsNeg;
            double[] ParamsPos = cp.ParamsPos;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);
            //Flux for negativ side
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsNeg[d] * cp.n[d];
                //flx *= U_Neg[0];
                //FlxNeg = flx;
                                
                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsNeg;
                inp.Parameters_OUT = ParamsNegFict;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = base.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;


                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
            }
            // Flux for positive side
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsPos[d] * cp.n[d];
                //flx *= U_Pos[0];
                //FlxPos = flx;

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsPosFict;
                inp.Parameters_OUT = ParamsPos;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = base.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }
        }
        */

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
            //if (this.MaterialInterface) {
                //if (varMode == EquationAndVarMode.mom_p)
                //    throw new NotImplementedException();

                U_NegFict = U_Pos;
                U_PosFict = U_Neg;

            //} else {
            //    throw new NotImplementedException();
            //}
        }


        /*
        public override void PrimalVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParamsLs cp, 
            double[] U_Neg, double[] U_Pos) {
            FlxNeg = 0;
            FlxPos = 0;
        }

        public override void FluxPotential(out double G, double[] U) {
            G = 0;
        }

        public override void Nu(out double NuNeg, out double NuPos,
            ref CommonParamsLs cp) {
            NuPos = 1.0;
            NuNeg = 1.0;
        }
        */

        public double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.ParamsNeg;
            double[] ParamsPos = cp.ParamsPos;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);
            //Flux for negativ side
            double FlxNeg;
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsNeg[d] * cp.n[d];
                //flx *= U_Neg[0];
                //FlxNeg = flx;

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsNeg;
                inp.Parameters_OUT = ParamsNegFict;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;


                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
            }
            // Flux for positive side
            double FlxPos;
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsPos[d] * cp.n[d];
                //flx *= U_Pos[0];
                //FlxPos = flx;

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsPosFict;
                inp.Parameters_OUT = ParamsPos;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }

            if (movingmesh)
                return 0.0;
            else 
                return FlxNeg * v_Neg - FlxPos * v_Pos;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));  
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }
    }


    class GeneralizedConvectionAtLevelSet_DissipativePart : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public GeneralizedConvectionAtLevelSet_DissipativePart(int _d, int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _LFFA, double _LFFB, IncompressibleMultiphaseBoundaryCondMap _bcmap, 
            double _kA, double _kB, double _hVapA, double _hVapB) {

            m_D = _D;
            m_d = _d;
            rhoA = _rhoA;
            rhoB = _rhoB;
            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.hVapB = _hVapB;
            //M = _M;
            m_LsTrk = LsTrk;

            NegFlux = new GeneralizedConvectionInBulk_DissipativePart(_D, _bcmap, _d, _rhoA, _rhoB, _LFFA, double.NaN, LsTrk, _kA, kB, _hVapA, _hVapB);
            NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            PosFlux = new GeneralizedConvectionInBulk_DissipativePart(_D, _bcmap, _d, _rhoA, _rhoB, double.NaN, _LFFB, LsTrk, _kA, kB, _hVapA, _hVapB);
            PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));

        }

        double rhoA;
        double rhoB;

        double kA;
        double kB;
        double hVapA;   // for the identification of the liquid phase
        double hVapB;   // for the identification of the liquid phase

        //double M;

        int m_D;
        int m_d;

        // Use Fluxes as in Bulk Convection
        GeneralizedConvectionInBulk_DissipativePart NegFlux;
        GeneralizedConvectionInBulk_DissipativePart PosFlux;


        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }

        public Double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.ParamsNeg;
            double[] ParamsPos = cp.ParamsPos;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);
            //Flux for negativ side
            double FlxNeg;
            {
                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsNeg;
                inp.Parameters_OUT = ParamsNegFict;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;

                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
            }
            // Flux for positive side
            double FlxPos;
            {
                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsPosFict;
                inp.Parameters_OUT = ParamsPos;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;

                FlxPos = -this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }

            return FlxNeg * v_Neg + FlxPos * v_Pos;
        }


        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D), new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, m_D));
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

    }


    class GeneralizedConvectionInBulk_DissipativePart : LinearFlux, IEquationComponentSpeciesNotification {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        IncompressibleBoundaryCondMap m_bcmap;

        /// <summary>
        /// Component index of the momentum equation.
        /// </summary>
        protected int m_component;

        public GeneralizedConvectionInBulk_DissipativePart(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, int _component, double _rhoA, double _rhoB, double _LFFA, double _LFFB, LevelSetTracker _lsTrk,
            double _kA, double _kB, double _hVapA, double _hVapB) {

            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;
            m_component = _component;

            rhoA = _rhoA;
            rhoB = _rhoB;

            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.hVapB = _hVapB;

            this.lsTrk = _lsTrk;
            this.LFFA = _LFFA;
            this.LFFB = _LFFB;

            //M = _M;

        }

        LevelSetTracker lsTrk;

        double LFFA;
        double LFFB;

        double rhoA;
        double rhoB;

        double kA;
        double kB;
        double hVapA;   // for the identification of the liquid phase
        double hVapB;   // for the identification of the liquid phase

        protected double rho_in;
        protected double rho_out;

        protected double k_in;
        protected double k_out;
        protected double hVap_in;
        protected double hVap_out;

        //double M;

        private double ComputeEvaporationMass(double[] GradT_In, double[] GradT_Out, double[] n) {

            double mEvap = 0.0;

            // for testing purposes
            double prescribedVolumeFlux = 0.1;
            if(hVap_in > 0) {
                mEvap = -rho_in * prescribedVolumeFlux;
            } else {
                mEvap = rho_out * prescribedVolumeFlux;
            }

            // TODO 

            return mEvap;
        }


        /// <summary>
        /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
        /// </summary>
        protected double LaxFriedrichsSchemeSwitch = 1.0;

        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": this.rho_in = this.rhoA; this.rho_out = this.rhoB; LaxFriedrichsSchemeSwitch = LFFA;
                    this.k_in = this.kA; this.k_out = this.kB; this.hVap_in = this.hVapA; this.hVap_out = this.hVapB; break;
                case "B": this.rho_in = this.rhoB; this.rho_out = this.rhoA; LaxFriedrichsSchemeSwitch = LFFB;
                    this.k_in = this.kB; this.k_out = this.kA; this.hVap_in = this.hVapB; this.hVap_out = this.hVapA; break;
                default: throw new ArgumentException("Unknown species.");
            }
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(SpcId).VolumeMask.GetBitMaskWithExternal();
        }

        protected System.Collections.BitArray SubGrdMask;


        internal double IEF(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return this.InnerEdgeFlux(ref inp, Uin, Uout);
        }

        protected override double InnerEdgeFlux(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {

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

            double flx = 0.0;

            // Calculate central part
            // ======================

            //double rhoIn = 1.0;
            //double rhoOut = 1.0;

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            //r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
            //r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
            //if(m_SpatialDimension == 3) {
            //    r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normale[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normale[2];
            //}

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_SpatialDimension];
            double[] VelocityMeanOut = new double[m_SpatialDimension];
            for(int d = 0; d < m_SpatialDimension; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
            }

            double[] GradTempIn = new double[m_SpatialDimension];
            double[] GradTempOut = new double[m_SpatialDimension];
            for(int d = 0; d < m_SpatialDimension; d++) {
                GradTempIn[d] = inp.Parameters_IN[2 * m_SpatialDimension + d];
                GradTempOut[d] = inp.Parameters_OUT[2 * m_SpatialDimension + d];
            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, true);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, true);

            double Lambda = Math.Max(LambdaIn, LambdaOut);

            double M = ComputeEvaporationMass(GradTempIn, GradTempOut, inp.Normale);
            double uJump = -M * ((1/rho_in) - (1/rho_out)) * inp.Normale[0];

            flx += Lambda * uJump * LaxFriedrichsSchemeSwitch;

            flx *= 0.5;

            flx *= rho_in;

            // cleanup mess and return
            // -----------------------

            Uout[0] = UoutBkUp;
            Uin[0] = UinBkUp;
            inp.Parameters_IN = InParamsBkup;
            inp.Parameters_OUT = OutParamsBkup;

            return flx;

        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            return 0.0;
        }


        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            output.Clear();
        }


        /// <summary>
        /// name of the <em>d</em>-th velocity component
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[] { }; }
        }


        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_SpatialDimension), VariableNames.Velocity0MeanVector(m_SpatialDimension), new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, m_SpatialDimension));
            }
        }


    }

}
