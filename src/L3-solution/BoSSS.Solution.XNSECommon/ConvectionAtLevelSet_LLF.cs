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
using System.Collections;
using BoSSS.Solution.XheatCommon;

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


    class ConvectionAtLevelSet_nonMaterialLLF : EvaporationAtLevelSet {


        public ConvectionAtLevelSet_nonMaterialLLF(int _d, int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.D = _D;
            this.m_d = _d;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.m_LsTrk = lsTrk;

            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.Rint = _Rint;

            this.Tsat = _Tsat;
            this.sigma = _sigma;
            this.pc = _pc;
        }

        int m_d;
        double rhoA;
        double rhoB;



        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double hVap = (hVapA > 0) ? hVapA : -hVapA;
            double M = qEvap / hVap;

            //Console.WriteLine("mEvap - ConvectionAtLevelSet_nonMaterialLLF: {0}", M);

            return M;

        }


        public override double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = ComputeEvaporationMass(cp.ParamsNeg.GetSubVector(2*D, D+3), cp.ParamsPos.GetSubVector(2 * D, D + 3), cp.n, cp.jCell);
            if (M == 0.0)
                return 0.0;


            double[] VelocityMeanIn = new double[D];
            double[] VelocityMeanOut = new double[D];
            for (int d = 0; d < D; d++) {
                VelocityMeanIn[d] = cp.ParamsNeg[D + d];
                VelocityMeanOut[d] = cp.ParamsPos[D + d];
            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.n, false);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.n, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);

            double uJump = -M * ((1 / rhoA) - (1 / rhoB)) * cp.n[m_d];

            double flx = Lambda * uJump * 0.8;

            return -flx * (rhoA * vA - rhoB * vB);
        }



        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D),
                    new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D), VariableNames.Temperature, "Curvature", "DisjoiningPressure");
            }
        }


    }


    class ConvectionAtLevelSet_Divergence : EvaporationAtLevelSet {


        public ConvectionAtLevelSet_Divergence(int _d, int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            double vorZeichen, bool RescaleConti, double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.D = _D;
            this.m_d = _d;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.m_LsTrk = lsTrk;

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= rhoA;
                scaleB /= rhoB;
            }

            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.Rint = _Rint;

            this.Tsat = _Tsat;
            this.sigma = _sigma;
            this.pc = _pc;
        }

        int m_d;
        double rhoA;
        double rhoB;

        double scaleA;
        double scaleB;



        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }


        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double hVap = (hVapA > 0) ? hVapA : -hVapA;
            double M = qEvap / hVap;

            //Console.WriteLine("mEvap - ConvectionAtLevelSet_Divergence: {0}", M);

            return M;

        }


        public override double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = ComputeEvaporationMass(cp.ParamsNeg.GetSubVector(2 * D, D + 3), cp.ParamsPos.GetSubVector(2 * D, D + 3), cp.n, cp.jCell);
            if (M == 0.0)
                return 0.0;


            double Ucentral = 0.0;
            for (int d = 0; d < D; d++) {
                Ucentral += 0.5 * (cp.ParamsNeg[d] + cp.ParamsPos[d]) * cp.n[d];
            }

            double uAxN = Ucentral * (-M * (1 / rhoA) * cp.n[m_d]);
            double uBxN = Ucentral * (-M * (1 / rhoB) * cp.n[m_d]);

            uAxN += -M * (1 / rhoA) * 0.5 * (U_Neg[0] + U_Pos[0]);
            uBxN += -M * (1 / rhoB) * 0.5 * (U_Neg[0] + U_Pos[0]);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict;
            //uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            uAxN_fict = uBxN;

            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict;
            //uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            uBxN_fict = uAxN;

            // compute the fluxes: note that for the continuity equation, we use not a real flux,
            // but some kind of penalization, therefore the fluxes have opposite signs!
            double FlxNeg = -Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = +Flux(uBxN_fict, uBxN);  // flux on B-side

            FlxNeg *= rhoA;
            FlxPos *= rhoB;

            double Ret = FlxNeg * vA - FlxPos * vB;

            return -Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }



        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }


        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D),
                    new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D), VariableNames.Temperature, "Curvature", "DisjoiningPressure");
            }
        }


    }




    //class ConvectionInBulk_weightedLLF : LinearFlux, IEquationComponentSpeciesNotification {

    //    /// <summary>
    //    /// Spatial dimension;
    //    /// </summary>
    //    protected int m_SpatialDimension;

    //    IncompressibleBoundaryCondMap m_bcmap;

    //    /// <summary>
    //    /// Component index of the momentum equation.
    //    /// </summary>
    //    protected int m_component;

    //    public ConvectionInBulk_weightedLLF(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, int _component, double _rhoA, double _rhoB, double _LFFA, double _LFFB, LevelSetTracker _lsTrk) {

    //        m_SpatialDimension = SpatDim;
    //        m_bcmap = _bcmap;
    //        m_component = _component;

    //        rhoA = _rhoA;
    //        rhoB = _rhoB;

    //        this.lsTrk = _lsTrk;
    //        this.LFFA = _LFFA;
    //        this.LFFB = _LFFB;

    //        //M = _M;

    //    }

    //    LevelSetTracker lsTrk;

    //    double LFFA;
    //    double LFFB;

    //    double rhoA;
    //    double rhoB;

    //    protected double rho_in;
    //    protected double rho_out;



    //    /// <summary>
    //    /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
    //    /// </summary>
    //    protected double LaxFriedrichsSchemeSwitch = 1.0;

    //    public void SetParameter(String speciesName, SpeciesId SpcId) {
    //        switch (speciesName) {
    //            case "A":
    //                this.rho_in = this.rhoA; this.rho_out = this.rhoB; LaxFriedrichsSchemeSwitch = LFFA; break;
    //            case "B":
    //                this.rho_in = this.rhoB; this.rho_out = this.rhoA; LaxFriedrichsSchemeSwitch = LFFB; break;
    //            default: throw new ArgumentException("Unknown species.");
    //        }
    //        SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(SpcId).VolumeMask.GetBitMaskWithExternal();
    //    }

    //    protected System.Collections.BitArray SubGrdMask;


    //    internal double IEF(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
    //        return this.InnerEdgeFlux(ref inp, Uin, Uout);
    //    }

    //    protected override double InnerEdgeFlux(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {

    //        double UinBkUp = Uin[0];
    //        double UoutBkUp = Uout[0];
    //        double[] InParamsBkup = inp.Parameters_IN;
    //        double[] OutParamsBkup = inp.Parameters_OUT;


    //        // subgrid boundary handling
    //        // -------------------------

    //        if (inp.iEdge >= 0 && inp.jCellOut >= 0) {

    //            bool CellIn = SubGrdMask[inp.jCellIn];
    //            bool CellOut = SubGrdMask[inp.jCellOut];
    //            Debug.Assert(CellIn || CellOut, "at least one cell must be in the subgrid!");

    //            if (CellOut == true && CellIn == false) {
    //                // IN-cell is outside of subgrid: extrapolate from OUT-cell!
    //                Uin[0] = Uout[0];
    //                inp.Parameters_IN = inp.Parameters_OUT.CloneAs();

    //            }
    //            if (CellIn == true && CellOut == false) {
    //                // ... and vice-versa
    //                Uout[0] = Uin[0];
    //                inp.Parameters_OUT = inp.Parameters_IN.CloneAs();
    //            }
    //        }

    //        // evaluate flux function
    //        // ----------------------

    //        double flx = 0.0;

    //        // Calculate central part
    //        // ======================

    //        //double rhoIn = 1.0;
    //        //double rhoOut = 1.0;

    //        //// 2 * {u_i * u_j} * n_j,
    //        //// resp. 2 * {rho * u_i * u_j} * n_j for variable density
    //        flx += rho_in * Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
    //        flx += rho_out * Uout[0] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
    //        if (m_SpatialDimension == 3) {
    //            flx += rho_in * Uin[0] * inp.Parameters_IN[2] * inp.Normale[2] - rho_out * Uout[0] * inp.Parameters_OUT[2] * inp.Normale[2];
    //        }

    //        // Calculate dissipative part
    //        // ==========================

    //        double[] VelocityMeanIn = new double[m_SpatialDimension];
    //        double[] VelocityMeanOut = new double[m_SpatialDimension];
    //        for (int d = 0; d < m_SpatialDimension; d++) {
    //            VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
    //            VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
    //        }

    //        double LambdaIn;
    //        double LambdaOut;

    //        LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, true);
    //        LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, true);

    //        LambdaIn *= rho_in;
    //        LambdaOut *= rho_out;

    //        double Lambda = Math.Max(LambdaIn, LambdaOut);

    //        double uJump = Uin[0] - Uout[0];

    //        flx += Lambda * uJump * LaxFriedrichsSchemeSwitch;

    //        flx *= 0.5;

    //        //flx *= rho_in;

    //        // cleanup mess and return
    //        // -----------------------

    //        Uout[0] = UoutBkUp;
    //        Uin[0] = UinBkUp;
    //        inp.Parameters_IN = InParamsBkup;
    //        inp.Parameters_OUT = OutParamsBkup;

    //        return flx;

    //    }

    //    protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
    //        return 0.0;
    //    }


    //    protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
    //        output.Clear();
    //    }


    //    /// <summary>
    //    /// name of the <em>d</em>-th velocity component
    //    /// </summary>
    //    public override IList<string> ArgumentOrdering {
    //        get { return new string[] { VariableNames.Velocity_d(m_component) }; }
    //    }


    //    public override IList<string> ParameterOrdering {
    //        get {
    //            return ArrayTools.Cat(VariableNames.Velocity0Vector(m_SpatialDimension), VariableNames.Velocity0MeanVector(m_SpatialDimension));
    //        }
    //    }


    //}


    class ConvectionAtLevelSet_weightedLLF : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public ConvectionAtLevelSet_weightedLLF(int _d, int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _LFFA, double _LFFB, IncompressibleMultiphaseBoundaryCondMap _bcmap, bool _movingmesh) {

            m_D = _D;
            m_d = _d;
            rhoA = _rhoA;
            rhoB = _rhoB;

            m_LsTrk = LsTrk;
            //MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;
            LFF = 0.5 * (_LFFA + _LFFB);

            //NegFlux = new ConvectionInBulk_weightedLLF(_D, _bcmap, _d, _rhoA, _rhoB, _LFFA, double.NaN, LsTrk);
            //NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            //PosFlux = new ConvectionInBulk_weightedLLF(_D, _bcmap, _d, _rhoA, _rhoB, double.NaN, _LFFB, LsTrk);
            //PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));

        }

        double rhoA;
        double rhoB;
        bool movingmesh;
        double LFF;

        int m_D;
        int m_d;

        // Use Fluxes as in Bulk Convection
        //ConvectionInBulk_weightedLLF NegFlux;
        //ConvectionInBulk_weightedLLF PosFlux;


        //void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
        //    U_NegFict = U_Pos; 
        //    U_PosFict = U_Neg; 
        //}

        public Double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {

            double UinBkUp = U_Neg[0];
            double UoutBkUp = U_Pos[0];
            double[] InParamsBkup = cp.ParamsNeg;
            double[] OutParamsBkup = cp.ParamsPos;


            // evaluate flux function
            // ----------------------

            double flx = 0.0;

            double[] Uint = new double[] { 0.0, 0.0 };

            // Calculate central part
            // ======================

            //// 2 * {u_i * u_j} * n_j,
            //// resp. 2 * {rho * u_i * u_j} * n_j for variable density
            flx += rhoA * U_Neg[0] * ((cp.ParamsNeg[0] - Uint[0]) * cp.n[0] + (cp.ParamsNeg[1] - Uint[1]) * cp.n[1]);
            flx += rhoB * U_Pos[0] * ((cp.ParamsPos[0] - Uint[0]) * cp.n[0] + (cp.ParamsPos[1] - Uint[1]) * cp.n[1]);
            //if (m_D == 3) {
            //    flx += rhoA * U_Neg[0] * cp.ParamsNeg[2] * cp.n[2] + rhoB * U_Pos[0] * cp.ParamsPos[2] * cp.n[2];
            //}

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = cp.ParamsNeg[m_D + d] - Uint[d];
                VelocityMeanOut[d] = cp.ParamsPos[m_D + d] - Uint[d];
            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.n, true);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.n, true);

            LambdaIn *= rhoA;
            LambdaOut *= rhoB;

            double Lambda = Math.Max(LambdaIn, LambdaOut);

            double uJump = U_Neg[0] - U_Pos[0];

            flx += Lambda * uJump * LFF;

            flx *= 0.5;

            //flx *= rho_in;

            // cleanup mess and return
            // -----------------------

            U_Pos[0] = UoutBkUp;
            U_Neg[0] = UinBkUp;
            cp.ParamsNeg = InParamsBkup;
            cp.ParamsPos = OutParamsBkup;

            // ====================================

            double Flx = flx * v_Neg - flx * v_Pos;

            return Flx;

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


    class ConvectionAtLevelSet_nonMaterial : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_lsTrk;

        public ConvectionAtLevelSet_nonMaterial(int _d, int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.D = _D;
            this.m_d = _d;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.m_lsTrk = lsTrk;

            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.Rint = _Rint;

            this.Tsat = _Tsat;
            this.sigma = _sigma;
            this.pc = _pc;
        }

        int D;
        int m_d;
        double rhoA;
        double rhoB;


        double kA;
        double kB;
        double hVapA;   // for the identification of the liquid phase
        double Rint;

        double Tsat;
        double sigma;
        double pc;



        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
        }


        private double ComputeEvaporationMass_Macro(double[] GradT_A, double[] GradT_B, double[] n) {

            double hVap = 0.0;
            double qEvap = 0.0;
            if (hVapA > 0) {
                hVap = hVapA;
                for (int d = 0; d < D; d++)
                    qEvap += (kA * GradT_A[d] - kB * GradT_B[d]) * n[d];
            } else {
                hVap = -hVapA;
                for (int d = 0; d < D; d++)
                    qEvap += (kB * GradT_B[d] - kA * GradT_A[d]) * n[d];
            }

            return qEvap / hVap;
        }

        private double ComputeEvaporationMass_Micro(double T_A, double T_B, double curv, double p_disp) {

            if (hVapA == 0.0)
                return 0.0;

            double pc0 = (pc < 0.0) ? sigma * curv + p_disp : pc;      // augmented capillary pressure (without nonlinear evaporative masss part)

            double TintMin = 0.0;
            double hVap = 0.0;
            double qEvap = 0.0;
            if (hVapA > 0) {
                hVap = hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rhoA)));
                if (T_A > TintMin)
                    qEvap = -(T_A - TintMin) / Rint;
            } else if (hVapA < 0) {
                hVap = -hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rhoB)));
                if (T_B > TintMin)
                    qEvap = (T_B - TintMin) / Rint;
            }

            return qEvap / hVap;
        }


        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, bool microRegion) {

            double M = 0.0;
            if (microRegion) {
                M = ComputeEvaporationMass_Micro(paramsNeg[D], paramsPos[D], paramsNeg[D + 1], paramsNeg[D + 2]);
            } else {
                M = ComputeEvaporationMass_Macro(paramsNeg.GetSubVector(0, D), paramsPos.GetSubVector(0, D), N);
            }

            return M;

        }


        public double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            //Debug.Assert(cp.ParamsPos[D + 1] == cp.ParamsNeg[D + 1], "curvature must be continuous across interface");
            //Debug.Assert(cp.ParamsPos[D + 2] == cp.ParamsNeg[D + 2], "disjoining pressure must be continuous across interface");

            //double M = ComputeEvaporationMass_Macro(cp.ParamsNeg.GetSubVector(0, D), cp.ParamsPos.GetSubVector(0, D), cp.n);
            //double M = ComputeEvaporationMass_Micro(cp.ParamsNeg[D], cp.ParamsPos[D], cp.ParamsNeg[D + 1], cp.ParamsNeg[D + 2]);
            double M = -0.1; // ComputeEvaporationMass(cp.ParamsNeg, cp.ParamsPos, cp.n, evapMicroRegion[cp.jCell]);

            double[] Uint = new double[] { 0.0, 0.0 };
            double UintxN = 0.0;

            double uAxN = 0.0;
            double uBxN = 0.0;

            // [[ rho(u*n) ]] {{u}} * {{v}}
            // ============================

            //for (int d = 0; d < D; d++) {
            //    uAxN += rhoA * cp.ParamsNeg[d] * cp.n[d];
            //    uBxN += rhoB * cp.ParamsPos[d] * cp.n[d];
            //}

            uAxN += -rhoA * UintxN;
            uBxN += -rhoB * UintxN;

            double Uaver = 0.5 * (U_Neg[0] + U_Pos[0]);

            uAxN *= Uaver;
            uBxN *= Uaver;


            // {{ rho(u*n) }} [[u]] * {{v}}
            // ============================

            double UnCentral = 0.0;
            for (int d = 0; d < D; d++) {
                UnCentral += 0.5 * (rhoA * (cp.ParamsNeg[d] - Uint[d]) + rhoB * (cp.ParamsPos[d] - Uint[d])) * cp.n[d];
            }

            uAxN += UnCentral * (0.0 - (-M * (1 / rhoA) * cp.n[m_d]));
            uBxN += UnCentral * (0.0 - (-M * (1 / rhoB) * cp.n[m_d]));


            // ====================================================================


            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict = uBxN;

            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict = uAxN;

            double FlxNeg = -Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = +Flux(uBxN_fict, uBxN);  // flux on B-side


            double Ret = FlxNeg * vA - FlxPos * vB;

            return Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }


        BitArray evapMicroRegion;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
        }


        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }


        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D),
                    new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D), VariableNames.Temperature, "Curvature", "DisjoiningPressure");
            }
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_lsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_lsTrk.GetSpeciesId("A"); }
        }

    }



}
