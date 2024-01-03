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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator.Viscosity;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution.EnergyCommon {

    public class DivergencePressureEnergyInSpeciesBulk : LinearFlux, ISpeciesFilter {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;


        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelocFunction;

        protected Func<double[], double, double>[] PressFunction;


        public DivergencePressureEnergyInSpeciesBulk(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId) {
            m_D = SpatDim;
            m_bcMap = _bcmap;
            m_spcId = spcId;
            ValidSpecies = spcName;

            VelocFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < m_D; d++)
                VelocFunction.SetColumn(m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            PressFunction = m_bcMap.bndFunction[VariableNames.Pressure + "#" + spcName];
        }


        IncompressibleMultiphaseBoundaryCondMap m_bcMap;


        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        public override IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Pressure0);
            }
        }


        protected override void Flux(ref CommonParamsVol inp, Double[] U, Double[] output) {

            double[] Vel = inp.Parameters.GetSubVector(0, m_D);
            double Press = inp.Parameters[m_D];

            for (int d = 0; d < m_D; d++) {
                output[d] = Press * Vel[d];        // pressure term
            }
            //output.ScaleV(-1.0);
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double Press_IN = inp.Parameters_IN[m_D];

            double acc = 0;

            IncompressibleBcType edgType = m_bcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Wall: {
                        //for (int d = 0; d < m_D; d++) {
                        //    acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                        //}
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        for (int d = 0; d < m_D; d++) {
                            double VelD = VelocFunction[inp.EdgeTag,d](inp.X, inp.time);
                            acc -= Press_IN * VelD * inp.Normal[d];
                        }
                        break;
                    }
                case IncompressibleBcType.Pressure_Outlet: {
                        double pD = PressFunction[inp.EdgeTag](inp.X, inp.time);
                        for (int d = 0; d < m_D; d++) {
                            acc -= pD * Vel_IN[d] * inp.Normal[d];
                        }
                        break;
                    }
                default: {
                        throw new NotImplementedException("ToDo");
                    }
            }

            return -acc;

        }

        protected override double InnerEdgeFlux(ref CommonParams inp, Double[] Uin, Double[] Uout) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double[] Vel_OUT = inp.Parameters_OUT.GetSubVector(0, m_D);
            double Press_IN = inp.Parameters_IN[m_D];
            double Press_OUT = inp.Parameters_OUT[m_D];

            double acc = 0;

            for (int d = 0; d < m_D; d++) {
                acc -= 0.5 * (Press_IN * Vel_IN[d] + Press_OUT * Vel_OUT[d]) * inp.Normal[d];
            }

            return -acc;
        }


        override public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        override public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        override public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }
    }


    public class DivergencePressureEnergyAtLevelSet : ILevelSetForm {

        //LevelSetTracker m_LsTrk;

        public DivergencePressureEnergyAtLevelSet(int SpatialDimension) {
            //this.m_LsTrk = lstrk;
            this.m_D = SpatialDimension;
        }

        int m_D;


        public double InnerEdgeForm(ref CommonParams inp, Double[] uA, Double[] uB, Double[,] Grad_uA, Double[,] Grad_uB, Double vA, Double vB, Double[] Grad_vA, Double[] Grad_vB) {

            int D = inp.D;
            double[] Vel_A = inp.Parameters_IN.GetSubVector(0, D);
            double[] Vel_B = inp.Parameters_OUT.GetSubVector(0, D);
            double p_A = inp.Parameters_IN[D];
            double p_B = inp.Parameters_OUT[D];

            double ret = 0.0;

            for (int d = 0; d < D; d++) {
                ret += 0.5 * (p_A * Vel_A[d] + p_B * Vel_B[d]) * inp.Normal[d];  // pressure
            }

            ret *= (vA - vB);

            return ret;
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { }; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Pressure0);
            }
        }
    }


    public class PressureGradientConvection : IVolumeForm, ISpeciesFilter {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;


        public PressureGradientConvection(int SpatDim, string spcNmn, SpeciesId spcId) {
            m_D = SpatDim;
            m_spcId = spcId;
            ValidSpecies = spcNmn;
        }


        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        public IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.PressureGradient(m_D));
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[] PressGrad = cpv.Parameters.GetSubVector(m_D, m_D);

            double ret = 0;

            for (int d = 0; d < m_D; d++) {
                ret += PressGrad[d] * Vel[d];
            }

            return ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }




    public class ConvectivePressureTerm_LLF : LinearFlux, ISpeciesFilter, IEquationComponentCoefficient {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        IncompressibleMultiphaseBoundaryCondMap m_bcmap;
        LevelSetTracker lsTrk;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelFunction;

        protected Func<double[], double, double>[] PressFunction;


        public ConvectivePressureTerm_LLF(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId, double _LFF, LevelSetTracker _lsTrk) {

            this.m_SpatialDimension = SpatDim;
            this.LaxFriedrichsSchemeSwitch = _LFF;

            this.lsTrk = _lsTrk;

            this.m_bcmap = _bcmap;
            this.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                this.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            PressFunction = m_bcmap.bndFunction[VariableNames.Pressure + "#" + spcName];

            this.m_spcId = spcId;
            this.ValidSpecies = spcName;
        }


        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        /// <summary>
        /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
        /// </summary>
        protected double LaxFriedrichsSchemeSwitch = 1.0;

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_SpatialDimension), VariableNames.Velocity0MeanVector(m_SpatialDimension), VariableNames.Pressure0);
            }
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            for (int d = 0; d < m_SpatialDimension; d++)
                output[d] = inp.Parameters[2 * m_SpatialDimension] * inp.Parameters[d];
        }


        internal double IEF(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return this.InnerEdgeFlux(ref inp, Uin, Uout);
        }

        protected double InnerEdgeFlux_BaseCall(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            r += inp.Parameters_IN[2 * m_SpatialDimension] * (inp.Parameters_IN[0] * inp.Normal[0] + inp.Parameters_IN[1] * inp.Normal[1]);
            r += inp.Parameters_OUT[2 * m_SpatialDimension] * (inp.Parameters_OUT[0] * inp.Normal[0] + inp.Parameters_OUT[1] * inp.Normal[1]);
            if (m_SpatialDimension == 3) {
                r += inp.Parameters_IN[2 * m_SpatialDimension] * inp.Parameters_IN[2] * inp.Normal[2] + inp.Parameters_OUT[2 * m_SpatialDimension] * inp.Parameters_OUT[2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_SpatialDimension];
            double[] VelocityMeanOut = new double[m_SpatialDimension];
            for (int d = 0; d < m_SpatialDimension; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, true);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, true);

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = inp.Parameters_IN[2 * m_SpatialDimension] - inp.Parameters_OUT[2 * m_SpatialDimension];

            r += Lambda * uJump * LaxFriedrichsSchemeSwitch;

            r *= 0.5;
            return r;
        }


        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {

            //double UinBkUp = Uin[0];
            //double UoutBkUp = Uout[0];
            double[] InParamsBkup = inp.Parameters_IN;
            double[] OutParamsBkup = inp.Parameters_OUT;


            // subgrid boundary handling
            // -------------------------

            if (inp.iEdge >= 0 && inp.jCellOut >= 0) {

                bool CellIn = SubGrdMask[inp.jCellIn];
                bool CellOut = SubGrdMask[inp.jCellOut];
                Debug.Assert(CellIn || CellOut, "at least one cell must be in the subgrid!");

                if (CellOut == true && CellIn == false) {
                    // IN-cell is outside of subgrid: extrapolate from OUT-cell!
                    Uin[0] = Uout[0];
                    inp.Parameters_IN = inp.Parameters_OUT.CloneAs();

                }
                if (CellIn == true && CellOut == false) {
                    // ... and vice-versa
                    Uout[0] = Uin[0];
                    inp.Parameters_OUT = inp.Parameters_IN.CloneAs();
                }
            }

            // evaluate flux function
            // ----------------------

            var flx = this.InnerEdgeFlux_BaseCall(ref inp, Uin, Uout);

            // cleanup mess and return
            // -----------------------

            //Uout[0] = UoutBkUp;
            //Uin[0] = UinBkUp;
            inp.Parameters_IN = InParamsBkup;
            inp.Parameters_OUT = OutParamsBkup;

            return flx;

        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet: {

                        double r = 0.0;

                        // Setup params
                        // ============
                        CommonParams inp2 = new CommonParams();
                        inp2.GridDat = inp.GridDat;
                        inp2.Normal = inp.Normal;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;

                        // Specify Parameters_OUT
                        // ======================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                        double Uout = inp.Parameters_IN[2 * m_SpatialDimension];

                        // Outer values for Velocity and VelocityMean
                        for (int j = 0; j < m_SpatialDimension; j++) {

                            inp2.Parameters_OUT[j] = inp2.Parameters_IN[j]; 

                            // VelocityMeanOut = VelocityMeanIn
                            inp2.Parameters_OUT[m_SpatialDimension + j] = inp.Parameters_IN[m_SpatialDimension + j];
                        }

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================
                        r = InnerEdgeFlux(ref inp2, Uin, new double[] { Uout });

                        return r;
                    }
                case IncompressibleBcType.Pressure_Outlet: {

                        double r = 0.0;
                        double u1, u2, u3 = 0, u_d;

                        u_d = this.PressFunction[inp.EdgeTag](inp.X, inp.time);
                        u1 = inp.Parameters_IN[0];
                        u2 = inp.Parameters_IN[1];
                        if (m_SpatialDimension == 3)
                            u3 = inp.Parameters_IN[2];

                        r += u_d * (u1 * inp.Normal[0] + u2 * inp.Normal[1]);
                        if (m_SpatialDimension == 3) {
                            r += u_d * u3 * inp.Normal[2];
                        }

                        return r;
                    }
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }

        }


        protected System.Collections.BitArray SubGrdMask;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(m_spcId).VolumeMask.GetBitMaskWithExternal();
        }


        override public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        override public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        override public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

    }


    public class ConvectivePressureTermAtLevelSet_LLF : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        bool movingmesh;

        public ConvectivePressureTermAtLevelSet_LLF(int _D, LevelSetTracker LsTrk, double _LFFA, double _LFFB,
            bool _MaterialInterface, IncompressibleMultiphaseBoundaryCondMap _bcmap, bool _movingmesh) {

            m_D = _D;

            m_LsTrk = LsTrk;

            MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new ConvectivePressureTerm_LLF(_D, _bcmap, "A", LsTrk.GetSpeciesId("A"), _LFFA, LsTrk);
            PosFlux = new ConvectivePressureTerm_LLF(_D, _bcmap, "B", LsTrk.GetSpeciesId("B"), _LFFB, LsTrk);

        }

        bool MaterialInterface;

        int m_D;

        // Use Fluxes as in Bulk Convection
        ConvectivePressureTerm_LLF NegFlux;
        ConvectivePressureTerm_LLF PosFlux;



        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
            if (this.MaterialInterface) {

                U_NegFict = U_Pos;
                U_PosFict = U_Neg;

            } else {
                throw new NotImplementedException();
            }
        }


        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.Parameters_IN;
            double[] ParamsPos = cp.Parameters_OUT;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);
            //Flux for negativ side
            double FlxNeg;
            {

                CommonParams inp = new CommonParams(); ; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsNeg;
                inp.Parameters_OUT = ParamsNegFict;
                inp.Normal = cp.Normal;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.X;
                inp.time = cp.time;

                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
            }
            // Flux for positive side
            double FlxPos;
            {

                CommonParams inp = new CommonParams(); ; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsPosFict;
                inp.Parameters_OUT = ParamsPos;
                inp.Normal = cp.Normal;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.X;
                inp.time = cp.time;

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }

            if (movingmesh)
                return 0.0;
            else
                return FlxNeg * v_Neg - FlxPos * v_Pos;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] {  };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D), VariableNames.Pressure0);
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
    }


    public class ConvectivePressureTerm_Upwind : LinearFlux, ISpeciesFilter {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        IncompressibleMultiphaseBoundaryCondMap m_bcmap;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelFunction;

        protected Func<double[], double, double>[] PressFunction;


        public ConvectivePressureTerm_Upwind(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId) {

            this.m_SpatialDimension = SpatDim;

            this.m_bcmap = _bcmap;
            this.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                this.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            PressFunction = m_bcmap.bndFunction[VariableNames.Pressure + "#" + spcName];

            this.m_spcId = spcId;
            this.ValidSpecies = spcName;
        }


        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_SpatialDimension), VariableNames.Pressure0);
            }
        }


        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {

            for (int d = 0; d < m_SpatialDimension; d++)
                output[d] = inp.Parameters[m_SpatialDimension] * inp.Parameters[d];

        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {

            double c = 0.0;
            for (int d = 0; d < m_SpatialDimension; d++)
                c += inp.Parameters_IN[d] * inp.Normal[d];

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet: {
                        return (c * inp.Parameters_IN[m_SpatialDimension]);
                    }
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Pressure_Dirichlet: {
                        return (c * this.PressFunction[inp.EdgeTag](inp.X, inp.time));
                    }
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }
        }


        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {

            double c = 0.0;
            for (int d = 0; d < m_SpatialDimension; d++)
                c += 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normal[d];

            if (c > 0)
                return (c * inp.Parameters_IN[m_SpatialDimension]);
            else
                return (c * inp.Parameters_OUT[m_SpatialDimension]);


        }


        override public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        override public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        override public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

    }


}
