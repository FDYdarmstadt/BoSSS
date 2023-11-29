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


    public abstract class LinearizedConvection : LinearFlux {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelFunction;

        //protected Func<double[], double, double>[] ScalarFunction;


        public LinearizedConvection(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap) {
            m_SpatialDimension = SpatDim;
            m_bcMap = _bcmap;

            VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < m_SpatialDimension; d++)
                VelFunction.SetColumn(m_bcMap.bndFunction[VariableNames.Velocity_d(d)], d);

            //ScalarFunction = m_bcMap.bndFunction["KineticEnergy"];

        }

        /// <summary>
        /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
        /// </summary>
        protected double LaxFriedrichsSchemeSwitch = 1.0;


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.KineticEnergy };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_SpatialDimension), VariableNames.Velocity0MeanVector(m_SpatialDimension));
            }
        }


        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            return InnerEdgeFlux_impl(ref inp, Uin, Uout);  
        }

        /// <summary>
        /// implemented in a separate function to prevent from overloading, i.e. make sure that <see cref="BorderEdgeFlux"/> calls this implementation, not an overloaded one
        /// </summary>
        protected double InnerEdgeFlux_impl(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            double rhoIn = 1.0;
            double rhoOut = 1.0;

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normal[0] + inp.Parameters_IN[1] * inp.Normal[1]);
            r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normal[0] + inp.Parameters_OUT[1] * inp.Normal[1]);
            if (m_SpatialDimension == 3) {
                r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normal[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normal[2];
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
            double uJump = Uin[0] - Uout[0];

            r += Lambda * uJump * LaxFriedrichsSchemeSwitch;

            r *= 0.5;
            return r;
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

            IncompressibleBcType edgeType = m_bcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet: {

                        double r = 0.0;

                        // Setup params
                        // ============
                        CommonParams inp2 = new CommonParams(); ;
                        inp2.GridDat = inp.GridDat;
                        inp2.Normal = inp.Normal;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;

                        // Specify Parameters_OUT
                        // ======================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];


                        // Dirichlet value for scalar (kinetic energy)
                        double kinE_Diri = 0.0;
                        for (int i = 0; i < m_SpatialDimension; i++) {
                            Func<double[], double, double> boundVel = this.VelFunction[inp.EdgeTag, i];
                            kinE_Diri += boundVel(inp.X, inp.time) * boundVel(inp.X, inp.time);
                        }

                        double Uout = kinE_Diri / 2.0;


                        // Outer values for Velocity and VelocityMean
                        for (int j = 0; j < m_SpatialDimension; j++) {

                            inp2.Parameters_OUT[j] = inp2.Parameters_IN[j]; // VelFunction[inp.EdgeTag, j](inp.X, inp.time);

                            // VelocityMeanOut = VelocityMeanIn
                            inp2.Parameters_OUT[m_SpatialDimension + j] = inp.Parameters_IN[m_SpatialDimension + j];
                        }

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================
                        r = InnerEdgeFlux_impl(ref inp2, Uin, new double[] { Uout });

                        return r;

                    }
                case IncompressibleBcType.Pressure_Outlet: {

                        double r = 0.0;
                        double u1, u2, u3 = 0, u_d;

                        u_d = Uin[0];
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


        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            for (int d = 0; d < m_SpatialDimension; d++)
                output[d] = U[0] * inp.Parameters[d];
        }

    }


    public class KineticEnergyConvectionInSpeciesBulk : LinearizedConvection, ISpeciesFilter, IEquationComponentCoefficient {


        public KineticEnergyConvectionInSpeciesBulk(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId, double _rho, double _LFF, LevelSetTracker _lsTrk)
            : base(SpatDim, _bcmap) {

            this.rho = _rho;
            base.LaxFriedrichsSchemeSwitch = _LFF;

            this.m_bcMap = _bcmap;

            base.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                base.VelFunction.SetColumn(m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            this.lsTrk = _lsTrk;
            this.ValidSpecies = spcName;
            m_spcId = spcId;
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(spcId).VolumeMask.GetBitMaskWithExternal();
        }

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;
        LevelSetTracker lsTrk;


        protected double rho;

        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        protected System.Collections.BitArray SubGrdMask;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(m_spcId).VolumeMask.GetBitMaskWithExternal();
        }


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.KineticEnergy };
            }
        }


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

            var flx = base.InnerEdgeFlux(ref inp, Uin, Uout);
            flx *= rho;

            // cleanup mess and return
            // -----------------------

            Uout[0] = UoutBkUp;
            Uin[0] = UinBkUp;
            inp.Parameters_IN = InParamsBkup;
            inp.Parameters_OUT = OutParamsBkup;

            return flx;
        }
        


        protected override double BorderEdgeFlux(ref BoSSS.Foundation.CommonParamsBnd inp, double[] Uin) {
            double flx = base.BorderEdgeFlux(ref inp, Uin);
            flx *= rho;
            return flx;
        }


        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            base.Flux(ref inp, U, output);
            output.ScaleV(rho);
        }

    }


    public class KineticEnergyConvectionAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        bool movingmesh;

        public KineticEnergyConvectionAtLevelSet(int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _LFFA, double _LFFB, 
            bool _MaterialInterface, IncompressibleMultiphaseBoundaryCondMap _bcmap, bool _movingmesh) {

            m_D = _D;

            m_LsTrk = LsTrk;

            MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new KineticEnergyConvectionInSpeciesBulk(_D, _bcmap, "A", LsTrk.GetSpeciesId("A"), _rhoA, _LFFA, LsTrk);
            //NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            PosFlux = new KineticEnergyConvectionInSpeciesBulk(_D, _bcmap, "B", LsTrk.GetSpeciesId("B"), _rhoB, _LFFB, LsTrk);
            //PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));

        }

        bool MaterialInterface;

        int m_D;

        // Use Fluxes as in Bulk Convection
        KineticEnergyConvectionInSpeciesBulk NegFlux;
        KineticEnergyConvectionInSpeciesBulk PosFlux;



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

                CommonParams inp = new CommonParams(); // = default(BoSSS.Foundation.InParams);
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

                BoSSS.Foundation.CommonParams inp = new CommonParams(); ; // = default(BoSSS.Foundation.InParams);
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
                return new string[] { VariableNames.KineticEnergy };
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



    public class KineticEnergyConvectionInSpeciesBulk_Upwind : LinearFlux, ISpeciesFilter {

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


        public KineticEnergyConvectionInSpeciesBulk_Upwind(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId, double _rho) {

            this.m_SpatialDimension = SpatDim;

            this.rho = _rho;
            this.m_spcId = spcId;
            this.ValidSpecies = spcName;
            this.m_bcmap = _bcmap;

            this.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                this.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

        }

        double rho;

        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.KineticEnergy };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_SpatialDimension)); //, VariableNames.Velocity0MeanVector(m_SpatialDimension));
            }
        }


        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {

            for (int d = 0; d < m_SpatialDimension; d++)
                output[d] = rho * U[0] * inp.Parameters[d];

        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {

            double c = 0.0;
            for (int d = 0; d < m_SpatialDimension; d++)
                c += inp.Parameters_IN[d] * inp.Normal[d];

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet:  {
                        double kinE_Diri = 0.0;
                        for (int i = 0; i < m_SpatialDimension; i++) {
                            Func<double[], double, double> boundVel = this.VelFunction[inp.EdgeTag, i];
                            kinE_Diri += 0.5 * (boundVel(inp.X, inp.time) * boundVel(inp.X, inp.time));
                        }
                        return (c * rho * kinE_Diri);
                    }
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Pressure_Dirichlet: {
                        return (c * rho * Uin[0]);
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
                return (c * rho * Uin[0]);
            else
                return (c * rho * Uout[0]);


        }
    }


}
