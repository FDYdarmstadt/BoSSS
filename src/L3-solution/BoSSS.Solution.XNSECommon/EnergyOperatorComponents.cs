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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.XNSECommon.Operator.Energy {



    public class EnergyBoundaryCondMap : BoundaryCondMap<IncompressibleBcType> {

        static string[] BndFunctions(IGridData g) {
            int D = g.SpatialDimension;

            return ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure, "KineticEnergy");

        }

        /// <summary>
        /// ctor
        /// </summary>
        protected EnergyBoundaryCondMap(IGridData f, IDictionary<string, AppControl.BoundaryValueCollection> b, string[] BndFuncName)
            : base(f, b, BndFuncName) {
        }

        /// <summary>
        /// ctor
        /// </summary>
        public EnergyBoundaryCondMap(IGridData f, IDictionary<string, AppControl.BoundaryValueCollection> b)
            : base(f, b, BndFunctions(f)) {
        }

    }


    public class EnergyMultiphaseBoundaryCondMap : EnergyBoundaryCondMap {

        static string[] BndFunctions(IGridData g, string[] SpeciesNames) {
            int D = g.SpatialDimension;
            List<string> scalarFields = new List<string>();

            foreach(var S in SpeciesNames) {
                for(int d = 0; d < D; d++) {
                    scalarFields.Add(VariableNames.Velocity_d(d) + "#" + S);
                }
                scalarFields.Add(VariableNames.Pressure + "#" + S);
                scalarFields.Add("KineticEnergy#" + S);
            }

            return scalarFields.ToArray();
        }


        public EnergyMultiphaseBoundaryCondMap(IGridData f, IDictionary<string, BoSSS.Solution.Control.AppControl.BoundaryValueCollection> b, string[] SpeciesNames)
           : base(f, b, BndFunctions(f, SpeciesNames)) //
        {
            string S0 = "#" + SpeciesNames[0];

            int D = f.SpatialDimension;
            for(int d = 0; d < D; d++) {
                base.bndFunction.Add(VariableNames.Velocity_d(d), base.bndFunction[VariableNames.Velocity_d(d) + S0]);
            }
            base.bndFunction.Add(VariableNames.Pressure, base.bndFunction[VariableNames.Pressure + S0]);
            base.bndFunction.Add("KineticEnergy", base.bndFunction["KineticEnergy" + S0]);
        }

    }


    // ================
    // convection terms
    // ================

    public abstract class LinearizedConvection : LinearFlux {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        EnergyBoundaryCondMap m_bcmap;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelFunction;

        protected Func<double[], double, double>[] ScalarFunction;


        public LinearizedConvection(int SpatDim, EnergyBoundaryCondMap _bcmap) {
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;

            VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for(int d = 0; d < m_SpatialDimension; d++)
                VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d)], d);

            //ScalarFunction = m_bcmap.bndFunction["KineticEnergy"];

        }

        /// <summary>
        /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
        /// </summary>
        protected double LaxFriedrichsSchemeSwitch = 1.0;


        //public override IList<string> ArgumentOrdering {
        //    get {
        //        return new string[] { "KineticEnergy" };
        //    }
        //}

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

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch(edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet: {

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

                        // Specify Parameters_OUT
                        // ======================
                        inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];


                        // Dirichlet value for scalar
                        double Uout = ScalarFunction[inp.EdgeTag](inp.X, inp.time);


                        // Outer values for Velocity and VelocityMean
                        for(int j = 0; j < m_SpatialDimension; j++) {

                            inp2.Parameters_OUT[j] = inp2.Parameters_IN[j]; // VelFunction[inp.EdgeTag, j](inp.X, inp.time);

                            // VelocityMeanOut = VelocityMeanIn
                            inp2.Parameters_OUT[m_SpatialDimension + j] = inp.Parameters_IN[m_SpatialDimension + j]; ;
                        }

                        // Calculate BorderEdgeFlux as InnerEdgeFlux
                        // =========================================
                        r = InnerEdgeFlux(ref inp2, Uin, new double[] { Uout });

                        return r;

                    }
                case IncompressibleBcType.Pressure_Outlet: {

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
            for(int d = 0; d < m_SpatialDimension; d++)
                output[d] = U[0] * inp.Parameters[d];
        }

    }


    public class KineticEnergyConvectionInBulk : LinearizedConvection, IEquationComponentSpeciesNotification {


        public KineticEnergyConvectionInBulk(int SpatDim, EnergyMultiphaseBoundaryCondMap _bcmap, double _rhoA, double _rhoB, double _LFFA, double _LFFB, LevelSetTracker _lsTrk)
            : base(SpatDim, _bcmap) {

            rhoA = _rhoA;
            rhoB = _rhoB;
            //varMode = _varMode;
            this.lsTrk = _lsTrk;
            this.LFFA = _LFFA;
            this.LFFB = _LFFB;
            this.m_bcmap = _bcmap;
            base.VelFunction = null;
            base.ScalarFunction = null;
        }

        EnergyMultiphaseBoundaryCondMap m_bcmap;
        LevelSetTracker lsTrk;

        double LFFA;
        double LFFB;

        double rhoA;
        double rhoB;
        protected double rho;

        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": this.rho = this.rhoA; base.LaxFriedrichsSchemeSwitch = LFFA; this.SetBndfunc("A"); break;
                case "B": this.rho = this.rhoB; base.LaxFriedrichsSchemeSwitch = LFFB; this.SetBndfunc("B"); break;
                default: throw new ArgumentException("Unknown species.");
            }
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(SpcId).VolumeMask.GetBitMaskWithExternal();
        }

        void SetBndfunc(string S) {
            int SpatDim = base.m_SpatialDimension;
            base.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for(int d = 0; d < SpatDim; d++)
                base.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + S], d);

            base.ScalarFunction = m_bcmap.bndFunction["KineticEnergy#" + S];
        }

        protected System.Collections.BitArray SubGrdMask;


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { "KineticEnergy" };
            }
        }


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
                flx *= rho;

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

        public KineticEnergyConvectionAtLevelSet(int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _LFFA, double _LFFB, bool _MaterialInterface, EnergyMultiphaseBoundaryCondMap _bcmap, bool _movingmesh) {
            m_D = _D;

            rhoA = _rhoA;
            rhoB = _rhoB;
            m_LsTrk = LsTrk;

            MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new KineticEnergyConvectionInBulk(_D, _bcmap, _rhoA, _rhoB, _LFFA, double.NaN, LsTrk);
            NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            PosFlux = new KineticEnergyConvectionInBulk(_D, _bcmap, _rhoA, _rhoB, double.NaN, _LFFB, LsTrk);
            PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));

        }

        bool MaterialInterface;
        double rhoA;
        double rhoB;
        int m_D;

        // Use Fluxes as in Bulk Convection
        KineticEnergyConvectionInBulk NegFlux;
        KineticEnergyConvectionInBulk PosFlux;



        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
            if(this.MaterialInterface) {

                U_NegFict = U_Pos;
                U_PosFict = U_Neg;

            } else {
                throw new NotImplementedException();
            }
        }


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

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }

            if(movingmesh)
                return 0.0;
            else
                return FlxNeg * v_Neg - FlxPos * v_Pos;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { "KineticEnergy" };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), (new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0, m_D));
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


    // =============
    // laplace terms
    // =============

    public class swipViscosity : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IEquationComponentCoefficient {


        /// <summary>
        /// a multiplier which is applied to everything but the penalty terms
        /// </summary>
        protected double m_alpha = 1.0;

        /// <summary>
        /// see <see cref="BoundaryCondMap{BCType}.EdgeTag2Type"/>;
        /// </summary>
        protected IncompressibleBcType[] EdgeTag2Type;

        /// <summary>
        /// spatial dimension
        /// </summary>
        protected int m_D;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: spatial dimension <br/>
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[][] kinFunction;


        public swipViscosity(double _penaltyBase, int D, EnergyBoundaryCondMap bcmap) {

            this.m_penalty_base = _penaltyBase;
            this.m_D = D;

            kinFunction = D.ForLoop(d => bcmap.bndFunction["KineticEnergy"]);
            EdgeTag2Type = bcmap.EdgeTag2Type;
        }


        /// <summary>
        /// in the case of constant conductivity, the value of the conductivity
        /// </summary>
        double m_constantViscosityValue = double.NaN;

        /// <summary>
        /// the thermal conductivity
        /// </summary>
        virtual protected double Viscosity(double[] Parameters) {
            return m_constantViscosityValue;
        }


        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { "KineticEnergy" };
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
                acc *= Viscosity(cpv.Parameters) * this.m_alpha;
            return -acc;
        }


        public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            double muB = this.Viscosity(inp.Parameters_OUT);


            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (muA * _Grad_uA[0, d] + muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            }
            Acc *= this.m_alpha;

            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term

            return -Acc;
        }

        public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            IncompressibleBcType edgType = this.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet: {
                        // inhom. Dirichlet b.c.
                        // +++++++++++++++++++++

                        double g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, 0);

                        for(int d = 0; d < inp.D; d++) {
                            double nd = inp.Normale[d];
                            Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                            Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Pressure_Dirichlet: {

                        for(int d = 0; d < inp.D; d++) {
                            double nd = inp.Normale[d];
                            Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                            //Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        //Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                default:
                    throw new NotImplementedException("ToDo");
            }

            return -Acc;
        }


        /// <summary>
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing,
        /// designed to match <see cref="BoSSS.Application.ipViscosity.TestSolution.U"/>
        /// </summary>
        public Func<int, double[], double> g_Diri_Override;

        /// <summary>
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing,
        /// designed to match <see cref="BoSSS.Application.ipViscosityTestSolution.dU"/>
        /// </summary>
        public Func<int, double[], int, double> g_Neu_Override;


        /// <summary>
        /// Dirichlet boundary value: the given temperature at the boundary.
        /// </summary>
        protected double g_Diri(double[] X, double time, int EdgeTag, int d) {
            if(this.g_Diri_Override == null) {
                Func<double[], double, double> boundVel = this.kinFunction[d][EdgeTag];
                double ret = boundVel(X, time);

                return ret;
            } else {
                return g_Diri_Override(d, X);
            }
        }

        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag, int d) {
            if(this.g_Neu_Override == null) {
                return 0.0;
            } else {
                double Acc = 0;
                for(int i = 0; i < this.m_D; i++) {
                    Acc += N[i] * g_Neu_Override(d, X, i);
                }
                return Acc;
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


    public class KineticEnergyLaplace : swipViscosity, IEquationComponentSpeciesNotification {


        public KineticEnergyLaplace(double penalty, double sw, EnergyMultiphaseBoundaryCondMap bcMap, int D, double _muA, double _muB)
            : base(penalty, D, bcMap) {
            muA = _muA;
            muB = _muB;
            base.m_alpha = sw;
            this.m_bcMap = bcMap;
            base.kinFunction = null;
            this.m_penalty = penalty;
        }


        double muA;
        double muB;

        double currentmu = double.NaN;
        double complementmu = double.NaN;


        EnergyMultiphaseBoundaryCondMap m_bcMap;

        /// <summary>
        /// multiplier for the penalty computation
        /// </summary>
        double m_penalty;


        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": currentmu = muA; complementmu = muB; SetBndfunction("A"); break;
                case "B": currentmu = muB; complementmu = muA; SetBndfunction("B"); break;
                default: throw new ArgumentException("Unknown species.");
            }

            double muFactor = Math.Max(currentmu, complementmu) / currentmu;
            base.m_penalty_base = this.m_penalty * muFactor;
        }

        void SetBndfunction(string S) {
            int D = base.m_D;
            base.kinFunction = D.ForLoop(d => this.m_bcMap.bndFunction["KineticEnergy#" + S]);
        }

        protected override double Viscosity(double[] Parameters) {
            return currentmu;
        }

    }

    
    public class KineticEnergylaplceAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public KineticEnergylaplceAtLevelSet(LevelSetTracker lstrk, double _muA, double _muB, double _penalty) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.penalty = _penalty;
            this.m_D = lstrk.GridDat.SpatialDimension;

        }

        double muA;
        double muB;

        double penalty;

        int m_D;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            int D = N.Length;
            //Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for(int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            Debug.Assert(!(double.IsInfinity(hCutCellMin) || double.IsNaN(hCutCellMin)));

            if(hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            double Ret = 0.0;

            Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret -= 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * (uA[0] - uB[0]);     // symmetry term
            Ret += (penalty / hCutCellMin) * (uA[0] - uB[0]) * (vA - vB) * (Math.Abs(muA) > Math.Abs(muB) ? muA : muB); // penalty term


            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }



        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { "KineticEnergy" }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
    }


    // ============
    // source terms
    // ============

    public class StressDivergence : LinearFlux, IEquationComponentSpeciesNotification {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;

        EnergyBoundaryCondMap bcmap;


        public StressDivergence(int SpatDim, EnergyMultiphaseBoundaryCondMap _bcmap, double _muA, double _muB) {
            m_D = SpatDim;
            bcmap = _bcmap;
            muA = _muA;
            muB = _muB;
        }

        double muA;
        double muB;

        double mu;

        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": this.mu = this.muA; break;
                case "B": this.mu = this.muB; break;
                default: throw new ArgumentException("Unknown species.");
            }
        }


        public override IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector()); //, VariableNames.Pressure);
            }
        }

        static double[,] VelociytGradient(double[] GradVelX, double[] GradVelY) {
            Debug.Assert(GradVelX.Length == 2);
            Debug.Assert(GradVelY.Length == 2);

            int D = GradVelX.Length;
            double[,] GradVel = new double[D, D];

            GradVel[0, 0] = GradVelX[0];
            GradVel[0, 1] = GradVelX[1];
            GradVel[1, 0] = GradVelY[0];
            GradVel[1, 1] = GradVelY[1];

            return GradVel;

        }


        protected override void Flux(ref CommonParamsVol inp, Double[] U, Double[] output) {

            double[] Vel = inp.Parameters.GetSubVector(0, m_D);
            double[,] GradVel = VelociytGradient(inp.Parameters.GetSubVector(m_D, m_D), inp.Parameters.GetSubVector(2 * m_D, m_D));
            //double Press = inp.Parameters[3 * m_D];

            for(int d = 0; d < m_D; d++) {
                output[d] = 0; // -Press * Vel[d];        // pressure term
                for(int dd = 0; dd < m_D; dd++) {       
                    output[d] += mu * GradVel[d, dd] * Vel[dd];      // velocity grad
                    //output[d] += mu * GradVel[dd, d] * Vel[dd];      // velocity grad transposed term
                }
            }
            output.ScaleV(-1.0);
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double[,] GradVel_IN = VelociytGradient(inp.Parameters_IN.GetSubVector(m_D, m_D), inp.Parameters_IN.GetSubVector(2 * m_D, m_D));
            //double Press_IN = inp.Parameters_IN[3 * m_D];

            double acc = 0;

            IncompressibleBcType edgType = bcmap.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Wall: {
                        for(int d = 0; d < m_D; d++) {
                            //acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                            for(int dd = 0; dd < m_D; dd++) {
                                acc += mu * (GradVel_IN[d, dd] * Vel_IN[dd]) * inp.Normale[d];
                                //acc += mu * (GradVel_IN[dd, d] * Vel_IN[dd]) * inp.Normale[d];  // transposed term
                            }
                        }
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        for(int d = 0; d < m_D; d++) {
                            //acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                            for(int dd = 0; dd < m_D; dd++) {
                                acc += mu * (GradVel_IN[d, dd] * Vel_IN[dd]) * inp.Normale[d];
                                //acc += mu * (GradVel_IN[dd, d] * Vel_IN[dd]) * inp.Normale[d];  // transposed term
                            }
                        }
                        break;
                    }
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Pressure_Dirichlet: {
                        for(int d = 0; d < m_D; d++) {
                            //acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                            for(int dd = 0; dd < m_D; dd++) {
                                acc += mu * (GradVel_IN[d, dd] * Vel_IN[dd]) * inp.Normale[d];
                                //acc += mu * (GradVel_IN[dd, d] * Vel_IN[dd]) * inp.Normale[d];  // transposed term
                            }
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
            double[,] GradVel_IN = VelociytGradient(inp.Parameters_IN.GetSubVector(m_D, m_D), inp.Parameters_IN.GetSubVector(2 * m_D, m_D));
            double[,] GradVel_OUT = VelociytGradient(inp.Parameters_OUT.GetSubVector(m_D, m_D), inp.Parameters_OUT.GetSubVector(2 * m_D, m_D));
            //double Press_IN = inp.Parameters_IN[3 * m_D];
            //double Press_OUT = inp.Parameters_OUT[3 * m_D];

            double acc = 0;

            for(int d = 0; d < m_D; d++) {
                //acc -= 0.5 * (Press_IN * Vel_IN[d] + Press_OUT * Vel_OUT[d]) * inp.Normale[d];
                for(int dd = 0; dd < m_D; dd++) {
                    acc += 0.5 * mu * (GradVel_IN[d, dd] * Vel_IN[dd] + GradVel_OUT[d, dd] * Vel_OUT[dd]) * inp.Normale[d];
                    //acc += 0.5 * mu * (GradVel_IN[dd, d] * Vel_IN[dd] + GradVel_OUT[dd, d] * Vel_OUT[dd]) * inp.Normale[d];  // transposed term
                }
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


    public class StressDivergenceAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public StressDivergenceAtLevelSet(LevelSetTracker lstrk, double _muA, double _muB) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.m_D = lstrk.GridDat.SpatialDimension;
        }

        double muA;
        double muB;

        int m_D;



        static double[,] VelociytGradient(double[] GradVelX, double[] GradVelY) {
            Debug.Assert(GradVelX.Length == 2);
            Debug.Assert(GradVelY.Length == 2);

            int D = GradVelX.Length;
            double[,] GradVel = new double[D, D];

            GradVel[0, 0] = GradVelX[0];
            GradVel[0, 1] = GradVelX[1];
            GradVel[1, 0] = GradVelY[0];
            GradVel[1, 1] = GradVelY[1];

            return GradVel;

        }


        public Double LevelSetForm(ref CommonParamsLs inp, Double[] uA, Double[] uB, Double[,] Grad_uA, Double[,] Grad_uB, Double vA, Double vB, Double[] Grad_vA, Double[] Grad_vB) {

            double[] Vel_A = inp.ParamsNeg.GetSubVector(0, m_D);
            double[] Vel_B = inp.ParamsPos.GetSubVector(0, m_D);
            double[,] GradVel_A = VelociytGradient(inp.ParamsNeg.GetSubVector(m_D, m_D), inp.ParamsNeg.GetSubVector(2 * m_D, m_D));
            double[,] GradVel_B = VelociytGradient(inp.ParamsPos.GetSubVector(m_D, m_D), inp.ParamsPos.GetSubVector(2 * m_D, m_D));
            //double p_A = inp.ParamsNeg[3 * m_D];
            //double p_B = inp.ParamsPos[3 * m_D];

            double ret = 0.0;

            for(int d = 0; d < m_D; d++) {
                //ret += 0.5 * (p_A * Vel_A[d] + p_B * Vel_B[d]) * inp.n[d];  // pressure
                for(int dd = 0; dd < m_D; dd++) {
                    ret -= 0.5 * (muA * GradVel_A[d, dd] * Vel_A[dd] + muB * GradVel_B[d, dd] * Vel_B[dd]) * inp.n[d];  // gradU
                    //ret -= 0.5 * (muA * GradVel_A[dd, d] * Vel_A[dd] + muB * GradVel_B[dd, d] * Vel_B[dd]) * inp.n[d];  // gradU transposed
                }
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

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector() ); //, VariableNames.Pressure);
            }
        }


    }


    public class SurfaceEnergy : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public SurfaceEnergy(int _D, LevelSetTracker LsTrk, double _sigma) {
            m_LsTrk = LsTrk;
            this.m_D = _D;
            this.sigma = _sigma;
        }

        int m_D;

        double sigma;


        public double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            Debug.Assert(cp.ParamsPos[0] == cp.ParamsNeg[0], "interface velocityX must be continuous across interface");
            Debug.Assert(cp.ParamsPos[1] == cp.ParamsNeg[1], "interface velocityY must be continuous across interface");
            Debug.Assert(cp.ParamsPos[2] == cp.ParamsNeg[2], "curvature must be continuous across interface");

            double curvature = cp.ParamsPos[m_D];
            double[] Vel = cp.ParamsPos.GetSubVector(0, m_D);
            double[] Normal = cp.n;

            double surfE = 0;
            for(int d = 0; d < m_D; d++) {
                surfE -= curvature * sigma * (Vel[d] * Normal[d]);
            }

            double FlxNeg = -0.5 * surfE;
            double FlxPos = +0.5 * surfE;

            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }


        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat((new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0, m_D), "Curvature");
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
            get { return TermActivationFlags.V; }
        }


    }


    public class DivergencePressureEnergy : LinearFlux {
        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;

        EnergyBoundaryCondMap bcmap;


        public DivergencePressureEnergy(int SpatDim, EnergyMultiphaseBoundaryCondMap _bcmap) {
            m_D = SpatDim;
            bcmap = _bcmap;
        }


        public override IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.Pressure);
            }
        }


        protected override void Flux(ref CommonParamsVol inp, Double[] U, Double[] output) {

            double[] Vel = inp.Parameters.GetSubVector(0, m_D);
            double Press = inp.Parameters[m_D];

            for(int d = 0; d < m_D; d++) {
                output[d] = -Press * Vel[d];        // pressure term
            }
            output.ScaleV(-1.0);
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double Press_IN = inp.Parameters_IN[m_D];

            double acc = 0;

            IncompressibleBcType edgType = bcmap.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet: {
                        for(int d = 0; d < m_D; d++) {
                            acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                        }
                        break;
                    }
                case IncompressibleBcType.Pressure_Outlet: {
                        for(int d = 0; d < m_D; d++) {
                            acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
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

            for(int d = 0; d < m_D; d++) {
                acc -= 0.5 * (Press_IN * Vel_IN[d] + Press_OUT * Vel_OUT[d]) * inp.Normale[d];
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

        LevelSetTracker m_LsTrk;

        public DivergencePressureEnergyAtLevelSet(LevelSetTracker lstrk) {
            this.m_LsTrk = lstrk;
            this.m_D = lstrk.GridDat.SpatialDimension;
        }

        int m_D;


        public Double LevelSetForm(ref CommonParamsLs inp, Double[] uA, Double[] uB, Double[,] Grad_uA, Double[,] Grad_uB, Double vA, Double vB, Double[] Grad_vA, Double[] Grad_vB) {

            double[] Vel_A = inp.ParamsNeg.GetSubVector(0, m_D);
            double[] Vel_B = inp.ParamsPos.GetSubVector(0, m_D);
            double p_A = inp.ParamsNeg[m_D];
            double p_B = inp.ParamsPos[m_D];

            double ret = 0.0;

            for(int d = 0; d < m_D; d++) {
                ret += 0.5 * (p_A * Vel_A[d] + p_B * Vel_B[d]) * inp.n[d];  // pressure
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

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.Pressure);
            }
        }
    }


    public class PressureGradientConvection : IVolumeForm {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;


        public PressureGradientConvection(int SpatDim) {
            m_D = SpatDim;
        }



        public IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), (new string[] { "PressureGradX", "PressureGradY", "PressureGradZ" }).GetSubVector(0, m_D) );
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[] PressGrad = cpv.Parameters.GetSubVector(m_D, m_D);

            double ret = 0;

            for(int d = 0; d < m_D; d++) {
                ret -= PressGrad[d] * Vel[d];
            }

            return -ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }


    public class Dissipation : IVolumeForm, IEquationComponentSpeciesNotification {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;


        public Dissipation(int SpatDim, double _muA, double _muB) {
            m_D = SpatDim;
            muA = _muA;
            muB = _muB;
        }

        double muA;
        double muB;

        double mu;

        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": this.mu = this.muA; break;
                case "B": this.mu = this.muB; break;
                default: throw new ArgumentException("Unknown species.");
            }
        }

        public IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector());
            }
        }


        static double[,] VelociytGradient(double[] GradVelX, double[] GradVelY) {
            Debug.Assert(GradVelX.Length == 2);
            Debug.Assert(GradVelY.Length == 2);

            int D = GradVelX.Length;
            double[,] GradVel = new double[D, D];

            GradVel[0, 0] = GradVelX[0];
            GradVel[0, 1] = GradVelX[1];
            GradVel[1, 0] = GradVelY[0];
            GradVel[1, 1] = GradVelY[1];

            return GradVel;

        }


        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[,] GradVel = VelociytGradient(cpv.Parameters.GetSubVector(m_D, m_D), cpv.Parameters.GetSubVector(2 * m_D, m_D));

            double ret = 0;

            for(int d = 0; d < m_D; d++) {
                for(int dd = 0; dd < m_D; dd++) {
                    ret += GradVel[d, dd] * GradVel[dd, d];
                    ret += GradVel[dd, d] * GradVel[dd, d];  // transposed term
                }
            }

            return mu * ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }


    public class PowerofGravity : IVolumeForm, IEquationComponentSpeciesNotification {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;


        public PowerofGravity(int SpatDim, double _rhoA, double _rhoB) {
            m_D = SpatDim;
            rhoA = _rhoA;
            rhoB = _rhoB;
        }

        double rhoA;
        double rhoB;

        double rho;

        public void SetParameter(String speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": this.rho = this.rhoA; break;
                case "B": this.rho = this.rhoB; break;
                default: throw new ArgumentException("Unknown species.");
            }
        }

        public IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.GravityVector(m_D));
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[] Grav = cpv.Parameters.GetSubVector(m_D, m_D);

            double ret = 0;

            for(int d = 0; d < m_D; d++) {
                ret += Grav[d] * Vel[d];
            }

            return -rho * ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }


    // =============
    // pressure term
    // =============


    //public abstract class LinearizedConvectionForParameters : LinearFlux {

    //    /// <summary>
    //    /// Spatial dimension;
    //    /// </summary>
    //    protected int m_SpatialDimension;

    //    EnergyBoundaryCondMap m_bcmap;


    //    public LinearizedConvectionForParameters(int SpatDim, EnergyBoundaryCondMap _bcmap) {
    //        m_SpatialDimension = SpatDim;
    //        m_bcmap = _bcmap;

    //    }

    //    /// <summary>
    //    /// set to 0.0 to turn the Lax-Friedrichs scheme into an central difference scheme.
    //    /// </summary>
    //    protected double LaxFriedrichsSchemeSwitch = 1.0;


    //    public override IList<string> ArgumentOrdering {
    //        get {
    //            return new string[] {  };
    //        }
    //    }

    //    public override IList<string> ParameterOrdering {
    //        get {
    //            return ArrayTools.Cat(VariableNames.VelocityVector(m_SpatialDimension), (new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0, m_SpatialDimension), VariableNames.Pressure);
    //        }
    //    }


    //    protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
    //        double r = 0.0;

    //        // Calculate central part
    //        // ======================

    //        // 2 * {u_i * u_j} * n_j,
    //        // resp. 2 * {rho * u_i * u_j} * n_j for variable density
    //        r += inp.Parameters_IN[2 * m_SpatialDimension] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
    //        r += inp.Parameters_OUT[2 * m_SpatialDimension] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
    //        if(m_SpatialDimension == 3) {
    //            r += inp.Parameters_IN[2 * m_SpatialDimension] * inp.Parameters_IN[2] * inp.Normale[2] + inp.Parameters_OUT[2 * m_SpatialDimension] * inp.Parameters_OUT[2] * inp.Normale[2];
    //        }

    //        // Calculate dissipative part
    //        // ==========================

    //        double[] VelocityMeanIn = new double[m_SpatialDimension];
    //        double[] VelocityMeanOut = new double[m_SpatialDimension];
    //        for(int d = 0; d < m_SpatialDimension; d++) {
    //            VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
    //            VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
    //        }

    //        double LambdaIn;
    //        double LambdaOut;

    //        LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, true);
    //        LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, true);

    //        double Lambda = Math.Max(LambdaIn, LambdaOut);
    //        double uJump = inp.Parameters_IN[2 * m_SpatialDimension] - inp.Parameters_OUT[2 * m_SpatialDimension];

    //        r += Lambda * uJump * LaxFriedrichsSchemeSwitch;

    //        r *= 0.5;
    //        return r;
    //    }


    //    protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

    //        IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

    //        switch(edgeType) {
    //            case IncompressibleBcType.Wall:
    //            case IncompressibleBcType.Velocity_Inlet: {

    //                    double r = 0.0;

    //                    // Setup params
    //                    // ============
    //                    Foundation.CommonParams inp2;
    //                    inp2.GridDat = inp.GridDat;
    //                    inp2.Normale = inp.Normale;
    //                    inp2.iEdge = inp.iEdge;
    //                    inp2.Parameters_IN = inp.Parameters_IN;
    //                    inp2.X = inp.X;
    //                    inp2.time = inp.time;

    //                    // Specify Parameters_OUT
    //                    // ======================
    //                    inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];


    //                    // Dirichlet value for scalar
    //                    double Uout = 0.0;


    //                    // Outer values for Velocity and VelocityMean
    //                    for(int j = 0; j < m_SpatialDimension; j++) {

    //                        inp2.Parameters_OUT[j] = inp2.Parameters_IN[j]; 

    //                        // VelocityMeanOut = VelocityMeanIn
    //                        inp2.Parameters_OUT[m_SpatialDimension + j] = inp.Parameters_IN[m_SpatialDimension + j];

    //                        // scalar parameter
    //                        inp2.Parameters_OUT[2*m_SpatialDimension] = inp2.Parameters_IN[2 * m_SpatialDimension]; 
    //                    }

    //                    // Calculate BorderEdgeFlux as InnerEdgeFlux
    //                    // =========================================
    //                    r = InnerEdgeFlux(ref inp2, Uin, new double[] { Uout });

    //                    return r;

    //                }
    //            case IncompressibleBcType.Pressure_Outlet: {

    //                    double r = 0.0;
    //                    double u1, u2, u3 = 0, u_d;

    //                    u_d = inp.Parameters_IN[2 * m_SpatialDimension]; // Uin[0];
    //                    u1 = inp.Parameters_IN[0];
    //                    u2 = inp.Parameters_IN[1];
    //                    if(m_SpatialDimension == 3)
    //                        u3 = inp.Parameters_IN[2];

    //                    r += u_d * (u1 * inp.Normale[0] + u2 * inp.Normale[1]);
    //                    if(m_SpatialDimension == 3) {
    //                        r += u_d * u3 * inp.Normale[2];
    //                    }

    //                    return r;
    //                }
    //            default:
    //                throw new NotImplementedException("Boundary condition not implemented!");
    //        }


    //    }


    //    protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
    //        for(int d = 0; d < m_SpatialDimension; d++)
    //            output[d] = inp.Parameters[2*m_SpatialDimension] * inp.Parameters[d];
    //    }


    //    override public TermActivationFlags VolTerms {
    //        get {
    //            return TermActivationFlags.GradV;
    //        }
    //    }

    //    override public TermActivationFlags BoundaryEdgeTerms {
    //        get {
    //            return TermActivationFlags.V;
    //        }
    //    }

    //    override public TermActivationFlags InnerEdgeTerms {
    //        get {
    //            return TermActivationFlags.V;
    //        }
    //    }

    //}


    //public class PressureConvectionInBulk : LinearizedConvectionForParameters, IEquationComponentSpeciesNotification {


    //    public PressureConvectionInBulk(int SpatDim, EnergyMultiphaseBoundaryCondMap _bcmap, double _LFFA, double _LFFB, LevelSetTracker _lsTrk)
    //        : base(SpatDim, _bcmap) {

    //        this.lsTrk = _lsTrk;
    //        this.LFFA = _LFFA;
    //        this.LFFB = _LFFB;
    //        this.m_bcmap = _bcmap;

    //    }

    //    EnergyMultiphaseBoundaryCondMap m_bcmap;
    //    LevelSetTracker lsTrk;

    //    double LFFA;
    //    double LFFB;

    //    public void SetParameter(String speciesName, SpeciesId SpcId) {
    //        switch(speciesName) {
    //            case "A": base.LaxFriedrichsSchemeSwitch = LFFA; break;
    //            case "B": base.LaxFriedrichsSchemeSwitch = LFFB; break;
    //            default: throw new ArgumentException("Unknown species.");
    //        }
    //        SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(SpcId).VolumeMask.GetBitMaskWithExternal();
    //    }

    //    //void SetBndfunc(string S) {
    //    //    int SpatDim = base.m_SpatialDimension;
    //    //    base.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
    //    //    for(int d = 0; d < SpatDim; d++)
    //    //        base.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + S], d);

    //    //    base.ScalarFunction = m_bcmap.bndFunction["KineticEnergy#" + S];
    //    //}

    //    protected System.Collections.BitArray SubGrdMask;


    //    //public override IList<string> ArgumentOrdering {
    //    //    get {
    //    //        return new string[] { "KineticEnergy" };
    //    //    }
    //    //}


    //    internal double IEF(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
    //        return this.InnerEdgeFlux(ref inp, Uin, Uout);
    //    }

    //    protected bool basecall = false;

    //    protected override double InnerEdgeFlux(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
    //        if(basecall) {
    //            return base.InnerEdgeFlux(ref inp, Uin, Uout);
    //        } else {

    //            //double UinBkUp = Uin[0];
    //            //double UoutBkUp = Uout[0];
    //            double[] InParamsBkup = inp.Parameters_IN;
    //            double[] OutParamsBkup = inp.Parameters_OUT;


    //            // subgrid boundary handling
    //            // -------------------------

    //            if(inp.iEdge >= 0 && inp.jCellOut >= 0) {

    //                bool CellIn = SubGrdMask[inp.jCellIn];
    //                bool CellOut = SubGrdMask[inp.jCellOut];
    //                Debug.Assert(CellIn || CellOut, "at least one cell must be in the subgrid!");

    //                if(CellOut == true && CellIn == false) {
    //                    // IN-cell is outside of subgrid: extrapolate from OUT-cell!
    //                    //Uin[0] = Uout[0];
    //                    inp.Parameters_IN = inp.Parameters_OUT.CloneAs();

    //                }
    //                if(CellIn == true && CellOut == false) {
    //                    // ... and vice-versa
    //                    //Uout[0] = Uin[0];
    //                    inp.Parameters_OUT = inp.Parameters_IN.CloneAs();
    //                }
    //            }

    //            // evaluate flux function
    //            // ----------------------

    //            var flx = base.InnerEdgeFlux(ref inp, Uin, Uout);

    //            // cleanup mess and return
    //            // -----------------------

    //            //Uout[0] = UoutBkUp;
    //            //Uin[0] = UinBkUp;
    //            inp.Parameters_IN = InParamsBkup;
    //            inp.Parameters_OUT = OutParamsBkup;

    //            return flx;
    //        }
    //    }


    //    protected override double BorderEdgeFlux(ref BoSSS.Foundation.CommonParamsBnd inp, double[] Uin) {

    //        this.basecall = true;
    //        double flx = base.BorderEdgeFlux(ref inp, Uin);
    //        this.basecall = false;

    //        return flx;
    //    }


    //    //protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
    //    //    base.Flux(ref inp, U, output);
    //    //}

    //}


}
