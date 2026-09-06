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


        public void SetParameter(String speciesName) {
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
            base.fluxFunction = D.ForLoop(d => m_bcMap.bndFunction[VariableNames.HeatFluxVectorComponent(d) + "#" + S]);
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
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: spatial dimension <br/>
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[][] fluxFunction;


        public swipConductivity(double _penaltyBase, int D, ThermalBoundaryCondMap bcmap) {

            this.m_penalty_base = _penaltyBase;
            this.m_D = D;

            tempFunction = bcmap.bndFunction[VariableNames.Temperature];
            fluxFunction = D.ForLoop(d => bcmap.bndFunction[VariableNames.HeatFluxVectorComponent(d)]);
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


        private static int _diagCondVolCount = 0;
        public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            // ABLATION TOGGLE (coordinator-requested per-piece ablation, this round): fully
            // disables THIS form's VolumeForm contribution for one specific target cell.
            // XNSFE_REPRO_DISABLE_VOL_CELL_COND is conduction-only (independent of the shared
            // XNSFE_REPRO_DISABLE_VOL_CELL toggle that also affects convection), used to
            // isolate conduction's own volume contribution via DIAG_FIXED_ABLATION.
            {
                string disCellStr = System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_VOL_CELL");
                string disCellStrCond = System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_VOL_CELL_COND");
                if ((disCellStr != null && int.TryParse(disCellStr, out int disCell) && cpv.jCell == disCell)
                    || (disCellStrCond != null && int.TryParse(disCellStrCond, out int disCellC) && cpv.jCell == disCellC))
                    return 0.0;
                // GLOBAL ablation (coordinator-requested cheap-reproduction-scan, this round):
                // unconditionally zero this form's contribution for EVERY cell, to check whether
                // the periodicity-causal M!=f mismatch survives with conduction's own volume term
                // completely removed from the operator (not just scaled down, and not just for one
                // target cell) -- distinct from XNSFE_REPRO_CONDUCTION_SCALE, which still leaves
                // conduction's stencil columns structurally present with a tiny nonzero value.
                if (System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_COND_ALL") == "1"
                    || System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_COND_VOL_ALL") == "1")
                    return 0.0;
            }
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[0, d] * GradV[d];

            if(acc != 0.0)
                acc *= Conductivity(cpv.Parameters) * this.m_alpha;
            double retVol = -acc;

            // DIAGNOSTIC (env-gated, zero cost unless set): dump every VolumeForm call touching a
            // target cell -- mirrors HeatConvectionInSpeciesBulk's DIAG_REALHEAT_VOL, for
            // independently verifying conduction's own volume (Laplacian) contribution.
            if (System.Environment.GetEnvironmentVariable("DIAG_REALCOND_VOL") == "1") {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_REALCOND_CELL");
                int targetCell = targetStr != null ? int.Parse(targetStr) : -1;
                if (targetCell < 0 || cpv.jCell == targetCell) {
                    int c = System.Threading.Interlocked.Increment(ref _diagCondVolCount);
                    if (c <= 2000000) {
                        System.Console.Error.WriteLine(string.Format(
                            "DIAG_REALCOND_VOL[{0}]: jCell={1} X=({2:0.######},{3:0.######},{4:0.######}) U={5:0.########} GradU=[{6:0.########},{7:0.########},{8:0.########}] V={9:0.########} GradV=[{10:0.########},{11:0.########},{12:0.########}] k={13:0.########} alpha={14:0.########} acc={15:0.##########} ret={16:0.##########}",
                            c, cpv.jCell, cpv.Xglobal[0], cpv.Xglobal[1], cpv.Xglobal[2],
                            U[0], GradU[0, 0], cpv.D > 1 ? GradU[0, 1] : double.NaN, cpv.D > 2 ? GradU[0, 2] : double.NaN,
                            V, GradV[0], cpv.D > 1 ? GradV[1] : double.NaN, cpv.D > 2 ? GradV[2] : double.NaN,
                            Conductivity(cpv.Parameters), this.m_alpha, acc, retVol));
                    }
                }
            }

            return retVol;
        }


        public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            // ABLATION TOGGLE (coordinator-requested per-edge ablation, this round): fully
            // disables THIS form's InnerEdgeForm contribution for one specific target edge.
            {
                string disEdgeStr = System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_EDGE_COND");
                if (disEdgeStr != null && int.TryParse(disEdgeStr, out int disEdge) && inp.iEdge == disEdge)
                    return 0.0;
                // GLOBAL ablation (coordinator-requested cheap-reproduction-scan, this round): see
                // matching comment in VolumeForm above.
                if (System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_COND_ALL") == "1"
                    || System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_COND_EDGE_ALL") == "1")
                    return 0.0;
            }
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double kA = this.Conductivity(inp.Parameters_IN);
            double kB = this.Conductivity(inp.Parameters_OUT);


            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (kA * _Grad_uA[0, d] + kB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normal[d];  // consistency term
                Acc += 0.5 * (kA * _Grad_vA[d] + kB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normal[d];  // symmetry term
            }
            Acc *= this.m_alpha;

            double muMax = (Math.Abs(kA) > Math.Abs(kB)) ? kA : kB;
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term

            double retCond = -Acc;

            // DIAGNOSTIC (env-gated, zero cost unless set): dump every InnerEdgeForm call
            // touching a target cell/edge -- traces the real production SIP conduction flux's
            // raw per-node/per-probe values, mirroring HeatConvectionInSpeciesBulk's
            // DIAG_REALHEAT_EDGE, so the two forms contributing to the shared HeatEq row can be
            // independently hand-verified and summed for a from-scratch check.
            if (System.Environment.GetEnvironmentVariable("DIAG_REALCOND_EDGE") == "1") {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_REALCOND_CELL");
                int targetCell = targetStr != null ? int.Parse(targetStr) : -1;
                string targetEdgeStr = System.Environment.GetEnvironmentVariable("DIAG_REALCOND_EDGENO");
                int targetEdgeNo = targetEdgeStr != null ? int.Parse(targetEdgeStr) : -1;
                bool edgeMatch = targetEdgeNo < 0 || inp.iEdge == targetEdgeNo;
                if (edgeMatch && (targetCell < 0 || inp.jCellIn == targetCell || inp.jCellOut == targetCell)) {
                    int c = System.Threading.Interlocked.Increment(ref _diagCondEdgeCount);
                    if (c <= 5000000) {
                        System.Console.Error.WriteLine(string.Format(
                            "DIAG_REALCOND_EDGE[{0}]: edge={1} jIn={2} jOt={3} X=({4:0.######},{5:0.######},{6:0.######}) N=({7:0.########},{8:0.########},{9:0.########}) uA={10:0.########} uB={11:0.########} GradUA=[{12:0.########},{13:0.########},{14:0.########}] GradUB=[{15:0.########},{16:0.########},{17:0.########}] vA={18:0.####} vB={19:0.####} GradVA=[{20:0.########},{21:0.########},{22:0.########}] GradVB=[{23:0.########},{24:0.########},{25:0.########}] kA={26:0.########} kB={27:0.########} pnlty={28:0.########} alpha={29:0.########} Acc={30:0.##########} ret={31:0.##########}",
                            c, inp.iEdge, inp.jCellIn, inp.jCellOut,
                            inp.X[0], inp.X[1], inp.X[2],
                            inp.Normal[0], inp.Normal[1], inp.Normal[2],
                            _uA[0], _uB[0],
                            _Grad_uA[0, 0], inp.D > 1 ? _Grad_uA[0, 1] : double.NaN, inp.D > 2 ? _Grad_uA[0, 2] : double.NaN,
                            _Grad_uB[0, 0], inp.D > 1 ? _Grad_uB[0, 1] : double.NaN, inp.D > 2 ? _Grad_uB[0, 2] : double.NaN,
                            _vA, _vB,
                            _Grad_vA[0], inp.D > 1 ? _Grad_vA[1] : double.NaN, inp.D > 2 ? _Grad_vA[2] : double.NaN,
                            _Grad_vB[0], inp.D > 1 ? _Grad_vB[1] : double.NaN, inp.D > 2 ? _Grad_vB[2] : double.NaN,
                            kA, kB, pnlty, this.m_alpha, Acc, retCond));
                    }
                }
            }

            return retCond;
        }
        private static int _diagCondEdgeCount = 0;

        public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            // GLOBAL ablation (coordinator-requested cheap-reproduction-scan, this round): see
            // matching comment in VolumeForm above -- fully removes conduction's boundary
            // (Dirichlet/Robin/etc.) contribution too, for every boundary face.
            if (System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_COND_ALL") == "1"
                || System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_COND_BND_ALL") == "1")
                return 0.0;

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
                            double nd = inp.Normal[d];
                            Acc += (kA * _Grad_uA[0, d]) * (_vA) * nd;
                            Acc += (kA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        Acc -= kA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                case ThermalBcType.ZeroGradient: {

                        for(int d = 0; d < inp.D; d++) {
                            double nd = inp.Normal[d];
                            //Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                            //Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        //Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                case ThermalBcType.ConstantHeatFlux: {

                        for(int d = 0; d < inp.D; d++) {
                            double g_D = this.g_Flux(inp.X, inp.time, d, inp.EdgeTag);

                            Acc += g_D * (_vA) * inp.Normal[d];
                            //Acc += (kA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        }
                        Acc *= this.m_alpha;

                        //Acc -= kA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                        break;
                    }
                case ThermalBcType.TemperatureSlip: {

                        double ls = Lslip[inp.jCellIn];

                        if (ls == 0.0)
                            goto case ThermalBcType.ConstantTemperature;

                        // robin b.c.
                        // +++++++++++++++++++++
                        double g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag);
                        Acc -= kA / ls * (_uA[0] - g_D) * _vA;
                        Acc *= this.m_alpha;
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
        double g_Flux(double[] X, double time, int d, int EdgeTag) {
            if(this.g_Flux_Override == null) {
                Func<double[], double, double> boundVel = this.fluxFunction[d][EdgeTag];
                double ret = boundVel(X, time);

                return ret;

            } else {
                return g_Flux_Override(X, time);
            }
        }

        /// <summary>
        /// slip-length for the thermal slip B.C.
        /// </summary>
        protected MultidimensionalArray Lslip;

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
            if (cs.UserDefinedValues != null) {
                if (cs.UserDefinedValues.Keys.Contains("ThermalSlipLengths"))
                    Lslip = (MultidimensionalArray)cs.UserDefinedValues["ThermalSlipLengths"];                
            }
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
            double penaltySizeFactor_B = jCellOut >= 0 && !cj[jCellOut].IsNaNorInf() ? 1.0 / cj[jCellOut] : 0;
            // 0.0 case:
            // there is no OUT-cell, i.e., current edge is at the boundary of the domain
            // -- or -- 
            // in 3D, with certain cut-cell quad rules (e.g. Algoim),
            // when the Level Set passes exactly through the corner of the cell,
            // it might happen that the species-volume of the OUT-cell is empty,
            // but the edge integral is only almost zero (slightly positive, e.g., weights around 10e-14 or so).
            // In such cases, the edge integral is evaluated (since it is a non-empty rule),
            // but the OUT-cell already has a NAN-length scale assigned.
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            var mu = penaltySizeFactor * m_penalty * m_penalty_base;
            if(mu.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");

            return mu;
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
