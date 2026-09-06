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

    public class HeatConvectionInSpeciesBulk_LLF : LinearizedHeatConvection, ISpeciesFilter, IEquationComponentCoefficient {

        

        public HeatConvectionInSpeciesBulk_LLF(int SpatDim, ThermalMultiphaseBoundaryCondMap _bcmap, string spcName, double _cap, double _LFF)
            : base(SpatDim, _bcmap) {

            this.cap = _cap;
            //this.m_spcId = spcId;
            ValidSpecies = spcName;

            //this.lsTrk = _lsTrk;
            this.LFF = _LFF;
            this.m_bcmap = _bcmap;


            //base.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            //for (int d = 0; d < SpatDim; d++)
            //    base.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            base.TempFunction = m_bcmap.bndFunction[VariableNames.Temperature + "#" + spcName];

            //SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(spcId).VolumeMask.GetBitMaskWithExternal();
        }

        ThermalMultiphaseBoundaryCondMap m_bcmap;
        //LevelSetTracker lsTrk;

        double LFF;

        double cap;

        //SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
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
            flx *= cap;

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
            flx *= cap;
            return flx;
        }


        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            base.Flux(ref inp, U, output);
            output.ScaleV(cap);
        }


        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            //SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(m_spcId).VolumeMask.GetBitMaskWithExternal();
            SubGrdMask = cs.SpeciesSubGrdMask;
        }


    }

    public class HeatConvectionInSpeciesBulk_LLF_Newton : LinearizedHeatConvectionJacobi, ISpeciesFilter, IEquationComponentCoefficient {



        public HeatConvectionInSpeciesBulk_LLF_Newton(int SpatDim, ThermalMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId,
            double _cap, double _LFF, LevelSetTracker _lsTrk)
            : base(SpatDim, _bcmap) {

            this.cap = _cap;
            this.m_spcId = spcId;
            ValidSpecies = spcName;

            this.lsTrk = _lsTrk;
            this.LFF = _LFF;
            this.m_bcmap = _bcmap;


            //base.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            //for (int d = 0; d < SpatDim; d++)
            //    base.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            base.TempFunction = m_bcmap.bndFunction[VariableNames.Temperature + "#" + spcName];

            //SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(spcId).VolumeMask.GetBitMaskWithExternal();
        }

        ThermalMultiphaseBoundaryCondMap m_bcmap;
        LevelSetTracker lsTrk;

        double LFF;

        double cap;

        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        protected System.Collections.BitArray SubGrdMask;



        //internal double IEF(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
        //    return this.InnerEdgeFlux(ref inp, Uin, Uout);
        //}


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
            flx *= cap;

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
         
            flx *= cap;

            return flx;
        }


        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            base.Flux(ref inp, U, output);
            output.ScaleV(cap);
        }


        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(m_spcId).VolumeMask.GetBitMaskWithExternal();
            Scale = cs.HomotopyValue;
        }
    }

    /// <summary>
    /// Discretization of the Convective Form based on the discretization of the Hamiltonian, using a Roe-Type Flux, see (Cheng, Shu; 2006)
    /// </summary>
    public class HeatConvectionInSpeciesBulk_Hamiltonian_Newton : IVolumeForm, IEdgeForm, ISupportsJacobianComponent, ISpeciesFilter, IEquationComponentCoefficient {



        public HeatConvectionInSpeciesBulk_Hamiltonian_Newton(int SpatDim, ThermalMultiphaseBoundaryCondMap _bcmap, string spcName, double _cap, double _LFF) { 

            this.cap = _cap;
            //this.m_spcId = spcId;
            ValidSpecies = spcName;
            this.m_D = SpatDim;
            //this.lsTrk = _lsTrk;
            this.LFF = _LFF;
            this.m_bcmap = _bcmap;


            //base.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            //for (int d = 0; d < SpatDim; d++)
            //    base.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            this.TempFunction = m_bcmap.bndFunction[VariableNames.Temperature + "#" + spcName];

            //SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(spcId).VolumeMask.GetBitMaskWithExternal();
        }

        protected Func<double[], double, double>[] TempFunction;
        ThermalMultiphaseBoundaryCondMap m_bcmap;
        //LevelSetTracker lsTrk;

        double LFF;
        int m_D;
        double cap;
        double Scale = 1.0;

        //SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.UxV;

        public IList<string> ArgumentOrdering => new string[] { VariableNames.Temperature }.Cat(VariableNames.VelocityVector(m_D));

        public IList<string> ParameterOrdering => new string[] { };

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.GradUxV;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public static int _diagCoeffUpdateCount = 0;
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
            // DIAGNOSTIC (env-gated, zero cost unless set): dump every CoefficientUpdate call's
            // HomotopyValue, to check whether it's ever non-1 or updated inconsistently between
            // the evaluate-mode and linearize-mode m_AssembleMatrix calls in NonlinearSolver.cs.
            if (System.Environment.GetEnvironmentVariable("DIAG_HOMOTOPY_DUMP") == "1") {
                int c = System.Threading.Interlocked.Increment(ref _diagCoeffUpdateCount);
                System.Console.Error.WriteLine($"DIAG_HOMOTOPY[HeatConvectionInSpeciesBulk_Hamiltonian_Newton]: call#{c} species={ValidSpecies} HomotopyValue={cs.HomotopyValue:0.################} Scale={Scale:0.################}");
            }
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivVol, DivergenceDerivEdg };
        }

        public static int _diagVolCount = 0;
        // DIAGNOSTIC (coordinator-requested check-2, this round): shared output cap for
        // DIAG_METRIC_CONSISTENCY[HEAT-FORM] prints in InnerEdgeForm below.
        private static int s_metricConsistencyCountHeat = 0;
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            // ABLATION TOGGLE (coordinator-requested per-piece ablation, this round): fully
            // disables THIS form's VolumeForm contribution for one specific target cell.
            {
                string disCellStr = System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_VOL_CELL");
                if (disCellStr != null && int.TryParse(disCellStr, out int disCell) && cpv.jCell == disCell)
                    return 0.0;
            }
            double acc = 0;
            for (int d = 0; d < m_D; d++) {
                // starke Form
                acc += U[1 + d] * GradU[0, d] * V;
            }
            double ret = Scale * cap * acc;

            // DIAGNOSTIC (env-gated, zero cost unless set): dump every VolumeForm call touching a
            // target cell -- traces the REAL production convective term's raw per-node/per-probe
            // values (U/GradU/V, Scale, cap, species) directly, distinguishing the
            // linearize-mode/differentiator-driven calls (small perturbation deltas around Utmp)
            // from the raw nonlinear-evaluator calls (actual field values), by call-count and
            // magnitude, on an ACTUAL XnsfeReproTest.cs production run rather than an isolated
            // analog.
            if (System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_VOL") == "1") {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_CELL");
                int targetCell = targetStr != null ? int.Parse(targetStr) : -1;
                if (targetCell < 0 || cpv.jCell == targetCell) {
                    int c = System.Threading.Interlocked.Increment(ref _diagVolCount);
                    if (c <= 200000000) {
                        System.Console.Error.WriteLine(string.Format(
                            "DIAG_REALHEAT_VOL[{0}]: species={1} jCell={2} X=({3:0.######},{4:0.######},{5:0.######}) U=[{6}] GradU[0,:]=[{7:0.########},{8:0.########},{9:0.########}] V={10:0.########} Scale={11:0.####} cap={12:0.####} acc={13:0.##########} ret={14:0.##########}",
                            c, ValidSpecies, cpv.jCell, cpv.Xglobal[0], cpv.Xglobal[1], cpv.Xglobal[2],
                            string.Join(",", U.Select(v => v.ToString("0.########"))),
                            GradU[0, 0], m_D > 1 ? GradU[0, 1] : double.NaN, m_D > 2 ? GradU[0, 2] : double.NaN,
                            V, Scale, cap, acc, ret));
                    }
                }
            }

            return ret;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            // ABLATION TOGGLE (coordinator-requested per-edge ablation, this round): fully
            // disables THIS form's InnerEdgeForm contribution for one specific target edge,
            // same env-gated pattern used for the boundary-form ablations. Used to test edges
            // 35/38 individually.
            {
                string disEdgeStr = System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_EDGE");
                if (disEdgeStr != null && int.TryParse(disEdgeStr, out int disEdge) && inp.iEdge == disEdge)
                    return 0.0;
            }

            // DIAGNOSTIC (coordinator-requested check-2, this round): this InnerEdgeForm
            // implementation is called from BOTH the M-path (EdgeFormDifferentiator wrapping
            // it for f0/every f1 probe -- see FormDifferentiator.cs's own DIAG_METRIC_CONSISTENCY
            // dump) AND the F-path (NECQuadratureEdge.cs calling it directly, undifferentiated,
            // for f(u0)) -- so dumping inp.Normal/inp.X HERE, at G17 precision, from every call
            // touching the target cell, lets both paths' raw geometric input be compared by
            // matching physical X, without needing separate instrumentation in the two
            // orchestrator files.
            if (System.Environment.GetEnvironmentVariable("DIAG_METRIC_CONSISTENCY") == "1") {
                string mcCellStrH = System.Environment.GetEnvironmentVariable("DIAG_METRIC_CONSISTENCY_CELL");
                if (mcCellStrH != null && int.TryParse(mcCellStrH, out int mcCellH)
                    && (inp.jCellIn == mcCellH || inp.jCellOut == mcCellH)
                    && s_metricConsistencyCountHeat < 400) {
                    System.Threading.Interlocked.Increment(ref s_metricConsistencyCountHeat);
                    var xLocalHeat = inp.X;
                    Console.Error.WriteLine(string.Format(
                        "DIAG_METRIC_CONSISTENCY[HEAT-FORM]: edge={0} jIn={1} jOt={2} X=({3:G17}) Normal=({4:G17},{5:G17},{6:G17})",
                        inp.iEdge, inp.jCellIn, inp.jCellOut,
                        string.Join(",", Enumerable.Range(0, xLocalHeat.Dim).Select(d => xLocalHeat[d].ToString("G17"))),
                        inp.Normal[0], inp.Normal[1], inp.Normal[2]));
                }
            }

            double flx = 0.0;

            //===========================================================================================================
            //===========================================================================================================
            // First variant, using central flux for temperature            
            /*
            double FlxNeg = _uIN[0] * (_uIN[1] * inp.Normal[0] + _uIN[2] * inp.Normal[1]);
            double FlxPos = _uOUT[0] * (_uOUT[1] * inp.Normal[0] + _uOUT[2] * inp.Normal[1]);
            if (m_D == 3) {
                FlxNeg += _uIN[0] * _uIN[3] * inp.Normal[2];
                FlxPos += _uOUT[0] * _uOUT[3] * inp.Normal[2];
            }

            // Term from partial integration back to strong form
            double sflx = FlxNeg * _vIN - FlxPos * _vOUT;
            //funktioniert mit starker Form
            flx = (0.5 * (_uIN[0] + _uOUT[0])) * ((_uIN[1] * inp.Normal[0] + _uIN[2] * inp.Normal[1]) * _vIN - (_uOUT[1] * inp.Normal[0] + _uOUT[2] * inp.Normal[1]) * _vOUT) - sflx;
            */
            //===========================================================================================================
            //===========================================================================================================



            //===========================================================================================================
            //===========================================================================================================
            // Second variant using Roe-Type Scheme
            // Normal velocities
            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            double vINxN = 0.0, vOUTxN = 0.0;
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = _uIN[1 + d];
                vINxN += VelocityMeanIn[d] * inp.Normal[d];
                VelocityMeanOut[d] = _uOUT[1 + d];
                vOUTxN += VelocityMeanOut[d] * inp.Normal[d];
            }
            
            double uJump = _uIN[0] - _uOUT[0];

            flx = 0.5 * (Math.Min(vINxN, vOUTxN) - Math.Abs(Math.Min(vINxN, vOUTxN))) * -uJump * _vIN;
            flx += 0.5 * (Math.Max(vINxN, vOUTxN) + Math.Abs(Math.Max(vINxN, vOUTxN))) * -uJump * _vOUT;
            //===========================================================================================================
            //===========================================================================================================

            double ret = Scale * cap * flx;

            // DIAGNOSTIC (env-gated, zero cost unless set): dump every InnerEdgeForm call
            // touching a target cell -- traces the real production Roe-type flux's raw
            // per-node/per-probe values (including the Math.Min/Max/Abs kink-sensitive
            // vINxN/vOUTxN), on an ACTUAL XnsfeReproTest.cs production run.
            if (System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_EDGE") == "1") {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_CELL");
                int targetCell = targetStr != null ? int.Parse(targetStr) : -1;
                string targetEdgeStr = System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_EDGENO");
                int targetEdgeNo = targetEdgeStr != null ? int.Parse(targetEdgeStr) : -1;
                bool edgeMatch = targetEdgeNo < 0 || inp.iEdge == targetEdgeNo;
                bool speciesFilterOk = ValidSpecies == (System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_SPECIES") ?? ValidSpecies);
                if (edgeMatch && speciesFilterOk && (targetCell < 0 || inp.jCellIn == targetCell || inp.jCellOut == targetCell)) {
                    int c = System.Threading.Interlocked.Increment(ref _diagEdgeCount);
                    if (c <= 5000000) {
                        System.Console.Error.WriteLine(string.Format(
                            "DIAG_REALHEAT_EDGE[{0}]: species={1} edge={2} jIn={3} jOt={4} X=({5:0.######},{6:0.######},{7:0.######}) N=({8:0.########},{9:0.########},{10:0.########}) uIN=[{11}] uOT=[{12}] vIN={13:0.####} vOT={14:0.####} vINxN={15:0.########} vOUTxN={16:0.########} uJump={17:0.########} flx={18:0.##########} Scale={19:0.####} cap={20:0.####} ret={21:0.##########}",
                            c, ValidSpecies, inp.iEdge, inp.jCellIn, inp.jCellOut,
                            inp.X[0], inp.X[1], inp.X[2],
                            inp.Normal[0], inp.Normal[1], inp.Normal[2],
                            string.Join(",", _uIN.Select(v => v.ToString("0.####"))),
                            string.Join(",", _uOUT.Select(v => v.ToString("0.####"))),
                            _vIN, _vOUT, vINxN, vOUTxN, uJump, flx, Scale, cap, ret));
                    }
                }
            }

            return ret;
        }

        static int _diagEdgeCount = 0;
        static int _diagBndCount = 0;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            // ABLATION TOGGLE (coordinator-requested, this round): fully disables this form's
            // ConstantTemperature-wall contribution, to test on the REAL running solver whether
            // this newly-discovered contributor (never checked before this round) is implicated
            // in the periodic "Mismatch" defect. Env-gated, zero cost / identical to normal
            // behavior unless set.
            if (System.Environment.GetEnvironmentVariable("XNSFE_REPRO_DISABLE_BNDCONV") == "1")
                return 0.0;

            ThermalBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case ThermalBcType.ConstantTemperature: {
                    double flx = 0.0;



                    // Dirichlet value for temperature
                    double[] _uOUT = new double[m_D + 1];
                    _uOUT[0] = TempFunction[inp.EdgeTag](inp.X, inp.time);

                    // Outer values for Velocity equal inner values in this term
                    for (int j = 0; j < m_D; j++) {
                        _uOUT[j + 1] = _uIN[j + 1]; //velFunction[inp.EdgeTag, j](inp.X, inp.time);
                    }

                    //normal velocities
                    double[] VelocityMeanIn = new double[m_D];
                    double vINxN = 0.0;
                    for (int d = 0; d < m_D; d++) {
                        VelocityMeanIn[d] = _uIN[1 + d];
                        vINxN += VelocityMeanIn[d] * inp.Normal[d];
                    }

                    //===========================================================================================================
                    //===========================================================================================================
                    // First variant, using central flux for temperature
                    /*
                    double FlxNeg = _uIN[0] * (_uIN[1] * inp.Normal[0] + _uIN[2] * inp.Normal[1]);
                    if (m_D == 3) {
                        FlxNeg += _uIN[0] * _uIN[3] * inp.Normal[2];
                    }

                    // Term from partial integration back to strong form
                    double sflx = FlxNeg * _vIN;

                    // Flux with boundary value - Flux with inner value
                    flx = _uOUT[0] * (vINxN * _vIN ) - sflx;
                    */
                    //===========================================================================================================
                    //===========================================================================================================


                    //===========================================================================================================
                    //===========================================================================================================
                    // Second variant using Roe-Type Scheme
                    // if VxN < 0 (Inflow) enforce Dirichlet condition, if > 0 (outflow) we still enforce the boundary value?
                    double uJump = _uOUT[0] - _uIN[0];
                    flx = 0.5 * (vINxN - Math.Abs(vINxN)) * uJump * _vIN;
                    flx -= 0.5 * (vINxN + Math.Abs(vINxN)) * uJump * _vIN;
                    //===========================================================================================================
                    //===========================================================================================================

                    double retBnd = Scale * cap * flx;

                    // DIAGNOSTIC (env-gated, zero cost unless set): dump every BoundaryEdgeForm
                    // call touching a target cell -- coordinator-requested full M-row
                    // reconstruction found cell 8 has 3 wall (ConstantTemperature) boundary
                    // faces that were previously omitted from the interior-edges+volume-only
                    // reconstruction, and this form's velocity dependence (via vINxN) means
                    // they can contribute to Velocity-Temperature Jacobian columns too.
                    if (System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_BND") == "1") {
                        string targetStr = System.Environment.GetEnvironmentVariable("DIAG_REALHEAT_CELL");
                        int targetCell = targetStr != null ? int.Parse(targetStr) : -1;
                        if (targetCell < 0 || inp.jCellIn == targetCell) {
                            int c = System.Threading.Interlocked.Increment(ref _diagBndCount);
                            if (c <= 2000000) {
                                System.Console.Error.WriteLine(string.Format(
                                    "DIAG_REALHEAT_BND[{0}]: edge={1} jCellIn={2} EdgeTag={3} X=({4:0.######},{5:0.######},{6:0.######}) N=({7:0.########},{8:0.########},{9:0.########}) uIN=[{10}] uOUT_Dirichlet=[{11}] vIN={12:0.####} vINxN={13:0.########} uJump={14:0.########} flx={15:0.##########} ret={16:0.##########}",
                                    c, inp.iEdge, inp.jCellIn, inp.EdgeTag,
                                    inp.X[0], inp.X[1], inp.X[2],
                                    inp.Normal[0], inp.Normal[1], inp.Normal[2],
                                    string.Join(",", _uIN.Select(v => v.ToString("0.########"))),
                                    string.Join(",", _uOUT.Select(v => v.ToString("0.########"))),
                                    _vIN, vINxN, uJump, flx, retBnd));
                            }
                        }
                    }

                    return retBnd;
                }
                case ThermalBcType.ZeroGradient:
                case ThermalBcType.ConstantHeatFlux: { 
                    return 0.0;
                }
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }
        }
    }    

    public class HeatConvectionInSpeciesBulk_Upwind : LinearFlux, ISpeciesFilter {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        ThermalMultiphaseBoundaryCondMap m_bcmap;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelFunction;

        protected Func<double[], double, double>[] TempFunction;

        public HeatConvectionInSpeciesBulk_Upwind(int SpatDim, ThermalMultiphaseBoundaryCondMap _bcmap, string spcName, double _cap) {

            this.m_SpatialDimension = SpatDim;

            this.cap = _cap;
            //this.m_spcId = spcId;
            this.ValidSpecies = spcName;
            this.m_bcmap = _bcmap;

            this.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                this.VelFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            this.TempFunction = m_bcmap.bndFunction[VariableNames.Temperature + "#" + spcName];
        }

        double cap;

        //SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_SpatialDimension)); //, VariableNames.Velocity0MeanVector(m_SpatialDimension));
            }
        }


        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {

            for(int d = 0; d < m_SpatialDimension; d++)
                output[d] = cap * U[0] * inp.Parameters[d];

        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {

            double c = 0.0;
            for (int d = 0; d < m_SpatialDimension; d++)
                c += inp.Parameters_IN[d] * inp.Normal[d];

            ThermalBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case ThermalBcType.ConstantTemperature: {
                        return (c * cap * TempFunction[inp.EdgeTag](inp.X, inp.time));
                    }
                case ThermalBcType.ZeroGradient:
                case ThermalBcType.ConstantHeatFlux: {
                        return (c * cap * Uin[0]);
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
                return (c * cap * Uin[0]);
            else
                return (c * cap * Uout[0]);


        }
    }
}
