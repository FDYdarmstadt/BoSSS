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

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            Scale = cs.HomotopyValue;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivVol, DivergenceDerivEdg };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for (int d = 0; d < m_D; d++) {
                // starke Form
                acc += U[1 + d] * GradU[0, d] * V;
            }
            return Scale * cap * acc;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
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

            return Scale * cap * flx;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
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

                    return Scale * cap * flx;
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
