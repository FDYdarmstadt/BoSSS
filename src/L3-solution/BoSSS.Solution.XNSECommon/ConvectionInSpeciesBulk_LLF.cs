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
using System.Diagnostics;

using ilPSP;
using ilPSP.Utils;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;

using BoSSS.Solution.NSECommon;


namespace BoSSS.Solution.XNSECommon.Operator.Convection {

    public class ConvectionInSpeciesBulk_LLF : LinearizedConvection, ISpeciesFilter, IEquationComponentCoefficient {

        public ConvectionInSpeciesBulk_LLF(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId, int _component, 
            double _rho, double _LFF, LevelSetTracker _lsTrk) :
            base(SpatDim, _bcmap, _component, false) {
            //
            rho = _rho;
            m_spcId = spcId;

            //varMode = _varMode;
            this.lsTrk = _lsTrk;
            this.LFF = _LFF;

            this.m_bcmap = _bcmap;

            int dir = base.m_component;
            base.velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for(int d = 0; d < SpatDim; d++)
                base.velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            //SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(spcId).VolumeMask.GetBitMaskWithExternal();
        }

        IncompressibleMultiphaseBoundaryCondMap m_bcmap;
        LevelSetTracker lsTrk;

        double LFF;


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        protected System.Collections.BitArray SubGrdMask;


        protected double rho;

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

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(m_spcId).VolumeMask.GetBitMaskWithExternal();
        }
    }
}
