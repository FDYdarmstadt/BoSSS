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
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.XNSECommon.Operator.Convection {

    public class ConvectionInBulk_LLF : LinearizedConvection, IEquationComponentSpeciesNotification {

        public ConvectionInBulk_LLF(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, int _component, double _rhoA, double _rhoB, double _LFFA, double _LFFB, LevelSetTracker _lsTrk) :
            base(SpatDim, _bcmap, _component, false) 
        {
            //
            rhoA = _rhoA;
            rhoB = _rhoB;
            //varMode = _varMode;
            this.lsTrk = _lsTrk;
            this.LFFA = _LFFA;
            this.LFFB = _LFFB;
            this.m_bcmap = _bcmap;
            base.velFunction = null;
        }

        IncompressibleMultiphaseBoundaryCondMap m_bcmap;
        LevelSetTracker lsTrk;

        double LFFA;
        double LFFB;


        public void SetParameter(string speciesName, SpeciesId SpcId) {
            switch (speciesName) {
                case "A": this.rho = this.rhoA; base.LaxFriedrichsSchemeSwitch = LFFA; this.SetBndfunc("A"); break;
                case "B": this.rho = this.rhoB; base.LaxFriedrichsSchemeSwitch = LFFB; this.SetBndfunc("B"); break;
                default: throw new ArgumentException("Unknown species.");
            }
            SubGrdMask = lsTrk.Regions.GetSpeciesSubGrid(SpcId).VolumeMask.GetBitMaskWithExternal();
        }

        void SetBndfunc(string S) {
            int SpatDim = base.m_SpatialDimension;
            int dir = base.m_component;
            base.velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for(int d = 0; d < SpatDim; d++)
                base.velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + S], d);
        }

        protected System.Collections.BitArray SubGrdMask;


        double rhoA;
        double rhoB;
        protected double rho;
        
        internal double IEF(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return this.InnerEdgeFlux(ref inp, Uin, Uout);
        }

        protected bool basecall = false;

        protected override double InnerEdgeFlux(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            if (basecall) {
                return base.InnerEdgeFlux(ref inp, Uin, Uout);
            } else {

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
}
