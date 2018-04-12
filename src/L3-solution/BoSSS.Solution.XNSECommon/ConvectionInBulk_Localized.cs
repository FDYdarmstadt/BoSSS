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
using ilPSP;

namespace BoSSS.Solution.XNSECommon.Operator.Convection {
    class ConvectionInBulk_Localized : ConvectionInBulk_LLF {

        public ConvectionInBulk_Localized(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, int _component, double _rhoA, double _rhoB, double _LFFA, double _LFFB, LevelSetTracker _lsTrk) : base(SpatDim, _bcmap, _component, _rhoA, _rhoB, _LFFA, _LFFB, _lsTrk) { }

        protected override double InnerEdgeFlux(ref BoSSS.Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            if(basecall) {
                return base.InnerEdgeFlux(ref inp, Uin, Uout);
            } else {

                double UinBkUp = Uin[0];
                double UoutBkUp = Uout[0];
                double[] InParamsBkup = inp.Parameters_IN;
                double[] OutParamsBkup = inp.Parameters_OUT;
                int m_SpatialDimension = Uin.Length;

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

                double rhoIn = 1.0;


                // 2 * {u_i * u_j} * n_j,
                // same as standard flux without outer values
                flx += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
                if(m_SpatialDimension == 3) {
                    flx += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normale[2];
                }

                // Calculate dissipative part
                // ==========================

                IList<double> _VelocityMeanIn = new List<double>();
                IList<double> _VelocityMeanOut = new List<double>();
                for(int d = 0; d < m_SpatialDimension; d++) {
                    _VelocityMeanIn.Add(inp.Parameters_IN[m_SpatialDimension + d]);
                    _VelocityMeanOut.Add(inp.Parameters_OUT[m_SpatialDimension + d]);
                }
                double[] VelocityMeanIn = _VelocityMeanIn.ToArray();
                double[] VelocityMeanOut = _VelocityMeanOut.ToArray();

                flx *= base.rho;

                // cleanup mess and return
                // -----------------------

                Uout[0] = UoutBkUp;
                Uin[0] = UinBkUp;
                inp.Parameters_IN = InParamsBkup;
                inp.Parameters_OUT = OutParamsBkup;

                return flx;
            }
        }

    }
}
