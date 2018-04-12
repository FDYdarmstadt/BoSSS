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

namespace BoSSS.Solution.XNSECommon.Operator.Viscosity {
    class ViscosityInBulk_GradUTerm_Localized : ViscosityInBulk_GradUTerm {

        public ViscosityInBulk_GradUTerm_Localized(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, int d, int D, double _muA, double _muB, ViscosityImplementation _ViscosityImplementation)
            : base(penalty, sw, bcMap, d, D, _muA, _muB, _ViscosityImplementation) { }



        public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);
            double muA = this.Viscosity(inp.Parameters_IN);
            double muB = this.Viscosity(inp.Parameters_OUT);


            switch (base.m_implMode) {
                case ViscosityImplementation.H: {
                        //only inner edges
                        for (int d = 0; d < inp.D; d++) {
                            //Acc += 0.5 * (muA * _Grad_uA[0, d] + muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                            //Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
                            Acc += (muA * _Grad_uA[0, d]) * (_vA) * inp.Normale[d];  // consistency term
                            //Acc += (muA * _Grad_vA[d]) * (_uA[0]) * inp.Normale[d];  // symmetry term
                        }
                        Acc *= base.m_alpha;

                        double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
                        //No Penalty Term for Localized Version -> no dependency on jump


                        return -Acc;

                    }

                case ViscosityImplementation.SWIP: {
                        throw new NotImplementedException();
                    }
                default: throw new NotImplementedException();
            }
        }
    }

    class ViscosityInBulk_GradUtranspTerm_Localized : ViscosityInBulk_GradUtranspTerm {

        public ViscosityInBulk_GradUtranspTerm_Localized(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, int d, int D, double _muA, double _muB, ViscosityImplementation _ViscosityImplementation) : base(penalty, sw, bcMap, d, D, _muA, _muB, _ViscosityImplementation) { }

            public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
                return 0;
            }
        }

        class ViscosityInBulk_divTerm_Localized : ViscosityInBulk_divTerm {

            public ViscosityInBulk_divTerm_Localized(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, int d, int D, double _muA, double _muB, ViscosityImplementation _ViscosityImplementation) : base(penalty, sw, bcMap, d, D, _muA, _muB, _ViscosityImplementation) { }

            public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
                return 0;
            }
        }
}