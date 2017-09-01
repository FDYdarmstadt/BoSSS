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
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Matrix representation of density.
    /// </summary>
    public class QuadratureMatrix_Rho : QuadratureMatrix {

        MaterialLaw EoS;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Basis"></param>s        
        /// <param name="GridDat"></param>
        /// <param name="EoS"></param>
        /// <param name="Scalar"></param>
        public QuadratureMatrix_Rho(Basis Basis, GridData GridDat, MaterialLaw EoS, params SinglePhaseField[] Scalar)
            : base(Basis, GridDat, Scalar) {
                this.EoS = EoS;
        }

        protected override double FieldFunc(params double[] FieldVal) {            
            return EoS.GetDensity(FieldVal);
        }
    }
}
