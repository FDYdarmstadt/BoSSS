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
    /// Matrix representation of some variables for low Mach number flows.
    /// </summary>
    public class VariableMatrices {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="GridDat"></param>
        /// <param name="VelocityBasis"></param>        
        /// <param name="PressureBasis"></param>        
        /// <param name="EoS"></param>
        /// <param name="Scalar"></param>        
        public VariableMatrices(GridData GridDat, Basis VelocityBasis, Basis PressureBasis, MaterialLaw EoS, params SinglePhaseField[] Scalar) {

            // Construct matrices            
            m_Rho = new QuadratureMatrix_Rho(VelocityBasis, GridDat, EoS, Scalar);            

            // Initialize matrices
            m_Rho.Update();            
        }

        QuadratureMatrix m_Rho;

        /// <summary>
        /// Quadrature matrix for density.        
        /// </summary>
        public QuadratureMatrix Rho {
            get {
                return m_Rho;
            }
        }
    }
}
