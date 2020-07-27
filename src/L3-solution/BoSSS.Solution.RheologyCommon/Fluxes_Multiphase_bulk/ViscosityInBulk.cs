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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;

namespace BoSSS.Solution.RheologyCommon {

    /// <summary>
    /// Viscosity part of constitutive equations in bulk for multiphase.
    /// </summary>
    public class ViscosityInBulk : ConstitutiveEqns_Viscosity, ISpeciesFilter {

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

        /// <summary>
        /// Initialize Viscosity
        /// </summary>
        public ViscosityInBulk(int _Component, IncompressibleMultiphaseBoundaryCondMap _BcMap, double _beta, double[] _Penalty1, string spcName, SpeciesId spcId) : base(_Component, _BcMap, _beta, _Penalty1) {
            this.validSpeciesId = spcName;
            this.m_bcMap = _BcMap;

            base.VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2];
            base.VelFunction.SetColumn(m_bcMap.bndFunction[VariableNames.VelocityX + "#" + spcName], 0);
            base.VelFunction.SetColumn(m_bcMap.bndFunction[VariableNames.VelocityY + "#" + spcName], 1);
        }

        public string validSpeciesId { get; private set; }
    }
}
