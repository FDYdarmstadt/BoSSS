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
    /// <summary>
    /// Convective term of the momentum equation. Compatible with newton solver
    /// </summary>
    public class LowMachCombustionConvectionInSpeciesBulk_LLF_Newton : LinearizedConvectionJacobi, ISpeciesFilter {
        /// <summary>
        /// Constructor for 
        /// </summary>
        public LowMachCombustionConvectionInSpeciesBulk_LLF_Newton(string spcName, int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, MaterialLaw EoS, int NumberOfComponents = -1) : base(SpatDim, _bcmap, _component, EoS, NumberOfComponents) {
            ValidSpecies = spcName;
        }

        /// <summary>
        /// Constructor for incompressible flow. Just for testing
        /// </summary>
        public LowMachCombustionConvectionInSpeciesBulk_LLF_Newton(string spcName, int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component) : base(SpatDim, _bcmap, _component) {
            ValidSpecies = spcName;
        }
        public string ValidSpecies {
            get;
            private set;
        }
    }
    /// <summary>
    /// Convective term of the momentum equation.
    /// </summary>
    public class LowMachCombustionConvectionInSpeciesBulk_LLF : LinearizedConvection, ISpeciesFilter {
        /// <summary>
        /// Constructor for 
        /// </summary>
        public LowMachCombustionConvectionInSpeciesBulk_LLF(string spcName, int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, MaterialLaw EoS, int NumberOfComponents = -1) : base(SpatDim, _bcmap, _component, EoS, NumberOfComponents) {
            ValidSpecies = spcName;
        }

        /// <summary>
        /// Constructor for incompressible flow. Just for testing
        /// </summary>
        public LowMachCombustionConvectionInSpeciesBulk_LLF(string spcName, int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component) : base(SpatDim, _bcmap, _component) {
            ValidSpecies = spcName;
        }
        public string ValidSpecies {
            get;
            private set;
        }
    }




    /// <summary>
    /// Scalar convection in the bulk phase for the low mach combustion solver.
    /// </summary>
    public class LowMachCombustion_ScalarConvectionInSpeciesBulk_LLF : LinearizedScalarConvection2Jacobi, ISpeciesFilter {
        public LowMachCombustion_ScalarConvectionInSpeciesBulk_LLF(string spcName, int SpatDim, int NumberOfReactants,
                                                                   IncompressibleBoundaryCondMap BcMap, MaterialLaw EoS,
                                                                   int idx) : base(SpatDim, NumberOfReactants, BcMap, EoS, idx) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }



}