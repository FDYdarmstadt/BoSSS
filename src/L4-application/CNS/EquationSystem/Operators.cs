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

namespace CNS.EquationSystem {
    
    /// <summary>
    /// The different operators of the compressible Navier-Stokes equations
    /// </summary>
    [Flags]
    public enum Operators {

        /// <summary>
        /// Empty equation
        /// </summary>
        None = 0,

        /// <summary>
        /// The convective operator
        /// </summary>
        Convection = 1,

        /// <summary>
        /// The diffusive operator
        /// </summary>
        Diffusion = 2,

        /// <summary>
        /// Optional gravity volume force
        /// </summary>
        Gravity = 4,
        
        /// <summary>
        /// Custom sources defined by the user
        /// </summary>
        CustomSource = 8,

        /// <summary>
        /// Sponge layer acting as damper near boundaries
        /// </summary>
        SpongeLayer = 16,

        /// <summary>
        /// Artificial viscosity to smear out shocks
        /// </summary>
        ArtificialViscosity = 32
    }
}
