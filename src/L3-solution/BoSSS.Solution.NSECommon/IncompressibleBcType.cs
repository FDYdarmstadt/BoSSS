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

namespace BoSSS.Solution.NSECommon {
    
    
    /// <summary>
    /// types of boundary conditions for incompressible Navier-Stokes solvers;
    /// </summary>
    public enum IncompressibleBcType {

        /// <summary>
        /// fixed or moving wall with no-slip condition
        /// </summary>
        Wall = 0,

        /// <summary>
        /// inlet with predefined velocity: Dirichlet b.c. for transport and viscous operator
        /// </summary>
        Velocity_Inlet = 1,

        /// <summary>
        /// \[
        ///     \left( \psi \boldsymbol{I} - \frac{1}{\textrm{Re}} \nabla \underline{u} \right) \cdot \underline{n}_{\partial \Omega} = 0
        /// \]
        /// </summary>
        Pressure_Outlet = 2,

        /// <summary>
        /// Outflow with Neumann pressure boundary condition
        /// Within the SIMPLE algorithm this boundary condition implements
        /// \[
        ///     \left( \psi \boldsymbol{I} - \frac{1}{\textrm{Re}} \nabla \underline{u} \right) \cdot \underline{n}_{\partial \Omega} = 0
        /// \]
        /// </summary>
        SIMPLE_Outflow = 3,

        /// <summary>
        /// Dirichlet boundary condition for pressure.
        /// (no restriction for velocity / velocity gradient).
        /// </summary>
        Pressure_Dirichlet = 4,

        /// <summary>
        /// Wall with no-slip condition for velocity
        /// and Neumann boundary condition for scalar variable
        /// (i.e. Level-Set for multiphase and temperature for Low-Mach).
        /// </summary>
        NoSlipNeumann = 5,

        /// <summary>
        /// Wall with free slip condition:
        /// \[
        /// \underline{u} \cdot \underline{n} = 0 \text{ and } \underline{t} \cdot (\operatorname{grad}u + (\operatorname{grad}u)^T ) \cdot \underline{n} = 0 
        /// \]
        /// </summary>
        FreeSlip = 6,

        /// <summary>
        /// symmetry boundary condition with contact angle (90° contact angle) and free slip
        /// </summary>
        SlipSymmetry = 7,

        /// <summary>
        /// Generalized Navier Boundary condition with linear effective forces at wall and contact line
        /// \[
        /// \underline{u} \cdot \underline{n}_S = 0 \text{ and } 
        /// \boldsymbol{P}_S (\operatorname{grad}u + (\operatorname{grad}u)^T ) \underline{n}_S = -\beta_S \boldsymbol{P}_S \underline{u} \text{ on } \partial \Omega_{S}
        /// \sigma \boldsymbol{P}_S \tau_L = -\beta_L (\underline{u} \cdot \underline{n}_L) + \sigma \cos(\Theta_e) \underline{n}_L \text{ on } L
        /// \]
        /// </summary>
        NavierSlip_Linear = 8,

        /// <summary>
        /// Pressure outlet with Dirichlet values for the scalars.
        /// </summary>
        ScalarDirichlet_PressureOutlet = 9,

        /// <summary>
        /// outflow boundary condition according to Dong et. al.
        /// \[
        ///     \left( - \psi \boldsymbol{I} + \frac{1}{\textrm{Re}} \nabla \underline{u} - \frac{1}{2} \left\lvert{\underline{u}}\right\rvert^2 S_{out} (\underline{n} \cdot \underline{u}) \right) \cdot \underline{n}_{\partial \Omega} = 0
        /// \]
        /// where 
        /// \[
        ///     S_{out}(-) = \frac{1}{2}(1 - tanh(\frac{(-)}{U_0 \delta}))
        /// \]
        /// with U_0 is the characteristic velocity scale, and δ > 0 is a chosen non-dimensional positive constant that is sufficiently small
        /// </summary>
        Dong_OutFlow = 10

    }
}
