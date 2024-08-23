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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Application.XNSE_Solver.Tests;

namespace BoSSS.Application.XNSFE_Solver.Tests {


    interface IXHeatTest : ITest {
        /// <summary>
        /// heat capacity of fluid A
        /// </summary>
        double c_A { get; }

        /// <summary>
        /// heat capacity of fluid A
        /// </summary>
        double c_B { get; }

        /// <summary> heat conductivity fluid A </summary>
        double k_A { get; }

        /// <summary> heat conductivity fluid B </summary>
        double k_B { get; }

        /// <summary> saturation temperature </summary>
        double T_sat { get; }

        /// <summary> latent heat of evaporation </summary>
        double h_vap { get; }

        /// <summary>
        /// Exact solution/Initial value for Temperature, for species <paramref name="species"/>.
        /// </summary>
        Func<double[], double, double> GetT(string species);

        /// <summary>
        /// Exact solution for total thermal energy.
        /// </summary>
        Func<double, double> GetE();

        bool CheckT { get; }
        bool CheckE { get; }
    }

    interface IXNSFETest : IXNSETest, IXHeatTest {

        /// <summary> 
        /// Some volumetric heat source.
        /// </summary>
        Func<double[], double> GetQ(string species);
        ///// <summary>
        ///// contact angle
        ///// </summary>
        //double theta_e { get; }
    }
}
