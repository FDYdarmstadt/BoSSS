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

using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.XNSEC{
 

    public interface IXNSECTest_Heat: IXNSECTest {

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


        /// <summary> diffusivity factor fluid A </summary>
        double rhoD_A { get; }

        /// <summary> diffusivity factor fluid B  </summary>
        double rhoD_B { get; }


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
    public interface IXNSECTest : IXNSETest {

        /// <summary>
        /// Exact solution/Initial value for Temperature, for species <paramref name="species"/>.
        /// </summary>
        Func<double[], double, double> GetTemperature(string species);

        /// <summary>
        /// Exact solution/Initial value for mass fractions, for species <paramref name="species"/>, component <paramref name="comp"/>.
        /// </summary>
        Func<double[], double, double> GetMassFractions(string species, int comp);

        /// <summary>
        /// Total number of chemical components involved in the solution
        /// </summary>
        int NumberOfChemicalComponents { get; }

        /// <summary>
        /// Activate chemical reaction related terms in the energy and species equations
        /// </summary>
        bool ChemicalReactionTermsActive { get; }

        /// <summary>
        /// Activate MassFraction equations
        /// </summary>
        bool EnableMassFractions { get; }

        /// <summary>
        /// Activate temperature equation
        /// </summary>
        bool EnableTemperature { get; }

        /// <summary>
        /// Directional vector of gravity
        /// </summary>
        double[] GravityDirection { get; }
    }


    public interface IXNSECTest_MixtureFraction : IXNSETest {

        /// <summary>
        /// Exact solution/Initial value for Temperature, for species <paramref name="species"/>.
        /// </summary>
        Func<double[], double, double> GetMixtureFraction(string species);


        /// <summary>
        /// Directional vector of gravity
        /// </summary>
        double[] GravityDirection { get; }
    }



    public interface IPrescribedMass : IXNSETest {
        /// <summary>
        /// only available if no heat equation is solved
        /// </summary>
        Func<double[], double, double> GetPrescribedMassflux_Evaluator();
    }
}