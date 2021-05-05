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

namespace BoSSS.Application.XNSE_Solver.Tests {


    /// <summary>
    /// Interface for tests (historical stuff). 
    /// </summary>
    public interface ITest {

        /// <summary>
        /// Gust mat the grid.
        /// </summary>
        int SpatialDimension {
            get;
        }

        /// <summary>
        /// Level set field in dependence of time.
        /// </summary>
        Func<double[], double, double> GetPhi();

        ///// <summary>
        ///// if known the exact solution will be set
        ///// </summary>
        //bool setExtSol { get; }

        /// <summary>
        /// Time step size.
        /// </summary>
        double dt { get; }

        /// <summary>
        /// Grid creation function.
        /// </summary>
        /// <param name="Resolution">1,2,3, etc</param>
        GridCommons CreateGrid(int Resolution);

        /// <summary>
        /// Boundary conditions and values.
        /// </summary>
        IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig();

        /// <summary>
        /// density of fluid A
        /// </summary>
        double rho_A { get; }

        /// <summary>
        /// density of fluid A
        /// </summary>
        double rho_B { get; }

        ///// <summary> interface speed in normal direction at time <paramref name="t"/>. </summary>
        //BoSSS.Foundation.ScalarFunction GetS(double time);


        ///// <summary> some external surface force (usually, only of use for manufactured solutions)</summary>
        //ScalarFunction GetSF(double time, int d);

        /// <summary>
        /// is the interface a material one or is it non-material?
        /// </summary>
        bool Material { get; }

        /// <summary>
        /// steady or in-stationary testcase?
        /// </summary>
        bool steady { get; }

        /// <summary>
        /// Stokes (false) or Navier-Stokes (true);
        /// </summary>
        bool IncludeConvection { get; }

        /// <summary>
        /// required polynomial degree for the level-set function
        /// </summary>
        int LevelsetPolynomialDegree { get; }

        /// <summary>
        /// the acceptable error of the solution, in comparison to the exact solution, for velocity and pressure.
        /// </summary>
        double[] AcceptableL2Error {
            get;
        }

        /// <summary>
        /// the acceptable Residual, for momentum and continuity equation.
        /// </summary>
        double[] AcceptableResidual {
            get;
        }
    }

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

    public interface IXNSETest : ITest {

        /// <summary> 
        /// Some external volume force, e.g. gravity.
        /// </summary>
        Func<double[], double> GetF(string species, int d);

        /// <summary>
        /// Exact solution/Initial value for velocity, for species <paramref name="species"/>, vector component <paramref name="d"/>.
        /// </summary>
        Func<double[], double, double> GetU(string species, int d);

        /// <summary>
        /// Exact solution/Initial value for pressure, for species <paramref name="species"/>.
        /// </summary>
        Func<double[], double, double> GetPress(string species);

        /// <summary> dynamic viscosity fluid A </summary>
        double mu_A { get; }

        /// <summary> dynamic viscosity fluid B </summary>
        double mu_B { get; }

        /// <summary>
        /// surface tension
        /// </summary>
        double Sigma { get; }

        /// <summary>
        /// <see cref="XNSE_Control.UseImmersedBoundary"/>
        /// </summary>
        bool TestImmersedBoundary { get; }

        /// <summary>
        /// Second level-set for the immersed boundary, only called if <see cref="TestImmersedBoundary"/> is true;
        /// </summary>
        Func<double[], double, double> GetPhi2();

        /// <summary>
        /// Velocity for the immersed boundary, only called if <see cref="TestImmersedBoundary"/> is true;
        /// </summary>
        Func<double[], double, double> GetPhi2U(int d);
    }

    interface IXNSFETest : IXNSETest, IXHeatTest {

    }



    public interface IXNSECTest : IXNSETest {
        
        /// <summary>
        /// Total number of chemical components involved in the solution
        /// </summary>
        int NumberOfChemicalComponents { get; }
        
        
        /// <summary>
        /// Activate chemical reaction related terms in the energy and species equations
        /// </summary>
        bool ChemicalReactionTermsActive { get; }

        /// <summary>
        /// Directional vector of gravity 
        /// </summary>
        double[] GravityDirection { get; }
    }
}
