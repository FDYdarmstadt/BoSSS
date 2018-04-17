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

        /// <summary>
        /// Exact solution/Initial value for velocity, for species <paramref name="species"/>, vector component <paramref name="d"/>.
        /// </summary>
        Func<double[],double,double> GetU(string species, int d);
       
        /// <summary>
        /// Exact solution/Initial value for pressure, for species <paramref name="species"/>.
        /// </summary>
        Func<double[], double, double> GetPress(string species);

        /// <summary>
        /// Time step size.
        /// </summary>
        double dt { get; }

        /// <summary>
        /// Grid creation function.
        /// </summary>
        GridCommons CreateGrid();

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

        /// <summary> dynamic viscosity fluid A </summary>
        double mu_A { get; }

        /// <summary> dynamic viscosity fluid B </summary>
        double mu_B { get; }

        ///// <summary> interface speed in normal direction at time <paramref name="t"/>. </summary>
        //BoSSS.Foundation.ScalarFunction GetS(double time);

        /// <summary> 
        /// Some external volume force, e.g. gravity.
        /// </summary>
        Func<double[],double> GetF(string species, int d);
        
        ///// <summary> some external surface force (usually, only of use for manufactured solutions)</summary>
        //ScalarFunction GetSF(double time, int d);

        /// <summary>
        /// surface tension
        /// </summary>
        double Sigma { get; }

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

    
    
}
