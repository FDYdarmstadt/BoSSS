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
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using System.Globalization;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;


namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// Basic test class for the level-set unit tests
    /// </summary>
    abstract class LevelSetBaseTest : IXNSElsTest {

        public bool TestImmersedBoundary => false;

        /// <summary>
        /// nix
        /// </summary>
        public Func<double[], double, double> GetPhi2() {
            throw new NotImplementedException(); // will never be called, as long as 'TestImmersedBoundary' == false;
        }

        public Func<double[], double, double> GetPhi2U(int d) {
            throw new NotImplementedException();
        }


        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetBaseTest(int spatDim, int LevelSetDegree) {
            this.SpatialDimension = spatDim;
            this.LevelsetPolynomialDegree = LevelSetDegree;
        }


        public double dt {
            get {
                return -1.0;    // will be set in LevelSetTest() according to level set cfl 
            }
        }

        /// <summary>
        /// computes the timestep size according to the level-set CFL condition
        /// </summary>
        /// <param name="Resolution"></param>
        /// <param name="LSdegree"></param>
        /// <returns></returns>
        public abstract double ComputeTimestep(int Resolution, int LSdegree);

       
        /// <summary>
        /// creates grid according to set resolution
        /// </summary>
        /// <param name="Resolution"></param>
        /// <returns></returns>
        public abstract GridCommons CreateGrid(int Resolution);

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public abstract IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig();

        /// <summary>
        /// Level-Set function
        /// </summary>
        public abstract Func<double[], double, double> GetPhi();

        /// <summary>
        /// velocity field
        /// </summary>
        public abstract Func<double[], double, double> GetU(string species, int d);

        /// <summary>
        /// no pressure field 
        /// </summary>
        public Func<double[], double, double> GetPress(string species) {
            switch (SpatialDimension) {
                case 2: {
                    switch (species) {
                        case "A": return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                        case "B": return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                        default: throw new ArgumentException();
                    }
                }
                case 3: {
                    switch (species) {
                        case "A": return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                        case "B": return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                        default: throw new ArgumentException();
                    }
                }
                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// No bulk force.
        /// </summary>
        public Func<double[], double> GetF(string species, int d) {
            return (X => 0.0);
        }


        /// <summary>
        /// no surface tension effects during level set motion 
        /// </summary>
        public double Sigma {
            get {
                return 0.0;
            }
        }

        /// <summary>
        /// considering level set movement in single fluid (A=B)
        /// </summary>
        public double rho_A {
            get { return 1.0; }
        }

        /// <summary>
        ///  considering level set movement in single fluid (A=B)
        /// </summary>
        public double rho_B {
            get { return 1.0; }
        }

        /// <summary>
        ///  considering level set movement in single fluid (A=B)
        /// </summary>
        public double mu_A {
            get { return 1.0; }
        }

        /// <summary>
        ///  considering level set movement in single fluid (A=B)
        /// </summary>
        public double mu_B {
            get { return 1.0; }
        }


        public bool Material {
            get { return true; }
        }

        public bool steady {
            get { return false; }
        }

        public bool IncludeConvection {
            get { return false; }
        }

        public int LevelsetPolynomialDegree {
            get;
            private set;
        }

        public double[] AcceptableL2Error {
            get {
                  return new double[] { };
            }
        }

        public double[] AcceptableResidual {
            get {
                return  new double[] { };
            }
        }


        /// <summary>
        /// returns spatial dimension
        /// </summary>
        public int SpatialDimension {
            get;
            private set;
        }

    }

}
