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


namespace BoSSS.Application.LsTest {

    /// <summary>
    /// Basic test case for advecting the level-set field
    /// in a constant velocity field
    /// </summary>
    public class LevelSetAdvectionTest : ILevelSetTest {

        protected double L = 1.0;
        protected double Radius = 0.4;

        protected double Ux = 1;
        protected double Uy = 0.5;
        protected double Uz = 0.7;

        protected double T = 0.4; // (Radius / Ux)

        protected double sign;

        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetAdvectionTest(int spatDim, int LevelSetDegree, bool reversed) {
            this.SpatialDimension = spatDim;
            this.LevelsetPolynomialDegree = LevelSetDegree;
            sign = (reversed) ? -1 : 1;
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
        public double ComputeTimestep(int Resolution, int LSdegree, int AMRlevel) {
            int gridCells1D = (9 * Resolution) * (AMRlevel + 1);
            double h = 1.0 * L / (double)gridCells1D;
            double dt = h / Ux;
            dt /= (double)(LSdegree * LSdegree);

            int timesteps = (int)Math.Ceiling(T / dt); // make sure that the singular point in time is exactly on one timestep

            return (T/(double)timesteps);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double getEndTime() {
            return (Radius / Ux);
        }
       
        /// <summary>
        /// creates a square in 2D and a cube in 3D
        /// </summary>
        /// <param name="Resolution"></param>
        /// <returns></returns>
        public virtual GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();
            
            //double L = 1.0;
            var xNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);
            var yNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);
            var zNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);

            GridCommons grd;
            switch (SpatialDimension) {
                case 2:
                grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                break;
                case 3:
                grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                break;
                default:
                throw new ArgumentOutOfRangeException();
            }

            grd.EdgeTagNames.Add(1, "wall");


            grd.DefineEdgeTags(delegate (double[] X) {
                return 1;
            });

            return grd;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public virtual IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("wall", new AppControl.BoundaryValueCollection()); 

            return config;
        }

        bool quadratic = false;

        /// <summary>
        /// Level-Set movement.
        /// </summary>
        public Func<double[], double, double>[] GetPhi() {
            return new Func<double[], double, double>[] { delegate (double[] X, double time) {

                double x0 = -sign * (Radius / 2.0);
                double y0 = -sign * (Radius / 4.0);
                double x = X[0] - x0, y = X[1] - y0;
                x -= Math.Sign((Radius / Ux) - time) * time * sign * Ux;
                y -= Math.Sign((Radius / Ux) - time) * time * sign * Uy;

                double dist;
                switch (SpatialDimension) {
                    case 2: { dist = Math.Sqrt(x * x + y * y); break; }
                    case 3: { double z0 = -sign * (Radius / 3.0);
                        double z = X[2]; z -= Math.Sign((Radius / Ux) - time) * time * sign * Uz;
                        dist = Math.Sqrt(x * x + y * y + z * z); break; }
                    default:
                    throw new ArgumentOutOfRangeException();
                }

                return quadratic ? (Radius * Radius - dist * dist) : Radius - dist;

            }};
        }

        /// <summary>
        /// velocity: constant velocity field with Ux in x-direction.
        /// </summary>
        public Func<double[], double, double>[][] GetU() {
            var Ret = new Func<double[], double, double>[this.NoOfLevelsets][];
            switch (SpatialDimension) {
                case 2:
                    Ret[0] = new Func<double[], double, double>[] { 
                        (X, t) => Math.Sign((Radius / Ux) - t) * sign * Ux,
                        (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uy
                    };
                    break;
                case 3:
                    Ret[0] = new Func<double[], double, double>[] {
                        (X, t) => Math.Sign((Radius / Ux) - t) * sign * Ux,
                        (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uy,
                        (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uz
                    };
                    break;
                default:
                    throw new ArgumentOutOfRangeException();
            }
            return Ret;
        }
       
        public int LevelsetPolynomialDegree {
            get;
            private set;
        }

        public double[,] AcceptableError {
            get {
                return new double[,] { { 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3 } };
            }
        }

        /// <summary>
        /// returns spatial dimension
        /// </summary>
        public int SpatialDimension {
            get;
            private set;
        }

        public int NoOfLevelsets => 1;
    }


    /// <summary>
    /// Basic test case for advecting the level-set field
    /// in a constant velocity field, using slip walls (does make a difference in computation of StokesExtension)
    /// </summary>
    class LevelSetAdvectionOnSlipWallTest : LevelSetAdvectionTest {

        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetAdvectionOnSlipWallTest(int spatDim, int LevelSetDegree, bool reversed) : base(spatDim, LevelSetDegree, reversed) {           
        }

        /// <summary>
        /// creates a square in 2D and a cube in 3D
        /// </summary>
        /// <param name="Resolution"></param>
        /// <returns></returns>
        public override GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();

            //double L = 1.0;
            var xNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);
            var yNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);
            var zNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);

            GridCommons grd;
            switch (SpatialDimension) {
                case 2:
                    grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                    break;
                case 3:
                    grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                    break;
                default:
                    throw new ArgumentOutOfRangeException();
            }

            grd.EdgeTagNames.Add(1, "NavierSlip_linear");


            grd.DefineEdgeTags(delegate (double[] X) {
                return 1;
            });

            return grd;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("NavierSlip_linear", new AppControl.BoundaryValueCollection());

            return config;
        }
    }


}
