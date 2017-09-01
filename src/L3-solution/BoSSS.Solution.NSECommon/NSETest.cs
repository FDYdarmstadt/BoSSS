using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Test for the Navier-Stokes Operator
    /// </summary>
    [TestFixture]
    class NSETest {

        /// <summary>
        /// template for different tests (exact solutions of the NSE)
        /// </summary>
        interface ITest {

            /// <summary>
            /// velocity components (exact solutions)
            /// </summary>
            /// <param name="d">spatial direction</param>
            ScalarFunction GetU(int d);
            
            /// <summary>
            /// creates an appropriate grid for the test case. 
            /// </summary>
            void CreateGrid(Context ctx, out GridCommons grd, out IncompressibleBoundaryCondMap bcMap);

            /// <summary>
            /// pressure (exact solution).
            /// </summary>
            ScalarFunction GetPress(double time, string species);
            
            /// <summary> Reynolds number. </summary>
            double mu_B { get; }
        }


        /// <summary>
        /// stream from the origin to infinity
        /// </summary>
        class ExpandingFlow : ITest {

            const double c = 1.0;

            /// <summary>
            /// <latex mode="display">\vec{u}(\vec{x}) =  \vec{x} \frac{c}{ | \vec{x} |^2 }</latex>
            /// </summary>
            public ScalarFunction GetU(int d) {
                return ((_2D)((x, y) =>
                        (((new double[] { x, y })[d] / (x * x + y * y)) // a very inefficient way to select the d-th component!
                        * c))).Vectorize();
            }

            /// <summary>
            /// Cartesian grid, nothing fancy.
            /// </summary>
            public void CreateGrid(Context _Context, out GridCommons grd, out IncompressibleBoundaryCondMap bcMap) {
                grd = new Cartesian2DGrid(_Context.CommMaster, Grid1D.Linspace(-2, 2, 10), Grid1D.Linspace(-2, 2, 10));
                grd.EdgeTagsNames.Add(1, "Pressure_Outlet");
                grd.DefineEdgeTags(x => 1);

                bcMap = null;
            }

            public ScalarFunction GetPress(double time, string species) {
                throw new NotImplementedException();
            }

            public double mu_B {
                get { throw new NotImplementedException(); }
            }
        }




        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public void Init() {
            bool dummy;
            ilPSP.Enviroment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
        }


        /// <summary>
        /// MPI teardown 
        /// </summary>
        [TestFixtureTearDown]
        public void Cleanup() {
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }


    }
}
