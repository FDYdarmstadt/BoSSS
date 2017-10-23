/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled or executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using BoSSS.Solution;
using NUnit.Framework;
using System.Globalization;
using System.Threading;

namespace ALTSTests {
    /// <summary>
    /// NUnit test class for ALTS
    /// </summary>
    [TestFixture]
    class NUnitTests : Program {

        public static void Main(string[] args) {
            Application._Main(
                args,
                true,
                "",
                () => new NUnitTests());
        }

        [TestFixtureSetUp]
        static public void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out dummy);
            Thread.CurrentThread.CurrentCulture = CultureInfo.CurrentCulture;
        }

        [TestFixtureTearDown]
        static public void Cleanup() {
            //Console.Out.Dispose();
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        public static void ALTSDynClustering(int order, int subGrids, double maxEnergyNorm) {
            NUnitTests test = null;

            Application._Main(
                new string[0],
                true,
                "",
                delegate () {
                    test = new NUnitTests() {
                        ABOrder = order,
                        numOfSubgrids = subGrids,
                    };
                    test.m_GridPartitioningType = BoSSS.Foundation.Grid.GridPartType.none;
                    return test;
                });

            double energyNorm = test.energyNorm;

            Assert.IsTrue(energyNorm < maxEnergyNorm + 1e-15);
        }

        // Call tests
        [Test]
        // Here, A-LTS gives the same result as LTS because AB order 1 equals Explicit Euler.
        // In this case, restarting a LTS simulation with another clustering is possible
        // because a history is not needed. 
        public static void ALTSDynClust_order1_subgrids3() {
            //3 cells
            //ALTSDynClustering(order: 1, subGrids: 3, maxEnergyNorm: 7.772253056189100E-01);

            //4 cells
            ALTSDynClustering(order: 1, subGrids: 4, maxEnergyNorm: 7.905061733461980E-01);
        }

        [Test]
        public static void ALTSDynClust_order2_subgrids3() {
            //3 cells
            //ALTSDynClustering(order: 2, subGrids: 3, maxEnergyNorm: 7.772253058420590E-01);

            //4 cells
            ALTSDynClustering(order: 2, subGrids: 4, maxEnergyNorm: 7.905061732830720E-01);
        }

        [Test]
        public static void ALTSDynClust_order3_subgrids3() {
            //3 cells
            //ALTSDynClustering(order: 3, subGrids: 3, maxEnergyNorm: 7.772253058420650E-01);

            //4 cells
            ALTSDynClustering(order: 3, subGrids: 4, maxEnergyNorm: 7.905061732830850E-01);
        }
    }
}
