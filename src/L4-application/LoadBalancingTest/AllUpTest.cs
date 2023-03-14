﻿using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.LoadBalancing;
using MPI.Wrappers;
using NUnit.Framework;
using System;

namespace BoSSS.Application.LoadBalancingTest {

    /// <summary>
    /// Complete Test for the load balancing.
    /// </summary>
    [TestFixture]
    static public class AllUpTest {

        /*
        /// <summary>
        /// MPI init
        /// </summary>
        [OneTimeSetUp]
        public static void SetUp() {
            BoSSS.Solution.Application.InitMPI();
        }
        /*
        /// <summary>
        /// Da Test!
        /// </summary>
        [Test]
        static public void NoDynamicBalanceTest([Values(1, 2)] int DGdegree) {
            LoadBalancingTestMain p = null;

            BoSSS.Solution.Application._Main(
                new string[0],
                true,
                null,
                delegate() {
                    p = new LoadBalancingTestMain();
                    p.DynamicBalance = false;
                    p.DEGREE = DGdegree;
                    p.cellCostEstimatorFactory = CellCostEstimatorLibrary.OperatorAssemblyAndCutCellQuadrules;
                    return p;
                });
        }
        */

        /// <summary>
        /// Da Test!
        /// </summary>
        [Test]
        static public void RuntimeCostDynamicBalanceTest(
            [Values(1, 2)] int DGdegree) {
            LoadBalancingTestMain p = null;
            
            BoSSS.Solution.Application<AppControlSolver>._Main(
                new string[0],
                true,
                delegate () {
                    p = new LoadBalancingTestMain();
                    p.DynamicBalance = true;
                    p.DEGREE = DGdegree;
                    p.cellCostEstimatorFactory = () => CellCostEstimatorLibrary.OperatorAssemblyAndCutCellQuadrules;
                    return p;
                });
        }

        /*
        /// <summary>
        /// Da Test!
        /// </summary>
        [Test]
        static public void StaticCostDynamicBalanceTest(
            [Values(1, 2)] int DGdegree) {
            LoadBalancingTestMain p = null;

            BoSSS.Solution.Application._Main(
                new string[0],
                true,
                null,
                delegate () {
                    p = new LoadBalancingTestMain();
                    p.DynamicBalance = true;
                    p.DEGREE = DGdegree;
                    p.cellCostEstimatorFactory = delegate(IApplication<AppControl> app, int performanceClassCount) {
                        if (performanceClassCount != 2) {
                            throw new Exception();
                        }

                        int[] performanceClassToCostMap = new int[] {
                            1, // normal cells
                            10 // assume cut cells cost 10 times more
                        };
                        return new StaticCellCostEstimator(performanceClassToCostMap);
                    };
                    return p;
                });
        }
        */
       
    }
}
