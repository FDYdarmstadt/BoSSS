using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using XESF.Fluxes;
using ilPSP.LinSolvers;
using ApplicationWithIDT;
using ApplicationWithIDT.OptiLevelSets;
using MathNet.Numerics.Interpolation;
using System.Linq;
using NUnit.Framework;
using System.Diagnostics.Metrics;
using XESTSF;

namespace XESTSF.Tests
{
    [TestFixture]
    public static class XESTSFTestProgram
    {
        #region NUnit stuff
        [OneTimeSetUp]
        static public void Init()
        {
            BoSSS.Solution.Application.InitMPI();
        }

        [OneTimeTearDown]
        static public void Cleanup()
        {
        }
        #endregion
        [Test]
        #region Stationary Shock Wave
        public static void XDG_SSW()
        {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using (var p = new XESTSFMain())
            {
                var C = XESTSFHardCodedControl.StationaryShockWave(
                    MaxIterations: 100,
                    dgDegree: 0,
                    numOfCellsX: 10,
                    numOfCellsY: 10,
                    lsDegree: 1
                    );

                p.Init(C);
                p.RunSolverMode();
                var tol = 1e-07;
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol && p.ResidualVector.MPI_L2Norm() < tol), $"the L2 Error is greater than {tol} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");

            }
        }

        [Test]
        public static void XDG_Shock_Acoustic_Interaction()
        {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using (var p = new XESTSFMain())
            {
                var C = XESTSFHardCodedControl.AcousticWave1D(
                    MachL:1.5,
                    p_amp_neg:0.0,
                    p_amp_pos:0.0001,
                    waveLength:1.0,
                    wavePosition:0.0,
                    shockPosition:1.5,
                    MaxIterations: 100,
                    dgDegree: 2,
                    numOfCellsX: 20,
                    numOfCellsT: 20,
                    lsDegree: 3                    
                    );

                p.Init(C);
                p.RunSolverMode();

                //check if converged
                var tol = 1e-04;
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol && p.ResidualVector.MPI_L2Norm() < tol), $"the L2 Error is greater than {tol} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");
                //check if interaction is correct
                p.DerivedVariableToXDGFieldMap.TryGetValue(Variables.XESTSFVariables.PertubationPressure, out XDGField p_per);

                //CellMask AllCells = CellMask.GetFullMask(p.GridData);
                //double[] mins = new double[AllCells.NoOfItemsLocally];
                //double[] maxs = new double[AllCells.NoOfItemsLocally];
                //p_per.GetCellwiseExtremalValues(mins,maxs);
                //Assert.IsTrue(Math.Abs(maxs.Max()-p_per_max_theo)<tol,"value is not right");


            }
        }
        #endregion

    }
}
