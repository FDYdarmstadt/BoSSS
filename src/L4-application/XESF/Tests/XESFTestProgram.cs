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

namespace XESF.Tests
{
    [TestFixture]
    public static class XESFTestProgram
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

        #region SupersonicWedgeFlow using two LS on a Cartesian mesh
        [Test]
        public static void XDG_SWF_TwoLs()
        {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using (var p = new XESFMain())
            {
                var C = XESFHardCodedControl.XDGWedgeFlow_TwoLs_Base(
                    optiLSDegree: 1,
                    lsDegree: 1,
                    shocksetup: ApplicationWithIDT.GetLevelSet.FromParams,
                    optiLevelSetType: OptiLevelSetType.GlobalLevelSet,
                    initialValue: ApplicationWithIDT.GetInitialValue.FromFunctionPerSpecies,
                    MaxIterations: 200,
                    dgDegree: 0,
                    numOfCellsX: 15,
                    numOfCellsY: 10,
                    initialAngle_shockLS: 32,
                    PlotInterval: -1,
                    interfaceFluxLS2: ConvectiveInterfaceFluxes.GodunovInterface,
                    bulkFlux: ConvectiveBulkFluxes.OptimizedHLLC,
                    FluxVersion: Fluxes.FluxVersion.Optimized,
                    agg: 0.4,
                    globalization: ApplicationWithIDT.GlobalizationStrategy.LineSearch
                    );
                p.Init(C);
                p.RunSolverMode();
                var tol1 = 1e-02;
                var tol2 = 2.5;
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol2 && p.ResidualVector.MPI_L2Norm() < tol1), $"the L2 Error is greater than {tol1} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");
            }
        }
        #endregion
        #region SupersonicWedgeFlow using one LS on a Cartesian rotated mesh
        [Test]
        public static void XDG_SWF_OneLs_Cart()
        {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using (var p = new XESFMain())
            {
                var C = XESFHardCodedControl.XDGWedgeFlow_OneLs_Rotation(
                    MaxIterations: 100,
                    dgDegree: 0,
                    numOfCellsX: 10,
                    numOfCellsY: 15,
                    lsDegree: 1
                    );

                p.Init(C);
                p.RunSolverMode();
                var tol = 1e-07;
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol && p.ResidualVector.MPI_L2Norm() < tol), $"the L2 Error is greater than {tol} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");

            }
        }
        #endregion

        /// <summary>
        /// Helper function to store interpolating points of an SplineLevelSet
        /// </summary>
        /// <param name="p"></param>
        /// <param name="filename"></param>
        public static void SaveIsoContourToTextFile(XESFMain p, string filename)
        {
            if (p.LevelSetOpti is SplineOptiLevelSet spliny)
            {
                spliny.GetSpline();
                if (spliny.Spline is CubicSpline cSpliny)
                {
                    var yMax = ((GridData)p.Grid.iGridData).Vertices.Coordinates.ExtractSubArrayShallow(-1, 1).To1DArray().Max();
                    var yMin = ((GridData)p.Grid.iGridData).Vertices.Coordinates.ExtractSubArrayShallow(-1, 1).To1DArray().Min();
                    var yPoints = GenericBlas.Linspace(yMin, yMax, 200);
                    var xPoints = new double[100];
                    var cPoints = new double[100];
                    var allPoints = MultidimensionalArray.Create(100, 3);
                    for (int i = 0; i < 100; i++)
                    {
                        allPoints[i, 0] = yPoints[i];

                        xPoints[i] = cSpliny.Interpolate(yPoints[i]);
                        allPoints[i, 1] = xPoints[i];

                        cPoints[i] = cSpliny.Differentiate(yPoints[i]);
                        allPoints[i, 2] = cPoints[i];
                    }
                    allPoints.SaveToTextFile(filename);

                }
                else
                {
                    throw new NotSupportedException("not supported - but one could implement a Newton Root finding along lines of type x -> (x,y_i) for that");
                }
            }
        }

    }
}
