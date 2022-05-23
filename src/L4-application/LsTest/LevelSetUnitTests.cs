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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using ilPSP;
using BoSSS.Solution.AdvancedSolvers.Testing;
using ilPSP.Connectors.Matlab;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation;
using BoSSS.Solution.Tecplot;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Application.LsTest {


    public interface ILevelSetTest {


        /// <summary>
        /// Gust mat the grid.
        /// </summary>
        int SpatialDimension {
            get;
        }

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
        /// Level set field(s) in dependence of time, typically given as signed distance function
        /// </summary>
        Func<double[], double, double>[] GetPhi();

        /// <summary>
        /// Advection velocity for each level set
        /// 1st index: level set
        /// 2nd index: spatial direction
        /// </summary>
        Func<double[], double, double>[][] GetU();

        /// <summary>
        /// required polynomial degree for the level-set function, for now all level sets have the same degree
        /// </summary>
        int LevelsetPolynomialDegree { get; }

        /// <summary>
        /// number of level sets
        /// </summary>
        int NoOfLevelsets { get; }

        /// <summary>
        /// the acceptable error of L2-CG, L2-DG, volume, surface and contour.
        /// </summary>
        double[,] AcceptableError {
            get;
        }

        /// <summary>
        /// compute a suitable timestep for various combinations of grid resolutions and level-set degrees
        /// </summary>
        double ComputeTimestep(int gridResolution, int lsDegree, int AMRlevel);

        /// <summary>
        /// 
        /// </summary>
        double getEndTime();

    }

    /// <summary>
    /// A collection of all-up NUnit tests for the XNSE solver regarding various level set handling
    /// </summary>
    [TestFixture]
    static public partial class LevelSetUnitTests {

        /// <summary>
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        [Test]
        public static void LevelSetAdvectionTest2D_fwd(
           [Values(2, 3, 4)] int LSdegree,
           [Values(0, 1, 2)] int AMRlevel,
           [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
           [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {
            LevelSetAdvectionTest2D(LSdegree, AMRlevel, levelSetEvolution, levelSetHandling, false);

        }


        /// <summary>
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        [Test]
        public static void LevelSetAdvectionTest2D_reverse(
           [Values(2, 3, 4)] int LSdegree,
           [Values(0, 1, 2)] int AMRlevel,
           [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
           [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {
            LevelSetAdvectionTest2D(LSdegree, AMRlevel, levelSetEvolution, levelSetHandling, true);
        }



        /// <summary>
        /// <see cref="BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        public static void LevelSetAdvectionTest2D(

            [Values(2, 3, 4)] int LSdegree,
            [Values(0, 1, 2)] int AMRlevel,
            [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
            [Values(false, true)] bool reversed) {

            int gridResolution;
            switch (LSdegree) {
                case 2: gridResolution = 3; break;
                case 3: gridResolution = 2; break;
                //case 4: gridResolution = 1; break;
                default:
                    gridResolution = 1; break;
            }

            if (LSdegree == 4 && AMRlevel > 0 && levelSetEvolution == LevelSetEvolution.StokesExtension)
                return;

            var Tst = new LevelSetAdvectionTest(2, LSdegree, reversed);
            var C = LSTstObj2CtrlObj(Tst, 40, levelSetEvolution, levelSetHandling, gridResolution, AMRlevel);

            //string IO = $"LSAdvectionTest2D-deg{LSdegree}-amrLvl{AMRlevel}-lsEvo{levelSetEvolution}-rev{reversed}-grdRes{gridResolution}";
            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 3;

            LevelSetTest(Tst, C);

            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!

        }

        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetAdvectionTest3D(
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(0, 1, 2)] int AMRlevel,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling, 
        //    [Values(false, true)] bool reversed)
        //    {

        //    var Tst = new LevelSetAdvectionTest(3, LSdegree, reversed);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 20, levelSetEvolution, levelSetHandling, 1);
        //    C.SkipSolveAndEvaluateResidual = true;

        //    LevelSetTest(Tst, C);

        //}


        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetAdvectionOnWallTest2D(
        //    [Values(Math.PI/4, Math.PI/3, Math.PI/6)] double contactAngle,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(0, 1, 2)] int AMRlevel,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {

        //    var Tst = new LevelSetAdvectionOnWallTest(2, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling, 1, AMRlevel);
        //    C.SkipSolveAndEvaluateResidual = true;

        //    LevelSetTest(Tst, C);

        //}

        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetAdvectionOnWallTest3D(
        //    [Values(Math.PI / 4, Math.PI / 3, Math.PI / 6)] double contactAngle,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(0, 1, 2)] int AMRlevel,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {

        //    var Tst = new LevelSetAdvectionOnWallTest(3, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling, 1, AMRlevel);
        //    C.SkipSolveAndEvaluateResidual = true;

        //    LevelSetTest(Tst, C);

        //}

        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetScalingTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetScalingTest(
        //    [Values(2)] int spatialDimension,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling)
        //    //[Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timeSteppingScheme) 
        //    {
        //    // Todo: singleInit/multiInit, 

        //    var Tst = new LevelSetScalingTest(spatialDimension, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling);

        //    LevelSetTest(Tst, C);

        //}


        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetRotationTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetRotationTest(
        //    [Values(2)] int spatialDimension,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling)
        //    //[Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timeSteppingScheme) 
        //    {
        //    // Todo: singleInit/multiInit, 

        //    var Tst = new LevelSetRotationTest(spatialDimension, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling);

        //    LevelSetTest(Tst, C);

        //}




        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        public static void LevelSetTest(ILevelSetTest Tst, SolverWithLevelSetUpdaterTestControl C, string IO = null) {

            using (var solver = new SolverWithLevelSetUpdaterTestCenter()) {

                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;

                AgglomerationAlgorithm.Katastrophenplot = delegate (DGField[] plotFields) {

                    List<DGField> allfields = new();
                    allfields.AddRange(plotFields);
                    
                    foreach(var f in solver.RegisteredFields) {
                        if(!allfields.Contains(f, (a, b) => object.ReferenceEquals(a, b)))
                            allfields.Add(f);
                    }

                    Tecplot.PlotFields(allfields, "AgglomFail", 0.0, 4);
                };

                solver.Init(C);
                solver.RunSolverMode();
                //solver.OperatorAnalysis();

                //-------------------Evaluate Error ---------------------------------------- 
                var evaluator = new LevelSetErrorEvaluator<SolverWithLevelSetUpdaterTestControl>(solver);
                var LastErrors = evaluator.ComputeErrors(C.Endtime);

                // output
                var ErrThresh = Tst.AcceptableError;
                for (int iLevSet = 0; iLevSet < solver.Control.NoOfLevelSets; iLevSet++) {
                    string levSetName = VariableNames.LevelSetCGidx(iLevSet);

                    var spcIds = solver.LsTrk.SpeciesNames;
                    string[] ErrNames = new string[spcIds.Count + 2 + 2];
                    int j = 0;
                    ErrNames[j++] = levSetName;
                    ErrNames[j++] = VariableNames.LevelSetDGidx(iLevSet);
                    foreach(var spc in spcIds) {
                        ErrNames[j++] = "Volume#" + spc;
                    }
                    ErrNames[j++] = VariableNames.AsLevelSetVariable(levSetName, "Surface-Error");
                    ErrNames[j++] = VariableNames.AsLevelSetVariable(levSetName, "Contour-Error");

                    for (int i = 0; i < ErrThresh.GetLength(1); i++) { 
                        Console.WriteLine("L2 error, '{0}': \t{1}", ErrNames[i], LastErrors[iLevSet, i]);
                        
                    }
                }

                // Assertions
                for (int iLevSet = 0; iLevSet < solver.Control.NoOfLevelSets; iLevSet++) {
                    string levSetName = VariableNames.LevelSetCGidx(iLevSet);

                    var spcIds = solver.LsTrk.SpeciesNames;
                    string[] ErrNames = new string[spcIds.Count + 2 + 2];
                    int j = 0;
                    ErrNames[j++] = levSetName;
                    ErrNames[j++] = VariableNames.LevelSetDGidx(iLevSet);
                    foreach (var spc in spcIds) {
                        ErrNames[j++] = "Volume#" + spc;
                    }
                    ErrNames[j++] = VariableNames.AsLevelSetVariable(levSetName, "Surface-Error");
                    ErrNames[j++] = VariableNames.AsLevelSetVariable(levSetName, "Contour-Error");

                    for (int i = 0; i < ErrThresh.GetLength(1); i++) {
                        Assert.LessOrEqual(LastErrors[iLevSet, i], ErrThresh[iLevSet, i]);
                    }
                }

                //-------------------check parallel simulations ---------------------------------------- 
                int RefMPIsize = 1;
                if (IO != null) {


                    var projCheck = new TestingIO(solver.GridData, $"{IO}.csv", true, RefMPIsize);
                    for (int iLevSet = 0; iLevSet < solver.Control.NoOfLevelSets; iLevSet++) {
                        LevelSet PhiDG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCGidx(iLevSet)].DGLevelSet;
                        LevelSet PhiCG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCGidx(iLevSet)].CGLevelSet;
                        projCheck.AddDGField(PhiDG);
                        projCheck.AddDGField(PhiCG);
                    }
                    projCheck.DoIOnow();
                    for (int iLevSet = 0; iLevSet < solver.Control.NoOfLevelSets; iLevSet++) {
                        LevelSet PhiDG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCGidx(iLevSet)].DGLevelSet;
                        LevelSet PhiCG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCGidx(iLevSet)].CGLevelSet;
                        Assert.Less(projCheck.AbsError(PhiDG), 1.0e-15, "Mismatch in projected PhiDG between single-core and parallel run.");
                        Assert.Less(projCheck.AbsError(PhiCG), 1.0e-15, "Mismatch in projected PhiCG between single-core and parallel run.");
                    }
                }

            }
        }


        static SolverWithLevelSetUpdaterTestControl LSTstObj2CtrlObj(ILevelSetTest tst, int maxNoTimesteps,
            LevelSetEvolution levelSetEvolution, LevelSetHandling levelSetHandling, 
            int GridResolution = 1, int AMRlevel = 0) {

            SolverWithLevelSetUpdaterTestControl C = new SolverWithLevelSetUpdaterTestControl();
            int D = tst.SpatialDimension;

            // database setup
            // ==============
            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "LevelSet/" + tst.GetType().Name;
            C.ProjectDescription = "Test";


            // DG degree
            // =========
            C.NoOfLevelSets = tst.NoOfLevelsets;
            C.DegreeOfLevelSets = tst.LevelsetPolynomialDegree;
            C.SetDGdegree(tst.LevelsetPolynomialDegree);


            // grid
            // ====
            C.GridFunc = () => tst.CreateGrid(GridResolution);


            // boundary conditions
            // ===================
            foreach (var kv in tst.GetBoundaryConfig()) {
                C.BoundaryValues.Add(kv);
            }

            // initial values and exact solution
            // =================================
            for (int i = 0; i < tst.NoOfLevelsets; i++) {
                C.SetAdvectionVelocity(i, tst.GetU()[i]);
                C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(i), tst.GetPhi()[i]);
            }


            // timestepping and solver
            // =======================
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.Option_LevelSetEvolution = levelSetEvolution;
            C.Timestepper_LevelSetHandling = levelSetHandling;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler; // should not really matter, we are only projecting the exact underlying advection velocity in each timestep

            // adaptive mesh refinement
            if (AMRlevel > 0) {
                C.AdaptiveMeshRefinement = true;
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = AMRlevel });
                C.AMR_startUpSweeps = AMRlevel;
            }

            double dt = tst.ComputeTimestep(GridResolution, tst.LevelsetPolynomialDegree, AMRlevel);
            C.dtFixed = dt;
            double T = tst.getEndTime();
            int timesteps = (int)(T / dt);
            C.NoOfTimesteps = (timesteps > maxNoTimesteps) ? maxNoTimesteps : timesteps;
            C.Endtime = dt * C.NoOfTimesteps;

            return C;
        }


    }

    /// <summary>
    /// error evaluator for LevelSet-based tests 
    /// computes error fields for: Phi, PhiDG, gradient of PhiDG 
    /// integral values: area, length
    /// </summary>
    public class LevelSetErrorEvaluator<T> where T : SolverWithLevelSetUpdaterTestControl, new() {

        protected SolverWithLevelSetUpdater<T> solver;
        public LevelSetErrorEvaluator(SolverWithLevelSetUpdater<T> solver) {
            this.solver = solver;
        }

        LevelSet[] exactPhi;
        LevelSetTracker exactLsTrk;

        void SetExactPhiAndLevelSetTracker(double time) {
            int NoLS = solver.Control.NoOfLevelSets;
            exactPhi = new LevelSet[NoLS];
            // exact level-set fields
            for (int iLevSet = 0; iLevSet < NoLS; iLevSet++) {
                Func<double[], double, double> phiExactFunc = solver.Control.InitialValues_Evaluators_TimeDep[VariableNames.LevelSetCGidx(iLevSet)];
                exactPhi[iLevSet] = new LevelSet(new Basis(solver.GridData, (solver.Control.DegreeOfLevelSets + 1) * 2), "exactLevelSet");
                exactPhi[iLevSet].Clear();
                exactPhi[iLevSet].ProjectField(NonVectorizedScalarFunction.Vectorize(phiExactFunc, time));
            }
            // exact level-set tracker
            exactLsTrk = new LevelSetTracker((GridData)solver.GridData,
            XQuadFactoryHelper.MomentFittingVariants.Saye, 1, solver.LsTrk.SpeciesNames.ToArray(), exactPhi);
            exactLsTrk.UpdateTracker(time);
        }


        /// <summary>
        /// computes the error against the continuous level set field "Phi"
        /// </summary>
        /// <param name="exactLevelSetFunc"></param>
        /// <param name="time"></param>
        /// <param name="cm"></param>
        /// <returns></returns>
        public double ComputeLevelSetError(int iLevSet, CellMask cm) {
            SinglePhaseField PhiCG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCGidx(iLevSet)].CGLevelSet;
            string name = PhiCG.Identification;

            double L2Error = PhiCG.L2Error(exactPhi[iLevSet], cm);

            solver.QueryHandler.ValueQuery("L2err_" + name, L2Error, true);

            return L2Error;
        }

        /// <summary>
        /// computes the error against the discontinuous level set field "PhiDG"
        /// </summary>
        /// <param name="exactLevelSetFunc"></param>
        /// <param name="time"></param>
        /// <param name="cm"></param>
        /// <returns></returns>
        public double ComputeDGLevelSetError(int iLevSet, CellMask cm) {
            SinglePhaseField PhiDG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCGidx(iLevSet)].DGLevelSet;
            string name = PhiDG.Identification;

            double L2Error = PhiDG.L2Error(exactPhi[iLevSet], cm);

            solver.QueryHandler.ValueQuery("L2err_" + name, L2Error, true);

            return L2Error;
        }

        /// <summary>
        /// computes the error of the interface points with respect to the given interface form
        /// </summary>
        /// <returns></returns>
        public double ComputeContourError(int iLevSet) {

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), solver.QuadOrder(), 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet, solver.LsTrk.Regions.GetCutCellMask4LevSet(iLevSet));
            int Norm = 2;
            var t_exactPhi = exactPhi[iLevSet].CloneAs();
            t_exactPhi.AccLaidBack(-1.0, solver.LsUpdater.LevelSets[VariableNames.LevelSetCGidx(iLevSet)].CGLevelSet);
            double conterr = DGField.IntegralOverEx(cqs, (X,U,j) => U[0], 2, t_exactPhi);

            return conterr;
        }


        /// <summary>
        /// computes the interface length in 2D and area in 3D
        /// </summary>
        /// <returns></returns>
        public double ComputeSurfaceError(int iLevSet) {

            double exactInterfaceSize = LevelSetUtils.GetInterfaceLength(iLevSet, exactLsTrk, solver.QuadOrder());
            double interfaceSize = LevelSetUtils.GetInterfaceLength(iLevSet, solver.LsTrk, solver.QuadOrder());

            return Math.Abs(interfaceSize - exactInterfaceSize);
        }

        /// <summary>
        /// computes the species area
        /// </summary>
        /// <returns></returns>
        public Dictionary<SpeciesId, double> ComputeSpeciesVolumeError() {

            Dictionary<SpeciesId, double> spcArea = new Dictionary<SpeciesId, double>();

            foreach (SpeciesId spcId in solver.LsTrk.SpeciesIdS) {
                double exactArea = LevelSetUtils.GetSpeciesArea(exactLsTrk, spcId, solver.QuadOrder());
                double area = LevelSetUtils.GetSpeciesArea(solver.LsTrk, spcId, solver.QuadOrder());
                double error = Math.Abs(area - exactArea);
                spcArea.Add(spcId, error);
            }

            return spcArea;
        }


        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object, given by the time dependent initial evaluators
        /// </summary>
        public double[,] ComputeErrors(double time) {

            SetExactPhiAndLevelSetTracker(time);

            var spcIds = solver.LsTrk.SpeciesIdS;

            double[,] Ret = new double[solver.Control.NoOfLevelSets, spcIds.Count + 2 + 2];
            var vol_errs = ComputeSpeciesVolumeError();
            for (int iLevSet = 0; iLevSet < solver.Control.NoOfLevelSets; iLevSet++) {
                int j = 0;

                CellMask cm = solver.LsTrk.Regions.GetCutCellMask4LevSet(iLevSet);
                Ret[iLevSet, j++] = ComputeLevelSetError(iLevSet, cm);
                Ret[iLevSet, j++] = ComputeDGLevelSetError(iLevSet, cm);

                for (int i = 0; i < spcIds.Count; i++) {
                    Ret[iLevSet, j++] = vol_errs[spcIds[i]];
                }

                Ret[iLevSet, j++] = ComputeSurfaceError(iLevSet);
                Ret[iLevSet, j++] = ComputeContourError(iLevSet);
            }

            return Ret;
        }

    }


}
