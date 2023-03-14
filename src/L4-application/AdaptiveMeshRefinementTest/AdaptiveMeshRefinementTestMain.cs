﻿using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Tecplot;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Collections;
using NUnit.Framework;
using MPI.Wrappers;

namespace BoSSS.Application.AdaptiveMeshRefinementTest {

    /// <summary>
    /// App which performs basic tests on adaptive mesh refinement and dynamic load balancing.
    /// </summary>
    class AdaptiveMeshRefinementTestMain : BoSSS.Solution.Application {
        
        static void Main(string[] args) {
            BoSSS.Solution.Application._Main(
                args,
                true,
                () => new AdaptiveMeshRefinementTestMain());
            //AllUpTest.SetUp();
            //AllUpTest.RuntimeCostDynamicBalanceTest(2);
            //AllUpTest.OneTimeTearDown();
        }

        public int TestCase = 2;

        protected override IGrid CreateOrLoadGrid() {
            GridCommons grd;
            switch (TestCase) {
                case 1: {
                    // ++++++++++++++++++
                    // affine-linear mesh
                    // ++++++++++++++++++
                    double[] xNodes = GenericBlas.Linspace(-5, 5, 21);
                    double[] yNodes = GenericBlas.Linspace(-5, 5, 21);
                    grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                    break;
                }

                case 2: {
                    // +++++++++++
                    // curved mesh
                    // +++++++++++
                    double[] xNodes = GenericBlas.Linspace(-5, 5, 21);
                    double[] yNodes = GenericBlas.Linspace(-5, 5, 21);
                    grd = Grid2D.BilinearSquareGrid(xNodes, yNodes, factor: 0.8);
                    break;
                }

                default:
                throw new ArgumentException("unknown test-case index");
            }

            return grd;
        }

        public override void Init(Solution.Control.AppControl control) {
            //control.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;
            control.NoOfMultigridLevels = 1;
            base.Init(control);
        }

        SinglePhaseField TestData;
        SinglePhaseField u;
        VectorField<SinglePhaseField> Grad_u;
        SinglePhaseField MagGrad_u;

        XDGField uX;
        XDGField uXResidual;
        XDGField uXEx;
        LevelSet LevSet;

        protected override void CreateFields() {
            var Basis = new Basis(base.GridData, DEGREE);
            u = new SinglePhaseField(Basis, "u");
            Grad_u = new VectorField<SinglePhaseField>(base.GridData.SpatialDimension, Basis, "Grad_u", (bs, nmn) => new SinglePhaseField(bs, nmn));
            MagGrad_u = new SinglePhaseField(new Basis(base.GridData, 0), "Magnitude_Grad_u");
            TestData = new SinglePhaseField(Basis, "TestData");
            base.m_RegisteredFields.Add(u);
            base.m_RegisteredFields.AddRange(Grad_u);
            base.m_RegisteredFields.Add(MagGrad_u);
            base.m_RegisteredFields.Add(TestData);

            LevSet = new LevelSet(new Basis(this.GridData, 2), "LevelSet");
            base.LsTrk = new LevelSetTracker((BoSSS.Foundation.Grid.Classic.GridData) this.GridData, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 1, new string[] { "A", "B" }, LevSet);
            base.m_RegisteredFields.Add(LevSet);

            var xBasis = new XDGBasis(base.LsTrk, DEGREE);

            uX = new XDGField(xBasis, "uX");
            uXResidual = new XDGField(xBasis, "ResX");
            uXEx = new XDGField(xBasis, "uXEx");

            uX.UpdateBehaviour = BehaveUnder_LevSetMoovement.JustReallocate;
            uXResidual.UpdateBehaviour = BehaveUnder_LevSetMoovement.JustReallocate;
            uXEx.UpdateBehaviour = BehaveUnder_LevSetMoovement.JustReallocate;

            base.m_RegisteredFields.Add(uX);
            base.m_RegisteredFields.Add(uXResidual);
            base.m_RegisteredFields.Add(uXEx);
        }


        /// <summary>
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE = 4;

        /// <summary>
        /// Setting initial value.
        /// </summary>
        protected override void SetInitial(double t) {
            TestData.ProjectField(TestDataFunc);
            DelUpdateLevelset(new DGField[] { uX }, 0.0, 0.0, 1.0, false);

            uX.ProjectField((x, y) => 1.0);
            
            // check error
            double L2err = TestData.L2Error(TestDataFunc);
            Console.WriteLine("Projection error from old to new grid: " + L2err);
            Assert.LessOrEqual(L2err, 1.0e-8, "Projection error of test field to high.");


          
        }


        /// <summary>
        /// Equation coefficient, species A.
        /// </summary>
        const double alpha_A = 0.1;

        /// <summary>
        /// Equation coefficient, species B.
        /// </summary>
        const double alpha_B = 3.5;


        /// <summary>
        /// Sets level-set and solution at time (<paramref name="time"/> + <paramref name="dt"/>).
        /// </summary>
        double DelUpdateLevelset(DGField[] CurrentState, double time, double dt, double UnderRelax, bool _incremental) {

            // new time
            double t = time + dt;

            // project new level-set
            double s = 1.0;
            LevSet.ProjectField((x, y) => -(x - s * t).Pow2() - y.Pow2() + (2.4).Pow2());
            LsTrk.UpdateTracker(t, incremental: _incremental);
            LsTrk.PushStacks();

            // exact solution for new timestep
            uXEx.GetSpeciesShadowField("A").ProjectField((x, y) => x + alpha_A * t);
            uXEx.GetSpeciesShadowField("B").ProjectField((x, y) => x + alpha_B * t);

           // uX.Clear();
            //uX.Acc(1.0, uXEx);

            // single-phase-properties
            u.Clear();
            u.ProjectField(NonVectorizedScalarFunction.Vectorize(uEx, t));
            
            Grad_u.Clear();
            Grad_u.Gradient(1.0, u);

            MagGrad_u.Clear();
            MagGrad_u.ProjectFunction(1.0,
                (ilPSP.Vector X, double[] U, int jCell) => Math.Sqrt(U[0].Pow2() + U[1].Pow2()),
                new Foundation.Quadrature.CellQuadratureScheme(),
                Grad_u.ToArray());

            // return level-set residual
            return 0.0;
        }


        /// <summary>
        /// A polynomial expression which is projected onto <see cref="TestData"/>; since it is polynomial,
        /// the data should remain constant under refinement and coarsening.
        /// </summary>
        static double TestDataFunc(double x, double y) {
            return (x * x + y * y);
        }

        /// <summary>
        /// Function which is projected onto <see cref="u"/>;
        /// </summary>
        static double uEx(double[] X, double time) {
            double x = X[0];
            double y = X[1];
            
            double r = Math.Sqrt(x * x + y * y);
            double alpha = Math.Atan2(y, x);

            double alpha0 = alpha - time;

            double x0 = Math.Cos(alpha0) * r;
            double y0 = Math.Sin(alpha0) * r;
            
            double val = Math.Exp(- x0.Pow2() - (y0 - 2.5).Pow2()); 

            return val;
        }


        


        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            
            if (L == null) {
                
            } else {
                //throw new NotImplementedException("todo");
                //PlotCurrentState(100, new TimestepNumber(100, 1), 2);
            }
        }

        //bool LevelReminder = false;

        /// <summary>
        /// Very primitive refinement indicator, works on a gradient criterion.
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel) {
            double GradMag = MagGrad_u.GetMeanValue(j);

            int DesiredLevel_j = 0;
            if(GradMag > 0.6)
                DesiredLevel_j = 1;
            if(GradMag > 0.81)
                DesiredLevel_j = 2;

            //if(DesiredLevel_j < CurrentLevel) {
            //    DesiredLevel_j = CurrentLevel;
            //    if(!LevelReminder) {
            //        Console.WriteLine("Reminder: coarsening disabled");
            //        LevelReminder = true;
            //    }
            //}

            return DesiredLevel_j;
            
        }

        /// <summary>
        /// Creates the cellmask which should be refined.
        /// </summary>
        private int[] GetCellsToRefine() {
            int refinementLevel = 2;
            int coarseRefinementLevel = 1;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellsToRefine = new int[noOfLocalCells];
            for (int j = 0; j < noOfLocalCells; j++) {
                double GradMag = MagGrad_u.GetMeanValue(j);
                if (GradMag > 0.81)
                    cellsToRefine[j] = refinementLevel;
                else if (GradMag > 0.6)
                    cellsToRefine[j] = coarseRefinementLevel;
            }
            return cellsToRefine;
        }

        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            if(TimestepNo > 3 && TimestepNo % 3 != 0) {
            //{ 
                newGrid = null;
                old2NewGrid = null;
                return;
            }

            // Check grid changes
            // ==================
            bool AnyChange;
            List<int> CellsToRefineList;
            List<int[]> Coarsening;
            GridRefinementController gridRefinementController = new GridRefinementController((GridData)this.GridData, LsTrk.Regions.GetCutCellMask());
            if (MPISize > 1) {
                AnyChange = gridRefinementController.ComputeGridChange(GetCellsToRefine(), out CellsToRefineList, out Coarsening);
            }
            else
                AnyChange = gridRefinementController.ComputeGridChange(LevelIndicator, out CellsToRefineList, out Coarsening);
            int NoOfCellsToRefine = 0;
            int NoOfCellsToCoarsen = 0;
            if(AnyChange) {
                int[] glb = (new int[] {
                    CellsToRefineList.Count,
                    Coarsening.Sum(L => L.Length),
                }).MPISum();
                NoOfCellsToRefine = glb[0];
                NoOfCellsToCoarsen = glb[1];
            }
            long oldJ = this.GridData.CellPartitioning.TotalLength;

            // Update Grid
            // ===========
           
            if(AnyChange) {

                

                Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");


                newGrid = ((GridData)this.GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);
                
            } else {
                newGrid = null;
                old2NewGrid = null;
            }
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                // Elo
                Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime + " ... ");

                // timestepping params
                base.NoOfTimesteps = 100;
                dt = Math.PI*2 / base.NoOfTimesteps;
                //base.NoOfTimesteps = 20;

                //UpdateBaseGrid(phystime + dt);

                //UpdateRefinedGrid(phystime + dt, TimestepNo);
                DelUpdateLevelset(new DGField[] { uX }, phystime, dt, 1.0, false);

                // check error
                double L2err = TestData.L2Error(TestDataFunc);

                Console.WriteLine("Projection error from old to new grid: " + L2err);
                Assert.LessOrEqual(L2err, 1.0e-8, "Projection error from old to new grid to high.");

                // return
                //base.TerminationKey = true;
                return dt;
            }
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "AdaptiveMeshRefinementTest." + timestepNo;
            Tecplot.PlotFields(base.m_RegisteredFields.ToArray(), filename, physTime, superSampling);


            //DGField[] RefinedFields = new[] { Refined_u, Refined_TestData, Refined_Grad_u[0], Refined_Grad_u[1], Refined_MagGrad_u };
            //string filename2 = "RefinedGrid." + timestepNo;
            //Tecplot.PlotFields(RefinedFields, filename2, physTime, superSampling);
        }
    }
}
