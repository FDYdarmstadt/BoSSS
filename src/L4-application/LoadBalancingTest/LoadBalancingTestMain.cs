using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.LoadBalancingTest {

    /// <summary>
    /// App which performs basic tests on dynamic load balancing.
    /// </summary>
    class LoadBalancingTestMain : BoSSS.Solution.Application {

        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

            BoSSS.Solution.Application._Main(
                args,
                true,
                () => new LoadBalancingTestMain());
        }

        public override void Init(BoSSS.Solution.Control.AppControl control) {
            control.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;
            base.Init(control);
        }

        protected override GridCommons CreateOrLoadGrid() {
            double[] nodes = GenericBlas.Linspace(-5, 5, 21);
            var grd = Grid2D.Cartesian2DGrid(nodes, nodes);
            this.Control.NoOfMultigridLevels = 1; // required for XDG-BDF timestepper
            return grd;
        }
              

        XDGField u;
        XDGField uResidual;
        XDGField uEx;

        SinglePhaseField Amarker;
        SinglePhaseField Bmarker;

        SinglePhaseField MPICellRank;

        LevelSet LevSet;

        protected override void CreateFields() {
            LevSet = new LevelSet(new Basis(this.GridData, 2), "LevelSet");
            base.LsTrk = new LevelSetTracker(this.GridData, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 1, new string[] { "A", "B" }, LevSet);

            var xBasis = new XDGBasis(base.LsTrk, DEGREE);
            u = new XDGField(xBasis, "u");
            uResidual = new XDGField(xBasis, "Res");
            uEx = new XDGField(xBasis, "uEx");

            Amarker = new SinglePhaseField(new Basis(this.GridData, 0), "Amarker");
            Bmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Bmarker");
            MPICellRank = new SinglePhaseField(new Basis(this.GridData, 0), "MPIRank");
            base.m_RegisteredFields.Add(LevSet);
            base.m_RegisteredFields.Add(u);
            base.m_RegisteredFields.Add(uEx);
            base.m_RegisteredFields.Add(uResidual);
            base.m_RegisteredFields.Add(Amarker);
            base.m_RegisteredFields.Add(Bmarker);
            base.m_RegisteredFields.Add(MPICellRank);
        }

        /// <summary>
        /// turn dynamic balancing on/off
        /// </summary>
        internal bool DynamicBalance = true;

        /// <summary>
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE = 3;

        
        /// <summary>
        /// Cell Agglomeration threshold
        /// </summary>
        internal double THRESHOLD = 0.1;

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
            LsTrk.UpdateTracker(incremental: _incremental);

            // exact solution for new timestep
            uEx.GetSpeciesShadowField("A").ProjectField((x, y) => x + alpha_A * t);
            uEx.GetSpeciesShadowField("B").ProjectField((x, y) => x + alpha_B * t);

            u.Clear();
            u.Acc(1.0, uEx);

            // markieren, wo ueberhaupt A und B sind
            Amarker.Clear();
            Bmarker.Clear();
            Amarker.AccConstant(+1.0, LsTrk.Regions.GetSpeciesSubGrid("A").VolumeMask);
            Bmarker.AccConstant(1.0, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);

            // MPI rank
            MPICellRank.Clear();
            MPICellRank.AccConstant(base.MPIRank);

            // return level-set residual
            return 0.0;
        }

        /// <summary>
        /// Setting initial value.
        /// </summary>
        protected override void SetInitial() {
            this.DelUpdateLevelset(null, 0.0, 0.0, 0.0, false);
        }

        /// <summary>
        /// The operator d/dx
        /// </summary>
        XSpatialOperator Op;

        /// <summary>
        /// The BDF time integrator - makes load balancing challenging.
        /// </summary>
        XdgBDFTimestepping TimeIntegration;

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            int quadorder = this.u.Basis.Degree * 2 + 1;

            Op = new XSpatialOperator(1, 0, 1, (A, B, C) => quadorder, "u", "c1");

            var blkFlux = new DxFlux(this.LsTrk, alpha_A, alpha_B);
            Op.EquationComponents["c1"].Add(blkFlux); // Flux in Bulk Phase;
            Op.EquationComponents["c1"].Add(new LevSetFlx(this.LsTrk, alpha_A, alpha_B)); // flux am lev-set 0
            
            Op.Commit();

            if (L == null) {
                TimeIntegration = new XdgBDFTimestepping(
                    new DGField[] { u }, new DGField[] { uResidual }, base.LsTrk,
                    true,
                    DelComputeOperatorMatrix, DelUpdateLevelset,
                    3, // BDF3
                       //-1, // Crank-Nicolson
                       //0, // Explicit Euler
                    LevelSetHandling.LieSplitting,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    SpatialOperatorType.LinearTimeDependent,
                    MassScale,
                    MultigridOperatorConfig,
                    this.MultigridSequence,
                    this.LsTrk.SpeciesIdS.ToArray(),
                    quadorder,
                    this.THRESHOLD,
                    true);
            } else {
                Debug.Assert(object.ReferenceEquals(this.MultigridSequence[0].ParentGrid, this.GridData));
                TimeIntegration.DataRestoreAfterBalancing(L, new DGField[] { u }, new DGField[] { uResidual }, base.LsTrk, this.MultigridSequence);
            }
        }

        /// <summary>
        /// Block scaling of the mass matrix.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double[] _rho = new double[] { 1 };
                IDictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _rho);
                R.Add(this.LsTrk.GetSpeciesId("B"), _rho);
                return R;
            }
        }

        MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int pu = this.u.Basis.Degree;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration entry will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[this.MultigridSequence.Length][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[1];


                    // configuration for pressure
                    configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(1, pu - iLevel),
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                        VarIndex = new int[] { 0 }
                    };

                }

                return configs;
            }
        }


        const XQuadFactoryHelper.MomentFittingVariants HMF = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

        protected virtual void DelComputeOperatorMatrix(BlockMsrMatrix OpMatrix, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {
            OpMatrix.Clear();
            OpAffine.ClearEntries();

            Op.ComputeMatrixEx(base.LsTrk,
                u.Mapping, null, uResidual.Mapping,
                OpMatrix, OpAffine, false,
                phystime,
                false,
                base.LsTrk.SpeciesIdS.ToArray());
        }

        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
         
            TimeIntegration.DataBackupBeforeBalancing(L);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                // Elo
                Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime + " ... ");

                // timestepping params
                base.NoOfTimesteps = 20;
                dt = 0.1;

                // initial value
                if (TimestepNo == 1)
                    TimeIntegration.SingleInit();

                // compute one time-step
                TimeIntegration.Solve(phystime, dt, ComputeOnlyResidual: true);
                double ResidualNorm = this.uResidual.L2Norm();
                Assert.LessOrEqual(ResidualNorm, 1.0e-8, "Unusually large and ugly residual norm detected.");

                // done.
                Console.WriteLine("    done.");

                // return
                //base.TerminationKey = true;
                return dt;
            }
        }

        int m_NoOfReparts = 0;

        /// <summary>
        /// Check is load-balancing was actually tested.
        /// </summary>
        protected override void Bye() {
            if(MPISize > 1) {
                int TotalNoOfReparts = m_NoOfReparts.MPISum();
                Assert.Greater(TotalNoOfReparts, 0, "Load-balancing was not tested at all.");
            }
        }


        internal Func<IApplication, int, ICellCostEstimator> cellCostEstimatorFactory = CellCostEstimatorLibrary.OperatorAssemblyAndCutCellQuadrules;


        /// <summary>
        /// A dummy routine in order to test cell dynamic load balancing 
        /// (**not** a good balancing, but triggers redistribution).
        /// </summary>
        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            if (!DynamicBalance || MPISize <= 1)
                return null;

            //if(MPIRank == 0)
            //    Debugger.Launch();
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] NewPart;
            
            int[] PerformanceClasses = new int[J];
            var CC = this.LsTrk.Regions.GetCutCellMask();
            foreach(int j in CC.ItemEnum) {
                PerformanceClasses[j] = 1;
            }
            int NoCutCells = CC.NoOfItemsLocally.MPISum();
            Console.WriteLine("Number of cut cells: " + NoCutCells);

            if (balancer == null) {
                balancer = new LoadBalancer(
                    new List<Func<IApplication, int, ICellCostEstimator>>() { cellCostEstimatorFactory }
                    );
            }

            NewPart = balancer.GetNewPartitioning(
                this,
                2,
                PerformanceClasses,
                TimeStepNo,
                GridPartType.none,
                "",
                imbalanceThreshold: 0.0,
                Period: 3,
                redistributeAtStartup: false);

            

            /*

            if (MPISize == 4 && TimeStepNo > 5) {
                NewPart = new int[J];

                double speed = this.GridData.Cells.h_maxGlobal / 0.3;
                for (int j = 0; j < J; j++) {
                    double[] X = this.GridData.Cells.CellCenter.GetRow(j);
                    double x = X[0];
                    double y = X[1];

                    int px = x < Math.Min(physTime * speed, 2) ? 0 : 1;
                    int py = y < 0 ? 0 : 1;

                    NewPart[j] = px * 2 + py;
                }

                //return Part;
            } else if (MPISize == 2) {
                NewPart = new int[J];

                double speed = this.GridData.Cells.h_maxGlobal / 0.3;

                for (int j = 0; j < J; j++) {
                    double[] X = this.GridData.Cells.CellCenter.GetRow(j);
                    double x = X[0];
                    double y = X[1];

                    int px = x < Math.Min(physTime * speed, 2) ? 0 : 1;
                    int py = y < 0 ? 0 : 1;

                    NewPart[j] = px;
                }

                //return null;
                //return Part;
            } else if (MPISize == 1) {
                
                return null;
            } else {
                return null;
            }

            //*/


            if(NewPart != null) {
                int myRank = this.MPIRank;
                foreach(int tr in NewPart) {
                    if(tr != myRank)
                        m_NoOfReparts++;
                }
            }
            
            return NewPart;
        }

        LoadBalancer balancer;

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "LoadBalancingTest." + timestepNo;
            Tecplot.PlotFields(base.m_RegisteredFields.ToArray(), filename, physTime, superSampling);
        }
    }
}
