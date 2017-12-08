using System;
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System.Linq;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution;
using BoSSS.Solution.XdgTimestepping;
using System.Collections.Generic;
using BoSSS.Solution.Multigrid;
using System.Collections;
using ilPSP.Tracing;

namespace BoSSS.Application.RefineAndLoadBal {

    /// <summary>
    /// App which performs basic tests on dynamic load balancing.
    /// </summary>
    class RefineAndLoadBalTestMain : BoSSS.Solution.Application<RefineAndLoadBalControl> {

        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

            BoSSS.Solution.Application<RefineAndLoadBalControl>._Main(
                args,
                true,
                null,
                () => new RefineAndLoadBalTestMain());
        }

        protected override GridCommons CreateOrLoadGrid() {
            double[] nodes = GenericBlas.Linspace(-5, 5, 21);
            var grd = Grid2D.Cartesian2DGrid(nodes, nodes);
            base.m_GridPartitioningType = GridPartType.none;
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
            base.LsTrk = new LevelSetTracker(this.GridData, 1, new string[] { "A", "B" }, LevSet);

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
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE = 3;

        /// <summary>
        /// Cell Agglomeration threshold
        /// </summary>
        internal double THRESHOLD = 0.1;


        /// <summary>
        /// Sets level-set and solution at time (<paramref name="time"/> + <paramref name="dt"/>).
        /// </summary>
        double DelUpdateLevelset(DGField[] CurrentState, double time, double dt, double UnderRelax, bool _incremental) {

            // new time
            double t = time + dt;

            // project new level-set
            LevSet.ProjectField(this.Control.LevelSet.Vectorize(t));
            LsTrk.UpdateTracker(incremental: _incremental);

            // exact solution for new timestep
            uEx.GetSpeciesShadowField("A").ProjectField(NonVectorizedScalarFunction.Vectorize(this.Control.uEx_A, t));
            uEx.GetSpeciesShadowField("B").ProjectField(NonVectorizedScalarFunction.Vectorize(this.Control.uEx_B, t));

            u.Clear();
            u.Acc(1.0, uEx);
            
            // markieren, wo ueberhaupt A und B sind
            Amarker.Clear();
            Bmarker.Clear();
            Amarker.AccConstant(+1.0, LsTrk._Regions.GetSpeciesSubGrid("A").VolumeMask);
            Bmarker.AccConstant(1.0, LsTrk._Regions.GetSpeciesSubGrid("B").VolumeMask);

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

        XSpatialOperator Op;

        XdgBDFTimestepping TimeIntegration; 

        protected override void CreateEquationsAndSolvers(LoadbalData L) {
            Op = new XSpatialOperator(1, 0, 1, QuadOrderFunc.SumOfMaxDegrees(RoundUp: true), "u", "c1");

            var blkFlux = new DxFlux(this.LsTrk, Control.alpha_A, Control.alpha_B);
            Op.EquationComponents["c1"].Add(blkFlux); // Flux in Bulk Phase;
            Op.EquationComponents["c1"].Add(new LevSetFlx(this.LsTrk, Control.alpha_A, Control.alpha_B)); // flux am lev-set 0
            Op.OnIntegratingBulk += blkFlux.NowIntegratingBulk;
            
            Op.Commit();

            if (L == null) {
                TimeIntegration = new XdgBDFTimestepping(
                    new DGField[] { u }, new DGField[] { uResidual }, base.LsTrk,
                    true,
                    DelComputeOperatorMatrix, DelUpdateLevelset, DelUpdateCutCellMetrics,
                    3, // BDF3
                    //-1, // Crank-Nicolson
                    //0, // Explicit Euler
                    LevelSetHandling.LieSplitting,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    SpatialOperatorType.LinearTimeDependent,
                    MassScale,
                    MultigridOperatorConfig,
                    this.MultigridSequence,
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

        CutCellMetrics DelUpdateCutCellMetrics() {
            return new CutCellMetrics(
                HMF, 
                this.Op.QuadOrderFunction(new int[] { u.Basis.Degree }, new int[0], new int[] { uResidual.Basis.Degree }), 
                this.LsTrk, 
                this.LsTrk.SpeciesIdS.ToArray());
        }

        const XQuadFactoryHelper.MomentFittingVariants HMF = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

        protected virtual void DelComputeOperatorMatrix(BlockMsrMatrix OpMatrix, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, MultiphaseCellAgglomerator Agglomerator, double phystime) {
            OpMatrix.Clear();
            OpAffine.ClearEntries();

            Op.ComputeMatrixEx(base.LsTrk,
                u.Mapping, null, uResidual.Mapping,
                OpMatrix, OpAffine, false, 
                phystime,
                false,
                XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 
                base.LsTrk.SpeciesIdS.ToArray());
        }

        public override void DataBackupBeforeBalancing(LoadbalData L) {
            TimeIntegration.DataBackupBeforeBalancing(L);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                // Elo
                Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime + " ... ");

                // timestepping params
                dt = Control.dtFixed;

                // initial value
                if (TimestepNo == 1)
                    TimeIntegration.SingleInit();

                // compute one time-step
                TimeIntegration.Solve(phystime, dt, ComputeOnlyResidual: true);
                double ResidualNorm = this.uResidual.L2Norm();

                // done.
                Console.WriteLine("    done.");

                // return
                //base.TerminationKey = true;
                return dt;
            }
        }

        /// <summary>
        /// Assigns two performance classes (of cells)
        /// - normal, non-cut cells: 0
        /// - cut-cells: 1
        /// </summary>
        override protected void GetCellPerformanceClasses(out int NoOfClasses, out int[] CellPerfomanceClasses) {
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;

            CellPerfomanceClasses = new int[J];
            var CC = this.LsTrk._Regions.GetCutCellMask();
            foreach (int j in CC.ItemEnum)
                CellPerfomanceClasses[j] = 1;
            NoOfClasses = 2;
        }

 
        /// <summary>
        /// Ususal plotting.
        /// </summary>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "LoadBalancingTest." + timestepNo;
            Tecplot.PlotFields(base.m_RegisteredFields.ToArray(), filename, physTime, superSampling);
        }
    }
}
