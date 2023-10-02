using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
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
using BoSSS.Solution.LoadBalancing;

namespace BoSSS.Application.LoadBalancingTest {

    /// <summary>
    /// App which performs basic tests on dynamic load balancing.
    /// </summary>
    public class LoadBalancingTestMain : BoSSS.Solution.Application<AppControlSolver> {

        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

            //MultiphaseCellAgglomerator.Katastrophenplot = KatastrophenPlot;
            //InitMPI();
            // dbg_launch();
            //BoSSS.Application.LoadBalancingTest.AllUpTest.RuntimeCostDynamicBalanceTest(1);
            //throw new NotImplementedException("remove me");


            BoSSS.Solution.Application<AppControlSolver>._Main(
                args,
                true,
                () => new LoadBalancingTestMain());
        }


        public override void Init(BoSSS.Solution.Control.AppControl control) {
            control.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;
            base.Init(control);
        }

        protected override IGrid CreateOrLoadGrid() {
            double[] nodes = GenericBlas.Linspace(-5, 5, 21);
            var grd = Grid2D.Cartesian2DGrid(nodes, nodes);
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
            base.LsTrk = new LevelSetTracker((GridData)this.GridData, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 1, new string[] { "A", "B" }, LevSet);

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
        /// update of level-set
        /// </summary>
        class LevelSetTimestepping : ISlaveTimeIntegrator {

            public LevelSetTimestepping(LoadBalancingTestMain __owner) {
                m_owner = __owner;
            }
            LoadBalancingTestMain m_owner;

            public void Pop() {
                throw new NotImplementedException();
            }

            public void Push() {
                throw new NotImplementedException();
            }

            public double Update(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {

                return m_owner.DelUpdateLevelset(CurrentState, phystime, dt, UnderRelax, incremental);
            }
        }

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
        protected override void SetInitial(double t) {
            this.DelUpdateLevelset(null, t, 0.0, 0.0, false);
        }

        /// <summary>
        /// The operator d/dx
        /// </summary>
        XDifferentialOperatorMk2 Op;

        /// <summary>
        /// The BDF time integrator - makes load balancing challenging.
        /// </summary>
        XdgTimestepping TimeIntegration;

        //XdgBDFTimestepping AltTimeIntegration;

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            int quadorder = this.u.Basis.Degree * 2 + 1;

            Op = new XDifferentialOperatorMk2(1, 0, 1, (A, B, C) => quadorder, LsTrk.SpeciesNames, "u", "c1");

            var blkFlux = new DxFlux(this.LsTrk, alpha_A, alpha_B);
            Op.EquationComponents["c1"].Add(blkFlux); // Flux in Bulk Phase;
            Op.EquationComponents["c1"].Add(new LevSetFlx( alpha_A, alpha_B)); // flux am lev-set 0

            Op.LinearizationHint = LinearizationHint.AdHoc;
            Op.AgglomerationThreshold = this.THRESHOLD;
            Op.TemporalOperator = new ConstantXTemporalOperator(Op, 1.0);

            Op.Commit();

            if (L == null) {
                
                
                TimeIntegration = new XdgTimestepping(
                    Op,
                    new DGField[] { u }, new DGField[] { uResidual },
                    TimeSteppingScheme.BDF3,
                    () => new LevelSetTimestepping(this), LevelSetHandling.LieSplitting,
                    MultigridOperatorConfig, 
                    _AgglomerationThreshold: this.THRESHOLD,
                    LinearSolver: this.Control.LinearSolver, NonLinearSolver: this.Control.NonLinearSolver);
                


            } else {
                Debug.Assert(object.ReferenceEquals(this.MultigridSequence[0].ParentGrid, this.GridData));
                TimeIntegration.DataRestoreAfterBalancing(L, new DGField[] { u }, new DGField[] { uResidual }, base.LsTrk, this.MultigridSequence, this.Op);
                //AltTimeIntegration.DataRestoreAfterBalancing(L, new DGField[] { u }, new DGField[] { uResidual }, base.LsTrk, this.MultigridSequence);
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
                        DegreeS = new int[] { Math.Max(1, pu - iLevel) },
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                        VarIndex = new int[] { 0 }
                    };

                }

                return configs;
            }
        }

       
        protected virtual void DelComputeOperatorMatrix(BlockMsrMatrix OpMatrix, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {
            
            OpAffine.ClearEntries();

            bool Eval = false;
            if(OpMatrix == null) {
                Eval = true;
                OpMatrix = new BlockMsrMatrix(uResidual.Mapping, u.Mapping);
            } else {
                OpMatrix.Clear();
            }

            //Op.ComputeMatrixEx(base.LsTrk,
            //    u.Mapping, null, uResidual.Mapping,
            //    OpMatrix, OpAffine, false,
            //    phystime,
            //    false,
            //    base.LsTrk.SpeciesIdS.ToArray());
            XDifferentialOperatorMk2.XEvaluatorLinear mtxBuilder = Op.GetMatrixBuilder(base.LsTrk, u.Mapping, null, uResidual.Mapping);

            mtxBuilder.CellLengthScales.AddRange(AgglomeratedCellLengthScales);

            mtxBuilder.time = phystime;
            mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);

            if(Eval) {
                OpMatrix.SpMV(1.0, new CoordinateVector(CurrentState), 1.0, OpAffine);
            }

        }
        


        public override void DataBackupBeforeBalancing(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
         
            TimeIntegration.DataBackupBeforeBalancing(L);
            //AltTimeIntegration.DataBackupBeforeBalancing(L);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                // Elo
                Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime + " ... ");

                // timestepping params
                base.NoOfTimesteps = 20;
                dt = 0.1;


                // compute one time-step
                TimeIntegration.Solve(phystime, dt, SkipSolveAndEvaluateResidual: true);
                //AltTimeIntegration.Solve(phystime, dt, ComputeOnlyResidual: true);
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


        internal Func<ICellCostEstimator[]> cellCostEstimatorFactory = () => CellCostEstimatorLibrary.OperatorAssemblyAndCutCellQuadrules;


        /// <summary>
        /// A dummy routine in order to test cell dynamic load balancing 
        /// (**not** a good balancing, but triggers redistribution).
        /// </summary>
        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            if (!DynamicBalance || MPISize <= 1)
                return null;

            //if(MPIRank == 0)
            //     dbg_launch();
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
                balancer = new LoadBalancer(cellCostEstimatorFactory(), this);
            }

            NewPart = balancer.GetNewPartitioning(
                this,
                TimeStepNo,
                GridPartType.none,
                "",
                imbalanceThreshold: 0.0,
                Period: 3,
                redistributeAtStartup: false,
                TimestepNoRestart: null);

            

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
