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

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Timestepping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.LinSolvers;
using System.Diagnostics;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Foundation.Grid.Aggregation;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// Explicit or implicit timestepping using Runge-Kutta formulas,
    /// specialized for XDG applications.
    /// </summary>
    public class XdgRKTimestepping : XdgTimesteppingBase {

        public RungeKuttaScheme m_RKscheme {
            get;
            private set;
        }

        /// <summary>
        /// Constructor;
        /// </summary>
        /// <param name="Fields"></param>
        /// <param name="IterationResiduals"></param>
        /// <param name="LsTrk"></param>
        /// <param name="_ComputeOperatorMatrix">See <see cref="ComputeOperatorMatrix"/>.</param>
        /// <param name="_UpdateLevelset">See <see cref="UpdateLevelset"/>.</param>
        /// <param name="BDForder">
        /// The order of the BDF scheme from 1 to 6; in addition, 0 encodes Explicit Euler and -1 encodes Crank-Nicolson.
        /// </param>
        /// <param name="_LevelSetHandling"></param>
        /// <param name="_MassMatrixShapeandDependence"></param>
        /// <param name="_SpatialOperatorType"></param>
        /// <param name="_MassScale"></param>
        /// <param name="_AgglomerationThreshold"></param>
        /// <param name="_RKscheme"></param>
        /// <param name="useX">
        /// U dont want to know!
        /// </param>
        /// <param name="_MultigridSequence"></param>
        /// <param name="_CutCellQuadOrder">Order of quadrature in cut cells, required e.g. for <see cref="LevelSetTracker.GetXDGSpaceMetrics(SpeciesId[], int, int)"/></param>
        /// <param name="_SpId">Species to compute, actually a subset of <see cref="LevelSetTracker.SpeciesIdS"/></param>
        /// <param name="_MultigridOperatorConfig">
        /// Configuration of block-preconditioner, if null a default value is chosen.
        /// </param>
        public XdgRKTimestepping(DGField[] Fields,
            DGField[] IterationResiduals,
            LevelSetTracker LsTrk,
            DelComputeOperatorMatrix _ComputeOperatorMatrix,
            DelUpdateLevelset _UpdateLevelset,
            RungeKuttaScheme _RKscheme,
            LevelSetHandling _LevelSetHandling,
            MassMatrixShapeandDependence _MassMatrixShapeandDependence,
            SpatialOperatorType _SpatialOperatorType,
            IDictionary<SpeciesId, IEnumerable<double>> _MassScale,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
            AggregationGridData[] _MultigridSequence,
            SpeciesId[] _SpId,
            int _CutCellQuadOrder,
            double _AgglomerationThreshold, bool useX,
            Control.NonLinearSolverConfig nonlinconfig,
            Control.LinearSolverConfig linearconfig) : base (nonlinconfig, linearconfig) {

            // check args, set internals
            // -------------------------

            if (Fields.Length != IterationResiduals.Length)
                throw new ArgumentException("Expecting the same number of fields and residuals.");
            for (int iFld = 0; iFld < Fields.Length; iFld++) {
                if (!Fields[iFld].Basis.Equals(IterationResiduals[iFld].Basis))
                    throw new ArgumentException(string.Format("Mismatch between {0}-th basis of fields and residuals.", iFld));
            }

            if (_MassScale != null) {
                if (!IEnumerableExtensions.SetEquals(_SpId, _MassScale.Keys))
                    throw new ArgumentException();
            }

            base.Residuals = new CoordinateVector(IterationResiduals);

            if (!(_RKscheme.IsExplicit || _RKscheme.IsDiagonallyImplicit)) {
                throw new NotSupportedException("Only supporting explicit or diagonally implicit schemes.");
            }

            base.m_LsTrk = LsTrk;
            base.Config_LevelSetHandling = _LevelSetHandling;
            base.Config_MassMatrixShapeandDependence = _MassMatrixShapeandDependence;
            base.Config_SpatialOperatorType = _SpatialOperatorType;
            base.ComputeOperatorMatrix = _ComputeOperatorMatrix;
            base.UpdateLevelset = _UpdateLevelset;
            base.Config_MassScale = _MassScale;
            base.Config_AgglomerationThreshold = _AgglomerationThreshold;
            this.m_RKscheme = _RKscheme.CloneAs();
            base.MultigridSequence = _MultigridSequence;
            base.Config_SpeciesToCompute = _SpId;
            base.Config_CutCellQuadratureOrder = _CutCellQuadOrder;
            if (_MultigridSequence == null || _MultigridSequence.Length < 1)
                throw new ArgumentException("At least one grid level is required.");

            m_CurrentState = new CoordinateVector(Fields);

            if (_MultigridOperatorConfig != null) {
                Config_MultigridOperator = _MultigridOperatorConfig;
            } else {
                SetConfig_MultigridOperator_Default(Fields);
            }

            base.CommonConfigurationChecks();

            // configure stack of level-set-tracker
            // ------------------------------------

            if (Config_LevelSetHandling == LevelSetHandling.None) {
                m_LsTrk.IncreaseHistoryLength(0);
            } else if (Config_LevelSetHandling == LevelSetHandling.LieSplitting
                  || Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                m_LsTrk.IncreaseHistoryLength(1);
            } else {
                m_LsTrk.IncreaseHistoryLength(m_RKscheme.Stages);
            }
            //m_LsTrk.IncreaseHistoryLength(1);

            // multigrid - init
            // ----------------

            InitMultigrid(Fields, useX);
        }

        /// <summary>
        /// Actual solution.
        /// </summary>
        CoordinateVector m_CurrentState;

        protected override CoordinateMapping CurrentStateMapping {
            get {
                return m_CurrentState.Mapping;
            }
        }

        //double m_LastLevelSetResidual;


        private void MoveLevelSetAndRelatedStuff(DGField[] locCurSt, double PhysTime, double dt, double UnderRelax,
            BlockMsrMatrix[] MassMatrixStack,
            double[][] k) {


            // level-set evolution
            int oldPushCount = m_LsTrk.PushCount;
            int oldVersion = m_LsTrk.VersionCnt;
            m_LastLevelSetResidual = this.UpdateLevelset(locCurSt, PhysTime, dt, UnderRelax, (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting));
            int newVersion = m_LsTrk.VersionCnt;
            int newPushCount = m_LsTrk.PushCount;

            if ((newVersion - oldVersion) != 1)
                throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");
            if (oldPushCount != newPushCount) {
                throw new ApplicationException("Phushing the history stacks of the level-set tracker is reserved to the timestepper (during one timestep).");
            }


            if (MassMatrixStack != null && MassMatrixStack.Length > 1) {
                Debug.Assert(base.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative
                         || base.Config_LevelSetHandling == LevelSetHandling.Coupled_Once);

                var MMS0 = MassMatrixStack.GetSubVector(0, 1);
                TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, MMS0, CurrentStateMapping, CurrentStateMapping);
                MassMatrixStack[0] = MMS0[0];
                MassMatrixStack[1] = null; // must be re-computed according to new levset pos.
            }

            if (k != null) {
                Debug.Assert(base.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative
                         || base.Config_LevelSetHandling == LevelSetHandling.Coupled_Once);

                for (int i = 0; i < k.Length; i++) {
                    if (k[i] != null) {
                        TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, k[i], CurrentStateMapping);
                    }
                }
            }

            m_CurrentAgglomeration = null; // also invalid now
            m_PrecondMassMatrix = null; // also invalid now

        }

        //#if DEBUG
        //        int[][] m_OldSources;
        //        double[][][] m_OldVolumes;
        //#endif
        List<int> m_Versions = new List<int>();
        int m_RequiredTimeLevels = 0;

        void UpdateAgglom(bool ReplaceTop) {

            if (m_RequiredTimeLevels == 0)
                m_Versions.Clear();


            if (m_RequiredTimeLevels == 0 && ReplaceTop == true)
                throw new NotSupportedException();
            if (!ReplaceTop) {
                m_RequiredTimeLevels++;
                m_Versions.Add(m_LsTrk.Regions.Version);
            } else {
                m_Versions[m_Versions.Count - 1] = m_LsTrk.Regions.Version;
            }
            Debug.Assert(m_RequiredTimeLevels == m_Versions.Count);

            for (int i = 0; i < m_Versions.Count; i++) {
                if (m_Versions[m_Versions.Count - 1 - i] != m_LsTrk.RegionsHistory[1 - i].Version)
                    throw new ApplicationException("Internal Error, level-set-tracker history stack messed up."); // cheap test, also affordable in release
            }


            double[] oldAggTrsh;
            if (m_RequiredTimeLevels > 1) {
                oldAggTrsh = new double[m_RequiredTimeLevels - 1];
                ArrayTools.SetAll(oldAggTrsh, this.Config_AgglomerationThreshold);
            } else {
                oldAggTrsh = null;
            }
            Debug.Assert(m_LsTrk.PopulatedHistoryLength >= m_RequiredTimeLevels - 1);

            //#if DEBUG
            //            if(m_RequiredTimeLevels > 1) {
            //                Debug.Assert(m_OldSources != null);
            //            } else {
            //                m_OldSources = null;
            //                m_OldVolumes = null;
            //            }
            //#endif
            m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(Config_SpeciesToCompute, Config_CutCellQuadratureOrder,
                this.Config_AgglomerationThreshold,
                AgglomerateNewborn: oldAggTrsh != null, AgglomerateDecased: (oldAggTrsh != null), ExceptionOnFailedAgglomeration: true,
                oldTs__AgglomerationTreshold: oldAggTrsh);
            //#if DEBUG
            //            int[][] NewSources = new int[Config_SpeciesToCompute.Length][];
            //            {
            //                for(int iSpc = 0; iSpc < NewSources.Length; iSpc++) {
            //                    var Spc = Config_SpeciesToCompute[iSpc];
            //                    NewSources[iSpc] = m_CurrentAgglomeration.GetAgglomerator(Spc).AggInfo.SourceCells.ItemEnum.ToArray();
            //                }
            //            }
            //            double[][][] NewVolumes = new double[m_RequiredTimeLevels][][];
            //            for(int iTs = 0; iTs < m_RequiredTimeLevels; iTs++) {
            //                NewVolumes[iTs] = new double[Config_SpeciesToCompute.Length][];
            //                for(int iSpc = 0; iSpc < NewSources.Length; iSpc++) {
            //                    var Spc = Config_SpeciesToCompute[iSpc];
            //                    NewVolumes[iTs][iSpc] = m_LsTrk.GetXDGSpaceMetrics(Config_SpeciesToCompute, Config_CutCellQuadratureOrder, 1 - iTs).CutCellMetrics.CutCellVolumes[Spc].To1DArray();
            //                }
            //            }


            //            if(m_RequiredTimeLevels > 1) {
            //                for(int iSpc = 0; iSpc < NewSources.Length; iSpc++) {
            //                    var Spc = Config_SpeciesToCompute[iSpc];
            //                    int[] _OldSources_spc = m_OldSources[iSpc];
            //                    int[] _NewSources_spc = NewSources[iSpc];
            //                    Debug.Assert(_OldSources_spc.IsSubsetOf(_NewSources_spc));
            //                }
            //            }
            //            m_OldSources = NewSources;
            //            m_OldVolumes = NewVolumes;
            //#endif

            //foreach (var spId in m_CurrentAgglomeration.SpeciesList) {
            //    string SpName = m_LsTrk.GetSpeciesName(spId);
            //    int NoOfAgg = m_CurrentAgglomeration.GetAgglomerator(spId).AggInfo.AgglomerationPairs.Length;
            //    Console.WriteLine("Species {0}, time {2}, number of agglomerations: {1}", SpName, NoOfAgg, m_RequiredTimeLevels);
            //}
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="PrecondMassMatrix">
        /// No scaling, but with agglomeration -- used to compute the DG basis of the preconditioner
        /// </param>
        /// <param name="ScaledMassMatrix">
        /// No agglomeration, but with <see cref="XdgTimesteppingBase.Config_MassScale"/> applied.
        /// </param>
        /// <param name="time"></param>
        void UpdateMassMatrix(out BlockMsrMatrix PrecondMassMatrix, out BlockMsrMatrix ScaledMassMatrix, double time) {
            if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                // i.e. if the Level-Set is outside the domain.

                ScaledMassMatrix = null;
                PrecondMassMatrix = null;
            } else {
                // checks:
                Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                    || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent
                    || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent,
                    "Something is not implemented here.");
                if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                    && Config_LevelSetHandling != LevelSetHandling.None)
                    throw new NotSupportedException("Illegal configuration;");


                if ((this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity)
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent)
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent)
                    ) {
                    //MassMatrixFactory MassFact = new MassMatrixFactory(CurrentStateMapping.BasisS.ElementAtMax(b => b.Degree), m_CurrentAgglomeration);
                    MassMatrixFactory MassFact = m_LsTrk.GetXDGSpaceMetrics(Config_SpeciesToCompute, Config_CutCellQuadratureOrder, 1).MassMatrixFactory;
                    PrecondMassMatrix = MassFact.GetMassMatrix(CurrentStateMapping, false);
                    m_CurrentAgglomeration.ManipulateMatrixAndRHS(PrecondMassMatrix, default(double[]), CurrentStateMapping, CurrentStateMapping);
                    ScaledMassMatrix = new BlockMsrMatrix(CurrentStateMapping);

                    int NF = this.CurrentStateMapping.Fields.Count;
                    //MassFact.AccMassMatrix(ScaledMassMatrix, CurrentStateMapping, _alpha: Config_MassScale);
                    base.ComputeMassMatrixImpl(ScaledMassMatrix, time);
                } else {
                    throw new NotSupportedException();
                }
            }
        }

        bool OneTimeMgInit = false;

        /// <summary>
        /// Performs one timestep, on the DG fields in <see cref="XdgTimesteppingBase.CurrentStateMapping"/>.
        /// </summary>
        public void Solve(double phystime, double dt) {

            Debug.Assert(m_RKscheme.c.Length == m_RKscheme.a.GetLength(0));
            Debug.Assert(m_RKscheme.b.Length == m_RKscheme.a.GetLength(1));

            // update multigrid basis _once_ in object lifetime for steady level set:
            if (this.Config_LevelSetHandling == LevelSetHandling.None && OneTimeMgInit == false) {
                base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
                OneTimeMgInit = true;
            }



            // ===========================================
            // update level-set (in the case of splitting)
            // ===========================================
            bool performSplitting;
            if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {

                double ls_dt = dt;
                if (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting)
                    ls_dt *= 0.5;

                // remember which old cells had values
                //var oldCCM = this.UpdateCutCellMetrics();

                // evolve the level set

                m_LsTrk.PushStacks();

                int oldPushCount = m_LsTrk.PushCount;
                int oldVersion = m_LsTrk.VersionCnt;

                this.MoveLevelSetAndRelatedStuff(m_CurrentState.Mapping.Fields.ToArray(), phystime, ls_dt, 1.0, null, null);

                int newPushCount = m_LsTrk.PushCount;
                int newVersion = m_LsTrk.VersionCnt;
                if ((newPushCount - oldPushCount) != 0)
                    throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");
                if ((newVersion - oldVersion) != 1)
                    throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");


                // in the case of splitting, the fields must be extrapolated 
                //var newCCM = this.UpdateCutCellMetrics();
                //var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
                Debug.Assert(m_LsTrk.HistoryLength >= 1);
                var SplittingAgg = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                    __AgglomerationTreshold: 0.0, AgglomerateNewborn: true, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true,
                    oldTs__AgglomerationTreshold: new double[] { 0.0 });
                SplittingAgg.Extrapolate(this.CurrentStateMapping);

                // yes, we use splitting (i.e. only one mass matrix is required)
                performSplitting = true;
            } else {
                performSplitting = false;
            }

            // ==============================================
            // solve main system
            // ==============================================

            // init mass matrix & cut-cell metrics 
            BlockMsrMatrix[] MassMatrix = new BlockMsrMatrix[performSplitting ? 1 : 2];
            //m_AllCCM.Clear();
            m_RequiredTimeLevels = 0;
            m_LsTrk.PushStacks();
            this.UpdateAgglom(false);
            base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
            BlockMsrMatrix PM, SM;
            UpdateMassMatrix(out PM, out SM, phystime);
            MassMatrix[0] = SM;
            m_PrecondMassMatrix = PM;

            // initial value
            CoordinateVector u0 = new CoordinateVector(this.CurrentStateMapping.Fields.Select(f => f.CloneAs()).ToArray());
            foreach (var f in u0.Mapping.Fields) {
                if (f is XDGField) {
                    ((XDGField)f).UpdateBehaviour = BehaveUnder_LevSetMoovement.PreserveMemory;
                }
            }

            // loop over Runge-Kutta stages...
            double[][] k = new double[m_RKscheme.Stages][];
            for (int s = 0; s < m_RKscheme.Stages; s++) {
                RKstage(phystime, dt, k, s, MassMatrix, u0, s > 0 ? m_RKscheme.c[s - 1] : 0.0);
                k[s] = new double[this.CurrentStateMapping.LocalLength];
                UpdateChangeRate(phystime + dt * m_RKscheme.c[s], k[s]);
            }

            // final stage
            RKstageExplicit(phystime, dt, k, m_RKscheme.Stages, MassMatrix, u0, m_RKscheme.c[m_RKscheme.Stages - 1], m_RKscheme.b, 1.0);

            // ===========================================
            // update level-set (in the case of splitting)
            // ===========================================

            if (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                // remember which old cells had values
                //var oldCCM = this.UpdateCutCellMetrics();

                // evolve the level set
                m_LsTrk.PushStacks();

                int oldPushCount = m_LsTrk.PushCount;
                int oldVersion = m_LsTrk.VersionCnt;

                this.MoveLevelSetAndRelatedStuff(m_CurrentState.Mapping.Fields.ToArray(), phystime + dt * 0.5, dt * 0.5, 1.0, null, null);

                int newPushCount = m_LsTrk.PushCount;
                int newVersion = m_LsTrk.VersionCnt;
                if ((newPushCount - oldPushCount) != 0)
                    throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");
                if ((newVersion - oldVersion) != 1)
                    throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");

                // in the case of splitting, the fields must be extrapolated 
                //var newCCM = this.UpdateCutCellMetrics();
                //var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
                Debug.Assert(m_LsTrk.HistoryLength >= 1);
                var SplittingAgg = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                    __AgglomerationTreshold: 0.0, AgglomerateNewborn: true, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true,
                    oldTs__AgglomerationTreshold: new double[] { 0.0 });
                SplittingAgg.Extrapolate(this.CurrentStateMapping);
            }
        }


        private void RKstage(double PhysTime, double dt, double[][] k, int s, BlockMsrMatrix[] Mass, CoordinateVector u0,
            double ActualLevSetRelTime) {

            // detect whether the stage s is explicit or not: (implicit schemes can have some explicit stages too)
            bool isExplicit = m_RKscheme.a[s, s] == 0;
            double[] RK_as = m_RKscheme.a.GetRow(s);

            double RelTime = m_RKscheme.c[s];


            // =========================
            // Compute intermediate step
            // =========================


            if (isExplicit) {
                // +++++++++++++++++++++
                // Explicit stage branch
                // +++++++++++++++++++++

                RKstageExplicit(PhysTime, dt, k, s, Mass, u0, ActualLevSetRelTime, RK_as, RelTime);
            } else {
                // +++++++++++++++++++++
                // Implicit stage branch
                // +++++++++++++++++++++

                RKstageImplicit(PhysTime, dt, k, s, Mass, u0, ActualLevSetRelTime, RK_as, RelTime);
            }


        }

        private void RKstageImplicit(double PhysTime, double dt, double[][] k, int s, BlockMsrMatrix[] Mass, CoordinateVector u0, double ActualLevSetRelTime, double[] RK_as, double RelTime) {
            Debug.Assert(s < m_RKscheme.Stages);
            Debug.Assert(m_RKscheme.c[s] > 0);
            Debug.Assert(RK_as[s] != 0);

            int Ndof = m_CurrentState.Count;

            // =========
            // RHS setup
            // =========

            m_ImplStParams = new ImplicitStage_AssiParams() {
                m_CurrentDt = dt,
                m_CurrentPhystime = PhysTime,
                m_IterationCounter = 0,
                m_ActualLevSetRelTime = ActualLevSetRelTime,
                m_RelTime = RelTime,
                m_k = k,
                m_u0 = u0,
                m_Mass = Mass,
                m_RK_as = RK_as,
                m_s = s
            };


            // ================
            // solve the system
            // ================

            NonlinearSolver nonlinSolver;
            ISolverSmootherTemplate linearSolver;
            GetSolver(out nonlinSolver, out linearSolver);


            if (RequiresNonlinearSolver) {

                // Nonlinear Solver (Navier-Stokes)
                // --------------------------------

                nonlinSolver.SolverDriver(m_CurrentState, default(double[])); // Note: the RHS is passed as the affine part via 'this.SolverCallback'

            } else {
                // Linear Solver (Stokes)
                // ----------------------


                // build the saddle-point matrix
                BlockMsrMatrix System, MaMa;
                double[] RHS;
                this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true);
                RHS.ScaleV(-1);

                // update the multigrid operator
                MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                    System, MaMa,
                    this.Config_MultigridOperator);

                // init linear solver
                linearSolver.Init(mgOperator);

                // try to solve the saddle-point system.
                mgOperator.UseSolver(linearSolver, m_CurrentState, RHS);

                // 'revert' agglomeration
                m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);
            }

            // ================
            // reset
            // ================
            m_ImplStParams = null;
        }

        BlockMsrMatrix m_PrecondMassMatrix;

        ImplicitStage_AssiParams m_ImplStParams = null;

        /// <summary>
        /// Helper/global variables for implicit RK-Stages, used by <see cref="AssembleMatrixCallback(out BlockMsrMatrix, out double[], out BlockMsrMatrix, DGField[])"/>.
        /// </summary>
        class ImplicitStage_AssiParams {
            public double m_CurrentPhystime;
            public double m_CurrentDt = -1;
            public int m_IterationCounter = 0;
            public double m_ActualLevSetRelTime;
            public double m_RelTime;
            public double[][] m_k;
            public CoordinateVector m_u0;
            public BlockMsrMatrix[] m_Mass;
            public double[] m_RK_as;
            public int m_s;
        }

        /// <summary>
        /// Matrix/Affine assembly in the case of an implicit RK stage.
        /// </summary>
        protected override void AssembleMatrixCallback(out BlockMsrMatrix System, out double[] Affine, out BlockMsrMatrix PcMassMatrix, DGField[] argCurSt, bool Linearization) {

            if (Linearization == false)
                throw new NotImplementedException("todo");

            int Ndof = m_CurrentState.Count;

            // copy data from 'argCurSt' to 'CurrentStateMapping', if necessary 
            // -----------------------------------------------------------
            DGField[] locCurSt = CurrentStateMapping.Fields.ToArray();
            if (locCurSt.Length != argCurSt.Length) {
                throw new ApplicationException();
            }
            int NF = locCurSt.Length;
            for (int iF = 0; iF < NF; iF++) {
                if (object.ReferenceEquals(locCurSt[iF], argCurSt[iF])) {
                    // nothing to do
                } else {
                    locCurSt[iF].Clear();
                    locCurSt[iF].Acc(1.0, argCurSt[iF]);
                }
            }

            // update of level-set
            // ----------------------

            bool updateAgglom = false;
            if (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once && m_ImplStParams.m_IterationCounter == 0
                || this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {

                //MoveLevelSetAndRelatedStuff(locCurSt, m_CurrentPhystime, m_CurrentDt, 1.0);
                if (Math.Abs(m_ImplStParams.m_ActualLevSetRelTime - m_ImplStParams.m_RelTime) > 1.0e-14) {
                    if (m_ImplStParams.m_IterationCounter <= 0)// only push tracker in the first iter
                        m_LsTrk.PushStacks();
                    MoveLevelSetAndRelatedStuff(locCurSt,
                        m_ImplStParams.m_CurrentPhystime, m_ImplStParams.m_CurrentDt * m_ImplStParams.m_RelTime, IterUnderrelax,
                        m_ImplStParams.m_Mass, m_ImplStParams.m_k);

                    // note that we need to update the agglomeration
                    updateAgglom = true;
                }
            }

            // update agglomeration
            // --------------------
#if DEBUG
            if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                Debug.Assert(m_CurrentAgglomeration != null);
                Debug.Assert(m_PrecondMassMatrix != null);
                // ensure, that, when splitting is used we update the agglomerator in the very first iteration.
            }
#endif
            if (updateAgglom || m_CurrentAgglomeration == null) {
                this.UpdateAgglom(m_ImplStParams.m_IterationCounter > 0);



                // update Multigrid-XDG basis
                base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
            }


            // mass matrix update
            // ------------------

            if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                // i.e. if the Level-Set is outside the domain.

                PcMassMatrix = null;

            } else {
                // checks:
                Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                    || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent
                    || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent,
                    "Something is not implemented here.");
                if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                    && Config_LevelSetHandling != LevelSetHandling.None)
                    throw new NotSupportedException("Illegal configuration;");

                if (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once || this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {

                    if ((this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity && m_PrecondMassMatrix == null) // compute mass matrix (only once in application lifetime)
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent && m_ImplStParams.m_IterationCounter == 0) // compute mass matrix once per timestep
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent) // re-compute mass matrix in every iteration
                    ) {

                        BlockMsrMatrix PM, SM;
                        UpdateMassMatrix(out PM, out SM, m_ImplStParams.m_CurrentPhystime + m_ImplStParams.m_CurrentDt * m_ImplStParams.m_RelTime);
                        m_ImplStParams.m_Mass[1] = SM;
                        m_PrecondMassMatrix = PM;
                    }
                }

                PcMassMatrix = m_PrecondMassMatrix;
            }

            // operator matrix update
            // ----------------------

            // we perform the extrapolation to have valid parameters if
            // - the operator matrix depends on these values
            this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

            var OpMatrix = new BlockMsrMatrix(CurrentStateMapping);
            var OpAffine = new double[Ndof];

            this.ComputeOperatorMatrix(OpMatrix, OpAffine, CurrentStateMapping, locCurSt, base.GetAgglomeratedLengthScales(), m_ImplStParams.m_CurrentPhystime + m_ImplStParams.m_CurrentDt * m_ImplStParams.m_RelTime);


            // assemble system
            // ---------------

            double dt = m_ImplStParams.m_CurrentDt;

            // select mass matrix (and some checks)
            double[] RHS = new double[Ndof];
            BlockMsrMatrix MamaRHS, MamaLHS;
            if (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once
                || this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {

                Debug.Assert(m_ImplStParams.m_Mass.Length == 2);
                MamaRHS = m_ImplStParams.m_Mass[0];
                MamaLHS = m_ImplStParams.m_Mass[1];

            } else if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                || this.Config_LevelSetHandling == LevelSetHandling.None) {

                Debug.Assert(m_ImplStParams.m_Mass.Length == 1);
                MamaRHS = m_ImplStParams.m_Mass[0];
                MamaLHS = m_ImplStParams.m_Mass[0];

            } else {
                throw new NotImplementedException();
            }

            // right-hand-side, resp. affine vector
            if (MamaRHS != null) {
                MamaRHS.SpMV(1.0 / dt, m_ImplStParams.m_u0, 0.0, RHS);
            } else {
                Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                RHS.SetV(m_ImplStParams.m_u0, 1.0 / dt);
            }
            for (int l = 0; l < m_ImplStParams.m_s; l++) {
                if (m_ImplStParams.m_RK_as[l] != 0.0)
                    RHS.AccV(-m_ImplStParams.m_RK_as[l], m_ImplStParams.m_k[l]);
            }

            Affine = RHS;
            Affine.ScaleV(-1.0);
            Affine.AccV(m_ImplStParams.m_RK_as[m_ImplStParams.m_s], OpAffine);

            // left-hand-side
            System = OpMatrix.CloneAs();
            System.Scale(m_ImplStParams.m_RK_as[m_ImplStParams.m_s]);
            if (MamaLHS != null) {
                System.Acc(1.0 / dt, MamaLHS);
            } else {
                System.AccEyeSp(1.0 / dt);
            }

            // perform agglomeration
            // ---------------------
            Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
            m_CurrentAgglomeration.ManipulateMatrixAndRHS(System, Affine, CurrentStateMapping, CurrentStateMapping);

            // increase iteration counter         
            // --------------------------
            m_ImplStParams.m_IterationCounter++;
        }


        //protected override void TriggerLevelSetUpdate() {
        //    throw new NotImplementedException();
        //}



        private void RKstageExplicit(double PhysTime, double dt, double[][] k, int s, BlockMsrMatrix[] Mass, CoordinateVector u0, double ActualLevSetRelTime, double[] RK_as, double RelTime) {
            Debug.Assert(s <= m_RKscheme.Stages);
            for (int i = 0; i < s; i++) {
                Debug.Assert(k[i] != null);
            }


            int Ndof = CurrentStateMapping.LocalLength;
            BlockMsrMatrix System;
            double[] RHS = new double[Ndof];


            if (RelTime > 0) {
                //
                // some ordinary stage above stage 0
                //

                if (base.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative
                     || base.Config_LevelSetHandling == LevelSetHandling.Coupled_Once) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // Moving interface, RK stage 'i':
                    // (c[s]/dt)*(Ms*us - M0*u0) + sum(k[l]*a[s,l], 0 <= l < s)) = 0
                    // For explicit timestepping, both versions of level-set coupling 
                    // will lead to identical results.
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    // move level-set:
                    if (Math.Abs(ActualLevSetRelTime - RelTime) > 1.0e-14) {

                        this.m_LsTrk.PushStacks();
                        this.MoveLevelSetAndRelatedStuff(u0.Mapping.Fields.ToArray(), PhysTime, dt * RelTime, IterUnderrelax, Mass, k);

                        this.UpdateAgglom(false);
                        BlockMsrMatrix PM, SM;
                        UpdateMassMatrix(out PM, out SM, PhysTime + dt * RelTime);
                        Mass[1] = SM;
                    }

                    var M0 = Mass[0];
                    var Ms = Mass[1];

                    // left-hand-side
                    if (Ms != null) {
                        System = Ms.CloneAs();
                        System.Scale(1.0 / dt);
                    } else {
                        Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                        System = null;
                    }

                    // right-hand-side
                    if (M0 != null) {
                        M0.SpMV(1.0 / dt, u0, 0.0, RHS);
                    } else {
                        Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                        RHS.SetV(u0, 1.0 / dt);
                    }
                    for (int l = 0; l < s; l++) {
                        if (RK_as[l] != 0.0)
                            RHS.AccV(-RK_as[l], k[l]);
                    }

                } else if (base.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                     || base.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                     || base.Config_LevelSetHandling == LevelSetHandling.None) {

                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // Level-Set-Splitting: Mass matrix remains constant during Runge-Kutta
                    // Stage:
                    // (c[s]/dt)*M0*(us - u0) + sum(k[l]*a[s,l], l=0..(s-1)) = 0
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    Debug.Assert(Mass.Length == 1);
                    var Ms = Mass[0];

                    // left-hand-side
                    if (Ms != null) {
                        System = Ms.CloneAs();
                        System.Scale(1.0 / dt);
                    } else {
                        Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                        System = null;
                    }

                    // right-hand-side
                    if (Ms != null) {
                        Ms.SpMV(1.0 / dt, u0, 0.0, RHS);
                    } else {
                        Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                        RHS.SetV(u0, 1.0 / dt);
                    }
                    for (int l = 0; l < s; l++) {
                        if (RK_as[l] != 0.0)
                            RHS.AccV(-RK_as[l], k[l]);
                    }

                } else {
                    throw new NotImplementedException();
                }

                // solve system
                if (System != null) {
                    Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                    m_CurrentAgglomeration.ManipulateMatrixAndRHS(System, RHS, this.CurrentStateMapping, this.CurrentStateMapping);
                    BlockSol(System, m_CurrentState, RHS);
                    m_CurrentAgglomeration.Extrapolate(this.CurrentStateMapping);
                } else {
                    // system is diagonal: 1/dt
                    m_CurrentState.SetV(RHS, dt);
                }

            } else {
                // +++++++++++++++++++++++++++
                // this is probably stage zero
                // +++++++++++++++++++++++++++

                m_CurrentState.Clear();
                m_CurrentState.Acc(1.0, u0);
            }
        }



        private void UpdateChangeRate(double PhysTime, double[] k) {

            BlockMsrMatrix OpMtx = new BlockMsrMatrix(this.CurrentStateMapping);
            double[] OpAff = new double[this.CurrentStateMapping.LocalLength];
            base.ComputeOperatorMatrix(OpMtx, OpAff, this.CurrentStateMapping, this.CurrentStateMapping.Fields.ToArray(), base.GetAgglomeratedLengthScales(), PhysTime);

            k.SetV(OpAff);
            OpMtx.SpMV(1.0, this.m_CurrentState, 1.0, k);
        }

        void BlockSol<V1, V2>(BlockMsrMatrix M, V1 X, V2 B)
            where V1 : IList<double>
            where V2 : IList<double> //
        {
            int i0 = M.RowPartitioning.i0;
            int iE = M.RowPartitioning.iE;

            var Part = M.RowPartitioning;
            Debug.Assert(Part.EqualsPartition(this.CurrentStateMapping));

            int J = m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells;

            double[] MtxVals = null;
            int[] Indices = null;

            MultidimensionalArray Block = null;
            double[] x = null, b = null;
            for (int j = 0; j < J; j++) {
                int bS = this.CurrentStateMapping.LocalUniqueCoordinateIndex(0, j, 0);
                int Nj = this.CurrentStateMapping.GetTotalNoOfCoordinatesPerCell(j);

                if (Block == null || Block.NoOfRows != Nj) {
                    Block = MultidimensionalArray.Create(Nj, Nj);
                    x = new double[Nj];
                    b = new double[Nj];
                } else {
                    Block.Clear();
                }


                // extract block and part of RHS
                for (int iRow = 0; iRow < Nj; iRow++) {
                    bool ZeroRow = true;
                    //MsrMatrix.MatrixEntry[] row = M.GetRow(iRow + bS + i0);
                    int LR = M.GetRow(iRow + bS + i0, ref Indices, ref MtxVals);

                    //foreach (var entry in row) {
                    for (int lr = 0; lr < LR; lr++) {
                        int ColIndex = Indices[lr];
                        double Value = MtxVals[lr];

                        Block[iRow, ColIndex - (bS + i0)] = Value;
                        if (Value != 0.0)
                            ZeroRow = false;
                    }
                    b[iRow] = B[iRow + bS];

                    if (ZeroRow) {
                        if (b[iRow] != 0.0)
                            throw new ArithmeticException();
                        else
                            Block[iRow, iRow] = 1.0;
                    }
                }

                // solve
                Block.SolveSymmetric(x, b);

                // store solution
                for (int iRow = 0; iRow < Nj; iRow++) {
                    X[iRow + bS] = x[iRow];
                }
            }
        }

    }
}
