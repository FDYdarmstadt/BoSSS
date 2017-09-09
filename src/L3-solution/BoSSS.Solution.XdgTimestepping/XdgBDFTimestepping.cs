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
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform;
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.Timestepping;
using ilPSP;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using MPI.Wrappers;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// switch for the initialization of the <see cref="XdgBDFTimestepping"/> 
    /// </summary>
    public enum TimestepperInit {

        /// Initialization from a single timestep, i.e. if this time-stepper should use BDF4,
        /// it starts with BDF1, BDF2, BDF3 in the first, second and third time-step.
        SingleInit,

        /// same initialization for SingleInit, but the first timesteps 
        /// are computed with a smaller timestepsize
        IncrementInit,

        /// Initialization for a multi-step method, e.g. BDF4 requires 4 timesteps.
        /// can be used if an analytic solution is known or simulation is restarted form another session
        MultiInit
    }

    /// <summary>
    /// Implicit time-stepping using Backward-Differentiation-Formulas (BDF),
    /// specialized for XDG applications.
    /// </summary>
    public class XdgBDFTimestepping : XdgTimesteppingBase {

        /// <summary>
        /// Constructor;
        /// </summary>
        /// <param name="Fields"></param>
        /// <param name="IterationResiduals"></param>
        /// <param name="LsTrk"></param>
        /// <param name="_ComputeOperatorMatrix">See <see cref="ComputeOperatorMatrix"/>.</param>
        /// <param name="_UpdateLevelset">See <see cref="UpdateLevelset"/>.</param>
        /// <param name="_UpdateCutCellMetrics">See <see cref="UpdateCutCellMetrics"/>.</param>
        /// <param name="BDForder">
        /// The order of the BDF scheme from 1 to 6; in addition, 0 encodes Explicit Euler and -1 encodes Crank-Nicolson.
        /// </param>
        /// <param name="_LevelSetHandling"></param>
        /// <param name="_MassMatrixShapeandDependence"></param>
        /// <param name="_SpatialOperatorType"></param>
        /// <param name="_MassScale"></param>
        /// <param name="_AgglomerationThreshold"></param>
        /// <param name="DelayInit">
        /// Triggers a delayed initialization of the BDF scheme.
        /// If true, it is the users responsibility to call either <see cref="SingleInit"/> or <see cref="MultiInit(double, double, Action{int, double, DGField[]})"/>.
        /// </param>
        /// <param name="useX">
        /// Use XDG-Fields because of Level Set. -> DIRTY HACK
        /// </param>
        /// <param name="_MultigridOperatorConfig">
        /// Configuration of block-preconditioner, if null a default value is chosen.
        /// </param>
        public XdgBDFTimestepping(IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            LevelSetTracker LsTrk,
            bool DelayInit,
            DelComputeOperatorMatrix _ComputeOperatorMatrix,
            DelUpdateLevelset _UpdateLevelset,
            DelUpdateCutCellMetrics _UpdateCutCellMetrics,
            int BDForder,
            LevelSetHandling _LevelSetHandling,
            MassMatrixShapeandDependence _MassMatrixShapeandDependence,
            SpatialOperatorType _SpatialOperatorType,
            IDictionary<SpeciesId, IEnumerable<double>> _MassScale,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
            AggregationGrid[] _MultigridSequence,
            double _AgglomerationThreshold, bool _useX) {

            if (Fields.Count() != IterationResiduals.Count())
                throw new ArgumentException("Expecting the same number of fields and residuals.");
            for (int iFld = 0; iFld < Fields.Count(); iFld++) {
                if (!Fields.ElementAt(iFld).Basis.Equals(IterationResiduals.ElementAt(iFld).Basis))
                    throw new ArgumentException(string.Format("Mismatch between {0}-th basis of fields and residuals.", iFld));
            }
                       

            this.Config_LevelSetHandling = _LevelSetHandling;
            this.Config_MassMatrixShapeandDependence = _MassMatrixShapeandDependence;
            this.Config_SpatialOperatorType = _SpatialOperatorType;
            this.ComputeOperatorMatrix = _ComputeOperatorMatrix;
            this.UpdateLevelset = _UpdateLevelset;
            this.UpdateCutCellMetrics = _UpdateCutCellMetrics;
            this.Config_MassScale = _MassScale;
            this.Config_AgglomerationThreshold = _AgglomerationThreshold;
            this.useX = _useX;
            base.MultigridSequence = _MultigridSequence;
            if (_MultigridSequence == null || _MultigridSequence.Length < 1)
                throw new ArgumentException("At least one grid level is required.");

            base.Residuals = new CoordinateVector(IterationResiduals.ToArray());
            
            if (_MultigridOperatorConfig != null) {
                base.Config_MultigridOperator = _MultigridOperatorConfig;
            } else {
                SetConfig_MultigridOperator_Default(Fields);
            }

            base.CommonConfigurationChecks();

            switch (BDForder) {
                case -1:
                m_TSCchain = new BDFSchemeCoeffs[] { BDFSchemeCoeffs.CrankNicolson() };
                break;

                case 0:
                m_TSCchain = new BDFSchemeCoeffs[] { BDFSchemeCoeffs.ExplicitEuler() };
                break;

                default:
                m_TSCchain = new BDFSchemeCoeffs[BDForder];
                for (int i = BDForder; i >= 1; i--) {
                    m_TSCchain[BDForder - i] = BDFSchemeCoeffs.BDF(i);
                }
                break;
            }
            m_LsTrk = LsTrk;

            int S = m_TSCchain[0].S;
            Debug.Assert(S == m_TSCchain.Length);

            // cut-cell-metrics stack
            // ----------------------

            if (Config_LevelSetHandling == LevelSetHandling.None
                || Config_LevelSetHandling == LevelSetHandling.LieSplitting
                || Config_LevelSetHandling == LevelSetHandling.StrangSplitting)
                m_Stack_CutCellMetrics = new CutCellMetrics[1];
            else
                m_Stack_CutCellMetrics = new CutCellMetrics[S + 1];


            // mass-matrix stack
            // -----------------

            if (Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                m_Stack_MassMatrix = new BlockMsrMatrix[1];
            } else if (Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity) {
                m_Stack_MassMatrix = new BlockMsrMatrix[1];
            } else {
                if (Config_LevelSetHandling == LevelSetHandling.LieSplitting
                    || Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                    m_Stack_MassMatrix = new BlockMsrMatrix[1];
                } else {
                    m_Stack_MassMatrix = new BlockMsrMatrix[S + 1];
                }
            }

            // operator matrix - stack
            // -----------------------

            m_Stack_OpMatrix = new BlockMsrMatrix[m_TSCchain[0].theta0 != 0.0 ? 2 : 1]; // only required for Crank.Nic. and Expl. Euler,
            m_Stack_OpAffine = new double[m_Stack_OpMatrix.Length][]; //      in this case 'theta0' is unequal 0.0.

            // m_Stack_u
            // ---------

            m_Stack_u = new CoordinateVector[S + 1];

            m_Stack_u[0] = new CoordinateVector(Fields.ToArray());
            for (int s = 1; s <= S; s++) {
                m_Stack_u[s] = new CoordinateVector(Fields.Select(dgf => dgf.CloneAs()).ToArray());
            }
            for (int s = 0; s <= S; s++) {
                DGField[] DgFs = m_Stack_u[s].Mapping.Fields.ToArray();
                foreach (var Dgf in DgFs) {
                    if (Dgf is XDGField) {
                        ((XDGField)Dgf).UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
                    } else {
                        // should not really matter, 
                        // the agglomeration of newborn cells should take care of it!
                    }
                }
            }

            // multigrid - init
            // ----------------

            InitMultigrid(Fields.ToArray(), useX);


            // other stuff
            // -----------

            if (!DelayInit)
                InitTimestepping(true);
        }

        BDFSchemeCoeffs[] m_TSCchain;

        /// <summary>
        /// DG coefficient mapping for the test- and trial-space.
        /// </summary>
        protected override CoordinateMapping CurrentStateMapping {
            get {
                return m_Stack_u[0].Mapping;
            }
        }


        /// <summary>
        /// Switch for using XDG-Fields
        /// </summary>
        bool useX;

        //
        // stack of cut-cell metrics
        //
        CutCellMetrics[] m_Stack_CutCellMetrics;

        //
        // stack of mass matrices (matrices _without_ agglomeration)
        //
        BlockMsrMatrix[] m_Stack_MassMatrix;

        //
        // stack of operator matrices and affine vectors (matrices _without_ agglomeration)
        // (since we also do Crank-Nicolson, we may need one previous operator matrix)
        //
        BlockMsrMatrix[] m_Stack_OpMatrix;
        double[][] m_Stack_OpAffine;

        //
        // stack of solution vectors
        //
        CoordinateVector[] m_Stack_u;


        int m_PopulatedStackDepth = 0;

        /// <summary>
        /// Unscaled, agglomerated mass matrix used by the preconditioner.
        /// </summary>
        BlockMsrMatrix m_PrecondMassMatrix;


        /// <summary>
        /// 
        /// </summary>
        void PushStack(int increment) {

            m_PopulatedStackDepth++;
            if (m_PopulatedStackDepth > m_TSCchain[0].S)
                m_PopulatedStackDepth = m_TSCchain[0].S;

            // operator stack
            // --------------
            
            Debug.Assert(m_Stack_OpMatrix.Length <= 2);
            Debug.Assert(m_Stack_OpAffine.Length == m_Stack_OpMatrix.Length);
            if (m_Stack_OpMatrix.Length == 2) {
                var Mtmp = m_Stack_OpMatrix[0];
                var Atmp = m_Stack_OpAffine[0];
                m_Stack_OpMatrix[0] = m_Stack_OpMatrix[1];
                m_Stack_OpAffine[0] = m_Stack_OpAffine[1];
                m_Stack_OpMatrix[1] = Mtmp;
                m_Stack_OpAffine[1] = Atmp;
            }
            if (m_Stack_OpMatrix[0] != null)
                m_Stack_OpMatrix[0] = null;
            if (m_Stack_OpAffine[0] != null)
                m_Stack_OpAffine[0].ClearEntries();

            // Solution-Stack
            // --------------
            // entry 0 should remain the same object all the time
            var Cvtmp = m_Stack_u[m_Stack_u.Length - 1];
            for (int i = m_Stack_u.Length - 1; i >= 2; i--) {
                m_Stack_u[i] = m_Stack_u[i - 1];
            }
            m_Stack_u[1] = Cvtmp;
            m_Stack_u[1].Clear();
            m_Stack_u[1].Acc(1.0, m_Stack_u[0]);

            // mass-matrix stack
            // -----------------


            if (this.Config_MassMatrixShapeandDependence != MassMatrixShapeandDependence.IsNonIdentity) {
                for (int i = m_Stack_MassMatrix.Length - 1; i >= 1; i--) {
                    m_Stack_MassMatrix[i] = m_Stack_MassMatrix[i - 1];
                }
                m_Stack_MassMatrix[0] = null;
                m_PrecondMassMatrix = null;
            } else {
                // should already be initialized - see 'InitTimestepping'
                Debug.Assert(m_Stack_MassMatrix.Length == 1);
                Debug.Assert(m_Stack_MassMatrix[0] != null);
                Debug.Assert(m_PrecondMassMatrix != null);
            }

            // cut-cell metrics
            // ----------------
            if (this.Config_LevelSetHandling != LevelSetHandling.None) {
                // we expect the level-set to change in every timestep

                for (int i = m_Stack_CutCellMetrics.Length - 1; i >= 1; i--) {
                    m_Stack_CutCellMetrics[i] = m_Stack_CutCellMetrics[i - 1];
                }
                m_Stack_CutCellMetrics[0] = null;
            } else {
                // a level-set which is static over all timeteps

                Debug.Assert(m_Stack_CutCellMetrics[0] != null);
            }


            // saves first timesteps (actual timerstpe size dt) in  case of incrementInit
            if (incrementTimesteps > 1 && increment == 1 && m_CurrentPhystime > 0.0)
                PushIncrementStack();
        }

        /// <summary>
        /// saves the first necessary timesteps at the actual time level for a incremental initialization for higher BDF schemes 
        /// </summary>
        void PushIncrementStack() {

            incrementHist++;
            if (incrementHist == m_TSCchain[0].S - 1) {

                // Copy increment history to the actual stack with timestepsize dt
                // ---------------------------------------------------------------
                for (int i = m_TSCchain[0].S; i > 1; i--) {
                    // Solution-Stack
                    m_Stack_u[i].Clear();
                    m_Stack_u[i].Acc(1.0, m_Stack_u_incHist[m_TSCchain[0].S - i]);

                    // mass-matrix stack
                    if (m_Stack_MassMatrix_incHist != null)
                        m_Stack_MassMatrix[i] = m_Stack_MassMatrix_incHist[m_TSCchain[0].S - i];

                    // cut-cell metrics
                    if (m_Stack_CutCellMetrics_incHist != null)
                        m_Stack_CutCellMetrics[i] = m_Stack_CutCellMetrics_incHist[m_TSCchain[0].S - i];

                    // get rid off the incrementHistory
                    m_Stack_u_incHist = null;
                    if (m_Stack_MassMatrix_incHist != null)
                        m_Stack_MassMatrix_incHist = null;
                    if (m_Stack_CutCellMetrics_incHist != null)
                        m_Stack_CutCellMetrics_incHist = null;

                }

                // restore the actual timestepsize
                m_CurrentDt *= incrementTimesteps;
                incrementTimesteps = 1;

                return;
            }

            // save history for the actual timestep size

            // Solution-Stack
            m_Stack_u_incHist[incrementHist].Clear();
            m_Stack_u_incHist[incrementHist].Acc(1.0, m_Stack_u[0]);

            // mass-matrix stack
            if (m_Stack_MassMatrix_incHist != null) 
                m_Stack_MassMatrix_incHist[incrementHist] = m_Stack_MassMatrix[1];

            // cut-cell metrics
            if (m_Stack_CutCellMetrics_incHist != null)
                m_Stack_CutCellMetrics_incHist[incrementHist] = m_Stack_CutCellMetrics[1];

        }


        int m_IterationCounter = 0;

        /// <summary>
        /// Initialization from a single timestep, i.e. if this time-stepper should use BDF4,
        /// it starts with BDF1, BDF2, BDF3 in the first, second and third time-step.
        /// </summary>
        /// <remarks>
        /// This approach in general does not provide the desired convergence order, but is applicable 
        /// is the simulation does not depend on the initial value.
        /// </remarks>
        public void SingleInit() {

            using (new FuncTrace()) { }

                InitTimestepping(true);

                if (Timestepper_Init == TimestepperInit.IncrementInit) {
                    if (incrementTimesteps <= 1)
                        throw new ArgumentOutOfRangeException("incrementInit needs a number of increment timesteps larger than 1");

                    InitIncrementStack();
                }

        }

        /// <summary>
        /// Initialization for a multi-step method, e.g. BDF4 requires 4 timesteps.
        /// </summary>
        /// <param name="physTime"></param>
        /// <param name="dt"></param>
        /// <param name="SetTimestep"></param>
        public void MultiInit(double physTime, int TimestepNo, double dt, Action<int, double, DGField[]> SetTimestep) {
            using (new FuncTrace()) {
                if (dt <= 0)
                    throw new ArgumentOutOfRangeException();
                if (m_CurrentDt > 0 && Math.Abs(dt / m_CurrentDt - 1.0) > 1.0e-14)
                    throw new ArgumentOutOfRangeException();

                int S = m_TSCchain[0].S;
                Debug.Assert(S == m_TSCchain.Length);

                for (int iStage = 0; iStage < S; iStage++) {
                    int TimeIndex = -S + iStage + 1;
                    double time = physTime + TimeIndex * dt;

                    {
                        int oldVersion = m_LsTrk.VersionCnt;
                        SetTimestep(TimeIndex + TimestepNo, time, m_Stack_u[0].Mapping.Fields.ToArray()); // the push-operation will init the other steps.
                        int newVersion = m_LsTrk.VersionCnt;

                        if ((newVersion - oldVersion) != 1)
                            throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");
                    }

                    // re-sort mass matrices
                    {
                        var Mtx2Update = m_Stack_MassMatrix.GetSubVector(1, m_Stack_MassMatrix.Length - 1);
                        TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, Mtx2Update, CurrentStateMapping, CurrentStateMapping);
                        Array.Copy(Mtx2Update, 0, m_Stack_MassMatrix, 1, Mtx2Update.Length);
                    }

                    // all the other stuff (cut-cell-metrics, ...)
                    InitTimestepping(iStage == (S - 1));

                    if (iStage < (S - 1))
                        PushStack(TimestepNo);
                }
            }
        }


        /// <summary>
        /// final initialization for the BDF timestepper scheme, all necessary timesteps have to be initialized
        /// </summary>
        /// <param name="OpInit"></param>
        private void InitTimestepping(bool OpInit) {

            {
                int[] Jtot =
                    (new int[] { base.m_LsTrk._Regions.GetCutCellMask().NoOfItemsLocally.MPISum(), base.m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells })
                    .MPISum();
                //Console.WriteLine("No of cells {0}, No of cut cells {1}.", Jtot[1], Jtot[0]);
                if (Jtot[0] == Jtot[1])
                    throw new ArithmeticException("All cells are cut cells - check your settings!");
            }


            m_Stack_CutCellMetrics[0] = this.UpdateCutCellMetrics();
            if (!m_Stack_CutCellMetrics[0].SpeciesList.SetEquals(Config_MassScale.Keys))
                throw new ApplicationException("Mismatch between species lists.");

            m_CurrentAgglomeration = new MultiphaseCellAgglomerator(
                    m_Stack_CutCellMetrics[0],
                    this.Config_AgglomerationThreshold, false, false, true,
                    null, null);

            // update Multigrid-XDG basis
            Debug.Assert(object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
            base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);

            // mass matrix update
            // ------------------

            if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                // i.e. if the Level-Set is outside the domain.

                foreach(var kv in this.Config_MassScale) {
                    SpeciesId spId = kv.Key;
                    double[] scaleVec = kv.Value.ToArray();
                    for(int i = 0; i < scaleVec.Length; i++) {
                        //if (scaleVec[i] != 1.0)
                        //    throw new ArithmeticException(string.Format("XDG time-stepping: illegal mass scale, mass matrix option {0} is set, but scaling factor for species {1}, variable no. {2} ({3}) is set to {4} (expecting 1.0).", 
                        //        MassMatrixShapeandDependence.IsIdentity, this.m_LsTrk.GetSpeciesName(kv.Key), i, this.CurrentStateMapping.Fields[i].Identification, scaleVec[i]));
                    }
                }

            } else {
                // checks:
                Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                    || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent
                    || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent,
                    "Something is not implemented here.");
                if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                    && Config_LevelSetHandling != LevelSetHandling.None)
                    throw new NotSupportedException("Illegal configuration;");


                // compute mass matrix (only once in application lifetime)
                Debug.Assert((m_PrecondMassMatrix == null) == (m_Stack_MassMatrix[0] == null));
                if ((this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity && m_PrecondMassMatrix == null)
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent && m_IterationCounter == 0)
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent)
                    ) {
                    int NF = CurrentStateMapping.Fields.Count;
                    MassMatrixFactory MassFact = new MassMatrixFactory(CurrentStateMapping.BasisS.ElementAtMax(b => b.Degree), m_CurrentAgglomeration);
                    m_PrecondMassMatrix = MassFact.GetMassMatrix(CurrentStateMapping, false);
                    m_Stack_MassMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                    MassFact.AccMassMatrix(m_Stack_MassMatrix[0], CurrentStateMapping, _alpha: Config_MassScale, VariableAgglomerationSwitch: NF.ForLoop(i => false));
                }
            }


            // operator matrix update
            // ----------------------

            if (OpInit && m_TSCchain[0].theta0 != 0.0) {
                // we perform the extrapolation to have valid parameters if
                // - the operator matrix depends on these values
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

                if (m_Stack_OpMatrix[0] == null) {
                    m_Stack_OpMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                }
                if (m_Stack_OpAffine[0] == null) {
                    m_Stack_OpAffine[0] = new double[CurrentStateMapping.LocalLength];
                }

                Debug.Assert(m_Stack_OpMatrix[0].InfNorm() == 0);
                Debug.Assert(m_Stack_OpAffine[0].L2Norm() == 0);
                this.ComputeOperatorMatrix(m_Stack_OpMatrix[0], m_Stack_OpAffine[0], CurrentStateMapping, CurrentStateMapping.Fields.ToArray(), m_CurrentAgglomeration, m_CurrentPhystime + m_CurrentDt);
            }
        }


        private string GetName__Stack_u(int i, int iF) {
            return this.GetType().FullName + "::Stack_u[" + i + "," + iF + "]";
        }

        private string GetName__Stack_MassMatrix(int i) {
            return this.GetType().FullName + "::Stack_MassMatrix[" + i + "]";
        }

        private string GetName__Stack_OpAffine(int i) {
            return this.GetType().FullName + "::Stack_OpAffine[" + i + "]";
        }

        private string GetName__Stack_OpMatrix(int i) {
            return this.GetType().FullName + "::Stack_OpMatrix[" + i + "]";
        }

        private string GetName__CutCellMetrics(int i) {
            return this.GetType().FullName + "::CutCellMetrics[" + i + "]";
        }

        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public void DataBackupBeforeBalancing(LoadBalancingData L) {
            using (new FuncTrace()) {
                if (m_PrivateBalancingInfo != null)
                    throw new NotSupportedException();
                m_PrivateBalancingInfo = new PrivateBalancingInfo();
                m_PrivateBalancingInfo.NoOfFields = m_Stack_u[0].Mapping.Fields.Count;

                // backup solution stack
                m_PrivateBalancingInfo.m_Stack_u = new bool[m_Stack_u.Length];
                for (int i = 0; i < m_Stack_u.Length; i++) {
                    CoordinateVector U_i = m_Stack_u[i];
                    if (U_i != null) {
                        m_PrivateBalancingInfo.m_Stack_u[i] = true;

                        DGField[] Fields = U_i.Mapping.Fields.ToArray();
                        for (int iF = 0; iF < Fields.Length; iF++) {
                            L.BackupField(Fields[iF], GetName__Stack_u(i, iF));
                        }
                    }
                }
                var map = CurrentStateMapping;
                m_Stack_u = null;

                // backup mass matrices
                m_PrivateBalancingInfo.m_Stack_MassMatrix = new bool[m_Stack_MassMatrix.Length];
                for (int i = 0; i < m_Stack_MassMatrix.Length; i++) {
                    if (m_Stack_MassMatrix[i] != null) {
                        m_PrivateBalancingInfo.m_Stack_MassMatrix[i] = true;
                        L.BackupMatrix(m_Stack_MassMatrix[i], GetName__Stack_MassMatrix(i), map, map);
                    }
                }
                m_Stack_MassMatrix = null;

                // backup operator
                Debug.Assert(m_Stack_OpMatrix.Length == m_Stack_OpAffine.Length);
                m_PrivateBalancingInfo.m_Stack_OpMatrix = new bool[m_Stack_OpMatrix.Length];
                m_PrivateBalancingInfo.m_Stack_OpAffine = new bool[m_Stack_OpAffine.Length];
                for (int i = 0; i < m_Stack_OpMatrix.Length; i++) {
                    if (m_Stack_OpMatrix[i] != null) {
                        m_PrivateBalancingInfo.m_Stack_OpMatrix[i] = true;
                        L.BackupMatrix(m_Stack_OpMatrix[i], GetName__Stack_OpMatrix(i), map, map);
                    }

                    if (m_Stack_OpAffine[i] != null) {
                        m_PrivateBalancingInfo.m_Stack_OpAffine[i] = true;
                        L.BackupVector(m_Stack_OpAffine[i], GetName__Stack_OpAffine(i));
                    }
                }
                m_Stack_OpMatrix = null;
                m_Stack_OpAffine = null;

                // backup the cut-cell metrics
                m_PrivateBalancingInfo.m_Stack_CutCellMetrics = new bool[m_Stack_CutCellMetrics.Length];
                for (int i = 0; i < m_Stack_CutCellMetrics.Length; i++) {
                    if (m_Stack_CutCellMetrics[i] != null) {
                        m_PrivateBalancingInfo.m_Stack_CutCellMetrics[i] = true;
                        L.BackupCutCellMetrics(m_Stack_CutCellMetrics[i], GetName__CutCellMetrics(i));
                    }
                }
                m_Stack_CutCellMetrics = null;

                // Delete agglomeration
                m_CurrentAgglomeration = null;
                base.MultigridBasis = null;
                base.MultigridSequence = null;
            }
        }

        PrivateBalancingInfo m_PrivateBalancingInfo;

        class PrivateBalancingInfo {
            public int NoOfFields; 
            public bool[] m_Stack_u;
            public bool[] m_Stack_MassMatrix;
            public bool[] m_Stack_OpMatrix;
            public bool[] m_Stack_OpAffine;
            public bool[] m_Stack_CutCellMetrics;
        }



        /// <summary>
        /// Step 2 of 2 for dynamic load balancing: restore this objects 
        /// status after the grid has been re-distributed.
        /// </summary>
        public void DataRestoreAfterBalancing(LoadBalancingData L,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            LevelSetTracker LsTrk,
            AggregationGrid[] _MultigridSequence) //
        {
            using (new FuncTrace()) {
                if (m_PrivateBalancingInfo == null)
                    throw new NotSupportedException();

                base.m_LsTrk = LsTrk;

                if (Fields.Count() != m_PrivateBalancingInfo.NoOfFields)
                    throw new ArgumentException();
                if (IterationResiduals.Count() != m_PrivateBalancingInfo.NoOfFields)
                    throw new ArgumentException();

                // restore solution stack
                m_Stack_u = new CoordinateVector[m_PrivateBalancingInfo.m_Stack_u.Length];
                m_Stack_u[0] = new CoordinateVector(Fields.ToArray());
                for (int s = 1; s < m_Stack_u.Length; s++) {
                    m_Stack_u[s] = new CoordinateVector(Fields.Select(dgf => dgf.CloneAs()).ToArray());
                }

                base.Residuals = new CoordinateVector(IterationResiduals.ToArray());

                for (int i = 0; i < m_Stack_u.Length; i++) {
                    if (m_PrivateBalancingInfo.m_Stack_u[i]) {
                        DGField[] _Fields = m_Stack_u[i].Mapping.Fields.ToArray();
                        for (int iF = 0; iF < _Fields.Length; iF++) {
                            L.RestoreDGField(_Fields[iF], GetName__Stack_u(i, iF));
                        }
                    }
                }

                // restore mass matrix
                m_Stack_MassMatrix = new BlockMsrMatrix[m_PrivateBalancingInfo.m_Stack_MassMatrix.Length];
                for (int i = 0; i < m_Stack_MassMatrix.Length; i++) {
                    if (m_PrivateBalancingInfo.m_Stack_MassMatrix[i]) {
                        m_Stack_MassMatrix[i] = new BlockMsrMatrix(this.CurrentStateMapping);
                        L.RestoreMatrix(m_Stack_MassMatrix[i], GetName__Stack_MassMatrix(i), CurrentStateMapping, CurrentStateMapping);
                    }
                }

                // restore operator matrix
                Debug.Assert(m_PrivateBalancingInfo.m_Stack_OpMatrix.Length == m_PrivateBalancingInfo.m_Stack_OpAffine.Length);
                m_Stack_OpMatrix = new BlockMsrMatrix[m_PrivateBalancingInfo.m_Stack_OpMatrix.Length];
                m_Stack_OpAffine = new double[m_PrivateBalancingInfo.m_Stack_OpAffine.Length][];
                for (int i = 0; i < m_Stack_OpMatrix.Length; i++) {
                    if (m_PrivateBalancingInfo.m_Stack_OpMatrix[i]) {
                        m_Stack_OpMatrix[i] = new BlockMsrMatrix(this.CurrentStateMapping);
                        L.RestoreMatrix(m_Stack_OpMatrix[i], GetName__Stack_OpMatrix(i), CurrentStateMapping, CurrentStateMapping);
                    }

                    if (m_PrivateBalancingInfo.m_Stack_OpAffine[i]) {
                        m_Stack_OpAffine[i] = new double[CurrentStateMapping.LocalLength];
                        L.RestoreVector(m_Stack_OpAffine[i], GetName__Stack_OpAffine(i));
                    }
                }

                // backup the cut-cell metrics
                m_Stack_CutCellMetrics = new CutCellMetrics[m_PrivateBalancingInfo.m_Stack_CutCellMetrics.Length];
                for (int i = 0; i < m_Stack_CutCellMetrics.Length; i++) {
                    if (m_PrivateBalancingInfo.m_Stack_CutCellMetrics[i]) {
                        L.RestoreCutCellMetrics(out m_Stack_CutCellMetrics[i], GetName__CutCellMetrics(i));
                    }
                }

                // Agglomerator
                if (m_Stack_CutCellMetrics[0] != null) {
                    CutCellMetrics[] prevCCM;
                    double[] oldAggTrsh;

                    prevCCM = m_Stack_CutCellMetrics.Skip(1).Where(ccm => ccm != null).ToArray();
                    oldAggTrsh = new double[prevCCM.Length];
                    ArrayTools.SetAll(oldAggTrsh, this.Config_AgglomerationThreshold);
                    if (prevCCM.Length <= 0) {
                        prevCCM = null;
                        oldAggTrsh = null;
                    }

                    if (!m_Stack_CutCellMetrics[0].SpeciesList.SetEquals(Config_MassScale.Keys))
                        throw new ApplicationException("Mismatch between species lists.");

                    m_CurrentAgglomeration = new MultiphaseCellAgglomerator(
                        m_Stack_CutCellMetrics[0],
                        this.Config_AgglomerationThreshold, prevCCM != null, prevCCM != null, true,
                        prevCCM, oldAggTrsh);
                }

                // finished
                m_PrivateBalancingInfo = null;
                base.MultigridSequence = _MultigridSequence;
                InitMultigrid(Fields.ToArray(), this.useX);
            }
        }


        public TimestepperInit Timestepper_Init;

        /// <summary>
        /// In case of a delayed initialization of <see cref="XdgBDFTimestepping"/>  
        /// is calling either <see cref="SingleInit"/> or <see cref="MultiInit(double, double, Action{int, double, DGField[]})"/>
        /// </summary>
        /// <param name="phystime"></param>
        /// <param name="TimestepNo"></param>
        /// <param name="dt"></param>
        /// <param name="SetTimestep">application specific action to set the previous timesteps for MultiInit, different actions for SetInitial and LoadRestart</param>
        public void DelayedTimestepperInit(double phystime, int TimestepNo, double dt, Action<int, double, DGField[]> SetTimestep) {

            if (Timestepper_Init == TimestepperInit.MultiInit) {
                MultiInit(phystime, TimestepNo, dt, SetTimestep); 
            } else {
                SingleInit();
            }

            // no increment solve for SinlgeInit and MultiInit!!!
            if (Timestepper_Init != TimestepperInit.IncrementInit)
                incrementTimesteps = 1;

        }


        public int incrementTimesteps = 1;

        int incrementHist = 0;

        CutCellMetrics[] m_Stack_CutCellMetrics_incHist;

        BlockMsrMatrix[] m_Stack_MassMatrix_incHist;

        CoordinateVector[] m_Stack_u_incHist;

        /// <summary>
        /// 
        /// </summary>
        private void InitIncrementStack() {

            int S = m_TSCchain[0].S;
            Debug.Assert(S >= 2);

            // cut-cell-metrics stack
            // ----------------------
            if (Config_LevelSetHandling == LevelSetHandling.None
                || Config_LevelSetHandling == LevelSetHandling.LieSplitting
                || Config_LevelSetHandling == LevelSetHandling.StrangSplitting)
                m_Stack_CutCellMetrics_incHist = null;
            else {
                m_Stack_CutCellMetrics_incHist = new CutCellMetrics[S - 1];
                m_Stack_CutCellMetrics_incHist[0] = m_Stack_CutCellMetrics[0];
            }

            // mass-matrix stack
            // -----------------
            if (Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                m_Stack_MassMatrix_incHist = null;
            } else if (Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity) {
                m_Stack_MassMatrix_incHist = null;
            } else {
                if (Config_LevelSetHandling == LevelSetHandling.LieSplitting
                    || Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                    m_Stack_MassMatrix_incHist = null;
                } else {
                    m_Stack_MassMatrix_incHist = new BlockMsrMatrix[S - 1];
                    m_Stack_MassMatrix_incHist[0] = m_Stack_MassMatrix[0];
                }
            }

            // m_Stack_u
            // ---------
            m_Stack_u_incHist = new CoordinateVector[S - 1];
            for (int s = 0; s < S-1; s++) {
                m_Stack_u_incHist[s] = new CoordinateVector(m_Stack_u[s].Mapping);
            }

        }



        /// <summary>
        /// Callback-routine  to update the linear resp. linearized system, 
        /// see <see cref="AssembleMatrixDel"/> resp. <see cref="NonlinearSolver.m_AssembleMatrix"/>.
        /// </summary>
        /// <param name="argCurSt">Input, current state of solution.</param>
        /// <param name="System">Output.</param>
        /// <param name="Affine">Output.</param>
        /// <param name="MassMatrix">
        /// Mass matrix including agglomeration, without any scaling,
        /// required for block-precond.
        /// </param>
        protected override void AssembleMatrixCallback(out BlockMsrMatrix System, out double[] Affine, out BlockMsrMatrix MassMatrix, DGField[] argCurSt) {
            using (new FuncTrace()) {
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
                if (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once && m_IterationCounter == 0
                    || this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {

                    MoveLevelSetAndRelatedStuff(locCurSt, m_CurrentPhystime, m_CurrentDt, 1.0);

                    // note that we need to update the agglomeration
                    updateAgglom = true;
                }

                // update agglomeration
                // --------------------
#if DEBUG
                if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                    if (m_IterationCounter == 0)
                        Debug.Assert(m_CurrentAgglomeration == null);
                    // ensure, that, when splitting is used we update the agglomerator in the very first iteration.
                }
#endif
                if (updateAgglom || m_CurrentAgglomeration == null) {
                    CutCellMetrics[] prevCCM;
                    double[] oldAggTrsh;

                    prevCCM = m_Stack_CutCellMetrics.Skip(1).Where(ccm => ccm != null).ToArray();
                    oldAggTrsh = new double[prevCCM.Length];
                    ArrayTools.SetAll(oldAggTrsh, this.Config_AgglomerationThreshold);
                    if (prevCCM.Length <= 0) {
                        prevCCM = null;
                        oldAggTrsh = null;
                    }

                    Debug.Assert(m_Stack_CutCellMetrics[0] != null);
                    if (!m_Stack_CutCellMetrics[0].SpeciesList.IsSetEqual(Config_MassScale.Keys))
                        throw new ApplicationException("Mismatch between species lists.");

                    m_CurrentAgglomeration = new MultiphaseCellAgglomerator(
                        m_Stack_CutCellMetrics[0],
                        this.Config_AgglomerationThreshold, prevCCM != null, prevCCM != null, true,
                        prevCCM, oldAggTrsh);

                    // update Multigrid-XDG basis
                    Debug.Assert(object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
                    base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
                }

                // mass matrix update
                // ------------------

                if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                    // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                    // i.e. if the Level-Set is outside the domain.

                    MassMatrix = null;
                } else {
                    // checks:
                    Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                        || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent
                        || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent,
                        "Something is not implemented here.");
                    if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                        && Config_LevelSetHandling != LevelSetHandling.None)
                        throw new NotSupportedException("Illegal configuration;");

                    Debug.Assert((m_PrecondMassMatrix == null) == (m_Stack_MassMatrix[0] == null));
                    if ((this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity && m_PrecondMassMatrix == null) // compute mass matrix (only once in application lifetime)
                        || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent && m_IterationCounter == 0) // compute mass matrix once per timestep
                        || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent) // re-compute mass matrix in every iteration
                        ) {
                        Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                        MassMatrixFactory MassFact = new MassMatrixFactory(CurrentStateMapping.BasisS.ElementAtMax(b => b.Degree), m_CurrentAgglomeration);
                        m_PrecondMassMatrix = MassFact.GetMassMatrix(CurrentStateMapping, false);
                        m_Stack_MassMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                        MassFact.AccMassMatrix(m_Stack_MassMatrix[0], CurrentStateMapping, _alpha: Config_MassScale, VariableAgglomerationSwitch: NF.ForLoop(i => false));
                    }

                    MassMatrix = m_PrecondMassMatrix;
                }


                // operator matrix update
                // ----------------------

                // we perform the extrapolation to have valid parameters if
                // - the operator matrix depends on these values
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

                // clear operator matrix (clearing and re-alloc are pretty equal, i.e. 'Clear()' just releases all internal memory)
                m_Stack_OpMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);

                // clear affine part
                if (m_Stack_OpAffine[0] == null) {
                    m_Stack_OpAffine[0] = new double[CurrentStateMapping.LocalLength];
                } else {
                    m_Stack_OpAffine[0].ClearEntries();
                }

                // assemble matrix & affine part
                Debug.Assert(m_Stack_OpMatrix[0].InfNorm() == 0);
                Debug.Assert(m_Stack_OpAffine[0].L2Norm() == 0);
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.ComputeOperatorMatrix(m_Stack_OpMatrix[0], m_Stack_OpAffine[0], CurrentStateMapping, locCurSt, m_CurrentAgglomeration, m_CurrentPhystime + m_CurrentDt);

                // assemble system
                // ---------------

                double dt = m_CurrentDt;

                double[] CurrentAffine = m_Stack_OpAffine[0];
                BlockMsrMatrix CurrentOpMatrix = m_Stack_OpMatrix[0];
                BlockMsrMatrix CurrentMassMatrix = m_Stack_MassMatrix[0];

                // choose BDF scheme 
                BDFSchemeCoeffs Tsc;
                int Smax = m_TSCchain[0].S;
                Debug.Assert(Smax == m_TSCchain.Length);
                Tsc = m_TSCchain[Smax - m_PopulatedStackDepth];

                // right-hand-side, resp. affine vector
                double[] RHS = new double[CurrentAffine.Length];
                if (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once
                    || this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // MovingMesh:
                    // (1/dt)*(M1*u1 - M0*u0) + theta1*(Op1*u1 + b1) + theta0*(Op0*u0 + b0);
                    // RHS = (1/dt)*M0*u0 - theta1*b1 - theta0*Op0*u0 - theta0*b0
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    for (int s = 1; s <= Tsc.S; s++) { // loop over BDF stages
                        if (m_Stack_MassMatrix[s] != null) {
                            m_Stack_MassMatrix[s].SpMV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s], 1.0, RHS); //   (1/dt)*M0*u0
                        } else {
                            // mass matrix is identity
                            Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                            RHS.AccV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s]);
                        }
                    }

                    RHS.AccV(-Tsc.theta1, CurrentAffine); //                                     -theta1*b1
                    if (Tsc.theta0 != 0.0) {
                        m_Stack_OpMatrix[1].SpMV(-Tsc.theta0, m_Stack_u[1], 1.0, RHS); // -theta0*Op0*u0
                        RHS.AccV(-Tsc.theta0, m_Stack_OpAffine[1]); //                    -theta0*b0 
                    }

                } else if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                    || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                    || this.Config_LevelSetHandling == LevelSetHandling.None) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // Splitting:
                    // (1/dt)M1*(u1 - u0) + theta1*(Op1*u1 + b1) + theta0*(Op1*u0 + b1);
                    // Note: only the new operator is used!!!!!!!
                    // RHS = (1/dt)*M1*u0 - theta1*b1 - theta0*Op1*u0 - theta0*b1
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    for (int s = 1; s <= Tsc.S; s++) { // loop over BDF stages
                        if (CurrentOpMatrix != null) {
                            if (CurrentMassMatrix != null) {
                                CurrentMassMatrix.SpMV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s], 1.0, RHS); //   (1/dt)*M0*u0 
                            } else {
                                Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                                RHS.AccV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s]);
                            }
                        } else {
                            Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                            RHS.AccV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s]);
                        }
                    }

                    RHS.AccV(-Tsc.theta1, CurrentAffine); //                                     -theta1*b1
                    if (Tsc.theta0 != 0.0) {
                        // For XDG & splitting, we have to use the _actual_ operator matrix for the old timestep,
                        // since the old one was created for a different interface. (Maybe)?

                        m_Stack_OpMatrix[0].SpMV(-Tsc.theta0, m_Stack_u[1], 1.0, RHS); //  -theta0*Op1*u0
                        RHS.AccV(-Tsc.theta0, m_Stack_OpAffine[0]); //                     -theta0*b1
                    }

                } else {
                    throw new NotImplementedException();
                }
                Affine = RHS;
                Affine.ScaleV(-1.0);

                // left-hand-side
                System = CurrentOpMatrix.CloneAs();
                System.Scale(Tsc.theta1);
                if (CurrentMassMatrix != null) {
                    System.Acc(1.0 / dt, CurrentMassMatrix);
                } else {
                    Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                    System.AccEyeSp(1.0 / dt);
                }

                // perform agglomeration
                // ---------------------
                Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                m_CurrentAgglomeration.ManipulateMatrixAndRHS(System, RHS, CurrentStateMapping, CurrentStateMapping);

                // increase iteration counter         
                // --------------------------

                m_IterationCounter++;
            }
        }

        double m_CurrentPhystime;
        double m_CurrentDt = -1;
        //double m_LastLevelSetResidual;
        //int m_TimestepCounter = 0;

        static double MatrixDist(MsrMatrix _A, MsrMatrix B) {
            MsrMatrix A = _A.CloneAs();
            A.Acc(-1.0, B);
            double w = A.InfNorm();
            return w;
        }


        public DelPushLevelSetRelatedStuff PushLevelSet {
            get;
            set;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <param name="timestepNumber"></param>
        /// <param name="ComputeOnlyResidual"></param>
        public void Solve(double phystime, double dt, bool ComputeOnlyResidual = false) {
            if (dt <= 0)
                throw new ArgumentOutOfRangeException();
            if (m_CurrentDt_Timestep > 0 && Math.Abs(dt / m_CurrentDt_Timestep - 1.0) > 1.0e-14)
                throw new ArgumentOutOfRangeException();

            m_CurrentDt_Timestep = dt;

            for (int i = 1; i <= incrementTimesteps; i++) {
                // push levelsets for every incremental timestep
                if (i > 1) 
                    PushLevelSet();

                // solve timestep with incremental timestep size
                double incTimestepSize = dt / (double)incrementTimesteps;
                Solve_Increment(i, phystime, incTimestepSize, ComputeOnlyResidual);
                phystime += incTimestepSize;
            }

        }

        double m_CurrentDt_Timestep = -1;

        bool OneTimeMgInit = false;

        /// <summary>
        /// Solver;
        /// </summary>
        /// <param name="phystime"></param>
        /// <param name="dt">
        /// Time-step size, must be the same value for each call in the lifetime of this object.
        /// </param>
        /// <param name="ComputeOnlyResidual">
        /// If true, no solution is performed; only the residual of the actual solution is computed.
        /// </param>
        public void Solve_Increment(int increment, double phystime, double dt, bool ComputeOnlyResidual = false) {
            if (dt <= 0)
                throw new ArgumentOutOfRangeException();
            if (m_CurrentDt > 0 && Math.Abs(dt / m_CurrentDt - 1.0) > 1.0e-14)
                throw new ArgumentOutOfRangeException();

            m_CurrentPhystime = phystime;
            m_CurrentDt = dt;
            m_IterationCounter = 0;

            PushStack(increment);
            if (incrementTimesteps == 1)
                dt = m_CurrentDt;
            else
                Console.WriteLine("Increment solve, timestep #{0}, dt = {1} ...", increment, dt);


            // update multigrid basis _once_ in object lifetime for steady level set:
            if (this.Config_LevelSetHandling == LevelSetHandling.None && OneTimeMgInit == false) {
                Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                Debug.Assert(object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
                base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
                OneTimeMgInit = true;
            }


            // ===========================================
            // update level-set (in the case of splitting)
            // ===========================================
            if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {

                double ls_dt = dt;
                if (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting)
                    ls_dt *= 0.5;

                // remember which old cells had values
                var oldCCM = this.UpdateCutCellMetrics();

                // evolve the level set
                this.MoveLevelSetAndRelatedStuff(m_Stack_u[0].Mapping.Fields.ToArray(), phystime, ls_dt, 1.0);

                // in the case of splitting, the fields must be extrapolated 
                var newCCM = this.UpdateCutCellMetrics();
                var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
                for (int i = 0; i < this.m_Stack_u.Length; i++)
                    SplittingAgg.Extrapolate(this.m_Stack_u[i].Mapping);

                // clear the old agglom; in case of splitting, the agglom does not depend on previous time-steps
                m_CurrentAgglomeration = null;
                Debug.Assert(m_Stack_CutCellMetrics.Length == 1);
                m_Stack_CutCellMetrics[0] = base.UpdateCutCellMetrics();
            }

            // ==============================================
            // solve main system
            // ==============================================

            {
                int[] Jtot =
                    (new int[] { base.m_LsTrk._Regions.GetCutCellMask().NoOfItemsLocally.MPISum(), base.m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells })
                    .MPISum();
                //Console.WriteLine("No of cells {0}, No of cut cells {1}.", Jtot[1], Jtot[0]);
                if (Jtot[0] == Jtot[1])
                    throw new ArithmeticException("All cells are cut cells - check your settings!");
            }
            if (!ComputeOnlyResidual) {

                // ++++++++++++++++++
                // normal solver run 
                // ++++++++++++++++++

                NonlinearSolver nonlinSolver;
                ISolverSmootherTemplate linearSolver;
                GetSolver(out nonlinSolver, out linearSolver);


                if (RequiresNonlinearSolver) {

                    // Nonlinear Solver (Navier-Stokes)
                    // --------------------------------

                    // use solver
                    nonlinSolver.SolverDriver(m_Stack_u[0], default(double[])); // Note: the RHS is passed as the affine part via 'this.SolverCallback'
                    
                    // 'revert' agglomeration
                    Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                    m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

                } else {
                    // Linear Solver (Stokes)
                    // ----------------------

                    // build the saddle-point matrix
                    //AssembleMatrix(this.CurrentVel, dt, phystime + dt);
                    BlockMsrMatrix System, MaMa;
                    double[] RHS;
                    this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray());
                    RHS.ScaleV(-1);

                    // update the multigrid operator
                    MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                        System, MaMa,
                        this.Config_MultigridOperator);

                    // init linear solver
                    linearSolver.Init(mgOperator);

                    // try to solve the saddle-point system.
                    mgOperator.UseSolver(linearSolver, m_Stack_u[0], RHS);

                    // 'revert' agglomeration
                    Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                    m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);
                }

            } else {
                // ++++++++++++++++++++++++++++++++++++
                // compute residual of actual solution 
                // ++++++++++++++++++++++++++++++++++++

                BlockMsrMatrix System, MaMa;
                double[] Affine;
                this.AssembleMatrixCallback(out System, out Affine, out MaMa, CurrentStateMapping.Fields.ToArray());

                base.Residuals.Clear();
                base.Residuals.SetV(Affine);
                System.SpMV(-1.0, m_Stack_u[0], -1.0, base.Residuals);

                var ResidualFields = base.Residuals.Mapping.Fields.ToArray();

                
                for (int i = 0; i < ResidualFields.Length; i++) {
                    double L2Res = ResidualFields[i].L2Norm();
                    if (m_ResLogger != null) {
                        m_ResLogger.CustomValue(L2Res, m_ResidualNames[i]);
                    } else {
                        Console.WriteLine("Residual {0}: {1}", m_ResidualNames != null ? m_ResidualNames[i] : ResidualFields[i].Identification, L2Res);
                    }
                }
                
                if (Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && m_ResLogger != null) {
                    m_ResLogger.CustomValue(0.0, "LevelSet");
                }

                if (Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && m_ResLogger != null) {
                    m_ResLogger.CustomValue(0.0, "LevelSet");
                }

                //Tecplot.Tecplot.PlotFields(ArrayTools.Cat(m_Stack_u[0].Mapping.Fields, ResidualFields), "strange", 0.0, 4);

                if (m_ResLogger != null)
                    m_ResLogger.NextIteration(true);
            }

            // ===========================================
            // update level-set (in the case of splitting)
            // ===========================================

            if (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {

                // remember which old cells had values
                var oldCCM = this.UpdateCutCellMetrics();

                // evolve the level set
                this.MoveLevelSetAndRelatedStuff(m_Stack_u[0].Mapping.Fields.ToArray(), phystime + dt * 0.5, dt * 0.5, 1.0);

                // in the case of splitting, the fields must be extrapolated 
                var newCCM = this.UpdateCutCellMetrics();
                var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
                SplittingAgg.Extrapolate(this.CurrentStateMapping);
            }
        }


        /// <summary>
        /// Visualization of multigrid solver convergence. 
        /// </summary>
        public ConvergenceObserver TestSolverOnActualSolution(string TecOutBaseName) {

            // build the saddle-point matrix
            //AssembleMatrix(this.CurrentVel, dt, phystime + dt);
            BlockMsrMatrix System, MaMa;
            double[] RHS;
            this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray());
            RHS.ScaleV(-1);

            // update the multigrid operator
            MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                System, MaMa,
                this.Config_MultigridOperator);

            // create solver
            ISolverWithCallback linearSolver = new OrthonormalizationScheme() {
                MaxIter = 50000,
                PrecondS = new ISolverSmootherTemplate[] {
                        ClassicMultigrid.InitMultigridChain(mgOperator,
                        i => new Schwarz() {
                            // this creates the pre-smoother for each level
                            m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                Depth = 1
                            },
                            overlap = 0
                        },
                        i => new Schwarz() {
                            // this creates the post-smoother for each level
                            m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                Depth = 1
                            },
                            overlap = 0
                        },
                        (i, mg) => {
                            mg.Gamma = 1;
                            mg.m_MaxIterations = 1;
                        },
                        () => new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS }) },
                Tolerance = 1.0e-10
            }; 


            // set-up the convergence observer
            double[] uEx = this.m_Stack_u[0].ToArray();
            int L = uEx.Length;
            m_CurrentAgglomeration.ClearAgglomerated(uEx, this.m_Stack_u[0].Mapping);
            var CO = new ConvergenceObserver(mgOperator, MaMa, uEx);
            uEx = null;
            CO.TecplotOut = TecOutBaseName;
            //CO.PlotDecomposition(this.u.CoordinateVector.ToArray());
            ((ISolverWithCallback)linearSolver).IterationCallback = CO.IterationCallback;
            
            // init linear solver
            linearSolver.Init(mgOperator);

            // try to solve the saddle-point system.
            mgOperator.UseSolver(linearSolver, new double[L], RHS);


            // return
            return CO;
        }


        /// <summary>
        /// Performs:
        ///  - level-set evolution
        ///  - re-sorting of matrices and vectors if the ordering of DOFs has changed.
        /// </summary>
        private void MoveLevelSetAndRelatedStuff(DGField[] locCurSt, double PhysTime, double dt, double UnderRelax) {
            // perform extrapolation:
            // If we use Agglomeration, the extrapolation is
            // also necessary for SinglePhaseFields, since we have no valid values in cells which are agglomerated.
            Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
            this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

            // level-set evolution
            int oldVersion = m_LsTrk.VersionCnt;
            m_LastLevelSetResidual = this.UpdateLevelset(locCurSt, PhysTime, dt, UnderRelax, (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting));
            int newVersion = m_LsTrk.VersionCnt;

            if ((newVersion - oldVersion) != 1)
                throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");

            // new cut-cell metric
            m_Stack_CutCellMetrics[0] = this.UpdateCutCellMetrics();
            if (!m_Stack_CutCellMetrics[0].SpeciesList.SetEquals(Config_MassScale.Keys))
                throw new ApplicationException("Mismatch between species lists.");


            // re-sort operator matrix
            if (m_Stack_OpMatrix.Length > 1) {
                TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, m_Stack_OpAffine[1], CurrentStateMapping);
            }

            // re-sort mass matrices
            {
                var Mtx2Update = new BlockMsrMatrix[m_Stack_MassMatrix.Length + m_Stack_OpMatrix.Length];
                Debug.Assert(m_Stack_OpMatrix.Length <= 2);
                if (m_Stack_OpMatrix.Length > 1)
                    Mtx2Update[0] = m_Stack_OpMatrix[1];
                for (int i = 1; i < m_Stack_MassMatrix.Length; i++) {
                    Mtx2Update[1 + i] = m_Stack_MassMatrix[i];
                }

                TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, Mtx2Update, CurrentStateMapping, CurrentStateMapping);

                if (m_Stack_OpMatrix.Length > 1)
                    m_Stack_OpMatrix[1] = Mtx2Update[0];
                for (int i = 1; i < m_Stack_MassMatrix.Length; i++) {
                    m_Stack_MassMatrix[i] = Mtx2Update[1 + i];
                }
            }
        }
    }

}