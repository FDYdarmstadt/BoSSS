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
using NUnit.Framework;

namespace BoSSS.Solution.XdgTimestepping {

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
        /// <param name="_useX">
        /// Use XDG-Fields because of Level Set. -> DIRTY HACK
        /// </param>
        /// <param name="_MultigridOperatorConfig">
        /// Configuration of block-preconditioner, if null a default value is chosen.
        /// </param>
        /// <param name="_CutCellQuadOrder">Order of quadrature in cut cells, required e.g. for <see cref="LevelSetTracker.GetXDGSpaceMetrics(SpeciesId[], int, int)"/></param>
        /// <param name="_SpId">Species to compute, actually a subset of <see cref="LevelSetTracker.SpeciesIdS"/></param>
        /// <param name="_MultigridSequence"></param>
        public XdgBDFTimestepping(IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            LevelSetTracker LsTrk,
            bool DelayInit,
            DelComputeOperatorMatrix _ComputeOperatorMatrix,
            DelUpdateLevelset _UpdateLevelset,
            int BDForder,
            LevelSetHandling _LevelSetHandling,
            MassMatrixShapeandDependence _MassMatrixShapeandDependence,
            SpatialOperatorType _SpatialOperatorType,
            IDictionary<SpeciesId, IEnumerable<double>> _MassScale,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
            AggregationGrid[] _MultigridSequence,
            SpeciesId[] _SpId,
            int _CutCellQuadOrder,
            double _AgglomerationThreshold, bool _useX) {

            if (Fields.Count() != IterationResiduals.Count())
                throw new ArgumentException("Expecting the same number of fields and residuals.");
            for (int iFld = 0; iFld < Fields.Count(); iFld++) {
                if (!Fields.ElementAt(iFld).Basis.Equals(IterationResiduals.ElementAt(iFld).Basis))
                    throw new ArgumentException(string.Format("Mismatch between {0}-th basis of fields and residuals.", iFld));
            }

            if (_MassScale != null) {
                if (!IEnumerableExtensions.SetEquals(_SpId, _MassScale.Keys))
                    throw new ArgumentException();
            }

            this.Config_LevelSetHandling = _LevelSetHandling;
            this.Config_MassMatrixShapeandDependence = _MassMatrixShapeandDependence;
            this.Config_SpatialOperatorType = _SpatialOperatorType;
            this.ComputeOperatorMatrix = _ComputeOperatorMatrix;
            this.UpdateLevelset = _UpdateLevelset;
            this.Config_MassScale = _MassScale;
            this.Config_AgglomerationThreshold = _AgglomerationThreshold;
            this.useX = _useX;
            base.MultigridSequence = _MultigridSequence;
            base.Config_SpeciesToCompute = _SpId;
            base.Config_CutCellQuadratureOrder = _CutCellQuadOrder;
            if (_MultigridSequence == null || _MultigridSequence.Length < 1)
                throw new ArgumentException("At least one grid level is required.");

            base.Residuals = new CoordinateVector(IterationResiduals.ToArray());

            if (_MultigridOperatorConfig != null) {
                base.Config_MultigridOperator = _MultigridOperatorConfig;
            } else {
                SetConfig_MultigridOperator_Default(Fields);
            }

            base.CommonConfigurationChecks();

            m_TSCchain = BDFCommon.GetChain(BDForder);

            m_LsTrk = LsTrk;

            int S = m_TSCchain[0].S;
            Debug.Assert(S == m_TSCchain.Length);

            // cut-cell-metrics stack
            // ----------------------

            if (Config_LevelSetHandling == LevelSetHandling.None) {
                m_LsTrk.IncreaseHistoryLength(0);
            } else if (Config_LevelSetHandling == LevelSetHandling.LieSplitting
                  || Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                m_LsTrk.IncreaseHistoryLength(1);
            } else {
                m_LsTrk.IncreaseHistoryLength(S + 1);
            }

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

        public int GetNumberOfStages {
            get {
                return m_TSCchain[0].S;
            }
        }

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
        /// returns the solution fields (and level-set if necessary) for older stages in case of S > 1
        /// </summary>
        /// <returns></returns>
        public ICollection<DGField>[] GetRestartInfos() {

            if(m_PopulatedStackDepth < m_TSCchain[0].S)
                return null;

            Debug.Assert(m_PopulatedStackDepth == m_TSCchain[0].S);

            ICollection<DGField>[] restartInfo = new List<DGField>[m_PopulatedStackDepth - 1];

            for(int i = 1; i < m_TSCchain[0].S; i++) {
                restartInfo[i - 1] = new List<DGField>();

                if(this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once
                    || this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {

                    DGField phiField = (DGField)m_LsTrk.LevelSetHistories[0][1 - i];
                    restartInfo[i - 1].Add(phiField);
                }

                DGField[] solFields = m_Stack_u[i].Mapping.Fields.ToArray();
                foreach(DGField f in solFields) {
                    restartInfo[i - 1].Add(f);
                }

            }

            return restartInfo;
        }     


        /// <summary>
        /// 
        /// </summary>
        void PushStack(int increment) {

            m_PopulatedStackDepth++;
            if (m_PopulatedStackDepth > m_TSCchain[0].S)
                m_PopulatedStackDepth = m_TSCchain[0].S;

            // Level-Set tracker
            // -----------------

            switch (Config_LevelSetHandling) {
                case LevelSetHandling.None:
                    // noop
                    break;

                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.StrangSplitting:
                    // noop.
                    break;

                case LevelSetHandling.Coupled_Iterative:
                case LevelSetHandling.Coupled_Once:
                    m_LsTrk.IncreaseHistoryLength(m_TSCchain[0].S);
                    m_LsTrk.PushStacks();
                    break;

                default:
                    throw new NotImplementedException();
            }

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

            //// cut-cell metrics
            //// ----------------
            //if (this.Config_LevelSetHandling != LevelSetHandling.None) {
            //    // we expect the level-set to change in every timestep

            //    for (int i = m_Stack_CutCellMetrics.Length - 1; i >= 1; i--) {
            //        m_Stack_CutCellMetrics[i] = m_Stack_CutCellMetrics[i - 1];
            //    }
            //    m_Stack_CutCellMetrics[0] = null;
            //} else {
            //    // a level-set which is static over all timeteps

            //    Debug.Assert(m_Stack_CutCellMetrics[0] != null);
            //}


            // special hack: increment init
            // ----------------------------

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

                    //// cut-cell metrics
                    //if (m_Stack_CutCellMetrics_incHist != null)
                    //    m_Stack_CutCellMetrics[i] = m_Stack_CutCellMetrics_incHist[m_TSCchain[0].S - i];

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

            //// cut-cell metrics
            //if (m_Stack_CutCellMetrics_incHist != null)
            //    m_Stack_CutCellMetrics_incHist[incrementHist] = m_Stack_CutCellMetrics[1];

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

            if (Timestepper_Init == TimeStepperInit.IncrementInit) {
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
                    (new int[] { base.m_LsTrk.Regions.GetCutCellMask().NoOfItemsLocally.MPISum(), base.m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells })
                    .MPISum();
                //Console.WriteLine("No of cells {0}, No of cut cells {1}.", Jtot[1], Jtot[0]);
                if (Jtot[0] == Jtot[1])
                    throw new ArithmeticException("All cells are cut cells - check your settings!");
            }


            // update multigrid basis _once_ in object lifetime for steady level set:
            // ----------------------
            if (this.Config_LevelSetHandling == LevelSetHandling.None && OneTimeMgInit == false) {
                m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(Config_SpeciesToCompute, Config_CutCellQuadratureOrder, Config_AgglomerationThreshold,
                    AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);

                Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                Debug.Assert(object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
                base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
                OneTimeMgInit = true;


                // matrix used for precond (must be agglomerated)
                if (this.Config_MassMatrixShapeandDependence != MassMatrixShapeandDependence.IsIdentity) {
                    MassMatrixFactory MassFact = m_LsTrk.GetXDGSpaceMetrics(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder).MassMatrixFactory;
                    m_PrecondMassMatrix = MassFact.GetMassMatrix(CurrentStateMapping, false);
                    m_CurrentAgglomeration.ManipulateMatrixAndRHS(m_PrecondMassMatrix, default(double[]), CurrentStateMapping, CurrentStateMapping);
                }
            }

            // update Multigrid-XDG basis
            //Debug.Assert(object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
            //base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);

            // mass matrix update
            // ------------------

            if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                // i.e. if the Level-Set is outside the domain.

                foreach (var kv in this.Config_MassScale) {
                    SpeciesId spId = kv.Key;
                    double[] scaleVec = kv.Value.ToArray();
                    for (int i = 0; i < scaleVec.Length; i++) {
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
                //Debug.Assert((m_PrecondMassMatrix == null) == (m_Stack_MassMatrix[0] == null));
                if ((this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity)// && m_PrecondMassMatrix == null)
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent && m_IterationCounter == 0)
                    || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent)
                    ) {
                    int NF = CurrentStateMapping.Fields.Count;
                    //MassMatrixFactory MassFact = new MassMatrixFactory(CurrentStateMapping.BasisS.ElementAtMax(b => b.Degree), m_CurrentAgglomeration);

                    // matrix for time derivative
                    MassMatrixFactory MassFact = m_LsTrk.GetXDGSpaceMetrics(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder).MassMatrixFactory;
                    m_Stack_MassMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                    MassFact.AccMassMatrix(m_Stack_MassMatrix[0], CurrentStateMapping, _alpha: Config_MassScale);

                }
            }


            // operator matrix update
            // ----------------------

            if (OpInit && m_TSCchain[0].theta0 != 0.0) {
                // we perform the extrapolation to have valid parameters if
                // - the operator matrix depends on these values
                //if(this.m_CurrentAgglomeration != null) {
                //} else {
                //    Debug.Assert(m_IterationCounter == 0 || Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                //}

                //    Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                //    this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);
                Debug.Assert(m_CurrentAgglomeration == null);

                if (m_Stack_OpMatrix[0] == null) {
                    m_Stack_OpMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                }
                if (m_Stack_OpAffine[0] == null) {
                    m_Stack_OpAffine[0] = new double[CurrentStateMapping.LocalLength];
                }

                Debug.Assert(m_Stack_OpMatrix[0].InfNorm() == 0);
                Debug.Assert(m_Stack_OpAffine[0].L2Norm() == 0);
                this.ComputeOperatorMatrix(m_Stack_OpMatrix[0], m_Stack_OpAffine[0], CurrentStateMapping, CurrentStateMapping.Fields.ToArray(), base.GetAgglomeratedLengthScales(), m_CurrentPhystime + m_CurrentDt);
            }
        }


        private string GetName__Stack_u(int i, int iF) {
            return this.GetType().FullName + "::Stack_u[" + i + "," + iF + "]";
        }



        private string GetName__Stack_OpAffine(int i) {
            return this.GetType().FullName + "::Stack_OpAffine[" + i + "]";
        }

        private string GetName__Stack_OpMatrix(int i) {
            return this.GetType().FullName + "::Stack_OpMatrix[" + i + "]";
        }

        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
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
                    }
                }
                m_Stack_MassMatrix = null;

                // agglomeration
                if (m_CurrentAgglomeration != null) {
                    m_PrivateBalancingInfo.m_Agglomeration = true;
                    m_PrivateBalancingInfo.m_Agglomeration_oldTrsh = m_CurrentAgglomeration.AgglomerationThreshold_Oldtimesteps;
                }

                // backup operator
                Debug.Assert(m_Stack_OpMatrix.Length == m_Stack_OpAffine.Length);
                m_PrivateBalancingInfo.m_Stack_Operator = new bool[m_Stack_OpMatrix.Length];
                for (int i = 0; i < m_Stack_OpMatrix.Length; i++) {
                    Debug.Assert((m_Stack_OpMatrix[i] != null) == (m_Stack_OpAffine[i] != null));
                    if ((m_Stack_OpMatrix[i] != null) || (m_Stack_OpAffine[i] != null)) {
                        m_PrivateBalancingInfo.m_Stack_Operator[i] = true;
                    }
                }
                m_Stack_OpMatrix = null;
                m_Stack_OpAffine = null;

                // Delete agglomeration
                m_CurrentAgglomeration = null;
                base.MultigridBasis = null;
                base.MultigridSequence = null;
                OneTimeMgInit = false;
            }
        }

        PrivateBalancingInfo m_PrivateBalancingInfo;

        class PrivateBalancingInfo {
            public int NoOfFields;
            public bool[] m_Stack_u;
            public bool[] m_Stack_MassMatrix;
            public bool[] m_Stack_Operator;
            public bool m_Agglomeration;
            public double[] m_Agglomeration_oldTrsh;
        }



        /// <summary>
        /// Step 2 of 2 for dynamic load balancing: restore this objects 
        /// status after the grid has been re-distributed.
        /// </summary>
        public void DataRestoreAfterBalancing(GridUpdateDataVaultBase L,
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
                            _Fields[iF].Clear();
                            L.RestoreDGField(_Fields[iF], GetName__Stack_u(i, iF));
                        }
                    }
                }

                // restore mass matrix
                m_Stack_MassMatrix = new BlockMsrMatrix[m_PrivateBalancingInfo.m_Stack_MassMatrix.Length];
                for (int i = 0; i < m_Stack_MassMatrix.Length; i++) {
                    if (m_PrivateBalancingInfo.m_Stack_MassMatrix[i]) {
                        m_Stack_MassMatrix[i] = new BlockMsrMatrix(this.CurrentStateMapping);
                        //L.RestoreMatrix(m_Stack_MassMatrix[i], GetName__Stack_MassMatrix(i), CurrentStateMapping, CurrentStateMapping);
                        m_LsTrk.GetXDGSpaceMetrics(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder, 1 - i)
                            .MassMatrixFactory
                            .AccMassMatrix(m_Stack_MassMatrix[i], CurrentStateMapping, _alpha: Config_MassScale);
                    }
                }


                // Agglomerator
                if (m_PrivateBalancingInfo.m_Agglomeration) {
                    double[] oldAggTrsh = m_PrivateBalancingInfo.m_Agglomeration_oldTrsh;
                    if (oldAggTrsh != null && oldAggTrsh.Length <= 0) {
                        oldAggTrsh = null;
                    }

                    m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                        __AgglomerationTreshold: base.Config_AgglomerationThreshold,
                        AgglomerateNewborn: (oldAggTrsh != null), AgglomerateDecased: (oldAggTrsh != null),
                        ExceptionOnFailedAgglomeration: true,
                        oldTs__AgglomerationTreshold: oldAggTrsh);

                }

                // restore operator matrix
                m_Stack_OpMatrix = new BlockMsrMatrix[m_PrivateBalancingInfo.m_Stack_Operator.Length];
                m_Stack_OpAffine = new double[m_PrivateBalancingInfo.m_Stack_Operator.Length][];
                for (int i = 0; i < m_Stack_OpMatrix.Length; i++) {
                    //if (m_PrivateBalancingInfo.m_Stack_OpMatrix[i]) {
                    //    m_Stack_OpMatrix[i] = new BlockMsrMatrix(this.CurrentStateMapping);
                    //    L.RestoreMatrix(m_Stack_OpMatrix[i], GetName__Stack_OpMatrix(i), CurrentStateMapping, CurrentStateMapping);
                    //}

                    //if (m_PrivateBalancingInfo.m_Stack_OpAffine[i]) {
                    //    m_Stack_OpAffine[i] = new double[CurrentStateMapping.LocalLength];
                    //    L.RestoreVector(m_Stack_OpAffine[i], GetName__Stack_OpAffine(i));
                    //}

                    if (!m_PrivateBalancingInfo.m_Stack_Operator[i])
                        continue;

                    m_Stack_OpMatrix[i] = new BlockMsrMatrix(CurrentStateMapping);
                    m_Stack_OpAffine[i] = new double[CurrentStateMapping.LocalLength];
                    this.ComputeOperatorMatrix(m_Stack_OpMatrix[i], m_Stack_OpAffine[i],
                        m_Stack_u[i].Mapping, m_Stack_u[i].Mapping.Fields.ToArray(), base.GetAgglomeratedLengthScales(), m_CurrentPhystime + m_CurrentDt);
                }

                // finished
                m_PrivateBalancingInfo = null;
                base.MultigridSequence = _MultigridSequence;
                InitMultigrid(Fields.ToArray(), this.useX);
            }
        }


        public TimeStepperInit Timestepper_Init;

        /// <summary>
        /// In case of a delayed initialization of <see cref="XdgBDFTimestepping"/>  
        /// is calling either <see cref="SingleInit"/> or <see cref="MultiInit(double, double, Action{int, double, DGField[]})"/>
        /// </summary>
        /// <param name="phystime"></param>
        /// <param name="TimestepNo"></param>
        /// <param name="dt"></param>
        /// <param name="SetTimestep">application specific action to set the previous timesteps for MultiInit, different actions for SetInitial and LoadRestart</param>
        public void DelayedTimestepperInit(double phystime, int TimestepNo, double dt, Action<int, double, DGField[]> SetTimestep) {

            if (Timestepper_Init == TimeStepperInit.MultiInit) {
                MultiInit(phystime, TimestepNo, dt, SetTimestep);
            } else {
                SingleInit();
            }

            // no increment solve for SinlgeInit and MultiInit!!!
            if (Timestepper_Init != TimeStepperInit.IncrementInit)
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
                //m_Stack_CutCellMetrics_incHist[0] = m_Stack_CutCellMetrics[0];
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
            for (int s = 0; s < S - 1; s++) {
                m_Stack_u_incHist[s] = new CoordinateVector(m_Stack_u[s].Mapping);
            }

        }


        /// <summary>
        /// Callback-routine  to update the linear resp. linearized system, 
        /// see <see cref="OperatorEvalOrLin"/> resp. <see cref="NonlinearSolver.m_AssembleMatrix"/>.
        /// </summary>
        /// <param name="argCurSt">Input, current state of solution.</param>
        /// <param name="System">Output.</param>
        /// <param name="Affine">Output.</param>
        /// <param name="PrecondMassMatrix">
        /// Mass matrix including agglomeration, without any scaling,
        /// required for block-precond.
        /// </param>
        /// <param name="Linearization">
        /// - true: assemble matrix and affine vector
        /// - false: evaluate operator (<paramref name="System"/> will be null)
        /// </param>
        protected override void AssembleMatrixCallback(out BlockMsrMatrix System, out double[] Affine, out BlockMsrMatrix PrecondMassMatrix, DGField[] argCurSt, bool Linearization) {
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
                if ((this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once && m_IterationCounter == 0)
                    || (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && CoupledIteration)) {

                    m_CoupledIterations++;
                    if(this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative)
                        Console.WriteLine("Coupled Iteration {0}:", m_CoupledIterations);

                    MoveLevelSetAndRelatedStuff(locCurSt, m_CurrentPhystime, m_CurrentDt, IterUnderrelax);

                    // note that we need to update the agglomeration
                    updateAgglom = true;
                }

                if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                    if (m_IterationCounter == 0) {
                        Debug.Assert(m_CurrentAgglomeration == null);
                        updateAgglom = true;
                    } else {
                        Debug.Assert(m_CurrentAgglomeration != null);
                    }
                    // ensure, that, when splitting is used we update the agglomerator in the very first iteration.
                }


                // update agglomeration
                // --------------------

                if (updateAgglom || m_CurrentAgglomeration == null) {

                    if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {
                        // Agglomeration update in the case of splitting - agglomeration does **NOT** depend on previous time-steps

                        Debug.Assert(m_IterationCounter == 0);
                        m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder, this.Config_AgglomerationThreshold,
                            AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);

                    } else {
                        // Agglomeration update in the case of a moving interface - consider previous time-steps
                        double[] oldAggTrsh;
                        if (m_PopulatedStackDepth > 0) {
                            oldAggTrsh = new double[m_PopulatedStackDepth];
                            ArrayTools.SetAll(oldAggTrsh, this.Config_AgglomerationThreshold);
                        } else {
                            oldAggTrsh = null;
                        }
                        if (this.Config_LevelSetHandling != LevelSetHandling.Coupled_Iterative) {
                            Debug.Assert(m_Stack_MassMatrix.Where(mm => mm != null).Count() == m_PopulatedStackDepth);
                        }

                        //TS++;

                        //if(TS == 3)
                        //    Console.WriteLine("break");

                        m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                            __AgglomerationTreshold: base.Config_AgglomerationThreshold,
                            AgglomerateNewborn: (oldAggTrsh != null), AgglomerateDecased: (oldAggTrsh != null),
                            ExceptionOnFailedAgglomeration: true,
                            oldTs__AgglomerationTreshold: oldAggTrsh);


                        //m_CurrentAgglomeration.PlotAgglomerationPairs("agglom-" + TS );
                        //Console.WriteLine("internal ts" + TS);
                        //Console.WriteLine("  no of agg, A {0} ", m_CurrentAgglomeration.GetAgglomerator(m_LsTrk.GetSpeciesId("A")).TotalNumberOfAgglomerations);
                        //Console.WriteLine("  no of agg, B {0} ", m_CurrentAgglomeration.GetAgglomerator(m_LsTrk.GetSpeciesId("B")).TotalNumberOfAgglomerations);
                    }


                    // update Multigrid-XDG basis
                    Debug.Assert(object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
                    base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
                }

                // mass matrix update
                // ------------------

                if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                    // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                    // i.e. if the Level-Set is outside the domain.

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

                    Debug.Assert((m_PrecondMassMatrix == null) == (m_Stack_MassMatrix[0] == null));
                    if ((this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity && m_PrecondMassMatrix == null) // compute mass matrix (only once in application lifetime)
                        || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent && m_IterationCounter == 0) // compute mass matrix once per timestep
                        || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent && CoupledIteration) // re-compute mass matrix in every iteration
                        ) {
                        Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));

                        //MassMatrixFactory MassFact = new MassMatrixFactory(CurrentStateMapping.BasisS.ElementAtMax(b => b.Degree), m_CurrentAgglomeration);
                        MassMatrixFactory MassFact = m_LsTrk.GetXDGSpaceMetrics(Config_SpeciesToCompute, Config_CutCellQuadratureOrder, 1).MassMatrixFactory;

                        // matrix for preconditioner (Agglom required)
                        m_PrecondMassMatrix = MassFact.GetMassMatrix(CurrentStateMapping, false);
                        m_CurrentAgglomeration.ManipulateMatrixAndRHS(m_PrecondMassMatrix, default(double[]), CurrentStateMapping, CurrentStateMapping);

                        // mass matrix for time derivative
                        m_Stack_MassMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                        MassFact.AccMassMatrix(m_Stack_MassMatrix[0], CurrentStateMapping, _alpha: Config_MassScale);
                    }

                    PrecondMassMatrix = m_PrecondMassMatrix;

                    CoupledIteration = true;
                }


                // operator matrix update
                // ----------------------

                // we perform the extrapolation to have valid parameters if
                // - the operator matrix depends on these values
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

                // clear operator matrix (clearing and re-alloc are pretty equal, i.e. 'BlockMsrMatrix.Clear()' just releases all internal memory)
                if(Linearization)
                    m_Stack_OpMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                else
                    m_Stack_OpMatrix[0] = null;

                // clear affine part
                if (m_Stack_OpAffine[0] == null) {
                    m_Stack_OpAffine[0] = new double[CurrentStateMapping.LocalLength];
                } else {
                    m_Stack_OpAffine[0].ClearEntries();
                }

                // assemble matrix & affine part
                Debug.Assert(m_Stack_OpMatrix[0] == null || m_Stack_OpMatrix[0].InfNorm() == 0);
                Debug.Assert(m_Stack_OpAffine[0].L2Norm() == 0);
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.ComputeOperatorMatrix(m_Stack_OpMatrix[0], m_Stack_OpAffine[0], CurrentStateMapping, locCurSt, base.GetAgglomeratedLengthScales(), m_CurrentPhystime + m_CurrentDt);

                


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
                        if(Linearization == false)
                            throw new NotImplementedException();
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
                        if (CurrentAffine != null) {
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
                        if(Linearization == false)
                            throw new NotImplementedException();
                        m_Stack_OpMatrix[0].SpMV(-Tsc.theta0, m_Stack_u[1], 1.0, RHS); //  -theta0*Op1*u0
                        RHS.AccV(-Tsc.theta0, m_Stack_OpAffine[0]); //                     -theta0*b1
                    }

                } else {
                    throw new NotImplementedException();
                }
                Affine = RHS;
                Affine.ScaleV(-1.0);

                // left-hand-side
                if(Linearization) {
                    System = CurrentOpMatrix.CloneAs();
                    if(Tsc.theta1 != 1.0)
                        System.Scale(Tsc.theta1);
                } else {
                    System = null;
                }
                if (CurrentMassMatrix != null) {
                    if(Linearization) {
                        System.Acc(1.0 / dt, CurrentMassMatrix);
                    } else {
                        CurrentMassMatrix.SpMV(1.0 / dt, new CoordinateVector(CurrentStateMapping), 1.0, Affine);
                    }
                } else {
                    Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                    if(Linearization) {
                        System.AccEyeSp(1.0 / dt);
                    } else {
                        Affine.AccV(1.0 / dt, new CoordinateVector(CurrentStateMapping));
                    }
                }

#if DEBUG
                if (Config_MassMatrixShapeandDependence != MassMatrixShapeandDependence.IsIdentity) {
                    // compare "private" and "official" mass matrix stack
                    // (private may be removed soon)

                    for (int i = 0; i < m_Stack_MassMatrix.Length; i++) {
                        var MM = m_Stack_MassMatrix[i];
                        if (MM == null)
                            continue;
                        if (i >= m_LsTrk.PopulatedHistoryLength)
                            continue;
                        var MMcomp = new BlockMsrMatrix(CurrentStateMapping);
                        m_LsTrk.GetXDGSpaceMetrics(Config_SpeciesToCompute, Config_CutCellQuadratureOrder, 1 - i)
                            .MassMatrixFactory
                            .AccMassMatrix(MMcomp, CurrentStateMapping, _alpha: Config_MassScale);
                        var sollNull = MMcomp.CloneAs();
                        sollNull.Acc(-1.0, MM);
                        double normMMcomp = sollNull.InfNorm();
                        Debug.Assert(normMMcomp == 0.0);
                    }
                }

#endif

                // perform agglomeration
                // ---------------------
                Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                m_CurrentAgglomeration.ManipulateMatrixAndRHS(System, RHS, CurrentStateMapping, CurrentStateMapping);


                // increase iteration counter         
                // --------------------------

                m_IterationCounter++;
            }
        }


        bool CoupledIteration = true;

        int m_CoupledIterations = 0;

        int m_InnerCoupledIterations = 0;


        /// <summary>
        /// callback routine for the handling of the coupled level-set iteration
        /// </summary>
        /// <param name="iterIndex"></param>
        /// <param name="currentSol"></param>
        /// <param name="currentRes"></param>
        /// <param name="Mgop"></param>
        protected void CoupledIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {

            double ResidualNorm = currentRes.L2NormPow2().MPISum().Sqrt();
            //Console.WriteLine("ResidualNorm in CoupledIterationCallback is {0}", ResidualNorm);
            // delay the update of the level-set until the flow solver converged
            if(ResidualNorm >= this.Config_SolverConvergenceCriterion) {
                this.CoupledIteration = false;
            } else {
                m_InnerCoupledIterations = 0;
            }

        }

        protected int CoupledIterationCounter(int NoIter, ref int coupledIter) {

            coupledIter = m_CoupledIterations;
            m_InnerCoupledIterations++;

            return m_InnerCoupledIterations;

        }



        double m_CurrentPhystime;
        double m_CurrentDt = -1;


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
        /// Solver.
        /// </summary>
        /// <param name="phystime">
        /// Physical time for the initial value.
        /// </param>
        /// <param name="dt">
        /// Time-step size, must be the same value for each call in the lifetime of this object.
        /// </param>
        /// <param name="ComputeOnlyResidual">
        /// If true, no solution is performed; only the residual of the actual solution is computed.
        /// </param>
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
        /// <param name="phystime">
        /// Physical time for the initial value.
        /// </param>
        /// <param name="dt">
        /// Time-step size, must be the same value for each call in the lifetime of this object.
        /// </param>
        /// <param name="ComputeOnlyResidual">
        /// If true, no solution is performed; only the residual of the actual solution is computed.
        /// </param>
        /// <param name="increment">
        /// Sub-timestep index (used during BDF startup).
        /// </param>
        void Solve_Increment(int increment, double phystime, double dt, bool ComputeOnlyResidual = false) {
            if (dt <= 0)
                throw new ArgumentOutOfRangeException();
            if (m_CurrentDt > 0 && Math.Abs(dt / m_CurrentDt - 1.0) > 1.0e-14)
                throw new ArgumentOutOfRangeException();

            m_CurrentPhystime = phystime;
            m_CurrentDt = dt;
            m_IterationCounter = 0;
            m_CoupledIterations = 0;
            m_InnerCoupledIterations = 0;

            PushStack(increment);
            if (incrementTimesteps == 1)
                dt = m_CurrentDt;
            else
                Console.WriteLine("Increment solve, timestep #{0}, dt = {1} ...", increment, dt);


            // ===========================================
            // update level-set (in the case of splitting)
            // ===========================================
            if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {

                Debug.Assert(m_CurrentAgglomeration == null);

                double ls_dt = dt;
                if (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting)
                    ls_dt *= 0.5;

                // remember which old cells had values
                //var oldCCM = this.UpdateCutCellMetrics();

                // evolve the level set
                m_LsTrk.IncreaseHistoryLength(1);
                m_LsTrk.PushStacks();

                int oldPushCount = m_LsTrk.PushCount;
                int oldVersion = m_LsTrk.VersionCnt;

                this.MoveLevelSetAndRelatedStuff(m_Stack_u[0].Mapping.Fields.ToArray(), phystime, ls_dt, 1.0);

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
                for (int i = 0; i < this.m_Stack_u.Length; i++)
                    SplittingAgg.Extrapolate(this.m_Stack_u[i].Mapping);

                // delete new agglomeration; in case of splitting, the agglomeration for the **bulk operator timestep** does not depend on previous time-steps
                m_CurrentAgglomeration = null;
            }

            // ==============================================
            // solve main system
            // ==============================================

            int oldLsTrkPushCount = m_LsTrk.PushCount;

            {
                int[] Jtot =
                    (new int[] { base.m_LsTrk.Regions.GetCutCellMask().NoOfItemsLocally.MPISum(), base.m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells })
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
                    this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true);
                    RHS.ScaleV(-1);

                    // update the multigrid operator
                    MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                        System, MaMa,
                        this.Config_MultigridOperator);

                    using (var tr = new FuncTrace()) {
                        // init linear solver
                        using (new BlockTrace("Slv Init", tr)) {
                            linearSolver.Init(mgOperator);
                        }

                        // try to solve the saddle-point system.
                        using (new BlockTrace("Slv Iter", tr)) {
                            mgOperator.UseSolver(linearSolver, m_Stack_u[0], RHS);
                        }
                    }

                    // 'revert' agglomeration
                    Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                    m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);
                }

            } else {
                // ++++++++++++++++++++++++++++++++++++
                // compute residual of actual solution 
                // ++++++++++++++++++++++++++++++++++++

                
                double[] Affine;
                this.AssembleMatrixCallback(out BlockMsrMatrix System, out Affine, out BlockMsrMatrix MaMa, CurrentStateMapping.Fields.ToArray(), false);
                Debug.Assert(System == null);

                base.Residuals.Clear();
                base.Residuals.SetV(Affine, -1.0);
                //System.SpMV(-1.0, m_Stack_u[0], +1.0, base.Residuals);

#if DEBUG
                {

                    this.AssembleMatrixCallback(out BlockMsrMatrix checkSystem, out double[] checkAffine, out BlockMsrMatrix MaMa1, CurrentStateMapping.Fields.ToArray(), true);

                    double[] checkResidual = new double[checkAffine.Length];
                    checkResidual.SetV(checkAffine, -1.0);
                    checkSystem.SpMV(-1.0, m_Stack_u[0], +1.0, checkResidual);

                    double distL2 = GenericBlas.L2DistPow2(checkResidual, base.Residuals).MPISum().Sqrt();
                    double refL2 = (new double[] { GenericBlas.L2NormPow2(m_Stack_u[0]), GenericBlas.L2NormPow2(checkResidual), GenericBlas.L2NormPow2(base.Residuals) }).MPISum().Max().Sqrt();

                    //Assert.Less(distL2, refL2 * 1.0e-5, "argh");

                }
#endif


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

                //if (Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && m_ResLogger != null) {
                //    m_ResLogger.CustomValue(0.0, "LevelSet");
                //}

                //Tecplot.Tecplot.PlotFields(ArrayTools.Cat(m_Stack_u[0].Mapping.Fields, ResidualFields), "strange", 0.0, 4);

                if (m_ResLogger != null)
                    m_ResLogger.NextIteration(true);
            }

            int newLsTrkPushCount = m_LsTrk.PushCount;
            if (newLsTrkPushCount != oldLsTrkPushCount)
                throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");

            if (Config_LevelSetHandling != LevelSetHandling.None)
                m_CurrentAgglomeration = null;

            // ===========================================
            // update level-set (in the case of splitting)
            // ===========================================

            if (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {

                Debug.Assert(m_CurrentAgglomeration == null);

                // evolve the level set
                m_LsTrk.IncreaseHistoryLength(1);
                m_LsTrk.PushStacks();

                int oldPushCount = m_LsTrk.PushCount;
                int oldVersion = m_LsTrk.VersionCnt;

                this.MoveLevelSetAndRelatedStuff(m_Stack_u[0].Mapping.Fields.ToArray(), phystime + dt * 0.5, dt * 0.5, 1.0);

                int newPushCount = m_LsTrk.PushCount;
                int newVersion = m_LsTrk.VersionCnt;
                if ((newPushCount - oldPushCount) != 0)
                    throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");
                if ((newVersion - oldVersion) != 1)
                    throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");

                // in the case of splitting, the fields must be extrapolated 
                //var newCCM = this.UpdateCutCellMetrics();
                //var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
                //SplittingAgg.Extrapolate(this.CurrentStateMapping);
                Debug.Assert(m_LsTrk.HistoryLength >= 1);
                var SplittingAgg = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                    __AgglomerationTreshold: 0.0, AgglomerateNewborn: true, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true,
                    oldTs__AgglomerationTreshold: new double[] { 0.0 });
                for (int i = 0; i < this.m_Stack_u.Length; i++)
                    SplittingAgg.Extrapolate(this.m_Stack_u[i].Mapping);

                Debug.Assert(m_CurrentAgglomeration == null);
            }


            // ====================
            // release end-of-stack
            // ====================
            int ie = m_Stack_OpMatrix.Length - 1;
            Debug.Assert(m_Stack_OpMatrix.Length == m_Stack_OpAffine.Length);
            //Debug.Assert((m_Stack_OpMatrix[ie] == null) == (m_Stack_OpAffine[ie] == null));
            m_Stack_OpMatrix[ie] = null;
            m_Stack_OpAffine[ie] = null;
            //m_Stack_MassMatrix[m_Stack_MassMatrix.Length - 1] = null;
        }


        /// <summary>
        /// Visualization of multigrid solver convergence. 
        /// </summary>
        public ConvergenceObserver TestSolverOnActualSolution(string TecOutBaseName) {

            // build the saddle-point matrix
            //AssembleMatrix(this.CurrentVel, dt, phystime + dt);
            BlockMsrMatrix System, MaMa;
            double[] RHS;
            this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true);
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
                            Overlap = 0
                        },
                        i => new Schwarz() {
                            // this creates the post-smoother for each level
                            m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                Depth = 1
                            },
                            Overlap = 0
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


        public event Action<int, double[], double[], MultigridOperator> CustomIterationCallback;


        protected override string GetSolver(out NonlinearSolver nonlinSolver, out ISolverSmootherTemplate linearSolver) {
            string description = base.GetSolver(out nonlinSolver, out linearSolver);

            if (RequiresNonlinearSolver) {

                if(this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                    nonlinSolver.IterationCallback += this.CoupledIterationCallback;
                    if(nonlinSolver is FixpointIterator)
                        ((FixpointIterator)nonlinSolver).Iteration_Count = this.CoupledIterationCounter;
                }

                nonlinSolver.IterationCallback += this.CustomIterationCallback;

            } else {
                if (linearSolver is ISolverWithCallback) {
                    ((ISolverWithCallback)linearSolver).IterationCallback += this.CustomIterationCallback;
                }
            }

            return description;
        }

        public bool coupledOperator = false;

        /// <summary>
        /// Performs:
        ///  - level-set evolution
        ///  - re-sorting of matrices and vectors if the ordering of DOFs has changed.
        /// </summary>
        private void MoveLevelSetAndRelatedStuff(DGField[] locCurSt, double PhysTime, double dt, double UnderRelax) {
            if (Config_LevelSetHandling == LevelSetHandling.None) {
                throw new ApplicationException("internal error");
            }
            if (m_IterationCounter > 0) {
                if (Config_LevelSetHandling != LevelSetHandling.Coupled_Iterative)
                    throw new ApplicationException("internal error");
            }

            // perform extrapolation:
            // If we use Agglomeration, the extrapolation is
            // also necessary for SinglePhaseFields, since we have no valid values in cells which are agglomerated.
            if (m_CurrentAgglomeration != null) {
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);
            } else {
                //if(m_IterationCounter > 1) {
                //    Debug.Assert(Config_LevelSetHandling == LevelSetHandling.LieSplitting || Config_LevelSetHandling == LevelSetHandling.StrangSplitting);
                //}
                Debug.Assert(m_IterationCounter == 0 || Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
            }

            // level-set evolution
            int oldVersion = m_LsTrk.VersionCnt;
            int oldPushCount = m_LsTrk.PushCount;

            m_LastLevelSetResidual = this.UpdateLevelset(locCurSt, PhysTime, dt, UnderRelax, (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting));

            int newVersion = m_LsTrk.VersionCnt;
            int newPushCount = m_LsTrk.PushCount;


            if ((newVersion - oldVersion) != 1 && !coupledOperator)
                throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");
            if((newVersion - oldVersion) != 0 && coupledOperator)
                throw new ApplicationException("Expecting exactly no call to 'UpdateTracker(...)' in 'UpdateLevelset(...)' for coupled Operators.");
            if ((newPushCount - oldPushCount) != 0)
                throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");


            //// new cut-cell metric
            //m_Stack_CutCellMetrics[0] = this.UpdateCutCellMetrics();
            //if (!m_Stack_CutCellMetrics[0].SpeciesList.SetEquals(Config_MassScale.Keys))
            //    throw new ApplicationException("Mismatch between species lists.");


            // re-sort mass matrices 
            // and operator matrix (Exp. Euler or Crank-Nicolson)
            {
                var Mtx2Update = new BlockMsrMatrix[m_Stack_MassMatrix.Length + m_Stack_OpMatrix.Length];
                Debug.Assert(m_Stack_OpMatrix.Length <= 2);
                if (m_Stack_OpMatrix.Length > 1)
                    Mtx2Update[0] = m_Stack_OpMatrix[1];
                for (int i = 1; i < m_Stack_MassMatrix.Length; i++) {
                    Mtx2Update[i + 1] = m_Stack_MassMatrix[i];
                }

                TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, Mtx2Update, CurrentStateMapping, CurrentStateMapping);

                if (m_Stack_OpMatrix.Length > 1)
                    m_Stack_OpMatrix[1] = Mtx2Update[0];
                for (int i = 1; i < m_Stack_MassMatrix.Length; i++) {
                    m_Stack_MassMatrix[i] = Mtx2Update[i + 1];
                }
            }


            // re-sort RHS  (Exp. Euler or Crank-Nicolson)
            if (m_Stack_OpMatrix.Length > 1) {
                TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, m_Stack_OpAffine[1], CurrentStateMapping);
            }


        }
    }

}