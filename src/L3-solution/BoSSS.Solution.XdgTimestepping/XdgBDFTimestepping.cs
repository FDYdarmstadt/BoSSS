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
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Timestepping;
using ilPSP;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using MPI.Wrappers;
using NUnit.Framework;
using BoSSS.Solution.Gnuplot;
using System.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Queries;
using NUnit.Framework.Constraints;

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
        /// <param name="BDForder">
        /// The order of the BDF scheme from 1 to 6; in addition, 0 encodes Explicit Euler and -1 encodes Crank-Nicolson.
        /// </param>
        /// <param name="_LevelSetHandling"></param>
        /// <param name="_MassMatrixShapeandDependence"></param>
        /// <param name="_SpatialOperatorType"></param>
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
        /// <param name="linearconfig"></param>
        /// <param name="nonlinconfig"></param>
        /// <param name="abstractOperator"></param>
        /// <param name="__Parameters"></param>
        public XdgBDFTimestepping(
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> __Parameters,
            IEnumerable<DGField> IterationResiduals,
            LevelSetTracker LsTrk,
            bool DelayInit,
            DelComputeOperatorMatrix _ComputeOperatorMatrix,
            ISpatialOperator abstractOperator,
            Func<ISlaveTimeIntegrator> _UpdateLevelset,
            int BDForder,
            LevelSetHandling _LevelSetHandling,
            MassMatrixShapeandDependence _MassMatrixShapeandDependence,
            SpatialOperatorType _SpatialOperatorType,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
            AggregationGridData[] _MultigridSequence,
            SpeciesId[] _SpId,
            int _CutCellQuadOrder,
            double _AgglomerationThreshold, bool _useX, Control.NonLinearSolverConfig nonlinconfig,
            ISolverFactory linearconfig) 
            : base(nonlinconfig, linearconfig) //
        {

            m_nonlinconfig = nonlinconfig;

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
            this.AbstractOperator = AbstractOperator;
            this.Config_AgglomerationThreshold = _AgglomerationThreshold;
            this.useX = _useX;
            base.MultigridSequence = _MultigridSequence;
            base.Config_SpeciesToCompute = _SpId;
            base.Config_CutCellQuadratureOrder = _CutCellQuadOrder;
            base.CurrentParameters = __Parameters.ToArray();
            base.AbstractOperator = abstractOperator;

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

            base.m_LsTrk = LsTrk;

            int S = m_TSCchain[0].S;
            Debug.Assert(S == m_TSCchain.Length);

            // cut-cell-metrics stack
            // ----------------------

            if (Config_LevelSetHandling == LevelSetHandling.None) {
                m_LsTrk.IncreaseHistoryLength(0);
            } else if (Config_LevelSetHandling == LevelSetHandling.LieSplitting
                  || Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                  || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {
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
                    || Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                    || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {
                    m_Stack_MassMatrix = new BlockMsrMatrix[1];
                } else {
                    m_Stack_MassMatrix = new BlockMsrMatrix[S + 1];
                }
            }

          

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

            //if (!DelayInit)
            InitTimestepping(true);
        }

        BDFSchemeCoeffs[] m_TSCchain;

        /// <summary>
        /// 1 for implicit/explicit Euler, Crank-Nicholson; 2 for BDF2, 3 for BDF3, etc. 
        /// </summary>
        public int GetNumberOfStages {
            get {
                return m_TSCchain[0].S;
            }
        }

        /// <summary>
        /// DG coefficient mapping for the test- and trial-space.
        /// </summary>
        public override CoordinateMapping CurrentStateMapping {
            get {
                return m_Stack_u[0].Mapping;
            }
        }


        /// <summary>
        /// Switch for using XDG-Fields
        /// </summary>
        internal bool useX;

        //
        // stack of mass matrices (matrices _without_ agglomeration)
        //
        BlockMsrMatrix[] m_Stack_MassMatrix;


        //
        // stack of solution vectors
        //
        CoordinateVector[] m_Stack_u;


        int m_PopulatedStackDepth = 0;




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
        void PushStack() {

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
                case LevelSetHandling.FSILieSplittingFullyCoupled:
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


            // Solution-Stack
            // --------------
            if(m_CurrentPhystime != fsiOldPhystime || Config_LevelSetHandling != LevelSetHandling.FSILieSplittingFullyCoupled) { // only true in case of fsi_splitting fully coupled
                // entry 0 should remain the same object all the time

                var Cvtmp = m_Stack_u[m_Stack_u.Length - 1];
                for (int i = m_Stack_u.Length - 1; i >= 2; i--) {
                    m_Stack_u[i] = m_Stack_u[i - 1];
                }
                m_Stack_u[1] = Cvtmp;
                m_Stack_u[1].Clear();
                m_Stack_u[1].Acc(1.0, m_Stack_u[0]);
            }

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

        }

        internal void PopStack() {
            //m_PopulatedStackDepth++;
            //if (m_PopulatedStackDepth > m_TSCchain[0].S)
            //    m_PopulatedStackDepth = m_TSCchain[0].S;

            m_PopulatedStackDepth--;


            /*
            // Level-Set tracker
            // -----------------

            switch (Config_LevelSetHandling) {
                case LevelSetHandling.None:
                    // noop
                    break;

                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.StrangSplitting:
                case LevelSetHandling.FSILieSplittingFullyCoupled:
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
            */



            // Solution-Stack
            // --------------
            if(m_CurrentPhystime != fsiOldPhystime || Config_LevelSetHandling != LevelSetHandling.FSILieSplittingFullyCoupled) { // only true in case of fsi_splitting fully coupled
                // entry 0 should remain the same object all the time

                m_Stack_u[0].Clear();
                m_Stack_u[0].Acc(1.0, m_Stack_u[1]);

                if(m_Stack_u.Length > 1) {
                    var Cvtmp = m_Stack_u[1];
                    for(int i = 1; i <= m_Stack_u.Length - 2; i++) {
                        m_Stack_u[i] = m_Stack_u[i + 1];
                    }
                    m_Stack_u[m_Stack_u.Length - 1] = Cvtmp;
                }
            }

            // mass-matrix stack
            // -----------------


            if (this.Config_MassMatrixShapeandDependence != MassMatrixShapeandDependence.IsNonIdentity) {
                for (int i = 0; i <= m_Stack_MassMatrix.Length - 2; i++) {
                    m_Stack_MassMatrix[i] = m_Stack_MassMatrix[i + 1];
                }
                m_Stack_MassMatrix[m_Stack_MassMatrix.Length - 1] = null;
                m_PrecondMassMatrix = null;
            } else {
                // should already be initialized - see 'InitTimestepping'
                Debug.Assert(m_Stack_MassMatrix.Length == 1);
                Debug.Assert(m_Stack_MassMatrix[0] != null);
                Debug.Assert(m_PrecondMassMatrix != null);
            }
        }


        int m_IterationCounter = 0;

        bool initialized = false;

        /// <summary>
        /// Initialization from a single timestep, i.e. if this time-stepper should use BDF4,
        /// it starts with BDF1, BDF2, BDF3 in the first, second and third time-step.
        /// </summary>
        /// <remarks>
        /// This approach in general does not provide the desired convergence order, but is applicable 
        /// is the simulation does not depend on the initial value.
        /// </remarks>
        public void SingleInit() {
            using (new FuncTrace()) {
                InitTimestepping(true);

               

                initialized = true;
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

                    if (iStage < (S - 1)) // push, but not the last time in the loop
                        PushStack();
                }

                initialized = true;
            }
        }


        /// <summary>
        /// final initialization for the BDF timestepper scheme, all necessary timesteps have to be initialized
        /// </summary>
        /// <param name="OpInit"></param>
        private void InitTimestepping(bool OpInit) {

            {
                int[] Jtot =
                    (new int[] { base.m_LsTrk.Regions.GetCutCellMask().NoOfItemsLocally, base.m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells })
                    .MPISum();
                //Console.WriteLine("No of cells {0}, No of cut cells {1}.", Jtot[1], Jtot[0]);
                if (Jtot[0] == Jtot[1]) {
                    
                    Console.Error.WriteLine($"MPI rank {this.m_LsTrk.GridDat.MpiRank}: NoOfItems = {Jtot[1]}, NoOfCell = {Jtot[1]}");
                    var CC = new SinglePhaseField(new Basis(this.m_LsTrk.GridDat, 0), "CutCells");
                    CC.AccConstant(1.0, base.m_LsTrk.Regions.GetCutCellMask());

                    var FieldsToPlot = new List<DGField>();
                    FieldsToPlot.Add(CC);

                    for (int iLs = 0; iLs < this.m_LsTrk.NoOfLevelSets; iLs++) {

                        var CC_iLs = new SinglePhaseField(new Basis(this.m_LsTrk.GridDat, iLs), $"CutCells-Ls{iLs}");
                        CC_iLs.AccConstant(1.0, base.m_LsTrk.Regions.GetCutCellMask4LevSet(iLs));
                        FieldsToPlot.Add(CC_iLs);

                        if (m_LsTrk.LevelSets[iLs] is DGField dgLs) {
                            FieldsToPlot.Add(dgLs);
                        }

                    }
                    Tecplot.Tecplot.PlotFields(FieldsToPlot.ToArray(), "Error", 0.0, 2);
                    
                    throw new ArithmeticException("All cells are cut cells - check your settings!");
                }
            }


            // update multigrid basis _once_ in object lifetime for steady level set:
            // ----------------------
            if (this.Config_LevelSetHandling == LevelSetHandling.None && OneTimeMgInit == false) {
                m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(Config_SpeciesToCompute, Config_CutCellQuadratureOrder, Config_AgglomerationThreshold,
                    AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true, Tag: "InitTimestepping");

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
                    //MassMatrixFactory MassFact = m_LsTrk.GetXDGSpaceMetrics(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder).MassMatrixFactory;
                    m_Stack_MassMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                    //MassFact.AccMassMatrix(m_Stack_MassMatrix[0], CurrentStateMapping, _alpha: Config_MassScale);
                    base.ComputeMassMatrixImpl(m_Stack_MassMatrix[0], m_LsTrk.RegionsHistory.Current.Time);

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

                /*
                if (m_Stack_OpMatrix[0] == null) {
                    m_Stack_OpMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                }
                if (m_Stack_OpAffine[0] == null) {
                    m_Stack_OpAffine[0] = new double[CurrentStateMapping.LocalLength];
                }

                Debug.Assert(m_Stack_OpMatrix[0].InfNorm() == 0);
                Debug.Assert(m_Stack_OpAffine[0].L2Norm() == 0);
                this.ComputeOperatorMatrix(m_Stack_OpMatrix[0], m_Stack_OpAffine[0], CurrentStateMapping, CurrentStateMapping.Fields.ToArray(), base.GetAgglomeratedLengthScales(), m_CurrentPhystime + m_CurrentDt);
                */
            }
        }


        private string GetName__Stack_u(int i, int iF) {
            if(!coupledOperator)
                return this.GetType().FullName + "::Stack_u[" + i + "," + iF + "]";
            else
                return this.GetType().FullName + "::CoupledStack_u[" + i + "," + iF + "]";
        }



        private string GetName__Stack_OpAffine(int i) {
            if(!coupledOperator)
                return this.GetType().FullName + "::Stack_OpAffine[" + i + "]";
            else
                return this.GetType().FullName + "::CoupledStack_OpAffine[" + i + "]";
        }

        private string GetName__Stack_OpMatrix(int i) {
            if(!coupledOperator)
                return this.GetType().FullName + "::Stack_OpMatrix[" + i + "]";
            else
                return this.GetType().FullName + "::CoupledStack_OpMatrix[" + i + "]";
        }

        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public void DataBackupBeforeBalancing(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            using (new FuncTrace()) {
                if (m_PrivateBalancingInfo != null)
                    throw new NotSupportedException("Method has already been called without matching call to `DataRestoreAfterBalancing`");
                
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

                /*
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
                */

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
            //public bool[] m_Stack_Operator;
            public bool m_Agglomeration;
            public double[] m_Agglomeration_oldTrsh;
        }



        /// <summary>
        /// Step 2 of 2 for dynamic load balancing: restore this objects 
        /// status after the grid has been re-distributed.
        /// </summary>
        public void DataRestoreAfterBalancing(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> Parameters,
            IEnumerable<DGField> IterationResiduals,
            LevelSetTracker LsTrk,
            AggregationGridData[] _MultigridSequence,
            ISpatialOperator abstractOperator) //
        {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = false;

                if (m_PrivateBalancingInfo == null)
                    throw new NotSupportedException();

                base.m_LsTrk = LsTrk;

                base.AbstractOperator = abstractOperator;

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

                base.CurrentParameters = Parameters != null ? Parameters.ToArray() : new DGField[0];

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

                        if(i == 0 // most recent time-step: mass matrix for all circumstances
                            || Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative
                            || Config_LevelSetHandling == LevelSetHandling.Coupled_Once) {
                            // in general, only for moving interface (i.e. Coupled_*)
                         
                            tr.Info($"Restoring Mass Matrix after AMR or load-balancing for time-index {1 - i}, should be physical time {m_CurrentPhystime - m_CurrentDt * i}");
                            m_Stack_MassMatrix[i] = new BlockMsrMatrix(this.CurrentStateMapping);
                            base.ComputeMassMatrixImpl(m_Stack_MassMatrix[i], LsTrk.RegionsHistory[1 - i].Time);
                        } else {
                            // for Splitting, None: only mass matrix for most recent timestep

                            m_Stack_MassMatrix[i] = m_Stack_MassMatrix[0];
                        }
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
                        oldTs__AgglomerationTreshold: oldAggTrsh, Tag: "DataRestoreAfterBalancing");

                }

                // finished
                m_PrivateBalancingInfo = null;
                base.MultigridSequence = _MultigridSequence;
                InitMultigrid(Fields.ToArray(), this.useX);

                // in case of steady level set the xdgAggBasis need to be updated
                if (this.Config_LevelSetHandling == LevelSetHandling.None && OneTimeMgInit == false) {
                    Debug.Assert(m_CurrentAgglomeration == null || object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                    Debug.Assert(m_CurrentAgglomeration == null || object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
                    base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);

                    // matrix used for precond (must be agglomerated)
                    if (this.Config_MassMatrixShapeandDependence != MassMatrixShapeandDependence.IsIdentity) {
                        MassMatrixFactory MassFact = m_LsTrk.GetXDGSpaceMetrics(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder).MassMatrixFactory;
                        m_PrecondMassMatrix = MassFact.GetMassMatrix(CurrentStateMapping, false);
                        if (this.m_CurrentAgglomeration != null) {
                            m_CurrentAgglomeration.ManipulateMatrixAndRHS(m_PrecondMassMatrix, default(double[]), CurrentStateMapping, CurrentStateMapping);
                        }
                    }
                }

            }
        }



        public void ResetDataAfterBalancing(IEnumerable<DGField> Fields)
            //IEnumerable<DGField> IterationResiduals
            //LevelSetTracker LsTrk,
            //AggregationGridData[] _MultigridSequence) //
        {
            using (new FuncTrace()) {
                m_Stack_u[0] = new CoordinateVector(Fields.ToArray());

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
        /// <param name="abstractOperator">
        ///  the original operator that somehow produced the matrix; yes, this API is convoluted piece-of-shit
        /// </param>
        internal protected override void AssembleMatrixCallback(out BlockMsrMatrix System, out double[] Affine, out BlockMsrMatrix PrecondMassMatrix, DGField[] argCurSt, bool Linearization, out ISpatialOperator abstractOperator) {
            using (var tr = new FuncTrace()) {

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
                    || (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && m_IterationCounter == 0)
                    || (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && CoupledIteration)) {

                    m_CoupledIterations++;
                    if(this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative)
                        Console.WriteLine("Coupled Iteration {0}:", m_CoupledIterations);

                    MoveLevelSetAndRelatedStuff(locCurSt, m_CurrentPhystime, m_CurrentDt, IterUnderrelax);

                    // note that we need to update the agglomeration
                    updateAgglom = true;
                }

                if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                    || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {
                    if (m_IterationCounter == 0) {
                        if(m_CurrentAgglomeration != null)
                            throw new ApplicationException();
                        updateAgglom = true;
                    } else {
                        if (m_CurrentAgglomeration == null)
                            Console.WriteLine("throw new ApplicationException();");
                            //throw new ApplicationException();
                    }
                    // ensure, that, when splitting is used we update the agglomerator in the very first iteration.
                }


                // update agglomeration
                // --------------------

                if (updateAgglom || m_CurrentAgglomeration == null) {

                    if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                        || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {
                        // Agglomeration update in the case of splitting - agglomeration does **NOT** depend on previous time-steps

                        Debug.Assert(m_IterationCounter == 0);
                        m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder, this.Config_AgglomerationThreshold,
                            AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true, Tag: "AssembleMatrixCallback");

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

                        m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                            __AgglomerationTreshold: base.Config_AgglomerationThreshold,
                            AgglomerateNewborn: (oldAggTrsh != null), AgglomerateDecased: (oldAggTrsh != null),
                            ExceptionOnFailedAgglomeration: true,
                            oldTs__AgglomerationTreshold: oldAggTrsh,
                            Tag: "AssembleMatrixCallback");
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
                        base.ComputeMassMatrixImpl(m_Stack_MassMatrix[0], m_LsTrk.RegionsHistory.Current.Time);
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
                BlockMsrMatrix OpMatrix;
                if(Linearization)
                    OpMatrix = new BlockMsrMatrix(CurrentStateMapping);
                else
                    OpMatrix = null;

                /*
                void FillMatrixWithRandomShit(BlockMsrMatrix OpMtx) {
                    Random rnd = new Random();
                    double[] buf = 10000.ForLoop(i => rnd.NextDouble());

                    int c = 0;
                    int J = this.m_LsTrk.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                    var bs = CurrentStateMapping.BasisS;
                    for(int row_j = 0; row_j < J; row_j++) {
                        this.m_LsTrk.GridDat.GetCellNeighbours(row_j, GetCellNeighbours_Mode.ViaEdges, out int[] Neighs, out _);
                        row_j.AddToArray(ref Neighs);

                        foreach(int col_j in Neighs) {
                            for(int rowVar = 0; rowVar < bs.Count; rowVar++) {
                                for(int colVar = 0; colVar < bs.Count; colVar++) {
                                    int N = bs[rowVar].GetLength(row_j);
                                    int M = bs[colVar].GetLength(col_j);

                                    for(int n = 0; n < N; n++) {
                                        for(int m = 0; m < M; m++) {
                                            long iRow = CurrentStateMapping.GlobalUniqueCoordinateIndex(rowVar, row_j, n);
                                            long iCol = CurrentStateMapping.GlobalUniqueCoordinateIndex(colVar, col_j, m);

                                            OpMtx[iRow, iCol] = buf[c];
                                            c++;
                                            if(c >= buf.Length)
                                                c = 0;
                                        }
                                    }
                                }
                            }
                        }

                    }
                }

                if(OpMatrix != null) {
                    FillMatrixWithRandomShit(OpMatrix);
                }
                //*/

                // clear affine part
                double[] OpAffine = new double[CurrentStateMapping.LocalLength];


                
                // assemble matrix & affine part
                Debug.Assert(OpMatrix == null || OpMatrix.InfNorm() == 0);
                Debug.Assert(OpAffine.L2Norm() == 0);
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.ComputeOperatorMatrix(OpMatrix, OpAffine, CurrentStateMapping, locCurSt, base.GetAgglomeratedLengthScales(), m_CurrentPhystime + m_CurrentDt, 1);

                


                // assemble system
                // ---------------

                double dt = m_CurrentDt;

                double[] CurrentAffine = OpAffine;
                BlockMsrMatrix CurrentOpMatrix = OpMatrix;
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
                        double[] evalBuffer = new double[RHS.Length];
                        this.ComputeOperatorMatrix(null, evalBuffer, m_Stack_u[1].Mapping, m_Stack_u[1].Fields.ToArray(), base.GetAgglomeratedLengthScales(), m_CurrentPhystime, 0);
                        RHS.AccV(-Tsc.theta0, evalBuffer); // RHS -= -theta0*Op0(u0)
                    }

                } else if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                    || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                    || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled
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
                        // For XDG & splitting, we have to evaluate the operator with 
                        // the _actual_ level set position, but 
                        // everything else (field state, etc.) from at the _old_ timestep.
                                                
                        double[] evalBuffer = new double[RHS.Length];
                        this.ComputeOperatorMatrix(null, evalBuffer, m_Stack_u[1].Mapping, m_Stack_u[1].Fields.ToArray(), base.GetAgglomeratedLengthScales(), m_CurrentPhystime, 1);
                        RHS.AccV(-Tsc.theta0, evalBuffer); // RHS -= -theta0*Op(u0) at current level-set
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

                // perform agglomeration
                // ---------------------
                Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                m_CurrentAgglomeration.ManipulateMatrixAndRHS(System, Affine, CurrentStateMapping, CurrentStateMapping);

                if(Linearization) {
                    m_LsTrk.CheckVectorZeroInEmptyCutCells(Affine, CurrentStateMapping, this.Config_SpeciesToCompute, m_CurrentAgglomeration, this.Config_CutCellQuadratureOrder);
                    m_LsTrk.CheckMatrixZeroInEmptyCutCells(System, CurrentStateMapping, this.Config_SpeciesToCompute, m_CurrentAgglomeration, this.Config_CutCellQuadratureOrder);
                } else {
                    m_LsTrk.CheckVectorZeroInEmptyCutCells(Affine, CurrentStateMapping, this.Config_SpeciesToCompute, m_CurrentAgglomeration, this.Config_CutCellQuadratureOrder);
                }

                // increase iteration counter         
                // --------------------------
                abstractOperator = AbstractOperator;
                m_IterationCounter++;
            }
        }

        /// <summary>
        /// The time associated with the current solution (<see cref="CurrentState"/>)
        /// </summary>
        public override double GetSimulationTime() {
            return m_CurrentPhystime;
        }

        bool CoupledIteration = true;

        int m_CoupledIterations = 0;

        int m_InnerCoupledIterations = 0;

        Control.NonLinearSolverConfig m_nonlinconfig;

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
            if(ResidualNorm >= this.m_nonlinconfig.ConvergenceCriterion) {
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


        //static double MatrixDist(MsrMatrix _A, MsrMatrix B) {
        //    MsrMatrix A = _A.CloneAs();
        //    A.Acc(-1.0, B);
        //    double w = A.InfNorm();
        //    return w;
        //}
        

        /// <summary>
        /// Perform temporal integration/implicit timestepping
        /// </summary>
        public override bool Solve(double phystime, double dt) {
            return Solve(phystime, dt, ComputeOnlyResidual: false);
        }


        /// <summary>
        /// Solver/Implicit Time Integrator
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
        /// <returns>
        /// - true: solver algorithm successfully converged
        /// - false: something went wrong
        /// </returns>
        public bool Solve(double phystime, double dt, bool ComputeOnlyResidual = false) {
            if(!initialized)
                SingleInit();

                        
            if (dt <= 0)
                throw new ArgumentOutOfRangeException();
            //if (m_CurrentDt_Timestep > 0 && Math.Abs(dt / m_CurrentDt_Timestep - 1.0) > 1.0e-14)
            //    throw new ArgumentOutOfRangeException();
            if (m_CurrentDt_Timestep > 0 && Math.Abs(m_CurrentDt_Timestep - dt) > 1e-14)
                AdaptToNewTimestep(dt, m_CurrentDt_Timestep);

            m_CurrentDt_Timestep = dt;

            bool success = true;
            //for (int i = 1; i <= incrementTimesteps; i++) {
            { 
                // push levelsets for every incremental timestep
                //if (i > 1)
                //    PushLevelSet();

                //// solve timestep with incremental timestep size
                //double incTimestepSize = dt / (double)incrementTimesteps;

                //success = success && Solve_Increment(i, phystime, incTimestepSize, ComputeOnlyResidual);
                success = success && Solve_Increment(phystime, dt, ComputeOnlyResidual);

                //phystime += incTimestepSize;
                phystime += dt;
            }

            


            return success;
        }

        private void AdaptToNewTimestep(double newTimestep, double oldTimestep) {
            if (Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative || Config_LevelSetHandling == LevelSetHandling.Coupled_Once)
                throw new NotImplementedException("Interpolation of mass_matrix_stack is not implemented");

            int timestepHistory = m_Stack_u.Length;
            CoordinateVector[] stuetzstelle = new CoordinateVector[m_Stack_u.Length];
            CoordinateVector[] newStackU = new CoordinateVector[m_Stack_u.Length];
            for (int i = 0; i < m_Stack_u.Length; i++)
            {
                stuetzstelle[i] = new CoordinateVector(m_Stack_u[i].Mapping);
                newStackU[i] = new CoordinateVector(m_Stack_u[i].Mapping);
                stuetzstelle[i].CopyEntries(m_Stack_u[i]);
                newStackU[i].CopyEntries(m_Stack_u[i]);
            }

            for (int i = 1; i < timestepHistory; i++) {
                double currentNewTimestep = -i * newTimestep;
                double[] langrangePoly = CalculateLangrangePolynom(currentNewTimestep, oldTimestep);
                newStackU[i].Scale(0);
                newStackU[i].CopyEntries(stuetzstelle[0]);
                newStackU[i].Scale(langrangePoly[0]);
                for (int j = 1; j < timestepHistory; j++) {
                    newStackU[i].Acc(langrangePoly[j], stuetzstelle[j]);
                }
            }
            m_Stack_u.Clear();
            m_Stack_u = newStackU.CloneAs();
        }

        private double[] CalculateLangrangePolynom(double time, double oldTimestep) {
            int timestepHistory = m_Stack_u.Length;
            double[] lPoly = new double[timestepHistory];
            for (int i = 0; i < timestepHistory; i++) {
                if(i != 0)
                    lPoly[i] = time / (-i * oldTimestep);
                for (int j = 1; j < timestepHistory; j++) {
                    if (j == i)
                        continue;
                    if (i == 0 && j == 1)
                        lPoly[i] = (time + j * oldTimestep) / ((j - i) * oldTimestep);
                    lPoly[i] *= (time + j * oldTimestep) / ((j - i) * oldTimestep);
                }
            }
            return lPoly;
        }

        double m_CurrentDt_Timestep = -1;

        bool OneTimeMgInit = false;

        double fsiOldPhystime = 0;




        /// <summary>
        /// Solver/Time Integrtor
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
        /// <returns>
        /// - true: solver algorithm successfully converged
        /// - false: something went wrong
        /// </returns>
        bool Solve_Increment(double phystime, double dt, bool ComputeOnlyResidual = false) {
            using (var tr = new FuncTrace()) {
                if (dt <= 0)
                    throw new ArgumentOutOfRangeException();
                //if (m_CurrentDt > 0 && Math.Abs(dt / m_CurrentDt - 1.0) > 1.0e-14)
                //    throw new ArgumentOutOfRangeException();

                m_CurrentPhystime = phystime;
                m_CurrentDt = dt;
                m_IterationCounter = 0;
                m_CoupledIterations = 0;
                m_InnerCoupledIterations = 0;
                PushStack();
                fsiOldPhystime = phystime;
                //if (incrementTimesteps == 1)
                //    dt = m_CurrentDt;
                //else
                //    Console.WriteLine("Increment solve, timestep #{0}, dt = {1} ...", increment, dt);
                dt = m_CurrentDt;


                // ===========================================
                // update level-set (in the case of splitting)
                // ===========================================
                if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                    || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                    || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {

                    Debug.Assert(m_CurrentAgglomeration == null);

                    double ls_dt = dt;
                    if (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting)
                        ls_dt *= 0.5;

                    // remember which old cells had values
                    //var oldCCM = this.UpdateCutCellMetrics();

                    // evolve the level set
                    if (Config_LevelSetHandling != LevelSetHandling.FSILieSplittingFullyCoupled && !this.coupledOperator) {
                        m_LsTrk.IncreaseHistoryLength(1);
                        m_LsTrk.PushStacks();
                    }

                    int oldPushCount = m_LsTrk.PushCount;
                    int oldVersion = m_LsTrk.VersionCnt;
                    this.MoveLevelSetAndRelatedStuff(m_Stack_u[0].Mapping.Fields.ToArray(), phystime, ls_dt, 1.0);

                    int newPushCount = m_LsTrk.PushCount;
                    int newVersion = m_LsTrk.VersionCnt;
                    if ((newPushCount - oldPushCount) != 0 && !coupledOperator)
                        throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");
                    if ((newVersion - oldVersion) != 1 && !coupledOperator)
                        throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");

                    // in the case of splitting, the fields must be extrapolated 
                    //var newCCM = this.UpdateCutCellMetrics();
                    //var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
                    Debug.Assert(m_LsTrk.HistoryLength >= 1);
                    var SplittingAgg = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                        __AgglomerationTreshold: 0.0, AgglomerateNewborn: true, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true,
                        oldTs__AgglomerationTreshold: new double[] { 0.0 }, Tag: "SolveIncrement");
                    for (int i = 0; i < this.m_Stack_u.Length; i++)
                        SplittingAgg.Extrapolate(this.m_Stack_u[i].Mapping);

                    // delete new agglomeration; in case of splitting, the agglomeration for the **bulk operator timestep** does not depend on previous time-steps
                    m_CurrentAgglomeration = null;
                }

                // ==============================================
                // solve main system
                // ==============================================
                bool success;

                int oldLsTrkPushCount = m_LsTrk.PushCount;

                {
                    int[] Jtot =
                        (new int[] { base.m_LsTrk.Regions.GetCutCellMask().NoOfItemsLocally, base.m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells })
                        .MPISum();
                    //Console.WriteLine("No of cells {0}, No of cut cells {1}.", Jtot[1], Jtot[0]);
                    if (Jtot[0] == Jtot[1])
                        throw new ArithmeticException("All cells are cut cells - check your settings!");
                }

                //{
                //    var rnd = new Random(1);
                //    var vec = new CoordinateVector(CurrentStateMapping);
                //    int L = CurrentStateMapping.LocalLength;
                //    for(int i = 0; i < L; i++)
                //        vec[i] = rnd.NextDouble();


                //    double[] Affine;
                //    BoSSS.Foundation.Quadrature.NonLin.Arsch.ShutTheFuckUp = true;
                //    this.AssembleMatrixCallback(out BlockMsrMatrix System, out Affine, out BlockMsrMatrix MaMa, CurrentStateMapping.Fields.ToArray(), false, out var dummy);
                //    Debug.Assert(System == null);

                //    base.Residuals.Clear();
                //    base.Residuals.SetV(Affine, -1.0);


                //}




                if (!ComputeOnlyResidual) {

                    // ++++++++++++++++++
                    // normal solver run 
                    // ++++++++++++++++++


                    if (RequiresNonlinearSolver) {
                        // ++++++++++++++++++++++++++++++++
                        // Nonlinear Solver (Navier-Stokes)
                        // ++++++++++++++++++++++++++++++++

                        // use solver
                        var nonlinSolver = GetNonlinSolver();
                        success = nonlinSolver.SolverDriver(m_Stack_u[0], default(double[])); // Note: the RHS is passed as the affine part via 'this.SolverCallback'

                        if (base.QueryHandler != null) {
                            base.QueryHandler.ValueQuery(QueryHandler.Conv, success ? 1.0 : 0.0, true);
                            base.QueryHandler.ValueQuery(QueryHandler.NonLinIter, nonlinSolver.NoOfNonlinearIter, true);
                            base.QueryHandler.ValueQuery(QueryHandler.NoOfCells, this.m_LsTrk.GridDat.CellPartitioning.TotalLength, true);
                            base.QueryHandler.ValueQuery(QueryHandler.DOFs, nonlinSolver.EssentialDOFs, true); // 'essential' DOF, in the XDG case less than cordinate mapping length 
                        }


                        // 'revert' agglomeration
                        Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                        m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

                    } else {
                        // ++++++++++++++++++++++++++++++++
                        // Linear Solver (Stokes)
                        // ++++++++++++++++++++++++++++++++
                        tr.Info("Using linear solver.");

                        // build the saddle-point matrix
                        BlockMsrMatrix System, MaMa;
                        double[] RHS;
                        this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true, out var dummy);
                        RHS.ScaleV(-1);

                        // update the multigrid operator
                        csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                        MultigridOperator mgOperator;
                        using (new BlockTrace("MultigridOperator setup", tr)) {
                            mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                                System, MaMa,
                                this.Config_MultigridOperator,
                                dummy);
                        }

                        using (var linearSolver = GetLinearSolver(mgOperator)) {



                            // try to solve the saddle-point system.
                            TimeSpan duration;
                            using (new BlockTrace("Solver_Run", tr)) {
                                var st = DateTime.Now;
                                mgOperator.UseSolver(linearSolver, m_Stack_u[0], RHS);
                                //mgOperator.ComputeResidual(this.Residuals, m_Stack_u[0], RHS);
                                duration = DateTime.Now - st;
                            }
                            Console.WriteLine("solver success: " + linearSolver.Converged + "; runtime: " + duration.TotalSeconds + " sec.");
                            success = linearSolver.Converged;


                            // 'revert' agglomeration
                            Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                            m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);


                            if (base.QueryHandler != null) {
                                base.QueryHandler.ValueQuery(QueryHandler.Conv, linearSolver.Converged ? 1.0 : 0.0, true);
                                base.QueryHandler.ValueQuery(QueryHandler.NoIter, linearSolver.ThisLevelIterations, true);
                                base.QueryHandler.ValueQuery(QueryHandler.NoOfCells, this.m_LsTrk.GridDat.CellPartitioning.TotalLength, true);
                                base.QueryHandler.ValueQuery(QueryHandler.DOFs, mgOperator.Mapping.TotalLength, true); // 'essential' DOF, in the XDG case less than cordinate mapping length 
                            }

                        }

                        //ExtractSomeSamplepoints("samples");
                    }



                } else {
                    // ++++++++++++++++++++++++++++++++++++
                    // compute residual of actual solution 
                    // ++++++++++++++++++++++++++++++++++++


                    double[] Affine;
                    this.AssembleMatrixCallback(out BlockMsrMatrix System, out Affine, out BlockMsrMatrix MaMa, CurrentStateMapping.Fields.ToArray(), false, out var dummy);
                    Debug.Assert(System == null);

                    base.Residuals.Clear();
                    base.Residuals.SetV(Affine, -1.0);

                    success = true;

#if DEBUG
                    {

                        this.AssembleMatrixCallback(out BlockMsrMatrix checkSystem, out double[] checkAffine, out BlockMsrMatrix MaMa1, CurrentStateMapping.Fields.ToArray(), true, out var dummy2);

                        double[] checkResidual = new double[checkAffine.Length];
                        checkResidual.SetV(checkAffine, -1.0);
                        checkSystem.SpMV(-1.0, m_Stack_u[0], +1.0, checkResidual);

                        Console.WriteLine("Norm of evaluated residual: " + base.Residuals.MPI_L2Norm());
                        Console.WriteLine("Norm of reference residual: " + checkResidual.MPI_L2Norm());


                        double distL2 = GenericBlas.L2DistPow2(checkResidual, base.Residuals).MPISum().Sqrt();
                        double refL2 = (new double[] { GenericBlas.L2NormPow2(m_Stack_u[0]), GenericBlas.L2NormPow2(checkResidual), GenericBlas.L2NormPow2(base.Residuals) }).MPISum().Max().Sqrt();

                        if (distL2 >= refL2 * 1.0e-5) {
                            double __distL2 = GenericBlas.L2DistPow2(checkAffine, base.Residuals).MPISum().Sqrt();
                        }

                        Tecplot.Tecplot.PlotFields(base.Residuals.Fields, "resi", 0.0, 2);

                        Assert.LessOrEqual(distL2, refL2 * 1.0e-5, "Significant difference between linearized and non-linear evaluation.");

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

                    if (m_ResLogger != null)
                        m_ResLogger.NextIteration(true);
                }

                bool calculateCondNumbers = false;
                
                if (calculateCondNumbers) {
                    var table = base.OperatorAnalysis(plotStencilCondNumViz: false, calculateStencils: false, calculateMassMatrix: true);
                    table.SaveToTextFileDebugUnsteady("CondEst", ".txt");
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
                        oldTs__AgglomerationTreshold: new double[] { 0.0 }, Tag: "SolveIncrement");
                    for (int i = 0; i < this.m_Stack_u.Length; i++)
                        SplittingAgg.Extrapolate(this.m_Stack_u[i].Mapping);

                    Debug.Assert(m_CurrentAgglomeration == null);
                }

                // ======
                // return 
                // ======
                //string path = Directory.GetCurrentDirectory();
                //var dinfo = Directory.CreateDirectory(path+@"\plots");
                //ExecuteWaterfallAnalysis(dinfo.FullName);
                //CreateFAMatrices(dinfo.FullName);


                m_CurrentPhystime = phystime + dt;

                if (Config_LevelSetHandling == LevelSetHandling.None) {
                    m_LsTrk.UpdateTracker(m_CurrentPhystime); // call is required to bring the internal time-stamp up-to-date;
                }

                return success;
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
            this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true, out var opi);
            RHS.ScaleV(-1);

            // update the multigrid operator
            MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                System, MaMa,
                this.Config_MultigridOperator,
                opi);            

            // set-up the convergence observer
            double[] uEx = this.m_Stack_u[0].ToArray();
            int L = uEx.Length;
            m_CurrentAgglomeration.ClearAgglomerated(uEx, this.m_Stack_u[0].Mapping);
            var CO = new ConvergenceObserver(mgOperator, MaMa, uEx);
            uEx = null;
            CO.TecplotOut = TecOutBaseName;
            //CO.PlotDecomposition(this.u.CoordinateVector.ToArray());

            // init linear solver
            using(var linearSolver = GetLinearSolver(mgOperator)) {
                ((ISolverWithCallback)linearSolver).IterationCallback = CO.IterationCallback;

                // try to solve the saddle-point system.
                mgOperator.UseSolver(linearSolver, new double[L], RHS);
            }

            // return
            return CO;
        }

        /*
        public void ExtractSomeSamplepoints(string OutputDir) {
            BlockMsrMatrix System, MaMa;
            double[] RHSsmall, RHSbig, RHS;
            double[] usmall, ubig;
            int Lsmall, Lbig;

            this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true, out var opi);
            RHS.ScaleV(-1);
            MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                System, MaMa,
                this.Config_MultigridOperator,
                opi.DomainVar.Select(varName => opi.FreeMeanValue[varName]).ToArray());

            Lbig = mgOperator.Mapping.ProblemMapping.LocalLength;
            Lsmall = mgOperator.Mapping.LocalLength;
            usmall = new double[Lsmall];
            ubig = new double[Lbig];
            RHSbig = RHS.CloneAs();
            RHSsmall = new double[Lsmall];

            var ReferenceSolver = new DirectSolver() {
                WhichSolver = DirectSolver._whichSolver.PARDISO
            };

            ReferenceSolver.Init(mgOperator);
            mgOperator.TransformRhsInto(RHSbig, RHSsmall, true);
            ReferenceSolver.Solve(usmall, RHSsmall);
            mgOperator.TransformSolFrom(ubig, usmall);

            m_CurrentAgglomeration.ClearAgglomerated(ubig, this.m_Stack_u[0].Mapping);
            var CO = new ConvergenceObserver(mgOperator, MaMa, ubig);
            CO.Resample(1, usmall, "samples");
        }
        */

        /*
        public void CreateFAMatrices(string OutputDir) {
            bool use_exact_solution = false; // switch between exact and zero solution

            // build the saddle-point matrix
            //AssembleMatrix(this.CurrentVel, dt, phystime + dt);
            BlockMsrMatrix System, MaMa;
            double[] RHSsmall, RHSbig, RHS;
            double[] usmall, ubig;
            double[] changeOfu;
            int Lsmall, Lbig;
            

            this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true, out var opi);
            RHS.ScaleV(-1);

            // update the multigrid operator
            MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                System, MaMa,
                this.Config_MultigridOperator,
                opi.DomainVar.Select(varName => opi.FreeMeanValue[varName]).ToArray());

            // Get Reference Solver
            Lbig = mgOperator.BaseGridProblemMapping.LocalLength;
            Lsmall = mgOperator.Mapping.LocalLength;

            if (!use_exact_solution) {
                RHS.Clear(); //zero solution
            }

            var StartSolution = GenericBlas.RandomVec(Lbig, 0);
            RHSbig = RHS.CloneAs();
            RHSsmall = new double[Lsmall];
            usmall = new double[Lsmall];
            ubig = new double[Lbig];

            if (use_exact_solution) {
                StartSolution = new double[Lbig];
                var ReferenceSolver = new DirectSolver() {
                    WhichSolver = DirectSolver._whichSolver.PARDISO
                };

                ReferenceSolver.Init(mgOperator);
                mgOperator.TransformRhsInto(RHSbig, RHSsmall, true);
                ReferenceSolver.Solve(usmall, RHSsmall);
                mgOperator.TransformSolFrom(ubig, usmall);
            }

            // Get the configured solvers
            ISolverSmootherTemplate linearSolver;
            NonlinearSolver NonlinearSolver;
            GetSolver(out NonlinearSolver, out linearSolver);

            // set-up the convergence observer
            m_CurrentAgglomeration.ClearAgglomerated(ubig, this.m_Stack_u[0].Mapping);
            var CO = new ConvergenceObserver(mgOperator, MaMa, ubig);
            CO.TecplotOut = OutputDir;
            changeOfu = new double[Lsmall];
            ((OrthonormalizationMultigrid)linearSolver).ExtractSamples = delegate (int iter, double[] u, string name){
                //if (iter % 5 != 0 && iter != 1)
                //    return;
                Debug.Assert(u.Length== usmall.Length);
                changeOfu.SetV(u);
                changeOfu.AccV(-1.0, usmall);
                CO.Resample(iter, changeOfu, name);
            };
            
            // init linear solver
            linearSolver.Init(mgOperator);

            // try to solve the saddle-point system.
            mgOperator.UseSolver(linearSolver, StartSolution, RHS);
        }
        */

        /// <summary>
        /// If an iterative linear solver is used:
        /// performs a waterfall analysis, <see cref="ConvergenceObserver.WaterfallAnalysis"/>,
        /// i.e. visualizes the decay of DG modes on different grid levels
        /// over the solver iterations
        /// </summary>
        /// <param name="OutputDir"></param>
        public void ExecuteWaterfallAnalysis(string OutputDir) {
            BlockMsrMatrix System, MaMa;
            double[] RHS;
            this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true, out var opi);
            RHS.ScaleV(-1);

            // update the multigrid operator
            MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                System, MaMa,
                this.Config_MultigridOperator,
                opi);

            if(LinearSolverConfig is IterativeSolverConfig ics)
                ics.MaxSolverIterations = 30;


            using(ISolverSmootherTemplate linearSolver = GetLinearSolver(mgOperator)) {

                var plots = ConvergenceObserver.WaterfallAnalysis((ISolverWithCallback)linearSolver, mgOperator, MaMa);
                // put this shit out
                foreach(var kv in plots) {
                    //var CL = kv.Value.ToGnuplot().PlotCairolatex(xSize: 14, ySize: 12);
                    //CL.WriteMinimalCompileableExample(Path.Combine(OutputDir, "plot_" + kv.Key + ".tex"), kv.Key + ".tex");
                    kv.Value.SavePgfplotsFile_WA(OutputDir + @"\" + kv.Key + ".tex");
                }
            }
        }

        
        /// <summary>
        /// Executes the linear solver with a random right-hand-side (RHS)
        /// </summary>
        /// <param name="Iter">on exit, the number of solver iterations used</param>
        public void ExecuteRandom(out int Iter) {
            // Get the configured solvers
            

            BlockMsrMatrix System, MaMa;
            double[] RHS;
            this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true, out var opi);
            RHS.ScaleV(-1);

            // update the multigrid operator
            MultigridOperator mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
                System, MaMa,
                this.Config_MultigridOperator,
                opi);

            int L = RHS.Length;
            var x0 = new double[L];

            // use a random init for intial guess.
            Random rnd = new Random();
            RHS = new double[L];
            for (int l = 0; l < L; l++) {
                RHS[l] = rnd.NextDouble();
            }

            using(var linearSolver = GetLinearSolver(mgOperator)) {
                // init linear solver
                linearSolver.Init(mgOperator);

                // try to solve the saddle-point system.
                mgOperator.UseSolver(linearSolver, x0, RHS);
                Iter = linearSolver.ThisLevelIterations;
            }
        }
        

        /// <summary>
        /// Bad design, should be removed; FK 18feb22
        /// </summary>
        public event Action<int, double[], double[], MultigridOperator> CustomIterationCallback;

        

        protected override NonlinearSolver GetNonlinSolver() {
            var nonlinSolver = base.GetNonlinSolver();

            if(this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                nonlinSolver.IterationCallback += this.CoupledIterationCallback;
                if(nonlinSolver is FixpointIterator)
                    ((FixpointIterator)nonlinSolver).Iteration_Count = this.CoupledIterationCounter;
            }

            nonlinSolver.IterationCallback += this.CustomIterationCallback;

            return nonlinSolver;
        }

        protected override ISolverSmootherTemplate GetLinearSolver(MultigridOperator op) {
            

            var linearSolver = base.GetLinearSolver(op);
            if(linearSolver is ISolverWithCallback swc) {
                swc.IterationCallback += this.CustomIterationCallback;
            }
            return linearSolver;
        }


        /// <summary>
        /// Bad, undocumented design! To be removed! Fk, 18feb22
        /// </summary>
        private bool coupledOperator = false;


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

            if(!ilPSP.DoubleExtensions.ApproxEqual(m_LsTrk.Regions.Time, PhysTime))
                throw new ApplicationException($"Before Level-Set update, mismatch in time between tracker (Regions.Time = {m_LsTrk.Regions.Time}) and physical time ({PhysTime}).");



            m_LastLevelSetResidual = this.UpdateLevelset().Update(locCurSt, PhysTime, dt, UnderRelax, (this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting));


            int newVersion = m_LsTrk.VersionCnt;
            int newPushCount = m_LsTrk.PushCount;


            if ((newVersion - oldVersion) != 1 && !coupledOperator)
                throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");
            if((newVersion - oldVersion) != 0 && coupledOperator)
                throw new ApplicationException("Expecting exactly no call to 'UpdateTracker(...)' in 'UpdateLevelset(...)' for coupled Operators.");
            if ((newPushCount - oldPushCount) != 0)
                throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");

            if(!ilPSP.DoubleExtensions.ApproxEqual(m_LsTrk.Regions.Time, PhysTime + dt))
                throw new ApplicationException($"After Level-Set update, mismatch in time between tracker (Regions.Time = {m_LsTrk.Regions.Time}) and physical time ({PhysTime + dt}).");
            

            //// new cut-cell metric
            //m_Stack_CutCellMetrics[0] = this.UpdateCutCellMetrics();
            //if (!m_Stack_CutCellMetrics[0].SpeciesList.SetEquals(Config_MassScale.Keys))
            //    throw new ApplicationException("Mismatch between species lists.");


            // re-sort mass matrices 
            {
                var Mtx2Update = m_Stack_MassMatrix.Skip(1).ToArray(); // shallow copy!
                TimeSteppingUtils.OperatorLevelSetUpdate(m_LsTrk, Mtx2Update, CurrentStateMapping, CurrentStateMapping);
                for (int i = 1; i < m_Stack_MassMatrix.Length; i++) {
                    m_Stack_MassMatrix[i] = Mtx2Update[i - 1];
                }
            }


            /*
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
            */

        }
    }

}
