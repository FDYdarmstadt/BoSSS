using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LoadBalancing;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// Not intended for direct usage, use <see cref="XdgApplicationWithSollver{T}"/> or <see cref="DgApplicationWithSollver{T}"/> instead.
    /// Base-class for applications with a monolithic operator and a time-integrator.
    /// </summary>
    abstract public class ApplicationWithSolver<T> : Application<T>
        where T : AppControlSolver, new() {



        /// <summary>
        /// Mapping for fields defined by <see cref="InstantiateSolutionFields"/>
        /// </summary>
        public virtual CoordinateMapping CurrentState {
            get {
                return CurrentStateVector.Mapping;
            }
           
        }
        CoordinateVector m_CurrentStateVector;

        /// <summary>
        /// DG coordinates of <see cref="CurrentState"/> in a single vector
        /// </summary>
        public virtual CoordinateVector CurrentStateVector {
            get {
                if(!object.ReferenceEquals(m_CurrentStateVector.Mapping.GridDat, this.GridData))
                    throw new ApplicationException("Grid data object mismatch - something gone terribly wrong, maybe during mesh adaptation.");

                return m_CurrentStateVector;
            }
            protected set {
                if (value != null) {
                    if (!object.ReferenceEquals(value.Mapping.GridDat, this.GridData))
                        throw new ArgumentException("Grid data object mismatch.");
                }
                m_CurrentStateVector = value;
            }
        }

        List<DGField> m_Parameters;

        /// <summary>
        /// List for parameter fields can be set by <see cref="CreateAdditionalFields"/>
        /// </summary>
        public virtual IList<DGField> Parameters {
            get {
                if(m_Parameters != null) {
                    foreach(var p in m_Parameters) {
                        if( p != null) {
                            if(!object.ReferenceEquals(p.GridDat, this.GridData))
                                throw new ApplicationException("Grid data object mismatch - something gone terribly wrong, maybe during mesh adaptation.");
                        }
                    }
                }

                return m_Parameters;
            }
            protected set {
                if(value != null) {
                    foreach(var p in value) {
                        if( p != null) {
                            if(!object.ReferenceEquals(p.GridDat, this.GridData))
                                throw new ArgumentException("Grid data object mismatch.");
                        }
                    }
                }


                if(m_Parameters == null)
                    m_Parameters = new List<DGField>();
                else
                    m_Parameters.Clear();
                m_Parameters.AddRange(value);
            }

        }

        /// <summary>
        /// Mapping for fields defined by <see cref="InstantiateResidualFields"/>
        /// </summary>
        public virtual CoordinateMapping CurrentResidual {
            get {
                return CurrentResidualVector.Mapping;
            }
        }

        CoordinateVector m_CurrentResidualVector;

        /// <summary>
        /// DG coordinates of <see cref="CurrentResidual"/> in a single vector
        /// </summary>
        public virtual CoordinateVector CurrentResidualVector {
            get {
                if(!object.ReferenceEquals(m_CurrentResidualVector.Mapping.GridDat, this.GridData))
                    throw new ApplicationException("Grid data object mismatch - something gone terribly wrong, maybe during mesh adaptation.");

                return m_CurrentResidualVector;
            }
            protected set {
                if (value != null) {
                    if (!object.ReferenceEquals(value.Mapping.GridDat, this.GridData))
                        throw new ArgumentException("Grid data object mismatch.");
                }
                m_CurrentResidualVector = value;
            }
        }

        /// <summary>
        /// makes direct use of <see cref="XdgTimesteppingBase.OperatorAnalysis"/>; aids the condition number scaling analysis
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis() {
            return this.Timestepping.OperatorAnalysis();
        }      

        //abstract protected void CreateTrackerHack();

        /// <summary>
        /// Called on startup and 
        /// also called after grid adaptation/MPI redistribution
        /// </summary>
        protected override void CreateFields() {

            // Before creating fields make sure the operator is cleared, so that the equations are instantiated with the correct LsTrk and stuff...
            // This is not that clear: This functionality should theoretically belong to <see cref="CreateEquationsAndSolvers(GridUpdateDataVaultBase)"/>,
            // However the Operator and its equations are instantiated on the first call of <see cref="Operator"/> when it is empty. 
            // These occur in this method e.g. in <see cref="InstantiateSolutionFields"/> in the derived classes <see cref="DgApplicationWithSolver{T}"/> and <see cref="XdgApplicationWithSolver{T}"/>
            // So nulling the Operator here ensures, that it is created again with the updated LsTrk and stuff...
            ClearOperator();

            base.CreateFields();

            if(this.LsTrk==null) CreateTracker();

            // solution:
            var solFields = InstantiateSolutionFields();
            CurrentStateVector = new CoordinateVector(solFields);
            foreach(var f in solFields) {
                base.RegisterField(f);
            }

            // residuals:
            var resFields = InstantiateResidualFields();
            CurrentResidualVector = new CoordinateVector(resFields);
            foreach(var f in resFields) {
                base.RegisterField(f, IOListOption.Always);
            }

            CreateAdditionalFields();
        }

        /// <summary>
        /// intended for the initialization of additional DG fields, which are not either solution or residuals;
        /// An example would be e.g. parameters such as Gravity which are set through the control file.
        /// </summary>
        protected virtual void CreateAdditionalFields() {

            // parameters
            var parameterFields = Operator.InvokeParameterFactory(CurrentState.Fields);
            Parameters = new List<DGField>();
            // Console.WriteLine(parameterFields.Count()); [Toprak]: I noticed that this is unnecessary, therefore I made it comment.

            foreach (var f in parameterFields) {
                this.Parameters.Add(f);
                if(f != null) 
                    // in some solvers, NULL parameters are present,
                    // i.e. the parameter is specified by the operator,
                    // but it is not needed in any equation component and therefore set to null.
                    base.RegisterField(f);
            }            

        }



        /// <summary>
        /// Override to instantiate the dependent variables of the PDE to solve;
        /// </summary>
        protected abstract IEnumerable<DGField> InstantiateSolutionFields();

        

        /// <summary>
        /// Override to instantiate fields to store the solver residuals
        /// </summary>
        public virtual IEnumerable<DGField> InstantiateResidualFields() {
            var SolFields = this.CurrentState.Fields;

            var CodNames = this.Operator.CodomainVar;
            if(CodNames.Count != SolFields.Count)
                throw new NotSupportedException("Unable to instantiate residual fields automatically - number of codomain variables is different than number of solution fields.");

            DGField[] R = new DGField[SolFields.Count];
            for(int i = 0; i < R.Length; i++) {
                R[i] = SolFields[i].CloneAs();
                R[i].Identification = "Residual-" + CodNames[i];
                R[i].Clear();
            }

            return R;
        }


        /// <summary>
        /// Main spatial operator
        /// </summary>
        public abstract IDifferentialOperator Operator {
            get;
        }


        /// <summary>
        /// initialization of <see cref="ApplicationWithSolver{T}.Timestepping"/>
        /// </summary>
        protected abstract void InitSolver();


        /// <summary>
        /// The time-integrator/solver
        /// </summary>
        public XdgTimestepping Timestepping {
            get;
            protected set;
        }


        // <summary>
        /// Number of time-steps required for restart, e.g. 1 for Runge-Kutta and implicit/explicit Euler, 2 for BDF2, etc.
        /// </summary>
        protected override int BurstSave {
            get {
                int Timestepping_bs;
                if(Timestepping != null) {
                    Timestepping_bs = Timestepping.BurstSave;
                } else {
                    string schStr = Control.TimeSteppingScheme.ToString().ToLower();
                    if(schStr.StartsWith("bdf")) {
                        Timestepping_bs = Convert.ToInt32(schStr.Substring(3));
                    } else {
                        Timestepping_bs = 1;
                    }
                }

                return Math.Max(Timestepping_bs, this.Control.BurstSave);
            }
        }


        /// <summary>
        /// call the solver
        /// </summary>
        public bool Solve(double phystime, double dt) {
            return this.Timestepping.Solve(phystime, dt);
        }

        /// <summary>
        /// setup AMR level indicators
        /// </summary>
        protected override void SetUpEnvironment() {
            base.SetUpEnvironment();

            // AMR level indactors
            //====================
            m_AMRLevelIndicators.Clear();
            if (this.Control != null && this.Control.activeAMRlevelIndicators != null) {
                m_AMRLevelIndicators.AddRange(this.Control.activeAMRlevelIndicators);
            }

            foreach (var lvlInd in ActiveAMRLevelIndicators) {
                lvlInd.Setup(this);
            }
        }


        /// <summary>
        /// 
        /// </summary>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {

            if (L == null) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                // Creation of time-integrator (initial, no balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                
                InitSolver();
                Timestepping.RegisterResidualLogger(this.ResLogger);
                if (Timestepping.m_BDF_Timestepper != null)
                    Timestepping.m_BDF_Timestepper.SingleInit();


            } else {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // restore time-stepper after grid redistribution (dynamic load balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                // currently, only supported for the BDF timestepper.
                
                Timestepping.DataRestoreAfterBalancing(L, CurrentState.Fields, CurrentResidual.Fields, base.LsTrk, base.MultigridSequence, this.Operator);

                if (!object.ReferenceEquals(this.Operator, Timestepping.Operator))
                    throw new ApplicationException();

            }
        }


        /// <summary>
        /// Step 2 of 2 for dynamic load balancing: restore this objects 
        /// status after the grid has been re-distributed.
        /// </summary>
        public override void DataBackupBeforeBalancing(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            if (Timestepping != null)   // in case of startUp backup
                Timestepping.DataBackupBeforeBalancing(L);
            CurrentStateVector = null;
            CurrentResidualVector = null;
            ClearOperator();
        }


        /// <summary>
        /// 
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            using (var tr = new FuncTrace()) {

                if (this.Control.AdaptiveMeshRefinement) {

                    // Check grid changes
                    // ==================

                    int[] desiredLevels = GetDesiredRefinementLevels();

                    GridRefinementController gridRefinementController = new GridRefinementController((GridData)this.GridData, null);
                    bool AnyChange = gridRefinementController.ComputeGridChange(desiredLevels, out List<int> CellsToRefineList, out List<int[]> Coarsening);

                    int NoOfCellsToRefine = 0;
                    int NoOfCellsToCoarsen = 0;
                    if (AnyChange.MPIOr()) {
                        int[] glb = (new int[] { CellsToRefineList.Count, Coarsening.Sum(L => L.Length) }).MPISum();
                        NoOfCellsToRefine = glb[0];
                        NoOfCellsToCoarsen = glb[1];
                    }
                    long oldJ = this.GridData.CellPartitioning.TotalLength;

                    // Update Grid
                    // ===========
                    if (AnyChange.MPIOr()) {

                        Console.WriteLine(" Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                        Console.WriteLine(" Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                        newGrid = ((GridData)this.GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                    } else {

                        newGrid = null;
                        old2NewGrid = null;
                    }

                } else {

                    newGrid = null;
                    old2NewGrid = null;
                }

            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private int[] GetDesiredRefinementLevels() {

            if(Control.AdaptiveMeshRefinement == true)
                if(ActiveAMRLevelIndicators == null || ActiveAMRLevelIndicators.Count <= 0) {
                    Console.Error.WriteLine("Control object configuration inconsistent: 'AdaptiveMeshRefinement == true', but no refinement indicators in 'activeAMRLevelIndicators' are set.");
                }

            int J = this.GridData.CellPartitioning.LocalLength;
            int[] levelChanges = new int[J];

            // combine all results of active level indicators
            foreach (var lvlInd in ActiveAMRLevelIndicators) {
                if(ActiveAMRLevelIndicators.IndexOf(lvlInd) == 0) {
                    levelChanges = lvlInd.DesiredCellChanges(); // levelChanges is instantiated to zero. Without this line, coarsening is impossible due to Max(a,b)
                }
                int[] lvls = lvlInd.DesiredCellChanges();
                //levelChanges = levelChanges.Zip(lvls, (a, b) => a + b).ToArray();
                levelChanges = levelChanges.Zip(lvls, (a, b) => Math.Max(a,b)).ToArray(); // keep finer level indicator, but don't double refine
            }

            // get desired level 
            int[] levels = new int[J];
            Cell[] cells = ((GridData)this.GridData).Grid.Cells;
            for (int j = 0; j < J; j++) {
                levels[j] = cells[j].RefinementLevel;
                if (levelChanges[j] > 0)
                    levels[j] += 1;
                else if (levelChanges[j] < 0)
                    levels[j] -= 1;
            }

            return levels;
        }


        List<AMRLevelIndicator> m_AMRLevelIndicators = new List<AMRLevelIndicator>();

        /// <summary>
        /// <see cref="Control.AppControl.activeAMRlevelIndicators"/>
        /// </summary>
        public IList<AMRLevelIndicator> ActiveAMRLevelIndicators {
            get {
                return m_AMRLevelIndicators;
            }
        }



        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>;
        /// temporary, to be moved into the spatial operator soon
        /// </summary>
        protected virtual MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int[] Degrees = this.CurrentState.BasisS.Select(b => b.Degree).ToArray();

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration entry will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[Degrees.Length];

                    // configurations for velocity
                    for (int d = 0; d < Degrees.Length; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Degrees[d] },
                            mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                            VarIndex = new int[] { d }
                        };
                    }
                }


                return configs;
            }

        }

        /// <summary>
        /// sets <see cref="Operator"/> to null.
        /// </summary>
        internal abstract void ClearOperator();

        /// <summary>
        /// Plot using Tecplot
        /// </summary>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            using (new FuncTrace()) {
                Tecplot.Tecplot.PlotFields(this.m_RegisteredFields, this.GetType().Name.Split('`').First() + "-" + timestepNo, physTime, superSampling);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="TimestepNo"></param>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetTimestep();
            Timestepping.Solve(phystime, dt);
            return dt;
        }

    }

    /// <summary>
    /// Base-class for XDG applications with a monolithic operator and a time-integrator.
    /// </summary>
    abstract public class XdgApplicationWithSolver<T> : ApplicationWithSolver<T>
         where T : AppControlSolver, new() //
    {

        /// <summary>
        /// Instantiation of fields according to the domain variable names in the spatial operator.
        /// </summary>
        protected override IEnumerable<DGField> InstantiateSolutionFields() {
            var DomNames = this.Operator.DomainVar;
            var ret = new DGField[DomNames.Count];
            for (int i = 0; i < DomNames.Count; i++) {
                string Name = DomNames[i];

                var fopts = this.Control.FieldOptions.Where(kv => kv.Key.WildcardMatch(Name)).SingleOrDefault().Value;
                if (fopts.Degree < 0) {
                    throw new ApplicationException($"Missing specification of DG degree for field {Name} in control object.");
                }

                var xb = new XDGBasis(this.LsTrk, fopts.Degree);
                ret[i] = new XDGField(xb, Name);
            }
            return ret;
        }

        /*
        ///// <summary>
        ///// Callback-template for level-set updates.
        ///// </summary>
        ///// <param name = "CurrentState" ></ param >
        ///// < param name= "time" >
        ///// Actual simulation time for the known value;
        ///// </param>
        ///// <param name = "dt" >
        ///// Timestep size.
        ///// </param>
        ///// <param name = "UnderRelax" >
        ///// </ param >
        ///// < param name="incremental">
        ///// true for Splitting schemes with subdivided level-set evolution(e.g.Strang-Splitting)
        ///// </param>
        ///// <returns>
        ///// Some kind of level-set-residual in order to check convergence in a fully coupled simulation
        ///// (see<see cref="LevelSetHandling.Coupled_Iterative"/>)
        ///// </returns>
        virtual public double UpdateLevelset(DGField[] CurrentState, double time, double dt, double UnderRelax, bool incremental) {
            LsTrk.UpdateTracker(time + dt);
            return 0.0;
        }
        */

        /// <summary>
        /// timestepper to update the all involved level-sets (typically one);
        /// can be null in case of a temporally constant level set.
        /// </summary>
        public virtual ISlaveTimeIntegrator GetLevelSetUpdater() {
            return null;
        }



        /// <summary>
        /// Instantiation of the spatial operator; cached in <see cref="XOperator"/>
        /// </summary>
        /// <returns></returns>
        /// <param name="D">spatial dimension</param>
        abstract protected XDifferentialOperatorMk2 GetOperatorInstance(int D);

        /// <summary>
        /// 
        /// </summary>
        public override IDifferentialOperator Operator {
            get {
                return XOperator;
            }
        }

        /// <summary>
        /// Instantiation of Level-Sets (<see cref="ILevelSet"/>, <see cref="LevelSet"/>) and the tracker.
        /// </summary>
        protected abstract LevelSetTracker InstantiateTracker();


        /// <summary>
        /// empty in the DG case
        /// </summary>
        protected override void CreateTracker() {
            var trk = InstantiateTracker();
            //var test = this.Operator;
            if (base.LsTrk == null) {
                base.LsTrk = trk;
            } else {
                if (!object.ReferenceEquals(trk, base.LsTrk))
                    throw new ApplicationException("It seems there is more then one Level-Set-Tracker in the application; not supported by the Application class.");
            }

            //foreach(var ls in LsTrk.LevelSetHistories.Select(history => history.Current)) {
            //    if(ls is DGField f) {
            //        base.RegisterField(f);
            //    }
            //}

        }

        private XDifferentialOperatorMk2 m_XOperator {
            get;
            set;
        }

        /// <summary>
        /// Cache for <see cref="GetOperatorInstance"/>
        /// </summary>
        virtual public XDifferentialOperatorMk2 XOperator {
            get {
                if (m_XOperator == null) {
                    m_XOperator = GetOperatorInstance(this.Grid.SpatialDimension);
                    if (!m_XOperator.IsCommitted)
                        throw new ApplicationException("Operator must be committed by user.");
                }
                return m_XOperator;
            }
        }



        protected abstract LevelSetHandling LevelSetHandling {
            get;
        }

        /// <summary>
        /// instantiation of <see cref="ApplicationWithSolver{T}.Timestepping"/>
        /// </summary>
        protected override void InitSolver() {
            if (base.Timestepping != null)
                return;

            XdgTimestepping solver = new XdgTimestepping(
                XOperator,
                CurrentState.Fields,
                CurrentResidual.Fields,
                Control.TimeSteppingScheme,
                this.GetLevelSetUpdater,
                LevelSetHandling,
                MultigridOperatorConfig,
                MultigridSequence,
                Control.AgglomerationThreshold,
                Control.LinearSolver, Control.NonLinearSolver,
                this.LsTrk,
                Parameters,
                this.QueryHandler);
            base.Timestepping = solver;
            Timestepping.RegisterResidualLogger(this.ResLogger);
            Timestepping.TimesteppingBase.Config_LevelSetConvergenceCriterion = Control.LevelSet_ConvergenceCriterion;
            if (!object.ReferenceEquals(base.LsTrk, solver.LsTrk))
                throw new ApplicationException();

        }

        /// <summary>
        /// releases the operator
        /// </summary>
        internal override void ClearOperator() {
            m_XOperator = null;
        }

        bool PlotShadowfields = false;
        /// <summary>
        /// Plot using Tecplot
        /// </summary>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            using (new FuncTrace()) {
                // Cells Numbers - Local
                var CellNumbers = this.m_RegisteredFields.Where(s => s.Identification == "CellNumbers").SingleOrDefault();
                if (CellNumbers == null) {
                    CellNumbers = new SinglePhaseField(new Basis(this.GridData, 0), "CellNumbers");
                    this.RegisterField(CellNumbers);
                }
                CellNumbers.Clear();
                CellNumbers.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = j0 + j;
                        }
                    }
                }, new CellQuadratureScheme());

                // Cells Numbers - Global
                var CellNumbersGlob = this.m_RegisteredFields.Where(s => s.Identification == "CellNumbersGlobal").SingleOrDefault();
                if (CellNumbersGlob == null) {
                    CellNumbersGlob = new SinglePhaseField(new Basis(this.GridData, 0), "CellNumbersGlobal");
                    this.RegisterField(CellNumbersGlob);
                }
                CellNumbersGlob.Clear();
                CellNumbersGlob.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    for (int j = 0; j < Len; j++) {
                        long gidx = this.GridData.CellPartitioning.i0 + j0 + j;
                        for (int k = 0; k < K; k++) {
                            result[j, k] = gidx;
                        }
                    }
                }, new CellQuadratureScheme());


                // GlobalID
                var GlobalID = this.m_RegisteredFields.Where(s => s.Identification == "GlobalID").SingleOrDefault();
                if (GlobalID == null) {
                    GlobalID = new SinglePhaseField(new Basis(this.GridData, 0), "GlobalID");
                    this.RegisterField(GlobalID);
                }
                GlobalID.Clear();
                GlobalID.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    for (int j = 0; j < Len; j++) {
                        long gid = this.GridData.iLogicalCells.GetGlobalID(j + j0); ;
                        for (int k = 0; k < K; k++) {
                            result[j, k] = gid;
                        }
                    }
                }, new CellQuadratureScheme());

                // MPI_rank
                int my_rank = ilPSP.Environment.MPIEnv.MPI_Rank;
                var MPI_rank = this.m_RegisteredFields.Where(s => s.Identification == "MPI_rank").SingleOrDefault();
                if (MPI_rank == null) {
                    MPI_rank = new SinglePhaseField(new Basis(this.GridData, 0), "MPI_rank");
                    this.RegisterField(MPI_rank);
                }
                MPI_rank.Clear();
                MPI_rank.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = my_rank;
                        }
                    }
                }, new CellQuadratureScheme());


                // CutCell
                var XNSE_classifier = new CutStateClassifier();
                var classifiedCells = XNSE_classifier.ClassifyCells(this);

                var cutCellClass = this.m_RegisteredFields.Where(s => s.Identification == "cutCellClass").SingleOrDefault();
                if (cutCellClass == null) {
                    cutCellClass = new SinglePhaseField(new Basis(this.GridData, 0), "cutCellClass");
                    this.RegisterField(cutCellClass);
                }
                cutCellClass.Clear();

                int J = cutCellClass.Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                for (int j = 0; j < J; j++) {
                    cutCellClass.SetMeanValue(j, (double)classifiedCells[j]);
                }

                // LvSetDist
                var LvSetDist = this.m_RegisteredFields.Where(s => s.Identification == "LvSetDist").SingleOrDefault();
                if (LvSetDist == null) {
                    LvSetDist = new SinglePhaseField(new Basis(this.GridData, 0), "LvSetDist");
                    this.RegisterField(LvSetDist);
                }
                LvSetDist.Clear();
                LvSetDist.AccLevelSetDist(1, this.LsTrk, 1);


                if (PlotShadowfields) {
                    List<DGField> Fields2Plot = new List<DGField>();
                    foreach (var field in this.m_RegisteredFields) {
                        if (field is XDGField xField) {
                            foreach (var spc in xField.Basis.Tracker.SpeciesNames) {
                                Fields2Plot.Add(xField.GetSpeciesShadowField(spc));
                            }
                        } else {
                            Fields2Plot.Add(field);
                        }
                    }
                    Tecplot.Tecplot.PlotFields(Fields2Plot, this.GetType().Name.Split('`').First() + "-" + timestepNo, physTime, superSampling);
                } else {
                    base.PlotCurrentState(physTime, timestepNo, superSampling);
                }
            }
        }

    }


    /// <summary>
    /// Base-class for XDG applications with a monolithic operator and a time-integrator.
    /// </summary>
    abstract public class DgApplicationWithSolver<T> : ApplicationWithSolver<T>
         where T : AppControlSolver, new() {


        /// <summary>
        /// Instantiation of fields according to the domain variable names in the spatial operator.
        /// </summary>
        protected override IEnumerable<DGField> InstantiateSolutionFields() {
            var DomNames = this.Operator.DomainVar;
            var ret = new DGField[DomNames.Count];
            for(int i = 0; i < DomNames.Count; i++) {
                string Name = DomNames[i];

                var fopts = this.Control.FieldOptions.Where(kv => kv.Key.WildcardMatch(Name)).SingleOrDefault().Value;
                if(fopts.Degree < 0) {
                    throw new ApplicationException($"Missing specification of DG degree for field {Name} in control object.");
                }

                var b = new Basis(this.GridData, fopts.Degree);
                ret[i] = new SinglePhaseField(b, Name);
            }
            return ret;
        }


        /// <summary>
        /// initialization of the main spatial operator
        /// </summary>
        /// <param name="D">spatial dimension</param>
        abstract protected DifferentialOperator GetOperatorInstance(int D);

        /// <summary>
        /// empty in the DG case
        /// </summary>
        protected override void CreateTracker() {
               
        }

        /// <summary>
        /// Contains hack 
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();
            if(Timestepping != null && Timestepping.LsTrk == null) {
                // re-create the internal Dummy-Tracker after dynamic load-balancing/mesh refinement
                Timestepping.RecreateDummyTracker(this.GridData);
                this.LsTrk = Timestepping.LsTrk;
            }
        }


        /// <summary>
        /// 
        /// </summary>
        public override IDifferentialOperator Operator {
            get {
                return SOperator;
            }
        }


        DifferentialOperator m_SOperator;

        /// <summary>
        /// Cache for <see cref="GetOperatorInstance"/>
        /// </summary>
        virtual public DifferentialOperator SOperator {
            get {
                if(m_SOperator == null) {
                    m_SOperator = GetOperatorInstance(this.Grid.SpatialDimension);
                    if(!m_SOperator.IsCommitted)
                        throw new ApplicationException("Operator must be comitted by user.");

                }
                return m_SOperator;
            }
        }

        /// <summary>
        /// releases the operator
        /// </summary>
        internal override void ClearOperator() {
            m_SOperator = null;
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void InitSolver() {
            if (base.Timestepping != null)
                return;

            XdgTimestepping solver = new XdgTimestepping(
                SOperator,
                CurrentState.Fields,
                CurrentResidual.Fields,
                Control.TimeSteppingScheme,
                MultigridOperatorConfig,
                MultigridSequence,
                Control.LinearSolver, Control.NonLinearSolver,
                Parameters,
                this.QueryHandler);

            LsTrk = solver.LsTrk; // register the dummy tracker which the solver created internally for the DG case
            base.Timestepping = solver;
            Timestepping.RegisterResidualLogger(this.ResLogger);
        }

    }

}
