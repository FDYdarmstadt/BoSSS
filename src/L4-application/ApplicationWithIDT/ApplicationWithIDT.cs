using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution.LevelSetTools;
using System.Collections;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.Quadrature;
using ApplicationWithIDT.OptiLevelSets;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Statistic;
using BoSSS.Solution.Gnuplot;
using BoSSS.Application.BoSSSpad;
using CNS.Source;
using System.Threading.Tasks;
using static BoSSS.Solution.GridImport.NASTRAN.NastranFile;
using static System.CommandLine.Help.DefaultHelpText;

namespace ApplicationWithIDT {
    /// <summary>
    /// Abstract class for implicit discontinuity tracking (IDT) via optimization
    /// - Solves optimization problem   f(U,phi) -> min!
    ///                                 subject to r(U,phi)=0
    /// Here:
    ///     - f is an objective function [Code: obj_f, obj_vec ...] 
    ///     - r is the discretized weak form of some stationary/space-time conservation law $\nabla \mathbf{f}(U)=0$  [Code: Residual,...] 
    ///     - U is the state vector [Code: ConservativeFields]
    ///     - phi a Level Set (only one Level Set to be optimized supported) [Code:LsTBO, OptiLevelSet]
    /// 
    /// Different types of Level Set Parametrization are supported (may not work for all examples) and defined in IOptiLevelSet.cs
    /// - SplineOptiLevelSet
    /// - GlobalOptiLevelSet
    /// - SpecFemOptiLevelSet
    /// - SinglePhaseFieldOptiLevekSet
    /// 
    /// Different Objective Functions are supported
    /// - Rankine Hugoniot (AllFaces or only Faces in CutCells)
    /// - Enriched Residual (AllFaces or only Faces in CutCells)
    /// 
    /// Main Method is RunSolverOneStep
    /// 
    /// TODO:
    /// - Remove legacy code
    /// - Introduce as many default parameters as possible in IDTControl
    /// - Clean up
    /// - More Documentation
    /// </summary>
    public abstract class ApplicationWithIDT<T> : Application<T>
        where T : IDTControl, new() {
        #region Member XDG Fields
        public XDGField[] ConservativeFields { get; set; }
        public XDGField[] ConservativeFieldsReinitialized { get; set; } // only used when in debug
        public XDGField personField { get; set; }
        public XDGField[] originalStep { get; set; }
        public XDGField[] LagMult { get; set; }
        public XDGField[] obj_f_fields { get; set; }
        public XDGField[] Residuals { get; set; }
        public XDGField[] UBackup { get; set; }
        public XDGField[] AgglomeratedUBackup { get; set; }
        #endregion
        #region LevelSets
        public LevelSet LevelSet { get; set; }
        public LevelSet LevelSetTwo { get; set; }
        public LevelSet LsTBO { get; set; } // This stands for LevelSet To Be Optimized and is a reference to either LevelSetOne or LevelSetTwo
        public LevelSet levelSetBackup { get; set; }
        public IOptiLevelSet LevelSetOpti { get; set; }
        public IOptiLevelSet LevelSetOptiBackup { get; set; }
        public IGrid LevelSetOptiGrid { get; set; }
        public LevelSet LevelSetStep { get; private set; }
        public LevelSetTracker LsTrkStep { get; private set; }
        public ContinuityProjection ContinuityProjection { get; set; }
        private SinglePhaseField ShockLevelSetField { get; set; }
        private SpeciesId[] _speciesToEvaluate_Ids { get; set; }
        public SpeciesId[] SpeciesToEvaluate_Ids {
            get {
                if(_speciesToEvaluate_Ids == null) {
                    if(this.Control.SpeciesToEvaluate != null) {
                        _speciesToEvaluate_Ids = new SpeciesId[Control.SpeciesToEvaluate.Length];
                        for(int i = 0; i < _speciesToEvaluate_Ids.Length; i++) {
                            _speciesToEvaluate_Ids[i] = LsTrk.GetSpeciesId(Control.SpeciesToEvaluate[i]);
                        }
                    } else {
                        throw new NotSupportedException("You did not specify any species to evaluate.");
                    }
                }
                return _speciesToEvaluate_Ids;
            }
        }
        public string[] SpeciesToEvaluate {
            get { return this.Control.SpeciesToEvaluate; }
        }
        #endregion
        #region Operator stuff
        public XDifferentialOperatorMk2 XSpatialOperator { get; set; }
        public XDifferentialOperatorMk2 Op_obj { get; set; }
        public IEvaluatorNonLin Eval_r { get; set; }
        public IOptProb Oproblem {get;set;}
        public IEvaluatorNonLin Eval_R { get; set; }
        public MultigridOperator.ChangeOfBasisConfig[][] MultiGridOperatorConfig { get; set; }
        public MultiphaseCellAgglomerator MultiphaseAgglomerator { get; set; }
        public IList<TimeStepConstraint> TimeStepConstraints { get; set; }
        public MultigridOperator multOp { get; set; }
        /// <summary>
        /// The Jacobi operators belonging to r and R, used when Linearization = JacobiOperator
        /// </summary>
        public IDifferentialOperator r_JacobiOperator { get; set; }
        public IDifferentialOperator R_JacobiOperator { get; set; }
        public PlotDriver plotDriver { get; set; }
        #endregion
        #region Coordinate Vectors and Mappings
        public CoordinateVector m_UnknownsVector { get; set; }
        public CoordinateVector UnknownsVector {
            get {
                if(m_UnknownsVector == null)
                    m_UnknownsVector = new CoordinateVector(ConservativeFields);
                return m_UnknownsVector;
            }
        }
        public CoordinateMapping UnknownsMap {
            get {
                return UnknownsVector.Mapping;
            }
        }
        public CoordinateVector ResidualVector { get; set; }
        public CoordinateMapping ResidualMap {
            get {
                return ResidualVector.Mapping;
            }
        }
        public CoordinateVector obj_f_vec { get; set; }
        public CoordinateMapping obj_f_map {
            get {
                return obj_f_vec.Mapping;
            }
        }
        #endregion
        #region Lists Tracking the History of OPT Parameters, functional values, Steps
        public List<double> Gammas { get; set; }
        public List<double> Kappas { get; set; }
        public List<double> Alphas { get; set; }
        public List<double[]> LevelSetOptiParams { get; set; }
        public List<double> Mus { get; set; }
        public List<double> ResNorms { get; set; }
        public List<double> obj_f_vals { get; set; }
        public List<double> StepCount { get; set; }
        public List<double[]> Steps { get; set; }
        #endregion
        #region Matrices
        public MsrMatrix LHS {get;set;}
        public MsrMatrix Jr_U { get; set; }
        public MsrMatrix Jr { get; set; }
        public MsrMatrix Jr_phi { get; set; }
        public MsrMatrix Jobj_U { get; set; }
        public MsrMatrix Jobj { get; set; }
        public MsrMatrix Jobj_phi { get; set; }
        BlockMsrMatrix LeftMul { get; set; }
        #endregion
        #region Arrays
        public double[] gradf_U { get; set; }
        public double[] gradf { get; set; }
        public double[] JacR0TimesStep { get; set; }
        public double[] JacR1TimesStep { get; set; }
        public double[] JacR0TimesR0 { get; set; }
        public double[] RHS { get; set; }
        public double[] stepIN { get; set; }
        public double[] stepIN_agg { get; set; }
        public double[] stepDL { get; set; }
        public double[] lambda { get; set; }
        public double[] merit_del { get; private set; }
        public double[] stepUphi { get; set; }

        #endregion
        #region Optimization Variables
        public double gamma { get; set; }
        public double mu { get; set; }
        public double kappa { get; set; }
        public double next_mu { get; set; }
        public double m_alpha { get; set; }
        double res_l1 { get; set; } = 0;
        /// <summary>
        /// computes the predicted value of the merit function for a step s (used for Globalization)
        /// $$ f_m(s) \approx f_m(0) + s^T \beta f_m'(0) $$
        /// 
        /// we have that f_m'(0)= 0.5* (J_R * R + mu* J_r *r)
        /// /// </summary>
        /// <param name="step"></param>
        /// <param name="beta"></param>
        /// <returns>the predicted value of the function</returns>
        //public double predictedMerit(double[] step,double alpha, double beta) {
        Func<double, double[], double, double> predictedMerit { get; set; }
        Func<double[], double, double> actualMerit { get; set; }
        public int CurMinIter { get; private set; }
        bool IncreaseDegreeNextTS { get; set; } = false;
        double NormalizingConstant { get; set; }
        /// <summary>
        /// Just the Current step number
        /// </summary>
        public int CurrentStepNo { get; set; }
        /// <summary>
        /// current l2 norm of constraint
        /// </summary>
        public double res_l2 { get; set; }
        /// <summary>
        /// current value of the objective function
        /// </summary>
        public double obj_f { get; set; }
        /// <summary>
        /// current value of the objective function
        /// </summary>
        public double f_phi { get; set; }
        /// <summary>
        /// Current Agglomeration Threshold
        /// </summary>
        public double CurrentAgglo { get; set; }
        /// <summary>
        /// Current L2 Norm of constraint (L2 means that Residual is regarded as XDG field)
        /// </summary>
        public double res_L2 { get; set; }
        //public double en_res_L2 { get; set; }
        /// <summary>
        /// Current Delta (for Dogleg)
        /// </summary>
        public double TrustRegionDelta { get; set; }
        /// <summary>
        /// initial Residual Norm
        /// </summary>
        public double InitResNorm { get; set; }
        /// <summary>
        /// Initial objective value
        /// </summary>
        public double Init_obj_f { get; set; }
        /// <summary>
        /// Maximum Number of SQP Iterations where Reinitialization() will be performed
        /// </summary>
        public int ReiniTMaxIter { get; private set; }

        #endregion
        #region Member Methods
        /// <summary>
        /// Override this to include Optimization variables
        /// </summary>
        /// <param name="timestepno"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        protected override TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t) {
            if(this.LevelSetOpti is SplineOptiLevelSet spliny) {
                return new IDTTimeStepInfo(t, this.CurrentSessionInfo, timestepno, this.IOFields,
                        this.gamma, this.kappa, this.m_alpha, this.mu, this.res_l2, this.obj_f, this.Alphas,
                        this.Gammas,this.Kappas, this.Mus, this.ResNorms, this.obj_f_vals, this.StepCount,
                        this.LevelSetOpti.GetParamsAsArray(),spliny.y );
            } else {
                return new IDTTimeStepInfo(t, this.CurrentSessionInfo, timestepno, this.IOFields,
                        this.gamma, this.kappa, this.m_alpha, this.mu, this.res_l2, this.obj_f, this.Alphas,
                        this.Gammas, this.Kappas, this.Mus, this.ResNorms, this.obj_f_vals, this.StepCount,
                        this.LevelSetOpti.GetParamsAsArray(), null);
            }
                   
        }
        /// <summary>
        /// Creates or Loads the grid
        /// </summary>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        protected override IGrid CreateOrLoadGrid() {
            IGrid grid = null;
            switch(Control.getGridFrom) {
                case GetGridFrom.DB:
                    if (Control.MeshPath != null)
                    {
                        grid = GridImporter.Import(Control.MeshPath);
                        var tmp = new byte[] { 1, 2, 3, 4, 5 };
                        for (int i = 0; i < Control.EdgTagNames.Length; i++)
                        {
                            grid.EdgeTagNames.Add(tmp[i], Control.EdgTagNames[i]);
                        }
                        grid.DefineEdgeTags(Control.EdgTagFunc);
                        return grid;
                    }
                    else
                    {
                        grid = base.CreateOrLoadGrid();
                        CompressibleEnvironment.Initialize(grid.SpatialDimension);
                        return grid;
                    }
                
                case GetGridFrom.GridFunc:
                default:
                // Optional: load grid from CNS calculation from which the shock level set has been reconstructed
                if(Control.ShockLevelSet_UseShockGrid) {
                    IDatabaseInfo dbi = DatabaseInfo.Open(Control.ShockLevelSet_Db);
                    ISessionInfo si = dbi.Controller.GetSessionInfo(Control.ShockLevelSet_Info.Item1);
                    ITimestepInfo tsi;

                    if(Control.ShockLevelSet_Info.Item2.MajorNumber < 0.0) {
                        tsi = si.Timesteps.Last();
                    } else {
                        throw new NotSupportedException("Currently, only the last time step can be used.");
                    }

                    grid = dbi.Controller.DBDriver.LoadGrid(tsi.GridID, dbi);

                    Console.WriteLine("Using grid from LEVEL SET RECONSTRUCTION");
                } else {
                    grid = base.CreateOrLoadGrid();
                }

                CompressibleEnvironment.Initialize(grid.SpatialDimension);

                return grid;
            }
        }
        /// <summary>
        /// This is called in CreateFields. Here the user must initialize the ConservativeFields and if 
        /// wanted also the DerivedFields. The DerivedFields also need to be registered.
        /// </summary>
        /// <returns></returns>
        public abstract void CreateConservativeFields(LevelSetTracker LsTrk,int DgDegree);

        /// <summary>
        /// This method creates the Fields used in this solver, also various lists are being created 
        /// as well as the data structure needed in the program (e.g. LHS,RHS, ...)
        /// It is the implementation of the Application member <see cref="Application{T}.CreateFields"/>
        /// </summary>
        /// <returns></returns>
        protected override void CreateFields() {

            CheckControl();

            InitializeLevelSets();

            //Initialize mandatory conservative fields and Residuals by abstract Function
            switch(Control.solRunType) {
                case SolverRunType.Standard:
                    CreateConservativeFields(this.LsTrk, Control.SolDegree);
                break;
                case SolverRunType.PContinuation:
                    CreateConservativeFields(this.LsTrk, Control.DgDegree_Start); //here we start with the lowest degree
                break;
                default:
                throw new NotImplementedException("this type of Solver run is undefined");
            }

            //Initialize all dependent Fields
            InitDependentFields();

            //Initialize the MultigridOperatorConfig
            InitializeMultiGridOpConfig();
            //Create history Lists to track the parameters
            Gammas = new List<double>();
            Gammas.Add(Control.Gamma_Start);

            Kappas = new List<double>();
            Alphas = new List<double>();
            Alphas.Add(Control.Alpha_Start);
            Mus = new List<double>();
            Mus.Add(Control.Mu_Start);
            ResNorms = new List<double>();
            obj_f_vals = new List<double>();
            StepCount = new List<double>();
            StepCount.Add(0);
            LevelSetOptiParams = new List<double[]>();
            LevelSetOptiParams.Add(LevelSetOpti.GetParamsAsArray());
            Steps = new List<double[]>();
            Steps.Add(stepIN);
            gamma = Control.Gamma_Start;
            CurMinIter = Control.minimalSQPIterations[0];
            ReiniTMaxIter = Control.ReiniTMaxIters[0];
            CurrentAgglo = Control.AgglomerationThreshold;
            //chose the Merit Function used in Globalization
            ChooseMeritFunction();
            ChooseFphiType();
            ChooseOptProblem();
            InitObjFields();
            //Register them (also registers other fields)
            RegisterFields();
        }
        /// <summary>
        /// Chooses the Optimization Problem that will be solved
        /// </summary>
        /// <exception cref="Exception"></exception>
        public void ChooseOptProblem() {
            switch(Control.optProblemType) {
                case OptProblemType.FullEnRes:
                    Oproblem = new SFFullEnRes(XSpatialOperator, pDiff: 1);
                break;
                case OptProblemType.EnResOnlyNearBand:
                    Oproblem = new SFNearBandEnRes(XSpatialOperator, pDiff: 1);
                break;
                case OptProblemType.EnResOnlyCutCells:
                    Oproblem = new SFCutCellEnRes(XSpatialOperator, pDiff: 1);
                break;
                case OptProblemType.RankineHugoniotFull:
                case OptProblemType.RankineHugoniotOnlyInterface:
                    Oproblem = new SFRankineHugoniotBase(XSpatialOperator,Op_obj, pDiff: 0);
                break;
                default: 
                throw new Exception("This OptProblem is not implemented yet");
            }

        }

        /// <summary>
        /// does some sanity checks on the control
        /// </summary>
        /// <exception cref="NotImplementedException"></exception>
        private void CheckControl() {
            var maxDeg = Control.SolDegree;
            if(Control.Alpha_Start <= 0 || Control.Alpha_Start > 1) {
                throw new NotSupportedException($"Control.Alpha_Start length must be between 0 and 1 but is {Control.Alpha_Start}");
            }
            switch(Control.solRunType) {
                case SolverRunType.PContinuation:
                    if(Control.ReiniTMaxIters.Length < maxDeg+1) {
                        throw new NotSupportedException($"Control.ReiniTMaxIters length mus bet at least {maxDeg + 1}");
                    }
                    if(Control.TerminationMinNs.Length < maxDeg + 1) {
                        throw new NotSupportedException($"Control.TerminationMinNs length must be at least {maxDeg+1}");
                    }
                    if(Control.tALNRs.Length < maxDeg + 1) {
                        throw new NotSupportedException($"Control.tALNRs length must be at least {maxDeg+1}");
                    }
                    if(Control.minimalSQPIterations.Length < maxDeg + 1) {
                        throw new NotSupportedException($"Control.MinPIter length must be at least {maxDeg+1}");
                    }

                break;
            }
            Console.WriteLine("*****************************************************************************************");
            Console.WriteLine("*  X   X   DDD   GGG   -   III  SSS  TTTTT  *");
            Console.WriteLine("*   X X    D  D G      -    I  S      T     *");
            Console.WriteLine("*    X     D  D G  GG  -    I   SSS   T     *      XDG Implicit Shock Tracking Solver    ");
            Console.WriteLine("*   X X    D  D G   G  -    I      S  T     *");
            Console.WriteLine("*  X   X   DDD   GGG   -   III  SSS   T     *");
            Console.WriteLine("*****************************************************************************************");
            Console.WriteLine($"*  ");
            Console.WriteLine($"*  Solver:                             {this.GetType().ToString()}");
            Console.WriteLine($"*  optimization problem:               {Control.optProblemType}");
            Console.WriteLine($"*  merit function:                     {Control.MeritFunctionType}");
            Console.WriteLine($"*  level set type:                     {Control.OptiLevelSetType}");
            Console.WriteLine($"*  solver type:                        {Control.solRunType}");
            Console.WriteLine($"*  f_phi type:                         {Control.fphiType}");
            Console.WriteLine($"*  Globalization Strategy:             {Control.GlobalizationStrategy}");
            Console.WriteLine($"*  ");
            Console.WriteLine("******************************* Params ************************************************");
            Console.WriteLine($"*  Solution Degree:                    {Control.SolDegree}");
            Console.WriteLine($"*  OptiLevelSet Degree:                {Control.OptiLevelSetDegree}");
            Console.WriteLine($"*  Level Set Degree:                   {Control.LevelSetDegree}");
            Console.WriteLine($"*  Level Set Two Degree:               {Control.LevelSetTwoDegree}");
            Console.WriteLine($"*  Max Iterations:                     {Control.MaxIterations}");
            //Console.WriteLine($"*  Linear Solver:                      {Control.LinearSolver?.GetType().ToString() ?? "None"}");
            Console.WriteLine($"*  Agglomeration Threshold:            {Control.AgglomerationThreshold}");
            //Console.WriteLine($"*  Is Far Config:                      {Control.isFarConfig}");
            Console.WriteLine($"*  Regularization Params:");
            Console.WriteLine($"*    Gamma Max:                        {Control.Gamma_Max}");
            Console.WriteLine($"*    Gamma Min:                        {Control.Gamma_Min}");
            Console.WriteLine($"*    Gamma Start:                      {Control.Gamma_Start}");
            Console.WriteLine($"*  Globalization Params:");
            Console.WriteLine($"*    Alpha Min:                        {Control.Alpha_Min}");
            Console.WriteLine($"*    Alpha Start:                      {Control.Alpha_Start}");
            Console.WriteLine($"*    Mu Start:                         {Control.Mu_Start}");
            Console.WriteLine($"*    Mu Max:                           {Control.Mu_Max}");
            Console.WriteLine($"*    Mu Omega:                         {Control.Mu_Omega}");
            Console.WriteLine($"*    Mu Rho:                           {Control.Mu_Rho}");
            Console.WriteLine($"*  Level Set penalty Params:");
            Console.WriteLine($"*    Kappa Xi:                         {Control.Kappa_Xi}");
            Console.WriteLine($"*    Kappa M:                          {Control.Kappa_M}");
            Console.WriteLine($"*    Kappa Min:                        {Control.Kappa_Min}");
            Console.WriteLine($"*    Kappa v:                          {Control.Kappa_v}");
            Console.WriteLine($"*  Reinitialization Params:");
            Console.WriteLine($"*    Apply ReiInit:                    {Control.ApplyReiInit}");
            Console.WriteLine($"*    c1:                               {Control.reInit_c1}");
            Console.WriteLine($"*    c2:                               {Control.reInit_c2}");
            Console.WriteLine($"*    c3:                               {Control.reInit_c3}");
            Console.WriteLine($"*  Adaptive Regularization Params:");
            Console.WriteLine($"*    L:                                {Control.L}");
            Console.WriteLine($"*    sigma_1:                          {Control.sigma_1}");
            Console.WriteLine($"*    sigma_2:                          {Control.sigma_2}");
            Console.WriteLine($"*    tauGamma:                         {Control.tauGamma}");
            //Console.WriteLine("****************************** Initialization & Level Set Config *************************");
            //Console.WriteLine($"*  Level Set One Initial Value:        {(Control.LevelSetOneInitialValue != null ? "Defined" : "None")}");
            //Console.WriteLine($"*  LsOne NegSpecies:                   {Control.LsOne_NegSpecies}");
            //Console.WriteLine($"*  LsOne PosSpecies:                   {Control.LsOne_PosSpecies}");
            //Console.WriteLine($"*  LsTwo NegSpecies:                   {Control.LsTwo_NegSpecies}");
            //Console.WriteLine($"*  LsTwo PosSpecies:                   {Control.LsTwo_PosSpecies}");
            //Console.WriteLine($"*  LsOne Species Pairs:                {Control.LsOne_SpeciesPairs.GetLength(0)} Pairs Defined");
            //Console.WriteLine($"*  LsTwo Species Pairs:                {Control.LsTwo_SpeciesPairs.GetLength(0)} Pairs Defined");
            //Console.WriteLine($"*  Write LS Coordinates:               {Control.WriteLSCoordinates}");
            //Console.WriteLine($"*  OptiLevelSet Parameters Defined:    {Control.OptiLevelSet_ParamNames?.Count ?? 0}");
            //Console.WriteLine($"*  Near Region Width:                  {Control.NearRegionWidth}");
            //Console.WriteLine($"*  Flux s_alpha:                       {Control.flux_s_alpha}");
            //Console.WriteLine("****************************** Data Bases ****************************************");

            //Console.WriteLine($"*  Mesh Path:                          {Control.MeshPath ?? "None"}");
            //Console.WriteLine($"*  Shock Level Set Db:                 {Control.ShockLevelSet_Db ?? "None"}");
            //Console.WriteLine($"*  Shock Level Set Use Shock Grid:     {Control.ShockLevelSet_UseShockGrid}");
            //Console.WriteLine($"*  Shock Level Set Seed From Db:       {Control.ShockLevelSet_SeedFromDb}");
            //Console.WriteLine($"*  Partially Fix Level Set Space Time: {Control.PartiallyFixLevelSetForSpaceTime}");
            //Console.WriteLine($"*  Save Matrices:                      {Control.SaveMatrices}");
            //Console.WriteLine("*****************************************************************************************");


        }


        /// <summary>
        /// 
        /// </summary>
        public void RegisterFields() {

            // Register fields for plotting
            this.m_IOFields.Add(LevelSet);
            if(Control.IsTwoLevelSetRun) {
                this.m_IOFields.Add(LevelSetTwo);
            }
            
            this.m_IOFields.AddRange(ConservativeFields);
            this.m_IOFields.AddRange(Residuals);
            this.m_IOFields.AddRange(obj_f_fields);
            this.m_IOFields.AddRange(LagMult);
            this.m_IOFields.Add(personField);

            // Register fields, e.g., for applying the initial conditions
            this.m_RegisteredFields.Add(LevelSet);
            if(Control.IsTwoLevelSetRun) {
                this.m_RegisteredFields.Add(LevelSetTwo);
            }
            this.m_RegisteredFields.AddRange(ConservativeFields);
            this.m_RegisteredFields.AddRange(Residuals);
            this.m_RegisteredFields.AddRange(obj_f_fields);
            this.m_RegisteredFields.AddRange(LagMult);
            this.m_RegisteredFields.Add(personField);


#if DEBUG
            //this.m_IOFields.AddRange(ConservativeFieldsReinitialized);
            //this.m_RegisteredFields.AddRange(ConservativeFieldsReinitialized);
#endif
        }

        /// <summary>
        /// Here the magic happens Divided in sub routines
        /// </summary>
        /// <param name="TimestepNo"></param>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //if(CurrentAgglo > 0 && TimestepNo % Control.ImmediatePlotPeriod ==0 && Control.ImmediatePlotPeriod!= -1) {
            //    MultiphaseAgglomerator.PlotAgglomerationPairs(Control.ProjectName + "_aggloMap" + "_" + CurrentStepNo);
            //}
            using (new FuncTrace())
            {
                if (TimestepNo == 1)
                {
                    //set initial kappa
                    //compute using kappa=1
                    kappa = 1;
                    var f_phi_vec = ComputeFphi();
                    f_phi = f_phi_vec.InnerProd(f_phi_vec);
                    if (f_phi > 1e-10)
                    {
                        kappa = Math.Sqrt(obj_f / f_phi_vec.InnerProd(f_phi_vec));
                    }
                    else
                    {
                        kappa = 1e3;
                    }
                    Kappas.Add(kappa);

                    Console.WriteLine($"Initial Value : l2: |r|= {string.Format("{0:#.#######E+00}", res_l2)}, |f_err|={string.Format("{0:#.#######E+00}", obj_f)}, | f_phi |={string.Format("{0:#.#######E+00}", f_phi)},");
                    UpdateAgglomerator();
                }

                // potentially change basis
                switch (Control.solRunType)
                {
                    case SolverRunType.PContinuation:
                        ChangeOfBasis();
                        break;
                }

                //reinitialization if p>0
                Reinitialization();

                Console.WriteLine("");
                Console.WriteLine($"Starting Opt-Iter No.{TimestepNo}: ");

                //backup of non-agglomerated state and Level Set
                createBackUp();

                //just to be safe, project spline Level Set on LsTBO
                if (LevelSetOpti is SplineOptiLevelSet spliny)
                {
                    LsTBO.Clear();
                    LsTBO.ProjectFromForeignGrid(1.0, spliny);
                }

                //Compute/Assemble linear System LHS*(stepIN;lambda)=RHS
                ComputeSystem();
                //solve it for step: (stepIN;lambda)
                SolveSystem();

                //save step to XDGField
                var StepOptiLevelSet = GetStepOptiLevelSet(stepIN);
                SaveStepToField(stepIN, StepOptiLevelSet);

                // makes the step continues if SinglePhaseField as OptiLevelSet is used
                //if(Control.OptiLevelSetType == OptiLevelSetType.SinglePhaseField && LevelSetOpti.GetGrid().iGeomCells.Count > 1) {
                //    // Does a continuity projection (hopefully) of the step
                //    MakeLevelSetContinous(StepOptiLevelSet);
                //    SaveStepFieldToStep(stepIN, StepOptiLevelSet);
                //    var StepOptiLevelSet_afterprojection = GetStepOptiLevelSet(stepIN);
                //    SaveStepToField(stepIN, StepOptiLevelSet_afterprojection);
                //}

                var orgStep = stepIN.CloneAs();
                // 
                try
                {
                    SQPStep();
                }
                catch
                {
                    // if something goes wrong here plot the failed state
                    try
                    {
                        UpdateDerivedVariables();
                        this.PlotCurrentState(0.0, 666, this.Control.SuperSampling);
                    }
                    catch
                    {
                        this.PlotCurrentState(0.0, 666, this.Control.SuperSampling);
                    }
                    throw;
                }

                ComputeAndWriteResiduals();
                WriteStuffAndAddToLists(TimestepNo);
                CmpRegParam();
                CurrentStepNo++;
                LsTrk.PushStacks();
                //this is called only due to plotting reasons
                if (ConservativeFields[0].Basis.Degree > 0)
                {
                    GetPerssonSensor(true);
                }
                //Check for termination or Increase of P degree
                CheckTermination();
                UpdateAgglomerator();
                // must return something (normally the time step)
                return 0;
            }
        }

        //Changes the Basis of the current solution, making it one DG higher 
        public void ChangeOfBasis() {
            using (new FuncTrace())
            {
                int CurrentDeg = ConservativeFields[0].Basis.Degree;
                if (IncreaseDegreeNextTS && CurrentDeg < Control.SolDegree)
                {
                    {
                        Console.WriteLine("################################################################################################");
                        Console.WriteLine("################################################################################################");
                        Console.WriteLine("###########################                                      #####################################################################");
                        Console.WriteLine("###########################      Changing polynomial degree      #########################################");
                        Console.WriteLine("###########################                                      #####################################################################");
                        Console.WriteLine($"###########################             From {CurrentDeg} to {CurrentDeg + 1}              #########################################");
                        Console.WriteLine("###########################                                      #####################################################################");
                        Console.WriteLine("################################################################################################");
                        Console.WriteLine("################################################################################################");

                        //change to Degree +1
                        CurrentDeg = CurrentDeg + 1;
                        ChangeTOBasisofDegree(CurrentDeg);

                        //increase Current Minimal iTerations by the respective Minimal iterations wanted for the polynomial degree

                        CurMinIter = CurrentStepNo + Control.minimalSQPIterations[CurrentDeg];
                        ReiniTMaxIter = CurrentStepNo + Control.ReiniTMaxIters[CurrentDeg];

                        IncreaseDegreeNextTS = false;
                        Console.WriteLine("");
                        Console.Write("OpEval after change:");
                        ComputeAndWriteResiduals();
                        Console.WriteLine("");
                        Console.WriteLine($" current DGdegree={CurrentDeg}, min It.:{Control.minimalSQPIterations[CurrentDeg]}, MaxReInit:{Control.ReiniTMaxIters[CurrentDeg]}, ReInitTol:{Control.reInitTols[CurrentDeg]}," +
                            $" TermN:{GetCurrentTermN()}, tALNR:{Control.tALNRs[CurrentDeg]} ");
                        Console.WriteLine("");
                        Console.WriteLine("################################################################################################");
                    }
                }
            }
        }

        public void Reinitialization() {
            int cDeg= ConservativeFields[0].Basis.Degree;
            switch(Control.reInitMode) {
                case ReInitMode.OneTolForEachP:
                    Control.reInit_c1 = Control.reInitTols[cDeg];
                break; 
            }
            if(ConservativeFields[0].Basis.Degree > 0 && Control.ApplyReiInit && CurrentStepNo < ReiniTMaxIter) {
                Reinitialization(false);
            }
        }

        /// <summary>
        /// simply updates the agglomerator
        /// </summary>
        public void UpdateAgglomerator() {
            if(CurrentAgglo > 0) {
                try {
                    var testAgg = LsTrk.GetAgglomerator(SpeciesToEvaluate_Ids, GetGlobalQuadOrder(), CurrentAgglo, ExceptionOnFailedAgglomeration: false);
                    MultiphaseAgglomerator = LsTrk.GetAgglomerator(SpeciesToEvaluate_Ids, GetGlobalQuadOrder(), CurrentAgglo, ExceptionOnFailedAgglomeration: false); 
                } catch {
                    Console.WriteLine("Failed to update Agglomerator");
                }
            }
        }
        #endregion
        #region Utility Methods
        /// <summary>
        /// A helper routine one can use in jupyter notebooks to obtain a Main object 
        /// where everything (e.g. SpatialOperator, Fields, etc.) is set up. 
        /// It is mainly used for testing.
        /// </summary>
        public void InitializeEverything() {
            this.SetUpEnvironment();
            this.CreateFields();
            this.SetInitial(0);
            this.CreateEquationsAndSolvers(null);
            PlotCurrentState(0.0, 0, Control.SuperSampling);
        }
        /// <summary>
        /// Function that Initializes LevelSets and Tracker
        /// </summary>
        /// <exception cref="NotImplementedException"></exception>
        public void InitializeLevelSets() {

            #region Create level set objects, and the level set tracker
            // Level set one (usually the geometry)
            this.LevelSet = new LevelSet(new Basis(this.GridData, this.Control.LevelSetDegree), "levelSet");
            int lsNumber;
            if(Control.IsTwoLevelSetRun) {
                this.LevelSetTwo = new LevelSet(new Basis(this.GridData, this.Control.LevelSetTwoDegree), "levelSetTwo");
                this.LevelSetStep = new LevelSet(new Basis(this.GridData, this.Control.LevelSetTwoDegree), "levelSetTwo_Step");
                base.LsTrk = new LevelSetTracker((GridData)this.GridData, Control.CutCellQuadratureType, 1, Control.SpeciesTable, this.LevelSet, this.LevelSetTwo);
                this.LsTrkStep = new LevelSetTracker((GridData)this.GridData, Control.CutCellQuadratureType, 6, Control.SpeciesTable, this.LevelSet, this.LevelSetStep);
                LsTBO = LevelSetTwo;
                lsNumber = 1;
            } else {
                this.LevelSetStep = new LevelSet(new Basis(this.GridData, this.Control.LevelSetDegree), "levelSet_Step");
                base.LsTrk = new LevelSetTracker((GridData)this.GridData, Control.CutCellQuadratureType, 1, new string[] { Control.LsOne_NegSpecies, Control.LsOne_PosSpecies }, this.LevelSet);
                this.LsTrkStep = new LevelSetTracker((GridData)this.GridData, Control.CutCellQuadratureType, 6, new string[] { Control.LsOne_NegSpecies, Control.LsOne_PosSpecies }, this.LevelSetStep);
                LsTBO = LevelSet;
                lsNumber = 0;
            }
            if(Control.LevelSetGridFunc == null) {
                Control.LevelSetGridFunc = Control.GridFunc;
                LevelSetOptiGrid = Grid;
            } else {
                LevelSetOptiGrid = Control.LevelSetGridFunc();
            }

            Basis LevelSetOptiBasis = new Basis(LevelSetOptiGrid, Control.OptiLevelSetDegree);
            //Debugger.Launch();
            switch(Control.OptiLevelSetType) {
                case OptiLevelSetType.GlobalLevelSet:
                this.LevelSetOpti = new GlobalOptiLevelSet(this.Grid, Control.OptiLevelSet_ParamNames, Control.OptiLevelSet_ParamValues, Control.OptiLevelSet_Param_Functions, Control.OptiLevelSetDegree, isOrthormal: Control.OptiLSIsOrthonormal);
                break;
                case OptiLevelSetType.SinglePhaseField:
                this.LevelSetOpti = new SinglePhaseFieldOptiLevelSet(LevelSetOptiBasis, "LevelSetOpti", LsTBO, LsTrk, lsNumber);
                LevelSetOpti.AssembleTracker();
                break;

                case OptiLevelSetType.SpecFemField:
                SpecFemBasis LevelSetOptiSpecBasis;
                switch(Control.getGridFrom) {
                    case GetGridFrom.GridFunc:
                    LevelSetOptiSpecBasis = new SpecFemBasis((GridData)LevelSetOptiGrid.iGridData, Control.OptiLevelSetDegree);
                    break;
                    case GetGridFrom.DB:
                    LevelSetOptiSpecBasis = new SpecFemBasis((GridData)this.Grid.iGridData, Control.OptiLevelSetDegree);
                    break;
                    default:
                    throw new NotImplementedException("not supported");
                }
                if(LsTBO.Basis.Degree < LevelSetOptiSpecBasis.ContainingDGBasis.Degree) {
                    Console.WriteLine("WARNING: The Basis of the DG LevelSet is " + LsTBO.Basis.Degree + " and therefore is lower than the SpecFem OptiLevelSet Degree " + LevelSetOptiSpecBasis.ContainingDGBasis.Degree + " --> DGLevelSet will eventually become discontinuous");
                }
                this.LevelSetOpti = new SpecFemOptiLevelSet(LevelSetOptiSpecBasis, LsTBO, true, LsTrk, lsNumber);
                break;
                case OptiLevelSetType.SplineLevelSet:
                this.LevelSetOpti = new SplineOptiLevelSet(LevelSetOptiBasis, LsTBO, LsTrk, lsNumber);
                break;
            }
            NormalizingConstant = LevelSetOpti.GetParam(0);
            #endregion
        }
        /// <summary>
        /// This method needs to be implemented by deriving classes and initialize the multi grid-operator config
        /// </summary>
        public abstract void InitializeMultiGridOpConfig();

        /// <summary>
        /// initialize fields dependent on ConVars (e.g. stepField, PersonFIeld, EnrichedField,Residuals)
        /// </summary>
        public void InitDependentFields() {
            m_UnknownsVector = new CoordinateVector(ConservativeFields);
            ConservativeFieldsReinitialized = new XDGField[ConservativeFields.Length];
            //Initialize Fields for the Step and Register
            originalStep = new XDGField[ConservativeFields.Length];
            LagMult = new XDGField[ConservativeFields.Length];
            for(int i = 0; i < ConservativeFields.Length; i++) {
#if DEBUG
                ConservativeFieldsReinitialized[i] = new XDGField(ConservativeFields[i].Basis, ConservativeFields[i].Identification + "_reInit");
#endif
                var stepBasis = new XDGBasis(LsTrkStep, ConservativeFields[i].Basis.Degree);
                originalStep[i] = new XDGField(stepBasis, ConservativeFields[i].Identification + "_step");
                LagMult[i] = new XDGField(ConservativeFields[i].Basis, ConservativeFields[i].Identification + "_lam");
            }
            //initialize Person-field, storing the sensor-values
            personField = new XDGField(new XDGBasis(LsTrk, 0), "personField");

            #region Create fields for residuals and the Enriched XDGFields
            Residuals = new XDGField[ConservativeFields.Length];
            for(int i = 0; i < ConservativeFields.Length; i++) {
                this.Residuals[i] = new XDGField(new XDGBasis(base.LsTrk, ConservativeFields[i].Basis.Degree), "res_" + ConservativeFields[i].Identification);
            }
            ResidualVector = new CoordinateVector(Residuals);
            #endregion
        }
        /// <summary>
        /// Made an extra method out of this, because it needs to be executed after Operator Initialization
        /// </summary>
        public void InitObjFields() {
            obj_f_fields = Oproblem.CreateObjField(Residuals);
            obj_f_vec = new CoordinateVector(obj_f_fields);
        }

        /// <summary>
        /// Initialized all MAtrices and Vectors used with the correct sizes (e.g dep. on polynomial degree)
        /// we want to call it at the beginning of a SolverRun an if the Polynomial Degree/LevelSet or Mesh is changed
        /// </summary>
        public void InitializeMatricesAndVectors() {

            int length_r = (int) new CoordinateMapping(ConservativeFields).TotalLength;
            int length_obj = Oproblem.GetObjLength(ConservativeFields);
            int length_phi = LevelSetOpti.GetLength();

            Jr = new MsrMatrix(length_r, length_r + length_phi, 1, 1);
            Jr_U = new MsrMatrix(length_r, length_r, 1, 1);
            Jr_phi = new MsrMatrix(length_r, length_phi, 1, 1);

            Jobj = new MsrMatrix(length_obj, length_r + length_phi, 1, 1);
            Jobj_U = new MsrMatrix(length_obj, length_r, 1, 1);
            Jobj_phi = new MsrMatrix(length_obj, length_phi, 1, 1);

            gradf = new double[length_r + length_phi];
            gradf_U = new double[length_r];
            stepUphi = new double[length_r + length_phi];
            merit_del = new double[length_r + length_phi];


            LHS = new MsrMatrix(2 * length_r + length_phi);
            LHS.AssumeSymmetric = true;
            RHS = new double[2 * length_r + length_phi];

            stepIN = new double[2 * length_r + length_phi];
            stepIN_agg = new double[2 * length_r + length_phi];
            stepDL = new double[2 * length_r + length_phi];
            lambda = new double[length_r];
        }

        /// <summary>
        /// This method contains a switch choosing the wished Merit Function
        /// </summary>
        /// <exception cref="ArgumentException"></exception>
        private void ChooseMeritFunction() {
            switch(Control.MeritFunctionType) {
                case MeritFunctionType.L1Merit:
                ///// <summary>
                ///// computes the merit function for a step s (used for Globalization)
                ///// $$ f_m(s) = 0.5*\Vert R(z_k + s) \Vert^2 + \mu \Vert r(z_k + s) \Vert_1$$
                ///// 
                ///// cf. (18.27) Numerical Optimization - Nocedal and Wright 2006
                ///// </summary>
                actualMerit = delegate (double[] step, double m_alpha) {
                    double[] objStep = new double[(int)obj_f_map.TotalLength];
                    double[] resStep = new double[(int)ResidualMap.TotalLength];
                    bool succes = AccumulateStep(step, m_alpha);

                    //TransformFromAggToSourceSpace();

                    ComputeResiduals(objStep, resStep);
                    double[] f_phi_vec = ComputeFphi();
                    double resStep_l1 = 0;
                    for(int i = 0; i < resStep.Length; i++) {
                        resStep_l1 += resStep[i].Abs();
                    }
                    resetStep(); // reset the step
                    return 0.5 * (objStep.InnerProd(objStep)+ f_phi_vec.InnerProd(f_phi_vec)) + mu * resStep_l1;
                };

                /// <summary>
                /// computes the predicted value of the merit function for a step s_k=\m_alpha * s (used for Globalization)
                /// $$ f_m(s) \approx f_m(0) + s^T \beta f_m'(0) $$
                /// 
                /// we have that 
                /// 
                /// $$ s^T * f_m'(0)= s^T * J_R * R  - mu* \vert r \vert_1  $$
                /// 
                /// so in total this should be 
                /// 
                /// $$ f_m(s)=  0.5*R^T*R + \mu * \vert r \vert_1 +  \beta \alpha_k (J_R * R *s - mu* \vert r \vert_1) $$#
                /// 
                /// which is the same as
                /// 
                /// $$ f_m(s)=  0.5*R^T*R + \mu * \vert r \vert_1 +  \beta \alpha_k (J_R * R *s - mu* \vert r \vert_1) $$
                /// /// </summary>
                predictedMerit = delegate (double m_alpha, double[] step, double beta) {

                    stepUphi.SetSubVector(step, 0, (int)UnknownsMap.TotalLength + LevelSetOpti.GetLength());
                    double gradfTimesStep = gradf.InnerProd(stepUphi);
                    //Jr.Transpose().SpMV(mu * beta, ResidualVector, 0, merit_del);
                    //merit_del.AccV(beta, gradf);
                    double resJacRTimesStep_l1 = 0;
                    res_l1 = 0;
                    for (int i = 0; i < JacR0TimesStep.Length; i++)
                    {
                        res_l1 += ResidualVector[i].Abs();
                        resJacRTimesStep_l1 += JacR0TimesStep[i].Abs();
                    }
                    return 0.5 * obj_f * obj_f + mu * res_l1 + beta * m_alpha * (gradfTimesStep + mu * resJacRTimesStep_l1);

                    //return 0.5 * obj_f * obj_f + mu * res_l1 + beta * m_alpha * (gradfTimesStep - mu * res_l1);
                };

                break;
                case MeritFunctionType.L2Merit:
                predictedMerit = delegate (double alpha, double[] step, double beta) {
                    stepUphi.SetSubVector(step, 0, (int)UnknownsMap.TotalLength + LevelSetOpti.GetLength());
                    double gradfTimesStep = gradf.InnerProd(stepUphi);
                    //Jr.Transpose().SpMV(mu * beta, ResidualVector, 0, merit_del);
                    //merit_del.AccV(beta, gradf);
                    return 0.5 * obj_f * obj_f + mu * res_l2 + beta * m_alpha * (gradfTimesStep + mu * JacR0TimesStep.L2Norm());
                    //return 0.5 * obj_f * obj_f + mu * res_l2 * res_l2 + beta * m_alpha * (gradfTimesStep + mu * JacR0TimesR0.InnerProd(stepUphi));
                };
                /// <summary>
                /// computes the l2 merit function for a step s and an alpha (used for Globalization)
                /// $$ f_m(s) = 0.5*(\mu \Vert r(z_k + s) \Vert^2 +\Vert R(z_k + s) \Vert^2)$$
                /// </summary>
                actualMerit = delegate (double[] step, double alpha) {
                    double[] objStep = new double[(int)obj_f_map.TotalLength];
                    double[] resStep = new double[(int)ResidualMap.TotalLength];
                    bool succes = AccumulateStep(step, alpha);
                    double[] f_phi_vec = ComputeFphi();
                    //TransformFromAggToSourceSpace();

                    ComputeResiduals(objStep, resStep);

                    resetStep(); // reset the step
                    return 0.5 * ( objStep.InnerProd(objStep) + f_phi_vec.InnerProd(f_phi_vec))+ mu * resStep.L2Norm();
                    //return (0.5 * (mu * resStep.InnerProd(resStep) + objStep.InnerProd(objStep) + f_phi_vec.InnerProd(f_phi_vec)));
                };
                break;
                default:
                throw new ArgumentException("You need to choose a merit type");
            }
        }
        /// <summary>
        /// deletes old text-files
        /// </summary>
        public static void DeleteOldTextFiles() {
            if(ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                Console.Write("rm");
                foreach(var pltFile in dir.GetFiles("*.txt")) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
            }
        }

        /// <summary>
        /// This method saves the coordinates of a OptiLevelSet onto the respective coordinates of the step 
        /// </summary>
        /// <param name="step">the step</param>
        /// <param name="stepOptiLevelSet">the OptiLevelSet </param>
        public void SaveStepFieldToStep(double[] step, IOptiLevelSet stepOptiLevelSet) {
            for(int iCord = 0; iCord < +stepOptiLevelSet.GetLength(); iCord++) {
                if(Math.Abs(step[(int)ResidualMap.TotalLength + iCord]) > 1e-14)
                    step[(int)ResidualMap.TotalLength + iCord] = stepOptiLevelSet.GetParam(iCord);
            }
        }
        /// <summary>
        /// Helper method to Write stuff to lists and Console
        /// </summary>
        /// <param name="TimestepNo"></param>
        public void WriteStuffAndAddToLists(int TimestepNo) {
            if(Control.WriteLSCoordinates) {
                for(int i = 0; i < LevelSetOpti.GetLength(); i++)
                    Console.Write(LevelSetOpti.GetParamName(i) + ": " + Math.Round(LevelSetOpti.GetParam(i), 4) + ", ");
            }
            Console.Write($" Alpha = {string.Format("{0:#.##E+00}", m_alpha)}, ");
            Console.Write($"Gamma = {string.Format("{0:#.##E+00}", gamma)}, ");
            Console.Write($"Kappa = {string.Format("{0:#.##E+00}", kappa)}, ");
            Console.Write($"Mu = {string.Format("{0:#.##E+00}", mu)}, ");
            Console.Write($"QuadOrder={GetGlobalQuadOrder()}");
            Alphas.Add(m_alpha);
            Gammas.Add(gamma);
            Kappas.Add(kappa);
            Mus.Add(mu);
            ResNorms.Add(ResidualVector.MPI_L2Norm());
            obj_f_vals.Add(obj_f_vec.MPI_L2Norm());
            LevelSetOptiParams.Add(LevelSetOpti.GetParamsAsArray());
            StepCount.Add(TimestepNo);
        }
        /// <summary>
        /// Computes both Residuals using the fields, that are also plotted and saves them to the variables res_l2 and obj_f
        /// </summary>
        public (double _res_l2, double _obj_f, double _res_L2) ComputeResiduals() {
            ComputeResiduals(obj_f_vec, ResidualVector);
            return (ResidualVector.MPI_L2Norm(), obj_f_vec.MPI_L2Norm(), L2Norm(Residuals));
            //en_res_L2 = L2Norm(obj_f_fields);
        }
        /// <summary>
        /// Gives the accumulated L2-Norm of an array of XDGFields
        /// </summary>
        /// <param name="residuals"></param>
        /// <returns></returns>
        public double L2Norm(XDGField[] residuals) {
            return L2NormPerField(residuals).Sum();
        }

        /// <summary>
        /// Gives the respective L2-Norms of an array of XDGFields
        /// </summary>
        /// <param name="residuals"></param>
        /// <returns></returns>
        public double[] L2NormPerField(XDGField[] residuals) {
            double[] l2norms = new double[residuals.Length];
            for(int iField = 0; iField < residuals.Length; iField++) {
                l2norms[iField] = l2norms[iField] + residuals[iField].L2NormAllSpecies();
            }
            return l2norms;
        }

        private void SaveStepToField(double[] step_t, IOptiLevelSet OptiLevelSetStep) {
            int iField = 0;
            int jCell = 0;
            int nMode = 0;
            int length_r = (int)UnknownsMap.TotalLength; // needs to be modified if more than one field is simulated

            for(int stepIndex = 0; stepIndex < length_r; stepIndex++) {
                UnknownsMap.LocalFieldCoordinateIndex(stepIndex, out iField, out jCell, out nMode);
                originalStep[iField].CoordinateVector[jCell * originalStep[iField].Basis.DOFperSpeciesPerCell * LsTrk.TotalNoOfSpecies + nMode] = step_t[stepIndex];
            }
            if(CurrentAgglo > 0) {
                if(MultiphaseAgglomerator.TotalNumberOfAgglomerations > 0) {
                    MultiphaseAgglomerator.Extrapolate(new CoordinateMapping(originalStep));
                }
            }
            OptiLevelSetStep.ProjectOntoLevelSet(LevelSetStep);
        }

        public void Normalize() {
            if(LevelSetOpti.GetParam(0) != 0) {
                for(int i = 0; i < LevelSetOpti.GetLength(); i++) {
                    LevelSetOpti.SetParam(i, LevelSetOpti.GetParam(i) / LevelSetOpti.GetParam(1));
                }

                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            }

        }
        #endregion

        /// <summary>
        /// Helper routine to compute and Output Residuals
        /// </summary>
        public void ComputeAndWriteResiduals() {

            TransformFromAggToSourceSpace();

            (res_l2, obj_f, res_L2) = ComputeResiduals();
            var kappa_copy = kappa;
            kappa = 1;
            var f_phi_vec = ComputeFphi();
            kappa = kappa_copy;
             Console.Write($"l2: ||R0||= {string.Format("{0:#.#######E+00}", res_l2)}, f_err={string.Format("{0:#.######E+00}", obj_f)}, f_phi={string.Format("{0:#.######E+00}", f_phi_vec.InnerProd(f_phi_vec))},");

        }
        /// <summary>
        /// evaluates both residuals
        /// </summary>
        public void ComputeResiduals<V>(V enRes_out, V res_out) where V : IList<double> {
            Oproblem.EvalConsAndObj(enRes_out, res_out,ConservativeFields);
            //Transform back into Agglomerated Space
            if(CurrentAgglo > 0) {
                MultiphaseAgglomerator.ManipulateMatrixAndRHS(default(MsrMatrix), enRes_out, obj_f_map, ResidualMap);
                MultiphaseAgglomerator.ManipulateMatrixAndRHS(default(MsrMatrix), res_out, ResidualMap, ResidualMap);
            }

        }


        #region Sub Methods of Optimization Algorithm
        /// <summary>
        /// simply sets a one one the diagonal for zero rows
        /// </summary>
        /// <param name="Mat"></param>
        static public void SetDiagonaForZeroRows(MsrMatrix Mat)
        {
            //Console.WriteLine("Check: Computed System ...");
            for (int rows = 0; rows < Mat.NoOfRows; rows++)
            {
                if (Mat.GetNoOfNonZerosPerRow(rows) == 0)
                {
                    Mat[rows, rows] = 1;
                }
            }
        }
        /// <summary>
        /// simply sets a one one the diagonal for zero rows
        /// </summary>
        /// <param name="Mat"></param>
        static public void SetDiagonaForZeroRows(BlockMsrMatrix Mat)
        {
            //Console.WriteLine("Check: Computed System ...");
            for (int rows = 0; rows < Mat.NoOfRows; rows++)
            {
                if (Mat.GetNoOfNonZerosPerRow(rows) == 0)
                {
                    Mat[rows, rows] = 1;
                }
            }
        }
        /// <summary>
        /// Solves the Linear System
        /// </summary>
        public void SolveSystem() {

            using (new FuncTrace())
            {
                //Adds ones on the diagonal to ghost rows (rows coincing with empty or agglomerated phase-cells) of the left hand side
                for (int rows = 0; rows < LHS.NoOfRows; rows++)
                {
                    if (LHS.GetNoOfNonZerosPerRow(rows) == 0)
                    {
                        LHS[rows, rows] = 1;
                        if (Math.Abs(RHS[rows]) > 1e-32)
                        {
                            Console.WriteLine("Nonzero Entry in ghost RHS row " + rows + "!!!!!! FIX Please");
                            RHS[rows] = 0;
                            Jobj_U.CheckForNanOrInfM();
                        }
                    }
                }
#if DEBUG
            //RHS.CheckForNanOrInfV();
            //LHS.CheckForNanOrInfM();
            //double minEig;
            //double[] VEig = new double[RHS.Length];
            //(minEig, VEig)=LHS.MinimalEigen();
            //var conditionNumber = LHS.condestArnoldi();
#endif

                try
                {
                    //save mat files 

                    SimpleSolversInterface.Solve_Direct(LHS, stepIN, RHS);
                }
                catch (Exception e)
                {
#if DEBUG
                //this.PlotCurrentState(0.0, 666, this.Control.SuperSampling);
                //Console.WriteLine("LHSInf=" + LHS.InfNorm());
                //Console.WriteLine("dR0dUInf=" + Jr_U.InfNorm());
                //Console.WriteLine("dR0dPhiInf=" + Jr_phi.InfNorm());
                //Console.WriteLine("dR1dUInf=" + Jobj_U.InfNorm());
                //Console.WriteLine("dR1dPhiInf=" + Jobj_phi.InfNorm());
                //try {
                //    Console.WriteLine("LHSMinEIg=" + LHS.MinimalEigen());
                //    Console.WriteLine("LHSMaxEIg=" + LHS.MaximalEigen());
                //    Console.WriteLine("dR0dUMinEIg=" + Jr_U.MinimalEigen());
                //    Console.WriteLine("dR0dUMaximalEIg=" + Jr_U.MaximalEigen());
                //} catch {
                //    Console.WriteLine("Exception thrown From MinEig, but we continue...");
                //}
                //Console.WriteLine("RHSNorm=" + RHS.L2Norm());
                //Console.WriteLine("gradfNorm=" + gradf.L2Norm());
                //Console.WriteLine("residualNorm=" + ResidualVector.L2Norm());
#endif
                    Console.WriteLine("Exception: " + e.ToString() + " From DirectSOlver Thrown, but we continue...");
                }
                CmpLineSearchResidualWeight();
                CmpFphiWeight();
                
                //agglomerated the solution vector
                //TransformFromSourceToAggSpace();
                //and make a backup
                //createAgglomeratedBackUp();
                
            }
        }
        /// <summary>
        /// Computes the weight kappa associated with fphi (such that f=f_err+kappa*fphi)
        /// </summary>
        public void CmpFphiWeight() {
            var f_phi_vec = ComputeFphi();
            var f_phi = f_phi_vec.InnerProd(f_phi_vec);
            if(CurrentStepNo< Control.Kappa_M) {
                if(obj_f < Control.Kappa_Xi*kappa*kappa*f_phi) {
                    kappa = Math.Max(Control.Kappa_Min, Control.Kappa_v * kappa);
                }
            }
            Kappas.Add(kappa);
        }     
        /// <summary>
        /// Computes the weigth (mu) associated with the line-search merit function
        /// </summary>
        public void CmpLineSearchResidualWeight() {

            //version 1
            //stepUphi = new double[(int)UnknownsMap.TotalLength + LevelSetOpti.GetLength()];
            //var BtimesS = new double[(int)UnknownsMap.TotalLength + LevelSetOpti.GetLength()];
            //stepUphi.SetSubVector(stepIN, 0, (int)UnknownsMap.TotalLength + LevelSetOpti.GetLength());
            //var B = Jobj.Transpose() * Jobj;
            //B.SpMV(0.5, stepUphi, 0.0, BtimesS);

            //res_l1 = 0;
            //for(int i = 0; i < ResidualVector.Length; i++) {
            //    res_l1 += ResidualVector[i].Abs();
            //}

            //next_mu = (gradf.InnerProd(stepUphi) + stepUphi.InnerProd(BtimesS)) / (1 - Control.Mu_Rho) / res_l1;
            //next_mu = gradf.InnerProd(stepUphi) / (1 - Control.Mu_Rho) / res_l1;
            //mu = Math.Min(Control.Mu_Max, Math.Max(Control.Mu_Omega * next_mu, mu));

            //Version 2
            //mu = 2 * lambda.MaximumAbsolute();

            //Version 3
            double max = 0;
            var LamCoord = new CoordinateVector(LagMult);
            for (int i=0; i < (int)UnknownsMap.TotalLength; i++) {
                LamCoord[i] = stepIN[i + (int)UnknownsMap.TotalLength + LevelSetOpti.GetLength()];
                if(LamCoord[i].Abs() > max) {
                    max = LamCoord[i].Abs();
                }
            }
            mu = 2 * max;
        }
        /// <summary>
        /// Computes the gamma associated with regularization of the linear system LHS*(stepIN;lambda)=RHS
        /// </summary>
        public void CmpRegParam() {

            SaveStepFieldToStep(stepIN, GetStepOptiLevelSet(stepIN));
            //double norm = LevelSetStep.L2Norm();
            double[] levelSetStepCoordinates = stepIN.GetSubVector((int)UnknownsMap.TotalLength, LevelSetOpti.GetLength());
            double norm = LevelSetOpti.Norm(levelSetStepCoordinates);


            if(norm < Control.sigma_1 * Control.L) {
                gamma = Math.Max(gamma / Control.tauGamma, Control.Gamma_Min);
            } else if(norm > Control.sigma_2 * Control.L) {
                gamma = Math.Max(gamma * Control.tauGamma, Control.Gamma_Min);
            } else {
                gamma = Math.Max(gamma, Control.Gamma_Min);
            }
            gamma = Math.Min(gamma, Control.Gamma_Max);

        }
#endregion
        
#region Plotting
        /// <summary>
        /// A public plot function (can't make Application member PlotCurrentState public)
        /// </summary>
        /// <param name="physTime"></param>
        /// <param name="timestepNo"></param>
        /// <param name="superSampling"></param>
        public void Plot(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            PlotCurrentState(physTime, timestepNo, superSampling);
        }
        /// <summary>
        /// a custom Plot-function (ignoring the physical time)
        /// </summary>
        /// <param name="timestepNo"></param>
        /// <param name="superSampling"></param>
        public void PlotCurrentState(TimestepNumber timestepNo, int superSampling = 0) {
            if(plotDriver == null) {
                plotDriver = new Tecplot(GridData, true, false, (uint)superSampling);
            }
            UpdateDerivedVariables();
            plotDriver.PlotFields(Control.ProjectName + "_" + timestepNo, 0.0, m_IOFields);
        }
        /// <summary>
        /// a custom Plot-function (ignoring the physical time), where the user sets a custom name
        /// </summary>
        /// <param name="name"></param>
        /// <param name="superSampling"></param>
        public void PlotCurrentState(string name, int superSampling = 0) {
            if(plotDriver == null) {
                plotDriver = new Tecplot(GridData, true, false, (uint)superSampling);
            }
            UpdateDerivedVariables();
            plotDriver.PlotFields(name, 0.0, m_IOFields);
        }
        /// <summary>
        /// Another plotting function with more options
        /// </summary>
        /// <param name="physTime"></param>
        /// <param name="timestepNo"></param>
        /// <param name="superSampling"></param>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            using (new FuncTrace())
            {
                if (plotDriver == null)
                {
                    plotDriver = new Tecplot(GridData, true, false, (uint)superSampling);
                }
                UpdateDerivedVariables();
                plotDriver.PlotFields(Control.ProjectName + "_" + timestepNo, physTime, m_IOFields);

                #region Debugging: write residuals to text file
                // Sample points
                //int noOfPoints = 16;
                //double[] nodes = GenericBlas.Linspace(-1.8, -0.3, noOfPoints);
                //MultidimensionalArray points = MultidimensionalArray.Create(noOfPoints, 2);
                //for (int i = 0; i < noOfPoints; i++) {
                //    points[i, 0] = nodes[i];
                //    points[i, 1] = -0.65;
                //}

                //// FieldEvaluation
                //MultidimensionalArray results = MultidimensionalArray.Create(noOfPoints, Residuals.Length);
                //for (int i = 0; i < Residuals.Length; i++) {
                //    FieldEvaluation fieldEvaluator = new FieldEvaluation((GridData)this.GridData);
                //    fieldEvaluator.Evaluate(1.0, Residuals, points, 0.0, results);
                //}

                //// StreamWriter
                //using (System.IO.StreamWriter sw = new System.IO.StreamWriter(String.Format("Residuals{0}.txt", timestepNo))) {
                //    //Console.WriteLine("x \t y \t result");
                //    sw.WriteLine("x \t y \t rho \t xMom \t yMom \t rhoE");
                //    string resultLine;
                //    for (int i = 0; i < noOfPoints; i++) {
                //        resultLine = points[i, 0] + "\t" + points[i, 1] + "\t" + results[i, 0] + "\t" + results[i, 1] + "\t" + results[i, 2] + "\t" + results[i, 3] + "\t";
                //        //Console.WriteLine(resultLine);
                //        sw.WriteLine(resultLine);
                //    }
                //    sw.Flush();
                //}
                #endregion

                #region Debugging: write DG fields to text file
                //// Sample points
                //int noOfPoints = 16;
                //double[] nodes = GenericBlas.Linspace(-1.8, -0.3, noOfPoints);
                //MultidimensionalArray points = MultidimensionalArray.Create(noOfPoints, 2);
                //for (int i = 0; i < noOfPoints; i++) {
                //    points[i, 0] = nodes[i];
                //    points[i, 1] = -0.65;
                //}

                //// FieldEvaluation
                //MultidimensionalArray resultsFields = MultidimensionalArray.Create(noOfPoints, ConservativeFields.Length);
                //for (int i = 0; i < ConservativeFields.Length; i++) {
                //    FieldEvaluation fieldEvaluator = new FieldEvaluation((GridData)this.GridData);
                //    fieldEvaluator.Evaluate(1.0, this.ConservativeFields, points, 0.0, resultsFields);
                //}

                //// StreamWriter
                //using (System.IO.StreamWriter sw = new System.IO.StreamWriter(String.Format("DGFields{0}.txt", timestepNo))) {
                //    //Console.WriteLine("x \t y \t result");
                //    //sw.WriteLine("x \t y \t rho \t xMom \t yMom \t rhoE");
                //    string resultLine;
                //    for (int i = 0; i < noOfPoints; i++) {
                //        resultLine = points[i, 0] + "\t" + points[i, 1] + "\t" + resultsFields[i, 0] + "\t" + resultsFields[i, 1] + "\t" + resultsFields[i, 2] + "\t" + resultsFields[i, 3] + "\t";
                //        //Console.WriteLine(resultLine);
                //        sw.WriteLine(resultLine);
                //    }
                //    sw.Flush();
                //}
                //Density.CoordinateVector.SaveToTextFile(String.Format("DensityCoordVec_{0}.txt", timestepNo));
                //Energy.CoordinateVector.SaveToTextFile(String.Format("EnergyCoordVec_{0}.txt", timestepNo));
                #endregion
            }
        }

#endregion


        /// <summary>
        /// saves everything relevant to database
        /// </summary>
        /// <param name="timestepno">time step number</param>
        /// <param name="t">physical time</param>
        protected override ITimestepInfo SaveToDatabase(TimestepNumber timestepno, double t) {
            using (new FuncTrace())
            {
                // Make sure that all derived variables are updated before saving
                UpdateDerivedVariables();
                return base.SaveToDatabase(timestepno, t);
            }
        }

        /// <summary>
        /// This function computes a Jacobian of an Operator with respect to the LevelSet Coordinates using Central FDs
        /// </summary>
        /// <param name="Eval"> Operator to be differentiated </param>
        public (MsrMatrix Jr, MsrMatrix Jobj) FD_LevelSet() {
            using(new FuncTrace()) {

                int nCol_obj = Oproblem.GetObjLength(ConservativeFields);
                int nCol_con = new CoordinateVector(ConservativeFields).Count;
                int nRow = LevelSetOpti.GetLength();

                MsrMatrix Jobj = new MsrMatrix(nCol_obj, nRow, 1, 1);
                MsrMatrix Jr = new MsrMatrix(nCol_con, nRow, 1, 1);

                var r_vec = new double[nCol_con];
                var r_vec_eps = new double[nCol_con];
                var r_vec_eps2 = new double[nCol_con];

                var obj_vec = new double[nCol_obj];
                var obj_vec_eps = new double[nCol_obj];
                var obj_vec_eps2 = new double[nCol_obj];


                double x;
                double val;
                double dx_right;
                double dx_left;

                // epsilon
                double eps = 1.0e-8;
                //project OptiLevelSet onto XDGLevelSet
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                LsTrk.UpdateTracker(CurrentStepNo);
                LevelSet phi0backup = new LevelSet(new Basis(LsTBO.GridDat.Grid, LsTBO.Basis.Degree), "LevelSetbackup");
                phi0backup.CopyFrom(LsTBO);

                //Evaluation of unperturbed
                Oproblem.EvalConsAndObj(obj_vec,r_vec,ConservativeFields);


                for(int n_param = 0; n_param < nRow; n_param++) {

                    //skip loop if param is non changeable
                    if (Control.PartiallyFixLevelSetForSpaceTime && LevelSetOpti is SplineOptiLevelSet splineLS)
                    {
                        double yMin = splineLS.y.Min(); //lower boundary of space time domain
                        if (n_param < splineLS.y.Length && Math.Abs(yMin - splineLS.y[n_param]) < 1e-14) //only accumalte if DOF if it is not on lower time boundary (here yMin=tMin)
                        {
                            continue;
                        }

                    }

                    //compute distortions
                    x = LevelSetOpti.GetParam(n_param);
                    dx_right = x + eps / 2;
                    dx_left = x - eps / 2;


                    try
                    {
                        // apply right distortion
                        LevelSetOpti.SetParam(n_param, dx_right);
                        //project
                        LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                        LsTrk.UpdateTracker(CurrentStepNo);
                        //compute Objective
                        Oproblem.EvalConsAndObj(obj_vec_eps, r_vec_eps, ConservativeFields);

                        // do the same on the left
                        LevelSetOpti.SetParam(n_param, dx_left);
                        LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                        LsTrk.UpdateTracker(CurrentStepNo);
                        Oproblem.EvalConsAndObj(obj_vec_eps2, r_vec_eps2, ConservativeFields);
                    }
                    catch
                    {
                        //reset
                        LevelSetOpti.SetParam(n_param, x);
                        LsTBO.CopyFrom(phi0backup);
                        LsTrk.UpdateTracker(CurrentStepNo);
                        //throw new Exception("Error in FD Computation");
                    }
                    
                    //Write value into the matrix
                    for(int iRow = 0; iRow < nCol_obj; iRow++) {
                        val = (obj_vec_eps[iRow] - obj_vec_eps2[iRow]) / eps;
                        if(val != 0 && !val.IsNaNorInf()) {
                            Jobj[iRow, n_param] = val;
                        }
                        
                    }
                    //IMutuableMatrixEx_Extensions.SaveToTextFile(Jobj,"Jobj.txt");
                    for(int iRow = 0; iRow < nCol_con; iRow++) {
                        val = (r_vec_eps[iRow] - r_vec_eps2[iRow]) / eps;
                        if(val != 0 && !val.IsNaNorInf()) {
                            Jr[iRow, n_param] = val;
                        }
                    }


                }
                return (Jr,Jobj);
            }
        }
        /// <summary>
        /// Creates an IOptiLevelSet for the step (\Delta \varphi), for potential plotting
        /// </summary>
        /// <param name="step"></param>
        /// <returns></returns>
        public IOptiLevelSet GetStepOptiLevelSet(double[] step) {
            var ret_OLS = LevelSetOpti.CloneAs();
            for(int iCord = 0; iCord < +ret_OLS.GetLength(); iCord++) {
                ret_OLS.SetParam(iCord, step[(int)ResidualMap.TotalLength + iCord]);
            }
            return ret_OLS;
        }
        /// <summary>
        /// does what it promises. Making the Level Set LSTBO continuous, Not really in Use
        /// </summary>
        /// <param name="targetOptiLevelSet"></param>
        public void MakeLevelSetContinous(IOptiLevelSet targetOptiLevelSet) {
            // Does a continuity projection (hopefully)
            IGridData gridData = targetOptiLevelSet.GetGrid();
            SinglePhaseField continuousLevelSet = new SinglePhaseField(new Basis(gridData, targetOptiLevelSet.GetDegree()), "LevelSet_cont");
            var tmp_LsTrk = new LevelSetTracker((GridData)gridData, Control.CutCellQuadratureType, 1, LsTrk.SpeciesTable, new LevelSet[] { targetOptiLevelSet.ToLevelSet(targetOptiLevelSet.GetDegree()) });
            CellMask nearBandMask;
            if(Control.IsTwoLevelSetRun) {
                nearBandMask = tmp_LsTrk.Regions.GetNearMask4LevSet(1, 1);
            } else {
                nearBandMask = tmp_LsTrk.Regions.GetNearMask4LevSet(0, 1);
            }

            var bitMask = nearBandMask.GetBitMask();
            BitArray bitMaskPos = bitMask.CloneAs();
            for(int i = 0; i < bitMask.Count; i++) {
                if(targetOptiLevelSet.GetMeanValue(i) > 0 && !bitMask[i]) {
                    bitMaskPos[i] = true;
                }
            }
            CellMask posMask = new CellMask(gridData, bitMaskPos);
            CellMask nearMask = new CellMask(gridData, bitMask);
            ContinuityProjection = new ContinuityProjection(continuousLevelSet.Basis, new Basis(targetOptiLevelSet.GetGrid(), targetOptiLevelSet.GetDegree()), gridData, ContinuityProjectionOption.ConstrainedDG);
            ContinuityProjection.MakeContinuous(targetOptiLevelSet.ToSinglePhaseField(targetOptiLevelSet.GetDegree()), continuousLevelSet, nearMask, posMask);
            targetOptiLevelSet.ProjectFromLevelSet(continuousLevelSet);
            //targetOptiLevelSet.ProjectOntoLevelSet(LsTBO);
        }

        ///// <summary>
        ///// New Version of Function 31.07.2023 
        ///// Here Agglomeration is first applied to the products of the Matrices, working worse then old version.
        ///// </summary>

        //public void ComputeSystem_new() {
        //    LHS.Clear();
        //    RHS.Clear();
        //    Jobj_U.Clear();
        //    Jr_U.Clear();
        //    Jobj.Clear();
        //    Jr.Clear();
        //    Jobj_phi.Clear();
        //    Jr_phi.Clear();
        //    gradf.Clear();
        //    gradf_U.Clear();
        //    lambda.Clear();

            
        //    int length_r = (int)UnknownsMap.TotalLength;
        //    int length_R = (int)obj_f_map.TotalLength;
        //    int length_l = LevelSetOpti.GetLength();

        //    int noOfRowsB = length_r + length_l;

        //    //***************************** */
        //    //// Compute Jr_U and Jobj_U
        //    //***************************** */
        //    (Jobj_U,Jr_U) = Oproblem.GetJacobians(ConservativeFields, XSpatialOperator.LinearizationHint);
        //    Oproblem.EvalConsAndObj(obj_f_vec, ResidualVector,ConservativeFields);

        //    //***************************** */
        //    //// Compute Jobj_phi & Jr_phi
        //    //***************************** */
        //    //obj_f_vec.SaveToTextFile("obj_f_vec_beforeFD.txt");
        //    //LsTBO.CoordinateVector.SaveToTextFile(@"C:\Users\sebastian\Documents\Forschung" + $"\\LSCoord_p{Control.SolDegree}_TS{CurrentStepNo}.txt");
        //    (Jr_phi, Jobj_phi) = FD_LevelSet();

        //    //Compute some B parts
        //    var BUU = Jobj_U.Transpose() * Jobj_U;
        //    var BUPhi = Jobj_U.Transpose() * Jobj_phi;

        //    //Compute gradf_U
        //    Jobj_U.Transpose().SpMV(1, obj_f_vec, 0, gradf_U);

        //    //***************************** */
        //    //// Agglomerate Jobj,Jr and obj_vec,r before multiplication
        //    //***************************** */
        //    if(CurrentAgglo > 0.0) {
        //        if(MultiphaseAgglomerator.TotalNumberOfAgglomerations > 0) {
        //            // Get the current MultigridOperator
        //            var MgSeq = CoarseningAlgorithms.CreateSequence(GridData);
        //            AggregationGridBasis[][] basisSeq = AggregationGridBasis.CreateSequence(MgSeq, obj_f_map.BasisS);
        //            basisSeq.UpdateXdgAggregationBasis(MultiphaseAgglomerator);
        //            var MassMatFactory = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, GetGlobalQuadOrder()).MassMatrixFactory;
        //            var MassMat = MassMatFactory.GetMassMatrix(ResidualMap, inverse: false);
        //            MultiphaseAgglomerator.ManipulateMatrixAndRHS(MassMat, default(double[]), ResidualMap, ResidualMap);
        //            var dummy = new BlockMsrMatrix(ResidualMap, ResidualMap);
        //            multOp = new MultigridOperator(basisSeq, ResidualMap, dummy, MassMat, MultiGridOperatorConfig, XSpatialOperator);

        //            // agglomeration for Jr_U, Jobj_U and residuals is done using the standard BoSSS function
        //            MultiphaseAgglomerator.ManipulateMatrixAndRHS(Jr_U, ResidualVector, new CoordinateMapping(ConservativeFields), new CoordinateMapping(ConservativeFields));
        //            MultiphaseAgglomerator.ManipulateMatrixAndRHS(BUU, gradf_U, new CoordinateMapping(ConservativeFields), new CoordinateMapping(ConservativeFields));

        //            // helper function to do the agglomeration for phi part 
        //            void AccRowManipulationMatrixForSpeciesAndFields(SpeciesId spc, XDGField[] fields, BlockMsrMatrix inputMatrix) {
        //                var AggRphi = MultiphaseAgglomerator.GetAgglomerator(spc);
        //                CellMask spcMask = LsTrk.Regions.GetSpeciesMask(spc);
        //                MiniMapping rowMini = new MiniMapping(new CoordinateMapping(fields), spc, LsTrk.Regions);
        //                inputMatrix.Acc(1.0, AggRphi.GetRowManipulationMatrix(new CoordinateMapping(fields), rowMini.MaxDeg, rowMini.NoOfVars, rowMini.i0Func, rowMini.NFunc, false, spcMask));
        //            }
        //            MsrMatrix AggMatrix;
        //            // ****Perform custom Agglomeration for BUPhi**
        //            LeftMul = new BlockMsrMatrix(new CoordinateMapping(ConservativeFields), new CoordinateMapping(ConservativeFields));
        //            foreach(SpeciesId spc in SpeciesToEvaluate_Ids) {
        //                AccRowManipulationMatrixForSpeciesAndFields(spc, ConservativeFields, LeftMul);
        //            }
                    
        //            //multiply
        //            AggMatrix = MsrMatrix.Multiply(LeftMul.ToMsrMatrix(), BUPhi);
        //            BUPhi = AggMatrix.CloneAs();


        //            // ****Perform custom Aggregation for Jr_phi***
        //            //multiply
        //            AggMatrix = MsrMatrix.Multiply(LeftMul.ToMsrMatrix(), Jr_phi);
        //            Jr_phi = AggMatrix.CloneAs();
        //        }
        //    }
        //    Assembler assembler = new Assembler(noOfRowsB + length_R + length_r);

        //    //*******************************/
        //    //// assemble R_{phi} (part of objective function corresponding to Level Set)
        //    ///*******************************/
        //    (double[] f_phi_vec, MsrMatrix Jf_phi) = ComputeFphiandJFphi();

        //    //***************************** */
        //    //// Compute BPhiPhi
        //    //***************************** */
        //    var BPhiPhi = Jobj_phi.Transpose() * Jobj_phi;
        //    //Add Fphi contribution
        //    assembler.AccBlock(Jf_phi.Transpose() * Jf_phi, BPhiPhi, 0, 0);
        //    //Add Regularization
        //    MsrMatrix D = this.LevelSetOpti.GetRegMatrix();
        //    D.Scale(gamma);
        //    assembler.AccBlock(D, BPhiPhi, 0, 0);

        //    //***************************** */
        //    //// Set the Blocks of B
        //    ///   B = (  BUU  ,     BUPhi^T 
        //    ///         BUPhi , BPhiPhi  )
        //    //***************************** */
        //    var B = new MsrMatrix( noOfRowsB, noOfRowsB,1, 1);
        //    assembler.AccBlock(BUU, B, 0, 0);
        //    assembler.AccBlock(BUPhi, B, 0,length_r);
        //    assembler.AccBlock(BUPhi.Transpose(), B, length_r,0);
        //    assembler.AccBlock(BPhiPhi, B, length_r, length_r);

        //    //***************************** */
        //    //// Set the Blocks of Jr
        //    ///  Jr = ( Jr_U  , Jr_Phi )
        //    //***************************** */
        //    assembler.AccBlock(Jr_U, Jr, 0, 0);
        //    assembler.AccBlock(Jr_phi, Jr, 0, length_r);

        //    //***************************** */
        //    //// Set the Blocks of LHS 
        //    ///  LHS = ( B  , Jr^T 
        //    ///          Jr ,   0   )
        //    //***************************** */
        //    //Add B
        //    assembler.AccBlock(B, LHS, 0, 0);
        //    //upper right block
        //    assembler.AccBlock(Jr.Transpose(), LHS, 0, noOfRowsB);
        //    //down left block
        //    assembler.AccBlock(Jr, LHS, noOfRowsB, 0);



        //    assembler.AccBlock(Jobj_U, Jobj, 0, 0);
        //    assembler.AccBlock(Jobj_phi, Jobj, 0, length_r);

        //    //IMutuableMatrixEx_Extensions.SaveToTextFile(J_R_comp, path + "\\J_R_comp_TT.txt");

        //    //***************************** */
        //    //// Set the Blocks of <see:cref and Assemble RHS
        //    ///  gradf = (              JRerr_U^T R_err(U)
        //    ///            JRerr_Phi^T R_err(U) + JfPhi_Phi^T R_phi(U)
        //    ///                                 r                       )          
        //    //***************************** */
        //    for(int i = 0; i < length_r; i++) {
        //        gradf[i] = gradf_U[i];
        //        RHS[i] = gradf[i];
        //    }

        //    //Compute JRerr_Phi^T R_err(U)
        //    double[] gradRerrPhi = new double[length_l];
        //    Jobj_phi.Transpose().SpMV(1, obj_f_vec, 0, gradRerrPhi);

        //    //Compute JfPhi_Phi^T R_phi(U)
        //    double[] gradfphi = new double[length_l];
        //    Jf_phi.Transpose().SpMV(1,f_phi_vec, 0, gradfphi);


        //    for(int i = length_r; i < gradf.Length; i++) {
        //        gradf[i] = gradRerrPhi[i - length_r] +gradfphi[i - length_r];
        //        RHS[i] = gradf[i];
        //    }

        //    //gradf_U.SetSubVector(gradf, 0, length_r);

        //    //***************************** */
        //    //// Assemble RHS
        //    //***************************** */
        //    double[] res = ResidualVector.ToArray();
        //    RHS.SetSubVector(res, gradf.Length, length_r);
        //    RHS.ScaleV(-1.0);
        //}
        /// <summary>
        /// Obtain the Block Jacobi Preconditioner of a BlockMsrMatrix
        /// </summary>
        /// <param name="M"></param>
        /// <param name="uMap"></param>
        /// <returns></returns>
        public static BlockMsrMatrix GetBlockJacobi(BlockMsrMatrix M,CoordinateMapping uMap)
        {
            var Diag = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
                   
            int Jloc = uMap.LocalNoOfBlocks;
            long j0 = uMap.FirstBlock;
            MultidimensionalArray temp = null;
            for (int j = 0; j < Jloc; j++)
            {
                long jBlock = j + j0;
                int Nblk = uMap.GetBlockLen(jBlock);
                long i0 = uMap.GetBlockI0(jBlock);

                if (temp == null || temp.NoOfCols != Nblk)
                    temp = MultidimensionalArray.Create(Nblk, Nblk);

                M.ReadBlock(i0, i0, temp);
                Diag.AccBlock(i0, i0, 1.0, temp, 0.0);
            }
            return Diag;
        }
        public void AgglomerateJacobiansAndResiduals()
        {
            using (new FuncTrace())
            {
                if (CurrentAgglo > 0.0)
                {
                    if (MultiphaseAgglomerator.TotalNumberOfAgglomerations > 0)
                    {
                        // Get the current MultigridOperator
                        var MgSeq = CoarseningAlgorithms.CreateSequence(GridData);
                        AggregationGridBasis[][] basisSeq = AggregationGridBasis.CreateSequence(MgSeq, obj_f_map.BasisS);
                        basisSeq.UpdateXdgAggregationBasis(MultiphaseAgglomerator);
                        var MassMatFactory = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, GetGlobalQuadOrder()).MassMatrixFactory;
                        var MassMat = MassMatFactory.GetMassMatrix(ResidualMap, inverse: false);
                        MultiphaseAgglomerator.ManipulateMatrixAndRHS(MassMat, default(double[]), ResidualMap, ResidualMap);
                        //XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = XSpatialOperator.GetMatrixBuilder(LsTrk, ResidualMap, null, obj_f_map);
                        //var OperatorMatrix = new BlockMsrMatrix(ResidualMap, obj_f_map);
                        //double[] Affine = new double[OperatorMatrix.RowPartitioning.LocalLength];
                        //mtxBuilder.ComputeMatrix(OperatorMatrix, Affine);
                        var dummy = new BlockMsrMatrix(ResidualMap, ResidualMap);
                        multOp = new MultigridOperator(basisSeq, ResidualMap, dummy, MassMat, MultiGridOperatorConfig, XSpatialOperator);
                        // agglomeration for Jr_U, Jobj_U and residuals is done using the standard BoSSS function
                        MultiphaseAgglomerator.ManipulateMatrixAndRHS(Jr_U, ResidualVector, new CoordinateMapping(ConservativeFields), new CoordinateMapping(ConservativeFields));
                        MultiphaseAgglomerator.ManipulateMatrixAndRHS(Jobj_U, obj_f_vec, new CoordinateMapping(obj_f_fields), new CoordinateMapping(ConservativeFields));
                        // helper function to do the agglomeration
                        void AccRowManipulationMatrixForSpeciesAndFields(SpeciesId spc, XDGField[] fields, BlockMsrMatrix inputMatrix)
                        {
                            var AggRphi = MultiphaseAgglomerator.GetAgglomerator(spc);
                            CellMask spcMask = LsTrk.Regions.GetSpeciesMask(spc);
                            MiniMapping rowMini = new MiniMapping(new CoordinateMapping(fields), spc, LsTrk.Regions);
                            inputMatrix.Acc(1.0, AggRphi.GetRowManipulationMatrix(new CoordinateMapping(fields), rowMini.MaxDeg, rowMini.NoOfVars, rowMini.i0Func, rowMini.NFunc, false, spcMask));
                        }


                        // ****Perform custom Aggregation for Jobj_phi***
                        LeftMul = new BlockMsrMatrix(new CoordinateMapping(obj_f_fields), new CoordinateMapping(obj_f_fields));
                        foreach (SpeciesId spc in SpeciesToEvaluate_Ids)
                        {
                            AccRowManipulationMatrixForSpeciesAndFields(spc, obj_f_fields, LeftMul);
                        }
                        MsrMatrix AggMatrix;
                        //multiply
                        AggMatrix = MsrMatrix.Multiply(LeftMul.ToMsrMatrix(), Jobj_phi);
                        Jobj_phi = AggMatrix.CloneAs();



                        // ****Perform custom Aggregation for Jr_phi***
                        LeftMul = new BlockMsrMatrix(new CoordinateMapping(ConservativeFields), new CoordinateMapping(ConservativeFields));
                        foreach (SpeciesId spc in SpeciesToEvaluate_Ids)
                        {
                            AccRowManipulationMatrixForSpeciesAndFields(spc, ConservativeFields, LeftMul);
                        }
                        //multiply
                        AggMatrix = MsrMatrix.Multiply(LeftMul.ToMsrMatrix(), Jr_phi);
                        Jr_phi = AggMatrix.CloneAs();
                    }
                }
            }
        }
        /// <summary>
        /// Assembles Left- and Right Hand side 
        /// To do so 
        ///  1. the Jacobians Jr_U,Jobj_U are computed with respect to the state DOFS -depending  on the Linearization chosen
        ///  2. the Jacobians Jobj_vec,Jobj_phi are computed with respect to the Level Set DOFS using the FD_LevelSet function
        ///  3. If agglomerate=true, agglomeration is performed on the Jacobians
        ///  4. everything is put together in the LHS variable
        ///  5. grad_f is computed and put together in the RHS variable
        /// </summary>
        public void ComputeSystem() {
            using (var f=new FuncTrace())
            {
                LHS.Clear();
                RHS.Clear();
                Jobj_U.Clear();
                Jr_U.Clear();
                Jobj.Clear();
                Jr.Clear();
                Jobj_phi.Clear();
                Jr_phi.Clear();
                gradf.Clear();
                gradf_U.Clear();
                lambda.Clear();


                int length_r = (int)UnknownsMap.TotalLength;
                int length_R = (int)obj_f_map.TotalLength;

                //***************************** */
                //// Compute Jr_U and Jobj_U
                //***************************** */
                (Jobj_U, Jr_U) = Oproblem.GetJacobians(ConservativeFields, XSpatialOperator.LinearizationHint);
                Oproblem.EvalConsAndObj(obj_f_vec, ResidualVector, ConservativeFields);
                res_l2 = ResidualVector.MPI_L2Norm();

                resetStep();
                (res_l2, obj_f, res_L2) = ComputeResiduals();

                //***************************** */
                //// Compute Jobj_phi & Jr_phi
                //***************************** */
                //obj_f_vec.SaveToTextFile("obj_f_vec_beforeFD.txt");
                //LsTBO.CoordinateVector.SaveToTextFile(@"C:\Users\sebastian\Documents\Forschung" + $"\\LSCoord_p{Control.SolDegree}_TS{CurrentStepNo}.txt");
                try
                {
                    (Jr_phi, Jobj_phi) = FD_LevelSet();
                }
                catch
                {
                    Console.WriteLine("Error in FD Level Set");
                    //AllthePossibleStepsPlot();
                    //throw new Exception("Error in FD Computation");
                }


                //***************************** */
                //// Agglomerate Jobj,Jr and obj_vec,r before multiplication
                //***************************** */
                AgglomerateJacobiansAndResiduals();

                //*******************************/
                //// assemble R_{phi} (part of objective function corresponding to Level Set)
                ///*******************************/
                (double[] f_phi_vec, MsrMatrix Jf_phi) = ComputeFphiandJFphi();
                
                using (new BlockTrace("Blocks Of LHS",f))
                {
                    //***************************** */
                    //// Set the Blocks of LHS
                    //***************************** */
                    int length_l = (int)Jobj_phi.NoOfCols;
                    int noOfRowsB = length_r + length_l;
                    // get helper object for assembling
                    MatrixAssembler assembler = new MatrixAssembler(noOfRowsB + length_R + length_r);
                    //Assemble J_r
                    assembler.AccBlock(Jr_U, Jr, 0, 0);
                    assembler.AccBlock(Jr_phi, Jr, 0, length_r);
                    //upper right block
                    assembler.AccBlock(Jr.Transpose(), LHS, 0, noOfRowsB);
                    //down left block
                    assembler.AccBlock(Jr, LHS, noOfRowsB, 0);
                    //assemble Jobj first
                    assembler.AccBlock(Jobj_U, Jobj, 0, 0);
                    assembler.AccBlock(Jobj_phi, Jobj, 0, length_r);
                    //upper left block
                    assembler.AccBlock(Jobj.Transpose() * Jobj, LHS, 0, 0);
                    //Regularization
                    MsrMatrix reg = this.LevelSetOpti.GetRegMatrix();
                    reg.Scale(gamma);
                    assembler.AccBlock(reg, LHS, length_r, length_r);
                    //Fphi contribution
                    assembler.AccBlock(Jf_phi.Transpose() * Jf_phi, LHS, length_r, length_r);

                }


                //***************************** */
                //// compute RHS
                //***************************** */
                using (new BlockTrace("RHS compuation", f))
                {
                    int length_l = (int)Jobj_phi.NoOfCols;
                    //calculate Grad f
                    Jobj.Transpose().SpMV(1, obj_f_vec, 0, gradf);

                    // add f_phi
                    double[] gradfphi = new double[length_l];
                    Jf_phi.Transpose().SpMV(1, f_phi_vec, 0, gradfphi);

                    // add gradient of f
                    for (int i = length_r; i < gradf.Length; i++)
                    {
                        gradf[i] = gradf[i] + gradfphi[i - length_r];
                    }
                    gradf_U.SetSubVector(gradf, 0, length_r);
                    for (int i = 0; i < gradf.Length; i++)
                    {
                        RHS[i] = gradf[i];
                    }

                    // add residual r(z)
                    RHS.SetSubVector(ResidualVector.ToArray(), gradf.Length, length_r);
                    RHS.ScaleV(-1.0);

                    //compute jacobian times residual (needed in line search)
                    JacR0TimesR0 = new double[gradf.Length];
                    Jr.Transpose().SpMV(1, ResidualVector, 0, JacR0TimesR0);
                }

                if (Control.SaveMatrices)
                {
                    // Save for preconditioner analysis
                    var uMap = new CoordinateMapping(ConservativeFields);

                    //create directory
                    var folder = "matrices";
                    if (!Directory.Exists(@".\" + folder))
                    {
                        Directory.CreateDirectory(@".\" + folder);
                    }
                    folder = @".\" + folder + @"\";

                    // get BlockMsrMatrix Version of Jacobian
                    var dR0dU = Jr_U.ToBlockMsrMatrix(uMap, uMap);
                    var dmy = new double[uMap.TotalLength];

                    //Save Matrices and RHS
                    Jr_phi.SaveToTextFileSparse(folder + "dR0dPhi" + CurrentStepNo + ".txt");
                    Jobj_phi.SaveToTextFileSparse(folder + "dR1dPhi" + CurrentStepNo + ".txt");
                    Jobj_U.SaveToTextFileSparse(folder + "dR1dU" + CurrentStepNo + ".txt");
                    RHS.SaveToTextFile(folder + "RHS" + CurrentStepNo + ".txt");
                    dR0dU.SaveToTextFileSparse(folder + "dR0dU" + CurrentStepNo + ".txt");

                    //get Ilu and save
                    var dR0dU_BILU = CellILU.ComputeILU(0, dR0dU.ToBlockMsrMatrix(uMap, uMap));
                    var dR0dU_BJ = GetBlockJacobi(dR0dU, uMap);
                    var dR0dU_BILU_MSR = dR0dU_BILU.ToMsrMatrix(); var dR0dU_BJ_MSR = dR0dU_BJ.ToMsrMatrix();
                    dR0dU_BILU_MSR.SaveToTextFileSparse(folder + "dR0dU_BILU" + CurrentStepNo + ".txt");
                    dR0dU_BJ_MSR.SaveToTextFileSparse(folder + "dR0dU_BJ" + CurrentStepNo + ".txt");

                    //P0 projection and Restriction operator
                    var incOp = GetInclusionOperator(ConservativeFields, 0);
                    incOp.SaveToTextFileSparse(folder + "Q0" + CurrentStepNo + ".txt");
                    var projOp = GetProjectionOperator(LsTrk, ConservativeFields, GetGlobalQuadOrder(), 0);
                    projOp.SaveToTextFileSparse(folder + "P0" + CurrentStepNo + ".txt");

                    //save some info about the matrices
                    var info_keys = new string[] { "p", "nCells", "gamma", };
                    var pDeg = ConservativeFields[0].Basis.Degree;
                    var nCells = ConservativeFields[0].Basis.GridDat.iGeomCells.Count;
                    var gamma = Gammas.Last();
                    var info = new double[] { pDeg, nCells, gamma };
                    info.SaveToTextFile(folder + "info" + CurrentStepNo + ".txt");
                    info_keys.SaveToTextFile(folder + "info_keys" + CurrentStepNo + ".txt");
                }
            }
        }
        public (double[] f_phi_vec, MsrMatrix Jf_phi_vec) ComputeFphiandJFphi() {
            using(new FuncTrace()) {

                int length_r = (int)UnknownsMap.TotalLength;
                int length_l = LevelSetOpti.GetLength();
                // compute unperturbed fphi
                var f_phi_vec = ComputeFphi();
                
                double[] f_phi_vec_eps = new double[f_phi_vec.Length];
                double[] f_phi_vec_eps2 = new double[f_phi_vec.Length];
                var Jf_phi_vec = new MsrMatrix(f_phi_vec.Length, length_l, 1, 1);

                double x;
                double val;
                double dx_right;
                double dx_left;

                double[] Epsilons = new double[length_l];
                for(int i = 0; i < length_l; i++) {
                    Epsilons[i] = 1.0e-8;
                }
                // create Backup of OptiLevelSet an project OptiLevelSet onto LevelSet
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                LevelSet phi0backup = new LevelSet(new Basis(LsTBO.GridDat.Grid, LsTBO.Basis.Degree), "LevelSetbackup");
                phi0backup.CopyFrom(LsTBO);


                for(int n_param = 0; n_param < length_l; n_param++) {

                    //compute distortions
                    x = LevelSetOpti.GetParam(n_param);
                    dx_right = x + Epsilons[n_param] / 2;
                    dx_left = x - Epsilons[n_param] / 2;


                    // apply right distortion
                    LevelSetOpti.SetParam(n_param, dx_right);
                    //project
                    LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                    //compute fphi
                    f_phi_vec_eps = ComputeFphi();

                    // do the same on the left
                    LevelSetOpti.SetParam(n_param, dx_left);
                    LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                    f_phi_vec_eps2 = ComputeFphi();

                    //Compute FD and write value into the matrix
                    for(int iRow = 0; iRow < f_phi_vec.Length; iRow++) {
                        val = (f_phi_vec_eps[iRow] - f_phi_vec_eps2[iRow]) / Epsilons[n_param];
                        if(val != 0) {
                            Jf_phi_vec[iRow, n_param] = val;
                        }

                    }

                    //reset
                    LevelSetOpti.SetParam(n_param, x);
                    LsTBO.CopyFrom(phi0backup);
                }
                return (f_phi_vec, Jf_phi_vec);
            }
        }

        Func<double[]> ComputeFphi;

        public void ChooseFphiType() {
            switch(Control.fphiType) {
                case FphiType.None:
                    ComputeFphi= delegate() { 
                        return new double[1]; 
                    };
                break;
                case FphiType.CurvatureAll:
                ComputeFphi = delegate () {
                    CellMask mask = CellMask.GetFullMask(this.LsTrk.GridDat);
                    return CurvFphi(mask);
                };
                break;
                case FphiType.CurvatureCut:
                ComputeFphi = delegate () {
                    CellMask mask = LsTrk.Regions.GetCutCellMask();
                    return CurvFphi(mask);
                };
                break;
                case FphiType.PerssonSensorAll:
                ComputeFphi = delegate () {
                    CellMask mask = CellMask.GetFullMask(this.LsTrk.GridDat);
                    return PerssonSensorFphi(mask);
                };
                break;
                case FphiType.PerssonSensorCut:
                ComputeFphi = delegate () {
                    CellMask mask = LsTrk.Regions.GetCutCellMask();
                    return PerssonSensorFphi(mask);
                };
                break;
            }
        }
        public double[] PerssonSensorFphi(CellMask mask) {
            double[] fphi = new double[mask.ItemEnum.Count()];
            var deg = LsTBO.Basis.Degree;
            if(LsTBO.Basis.Degree > 0) {
                var ba = new BitArray(Grid.NumberOfCells); //empty array
                int count = 0;
                foreach(int iCell in mask.ItemEnum) { // loop over all cells we consider
                    ba[iCell] = true; //set corresponding entry
                    var tmpmsk = new CellMask(GridData, ba);
                    // Copy only the coordinates that belong to the highest modes 
                    foreach(int coordinate in LsTBO.Basis.GetPolynomialIndicesForDegree(iCell, deg)) {
                        fphi[count] = fphi[count]+Math.Pow(LsTBO.Coordinates[iCell, coordinate],2);
                    }
                    fphi[count] = fphi[count] / Math.Pow(LsTBO.L2Norm(tmpmsk),2);
                    ba[iCell] = false; // reset
                    count++;
                }
            }
            return fphi;
        }
        public double[] GradFphi(CellMask mask) {
            double[] fphi = new double[mask.ItemEnum.Count()];
            if(LsTBO.Basis.Degree > 0) {
                Basis basisForGrad = new Basis(Grid, LsTBO.Basis.Degree - 1);
                SinglePhaseField grad;
                LsTBO.GradientNorm(out grad, mask);

                var ba = new BitArray(Grid.NumberOfCells); //empty array
                int count = 0;
                foreach(int iCell in mask.ItemEnum) { // loop over all cells we consider
                    ba[iCell] = true; //set corresponding entry
                    var tmpmsk = new CellMask(GridData, ba);
                    fphi[count] = grad.L2Norm(tmpmsk);
                    ba[iCell] = false; // reset
                    count++;
                }
            }
            fphi.ScaleV(kappa);
            return fphi;
        }

        public double[] CurvFphi(CellMask mask) {
            double[] fphi = new double[mask.ItemEnum.Count()];
            if(LsTBO.Basis.Degree > 1) {
                Basis basisForCurvature = new Basis(Grid, LsTBO.Basis.Degree - 2);
                SinglePhaseField totCurv = new SinglePhaseField(basisForCurvature);


                LsTBO.ComputeTotalCurvature(totCurv, mask);
                totCurv.Scale(0.5); // we want the mean curvature, we don't care about the sign

                var ba = new BitArray(Grid.NumberOfCells); //empty array
                int count = 0;
                foreach(int iCell in mask.ItemEnum) { // loop over all cells we consider
                    ba[iCell] = true; //set corresponding entry
                    var tmpmsk = new CellMask(GridData, ba);
                    fphi[count] = totCurv.L2Norm(tmpmsk);
                    ba[iCell] = false; // reset
                    count++;
                }
            }
            fphi.ScaleV(kappa);
            return fphi;
        }

        /// <summary>
        /// This method creates a backup of all the state variables
        /// </summary>
        public void createBackUp() {
            levelSetBackup = LsTBO.CloneAs();
            UBackup = new XDGField[ConservativeFields.Length];
            for(int i = 0; i < ConservativeFields.Length; i++) {
                UBackup[i] = new XDGField(ConservativeFields[i].Basis, ConservativeFields[i].Identification + "_copy");
                UBackup[i].CoordinateVector.SetV(ConservativeFields[i].CoordinateVector);
            }
            LevelSetOptiBackup = LevelSetOpti.CloneAs();
        }

        /// <summary>
        /// This method creates a backup of all the agglomerated state variables and the levelset
        /// </summary>
        public void createAgglomeratedBackUp()
        {
            levelSetBackup = LsTBO.CloneAs();
            AgglomeratedUBackup = new XDGField[ConservativeFields.Length];
            for (int i = 0; i < ConservativeFields.Length; i++)
            {
                AgglomeratedUBackup[i] = new XDGField(ConservativeFields[i].Basis, ConservativeFields[i].Identification + "_copy");
                AgglomeratedUBackup[i].CoordinateVector.SetV(ConservativeFields[i].CoordinateVector);
            }
            LevelSetOptiBackup = LevelSetOpti.CloneAs();
        }
        /// <summary>
        /// this method resets the actual state to the backup
        /// </summary>
        public void resetStep() {
            if(levelSetBackup == null || UBackup == null) {
                throw new NotSupportedException("cannot reset with not backup available");
            }
            LsTBO.CopyFrom(levelSetBackup);
            LsTrk.UpdateTracker(CurrentStepNo);  
            for(int i = 0; i < ConservativeFields.Length; i++) {
                ConservativeFields[i].CopyFrom(UBackup[i]);
            }
            LevelSetOpti.CopyParamsFrom(LevelSetOptiBackup);
        }
        /// <summary>
        /// this method resets the actual state to the agglomerated backup and the levelset
        /// </summary>
        public void resetStepAgglomerated()
        {
            if (levelSetBackup == null || AgglomeratedUBackup == null)
            {
                throw new NotSupportedException("cannot reset with not backup available");
            }
            LsTBO.CopyFrom(levelSetBackup);
            LsTrk.UpdateTracker(CurrentStepNo);
            for (int i = 0; i < ConservativeFields.Length; i++)
            {
                ConservativeFields[i].CopyFrom(AgglomeratedUBackup[i]);
            }
            LevelSetOpti.CopyParamsFrom(LevelSetOptiBackup);
        }
        /// <summary>
        /// this method computes a custom norm (basically the L2-Norm of the  coordinates ) where one can choose which entries in the double[] x are committed. 
        /// This is used here as the step in our method also contains DOFs corresponding to the Lagrange multiplier lambda which shall be omitted.
        /// </summary>
        /// <param name="x"> array containing the unknowns </param>
        /// <param name="length_r">first length should equal state DOFs</param>
        /// <param name="length_l">second length should equal LevelSet DOFs</param>
        /// <returns></returns>
        double StepNorm(double[] x, int length_r, int length_l) {
            double normU = 0;
            double normL = 0;
            if(x.Length < length_r + length_l)
                throw new ArgumentException();
            for(int i = 0; i < length_r; i++) {
                normU += x[i] * x[i];
            }
            for(int i = length_r; i < length_l; i++) {
                normL += x[i] * x[i];
            }
            return Math.Sqrt(normU + normL);
        }
        /// <summary>
        /// same as StepNorm but as a InnerProduct
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="length_r"></param>
        /// <param name="length_l"></param>
        /// <returns></returns>
        double StepInnerProd(double[] x, double[] y, int length_r, int length_l) {
            double normU = 0;
            double normL = 0;
            if(x.Length < length_r + length_l)
                throw new ArgumentException();
            for(int i = 0; i < length_r; i++) {
                normU += x[i] * y[i];
            }
            for(int i = length_r; i < length_l; i++) {
                normL += x[i] * y[i];
            }
            return normU + normL;
        }
        /// <summary>
        /// Transforms the state solution to another degree 
        /// </summary>
        /// <param name="degree"></param>
        public void ChangeTOBasisofDegree(int degree) {
            if(degree == ConservativeFields[0].Basis.Degree) {
                Console.WriteLine("ChangeTOBasisofDegree(): No need to do anything as degrees are the same");
            } else {

                //create a backup of the current field
                createBackUp();

                int oldDeg = ConservativeFields[0].Basis.Degree;
                this.m_IOFields.Clear();
                this.m_RegisteredFields.Clear();

                //Initiate fields of order degree
                CreateConservativeFields(LsTrk, degree);

                //initialize all dependent fields
                InitDependentFields();
                InitObjFields();

                //register the fields
                RegisterFields();

                //Recreate Matrices
                InitializeMatricesAndVectors();

                //Init MGridOpConfig
                InitializeMultiGridOpConfig();

                if(degree > oldDeg) {
                    //Copy DOFS from backup
                    for(int iField = 0; iField < ConservativeFields.Length; iField++) {
                        ConservativeFields[iField].AccLaidBack(1.0, UBackup[iField]);
                    }
                } else if(degree < oldDeg) {
                    //Do a projection
                    for(int iField = 0; iField < ConservativeFields.Length; iField++) {
                        var currentDeg = UBackup[iField].Basis.Degree - 1;
                        var proj = GetPMinus1Projection(LsTrk, UBackup[iField], UBackup[iField].Basis.Degree * currentDeg);
                        while(currentDeg > degree) {
                            currentDeg--;
                            proj = GetPMinus1Projection(LsTrk, proj, proj.Basis.Degree * currentDeg);
                        }
                        ConservativeFields[iField].Acc(1.0, proj);
                    }
                }
                gamma = Control.Gamma_Start;
                CurrentAgglo = Control.AgglomerationThreshold;
            }
            
            
        }
        public int GetCurrentTermN() {
            switch(Control.solRunType) {
                case SolverRunType.Standard:
                    return Control.TerminationMinNs[0];
                case SolverRunType.PContinuation:
                try {
                    return Control.TerminationMinNs[ConservativeFields[0].Basis.Degree];
                } catch {
                    return Control.TerminationMinNs.Last();
                }
                default:
                    return 5;
            }
        }
        public double GetCurrenttALNR() {
            switch(Control.solRunType) {
                case SolverRunType.Standard:
                return Control.tALNRs[0];
                case SolverRunType.PContinuation:
                try {
                    return Control.tALNRs[ConservativeFields[0].Basis.Degree];
                } catch {
                    return Control.tALNRs.Last();
                }
                default:
                    return 1.001;
            }
        }
        public void CheckTermination() {
            bool success = false;
            double tALNR = GetCurrenttALNR();
            int TerminationMinN = GetCurrentTermN();
            
            Console.WriteLine("");            
            if(CheckTermination(ref success, res_l2, InitResNorm, ResNorms, CurrentStepNo, TerminationMinN, tALNR) &&
                CheckTermination(ref success, obj_f, Init_obj_f, obj_f_vals, CurrentStepNo, TerminationMinN, tALNR)) {
                if(ConservativeFields[0].Basis.Degree >= Control.SolDegree) {
                    this.TerminationKey = true;
                } else {
                    IncreaseDegreeNextTS = true;
                }
            }
        }

        bool CheckTermination(ref bool success, double norm_CurRes, double fnorminit, List<double> normHistory, int itc, int N, double tALNR) {
            using(var tr = new FuncTrace()) {


                tr.InfoToConsole = false;
                tr.Info($"Checking termination criterion: Iter {itc}, Current residual norm {norm_CurRes}, Initial Residual Norm {fnorminit}");
                var MinIter = CurMinIter;

                if(Control.ConvCrit > 0) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // **if** the convergence criterion is set by the user, we should aim for it
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    tr.Info("NEWTON: going for convergence criterion " + Control.ConvCrit);


                    if(itc < MinIter) {
                        tr.Info("below minimal number of iterations");
                        return false;
                    }

                    if(norm_CurRes <= Control.ConvCrit * fnorminit + Control.ConvCrit) {
                        success = true;
                        tr.Info("Newton iterations SUCCESS!");
                        return true;
                    }

                    if(itc >= Control.MaxIterations) {
                        // run out of iterations
                        return true;
                    }

                    return false;
                } else {
                    bool terminateLoop = false;


                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // ... otherwise, we do the best we can and iterate until no further improvement can be made
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    tr.Info($"NEWTON: trying to converge as far as possible (ConvCrit is {Control.ConvCrit})");

                    //
                    // See:
                    // Kikker, Kummer, Oberlack:
                    // A fully coupled high‐order discontinuous {Galerkin} solver for viscoelastic fluid flow,
                    // IJNMF, No. 6, Vol. 93, 2021, 1736--1758,
                    // section 3.3.1
                    //

                    const double _MinConvCrit = 1e-6;

                    double LastAverageNormReduction() {
                        //const int N = 5; // look at latest N residuals
                        if(N < 1)
                            throw new ArgumentException();

                        if(normHistory.Count - N < 0)
                            return 1e100; // ignore if we have not at least 'N'  residuals so far.

                        int i0 = normHistory.Count - N;

                        // take the (minimum) skyline to be immune against residual oscillations...
                        double[] NormHistorySkyline = new double[normHistory.Count];
                        NormHistorySkyline[0] = normHistory[0];
                        for(int i = 1; i < normHistory.Count; i++) {
                            NormHistorySkyline[i] = Math.Min(NormHistorySkyline[i - 1], normHistory[i]);
                        }


                        double Avg = 0;
                        double Count = 0;
                        for(int i = i0; i < normHistory.Count - 1; i++) { // look at the last 'N' residual norms...
                            double ResNormReductionFactor = NormHistorySkyline[i] / Math.Max(NormHistorySkyline[i + 1], double.Epsilon);
                            if(ResNormReductionFactor < 1)
                                ResNormReductionFactor = 1; // should never happen anyway due to sky lining...
                            Count = Count + 1;
                            Avg = Avg + ResNormReductionFactor;
                        }
                        return Avg / Count;
                    }
                    double DistMinMax() {
                        if(N < 1)
                            throw new ArgumentException();
                        if(normHistory.Count - N < 0)
                            return 1e100; // ignore if we have not at least 'N'  residuals so far.
                        int i0 = normHistory.Count - N;

                        double[] lastNhistory = new double[N];

                        for(int i = i0; i < normHistory.Count; i++) {
                            lastNhistory[i - i0] = normHistory[i];
                        }
                        var min = lastNhistory.Min();
                        var max = lastNhistory.Max();

                        return Math.Abs(min - max) / Math.Abs(lastNhistory.Average());

                    }

                    if(itc >= MinIter) {

                        
                            
                            // only terminate if we reached the minimum number of iterations
                            //if(norm_CurRes <= _MinConvCrit * fnorminit + _MinConvCrit) {
                            // reached minimum convergence criterion

                            double ALNR = LastAverageNormReduction();
                            tr.Info($"LastAverageNormReduction is {ALNR}, {norm_CurRes} ...");
#if TEST
                        ALNR = 0;
#endif
                            // continue until the solver stalls
                            if(ALNR <= tALNR) {
                                Console.WriteLine();
                                Console.WriteLine($"ALNR={ALNR} <= {tALNR} so we terminate/increase degree");
                                tr.Info("... no sufficient further progress - terminating.");
                                success = true;
                                terminateLoop = true;
                            } else {
                                tr.Info("... but continue as long as we make progress;");
                                Console.WriteLine($"ALNR={ALNR} > {tALNR} so we continue");
                                //success = true;
                                //terminateLoop = true;
                            }
                            //} else {
                            //    tr.Info($"minimal build-in convergence criterion NOT reached (current residual is {norm_CurRes}, limit is {_MinConvCrit * fnorminit + _MinConvCrit}) yet.");
                            //}

                            if(itc >= Control.MaxIterations) {
                                // run out of iterations
                                tr.Info($"Maximum number of iterations reached ({Control.MaxIterations}) - terminating.");
                                terminateLoop = true;
                                TerminationKey = false;
                            }
                    } else {
                        Console.WriteLine($"**NOT** terminating now, minimal number of iterations not reached");
                        tr.Info($"**NOT** terminating now, minimal number of iterations (iteration {itc}, minimum is {MinIter}) not reached yet; CurRes = {norm_CurRes}, threshold is {_MinConvCrit * fnorminit + _MinConvCrit}, initial norm was {fnorminit}");
                    }

                    return terminateLoop;

                }

            }
        }



        /// <summary>
        /// The first step computed can be very bad and cause LevelSetCFL or the operator to crash.
        /// This method iteratively downscales the step ($s_k,\tau * s_k,\tau^2 * s_k,$) until it doesn't cause any of these problems
        /// </summary>
        /// <param name="step_t"> the step </param>
        /// <param name="tau"> the scaling factor </param>
        public void AccumulateInitialStep(double[] step_t, double tau) {
            if(tau >= 1 || tau <= 0) {
                throw new ArgumentException(" the scaling factor tau=" + tau + " must be between 0 and 1");
            }
            bool success = false;
            int counter1 = 0;
            int counter2 = 0;
            m_alpha = Control.Alpha_Start;
#if DEBUG
            //Create second Backup Debugging
            //double[] rho_backup = new double[UBackup[0].CoordinateVector.Length];
            //rho_backup.SetV(UBackup[0].CoordinateVector);
#endif

#region LevelSetCFL
            // we will check the LevelSetCFL not with the original LevelSetTrecker but with a nonShallow Copy
            // This is dirty but causing the LevelSetCFL with the original LevelSetTracker breaks everything 
            LevelSet LevelSet_copy;
            LevelSetTracker LsTrk_copy;

            int iField = 0;
            int jCell = 0;
            int nMode = 0;
            int length_r = (int)UnknownsMap.TotalLength; // needs to be modified if more than one field is simulated

            while(success == false && m_alpha >= Control.Alpha_Min) {
                success = true;

                //Get a copy from the LevelSet and the Tracker
                if(Control.IsTwoLevelSetRun) {
                    LevelSet_copy = LevelSetTwo.CloneAs();
                    LsTrk_copy = new LevelSetTracker((GridData)this.GridData, Control.CutCellQuadratureType, Control.NearRegionWidth, Control.SpeciesTable, this.LevelSet, LevelSet_copy);
                } else {
                    LevelSet_copy = LevelSet.CloneAs();
                    LsTrk_copy = new LevelSetTracker((GridData)this.GridData, Control.CutCellQuadratureType, Control.NearRegionWidth, new string[] { Control.LsOne_NegSpecies, Control.LsOne_PosSpecies }, LevelSet_copy);
                }

                // add the step to the state DOFs - scaled by m_alpha
                for(int stepIndex = 0; stepIndex < length_r; stepIndex++) {
                    UnknownsMap.LocalFieldCoordinateIndex(stepIndex, out iField, out jCell, out nMode);
                    ConservativeFields[iField].CoordinateVector[jCell * ConservativeFields[iField].Basis.DOFperSpeciesPerCell * LsTrk.TotalNoOfSpecies + nMode] += m_alpha * step_t[stepIndex];
                }

                // add the step to the LevelSet DOFs - scaled by m_alpha
                for(int i = 0; i < LevelSetOpti.GetLength(); i++) {
                    LevelSetOpti.AccToParam(i, m_alpha * step_t[i + length_r]);
                }
                //project onto LevelSet object 
                LevelSetOpti.ProjectOntoLevelSet(LevelSet_copy);

                //var tp = new Tecplot(GridData, 4);
                //tp = new Tecplot(GridData, 4);

                //List<DGField> list = new List<DGField>();
                //list.AddRange(UBackup);
                //list.Add(LevelSet_copy);

                //tp.PlotFields("BackUp" + "_" + 0, 0.0, list);

                try {
                    LsTrk_copy.UpdateTracker(CurrentStepNo, Control.NearRegionWidth);
                } catch { //if Exception is called no step is computed
                    Console.WriteLine($"alpha, ||R1|| {m_alpha} LevelSetCFL ");
                    success = false; //continue the loop
                    m_alpha *= tau;  //reduce m_alpha if Exception is thrown
                    counter1++;
                    //tp.PlotFields("BackUp" + "_" + 1, 0.0, list);
                }
                resetStep();
            }

            //if(counter1 > 0) {
            //    //Console.WriteLine("step induces LevelSetCFL and had to be shortened by" + Math.Pow(tau, counter1));
            //    Console.WriteLine("step induces LevelSetCFL and had to be shortened to alpha =" + m_alpha);
            //}

#if DEBUG
            //if(!rho_backup.SequenceEqual(UBackup[0].CoordinateVector)) {
            //    throw new Exception("WTF");
            //}
#endif
#endregion
            //so now we should have a step that doesn't throw any LevelSetCFLException so we can compute it
            success = !AccumulateStep(step_t, m_alpha);
#region Operator Exceptions
            //next we want to find a step which doesn't break our operator
            //So this loop gives us a Step which doesn't break the operator (e.g. negative rho/pressure etc...)
            
            while(success == false) {
                success = true;
                try {


                    TransformFromAggToSourceSpace();

                    double[] obj_f = new double[obj_f_vec.Length];
                    double[] res = new double[ResidualVector.Length];
                    //Eval_R = XSpatialOperator.GetEvaluatorEx(LsTrk, ConservativeFields, null, obj_f_map);
                    //Eval_R.Evaluate(1.0, 0.0, en_res);
                    ComputeResiduals(obj_f, res);
                    if(obj_f.MPI_L2Norm().IsNaN()) {
                        throw new NotFiniteNumberException("enriched residual norm is NaN");
                    }
                    //double[] res = new double[ResidualVector.Length];
                    //Eval_r = XSpatialOperator.GetEvaluatorEx(LsTrk, ConservativeFields, null, ResidualMap);
                    //Eval_r.Evaluate(1.0, 0.0, res);

                } catch (Exception e) {
                    if(e is NotFiniteNumberException) {
                        Console.WriteLine($"alpha, ||R1||  {m_alpha}, NaN ");
                    } else {
                        Console.WriteLine($"alpha, ||R1||  {m_alpha}, Error ");
                    }
                    resetStep();
                    counter2++;
                    m_alpha *= tau; //further reduce m_alpha if Exception is thrown
                    success = AccumulateStep(step_t, m_alpha);
                    success = false; //continue the loop
                    if(m_alpha < Control.Alpha_Min) {
                        success = true;
                    }
                }
                
            }
            // if alpha hat to be shortened to much throw exception
            if(m_alpha < Control.Alpha_Min) {
                AllthePossibleStepsPlot(step_t, tau);
                throw new Exception("step needed to be shortened to much (m_alpha < alpha_min)");
            }

            //if(counter2 > 0) {
            //    //Console.WriteLine("step breaks Operator and had to be shortened to alpha=" + Math.Pow(tau, counter1 + counter2));
            //    Console.WriteLine("step breaks Operator and had to be shortened to alpha=" + m_alpha);
            //}
#endregion
            resetStep();
#if DEBUG
            //if(!rho_backup.SequenceEqual(UBackup[0].CoordinateVector)) {
            //    throw new Exception("WTF");
            //}
#endif
            ////Finally the variables are reseted and the step is scaled with the m_alpha obtained 
            //if(m_alpha < 1) {
            //    step_t.ScaleV(m_alpha);
            //}

        }
        public void AllthePossibleStepsPlot(double eps=1e-8)
        {
            var length_r = (int)ResidualMap.TotalLength;
            //Plot the current state
            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            Plot(0.0, 1000);
            //Plot the possible steps
            int i = 0;

            int nCol_obj = Oproblem.GetObjLength(ConservativeFields);
            int nCol_con = new CoordinateVector(ConservativeFields).Count;
            int nRow = LevelSetOpti.GetLength();

            MsrMatrix Jobj = new MsrMatrix(nCol_obj, nRow, 1, 1);
            MsrMatrix Jr = new MsrMatrix(nCol_con, nRow, 1, 1);

            var r_vec = new double[nCol_con];
            var r_vec_eps = new double[nCol_con];
            var r_vec_eps2 = new double[nCol_con];

            var obj_vec = new double[nCol_obj];
            var obj_vec_eps = new double[nCol_obj];
            var obj_vec_eps2 = new double[nCol_obj];


            double x;
            double val;
            double dx_right;
            double dx_left;

            // epsilon
            //project OptiLevelSet onto XDGLevelSet
            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            LsTrk.UpdateTracker(CurrentStepNo);
            LevelSet phi0backup = new LevelSet(new Basis(LsTBO.GridDat.Grid, LsTBO.Basis.Degree), "LevelSetbackup");
            phi0backup.CopyFrom(LsTBO);

            //Evaluation of unperturbed
            Oproblem.EvalConsAndObj(obj_vec, r_vec, ConservativeFields);


            for (int n_param = 0; n_param < nRow; n_param++)
            {

                //skip loop if param is non changeable
                if (Control.PartiallyFixLevelSetForSpaceTime && LevelSetOpti is SplineOptiLevelSet splineLS)
                {
                    double yMin = splineLS.y.Min(); //lower boundary of space time domain
                    if (n_param < splineLS.y.Length && Math.Abs(yMin - splineLS.y[n_param]) < 1e-14) //only accumalte if DOF if it is not on lower time boundary (here yMin=tMin)
                    {
                        continue;
                    }

                }

                //compute distortions
                x = LevelSetOpti.GetParam(n_param);
                dx_right = x + eps / 2;
                dx_left = x - eps / 2;

                // apply right distortion
                LevelSetOpti.SetParam(n_param, dx_right);
                //project
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                
                //compute Objective

                try
                {
                    LsTrk.UpdateTracker(CurrentStepNo);

                    (double _resl2, double _obj, double _resL2) = ComputeResiduals();
                    resetStep();
                }
                catch (Exception e)
                {
                    Console.WriteLine("Exception at right fd of n_param = " +n_param + ": ");
                    Console.WriteLine(e.ToString());
                    Console.WriteLine("");
                    resetStep();
                }
                Plot(0.0, 1000 + (n_param + 1));

                // do the same on the left
                LevelSetOpti.SetParam(n_param, dx_left);
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                LsTrk.UpdateTracker(CurrentStepNo);

                try
                {
                    LsTrk.UpdateTracker(CurrentStepNo);
                    (double _resl2, double _obj, double _resL2) = ComputeResiduals();
                    resetStep();
                }
                catch (Exception e)
                {
                    Console.WriteLine("Exception at left fd of n_param = " + n_param + ": ");
                    Console.WriteLine(e.ToString());
                    Console.WriteLine("");
                    resetStep();
                }
                Plot(0.0, 1000 - (n_param + 1));

                //reset
                LevelSetOpti.SetParam(n_param, x);
                LsTBO.CopyFrom(phi0backup);
                LsTrk.UpdateTracker(CurrentStepNo);
            }
        }
        public void AllthePossibleStepsPlot<Tin>(Tin step_t,double tau ) where Tin: IList<double> {
            var length_r = (int) ResidualMap.TotalLength; 
            //Plot the current state
            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            Plot(0.0, 1000);
            //Plot the possible steps
            int i = 0;
            var alpha_p = Control.Alpha_Start;
            while(alpha_p >= Control.Alpha_Min) {
                alpha_p = Control.Alpha_Start * Math.Pow(tau, i);
                for(int stepIndex = 0; stepIndex < length_r; stepIndex++) {
                    UnknownsMap.LocalFieldCoordinateIndex(stepIndex, out int iField, out int jCell, out int nMode);
                    ConservativeFields[iField].CoordinateVector[jCell * ConservativeFields[iField].Basis.DOFperSpeciesPerCell * LsTrk.TotalNoOfSpecies + nMode] += alpha_p * step_t[stepIndex];
                }
                for(int j = 0; j < LevelSetOpti.GetLength(); j++) {
                    LevelSetOpti.AccToParam(j, alpha_p * step_t[j + length_r]);
                }
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                i++;
                try {
                    LsTrk.UpdateTracker(CurrentStepNo);

                    TransformFromAggToSourceSpace();


                    (double _resl2, double _obj, double _resL2) = ComputeResiduals();

                    Plot(0.0, 1000 + i + 1);
                    resetStep();
                } catch (Exception e){
                    Console.WriteLine("Exception at alpha = "+ alpha_p +": ");
                    Console.WriteLine(e.ToString());
                    Console.WriteLine("");
                    Plot(0.0, 1000 + i + 1);
                    resetStep();
                    m_alpha *= tau; //further reduce m_alpha if Exception is thrown
                }
            }
        }
        /// <summary>
        /// accumulates a step onto the current state and updates the LsTrk
        /// </summary>
        /// <param name="step_t">step to be accumulated</param>
        /// <param name="m_alpha">step size </param>
        /// <returns>if successful (no LevelSetCFL thrown) </returns>
        public bool AccumulateStep(double[] step_t, double m_alpha) {

            int iField = 0;
            int jCell = 0;
            int nMode = 0;
            int length_r = (int)UnknownsMap.TotalLength; // needs to be modified if more than one field is simulated

            for(int stepIndex = 0; stepIndex < length_r; stepIndex++) {
                UnknownsMap.LocalFieldCoordinateIndex(stepIndex, out iField, out jCell, out nMode);
                ConservativeFields[iField].CoordinateVector[jCell * ConservativeFields[iField].Basis.DOFperSpeciesPerCell * LsTrk.TotalNoOfSpecies + nMode] += m_alpha * step_t[stepIndex];
            }


            /// Level Set Update
            //special treatment for space time level sets
            if(Control.PartiallyFixLevelSetForSpaceTime && LevelSetOpti is SplineOptiLevelSet splineLS)
            {
                double yMin = splineLS.y.Min(); //lower boundary of space time domain
                for (int i = 0; i < LevelSetOpti.GetLength(); i++)
                {
                    if (i<splineLS.y.Length && Math.Abs(yMin - splineLS.y[i])> 1e-14) //only accumalte if DOF if it is not on lower time boundary (here yMin=tMin)
                    {
                        LevelSetOpti.AccToParam(i, m_alpha * step_t[i + length_r]); 
                    }
                }
            }
            else
            {
                for (int i = 0; i < LevelSetOpti.GetLength(); i++)
                {
                    LevelSetOpti.AccToParam(i, m_alpha * step_t[i + length_r]);
                }
            }

            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
            //var tp = new Tecplot(GridData, 4);
            //tp = new Tecplot(GridData, 4);

            //List<DGField> list = new List<DGField>();
            //list.AddRange(UBackup);
            //list.Add(LsTBO);
            //tp.PlotFields("BackUp" + "_" + 0, 0.0, list);

            try {
                LsTrk.UpdateTracker(CurrentStepNo);
                //UpdateAgglomerator();   
            } catch { //if Exception is called no step is computed
                resetStep();
                //tp.PlotFields("BackUp" + "_" + 1, 0.0,list);

                return false;
            }


            return true;

        }
        /// <summary>
        /// Transform the solutionVector back into the original Space (Non-agglomerated)
        /// </summary>
        public void TransformFromAggToSourceSpace() {

            //CoordinateVector SolutionVec = new CoordinateVector(ConservativeFields);
            //double[] v = new double[SolutionVec.Count];
            ////CoordinateVector SolutionVec_agg = new CoordinateVector(ConservativeFields);
            ////SolutionVec_agg.Clear();
            //LeftMul.Transpose().SpMV(1.0, SolutionVec, 0.0, v);
            //SolutionVec.SetV(v);
            //UpdateEnrichedFields();
            if(CurrentAgglo > 0) {
                if(MultiphaseAgglomerator.TotalNumberOfAgglomerations > 0) {
                    MultiphaseAgglomerator.Extrapolate(new CoordinateMapping(ConservativeFields));//reset ConservativeFields into non-agglomerated form
                }
            }
                   

            
        }
        public void TransformStepFromAggToSourceSpace(double[] vector)
        {

            CoordinateVector SolutionVec = new CoordinateVector(ConservativeFields);
            var dummy = ConservativeFields.CloneNonshallow();
            CoordinateVector vec= new CoordinateVector(dummy);
            vec.SetV(vector.GetSubVector(0, vec.Length));
            if (CurrentAgglo > 0)
            {
                if (MultiphaseAgglomerator.TotalNumberOfAgglomerations > 0)
                {
                    MultiphaseAgglomerator.Extrapolate(vec.Mapping);//reset ConservativeFields into non-agglomerated form
                }
            }
            vector.SetSubVector(vec.ToArray(), 0, vec.Length);
        }
        /// <summary>
        /// Transform the solutionVector into AggSpace (from the non-agglomerated into the agglomerated)
        /// </summary>
        public void TransformFromSourceToAggSpace() {
            if(CurrentAgglo > 0) {
                if(MultiphaseAgglomerator.TotalNumberOfAgglomerations > 0) {
                    CoordinateVector SolutionVec = new CoordinateVector(ConservativeFields);
                    var SolutionVec_agglo = new double[multOp.Mapping.TotalLength];
                    multOp.TransformSolInto(SolutionVec, SolutionVec_agglo); //this gives a shorter vector
                    multOp.TransformSolFrom(SolutionVec, SolutionVec_agglo); //this makes the vector have the right length

                }
            }

        }

        /// <summary>
        /// This method serves as a globalization to this solver.Starting with the computed step s^{IN} and depending on the globalization chosen a new step s
        /// is computed which either satisfies the condition of sufficient decrease
        /// $$ f_m(z_k +s) \leq f_m(z_k) + s^T f_M'(z_k) $$ 
        /// or was shortened by $alpha_min$.
        /// 
        /// The globalization strategies implemented so far:
        /// 
        /// 1. None : $s=s^{IN}$
        /// 2. Line Search: try $(s^{IN}, \tau s^{IN}, \tau^2 s^{IN}, ....)$
        /// 3. DogLeg: compute another search direction (the DogLeg-step) $s^{DL}$ and try some convex combinations $s_\lambda = \lambda s^{IN} + (1-\lambda) s^{DL)$ such that $\Vert s \Vert = \delta$
        /// 4. ZickZack LineSearch: apply different alphas for the state Coordinates and the LevelSet Coordinates
        /// 
        /// </summary>
        /// <returns></returns>
        public double SQPStep() {
            using (new FuncTrace())
            {
                double tau = Control.tauAlpha;
                double beta = 1e-4;
                double alpha_min = Control.Alpha_Min;
                bool succes = false;


                //Evaluate the Operator at the Current state (normally this evaluation is done in RunSolverOneStep)
                (res_l2, obj_f, res_L2) = ComputeResiduals();
                //resetStep();
                //TransformFromAggToSourceSpace();
                //(res_l2, obj_f, res_L2) = ComputeResiduals(); 
                //double merit_0=res_l2 + obj_f;
                //Goes Back to NonAggSpace
                //resetStepAgglomerated();

                /// pre calculate things later used to calculate the predicted merit function
                //switch (Control.MeritFunctionType)
                //{
                //    case MeritFunctionType.L1Merit:
                //        res_l1 = 0;
                //        for (int i = 0; i < ResidualVector.Length; i++)
                //        {
                //            res_l1 += ResidualVector[i].Abs();
                //        }
                //        break;
                //    case MeritFunctionType.L2Merit:
                //        Jr.Transpose().SpMV(mu * beta, ResidualVector, 0, merit_del);
                //        merit_del.AccV(beta, gradf);
                //        break;
                //    default: throw new Exception("you need to choose a MeritFunctionType");

                //}
                //transform everything back to non-agglomerated
                int length_R = (int)obj_f_map.TotalLength;
                int length_r = (int)UnknownsMap.TotalLength; // needs to be modified if more than one field is simulated
                int length_l = LevelSetOpti.GetLength();

                stepUphi.SetSubVector(stepIN, 0, (int)UnknownsMap.TotalLength + LevelSetOpti.GetLength());

                JacR0TimesStep = new double[length_r];
                Jr.SpMV(1.0, stepUphi, 0.0, JacR0TimesStep);

                TransformStepFromAggToSourceSpace(gradf);
                TransformStepFromAggToSourceSpace(JacR0TimesR0);
                TransformStepFromAggToSourceSpace(JacR0TimesStep);
                TransformStepFromAggToSourceSpace(stepIN);


                ///Safe Guard: see if step breaks operator or induces LevelSet CFL, if so, shorten the step 
                AccumulateInitialStep(stepIN, tau);

                double[] temp = new double[RHS.Length];
                double[] step = new double[stepIN.Length];




                switch (Control.GlobalizationStrategy)
                {
                    #region DogLeg
                    case GlobalizationStrategy.DogLeg: /// WARNING: May be deprecated
                        // compute Cauchy point
                        // ====================
                        Console.WriteLine("WARNING: DogLeg globalization was not used for a long time, and migth contain errors.");
                        double[] stepCP;
                        {
                            // step 1: calculate direction of Cauchy point / direction of steepest decent
                            if (RHS.Length != (int)LHS.NoOfRows)
                            {
                                throw new ArgumentException("No of RHS entries must equal LHS Rows, but LHS:" + (int)LHS.NoOfRows + ", RHS:" + RHS.Length);
                            }
                            double[] dk = new double[stepIN.Length];
                            MsrMatrix HTransp = LHS.Transpose();
                            HTransp.SpMV(1.0, RHS, 0.0, dk);
                            // step 2: computing Cauchy point
                            double[] w = new double[RHS.Length];
                            LHS.SpMV(1.0, dk, 0.0, w);

                            double lambda = RHS.InnerProd(w) / w.InnerProd(w);
                            stepCP = dk;
                            stepCP.ScaleV(lambda);
                            dk = null; // invalidate
                        }


                        // find point on Dogleg curve, within the trust region
                        // ===================================================

                        double l2_stepCP = StepNorm(stepCP, length_r, length_l);
                        double l2_stepIN = StepNorm(stepIN, length_r, length_l);
                        void PointOnDogleg(double _TrustRegionDelta, double[] step)
                        {
                            if (l2_stepIN <= _TrustRegionDelta)
                            {
                                // use Newton Step
                                step.SetV(stepIN);
                            }
                            else
                            {
                                if (l2_stepCP < _TrustRegionDelta)
                                {
                                    // interpolate between Cauchy-point and Newton-step
                                    double tau;
                                    {
                                        double A = l2_stepCP, B = l2_stepIN, C = C = StepInnerProd(stepCP, stepIN, length_r, length_l);
                                        tau = (A.Pow2() - C + Math.Sqrt((A.Pow2() + B.Pow2() - 2 * C) * _TrustRegionDelta.Pow2() - A.Pow2() * B.Pow2() + C.Pow2()))
                                            / (A.Pow2() + B.Pow2() - 2 * C);
                                        if (!(tau >= -0.00001 && tau <= 1.00001) || tau.IsNaN() || tau.IsInfinity())
                                            throw new ArithmeticException();
                                    }
                                    // do interpolation
                                    step.SetV(stepCP, (1 - tau));
                                    step.AccV(tau, stepIN);
                                    //Debug.Assert(Math.Abs((step.MPI_L2Norm() / _TrustRegionDelta) - 1.0) <= 1.0e-3, "interpolation step went wrong");
                                }
                                else
                                {
                                    // use reduced Cauchy point
                                    step.SetV(stepCP, (_TrustRegionDelta / l2_stepCP));
                                }
                            }
                        }

                        const double delta_min = 1e-6;
                        const double delta_max = 1e10;
                        double norm_step = StepNorm(stepIN, length_r, length_l);

                        if (norm_step < delta_min)
                            TrustRegionDelta = 2 * delta_min;
                        else
                            TrustRegionDelta = norm_step;

                        TrustRegionDelta = Math.Min(delta_max, TrustRegionDelta);

                        if (TrustRegionDelta < delta_min || TrustRegionDelta > delta_max)
                            throw new ArithmeticException("trust region width out of allowed range");
                        PointOnDogleg(TrustRegionDelta, step);



                        // check and adapt trust region
                        // ============================
                        double last_pred;
                        double last_ared;
                        try
                        {
                            last_pred = predictedMerit(1.0, step, beta);
                            last_ared = actualMerit(step, 1.0);
                        }
                        catch
                        {
                            last_ared = 1000000;
                            last_pred = -1;
                        }

                        //Console.WriteLine("delta start=" + TrustRegionDelta + " Ares:=" + last_ared +" Pres:=" + last_pred);
                        // trust region adaptation loop
                        while (last_ared >= last_pred)
                        {
                            double newTrustRegionDelta = TrustRegionDelta * 0.5;
                            if (newTrustRegionDelta <= delta_min)
                                break;
                            //compute new step         
                            PointOnDogleg(TrustRegionDelta, step);
                            TrustRegionDelta = Math.Max(delta_min, newTrustRegionDelta);
                            try
                            {
                                last_ared = actualMerit(step, 1.0);
                                last_pred = predictedMerit(1.0, step, beta);
                            }
                            catch
                            {
                                last_ared = 1000000;
                                last_pred = -1;
                                if (TrustRegionDelta == delta_min)
                                {
                                    throw new Exception("TrustregionDeltaNeeded to be shortened to much");
                                }
                            }
                            //Console.WriteLine("try delta=" + TrustRegionDelta + " Ares:=" + last_ared +" Pres:=" + last_pred);
                        }
                        succes = AccumulateStep(step, 1.0);
                        // //Final Trust-region update
                        // {
                        //     // taken from [Pawlovski et.al. 2006], section 2.4;
                        //     // originally from J. E. Dennis, Jr. and R. B. Schnabel. Numerical Methods for Unconstrained Optimization and Nonlinear Equations. Series in Automatic Computation. Prentice-Hall, Englewood Cliffs, NJ, 1983.

                        //     const double rho_s = 0.1;
                        //     const double rho_e = 0.75;
                        //     const double beta_s = 0.25;
                        //     const double beta_e = 4.0;

                        //     if(last_ared / last_pred < rho_s && l2_stepIN < TrustRegionDelta) {
                        //         TrustRegionDelta = Math.Max(l2_stepIN, delta_min);
                        //     } else if(last_ared / last_pred < rho_s) {
                        //         TrustRegionDelta = Math.Max(beta_s * TrustRegionDelta, delta_min); // shrinking
                        //     } else if(last_ared / last_pred > rho_e) {
                        //         TrustRegionDelta = Math.Min(beta_e * TrustRegionDelta, delta_max); // enhancing
                        //     }
                        // }
                        //stepIN.SetV(step);
                        stepIN.SetV(step);
                        return TrustRegionDelta;
                    #endregion
                    #region Line Search
                    case GlobalizationStrategy.LineSearch:

                        //Console.WriteLine("Starting Line Search:");
                        //m_alpha = Control.Alpha_Start;
                        while (m_alpha > alpha_min)
                        {

                            last_ared = actualMerit(stepIN, m_alpha);
                            last_pred = predictedMerit(m_alpha, stepIN, beta);
                            Console.WriteLine("alpha, M, predM  " + m_alpha + ", " + last_ared + ", " + last_pred);
                            if (last_ared > double.MaxValue)
                            {
                                Console.Write("");
                            }
#if DEBUG
                        //Plot and Delete the step
                        //var success = AccumulateStep(stepIN, m_alpha);
                        //if(Control.AgglomerationThreshold > 0) {
                        //    if(MultiphaseAgglomerator.TotalNumberOfAgglomerations > 0) {
                        //        TransformFromAggToSourceSpace();    //go back into non agglomerated Space
                        //    }
                        //}
                        //PlotCurrentState(Control.ProjectName + "_" + (CurrentStepNo + 1));
                        //var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                        //foreach(var pltFile in dir.GetFiles("*" + Control.ProjectName + "_" + (CurrentStepNo + 1) + ".plt")) {
                        //    pltFile.Delete();
                        //}
                        //resetStep(); // by resetting we also transform from source to agg
#endif
                            if (last_ared < last_pred)
                            {
                                succes = AccumulateStep(stepIN, m_alpha);
                                stepIN.ScaleV(m_alpha);
                                return m_alpha;
                            }
                            else
                            {
                                //reset (only if the loop continues)
                                resetStep();
                                if (m_alpha * tau > alpha_min)
                                {
                                    m_alpha = m_alpha * tau; //continues loop
                                }
                                else
                                {
                                    m_alpha = alpha_min;
                                    succes = AccumulateStep(stepIN, m_alpha);
                                    stepIN.ScaleV(m_alpha);
                                    return m_alpha;
                                }
                            }
                        }
                        //This is the case if while loop doesn't return (so if m_alpha < alpha_min)
                        m_alpha = alpha_min;
                        succes = AccumulateStep(stepIN, m_alpha);
                        stepIN.ScaleV(m_alpha);
                        return m_alpha;
                    #endregion
                    #region ZickZack Search
                    case GlobalizationStrategy.ZickZackSearch:
                        double FuncPhi(double m_alpha, bool biggerAlphaU)
                        {
                            double alphaPhi = m_alpha;
                            if (biggerAlphaU)
                            {
                                return Math.Pow(alphaPhi, 5);
                            }
                            else
                            {
                                return Math.Pow(alphaPhi, 0.2);
                            }
                        }
                        double FuncU(double m_alpha, bool biggerAlphaU)
                        {
                            double alphaU = m_alpha;
                            if (biggerAlphaU)
                            {
                                return Math.Pow(alphaU, 0.2);
                            }
                            else
                            {
                                return Math.Pow(alphaU, 5);
                            }
                        }

                        void ZickZack(double[] step, double m_alpha, bool biggerAlphaU)
                        {

                            double alpha_u = FuncU(m_alpha, biggerAlphaU);
                            double alpha_phi = FuncPhi(m_alpha, biggerAlphaU);
                            for (int i = 0; i < length_r; i++)
                            {
                                step[i] = step[i] * alpha_u;
                            }
                            for (int i = length_r; i < length_l; i++)
                            {
                                step[i] = step[i] * alpha_phi;
                            }
                        }
                        bool biggerAlphaU = false;

                        while (m_alpha > alpha_min)
                        {
                            biggerAlphaU = !biggerAlphaU;
                            //Compute Step
                            step.SetV(stepIN);
                            ZickZack(step, m_alpha, biggerAlphaU);
                            last_ared = actualMerit(step, 1.0);
                            last_pred = predictedMerit(1.0, step, beta);

                            if (last_ared < last_pred)
                            {
                                succes = AccumulateStep(step, 1.0);
                                stepIN.SetV(step);
                                return m_alpha;
                            }
                            else
                            {
                                //reset (only if the loop continues)
                                resetStep();
                                if (m_alpha * tau > alpha_min)
                                {
                                    m_alpha = m_alpha * tau; //continues loop
                                }
                                else
                                {
                                    m_alpha = alpha_min;
                                    succes = AccumulateStep(step, 1.0);
                                    stepIN.SetV(step);
                                    return m_alpha;
                                }
                            }
                        }

                        m_alpha = alpha_min;
                        step.SetV(stepIN, m_alpha);
                        //ZickZack(step,m_alpha,false);

                        succes = AccumulateStep(step, 1.0);
                        stepIN.SetV(step);
                        return m_alpha;
                    #endregion
                    #region NoGlob
                    case GlobalizationStrategy.None:
                        step.SetV(stepIN);
                        step.ScaleV(m_alpha);
                        succes = AccumulateStep(step, 1.0);
                        stepIN.SetV(step);
                        return m_alpha;
                    #endregion
                    default:
                        throw new NotSupportedException("Problem with Globalization Strategy");
                }
            }
        }

        /// <summary>
        /// This function shall update all derived Variables that are added. User has maximum flexibility here
        /// </summary>
        public abstract void UpdateDerivedVariables();
        
        /// <summary>
        /// This reinitializes all Cells prescribed
        /// </summary>
        /// <param name="isExcessiveLineSearch"></param>
        /// <param name="K_reinit">Cells to be Reinitialized</param>
        public void Reinitialize(List<Tuple<int, SpeciesId>> K_reinit) {

            List<Tuple<int, SpeciesId>> K_reinitAndNeighbours = new List<Tuple<int, SpeciesId>>();
            var ConservativeFieldsCopy = new XDGField[ConservativeFields.Length];
            for(int i = 0; i < ConservativeFields.Length; i++) {
                ConservativeFieldsCopy[i] = ConservativeFields[i].CloneAs();
            }

            // Next we compute for every cell and field a patch and then the average over that patch + we lastly set the 
            var fieldZero = ConservativeFields[0]; // on this one the jump norms are computed
            foreach(Tuple<int, SpeciesId> cell in K_reinit) {
                // gives the CutCellVolumes
                var metrics = LsTrk.GetXDGSpaceMetrics(new[] { cell.Item2 }, GetGlobalQuadOrder());
                MultidimensionalArray cutCellVolumes;
                metrics.CutCellMetrics.CutCellVolumes.TryGetValue(cell.Item2, out cutCellVolumes);

                // Compute The patch of a cell over witch we will average
                // will store the patch above which the average is computed
                double accPatchVolume = cutCellVolumes[cell.Item1];
                BitArray CellPatch = new BitArray(GridData.iGeomCells.Count);
                {
                    //cell itself belongs to CellPatch
                    CellPatch[cell.Item1] = true;
                    //unused connecting entities (output of getCellNeighbours)
                    int[] con;
                    int[] CellNeighbours;
                    var shadowZero = fieldZero.GetSpeciesShadowField(cell.Item2);
                    //gives the cell neighbors
                    this.GridData.GetCellNeighbours(cell.Item1, GetCellNeighbours_Mode.ViaEdges, out CellNeighbours, out con);
                    //in this loop over all neighbors the patch is created above which the average is computed
                    for(int iNeigh = 0; iNeigh < CellNeighbours.Length; iNeigh++) {
                        int neighbour = CellNeighbours[iNeigh];
                        //Here for each neighbor we check the jumpNorm, for this we need an edge mask
                        //var cellEdges = new int[] { GridData.CellToEdge(cell.Item1, 0), GridData.CellToEdge(cell.Item1, 1), GridData.CellToEdge(cell.Item1, 2), GridData.CellToEdge(cell.Item1, 3) };
                        var edge = new BitArray(GridData.iLogicalEdges.Count);
                        {
                            edge[con[iNeigh]] = true;
                        }
                        // ToDo: Check the JumpNorm
                        var jumpnorm = shadowZero.JumpNorm(new EdgeMask(GridData, edge));

                        // if jump-norm less than prescribed constant, add cell to patch and volume to accpatchvolume
                        if(jumpnorm <= Control.reInit_c3) {
                            CellPatch[neighbour] = true;
                            accPatchVolume += cutCellVolumes[neighbour];
                        }
                    }

                }
                CellMask PatchMask = new CellMask(GridData, CellPatch);

                // next for every field in ConFields we compute the average and reset with it the respective copy field
                for(int i = 0; i < ConservativeFields.Length; i++) {
                    var copyField = ConservativeFieldsCopy[i]; //here the new values are saved
                    var conservativeField = ConservativeFields[i];
                    var shadow = conservativeField.GetSpeciesShadowField(cell.Item2);
                    var shadowCopy = copyField.GetSpeciesShadowField(cell.Item2);

                    // Compute the average value on the patch in the original Field
                    double average = 0;
                    {
                        var SchemeHelper = metrics.XQuadSchemeHelper;
                        CellQuadratureScheme cqs = SchemeHelper.GetVolumeQuadScheme(cell.Item2, UseDefaultFactories: true, PatchMask);
                        var rule = cqs.Compile(this.GridData, GetGlobalQuadOrder());
                        CellQuadrature.GetQuadrature(new int[] { 1 },
                            this.GridData,
                            rule,
                            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                                shadow.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                            },
                            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                //for(int i = 0; i < ResultsOfIntegration.Lengths[0]; i++) {
                                average += ResultsOfIntegration.ExtractSubArrayShallow(-1, 0).Sum(); //sum over all cells
                                //}
                            }).Execute();

                        average = average / accPatchVolume; // finally divide by the average
                    }

                    //set DG coordinates of shadowCopy to zero in given cell
                    for(int iDeg = 0; iDeg <= shadow.Basis.Degree; iDeg++) {
                        foreach(int jIndex in shadow.Basis.GetPolynomialIndicesForDegree(cell.Item1, iDeg)) {
                            shadowCopy.Coordinates[cell.Item1, jIndex] = 0;
                        }
                    }
                    //set the average 
                    shadowCopy.SetMeanValue(cell.Item1, average);
#if DEBUG // Sanity check if average is higher or lower extremal value on Patch
                    //double[] minOverPatch = new double[PatchMask.NoOfItemsLocally];
                    //double[] maxOverPatch = new double[PatchMask.NoOfItemsLocally];
                    //shadow.GetCellwiseExtremalValues(minOverPatch, maxOverPatch, PatchMask, (x) => x);
                    //var max = maxOverPatch.Max();
                    //var min = minOverPatch.Min();
                    //if (average < min || average >max) {
                    //    Console.WriteLine("oh No");
                    //}
                    //Assert.GreaterOrEqual( average, min);
                    //Assert.LessOrEqual(average, max );
#endif
                }
            }

            //Overwrite the general solutionFIeld
            for(int i = 0; i < ConservativeFields.Length; i++) {
                ConservativeFields[i].CoordinateVector.SetV(ConservativeFieldsCopy[i].CoordinateVector, 1.0);

#if DEBUG
                ConservativeFieldsReinitialized[i].CoordinateVector.SetV(ConservativeFields[i].CoordinateVector, 1.0);
#endif
            }
        }
        /// <summary>
        /// Reinitializes all CutCells where the fraction $|K_s|/|K|$ is smaller than a desired size
        /// </summary>
        /// <param name="size">double between zero and 1</param>
        public void ReinitializeSmallCutCells(double size) {
            List<Tuple<int, SpeciesId>> K_reinit = new List<Tuple<int, SpeciesId>>();
            var metrics = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS, GetGlobalQuadOrder());
            var CutCellVolumes = metrics.CutCellMetrics.CutCellVolumes;

            var cutcellMask = LsTrk.Regions.GetCutCellMask();
            foreach(var iCell in cutcellMask.ItemEnum) {
                double cellsize = ((GridData)GridData).Cells.GetCellVolume(iCell);
                foreach(SpeciesId spc in SpeciesToEvaluate_Ids) {
                    if(CutCellVolumes[spc][iCell]/cellsize <= size) {
                        K_reinit.Add(new Tuple<int, SpeciesId>(iCell, spc));
                    }
                }
            }
            Reinitialize(K_reinit);
        }
        /// <summary>
        /// reinitializes all cut cells
        /// </summary>
        public void ReinitializeAllCutCells() {
            List<Tuple<int, SpeciesId>> K_reinit = new List<Tuple<int, SpeciesId>>();
            var cutcellMask = LsTrk.Regions.GetCutCellMask();
            foreach(var iCell in cutcellMask.ItemEnum) {
                foreach(SpeciesId spc in SpeciesToEvaluate_Ids) {
                    K_reinit.Add(new Tuple<int, SpeciesId>(iCell, spc));
                }
            }
            Reinitialize(K_reinit);
        }
        /// <summary>
        /// reinitializes all cut cells
        /// </summary>
        public void ReinitializeNearBand() {
            List<Tuple<int, SpeciesId>> K_reinit = new List<Tuple<int, SpeciesId>>();
            var cutcellMask = LsTrk.Regions.GetNearFieldMask(1);
            foreach(var iCell in cutcellMask.ItemEnum) {
                foreach(SpeciesId spc in SpeciesToEvaluate_Ids) {
                    K_reinit.Add(new Tuple<int, SpeciesId>(iCell, spc));
                }
            }
            Reinitialize(K_reinit);
        }
        /// <summary>
        /// Performs Reinitialization of the XDG-solution. That is setting the solution in oscillatory cells to constant values
        /// equaling the average above neighboring cells (that do not cross a shock).
        /// This function is used before each SQP-iteration and occasionally, if there is an excessive number
        /// of Line Search iterations 
        /// </summary>
        /// <param name="isExcessiveLineSearch">true if performed during line search</param>
        public void Reinitialization(bool isExcessiveLineSearch) {
            using (new FuncTrace())
            {
                //Compute the personField storing the values of the person sensor for each cell
                GetPerssonSensor(isExcessiveLineSearch);

                XDGField solutionField = ConservativeFields[0];
                CellMask fullMask = CellMask.GetFullMask(this.LsTrk.GridDat);

                double maxSensorValue = double.MinValue;
                {
                    //gives the maximum Sensor value, only used when performed during line-search

                    if (isExcessiveLineSearch == true)
                    {
                        maxSensorValue = personField.GetMeanValue(0);
                        foreach (SpeciesId speciesId in this.SpeciesToEvaluate_Ids)
                        {
                            var personField_shadow = personField.GetSpeciesShadowField(speciesId);
                            foreach (int cell in fullMask.ItemEnum)
                            {
                                if (maxSensorValue < personField_shadow.GetMeanValue(cell))
                                {
                                    maxSensorValue = personField_shadow.GetMeanValue(cell);
                                }
                            }
                        }
                    }
                }

                // Finds the set to be Reinitialized
                List<Tuple<int, SpeciesId>> K_reinit = new List<Tuple<int, SpeciesId>>();
                {
                    foreach (SpeciesId speciesId in this.SpeciesToEvaluate_Ids)
                    {
                        var personField_shadow = personField.GetSpeciesShadowField(speciesId);
                        foreach (int cell in fullMask.ItemEnum)
                        {
                            if (isExcessiveLineSearch == true)
                            {
                                //add cell if person sensor greater than prescribed fraction of maximum sensor value
                                if (personField_shadow.GetMeanValue(cell) >= Control.reInit_c2 * maxSensorValue)
                                {
                                    if (LsTrk.Regions.IsSpeciesPresentInCell(speciesId, cell))
                                    {
                                        K_reinit.Add(new Tuple<int, SpeciesId>(cell, speciesId));
                                    }
                                }
                            }
                            else
                            {
                                //add cell if person sensor greater than prescribed constant
                                if (personField_shadow.GetMeanValue(cell) > Control.reInit_c1)
                                {
                                    if (LsTrk.Regions.IsSpeciesPresentInCell(speciesId, cell))
                                    {
                                        K_reinit.Add(new Tuple<int, SpeciesId>(cell, speciesId));
                                    }
                                }
                            }
                        }
                    }
                }

                //print all cells that are being reinitialized
                List<int> CellsToReInit = new List<int>();
                if (K_reinit.Count > 0)
                {
                    Console.Write("Reinitialized Cells:");
                }

                foreach (Tuple<int, SpeciesId> cell in K_reinit)
                {
                    CellsToReInit.Add(cell.Item1);
                    Console.Write("(" + cell.Item1 + " " + LsTrk.GetSpeciesName(cell.Item2) + ") ,");
                }
                Console.WriteLine(" ");

                Reinitialize(K_reinit);
            }

        }
        /// <summary>
        /// This function computes the values of the PersonPerraire sensor
        /// </summary>/// <exception cref="NotSupportedException"> if the basis degree is 0 the sensor is not defined</exception>
        public void GetPerssonSensor(bool isCalledDuringLineSearch) {
            if(ConservativeFields[0].Basis.Degree == 0) {
                throw new NotSupportedException("Sensor not supported for DG degree 0");
            }
            
            XDGField xdgfieldToTest = ConservativeFields[0].CloneAs();

            //start by clearing the personField and setting all mean values to minus infinity
            personField.Clear();
            var fullMask = CellMask.GetFullMask(this.LsTrk.GridDat);
            foreach(int jCell in fullMask.ItemEnum) {
                foreach(SpeciesId speciesId in LsTrk.SpeciesIdS) {
                    personField.GetSpeciesShadowField(speciesId).SetMeanValue(jCell, -10);
                }
            }

            //Choose Quad Order
            int deg = xdgfieldToTest.Basis.Degree;
            var AvailOrders = LsTrk.GetCachedOrders().Where(order => order >= 2 * deg);
            int order2Pick = AvailOrders.Any() ? AvailOrders.Min() : 2 * deg;

            //get the p-1 Projection Up-1
            XDGField fieldPMinus1 = GetPMinus1Projection(LsTrk, xdgfieldToTest, order2Pick);
            var fieldPMinus1LaidBack = new XDGField(xdgfieldToTest.Basis, "fieldPMinus1LaidBack");
            fieldPMinus1LaidBack.AccLaidBack(1.0, fieldPMinus1);

            //do the substraction U - Up-1
            var diffField = xdgfieldToTest.CloneAs();
            diffField.Identification = "diffField";
            diffField.Acc(-1.0, fieldPMinus1LaidBack);
            int N = xdgfieldToTest.Basis.NonX_Basis.Length; // DOFs per cell per species.

            //helper functions
            double[] GetCoords(int j,SpeciesId spc,XDGField field) {
                int iSpc = LsTrk.Regions.GetSpeciesIndex(spc, j);
                if(iSpc < 0)
                    return null;
                else
                    return field.Coordinates.GetRowPart(j, iSpc * N, N);
            }
            void SetValuePersonField(int cell, SpeciesId speciesId, double nominator, double denominator) {
                if(nominator == 0) {
                    personField.GetSpeciesShadowField(speciesId).SetMeanValue(cell, -10);
                } else {
                    
                    if(denominator == 0) {
                        throw new ArgumentException("denominator zero but nominator not zero");
                    } else {
                        double PerssonSensor = (nominator / denominator).Sqrt();
                        double val = Math.Log10(PerssonSensor);
                        if(val == 0) {
                            Console.WriteLine("WTF");
                        }
                        personField.GetSpeciesShadowField(speciesId).SetMeanValue(cell, val);
                    }
                }
            }
            
            //compute the person field
            foreach(SpeciesId speciesId in this.SpeciesToEvaluate_Ids) {
                var speciesMask = LsTrk.Regions.GetSpeciesMask(speciesId);
                var cutCellMask = LsTrk.Regions.GetCutCellMask().Intersect(speciesMask);
                var UnCutSpeciesMask = speciesMask.Except(cutCellMask);

                var MMF = LsTrk.GetXDGSpaceMetrics(speciesId, order2Pick).MassMatrixFactory;
                var MMblox = MMF.GetMassMatrixBlocks(xdgfieldToTest.Basis.NonX_Basis, speciesId);

                // 1st, do all the cut cells with non-id mass matrix
                // =================================================
                double[] tmp = new double[N];
                double[] tmp2 = new double[N];
                for(int iSub = 0; iSub < MMblox.jSub2jCell.Length; iSub++) {
                    int jCell = MMblox.jSub2jCell[iSub];
                    double[] CoordsFTT = GetCoords(jCell, speciesId, xdgfieldToTest);
                    double[] CoordsDiff = GetCoords(jCell, speciesId, diffField);
                    if(CoordsFTT == null)
                       continue; // species not present in cell; no contribution.

                    var MM_j = MMblox.MassMatrixBlocks.ExtractSubArrayShallow(new int[] { iSub, 0, 0 }, new int[] { iSub - 1, N - 1, N - 1 });

                    MM_j.GEMV(1.0, CoordsFTT, 0.0, tmp);
                    double denominator = CoordsFTT.InnerProd(tmp);

                    MM_j.GEMV(1.0, CoordsDiff, 0.0, tmp2);
                    double nominator = CoordsDiff.InnerProd(tmp2);

                    SetValuePersonField(jCell, speciesId, nominator, denominator);
                    
                }

                //2nd we do the NonCutCells
                foreach(int iCell in UnCutSpeciesMask.ItemEnum) {
                    double[] CoordsFTT = GetCoords(iCell, speciesId, xdgfieldToTest);
                    double[] CoordsDiff = GetCoords(iCell, speciesId, diffField);

                    // in this iCell, we have an orthonormal basis, i.e. the mass matrix is the identity
                    double denominator = CoordsFTT.L2NormPow2();
                    double nominator = CoordsDiff.L2NormPow2();

                    SetValuePersonField(iCell, speciesId, nominator, denominator);
                }
            }
        }
        /// <summary>
        /// Gets an Operator Matrix that should project a high-order solution to a low-order solution in the L2 sense
        /// </summary>
        /// <param name="LsTrk"></param> LevelSetTracker
        /// <param name="fieldPhigh"></param> high order XDGField
        /// <param name="order2Pick"></param> Quad Order
        /// <param name="pOrder"></param> low porder polynomail degree
        /// <returns></returns>
        MsrMatrix GetProjectionOperator(LevelSetTracker LsTrk, XDGField[] fieldPhigh, int order2Pick, int pOrder)
        {
            //create a xdgfield of deg p and get some mappings
            var nVar = fieldPhigh.Length;
            var fieldP = new XDGField[fieldPhigh.Length];
            for (int i = 0; i < nVar; i++)
            {
                fieldP[i] = new XDGField(new XDGBasis(LsTrk, pOrder), "fieldP");
            }
            var highMap = new CoordinateMapping(fieldPhigh);
            var lowMap = new CoordinateMapping(fieldP);

            //get the MassMatrix
            MassMatrixFactory massMatrixFactory = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, order2Pick).MassMatrixFactory;
            BlockMsrMatrix massMatrix = massMatrixFactory.GetMassMatrix(highMap, inverse: false);
            BlockMsrMatrix MPlowPlow = massMatrixFactory.GetMassMatrix(lowMap, inverse: true);

            //get the subMassMatrices of deg pOrder
            MsrMatrix MPlowP;//Sub of Mass Mat with (rows pOrder ,column p)
            {
                //compute the lists of indices used by GetSubMatrix()
                var rowIndices = new List<long>();
                var colIndices = new List<long>();
                {
                    int MaxModeR = fieldPhigh[0].Mapping.MaxTotalNoOfCoordinatesPerCell / LsTrk.TotalNoOfSpecies;
                    int MaxModer = fieldP[0].Mapping.MaxTotalNoOfCoordinatesPerCell / LsTrk.TotalNoOfSpecies;

                    int iField;
                    int jCell;
                    int nMode;
                    for (int iRow = 0; iRow < lowMap.TotalLength; iRow++)
                    {
                        lowMap.LocalFieldCoordinateIndex(iRow, out iField, out jCell, out nMode);
                        double rMode = (double)nMode;
                        int fac = (int)Math.Floor(rMode / MaxModer);
                        int row = highMap.LocalUniqueCoordinateIndex(iField, jCell, nMode + fac * (MaxModeR - MaxModer));
                        rowIndices.Add(row);
                    }
                    for (int i = 0; i < massMatrix.NoOfCols; i++)
                    {
                        colIndices.Add(i);
                    }
                }
                //obtain the submatrix
                MPlowP = massMatrix.ToMsrMatrix().GetSubMatrix(rowIndices, colIndices);
            }

            // Here we compute M_{p-1,p-1}^{-1} * M_{p-1,p}
            return MsrMatrix.Multiply(MPlowPlow.ToMsrMatrix(), MPlowP);
            
        }

        /// <summary>
        /// Gives an inclusion operator that transforms a low order XDGField into a high-order XDGField by copying the coordinates
        /// A*c_low=c_high
        /// </summary>
        /// <param name="fieldPhigh"></param>
        /// <param name="pOrder">order of c_low</param>
        /// <returns></returns>
        public MsrMatrix GetInclusionOperator(XDGField[] fieldPhigh, int pOrder)
        {
            var nVar = fieldPhigh.Length;
            var fieldPlow = new XDGField[fieldPhigh.Length];
            var copyHigh = new XDGField[fieldPhigh.Length]; 
            for (int i = 0; i < nVar; i++)
            {
                fieldPlow[i] = new XDGField(new XDGBasis(LsTrk, pOrder), "fieldP");
                copyHigh[i] = fieldPhigh[i].CloneAs();
                copyHigh[i].Clear();
            }
            var cV_low = new CoordinateVector(fieldPlow);
            var highMap = new CoordinateMapping(copyHigh);
            var lowMap = new CoordinateMapping(fieldPlow);
            // numbering of entries
            for (int iRow = 0; iRow < lowMap.TotalLength; iRow++)
            {
                cV_low[iRow] = iRow + 1;
            }
            for (int i = 0; i < nVar; i++)
            {
                copyHigh[i].AccLaidBack(1.0, fieldPlow[i]);
            }

            var cV_high = new CoordinateVector(copyHigh);
            MsrMatrix Inclusion = new MsrMatrix((int) highMap.TotalLength, (int) lowMap.TotalLength,1,1) ;//Sub of Mass Mat with (rows pOrder ,column p)
            {
                for (int iRow = 0; iRow < highMap.TotalLength; iRow++)
                {
                    int lowRow = (int) cV_high[iRow];
                    if (lowRow> 0)
                    {
                        Inclusion[iRow, lowRow-1]=1;
                    }
                }
            }
            return Inclusion;

        }
        /// <summary>
        /// Computes the Projection onto the (p-1)-Polynomial space (relative to the input field)
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <param name="xdgfieldToTest"></param>
        /// <param name="order2Pick"></param>
        /// <returns>projection</returns>
        public static XDGField GetPMinus1Projection(LevelSetTracker LsTrk, XDGField xdgfieldToTest ,int order2Pick) {
            var fieldPMinus1 = new XDGField(new XDGBasis(LsTrk, xdgfieldToTest.Basis.Degree - 1), "fieldPMinus1");
            {
                //get the MassMatrix
                MassMatrixFactory massMatrixFactory = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, order2Pick).MassMatrixFactory;
                BlockMsrMatrix massMatrix = massMatrixFactory.GetMassMatrix(xdgfieldToTest.Mapping, inverse: false);
                BlockMsrMatrix MMPmin1Pmin1Inv = massMatrixFactory.GetMassMatrix(fieldPMinus1.Mapping, inverse: true);

                //get the subMassMatrices of deg p-1
                MsrMatrix MMPmin1P;//Sub of Mass Mat with (rows p-1 ,column p)
                {
                    //compute the lists of indices used by GetSubMatrix()
                    var rowIndices = new List<long>();
                    var colIndices = new List<long>();
                    {
                        int MaxModeR = xdgfieldToTest.Mapping.MaxTotalNoOfCoordinatesPerCell / LsTrk.TotalNoOfSpecies;
                        int MaxModer = fieldPMinus1.Mapping.MaxTotalNoOfCoordinatesPerCell  / LsTrk.TotalNoOfSpecies;

                        int iField;
                        int jCell;
                        int nMode;
                        for(int iRow = 0; iRow < fieldPMinus1.Mapping.TotalLength; iRow++) {
                            fieldPMinus1.Mapping.LocalFieldCoordinateIndex(iRow, out iField, out jCell, out nMode);
                            double rMode = (double)nMode;
                            int fac = (int)Math.Floor(rMode / MaxModer);
                            int row = xdgfieldToTest.Mapping.LocalUniqueCoordinateIndex(iField, jCell, nMode + fac * (MaxModeR - MaxModer));
                            rowIndices.Add(row);
                        }
                        for(int i = 0; i < massMatrix.NoOfCols; i++) {
                            colIndices.Add(i);
                        }
                    }
                    MMPmin1P = massMatrix.ToMsrMatrix().GetSubMatrix(rowIndices, colIndices);
                }

                // Here we compute M_{p-1,p-1}^{-1} * M_{p-1,p} * Coords(xdgfieldtoTest)
                var tmpVec = new double[fieldPMinus1.CoordinateVector.Length];
                MMPmin1P.SpMV(1.0, xdgfieldToTest.CoordinateVector, 0.0, tmpVec);
                MMPmin1Pmin1Inv.SpMV(1.0, tmpVec, 0.0, fieldPMinus1.CoordinateVector);


            }
            return fieldPMinus1;
        }

        /// <summary>
        /// Evolve the solution ConVars in Time using ExplicitTimeStepping, neither used or tested
        /// </summary>
        /// <param name="ConVars"></param>
        /// <param name="residuals"></param>
        void DoExplicitTimeStepping(XDGField[] ConVars, XDGField[] residuals) {
            //init an operator with temp operator
            var OpWithTemp = XSpatialOperator.CloneAs();
            var TempOp = new ConstantXTemporalOperator(OpWithTemp);
            OpWithTemp.TemporalOperator = TempOp;
            OpWithTemp.Commit();
            //init an appropriate timeStepper            
            var xdgtimestepping = new XdgTimestepping(
                                OpWithTemp,
                                ConVars, residuals,
                                TimeSteppingScheme.ExplicitEuler,
                                null, LevelSetHandling.None,
                                this.MultiGridOperatorConfig, this.Control.AgglomerationThreshold,
                                this.Control.LinearSolver, this.Control.NonLinearSolver);

            //get a starting Timestepsize
            double dt = Control.IG_dt_Start;
            //double dt = 0.1;
            //an adaptive function
            void AdaptTS(double nu) {
                //adaptive Time-stepping
                if(nu < Control.IG_nu_Min) {
                    dt = dt * Control.IG_alpha;
                } else if(nu < Control.IG_nu_Max) {
                    //do nothing with dt
                } else {
                    dt = dt * Control.IG_beta;
                }
            }

            //choose max iter and initialize counter
            int maxIt = 100;
            int it = 0;

            var SolVec = new CoordinateVector(ConVars);
            var u_n0 = new double[SolVec.Length];
            u_n0.SetV(SolVec);

            var ResVec = new CoordinateVector(residuals);
            Eval_r = XSpatialOperator.GetEvaluatorEx(LsTrk, ConVars, null, ResidualMap);
            Eval_r.Evaluate(1.0, 0.0, ResVec);


            var tp = new Tecplot(LsTrk.GridDat, 3);

            tp.PlotFields($"ComputeP0Solution_{it}", LsTrk.Regions.Time, new DGField[] { LevelSet, LevelSetTwo, ConVars[0], ConVars[1], ConVars[2], ConVars[3], residuals[0], residuals[1], residuals[2], residuals[3] });
            while(ResVec.MPI_L2Norm() > 1e-5 && it < maxIt && dt <= Math.Max(1e06, Control.IG_dt_Start)) {
                //xdgtimestepping.Solve(tmp_trk.Regions.Time, dt);

                try {

                    var succes = xdgtimestepping.Solve(LsTrk.Regions.Time, dt);
                } catch {
                    SolVec = new CoordinateVector(ConVars);
                    SolVec.SetV(u_n0);
                    Console.WriteLine($"time step dt ={dt} breaks operator and is shortened to {dt * Control.IG_beta}");
                    dt = dt * Control.IG_beta;
                    continue;
                }


                double[] delta_u = new double[SolVec.Count];
                SolVec = new CoordinateVector(ConVars);
                delta_u.AccV(1.0, SolVec);
                delta_u.AccV(-1.0, u_n0);
                var norm = delta_u.L2Norm();
                double MaxAbs(double[] u) {
                    for(int i = 0; i < u.Length; i++) {
                        u[i] = u[i].Abs();
                    }
                    return u.Max();
                }


                delta_u.ScaleV(1 / Math.Max(MaxAbs(u_n0), 1e-3));
                var nu = MaxAbs(delta_u);

                Eval_r = XSpatialOperator.GetEvaluatorEx(LsTrk, ConVars, null, ResidualMap);
                ResVec = new CoordinateVector(residuals);
                Eval_r.Evaluate(1.0, 0.0, ResVec);
                Console.WriteLine($"It {it}: ||r||={ResVec.MPI_L2Norm()}, ||du||={norm}, ||nu|| = {nu}, dt={dt}");


                AdaptTS(nu);
                u_n0.SetV(SolVec);
                it++;
                tp.PlotFields($"ComputeP0Solution_{it}", LsTrk.Regions.Time, new DGField[] { LevelSet, LevelSetTwo, ConVars[0], ConVars[1], ConVars[2], ConVars[3], residuals[0], residuals[1], residuals[2], residuals[3] });

            }

            Console.WriteLine("...Finished");

        }
        /// <summary>
        /// Evolves the solution in time using implicit Euler time stepping and adaptive time steps
        /// Goal: find steady state
        /// </summary>
        /// <param name="ConVars"> Conservative Variables </param>
        /// <param name="residuals">Residual fields</param>
        void DoImplicitTimeStepping(XDGField[] ConVars, XDGField[] residuals) {
            
            //init an operator with temp operator
            var OpWithTemp = XSpatialOperator.CloneAs();
            var TempOp = new ConstantXTemporalOperator(OpWithTemp);
            OpWithTemp.TemporalOperator = TempOp;
            OpWithTemp.Commit();

            //init an appropriate timeStepper
            var xdgtimestepping = new XdgTimestepping(OpWithTemp, ConVars, residuals,TimeSteppingScheme.RK_ImplicitEuler, _AgglomerationThreshold:CurrentAgglo);


            //get a starting Timestepsize
            double dt = Control.IG_dt_Start;

            //function to adapt timestep (using some predefined threshholds)
            void AdaptTS(double nu) {
                //adaptive Time stepping
                if(nu < Control.IG_nu_Min) {
                    dt = dt * Control.IG_alpha;
                } else if(nu < Control.IG_nu_Max) {
                    //do nothing with dt
                } else {
                    dt = dt * Control.IG_beta;
                }
            }
            
            //choose max iter and initialize counter
            int maxIt = 200;
            int it = 0; // counter

            //get resiudal and solution as vector
            var SolVec = new CoordinateVector(ConVars);
            var ResVec = new CoordinateVector(residuals);

            //backup for solution
            var u_n0 = new double[SolVec.Length];
            u_n0.SetV(SolVec);
            
            // eval inital residual
            Eval_r = XSpatialOperator.GetEvaluatorEx(LsTrk, ConVars, null, ResidualMap);
            Eval_r.Evaluate(1.0, 0.0, ResVec);

            //agglomerate residual
            MultiphaseAgglomerator = LsTrk.GetAgglomerator(SpeciesToEvaluate_Ids, GetGlobalQuadOrder(), CurrentAgglo, ExceptionOnFailedAgglomeration: false);
            MultiphaseAgglomerator.ManipulateMatrixAndRHS(default(MsrMatrix), ResVec, ResidualMap, new CoordinateMapping(ConVars));
            
            // while residual big, time step small and it < maxIt do timestepping
            while(ResVec.MPI_L2Norm() > 1e-12 && it < maxIt && dt <= Math.Max(1e06, Control.IG_dt_Start)) {
                
                // try implicit timestep, if fail make dt smaller
                try {
                    var succes = xdgtimestepping.Solve(LsTrk.Regions.Time, dt);
                } catch {
                    //reset solution
                    SolVec = new CoordinateVector(ConVars);
                    SolVec.SetV(u_n0);
                    Console.WriteLine($"time step dt ={dt} breaks operator and is shortened to {dt * Control.IG_beta}");
                    dt = dt * Control.IG_beta;
                    continue;
                }

                //obtain the step (delta u) and calc norm
                double[] delta_u = new double[SolVec.Count];
                SolVec = new CoordinateVector(ConVars);
                delta_u.AccV(1.0, SolVec);
                delta_u.AccV(-1.0, u_n0);
                var norm = delta_u.L2Norm();

                // infinity norm
                double MaxAbs(double[] u) {
                    for(int i = 0; i < u.Length; i++) {
                        u[i] = u[i].Abs();
                    }
                    return u.Max();
                }

                //normalize delta u by 1/|u|_infty or 1e-3
                delta_u.ScaleV(1 / Math.Max(MaxAbs(u_n0), 1e-3));
                var nu = MaxAbs(delta_u);

                //compute agglomerated residual
                Eval_r = XSpatialOperator.GetEvaluatorEx(LsTrk, ConVars, null, ResidualMap);
                ResVec = new CoordinateVector(residuals);
                Eval_r.Evaluate(1.0, 0.0, ResVec);
                MultiphaseAgglomerator.ManipulateMatrixAndRHS(default(MsrMatrix), ResVec, ResidualMap, new CoordinateMapping(ConVars));
                
                // print result
                Console.WriteLine($"It {it}: ||r||={ResVec.MPI_L2Norm()}, ||du||={norm}, ||nu|| = {nu}, dt={dt}");

                //adapt time step
                AdaptTS(nu);
                u_n0.SetV(SolVec);
                it++;
                //tp.PlotFields($"ComputeP0Solution_{it}", LsTrk.Regions.Time, new DGField[] { LevelSet, LevelSetTwo, ConVars[0], ConVars[1], ConVars[2], ConVars[3], residuals[0], residuals[1], residuals[2], residuals[3] });

            }

            Console.WriteLine("...Finished");
            
        }
        /// <summary>
        /// Compute a p0 initial guess for the sate given a fixed LevelSet
        /// </summary>
        public void ComputeP0Solution() {
            using (new FuncTrace())
            {
                //create a vector with the cons. fields with p=0 basis
                var P0Vars = new XDGField[ConservativeFields.Length];
                var P0Residuals = new XDGField[ConservativeFields.Length];
                //LevelSetTracker tmp_trk;
                //if(Control.IsTwoLevelSetRun) {
                //    tmp_trk = new LevelSetTracker(LsTrk.GridDat, LsTrk.CutCellQuadratureType, LsTrk.NearRegionWidth, LsTrk.SpeciesTable, this.LevelSet, this.LevelSetTwo);
                //} else {
                //    tmp_trk = new LevelSetTracker(LsTrk.GridDat, LsTrk.CutCellQuadratureType, LsTrk.NearRegionWidth, LsTrk.SpeciesTable, this.LevelSet);
                //}

                //var p0Basis = new XDGBasis(tmp_trk, 0);
                var p0Basis = new XDGBasis(LsTrk, 0);
                for (int iField = 0; iField < ConservativeFields.Length; iField++)
                {
                    P0Vars[iField] = new XDGField(p0Basis);
                    P0Vars[iField].Identification = ConservativeFields[iField].Identification;
                    P0Residuals[iField] = new XDGField(p0Basis);
                    P0Residuals[iField].Identification = Residuals[iField].Identification + "_res";
                    //Projects the Initial Value onto p0vars
                    foreach (SpeciesId spc in this.SpeciesToEvaluate_Ids)
                    {
                        P0Vars[iField].GetSpeciesShadowField(spc).ProjectFromForeignGrid(1.0, (ConventionalDGField)ConservativeFields[iField].GetSpeciesShadowField(spc)); ;
                        P0Residuals[iField].GetSpeciesShadowField(spc).ProjectFromForeignGrid(1.0, (ConventionalDGField)Residuals[iField].GetSpeciesShadowField(spc));
                    }
                }
                DoImplicitTimeStepping(P0Vars, P0Residuals);


                Console.WriteLine("...Finished");

                //Write p0Vars onto Conservative Fields
                for (int iField = 0; iField < ConservativeFields.Length; iField++)
                {
                    ConservativeFields[iField].Clear();
                    ConservativeFields[iField].AccLaidBack(1.0, P0Vars[iField]);
                }
            }
        }
        /// <summary>
        /// Assembles Jacobian and Residual dependent on which linearization is chosen. This routine used in P0Solution Computation
        /// </summary>
        /// <param name="DoVars"></param>
        /// <param name="CoDoVars"></param>
        /// <param name="Jacobian"></param>
        /// <param name="ResidualVec"></param>
        /// <exception cref="ArgumentException"></exception>
        public void ComputeJacobianAndRes(XDGField[] DoVars, XDGField[] CoDoVars, MsrMatrix Jacobian, CoordinateVector ResidualVec) {

            Jacobian = Jacobian==null? new MsrMatrix(new CoordinateMapping(DoVars), new CoordinateMapping(CoDoVars)) : Jacobian;
            ResidualVec = ResidualVec == null ? new CoordinateVector(CoDoVars):ResidualVec;
         
            
            switch(XSpatialOperator.LinearizationHint) {
                case LinearizationHint.FDJacobi:
                //var r_JacobianBuilder = XSpatialOperator.GetFDJacobianBuilder(ConservativeFields, null, ResidualMap);
                //r_JacobianBuilder.ComputeMatrix(Jr_U, ResidualVector);

                var R_JacobianBuilder = XSpatialOperator.GetFDJacobianBuilder(DoVars, null, new CoordinateMapping(CoDoVars));
                R_JacobianBuilder.ComputeMatrix(Jacobian, ResidualVec);
                break;
                case LinearizationHint.GetJacobiOperator:

                //var r_MatrixBuilder = r_JacobiOperator.GetMatrixBuilder(new CoordinateMapping(ConservativeFields), r_JacobiOperator.InvokeParameterFactory(ConservativeFields), ResidualMap);
                //r_MatrixBuilder.ComputeMatrix(Jr_U, ResidualVector);

                var R_MatrixBuilder = R_JacobiOperator.GetMatrixBuilder(new CoordinateMapping(DoVars), r_JacobiOperator.InvokeParameterFactory(DoVars), new CoordinateMapping(CoDoVars));
                R_MatrixBuilder.ComputeMatrix(Jacobian, ResidualVec);
                break;
                case LinearizationHint.AdHoc:
                //var r_MatrixBuilderAdhoc = XSpatialOperator.GetMatrixBuilder(new CoordinateMapping(ConservativeFields), r_JacobiOperator.InvokeParameterFactory(ConservativeFields), ResidualMap);
                //r_MatrixBuilderAdhoc.ComputeMatrix(Jr_U, ResidualVector);

                var R_MatrixBuilderAdhoc = XSpatialOperator.GetMatrixBuilder(new CoordinateMapping(DoVars), r_JacobiOperator.InvokeParameterFactory(DoVars), new CoordinateMapping(CoDoVars));
                R_MatrixBuilderAdhoc.ComputeMatrix(Jacobian, ResidualVec);
                break;
                default:
                throw new ArgumentException("this should not happen");
            }
        }
        /// <summary>
        /// evaluates the Quad order Function specified in the Control
        /// </summary>
        /// <returns></returns>
        public virtual int GetGlobalQuadOrder() {
            int[] DomainDegrees = ResidualMap.BasisS.Select(f => f.Degree).ToArray();
            int[] CodomainDegrees = obj_f_map.BasisS.Select(f => f.Degree).ToArray();
            var quadOrder = Control.quadOrderFunc(DomainDegrees,null,CodomainDegrees);
            return quadOrder;
        }
        /// <summary>
        /// Gives the Residual Plot using the current operator (so the test functions are the same for all evaluations), makes a difference in staggered setting 
        /// </summary>
        /// <param name="si"></param>
        /// <returns></returns>
        public Plot2Ddata GetUniformResPlot(ISessionInfo si) {

            List<double> residuals = new List<double>();
            List<double> enrichedResiduals = new List<double>();

            foreach(var timestep in si.Timesteps) {
                List<DGField> flds = new List<DGField>();
                var fields = timestep.Fields.ToList();
                LsTBO.Clear();
                if(Control.IsTwoLevelSetRun) {
                    LsTBO.CoordinateVector.SetV(timestep.GetField("levelSetTwo").CoordinateVector);
                } else {
                    LsTBO.CoordinateVector.SetV(timestep.GetField("levelSet").CoordinateVector);
                }
                LsTrk.UpdateTracker(0.0);

                //load the fields
                var l = ConservativeFields.Length;
                XDGField[] ConsFieldsTS = new XDGField[l];
                if(l==1) {
                    ConsFieldsTS[0] = new XDGField(new XDGBasis(LsTrk, timestep.GetField("c").Basis.Degree));
                    ConsFieldsTS[0].CoordinateVector.SetV(timestep.GetField("c").CoordinateVector);
                } else {
                    ConsFieldsTS[0] = new XDGField(new XDGBasis(LsTrk, timestep.GetField("rho").Basis.Degree));
                    ConsFieldsTS[1] = new XDGField(new XDGBasis(LsTrk, timestep.GetField("rho").Basis.Degree));
                    ConsFieldsTS[2] = new XDGField(new XDGBasis(LsTrk, timestep.GetField("rho").Basis.Degree));
                    ConsFieldsTS[3] = new XDGField(new XDGBasis(LsTrk, timestep.GetField("rho").Basis.Degree));
                    ConsFieldsTS[0].CoordinateVector.SetV(timestep.GetField("rho").CoordinateVector);
                    ConsFieldsTS[1].CoordinateVector.SetV(timestep.GetField("m0").CoordinateVector);
                    ConsFieldsTS[2].CoordinateVector.SetV(timestep.GetField("m1").CoordinateVector);
                    ConsFieldsTS[3].CoordinateVector.SetV(timestep.GetField("rhoE").CoordinateVector);
                }
                for(int i = 0; i < ConservativeFields.Length; i++) {
                    ConservativeFields[i].Clear();
                    ConservativeFields[i].AccLaidBack(1.0, ConsFieldsTS[i]);
                }

                (res_l2, obj_f, res_L2) = ComputeResiduals();

                residuals.Add(res_l2);
                enrichedResiduals.Add(obj_f);
            }
            return GetPlot(enrichedResiduals, "||R(z)||_2", residuals, "||r(z)||_2");
        }
        /// <summary>
        /// returns a plot against the iterations for variable lists
        /// </summary>
        /// <param name="y1"></param>
        /// <param name="l1"></param>
        /// <param name="y2"></param>
        /// <param name="l2"></param>
        /// <param name="y3"></param>
        /// <param name="l3"></param>
        /// <returns></returns>
        public Plot2Ddata GetPlot(List<double> y1, string l1, List<double> y2=null, string l2=null,List<double> y3 = null, string l3=null) {
            var plot = new Plot2Ddata();
            var Fmt = new PlotFormat();
            Fmt.PointType = PointTypes.OpenCircle;
            //Fmt.PointSize = 1.2;
            //Fmt.LineWidth = 1;
            Fmt.Style = Styles.LinesPoints;
            Fmt.LineColor = LineColors.Blue;
            Fmt.PointType = PointTypes.Diamond;

            var Fmt2 = Fmt.CloneAs();
            Fmt2.PointType = PointTypes.Asterisk;
            Fmt2.LineColor = LineColors.Red;

            var Fmt3 = Fmt.CloneAs();
            Fmt3.PointType = PointTypes.Plus;
            Fmt3.LineColor = LineColors.Black;

            plot.AddDataGroup(l1, StepCount, y1, Fmt);
            if(y2 != null && l2!=null) {
                plot.AddDataGroup(l2, StepCount, y2, Fmt2);
            }
            if(y3 != null && l3 != null) {
                plot.AddDataGroup(l3, StepCount, y3, Fmt3);
            }

            plot.Xlabel = "Iteration";
            plot.LogY = true;

            return plot;
        }
        /// <summary>
        /// Gives a PlotData where the Residual history is shown as computed by the method (test functions may be different if a basis change was done)
        /// </summary>
        /// <returns></returns>
        public Plot2Ddata GetResPlot() {
            GetPlotTable(obj_f_vals, "R1norms", ResNorms, "R0norms");
            return GetPlot(obj_f_vals, "||R(z)||", ResNorms, "||r(z)||");
        }
        public Plot2Ddata GetResEnthalpyPlot(ISessionInfo si, double EEN) {
            var EE = new List<double>();
            foreach(var timestep in si.Timesteps) {
                EE.Add(((XDGField)timestep.GetField("h_err")).L2NormAllSpecies() / EEN);
            }
            GetPlotTable(obj_f_vals, "obj_f_vals", ResNorms, "ResNorms", EE, "EnthalpyError");
            return GetPlot(obj_f_vals, "||R(z)||", ResNorms, "||r(z)||", EE, "||h_{err}||");
        }
        /// <summary>
        /// Used for Ticks
        /// </summary>
        /// <param name="y1"></param>
        /// <param name="l1"></param>
        /// <param name="y2"></param>
        /// <param name="l2"></param>
        /// <param name="y3"></param>
        /// <param name="l3"></param>
        public void GetPlotTable(List<double> y1, string l1, List<double> y2 = null, string l2 = null, List<double> y3 = null, string l3 = null) {
            SavePlotTable(StepCount, y1, l1, y2, l2, y3, l3);
        }

        /// <summary>
        /// Used for Ticks
        /// </summary>
        /// <param name="y1"></param>
        /// <param name="l1"></param>
        /// <param name="y2"></param>
        /// <param name="l2"></param>
        /// <param name="y3"></param>
        /// <param name="l3"></param>
        public static void SavePlotTable(List<double> StepCount,List<double> y1, string l1, List<double> y2 = null, string l2 = null, List<double> y3 = null, string l3 = null)
        {
            var table1 = MultidimensionalArray.Create(y1.Count, 2);
            table1.SetColumn(0, StepCount);
            table1.SetColumn(1, y1);
            table1.SaveToTextFile(l1 + ".txt");
            if (y2 != null && l2 != null)
            {
                var table2 = MultidimensionalArray.Create(y2.Count, 2);
                table2.SetColumn(0, StepCount);
                table2.SetColumn(1, y2);
                table2.SaveToTextFile(l2 + ".txt");
            }
            if (y3 != null && l3 != null)
            {
                var table3 = MultidimensionalArray.Create(y3.Count, 2);
                table3.SetColumn(0, StepCount);
                table3.SetColumn(1, y3);
                table3.SaveToTextFile(l3 + ".txt");
            }
        }



        //public Plot2Ddata GetResEnthalpyPlot(ISessionInfo si, double EEN) {
        //    var EE = new List<double>();
        //    foreach(var time step in si.Timesteps) {
        //        EE.Add(((XDGField)timestep.GetField("h_err")).L2NormAllSpecies() / EEN);
        //    }
        //    GetPlotTable(obj_f_vals, "obj_f_vals", ResNorms, "ResNorms", EE, "EnthalpyError");
        //}

        public void PlotShadowFields(ISessionInfo si) {
            var texplot = new Tecplot(ConservativeFields[0].GridDat, 3);
            foreach(var timestep in si.Timesteps) {
                List<DGField> flds = new List<DGField>();
                var fields = timestep.Fields.ToList();
                IDTTimeStepInfo lts = (IDTTimeStepInfo)((TimestepProxy)timestep).GetInternal();
                texplot.PlotFields($"{si.ProjectName}_SP_" + lts.TimeStepNumber, 0.0, lts.GetShadowFields());
            }
        }
    }
#region Helper Class to Assemble System
    /// <summary>
    /// This is a helper class serving the purpose to insert one MsrMatrix A into a bigger Matrix B.
    /// This is done by using AccSubMatrixTo - method.
    /// This method needs long - arrays (which command the row/column indices that are being set) to work. 
    /// These long arrays are initialized only once on construction of this assembler object
    /// </summary>
    public class MatrixAssembler {
        public int m_length;
        private long[] helper;
        public MatrixAssembler(int Mlength) {
            m_length = Mlength;
            this.helper = GetIndices(Mlength);
        }
        public static long[] GetIndices(int length) {
            var ret = new long[length];
            for(int i = 0; i < length; i++) {
                ret[i] = i;
            }
            return ret;
        }
        public MatrixAssembler(MsrMatrix source) {
            m_length = (int)source.NoOfCols;
            this.helper = new long[m_length];
            for(int i = 0; i < m_length; i++) {
                helper[i] = i;
            }
        }
        // method inserts source Matrix into target Matrix, at a specified position
        // iRow, jCol ---  Position in TargetMatrix where source is put
        public void AccBlock(MsrMatrix source, MsrMatrix target, int iRow, int jCol) {
            if(m_length < source.NoOfCols) {
                throw new ArgumentException("length of source bigger than target length");
            }

            var rows_source = new long[source.NoOfRows];
            Array.Copy(helper, 0, rows_source, 0, source.NoOfRows);
            var cols_source = new long[source.NoOfCols];
            Array.Copy(helper, 0, cols_source, 0, source.NoOfCols);

            var rows_target = new long[source.NoOfRows];
            Array.Copy(helper, iRow, rows_target, 0, source.NoOfRows);
            var cols_target = new long[source.NoOfCols];
            Array.Copy(helper, jCol, cols_target, 0, source.NoOfCols);

            source.AccSubMatrixTo(1, target, rows_source, rows_target, cols_source, cols_target);
        }


    }

#endregion
#region Helper Class for Agglomeration
    class MiniMapping {

        /// <summary>
        /// DOFs per cell, per variable, per species.
        /// index: variable.
        /// </summary>
        int[] NS;

        bool[] VarIsXdg;

        UnsetteledCoordinateMapping m_Map;
        SpeciesId m_spId;

        LevelSetTracker.LevelSetRegions m_LsRegion;

        public int MaxDeg = -1;

        public int NoOfVars;

        public MiniMapping(UnsetteledCoordinateMapping Map, SpeciesId spId, LevelSetTracker.LevelSetRegions r) {
            m_Map = Map;
            m_spId = spId;
            m_LsRegion = r;

            Basis[] BS = Map.BasisS.ToArray();
            NS = new int[BS.Length];
            VarIsXdg = new bool[BS.Length];
            NoOfVars = BS.Length;

            for(int iVar = 0; iVar < BS.Length; iVar++) {
                XDGBasis xBasis = BS[iVar] as XDGBasis;
                if(xBasis != null) {
                    NS[iVar] = xBasis.NonX_Basis.Length;
                    //m_LsTrk = xBasis.Tracker;
                    VarIsXdg[iVar] = true;
                } else {
                    NS[iVar] = BS[iVar].Length;
                    VarIsXdg[iVar] = false;
                }

                MaxDeg = Math.Max(MaxDeg, BS[iVar].Degree);
            }
        }

        public long i0Func(int jCell, int iVar) {
            if(VarIsXdg[iVar]) {
                int iSpc = m_LsRegion.GetSpeciesIndex(this.m_spId, jCell);
                return m_Map.GlobalUniqueCoordinateIndex(iVar, jCell, iSpc * NS[iVar]);
            } else {
                return m_Map.GlobalUniqueCoordinateIndex(iVar, jCell, 0);
            }
        }

        public int NFunc(int jCell, int iVar) {
            return NS[iVar];
        }


    }
#endregion
}

