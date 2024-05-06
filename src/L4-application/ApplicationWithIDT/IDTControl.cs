﻿
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;


namespace ApplicationWithIDT {
    /// <summary>
    /// Defines the base Control object class used in the ApplicationWithIDT solvers.
    /// </summary>
    public class IDTControl : CompressibleControl {

        /// <summary>
        /// Creates Default IDTControl
        /// </summary>
        public IDTControl() {
            InitialValueFunctionsPerSpecies = new Dictionary<string, Func<double[], double>>();
            InitialValueFunctionsPerSpecies.Add("L", x => 0);
            InitialValueFunctionsPerSpecies.Add("R", x => 0);
            quadOrderFunc = (int[] A, int[] B, int[] C) => Math.Abs(A.Max()) + Math.Abs(C.Max()) + Math.Max(this.LevelSetDegree, this.LevelSetTwoDegree);
            FixedSQPIterations = new int[] { 20, 20, 20, 20 };
            this.NonLinearSolver = new NonLinearSolverConfig();
            this.NonLinearSolver.MaxSolverIterations = 10;
            this.NonLinearSolver.ConvergenceCriterion = 1e-08;
            this.NonLinearSolver.verbose = false;
            this.SpeciesTable[0, 0] = "L";
            this.SpeciesTable[0, 1] = "R";
            this.SpeciesTable[1, 0] = "A";
            this.SpeciesTable[1, 1] = "B";
            base.NoOfMultigridLevels = 1;
        }

        //[ExclusiveLowerBound(0.0)]
        /// <summary>
        /// Threshold for agglomeration
        /// </summary>
        public double AgglomerationThreshold { get; set; } = 0.0;
        /// <summary>
        /// The configured polynomial DG degree of the level set
        /// </summary>
        public int LevelSetDegree { get; set; } = 1;

        /// <summary>
        /// here one can add a value that will be used to normalized the <field,spc> combination
        /// </summary>
        public Dictionary<Tuple<string, string>, double> BaseFlowPerSpeciesAndField { get; set; } = null;

        /// <summary>
        /// Function determining the QuadOrderdegree
        /// </summary>
        public Func<int[], int[], int[], int> quadOrderFunc { get; set; }


        /// <summary>
        /// path to mesh which has to be loaded
        /// </summary>
        public string MeshPath { get; set; } = null;

        /// <summary>
        /// The Dg Degree of the Solution
        /// </summary>
        public int SolDegree { get; set; }
        /// <summary>
        /// indicates if the initialLevelSet is far from the exact position. So every hard-coded control can have to options 
        /// </summary>
        public bool isFarConfig { get; set; } = true;
        /// <summary>
        /// The quadrature degree used for the bulk and level set terms
        /// </summary>
        public int NonlinearQuadratureDegree { get; set; } = 4;
        public bool ComputeOnlyResidual { get; set; } = false;

        public double flux_s_alpha { get; internal set; } = 10;
        [DataMember]
        //public LinearSolverConfig LinearSolver = new LinearSolverConfig();
        public BoSSS.Solution.AdvancedSolvers.ISolverFactory LinearSolver { get; set; } = new BoSSS.Solution.AdvancedSolvers.DirectSolver.Config();
        /// <summary>
        /// Configuration of 'primary' nonlinear solver, used for the initial value computation
        /// </summary>
        [DataMember]
        public NonLinearSolverConfig NonLinearSolver { get; set; }
        public MassMatrixShapeandDependence MassMatrixShapeandDependence { get; set; } = MassMatrixShapeandDependence.IsNonIdentity;


        #region Initial Value/Boundary Conditions
        public Dictionary<string, Func<double[], double>> InitialValueFunctionsPerSpecies { get; set; }
        public Func<double[], double> GetInitialValueFunc(string species) {
            Func<double[], double> ret;
            this.InitialValueFunctionsPerSpecies.TryGetValue(species, out ret);
            return ret;
        }
        // only used for Burgers and Scalar Advection, defines the Dirichlet Boundaries
        public Func<double[], double> DirichletBoundaryMap { get; set; }
        /// <summary>
        /// function defining EdgeTags for the Mesh
        /// </summary>
        public Func<double[], byte> EdgTagFunc { get; set; } = null;

        /// <summary>
        /// EdgeTagNames array
        /// </summary>
        public string[] EdgTagNames { get; set; } = null;
        public string SeedFromAV_Db { get; set; } = null;
        public Tuple<Guid, TimestepNumber> SeedFromAV_Db_Info { get; set; } = null;
        #endregion

        #region Level Set Stuff
        public bool IsTwoLevelSetRun { get; set; } = true;
        public Func<double[], double> LevelSetTwoInitialValue { get; set; } = x => 0.5 - x[0];
        public string[] SpeciesToEvaluate { get; set; } = null;
        public string[,] SpeciesTable { get; set; } = new string[2, 2];
        public string LsOne_NegSpecies { get; set; } = "V";
        public string LsOne_PosSpecies { get; set; } = "L";
        public string LsTwo_NegSpecies { get; set; } = "L";
        public string LsTwo_PosSpecies { get; set; } = "R";
        public string[,] LsOne_SpeciesPairs { get; set; } = new string[,] { { "L", "R" } };
        public string[,] LsTwo_SpeciesPairs { get; set; } = new string[,] { { "L", "R" } };
        public bool WriteLSCoordinates { get; set; } = false;
        public Func<IGrid> LevelSetGridFunc { get; set; }
        public bool OptiLSIsOrthonormal { get; set; } = false;
        public int OptiLevelSetDegree { get; set; } = 1;
        /// <summary>
        /// enum to choose where the Grid comes from
        /// </summary>
        public GetGridFrom getGridFrom { get; set; } = GetGridFrom.GridFunc;
        /// <summary>
        /// these are the Parameter Names of the OptiLevelSet
        /// </summary>
        public List<string> OptiLevelSet_ParamNames { get; set; }
        /// <summary>
        /// these are the initial ParamValues a_i the OptiLevelSet will consist of
        /// </summary>
        public List<double> OptiLevelSet_ParamValues { get; set; }
        /// <summary>
        /// These are the basis function the OptiLevelSet will consist of
        /// </summary>
        public List<Func<double[], double, double>> OptiLevelSet_Param_Functions { get; set; }
        /// <summary>
        /// LevelSetCFL is thrown if LevelSet moves more than this amount of cells (set this value higher if you want to allow more movement)
        /// </summary>
        public int NearRegionWidth { get; set; } = 1;
        /// <summary>
        /// The configured polynomial DG degree of the second level set
        /// </summary>
        public int LevelSetTwoDegree { get; set; } = 1;
        /// <summary>
        /// Position of LevelSetOne, only needed when two LevelSets are used
        /// </summary>
        public Func<double[], double> LevelSetOneInitialValue { get; set; } = null;
        public string ShockLevelSet_Db { get; set; } = null;
        public Tuple<Guid, TimestepNumber> ShockLevelSet_Info { get; set; } = null;
        public string ShockLevelSet_FieldName { get; set; } = null;
        public bool ShockLevelSet_UseShockGrid { get; set; } = false;
        public bool ShockLevelSet_SeedFromDb { get; set; } = false;
        public string ShockLevelSet_RestartDb { get; set; } = null;
        public Tuple<Guid, TimestepNumber> ShockLevelSet_RestartInfo { get; set; } = null;
        #endregion

        /// <summary>
        /// Adaptive regularization parameters
        /// </summary>
        public int L { get; set; } = 1;
        public double sigma_1 { get; set; } = 1e-2;
        public double sigma_2 { get; set; } = 1e-1;
        public double tauGamma { get; set; } = 1.5;
        internal double tauAlpha { get; set; } = 0.5;
        public double Gamma_Max { get; set; } = 1;
        public double Gamma_Min { get; set; } = 1e-4;
        public double Gamma_Start { get; set; } = 1;

        /// <summary>
        /// level set quality parameters
        /// </summary>
        public double Kappa_Xi { get; set; } = 1;
        public double Kappa_M { get; set; } = 100;
        public double Kappa_Min { get; set; } = 1e-10;
        public double Kappa_v { get; set; } = 0.5;

        /// <summary>
        /// Globalization parameters
        /// </summary>
        public double Alpha_Min { get; set; } = 1e-8;
        public double Alpha_Start { get; set; } = 1;
        public double Mu_Start { get; set; } = 1;
        public double Mu_Max { get; set; } = 1e6;
        public double Mu_Omega { get; internal set; } = 1.2;
        public double Mu_Rho { get; internal set; } = 0.95;

        /// <summary>
        /// Reinitialization parameters
        /// </summary>
        public bool ApplyReiInit { get; set; } = true;
        public double reInit_c1 { get; set; } = -2e-1; // c_6 close to eq. (65) taken from,  A robust, high-order implicit shock tracking method for simulation of complex, high-speed flows
        public double reInit_c2 { get; set; } = 1e-2; // c_7 close to eq. (67) taken from,  A robust, high-order implicit shock tracking method for simulation of complex, high-speed flows
        public double reInit_c3 { get; set; } = 1e-2; // c_8 close to eq. (70) taken from,  A robust, high-order implicit shock tracking method for simulation of complex, high-speed flows
        public double Kappa { get; set; } = 0.5;

        public double[] reInitTols { get; set; } = new double[] { -2e-1, -2e-1, -2e-1, -2e-1, -2e-1, 0 };

        /// <summary>
        /// Staggered solver
        /// </summary>
        public SolverRunType solRunType { get; set; } = SolverRunType.Standard;
        //public StaggerdRunConfig staggerdRunConfig { get; set; } = StaggerdRunConfig.ConstantTimesteps;
        public int DgDegree_Start { get; set; } = 0;

        public TerminationStrategy terStrat { get; set; } = TerminationStrategy.Skyline;
        /// <summary>
        /// linearization used for the State DOF Jacobian dR1dU,dR0dU
        /// </summary>
        public Linearization Linearization { get; set; }
        /// <summary>
        /// FDs used for differentiation of LevelSet DOFS
        /// </summary>
        public FDtype FDtype { get; set; } = FDtype.Central;

        /// <summary>
        /// type of IV initialization
        /// </summary>
        public GetInitialValue GetInitialValue { get; set; } = GetInitialValue.FromFunctionPerSpecies;

        /// <summary>
        /// type of LevelSet which will be optimized
        /// </summary>
        public OptiLevelSetType OptiLevelSetType { get; set; } = OptiLevelSetType.GlobalLevelSet;
        /// <summary>
        /// type of LevelSet initialization
        /// </summary>
        public GetLevelSet GetLevelSet { get; set; } = GetLevelSet.FromParams;

        /// <summary>
        /// globalization method used
        /// </summary>
        public GlobalizationStrategy GlobalizationStrategy { get; set; } = GlobalizationStrategy.LineSearch;

        /// <summary>
        /// Merit Function used in GLobalization algorithm
        /// </summary>
        public MeritFunctionType MeritFunctionType { get; set; } = MeritFunctionType.L1Merit;
        

        /// <summary>
        /// enum to choose optimization problem which is solved
        /// </summary>
        public OptProblemType optProblemType { get; set; } = OptProblemType.FullEnRes;
        /// <summary>
        /// enum to choose f_phi 
        /// </summary>
        public FphiType fphiType { get; set; } = FphiType.None;

        /// <summary>
        /// params needed if Initial guess is obtained via implicit timestepping 
        /// </summary>
        public double IG_dt_Start { get; set; } = 1e-03;
        public double IG_beta { get; set; } = 0.2;
        public double IG_nu_Min { get; set; } =0.05;
        public double IG_nu_Max { get; set; } = 0.1;
        public double IG_alpha { get; set; } = 2;
        public double ConvCrit { get; set; } = 0;

        /// <summary>
        /// enum to control different setups for reinitialization 
        /// </summary>
        public ReInitMode reInitMode { get; set; } = ReInitMode.OneTolForEachP;
        /// <summary>
        /// Prescribes fixed amount of iterations for every p-continuation level
        /// </summary>
        public int[] FixedSQPIterations { get; set; } = null;
        /// <summary>
        /// Prescribes minimum amount of iteration for every p-contiuation level
        /// </summary>
        public int[] minimalSQPIterations { get; set; } = new int[] { 30, 30, 10, 10, 10 ,10,10,10, 10 };
        /// <summary>
        /// Maximum number of total SQP Iterations
        /// </summary>
        public int MaxIterations { get; set; } = 200;
        /// <summary>
        /// maximum number of iterations for each p-continuation level which allow for reinitialization
        /// </summary>
        public int[] ReiniTMaxIters { get; set; } = new int[] { 30, 20, 5, 5, 5 ,5,5,5,5};

        /// <summary>
        /// params needed for termination
        /// </summary>
        public int[] TerminationMinNs { get; set; } = new int[] { 8, 8, 8, 8, 8 ,8,8,8,8};
        public double[] tALNRs { get; set; } = new double[] { 1.005, 1.005, 1.001, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01 };
        public bool PartiallyFixLevelSetForSpaceTime { get; set; } = false;
        public bool SaveMatrices { get; set; } = false;
    }
        /// <summary>
        /// this enum controls were we get the initial value from
        /// </summary>
        public enum GetInitialValue {
        FromFunctionPerSpecies,
        FromDBSinglePhase,
        FromDBXDG,
        FromDB_Partial_SeedInflowExactly,
        FromAVRun,
        FromP0Timestepping,
        OneFullNewtonStep,
    }
    /// <summary>
    /// this enum controls the Globalization Strategy
    /// </summary>
    public enum GlobalizationStrategy {
        None = 0,
        LineSearch,
        DogLeg,
        ZickZackSearch
    }


    public enum TerminationStrategy {
        Skyline =0,
        SkylineWithDifferntTermNs,
        MaxVsMin,
        Experimental,
    }
    /// <summary>
    /// this enum controls the FDtype we use for the LevelSetParams
    /// </summary>
    public enum FDtype {
        None = 0,
        Right,
        Left,
        Central
    }
    /// <summary>
    /// this enum controls which linearization we use
    /// </summary>
    public enum Linearization {
        FD,
        JacobiOperator,
        Adhoc
    }

    /// <summary>
    /// FromReconstruction: LevelSet is loaded from DB (from a previously run Reconstruction)
    /// FromOldSimulation: LevelSet is loaded from DB from an old simulation
    /// FromParams: User specifies basis Functions in the ControlFile
    /// FromFunction: User specifies a Func<double[],double> object witch is projected onto an ONB
    /// </summary>
    public enum GetLevelSet {
        FromReconstruction,
        FromOldSimulation,
        FromParams,
        FromFunction,
        FromReconstructionFromPoints,
        DirectyFromTimestep
    }
    /// <summary>
    /// GlobalLevelSet: LevelSet is defined by a set of basis functions and coordinates on the global domain
    /// SinglePhaseField: the BoSSS SinglePhaseField is used to obtain the OptiLevelSet
    /// </summary>
    public enum OptiLevelSetType {
        GlobalLevelSet,
        SinglePhaseField,
        SpecFemField,
        SplineLevelSet
    }
    /// <summary>
    /// the enum controls the Merit Function the globalization uses
    /// </summary>
    public enum MeritFunctionType {
        L1Merit,
        L2Merit
    }
    /// <summary>
    /// This enum controls where we get the grid from
    /// </summary>
    public enum GetGridFrom {
        GridFunc,
        DB
    }
    /// <summary>
    /// this enum tells us what kind of SolverRun we have
    /// </summary>
    public enum SolverRunType {
        Standard,
        PContinuation,
        Unsteady
    }
    /// <summary>
    /// For the staggered solver, controls when the Basis degree is changed
    /// </summary>
    //public enum StaggerdRunConfig {
    //    ConstantTimesteps,
    //    StagnatingLineSearch
    //}
    public enum OptProblemType {
        FullEnRes,
        EnResOnlyCutCells,
        EnResOnlyNearBand,
        RankineHugoniotFull,
        RankineHugoniotOnlyInterface
    }

    public enum FphiType {
        None,
        CurvatureAll,
        CurvatureCut,
        PerssonSensorAll,
        PerssonSensorCut
    }

    public enum ReInitMode {
        OneTolForAllP,
        OneTolForEachP,
    }
}