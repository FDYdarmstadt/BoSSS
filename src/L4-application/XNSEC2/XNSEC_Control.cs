using BoSSS.Application.XNSFE_Solver;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;
using static BoSSS.Solution.NSECommon.MaterialLaw_MultipleSpecies;
using static BoSSS.Solution.NSECommon.SIPDiffusionTemperature;

namespace BoSSS.Application.XNSEC {

    [DataContract]
    [Serializable]
    public class XNSEC_MF_Control : XNSEC_Control {

        public override Type GetSolverType() {
            return typeof(XNSEC_MixtureFraction);
        }


        /// <summary>
        /// Sets the DG polynomial degree
        /// </summary>
        /// <param name="DGp">Degree for velocity; pressure  will be one order lower.</param>
        public override void SetDGdegree(int DGp) {
            if (DGp < 1)
                throw new ArgumentOutOfRangeException("DG polynomial degree must be at least 1.");

            base.FieldOptions.Clear();
            //base.SetDGdegree(DGp);


            FieldOptions.Add("Velocity*", new FieldOpts() {
                Degree = DGp,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.Pressure, new FieldOpts() {
                Degree = DGp - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            FieldOptions.Add(VariableNames.LevelSetDG, new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetCG, new FieldOpts() {
                Degree = Math.Max(2, DGp),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetDGidx(1), new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetCGidx(1), new FieldOpts() {
                Degree = Math.Max(2, DGp),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            //FieldOptions.Add(VariableNames.Curvature, new FieldOpts() {
            //    Degree = Math.Max(2, p) * 2,
            //    SaveToDB = SaveCurvature
            //});

            FieldOptions.Add(VariableNames.MixtureFraction, new FieldOpts() { Degree = DGp, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.Rho, new FieldOpts() { Degree = DGp, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.cp, new FieldOpts() { Degree = DGp, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        }

    }

    /// <summary>
    /// Base control file
    /// </summary>
    [DataContract]
    [Serializable]
    public class XNSEC_Control : XNSFE_Control {

        public XNSEC_Control() {
            base.NoOfMultigridLevels = 1;
            base.SuperSampling = 2;
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            base.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
        }

        public XNSEC_Control(int dgDeg, double _pRef, double _uRef, double _TRef, bool useAdimensional, double _LRef = 0.0, double[] FuelYs = null, double[] OxidizerYs = null, double epsilon = 5, bool analyticalSolutionOK = false) {
            this.CC = new ChemicalConstants();
            this.AnalyticsolutionSwitch = analyticalSolutionOK;

            this.YFuelInlet = FuelYs[0];
            this.YOxInlet = OxidizerYs[1];
            this.TFuelInlet = 1.0;
            this.TOxInlet = 1.0;
            this.NumberOfChemicalSpecies = 4;
            SetDGdegree(dgDeg);
            this.physicsMode = PhysicsMode.Combustion;
            base.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            this.GravityDirection = new double[] { 0.0, 0.0, 0.0 }; //No gravity.
            base.NoOfMultigridLevels = 1;
            base.NonLinearSolver.ConvergenceCriterion = 1e-8;
            base.NonLinearSolver.verbose = true;
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            base.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            base.NonLinearSolver.MaxSolverIterations = 20;
            //base.NoOfTimesteps = int.MaxValue;
            this.BDFOrder = 1;
            base.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            base.saveperiod = 1;
            this.PenaltyViscMomentum = 1.0;
            this.PenaltyHeatConduction = 1.0;
            this.s = (CC.nu_O2 * CC.MW_O2) / (CC.nu_CH4 * CC.MW_CH4);
            this.phi = s * YFuelInlet / YOxInlet;
            this.zSt = 1.0 / (1.0 + phi);

            var MLC = new MaterialLawCombustion(273.11 * 0 + 300, new double[] { }, this.MatParamsMode, this.rhoOne, true, 1.0, 1, 1, YOxInlet, YFuelInlet, this.zSt, CC, 0.75);
            var ThermoProperties = new ThermodynamicalProperties();
            this.pRef = _pRef; // reference pressure
            this.TRef = _TRef;// reference temperature

            double[] MWs = new double[] { CC.MW_CH4, CC.MW_O2, CC.MW_CO2, CC.MW_H2O, CC.MW_N2 };

            this.MWRef = MLC.getAvgMW(MWs, OxidizerYs);
            this.rhoRef = pRef * MWRef / (8.314 * TRef * 1000); // Kg/m3. ok ;
            this.T_ref_Sutherland = 273;
            this.cpRef = 1.4;// ThermoProperties.Calculate_Cp_Mixture(new double[] { 0.23, 0.77 }, new string[] { "O2", "N2" }, 300); //
            Console.WriteLine("Cpref set to: " + this.cpRef);
            this.muRef = MLC.getViscosityDim(273);//revisar!! antes estaba en 300... pero supongo que si sutherland law ocupa 273 debo poner eso tambien en 273 para mantenerme consistente

            this.B = CC.PreExponentialFactor;
            this.DRef = MLC.get_LambdaCp_Term(TRef) / rhoRef; // lambda/cp = rho Di

            this.uRef = _uRef > 0 ? _uRef : Math.Pow(DRef * epsilon, 0.5); // Reference velocity
            this.LRef = _LRef > 0 ? _LRef : DRef / uRef; // reference length

            double Ta; // Activation temperature, K
            double heatRelease; //Mass heat release, KJ/kg

            Ta = CC.Ta;

            heatRelease = this.CC.HeatReleaseMass;

            double heatRelease_Ref = (TRef * cpRef);

            this.HeatRelease = heatRelease / heatRelease_Ref;

            this.MolarMasses = new double[] { this.CC.MW_CH4, this.CC.MW_O2, this.CC.MW_CO2, this.CC.MW_H2O, this.CC.MW_N2 };
            this.MolarMasses.ScaleV(1.0 / MWRef);//NonDimensionalized Molar masses
            this.StoichiometricCoefficients = new double[] { -1, -2, 1, 2, 0 };

            this.NumberOfProducts = 3;
            this.Damk = rhoRef * this.LRef * B / (uRef * MWRef);
            this.Reynolds = rhoRef * uRef * this.LRef / muRef;
            this.Prandtl = 0.75;////muRef * cpRef / lambdaRef; // Air prandtl number
            this.Schmidt = this.Prandtl; // Because Lewis number  is assumed as 1.0  (Le = Pr/Sc)

            //this.Lewis = new double[] { 0.97, 1.11, 1.39, 0.83, 1.0 }; // If i use this, N2 appears in the system?..........................
            this.Lewis = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };

            double g = 9.8; // m/s2
            this.Froude = Math.Sqrt(uRef * uRef / (this.LRef * g));

            double Ta_adim = Ta / (TRef);
            this.ReactionRateConstants = new double[] { this.Damk, Ta_adim, 1.0, 1.0 };

            if (!useAdimensional) { // Set up for a run with dimensional variables
                this.pRef = _pRef; // reference pressure
                this.TRef = 1;// reference temperature
                this.MWRef = MLC.getAvgMW(MWs, OxidizerYs);
                this.rhoRef = pRef * MWRef / (R_gas * TRef * 1000); // Kg/m3. ok ;
                this.cpRef = 1.31;// Representative value, KJ/Kg K
                this.muRef = MLC.getViscosityDim(300);
                this.DRef = 1.0;
                this.uRef = 1.0;
                this.LRef = 1.0;
                this.R_gas = 8.314; //KJ/KmolK should this be multipleid with 1000?????
                this.AmbientPressure = pRef;
                this.MolarMasses = new double[] { this.CC.MW_CH4, this.CC.MW_O2, this.CC.MW_CO2, this.CC.MW_H2O, this.CC.MW_N2 };

                this.StoichiometricCoefficients = new double[] { -1, -2, 1, 2, 0 };
                //this.ViscosityOption = ViscosityOption.VariableViscosity;

                this.Damk = 1.0;
                this.Froude = 1.0;
                this.Reynolds = 1.0 / muRef;
                this.Prandtl = 0.75;
                this.Schmidt = this.Prandtl;

                this.Lewis = new double[] { 0.97, 1.11, 1.39, 0.83, 1.0 };

                //this.GravityDirection = new double[] { 0.0, 9.8, 0.0 };

                this.T_ref_Sutherland = 300;

                this.HeatRelease = this.CC.HeatReleaseMass;
                this.ReactionRateConstants = new double[] { CC.PreExponentialFactor, CC.Ta, 1.0, 1.0 };
            }

            this.AdiabaticTemperature = zSt * TFuelInlet + (1 - zSt) * TOxInlet + this.HeatRelease * this.YFuelInlet * zSt;
            Console.WriteLine("Damköhler number is: {0}", ReactionRateConstants[0]);
            Console.WriteLine("Adimensional activation temperature is: {0}", Ta_adim);
            Console.WriteLine("Reynolds number is {0}", Reynolds);
            Console.WriteLine("Prandtl number is {0}", Prandtl);
            Console.WriteLine("The Maximum temperature reached should be {0}", this.AdiabaticTemperature);
            Console.WriteLine("Cp: {0}", cpRef);
            Console.WriteLine("Stoichiometric mixture fraction (z_st): {0}", this.zSt);
            //Console.WriteLine("Reference time in seconds: {0}", timeRef);
        }

        public override Type GetSolverType() {
            return typeof(XNSEC);
        }

        /// <summary>
        /// Sets the DG polynomial degree
        /// </summary>
        /// <param name="DGp">Degree for velocity; pressure  will be one order lower.</param>
        public override void SetDGdegree(int DGp) {
            if (DGp < 1)
                throw new ArgumentOutOfRangeException("DG polynomial degree must be at least 1.");

            base.FieldOptions.Clear();
            base.SetDGdegree(DGp);

            FieldOptions.Add(VariableNames.ThermodynamicPressure, new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            var bla = this.EnableTemperature ? DGp : 0;
            FieldOptions[VariableNames.Temperature] = new FieldOpts() { Degree = this.EnableTemperature ? DGp : 0, SaveToDB = FieldOpts.SaveToDBOpt.TRUE };
            for (int i = 0; i < this.NumberOfChemicalSpecies; i++) {
                FieldOptions.Add(VariableNames.MassFraction_n(i), new FieldOpts() { Degree = this.EnableMassFractions? DGp : 0, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            }







            FieldOptions.Add(VariableNames.MixtureFraction, new FieldOpts() { Degree = DGp, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.Rho, new FieldOpts() { Degree = DGp, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add("kReact", new FieldOpts() { Degree = DGp, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            FieldOptions.Add(VariableNames.cp, new FieldOpts() { Degree = DGp, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        }

        public void SetSaveOptions(string dataBasePath = null, int savePeriod = 1) {
            if (dataBasePath != null) {
                savetodb = true;
                DbPath = dataBasePath;
                saveperiod = savePeriod;
            } else
                savetodb = false;
        }

        public void SetTimeSteppingOptions(double dt, double endtime) {
            if (dt <= 0) {
                this.TimesteppingMode = _TimesteppingMode.Steady;
            } else {
                this.TimesteppingMode = _TimesteppingMode.Transient;
                this.dtFixed = dt;
                this.Endtime = endtime;
                this.NoOfTimesteps = (int)(endtime / dt);
            }
        }

        public void SetAdaptiveMeshRefinement(int amrLevel, int pseudoTimeStepsNo, int _AMR_startUpSweeps = -1) {
            NoOfTimesteps = pseudoTimeStepsNo;
            if (amrLevel == 0) {
                return;
            }
            AdaptiveMeshRefinement = amrLevel > 1 ? true : false;
            RefinementLevel = amrLevel;
            AMR_startUpSweeps = _AMR_startUpSweeps;

            Console.WriteLine("No of start up sweeps " + AMR_startUpSweeps);
        }

        /// <summary>
        /// DELETE!!
        ///  i added this in order to be able to use  SetEqualityBasedSessionJobControlCorrelation with the HPC cluster
        ///  SetNameBasedSessionJobControlCorrelation doesnt work..
        /// </summary>
        [DataMember]
        public int dummycounter = 1;

        // Solver configuration
        //---------------------
        /// <summary>
        /// Penalty factor for viscous term in the moment equation
        /// </summary>
        [DataMember]
        public double PenaltyViscMomentum = 1.0;

        /// <summary>
        /// Penalty factor for heat conduction term in the energy equation
        /// </summary>
        [DataMember]
        public double PenaltyHeatConduction = 1.0;

        ///// <summary>
        ///// Number of subdivisions of the homotopy algorithm
        ///// </summary>
        //[DataMember]
        //public int NumberOfHomotopyArraySubdivisions = 10;

        ///// <summary>
        ///// Terms activated in the SIP-Viscosity terms
        ///// </summary>
        //[DataMember]
        //public ViscosityTermsSwitch myviscosityTerms = (ViscosityTermsSwitch.grad_u | ViscosityTermsSwitch.grad_uT | ViscosityTermsSwitch.divU);

        ///<summary>
        /// Block-Preconditiond for the velocity/momentum-block of the saddle-point system
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;

        /// <summary>
        /// Block-Preconditiond for the pressure/continuity-block of the saddle-point system
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode PressureBlockPrecondMode = MultigridOperator.Mode.Eye; // no SymPart_Diag-Präcon, because there may be no zero on the diagonal!!!

        /// <summary>
        /// Unity spatial direction vector of gravity.
        /// </summary>
        [DataMember]
        public double[] GravityDirection = { 0.0, 1.0, 0.0 };

        /// <summary>
        ///
        /// </summary>
        [DataMember]
        public double phi;

        ///// <summary>
        ///// Switch for activating the selfmade Homotopie algorithm. The variable <see cref="homotopieVariableName"/>
        ///// is gradually incremented until the desired value <see cref="homotopieAimedValue"/> is reached.
        ///// </summary>
        //[DataMember]
        //public bool useHomotopie = false;

        /// <summary>
        /// Switch for using variable (equivalence ratio dependent) parameters of the one Step combustion model
        /// </summary>
        [DataMember]
        public bool VariableOneStepParameters = false;

        /// <summary>
        /// Switch for using the simplified geometry in the CounterFlow Configuration
        /// </summary>
        [DataMember]
        public bool usesimplifiedGeometry = false;

        /// <summary>
        ///
        /// </summary>
        [DataMember]
        public bool UseMixtureFractionsForCombustionInitialization = false;

        /// <summary>
        /// name of the variable going to be increased in the homotopie algorithm.
        /// </summary>
        [DataMember]
        public string homotopieVariableName;

        /// <summary>
        /// Aimed value of the homotopie variable.
        /// </summary>
        [DataMember]
        public double homotopieAimedValue;

        /// <summary>
        ///PlotNewtonIterations
        /// </summary>
        [DataMember]
        public bool PlotNewtonIterations = false;

        /// <summary>
        ///PlotNewtonIterations
        /// </summary>
        [DataMember]
        public bool PlotAdditionalParameters = true;

        /// <summary>
        /// Sensor Variable
        /// </summary>
        [DataMember]
        public string SensorVariable = VariableNames.ReactionRate;

        /// <summary>
        /// Reference point for pressure.
        /// Only needed if there is no Dirichlet boundary condition for the pressure.
        /// </summary>
        [DataMember]
        public double[] PressureReferencePoint;

        /// <summary>
        /// Desired BDFOrder for the timesteper
        /// </summary>
        [DataMember]
        public int BDFOrder = 1;

        /// <summary>
        /// Mode for the determination of the material parameters
        /// </summary>
        [DataMember]
        public MaterialParamsMode MatParamsMode = MaterialParamsMode.Constant;

        /// <summary>
        /// Switch for mode of heat capacity calculation
        /// </summary>
        [DataMember]
        public CpCalculationMode HeatCapacityMode = CpCalculationMode.constant;

        //}
        /// <summary>
        /// Exact solution for Mass Fractions, for each species (either A or B).
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public IDictionary<string, Func<double[], double, double>[]> ExactSolutionMassFractions;

        [DataMember]
        public double[] FuelInletConcentrations = new double[] { 0.2, 0.0, 0.0, 0.0, 0.8 };

        [DataMember]
        public double[] OxidizerInletConcentrations = new double[] { 0.0, 0.23, 0.0, 0.0, 0.77 };

        [DataMember]
        public double YFuelInlet { get; set; } = 1.0;

        [DataMember]
        public double YOxInlet { get; set; } = 1.0;

        [DataMember]
        public double TFuelInlet { get; set; } = 1.0;

        [DataMember]
        public double TOxInlet { get; set; } = 1.0;

        [DataMember]
        public double zSt { get; set; } = 0.05;

        [DataMember]
        public double s { get; set; } = 4.0;

        [DataMember]
        public double AdiabaticTemperature { get; set; }

        ///// <summary>
        ///// Viscosity, density and surface tension.
        ///// </summary>
        //[DataMember]
        //public PhysicalParameters PhysicalParameters = new PhysicalParametersCombustion() {
        //    Material = true,
        //    IncludeConvection = false,
        //    IncludeDiffusion = true,
        //    mu_A = 1.0,
        //    mu_B = 1.0,
        //    rho_A = 1.0,
        //    rho_B = 1.0,
        //    Sigma = 0.0,
        //    rhoD_A = 1.0,
        //    rhoD_B = 1.0

        //};

        ///// <summary>
        ///// Just a cast of <see cref="PhysicalParameters"/>
        ///// </summary>
        //public PhysicalParametersCombustion PhysicalParametersCombustion {
        //    get {
        //        return (PhysicalParametersCombustion)PhysicalParameters;
        //    }
        //}

        public ChemicalConstants CC {
            get; set;
        } = new ChemicalConstants();

        /// <summary>
        /// Sets rho equal to one, i.e. there is no dependency of mom-i equations and energy. Just for debugging purposes
        /// </summary>
        [DataMember]
        public bool rhoOne = false;

        /// <summary>
        /// Viscosity mode
        /// </summary>
        [DataMember]
        public ViscosityOption ViscosityOption = ViscosityOption.VariableViscosityDimensionless;

        /// <summary>
        /// Viscosity mode
        /// </summary>
        [DataMember]
        public ThermalWallType myThermalWallType = ThermalWallType.fixedTemperature;

        /// <summary>
        /// Enables the manufactured solution source terms in the operator
        /// </summary>
        [DataMember]
        public bool ManufacturedSolutionSwitch = false;

        /// <summary>
        /// Used for activating the calculation of the error based on the predefined solution
        /// </summary>
        [DataMember]
        public bool AnalyticsolutionSwitch = false;

        /// <summary>
        /// EdgeTags for boundaries, where the Nusselt number
        /// should be calculated.
        /// </summary>
        [DataMember]
        public string[] EdgeTagsNusselt = null;

        /// <summary>
        /// Physics mode used in the simulation. Low Mach solves u,v,p2,T and Combustion solves additionally for Y0,Y1, Y2 and implicitly Y3
        /// </summary>
        [DataMember]
        public PhysicsMode m_physicsMode = PhysicsMode.Combustion;

        public PhysicsMode physicsMode {
            get {
                return m_physicsMode;
            }
            set {
                m_physicsMode = value;
            }
        }

        /// <summary>
        ///
        /// </summary>
        [DataMember]
        public double[][] troubledPoints = null;

        /// <summary>
        /// desired minimum refinement level at interface
        /// </summary>
        [DataMember]
        public int BaseRefinementLevel = 1;

        /// <summary>
        /// maximum refinement level including additional refinement (contact line, curvature, etc.)
        /// </summary>
        [DataMember]
        public int RefinementLevel = 1;

        /// <summary>
        /// If true the Jacobi matrix will be calculated by Finite differences
        /// </summary>
        [DataMember]
        public bool UseFDJ = false;

        /// <summary>
        /// if true, the base vectorization will be used. Mainly used for performance measurements
        /// </summary>
        [DataMember]
        public bool IgnoreVectorizedFluxes = false;

        /// <summary>
        /// Mode for thermodynamic pressure.
        /// </summary>
        [DataMember]
        public ThermodynamicPressureMode ThermodynamicPressureMode = ThermodynamicPressureMode.Constant;

        /// <summary>
        /// Can to be prescribed by the user for <see cref="ThermodynamicPressureMode.MassDetermined"/>.
        /// The mass of the system is assumed by default to by 1.0
        /// </summary>
        [DataMember]
        public double InitialMass = 1.0;

        /// <summary>
        /// Value to initialize thermodynamic pressure
        /// </summary>
        [DataMember]
        public double AmbientPressure = 1.0;

        /// <summary>
        /// Flag for activating/deactivating the chemical reactions
        /// </summary>
        [DataMember]
        public bool ChemicalReactionActive = true;

        /// <summary>
        /// Flag for d rho/ dt term in the Continuity equation.
        /// Only for debugging purposes
        /// </summary>
        [DataMember]
        public bool timeDerivativeConti_OK = false;

        /// <summary>
        /// Flag for d p0/ dt term in the energy equation.
        /// Only for debugging purposes
        /// </summary>
        [DataMember]
        public bool timeDerivativeEnergyp0_OK = false;

        /// <summary>
        /// Enable mass fractions equations
        /// Only for debugging purposes
        /// </summary>
        [DataMember]
        public bool EnableMassFractions = true;

        /// <summary>
        /// Enable temperature equation
        /// Only for debugging purposes
        /// </summary>
        [DataMember]
        public bool EnableTemperature = true;

        /// <summary>
        /// Enable self made temporal operator of the LowMach solver
        /// Only for debugging purposes
        /// </summary>
        [DataMember]
        public bool UseSelfMadeTemporalOperator = true;

        /// <summary>
        /// Minimal and maximal values allowed for each variable
        /// </summary>
        //[NonSerialized]
        //[JsonIgnore]
        [DataMember]
        public Dictionary<string, Tuple<double, double>> VariableBounds = null;

        /// <summary>
        /// Homotopy values
        /// Allows gradual increase of the Reynolds number to improve convergence
        /// </summary>
        [DataMember]
        public double[] SelfDefinedHomotopyArray;


 

        public double[] HomotopyArray {
            get {
                //if (m_HomotopyArray == null && this.HomotopyApproach == HomotopyType.Manual) {
                //    SelfDefinedHomotopyArray
                //}
                return SelfDefinedHomotopyArray;
            }
            set {
                SelfDefinedHomotopyArray = value;
                this.NoOfTimesteps = SelfDefinedHomotopyArray.Length;
            }
        }


        [DataMember]
        public string NameFieldForSensor;

        /// <summary>
        /// Homotopy approach (None/manual/automatic)
        /// </summary>
        public enum HomotopyType {
            None = 0,
            Manual = 1,
            Automatic = 2
        }

        [DataMember]
        public HomotopyType HomotopyApproach = HomotopyType.None;

        /// <summary>
        /// Desired homotopy variable
        /// </summary>
        public enum HomotopyVariableEnum {
            None = 0,
            Reynolds = 1,
            VelocityInletMultiplier = 2,
            HeatOfCombustion = 3,
        }

        [DataMember]
        public HomotopyVariableEnum m_HomotopyVariable;

        public HomotopyVariableEnum HomotopyVariable {
            get {
                return m_HomotopyVariable;
            }
            set {
                if (value == HomotopyVariableEnum.Reynolds)
                    homotopieVariableName = HomotopieVariableNames.Reynolds;
                else if (value == HomotopyVariableEnum.VelocityInletMultiplier)
                    homotopieVariableName = HomotopieVariableNames.VelocityMultiplier;
                else if (value == HomotopyVariableEnum.HeatOfCombustion)
                    homotopieVariableName = HomotopieVariableNames.HeatOfReaction;
                else
                    throw new Exception("Wrong homotopyvariable");

                m_HomotopyVariable = value;
            }
        }

        // Physical Parameters
        //--------------------

        /// <summary>
        /// Reference temperature for sutherlands law for the viscosity in the GetViscosity function of the MaterialLawLowMach class.
        /// </summary>
        [DataMember]
        public double T_ref_Sutherland = 600;

        /// <summary>
        /// Reynolds number
        /// </summary>
        [DataMember]
        public double Reynolds = 1.0;

        /// <summary>
        /// Starting value used within the homotopy algorithm.
        /// Should be a value which makes the problem "easy" to solve.
        /// </summary>
        [DataMember]
        public double StartingHomotopyValue = 1.0;

        /// <summary>
        /// Rayleigh number. For the case of cavity natural convection Ra = Re*Re
        /// </summary>
        [DataMember]
        public double Rayleigh = 1.0;

        /// <summary>
        /// Heat capacity ratio
        /// </summary>
        [DataMember]
        [ExclusiveLowerBound(0.0)]
        public double HeatCapacityRatio = 1.0;

        /// <summary>
        /// Prandtl number
        /// </summary>
        [DataMember]
        public double Prandtl = 1.0;

        /// <summary>
        /// Schmidt number
        /// </summary>
        [DataMember]
        public double Schmidt = 1.0;

        /// <summary>
        /// Array if Lewis number for each species
        /// </summary>
        [DataMember]
        public double[] Lewis = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };

        /// <summary>
        /// Froude number.
        /// </summary>
        [DataMember]
        [InclusiveLowerBound(0.0)]
        public double Froude = 1.0;

        /// <summary>
        /// Dahmköhler number.
        /// </summary>
        [DataMember]
        [InclusiveLowerBound(0.0)]
        public double Damk = 1.0;

        //Combustion Parameters
        //---------------------

        /// <summary>
        /// The number of reaction products. Usually two (i.e. CO2 and H2O). Also inerts are considered here (N2 for example)
        /// </summary>
        [DataMember]
        public int NumberOfProducts { get; set; } = 3;

        /// <summary>
        /// The number of reaction products. Usually two (i.e. CO2 and H2O). Also inerts are considered here (N2 for example)
        /// </summary>
        [DataMember]
        public int NumberOfChemicalSpecies { get; set; } = 1;

        /// <summary>
        /// 1. 0 molar mass, 2. 1 molar mass, 3. molar mass 2, molar mass 3
        /// </summary>
        [DataMember]
        public double[] MolarMasses = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };

        /// <summary>
        /// 1. PreExpFactor/Damköhler number, 2. ActivationTemperature, 3. MassFraction0Exponent, 4. MassFraction1Exponent
        /// </summary>
        [DataMember]
        public double[] ReactionRateConstants = new double[] { 100.0, 3.0, 1.0, 1.0 };

        /// <summary>
        /// 1. 0 (fuel) coefficient, 2. 1 (oxidizer) coefficient, 3. 2 (CO2) coefficient, 3 (H2O) coefficient
        /// </summary>
        [DataMember]
        public double[] StoichiometricCoefficients = new double[] { -1, -2, 1, 2, 0 };

        /// <summary>
        /// Reference for the MolarMasses (i.e. value of M_\infty)
        /// </summary>
        [DataMember]
        public double M_ref = 1.0;

        /// <summary>
        /// Partial heat capacities of the species 1 to ns. Measured relative to the mean heat capacity.
        /// No good solution for non-constant heat capacity. Needs to be changed for non-constant heat capacities.
        /// </summary>
        [DataMember]
        public double[] PartialHeatCapacities = new double[] { 1.5, 1.0, 1.0, 1.0, 1.0 };

        /// <summary>
        /// Stoichiometric coefficient S
        /// </summary>
        [DataMember]
        public double S_coef = 4.0;

        /// <summary>
        /// ID of session to load as initial conditions.
        /// </summary>
        [DataMember]
        public string ID;

        /// <summary>
        /// Analytic solution temperature.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticTemperature = null;

        /// <summary>
        /// Analytic solution VelocityX.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticVelocityX = null;

        /// <summary>
        /// Analytic solution VelocityY.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticVelocityY = null;

        /// <summary>
        /// Analytic solution VelocityZ.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticVelocityZ = null;

        /// <summary>
        /// Analytic solution Pressure.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticPressure = null;

        /// <summary>
        /// Analytic solution Y0.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticY0 = null;

        /// <summary>
        /// Analytic solution Y1.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticY1 = null;

        /// <summary>
        /// Analytic solution Y2.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticY2 = null;

        /// <summary>
        /// Analytic solution Y3.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticY3 = null;

        /// <summary>
        /// Analytic solution Y4.
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double> AnalyticY4 = null;

        /// <summary>
        /// Manufactured solution Momentum X
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_MomentumX = null;

        /// <summary>
        /// Manufactured solution Momentum X
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_MomentumY = null;

        /// <summary>
        /// Manufactured solution Momentum X
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_MomentumZ = null;

        /// <summary>
        /// Manufactured solution continuity
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_Continuity = null;

        /// <summary>
        /// Manufactured solution Energy
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_Energy = null;

        /// <summary>
        /// Manufactured solution Species0
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_Species0 = null;

        /// <summary>
        /// Manufactured solution Species1
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_Species1 = null;

        /// <summary>
        /// Manufactured solution Species2
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_Species2 = null;

        /// <summary>
        /// Manufactured solution Species3
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> ManufacturedSolution_Species3 = null;

        /// <summary>
        /// Number of time steps stored in ScalarFieldHistory
        /// </summary>
        [DataMember]
        [InclusiveLowerBound(1)]
        [InclusiveUpperBound(4)]
        public int TimeOrder = 2;

        /// <summary>
        /// Enum that defines a variable, for which the pearson sensor will be calculated.
        /// </summary>
        public enum PearsonSensorVariable {
            None,

            Pressure,

            PressureGradient,

            Temperature,

            TemperatureGradient,

            VelocityX,

            VelocityXGradient,

            MassFraction2,

            kReact
        }

        /// <summary>
        /// Option for calculating thermodynamic pressure in Low-Mach approximation, i.e. p0.
        /// </summary>
        public enum MyThermodynamicPressureMode {

            /// <summary>
            /// Constant in time.
            /// This option is usually used for open systems,
            /// where p0 is set to the ambient pressure.
            /// </summary>
            Constant,

            /// <summary>
            /// This option can be used for closed systems,
            /// where p0 is determined to conserve mass, i.e.
            /// \f$ p_0 = \frac{m(t=0)}{\int 1/T dV} \f$ .
            /// </summary>
            MassDetermined
        }

        /// <summary>
        ///
        /// </summary>
        public static class HomotopieVariableNames {

            /// <summary>
            ///
            /// </summary>
            public const string Reynolds = "Reynolds";

            /// <summary>
            ///
            /// </summary>
            public const string Damkoehler = "Damkoehler";

            /// <summary>
            ///
            /// </summary>
            public const string VelocityMultiplier = "VelocityMultiplier";

            /// <summary>
            ///
            /// </summary>
            public const string HeatOfReaction = "HeatOfReaction";
        }

        /// <summary>
        /// Modus for thermodynamic pressure.
        /// </summary>
        [DataMember]
        public MyThermodynamicPressureMode _MyThermodynamicPressureMode;

        /// <summary>
        /// Treatment of the level-set.
        /// </summary>
        public enum ReactionModel {

            /// <summary>
            /// No Reaction
            /// </summary>
            None = 0,

            /// <summary>
            /// One step methane combustion CH4 + 2 O2 -> CO2 + 2 H2O
            /// </summary>
            OneStepMethane = 1,
        }

        /// <summary>
        /// Modus for thermodynamic pressure.
        /// </summary>
        [DataMember]
        public PearsonSensorVariable myPearsonSensorVariable;

        /// <summary>
        /// Use Persson Sensor to detect high energy modes of singularities
        /// </summary>
        [DataMember]
        public bool UsePearssonSensor = false;

        /// <summary>
        /// bound for perssonsensor
        /// </summary>
        [DataMember]
        public double SensorLimit = 1e-7;

        [DataMember]
        public double LRef { get; set; }

        [DataMember]
        public double pRef { get; set; }

        [DataMember]
        public double uRef { get; set; }

        [DataMember]
        public double MWRef { get; set; }

        [DataMember]
        public double muRef { get; set; } //Kg/(m.s)

        [DataMember]
        public double B { get; set; }

        [DataMember]
        public double DRef { get; set; }

        [DataMember]
        public double cpRef { get; set; } = 1.0;

        [DataMember]
        public double lambdaRef { get; set; }

        [DataMember]
        public double TRef { get; set; } = 300;

        [DataMember]
        public double rhoRef { get; set; }

        [DataMember]
        public double R_gas { get; set; } = 1.0;

        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> SourceTermScalar;

        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> SourceTermConti;

        [NonSerialized]
        [JsonIgnore]
        public Func<double[], double, double> SourceTermMomentum;

        [DataMember]
        private double m_HeatRelease = -1;

        /// <summary>
        /// Factor for smoothing the flame sheet
        /// </summary>
        [DataMember]
        public double smoothingFactor = 10;

        /// <summary>
        /// HeatReleaseFactor per kilogram of fuel (The sum M_alpha cp_alpha v_alpha for alpha = 1 to ns accounting for the enthalphy of reaction).
        /// Note that for rich flames the heat release is smaller than the enthalpy of reaction.
        /// See "A simple one-step chemistry model for partially premixed hydrocarbon combustion" (2006) by Eduardo Fernandez-Tarrazo et. al..
        [DataMember]
        public double HeatRelease {
            get {
                if (m_HeatRelease == -1 && physicsMode == PhysicsMode.Combustion) {
                    Console.WriteLine("Warning!! Heat release should be set by the user for combustion applications. Setting it automatically to zero");
                    m_HeatRelease = 0.0;
                }
                return m_HeatRelease;
            }
            set {
                m_HeatRelease = value;
            }
        }
    }
}