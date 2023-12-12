using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSFE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using CommandLine;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.XNSEC {
    /// <summary>
    /// Low Mach number flow solver. Supports temperature dependeant density and transport parameters.
    /// The mixture fraction solver <see cref="XNSEC_MixtureFraction"/> can be used for finding estimates for combustion applications
    /// </summary>
    public partial class XNSEC : SolverWithLevelSetUpdater<XNSEC_Control> {

        //===========
        // Main file
        //===========
        private static void Main(string[] args) {
            //-n 4 ./XNSEC.exe -c "cs:BoSSS.Application.XNSEC.FullNSEControlExamples.BackwardFacingStep()"

            //InitMPI();
            //BoSSS.Application.XNSEC.NUnitTest.ViscosityJumpTest(2, 1, 0.0d, ViscosityMode.FullySymmetric, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local);
            //BoSSS.Application.XNSEC.NUnitTest.XDG_PSEUDO1D_EVAPORATION_TEST();

            //DeleteOldPlotFiles();
            //DeleteOldTextFiles();
            //Debugger.Launch();

            //NUnitTest.XDG_PSEUDO1D_EVAPORATION_TEST();
            //NUnitTest.CavityNaturalConvection();
            //NUnitTest.XDG_DROPLET_COMBUSTION_TEST();
            //NUnit.Framework.Assert.AreEqual(true, false, "remove me");

            //BoSSS.Solution.Application<XNSEC_Control>._Main(new string[] { "--control", "cs:BoSSS.Application.XNSEC.FullNSEControlExamples.XDG_pseudo2dCombustion_MixtureFraction()", "--delplt" }, false, delegate () {
            //    var p = new XNSEC_MixtureFraction();
            //    return p;
            //});
            //NUnit.Framework.Assert.AreEqual(true, false, "remove me");

            //-n 8 ./XNSEC.exe -c "cs:BoSSS.Application.XNSEC.FullNSEControlExamples.XDG_DropletCombustion()"
            //System.Environment.Exit(111);

            bool MixtureFractionCalculation = false;
            try {
                // peek at control file and select correct solver depending on controlfile type
                // parse arguments
                args = ArgsFromEnvironmentVars(args);
                var CmdlineParseRes = Parser.Default.ParseArguments<CommandLineOptions>(args);
                bool argsParseSuccess = !CmdlineParseRes.Errors.Any();
                CommandLineOptions opt = CmdlineParseRes.Value;

                if (!argsParseSuccess) {
                    System.Environment.Exit(-1);
                }

                if (opt.ControlfilePath != null) {
                    opt.ControlfilePath = opt.ControlfilePath.Trim();
                }

                LoadControlFile(opt.ControlfilePath, out XNSEC_Control ctrlV2, out XNSEC_Control[] ctrlV2_ParameterStudy);
                MixtureFractionCalculation = ctrlV2 is XNSEC_MF_Control | ctrlV2_ParameterStudy is XNSEC_MF_Control[];
            } catch {
                Console.WriteLine("Error while determining control type, using default behavior for 'XNSEC_Control'");
            }

            if (MixtureFractionCalculation) {
                _Main(args, false, delegate () {
                    var p = new XNSEC_MixtureFraction();
                    return p;
                });
            } else {
                _Main(args, false, delegate () {
                    var p = new XNSEC();
                    return p;
                });
            }

            //MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        #region Operator configuration

        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            int D = GridData.SpatialDimension;

            if (iLevSet == 0) {
                // Main Difference to base implementation:
                //var levelSetVelocity = new LevelSetVelocityEvaporative("Phi", GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters, new XNSFE_OperatorConfiguration(Control);

                var config = new XNSEC_OperatorConfiguration(Control);
                //var levelSetVelocity = config.isEvaporation ?
                //    new LevelSetVelocityGeneralNonMaterial(VariableNames.LevelSetCG, GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters, config) :
                //    new LevelSetVelocity(VariableNames.LevelSetCG, GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                var levelSetVelocity =                     new LevelSetVelocity(VariableNames.LevelSetCG, GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                return levelSetVelocity;
            } else if (iLevSet == 1) {
                // +++++++++++++++++++++
                // the immersed boundary
                // +++++++++++++++++++++

                string[] VelocityNames = VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(1), VariableNames.VelocityVector(D)).ToArray();
                ScalarFunctionTimeDep[] VelFuncs = new ScalarFunctionTimeDep[D];
                for (int d = 0; d < D; d++) {
                    Control.InitialValues_EvaluatorsVec.TryGetValue(VelocityNames[d], out VelFuncs[d]);
                }

                ILevelSetParameter levelSetVelocity = new ExplicitLevelSetVelocity(VariableNames.LevelSetCGidx(1), VelFuncs);
                return levelSetVelocity;
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// - 4x the velocity degree if convection is included (quadratic term in convection times
        ///   density times test function yields quadruple order)
        /// - 2x the velocity degree in the Stokes case
        /// </summary>
        /// <remarks>
        /// Note: Sayes algorithm can be regarded as a nonlinear transformation to the [-1,1]
        /// reference Element. We transform $`\int f dx $` to the reference Element, $`\int f dx =
        /// \int f(T) |det D(T)| d\hat{x} $` Suppose $`f$` has degree $`n$` and suppose the
        /// transformation $`T$` has degree $`p$`, then the integrand in reference space has
        /// approximately degree $`\leq n * p + (p - 1) $` This is problematic, because we need to
        /// find $`\sqrt(n * p + (p - 1))$` roots of the level set function, if we want to integrate
        /// $`f$` exactly. This goes unnoticed when verifying the quadrature method via
        /// volume/surface integrals with constant $`f = 1$`. When evaluating a constant function,
        /// $`n = 0$`, the degree of the integrand immensely simplifies to $`(p - 1)$`.
        /// </remarks>
        public override int QuadOrder() {
            if (Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye
               && Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes) {
                throw new ArgumentException($"The XNSE solver is only verified for cut-cell quadrature rules " +
                    $"{XQuadFactoryHelper.MomentFittingVariants.Saye} and {XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes}; " +
                    $"you have set {Control.CutCellQuadratureType}, so you are notified that you reach into unknown territory; " +
                    $"If you do not know how to remove this exception, you should better return now!");
            }

            //QuadOrder

            int degU = VelocityDegree();
            int quadOrder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 4 : 2);
            if (this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                //See remarks
                quadOrder *= 2;
                quadOrder += 1;
            }
            return quadOrder;
        }


        /// <summary>
        /// Either fluids A and B; or A, B and solid C.
        /// </summary>
        protected override Array SpeciesTable {
            get {
                if (Control.UseImmersedBoundary) {
                    var r = new string[2, 2];
                    r[0, 0] = "A";
                    r[0, 1] = "C"; // solid
                    r[1, 0] = "B";
                    r[1, 1] = "C"; // also solid
                    return r;
                } else {
                    return new[] { "A", "B" };
                }
            }
        }

        /// <summary>
        /// Usually, the term "DG order of the calculation" means the velocity degree.
        /// </summary>
        protected int VelocityDegree() {
            int pVel;
            if (this.Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v)) {
                pVel = v.Degree;
            } else if (this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1)) {
                pVel = v1.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pVel;
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            int D = this.GridData.SpatialDimension;
            int pVel = VelocityDegree();
            int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
            int pTemp = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Temperature].Degree;
            int pMassFraction = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.MassFraction0].Degree;

            // configurations for velocity
            for (int d = 0; d < D; d++) {
                var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pVel },
                    mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                    //mode = MultigridOperator.Mode.Eye,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.VelocityVector(D)[d]) }
                };
                configsLevel.Add(configVel_d);
            }
            // configuration for pressure
            var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pPrs },
                mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                //mode = MultigridOperator.Mode.Eye,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Pressure) }
            };
            configsLevel.Add(configPres);

            // configuration for Temperature
            var confTemp = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pTemp },
                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                //mode = MultigridOperator.Mode.Eye,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Temperature) }
            };
            configsLevel.Add(confTemp);

            // configurations for Mass fractions
            int NumberOfSpecies = Control.NumberOfChemicalSpecies;
            var massFractionNames = VariableNames.MassFractions(NumberOfSpecies);
            for (int i = 0; i < (NumberOfSpecies); i++) {
                var configMF = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pMassFraction },
                    //mode = MultigridOperator.Mode.Eye,
                    mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(massFractionNames[i]) }
                };
                configsLevel.Add(configMF);
            }
        }

        #endregion Operator configuration

        #region Operator definition


        //protected LowMachCombustionMultiphaseBoundaryCondMap boundaryMap;
        protected IncompressibleBoundaryCondMap boundaryMap;
        
        private ThermalMultiphaseBoundaryCondMap m_thermBoundaryMap;

        protected override IncompressibleBoundaryCondMap GetBcMap() {
            if (boundaryMap == null)
                boundaryMap = new LowMachCombustionMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, new string[] { "A", "B" }, Control.NumberOfChemicalSpecies);



            return boundaryMap;
        }

        protected virtual void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            int quadOrder = QuadOrder();
            XNSEC_OperatorConfiguration config = new XNSEC_OperatorConfiguration(this.Control);

            GetBcMap();

            #region Equations of state

            var ChemicalModel = new OneStepChemicalModel(Control.VariableOneStepParameters, Control.YFuelInlet, Control.YOxInlet);

            if (boundaryMap.PhysMode == PhysicsMode.Combustion) {
                EoS_A = new MaterialLaw_MultipleSpecies(Control.MolarMasses, Control.MatParamsMode, Control.rhoOne, Control.R_gas, Control.T_ref_Sutherland, ChemicalModel, Control.cpRef, Control.HeatCapacityMode);
                EoS_B = new MaterialLaw_MultipleSpecies(Control.MolarMasses, Control.MatParamsMode, Control.rhoOne, Control.R_gas, Control.T_ref_Sutherland, ChemicalModel, Control.cpRef, Control.HeatCapacityMode);
            } else if (boundaryMap.PhysMode == PhysicsMode.MixtureFraction) {
                EoS_A = new MaterialLawMixtureFractionNew(Control.T_ref_Sutherland, Control.MolarMasses, Control.MatParamsMode, Control.rhoOne, Control.R_gas, Control.HeatRelease, Control.TOxInlet, Control.TFuelInlet, Control.YOxInlet, Control.YFuelInlet, Control.zSt, Control.CC, ChemicalModel, Control.cpRef, Control.smoothingFactor);
                EoS_B = new MaterialLawMixtureFractionNew(Control.T_ref_Sutherland, Control.MolarMasses, Control.MatParamsMode, Control.rhoOne, Control.R_gas, Control.HeatRelease, Control.TOxInlet, Control.TFuelInlet, Control.YOxInlet, Control.YFuelInlet, Control.zSt, Control.CC, ChemicalModel, Control.cpRef, Control.smoothingFactor);
            } else {
                throw new Exception("Wrong configuration");
            }

            //initialize EoS
            EoS_A.Initialize(Control.AmbientPressure);
            EoS_B.Initialize(Control.AmbientPressure);

            EoS_A.ConstantDensityValue = config.physParams.rho_A;
            EoS_B.ConstantDensityValue = config.physParams.rho_B;
            EoS_A.ConstantViscosityValue = config.physParams.mu_A;
            EoS_B.ConstantViscosityValue = config.physParams.mu_B;
            EoS_A.ConstantHeatConductivityValue = config.thermParams.k_A;
            EoS_B.ConstantHeatConductivityValue = config.thermParams.k_B;
            EoS_A.ConstantHeatCapacityValue = config.thermParams.c_A;
            EoS_B.ConstantHeatCapacityValue = config.thermParams.c_B;
            EoS_A.ConstantDiffusivityFactorVal = config.thermParams.k_A; // TODO
            EoS_B.ConstantDiffusivityFactorVal = config.thermParams.k_B; // TODO

            #endregion Equations of state

            // ============================
            // Momentum
            // ============================
            for (int d = 0; d < D; ++d) {
                DefineMomentumEquations(opFactory, config, d, D, lsUpdater);
                //Add Gravitation
                if (config.isGravity) {
                    var GravA = Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A, Control.GetGravity("A", d));
                    opFactory.AddParameter(GravA);
                    var GravB = Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B, Control.GetGravity("B", d));
                    opFactory.AddParameter(GravB);
                }
            }

            // ============================
            // Other Parameters
            // ============================
            DefineAditionalParameters(opFactory, config, D, lsUpdater, quadOrder);

            // ============================
            // Continuity
            // ============================
            if (config.isContinuity) {
                DefineContinuityEquation(opFactory, config, D, lsUpdater);
            }

            // ============================
            // Scalar Equations
            // ============================
            DefineScalarEquations(opFactory, config, D, lsUpdater);


            // ============================
            // Immersed BoundaryEquations
            // ============================
            if (Control.UseImmersedBoundary)
                DefineSystemImmersedBoundary(D, opFactory, config, lsUpdater);
        }

        protected virtual void DefineAditionalParameters(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater, int quadOrder) {
            // ============================== //
            // === additional parameters === //
            // ============================= //
            if (config.PlotAdditionalParameters) {
                opFactory.AddParameter(new Density(EoS_A, EoS_B, config.NoOfChemicalSpecies));
                opFactory.AddParameter(new Viscosity(EoS_A, EoS_B));
                opFactory.AddParameter(new HeatCapacity(EoS_A, EoS_B));
            }
            opFactory.AddCoefficient(new Solution.XNSECommon.SlipLengths(config, VelocityDegree()));




            if (config.isEvaporation) {
                var MassFluxExt = new MassFluxExtension_Evaporation(config);
                lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, MassFluxExt);
            }
            // ==================================== //
            // === level set related parameters === //
            // ==================================== //

            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);

            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, normalsParameter);
            var v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
            if (config.physParams.IncludeConvection && config.isTransport) {
                opFactory.AddParameter(v0Mean);
            }
            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, v0Mean);

            #region SurfaceTension

            switch (Control.AdvancedDiscretizationOptions.SST_isotropicMode) {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine:
                    MaxSigma maxSigmaParameter = new MaxSigma(Control.PhysicalParameters, Control.AdvancedDiscretizationOptions, QuadOrder(), Control.dtFixed);
                    opFactory.AddParameter(maxSigmaParameter);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, maxSigmaParameter);
                    BeltramiGradient lsBGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsBGradient);
                    break;

                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                    BeltramiGradient lsGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsGradient);
                    break;

                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                    BeltramiGradientAndCurvature lsGradientAndCurvature =
                        FromControl.BeltramiGradientAndCurvature(Control, "Phi", quadOrder, D);
                    opFactory.AddParameter(lsGradientAndCurvature);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsGradientAndCurvature);
                    break;

                case SurfaceStressTensor_IsotropicMode.Curvature_Fourier:
                    FourierLevelSet ls = (FourierLevelSet)lsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet;
                    var fourier = new FourierEvolver(
                        VariableNames.LevelSetCG,
                        ls,
                        Control.FourierLevSetControl,
                        Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Curvature].Degree);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, fourier);
                    //lsUpdater.AddEvolver(VariableNames.LevelSetCG, fourier);
                    opFactory.AddParameter(fourier);
                    break;

                default:
                    throw new NotImplementedException($"option {Control.AdvancedDiscretizationOptions.SST_isotropicMode} is not handled.");
            }

            #endregion SurfaceTension

            // ==================== //
            // === Coefficients === //
            // ==================== //
            opFactory.AddCoefficient(new ReynoldsNumber(config));

            opFactory.AddCoefficient(new EvapMicroRegion());
            if (config.prescribedMassflux != null)
                opFactory.AddCoefficient(new PrescribedMassFlux(config));
        }

        protected virtual void DefineContinuityEquation(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater) {
            opFactory.AddEquation(new LowMachContinuity(D, "A", config, boundaryMap, EoS_A, Control.dtFixed, Control.ManufacturedSolution_Continuity, Control.NonLinearSolver.SolverCode));
            opFactory.AddEquation(new LowMachContinuity(D, "B", config, boundaryMap, EoS_B, Control.dtFixed, Control.ManufacturedSolution_Continuity, Control.NonLinearSolver.SolverCode));

            //=== evaporation extension === //
            if (config.isEvaporation) {
                opFactory.AddEquation(new InterfaceContinuity_Evaporation_Newton_LowMach("A", "B", D, config));
            } else {
                opFactory.AddEquation(new InterfaceContinuityLowMach(config, D, LsTrk, config.isMatInt));

            }


            if (Control.timeDerivativeConti_OK) {
                var rho0 = new Density_t0(config.NoOfChemicalSpecies, (MaterialLaw_MultipleSpecies)EoS_A);
                opFactory.AddParameter(rho0);

                var rho00 = new Density_t00(config.NoOfChemicalSpecies, (MaterialLaw_MultipleSpecies)EoS_A);
                opFactory.AddParameter(rho00);
            }
            //lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, rho0);
        }

        protected virtual void DefineMomentumEquations(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int d, int D, LevelSetUpdater lsUpdater) {
            Func<double[], double, double> ManSol = d == 0 ? Control.ManufacturedSolution_MomentumX : Control.ManufacturedSolution_MomentumY;
            opFactory.AddEquation(new LowMachMomentumEquations("A", d, D, boundaryMap, config, EoS_A, ManSol, Control.NonLinearSolver.SolverCode));
            opFactory.AddEquation(new LowMachMomentumEquations("B", d, D, boundaryMap, config, EoS_B, ManSol, Control.NonLinearSolver.SolverCode));
            opFactory.AddEquation(new NSEInterface_LowMach("A", "B", d, D, boundaryMap, config, EoS_A, EoS_B, config.isMovingMesh));
            // opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, LsTrk, config)); // Maybe later...


             //=== evaporation extension === //
            if (config.isEvaporation) {
                //opFactory.AddEquation(new InterfaceNSE_Evaporation_Newton("A", "B", D, d, config));
                opFactory.AddEquation(new InterfaceNSE_Evaporation_LowMach("A", "B", D, d, config));

            }
        }

        /// <summary>
        /// Scalar equations used in the lowMach equations, namely Temperature equation and Mass Fractions equations
        /// </summary>
        public virtual void DefineScalarEquations(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater) {
            //================================
            // Energy equations (Temperature)
            //================================

            if (config.TemperatureEquationOK) {
                opFactory.AddEquation(new LowMachEnergy("A", D, boundaryMap, config, EoS_A, Control.HeatRelease, Control.ReactionRateConstants, Control.MolarMasses, Control.TRef, Control.cpRef, Control.dtFixed, Control.myThermalWallType, Control.ManufacturedSolution_Energy));
                opFactory.AddEquation(new LowMachEnergy("B", D, boundaryMap, config, EoS_B, Control.HeatRelease, Control.ReactionRateConstants, Control.MolarMasses, Control.TRef, Control.cpRef, Control.dtFixed, Control.myThermalWallType, Control.ManufacturedSolution_Energy));
                if (config.isEvaporation) {
                    opFactory.AddEquation(new HeatInterface_Evaporation_LowMach("A", "B", D, m_thermBoundaryMap, config, EoS_A,EoS_B, config.Reynolds,config.Prandtl, config.thermParams.T_sat ));
                } else {
                    opFactory.AddEquation(new HeatInterface_LowMach("A", "B", D,  boundaryMap, config, EoS_A, EoS_B,config.Reynolds,config.Prandtl));
                }

                opFactory.AddParameter(new dp0dt(EoS_A, Control.Reynolds, Control.Prandtl));
            } else {
                opFactory.AddEquation(new IdentityEquation("A", VariableNames.Temperature, EquationNames.HeatEquation));
                opFactory.AddEquation(new IdentityEquation("B", VariableNames.Temperature, EquationNames.HeatEquation));
            }
            opFactory.AddParameter(new ThermodynamicPressure(Control.InitialMass, Control.ThermodynamicPressureMode, EoS_A));

            //================================
            // Mass Fractions equations
            //================================
            
            for (int s = 0; s < config.NoOfChemicalSpecies; s++) {
                if (config.MassFractionEquationsOK) {
                    int chemicalSpeciesCounter = s;

                    opFactory.AddEquation(new LowMachMassFraction("A", D, boundaryMap, config, EoS_A, chemicalSpeciesCounter, Control.ReactionRateConstants, Control.StoichiometricCoefficients, Control.MolarMasses, Control));
                    opFactory.AddEquation(new LowMachMassFraction("B", D, boundaryMap, config, EoS_B, chemicalSpeciesCounter, Control.ReactionRateConstants, Control.StoichiometricCoefficients, Control.MolarMasses, Control));
                    if (Control.ChemicalReactionActive) {
                        opFactory.AddParameter(new ReactionRate(EoS_A, Control.ReactionRateConstants, config.NoOfChemicalSpecies));
                    }

                    if (config.isEvaporation) {
                        double[] valsAtInterface = new double[config.NoOfChemicalSpecies];
                        valsAtInterface[0] = 1.0;
                        valsAtInterface[1] = 0.0;

                        opFactory.AddEquation(new SpeciesMassTransferInterface_Evaporation_LowMach("A", "B", D, m_thermBoundaryMap, config, EoS_A, EoS_B, config.Reynolds, config.Prandtl, valsAtInterface[s],s));
                    } else {
                        //TODO
                        //opFactory.AddEquation(new HeatInterface_LowMach("A", "B", D, boundaryMap, config));
                    }
                } else {// Add identity equation for each MF
                    opFactory.AddEquation(new IdentityEquation("A", VariableNames.MassFractions(config.NoOfChemicalSpecies)[s], EquationNames.SpeciesMassBalanceName(s)));
                    opFactory.AddEquation(new IdentityEquation("B", VariableNames.MassFractions(config.NoOfChemicalSpecies)[s], EquationNames.SpeciesMassBalanceName(s)));
                }
            }

            //var p0_old = new ThermodynamicPressure_Oldtimestep(1.0, Control.ThermodynamicPressureMode, EoS_A);
            //opFactory.AddParameter(p0_old);
            //lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, p0_old);
        }

        /// <summary>
        /// 
        /// </summary>
        protected override int NoOfLevelSets {
            get {
                if (Control.UseImmersedBoundary)
                    return 2;
                else
                    return 1;
            }
        }


        /// <summary>
        /// Definition of the boundary condition on the immersed boundary (fluid-solid boundary, level-set 1), 
        /// <see cref="XNSE_Control.UseImmersedBoundary"/>;
        /// Override to customize.
        /// </summary>
        protected virtual void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, XNSEC_OperatorConfiguration config, LevelSetUpdater lsUpdater) {

            ////////////////////////////////////////////////
            ///     Momentum
            ////////////////////////////////////////////////
            for (int d = 0; d < D; ++d) {
                // so far only no slip!
                opFactory.AddEquation(new NSEimmersedBoundary_Newton_LowMach("A", "C", 1, d, D, boundaryMap, config, EoS_A, config.isMovingMesh,Control.physicsMode));
                opFactory.AddEquation(new NSEimmersedBoundary_Newton_LowMach("B", "C", 1, d, D, boundaryMap, config, EoS_B, config.isMovingMesh, Control.physicsMode));
            }
            ////////////////////////////////////////////////
            ///     Conti
            ////////////////////////////////////////////////
            opFactory.AddEquation(new ImmersedBoundaryContinuity_LowMach("A", "C", 1, config, D));
            opFactory.AddEquation(new ImmersedBoundaryContinuity_LowMach("B", "C", 1, config, D));
            ////////////////////////////////////////////////
            ///     Temperature
            ////////////////////////////////////////////////
            if (config.TemperatureEquationOK) {
                opFactory.AddEquation(new ImmersedBoundaryHeat_LowMach("A", "C", 1, D, config,EoS_A));
                opFactory.AddEquation(new ImmersedBoundaryHeat_LowMach("B", "C", 1, D, config, EoS_B));
            }
            ////////////////////////////////////////////////
            ///     MassFractions
            ////////////////////////////////////////////////
          

            for (int s = 0; s < config.NoOfChemicalSpecies; s++) {
                if (config.MassFractionEquationsOK) {
                    int chemicalSpeciesCounter = s;
                    opFactory.AddEquation(new ImmersedBoundaryMF_LowMach("A", "C", 1, D, config, EoS_A, chemicalSpeciesCounter,config.NoOfChemicalSpecies));
                    opFactory.AddEquation(new ImmersedBoundaryMF_LowMach("B", "C", 1, D, config, EoS_B, chemicalSpeciesCounter, config.NoOfChemicalSpecies));

                }  
            }
            opFactory.AddParameter((ParameterS)GetLevelSetVelocity(1));
        }


        /// <summary>
        /// Low-Mach unsteady part definition
        /// </summary>
        /// <param name="D"></param>
        /// <param name="opFactory"></param>
        /// <param name="lsUpdater"></param>
        protected virtual void DefineTemporalTerm(int D, OperatorFactory opFactory) {
            //  var EoS = base.Control.EoS;
            int NoOfChemSpecies = Control.NumberOfChemicalSpecies;
            // dbg_launch();

            if (boundaryMap.PhysMode == PhysicsMode.Combustion) {

                // Momentum
                // ============================
                for (int d = 0; d < D; d++) {
                    opFactory.AddEquation(new LowMachUnsteadyEquationPart("A", D, VariableNames.VelocityVector(D)[d], EquationNames.MomentumEquationComponent(d), NoOfChemSpecies, EoS_A));
                    opFactory.AddEquation(new LowMachUnsteadyEquationPart("B", D, VariableNames.VelocityVector(D)[d], EquationNames.MomentumEquationComponent(d), NoOfChemSpecies, EoS_B));
                }

                // Continuity
                // ============================
                opFactory.AddEquation(new LowMachUnsteadyEquationPart("A", D, VariableNames.Pressure, EquationNames.ContinuityEquation, NoOfChemSpecies, EoS_A, massScale: 0.0));
                opFactory.AddEquation(new LowMachUnsteadyEquationPart("B", D, VariableNames.Pressure, EquationNames.ContinuityEquation, NoOfChemSpecies, EoS_B, massScale: 0.0));


                // Energy (Temperature)
                // ============================
                opFactory.AddEquation(new LowMachUnsteadyEquationPart("A", D, VariableNames.Temperature, EquationNames.HeatEquation, NoOfChemSpecies, EoS_A, massScale: 1.0 / Control.HeatCapacityRatio));
                opFactory.AddEquation(new LowMachUnsteadyEquationPart("B", D, VariableNames.Temperature, EquationNames.HeatEquation, NoOfChemSpecies, EoS_B, massScale: 1.0 / Control.HeatCapacityRatio));
                // Mass Fractions
                // ============================
                for (int s = 0; s < NoOfChemSpecies; s++) {
                    opFactory.AddEquation(new LowMachUnsteadyEquationPart("A", D, VariableNames.MassFractions(NoOfChemSpecies)[s], EquationNames.SpeciesMassBalanceName(s), NoOfChemSpecies, EoS_A));
                    opFactory.AddEquation(new LowMachUnsteadyEquationPart("B", D, VariableNames.MassFractions(NoOfChemSpecies)[s], EquationNames.SpeciesMassBalanceName(s), NoOfChemSpecies, EoS_B));
                }
            } else {

                // Momentum
                // ============================
                for (int d = 0; d < D; d++) {
                    opFactory.AddEquation(new LowMachUnsteadyEquationPart_MF("A", D, VariableNames.VelocityVector(D)[d], EquationNames.MomentumEquationComponent(d), NoOfChemSpecies, EoS_A));
                    opFactory.AddEquation(new LowMachUnsteadyEquationPart_MF("B", D, VariableNames.VelocityVector(D)[d], EquationNames.MomentumEquationComponent(d), NoOfChemSpecies, EoS_B));
                }

                // Continuity
                // ============================
                opFactory.AddEquation(new LowMachUnsteadyEquationPart_MF("A", D, VariableNames.Pressure, EquationNames.ContinuityEquation, NoOfChemSpecies, EoS_A, massScale: 0.0));
                opFactory.AddEquation(new LowMachUnsteadyEquationPart_MF("B", D, VariableNames.Pressure, EquationNames.ContinuityEquation, NoOfChemSpecies, EoS_B, massScale: 0.0));



                // Mixture Fraction
                // ============================
                opFactory.AddEquation(new LowMachUnsteadyEquationPart_MF("A", D, VariableNames.MixtureFraction, EquationNames.MixtureFractionEquation, NoOfChemSpecies, EoS_A, massScale: 1.0 ));
                opFactory.AddEquation(new LowMachUnsteadyEquationPart_MF("B", D, VariableNames.MixtureFraction, EquationNames.MixtureFractionEquation, NoOfChemSpecies, EoS_B, massScale: 1.0 ));
            }
        }

        #endregion Operator definition

        public MaterialLaw_MultipleSpecies EoS_A;
        public MaterialLaw_MultipleSpecies EoS_B;

        private XDifferentialOperatorMk2 XOP;
        /// <summary>
        /// Low-Mach system of equations definition
        /// </summary>
        /// <param name="D"></param>
        /// <param name="opFactory"></param>
        /// <param name="lsUpdater"></param>

        protected override XDifferentialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {
            OperatorFactory opFactory = new OperatorFactory();

            DefineSystem(D, opFactory, levelSetUpdater);

            /*XSpatialOperatorMk2*/
            XOP = opFactory.GetSpatialOperator(QuadOrder());
            //final settings
            XOP.FreeMeanValue[VariableNames.Pressure] = !GetBcMap().DirichletPressureBoundary;



            if (Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                Console.WriteLine("Linearization Hint:" + LinearizationHint.GetJacobiOperator.ToString());
                XOP.LinearizationHint = LinearizationHint.GetJacobiOperator;
            } else if (Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard){
                Console.WriteLine("Linearization Hint:" + LinearizationHint.GetJacobiOperator.ToString());

                //throw new NotImplementedException("LowMach solver supports only Newton as NonLinearSolver");
            }

            XOP.ParameterUpdates.Add(PlotNewtonIterationsHack);

            XOP.IsLinear = false;
            XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;

            // ============================
            // Self made temporal operator
            // ============================
            if (Control.UseSelfMadeTemporalOperator && (base.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient)) {
                Console.WriteLine("Using low Mach temporal operator");
                OperatorFactory temporalOperatorFactory = new OperatorFactory();
                DefineTemporalTerm(D, temporalOperatorFactory);
                XDifferentialOperatorMk2 temporalXOP = temporalOperatorFactory.GetSpatialOperator(QuadOrder());
                temporalXOP.Commit();

                var DependentTemporalOp = new DependentXTemporalOperator(XOP);

                foreach (var c in temporalXOP.CodomainVar) {
                    foreach (var d in temporalXOP.EquationComponents[c]) {
                        DependentTemporalOp.EquationComponents[c].Add(d);
                    }
                }

                XOP.TemporalOperator = DependentTemporalOp;
            }

            ////============================
            //// Solver-Controlled Homotopy<
            ////============================

            if (this.Control.HomotopyApproach == XNSEC_Control.HomotopyType.Automatic) {
                if (Control.HomotopyVariable == XNSEC_Control.HomotopyVariableEnum.Reynolds) {
                    this.CurrentHomotopyValue = Control.Reynolds;
                    XOP.HomotopyUpdate.Add(delegate (double HomotopyScalar) {
                        if (HomotopyScalar < 0.0)
                            throw new ArgumentOutOfRangeException();
                        if (HomotopyScalar > 1.0)
                            throw new ArgumentOutOfRangeException();

                        // Using a linear function to prescribe the homotopy  path
                        // If HomotopyScalar = 0 => Reynolds = StartingValue
                        // If HomotopyScalar = 1 => Reynolds = Control.Reynolds
                        double StartingValue = Control.StartingHomotopyValue; // this should be an "easy" value for finding a solution

                        //double StartingValue = Control.Reynolds/10; // this should be an "easy" value for finding a solution
                        double AimedValue = Control.Reynolds;

                        //Linear
                        double slope = (AimedValue - StartingValue) / (1 - 0);
                        double val = slope * (HomotopyScalar - 0) + StartingValue;
                        this.CurrentHomotopyValue = val;

                        //////Exponential
                        //Console.WriteLine("Updating the homotopy value using a Exponential function ");
                        //double slope = (Math.Log10(AimedValue) - Math.Log10(StartingValue)) / (1 - 0);
                        //double reExponent = slope * (HomotopyScalar - 0) + Math.Log10(StartingValue);
                        //this.CurrentHomotopyValue = Math.Pow(10, reExponent);

                        Console.WriteLine("HomotopyScalar:" + HomotopyScalar);
                        Console.WriteLine("HomotopyValue:" + CurrentHomotopyValue);
                    });
                    var defaultcoefficients = XOP.OperatorCoefficientsProvider;
                    XOP.OperatorCoefficientsProvider = delegate (LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
                        CoefficientSet cs = defaultcoefficients(lstrk, spc, quadOrder, TrackerHistoryIdx, time);
                        cs.UserDefinedValues["Reynolds"] = this.CurrentHomotopyValue;
                        return cs;
                    };
                }

                if (Control.HomotopyVariable == XNSEC_Control.HomotopyVariableEnum.VelocityInletMultiplier) {
                    this.CurrentHomotopyValue = Control.homotopieAimedValue;
                    XOP.HomotopyUpdate.Add(delegate (double HomotopyScalar) {
                        if (HomotopyScalar < 0.0)
                            throw new ArgumentOutOfRangeException();
                        if (HomotopyScalar > 1.0)
                            throw new ArgumentOutOfRangeException();

                        ////Linear
                        //double StartingValue = 1.0; // this should be an "easy" value for finding a solution
                        //double AimedValue = Control.homotopieAimedValue;
                        //double slope = (AimedValue - StartingValue) / (1 - 0);
                        //double val = slope * (HomotopyScalar - 0) + StartingValue;

                        double StartingValue = 1.0 / Control.homotopieAimedValue; // this should be an "easy" value for finding a solution
                        double AimedValue = 1.0;
                        double slope = (AimedValue - StartingValue) / (1 - 0);
                        double val = slope * (HomotopyScalar - 0) + StartingValue;

                        this.CurrentHomotopyValue = val;

                        Console.WriteLine("HomotopyScalar:" + HomotopyScalar);
                        Console.WriteLine("HomotopyValue:" + CurrentHomotopyValue);
                    });
                    var defaultcoefficients = XOP.OperatorCoefficientsProvider;
                    XOP.OperatorCoefficientsProvider = delegate (LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
                        CoefficientSet cs = defaultcoefficients(lstrk, spc, quadOrder, TrackerHistoryIdx, time);
                        cs.UserDefinedValues["VelocityMultiplier"] = this.CurrentHomotopyValue;
                        return cs;
                    };
                }

                if (Control.HomotopyVariable == XNSEC_Control.HomotopyVariableEnum.HeatOfCombustion) {
                    this.CurrentHomotopyValue = Control.homotopieAimedValue;
                    XOP.HomotopyUpdate.Add(delegate (double HomotopyScalar) {
                        if (HomotopyScalar < 0.0)
                            throw new ArgumentOutOfRangeException();
                        if (HomotopyScalar > 1.0)
                            throw new ArgumentOutOfRangeException();

                        ////Linear

                        double StartingValue = 0.0; // this should be an "easy" value for finding a solution
                        double AimedValue = Control.HeatRelease;
                        double slope = (AimedValue - StartingValue) / (1 - 0);
                        double val = slope * (HomotopyScalar - 0) + StartingValue;

                        this.CurrentHomotopyValue = val;
                        ((MaterialLawMixtureFractionNew)EoS_A).Q = val;
                        Console.WriteLine("HomotopyScalar:" + HomotopyScalar);
                        Console.WriteLine("HomotopyValue:" + CurrentHomotopyValue);
                    });
                    var defaultcoefficients = XOP.OperatorCoefficientsProvider;
                    XOP.OperatorCoefficientsProvider = delegate (LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
                        CoefficientSet cs = defaultcoefficients(lstrk, spc, quadOrder, TrackerHistoryIdx, time);
                        cs.UserDefinedValues["HeatOfReaction"] = this.CurrentHomotopyValue;
                        return cs;
                    };
                }
            }

            ////============================
            //// Solver safe guard 
            ////============================
            if (Control.VariableBounds != null) {
                Console.WriteLine("Using solver safe guard!");
                XOP.SolverSafeguard = DelValidationCombustion;
            }
            XOP.Commit();

            PrintConfiguration();
            return XOP;
        }

        /// <summary>
        /// duh
        /// </summary>
        public double CurrentHomotopyValue {
            get {
                return (double)XOP.UserDefinedValues["A"][Control.homotopieVariableName];
            }
            set {
                double oldVal;
                if (XOP.UserDefinedValues["A"].ContainsKey(Control.homotopieVariableName))
                    oldVal = CurrentHomotopyValue;
                else
                    oldVal = double.NegativeInfinity;

                if (oldVal != value)
                    Console.WriteLine("setting" + Control.homotopieVariableName + " to " + value);
                XOP.UserDefinedValues["A"][Control.homotopieVariableName] = value;
            }
        }

        int homotopyStep = 0;
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls


           

            dt = GetTimestep();
            // Convert mixture fraction into temperature and mass fractions
            if (Control is XNSEC_MF_Control) {
                string[] names = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(Control.NumberOfChemicalSpecies));
                XDGField MixtureFraction = (XDGField)this.m_IOFields.Where(f => f.Identification == VariableNames.MixtureFraction).Single();
                Console.WriteLine("transforming back variables for MF calculation");


                foreach (var id in names) {
                    var field = (XDGField)(m_IOFields.Where(f => f.Identification == id).SingleOrDefault());

                    field.Clear();
                    var fieldToTransform_A = ((XDGField)field).GetSpeciesShadowField("A");
                    fieldToTransform_A.ProjectField(1.0,
                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                        int K = result.GetLength(1);
                        MultidimensionalArray ZArr = MultidimensionalArray.Create(Len, K);
                        MixtureFraction.Evaluate(j0, Len, NS, ZArr);
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                result[j, k] = EoS_A.getVariableFromZ(ZArr[j, k], id);
                            }
                        }
                    }, new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, null));


                    var fieldToTransform_B = ((XDGField)field).GetSpeciesShadowField("B");
                    fieldToTransform_B.ProjectField(1.0,
                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                        int K = result.GetLength(1);
                        MultidimensionalArray ZArr = MultidimensionalArray.Create(Len, K);
                        MixtureFraction.Evaluate(j0, Len, NS, ZArr);
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                result[j, k] = EoS_A.getVariableFromZ(ZArr[j, k], id);
                            }
                        }
                    }, new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, null));


                }


            }

            if (Control.timeDerivativeEnergyp0_OK) {
                //    var p0_old = this.Parameters.Where(f => f.Identification == VariableNames.ThermodynamicPressure + "_t0").Single();
                //    var p0 = this.Parameters.Where(f => f.Identification == VariableNames.ThermodynamicPressure).Single();
                //    p0_old.Clear();
                //    p0_old.Acc(1.0, p0);
            }

        
            if (Control.timeDerivativeConti_OK && base.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
                XOP.InvokeParameterUpdate(phystime, CurrentStateVector.Fields.ToArray(), Parameters.ToArray());
                var rho0_oldold = this.Parameters.Where(f => f.Identification == VariableNames.Rho + "_t00").Single();
                var rho0_old = this.Parameters.Where(f => f.Identification == VariableNames.Rho + "_t0").Single();
                var rho = this.Parameters.Where(f => f.Identification == VariableNames.Rho).Single();


                rho0_oldold.Clear();
                rho0_oldold.Acc(1.0, rho0_old);

                rho0_old.Clear();
                rho0_old.Acc(1.0, rho);
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!");
                Console.WriteLine("Rho mean value: " + rho.GetMeanValueTotal(null));
                Console.WriteLine("Rho old mean value: " + rho0_old.GetMeanValueTotal(null));
                Console.WriteLine("Rho old old mean value: " + rho0_oldold.GetMeanValueTotal(null));
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!");
            }



            // Selfmade homotopy
            if (this.Control.HomotopyApproach == XNSEC_Control.HomotopyType.Manual) {
                this.CurrentHomotopyValue= this.Control.SelfDefinedHomotopyArray[homotopyStep];
                Console.WriteLine("Setting reynolds number to " + this.CurrentHomotopyValue);
                var defaultcoefficients = XOP.OperatorCoefficientsProvider;
                XOP.OperatorCoefficientsProvider = delegate (LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
                    CoefficientSet cs = defaultcoefficients(lstrk, spc, quadOrder, TrackerHistoryIdx, time);
                    cs.UserDefinedValues[Control.homotopieVariableName] = this.CurrentHomotopyValue;
                    return cs;
                };
                homotopyStep++;
            }

            var overallstart = DateTime.Now;
            Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt}");
            bool SolverSuccess = Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            var overallstop = DateTime.Now;
            var overallduration = overallstop - overallstart;

            //int D = this.Grid.SpatialDimension;
            //double dtCFL = this.GridData.ComputeCFLTime(this.CurrentStateVector.Fields.Take(D), double.MaxValue);
            //Console.WriteLine("CFL time is " + dtCFL);

            Console.WriteLine("Duration of this timestep: " + overallduration);
            Console.WriteLine($"done with time step {TimestepNo}");

            if (Control.AnalyticsolutionSwitch || !Control.ExactSolutionVelocity.IsNullOrEmpty()) {
                CalcErrors();
            }

            //Calculate nusselt number
            if ((Control.EdgeTagsNusselt != null) && Control.TimesteppingMode == AppControl._TimesteppingMode.Steady) {
                Console.WriteLine("Calculating nusselt numbers!");
                var temperatureXdg = (XDGField)(CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.Temperature).SingleOrDefault());
                var temp = temperatureXdg.ProjectToSinglePhaseField(4);
                var NusseltResults = CalculateNusselt(TimestepNo, base.GridData, temp, Control);
                this.CurrentSessionInfo.KeysAndQueries.Add("NusseltNumber0", NusseltResults[0]);
                this.CurrentSessionInfo.KeysAndQueries.Add("NusseltNumber1", NusseltResults[1]);
                this.CurrentSessionInfo.KeysAndQueries.Add("NusseltNumber2", NusseltResults[2]);

                Console.WriteLine("Nusselt0:" + NusseltResults[0]);
                Console.WriteLine("Nusselt1:" + NusseltResults[1]);
                Console.WriteLine("Nusselt2:" + NusseltResults[2]);
            }

            double minTemperature = 0; double maxTemperature = 0;
            try {
                CurrentState.Fields.Where(f => f.Identification == VariableNames.Temperature).Single().GetExtremalValues(out minTemperature, out maxTemperature);
                Console.WriteLine("Min Temperature in this timestep: " + minTemperature);
                Console.WriteLine("Max Temperature in this timestep: " + maxTemperature);
            } catch (Exception e) {

            }
            //sensor.Update(CurrentState.Fields.Where(f => f.Identification == VariableNames.Temperature).Single());
            return dt;
        }

       

        private int hack_TimestepIndex = 0;
        //private PerssonSensor sensor;

        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            if (Control.UseMixtureFractionsForCombustionInitialization) {
                var ChemicalModel = new OneStepChemicalModel(Control.VariableOneStepParameters, Control.YFuelInlet, Control.YOxInlet);
                var m_EoS = new MaterialLawMixtureFractionNew(Control.T_ref_Sutherland, Control.MolarMasses, Control.MatParamsMode, Control.rhoOne, Control.R_gas, Control.HeatRelease, Control.TOxInlet, Control.TFuelInlet, Control.YOxInlet, Control.YFuelInlet, Control.zSt, Control.CC, ChemicalModel, Control.cpRef, Control.smoothingFactor);

                string[] names = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(Control.NumberOfChemicalSpecies));

                // Start the combustion calculation with the mixture fraction
                var MixtureFraction = this.m_IOFields.Where(f => f.Identification == VariableNames.MixtureFraction).Single();

                //var basis = this.m_RegisteredFields.Where(f => f.Identification == VariableNames.Temperature).First().Basis;
                //var MixtureFraction = new XDGField((XDGBasis)basis, VariableNames.MixtureFraction);
                base.RegisterField(MixtureFraction, IOListOption.Always);

                foreach (var id in names) {
                    var fieldToTransform = FindField(this.m_RegisteredFields.ToArray(), id);
                    var fieldToTransform_A = ((XDGField)fieldToTransform).GetSpeciesShadowField("A");
                    fieldToTransform_A.Clear();
                    fieldToTransform_A.ProjectField(1.0,
                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                        int K = result.GetLength(1);
                        MultidimensionalArray ZArr = MultidimensionalArray.Create(Len, K);
                        MixtureFraction.Evaluate(j0, Len, NS, ZArr);
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                result[j, k] = m_EoS.getVariableFromZ(ZArr[j, k], id);
                            }
                        }
                    }, new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, null));
                }
            }
       
        
        }

        protected override void CreateFields() {
     
            base.CreateFields();

            if (Control.UseMixtureFractionsForCombustionInitialization /*|| Control is XNSEC_MF_Control*/) {
                base.RegisterField(new XDGField((XDGBasis)(this.m_RegisteredFields.Where(f => f.Identification == VariableNames.Temperature).First().Basis), VariableNames.MixtureFraction), IOListOption.Always);
            }

            // Create fields for analytical solution and errors
            if (Control.AnalyticsolutionSwitch || !Control.ExactSolutionVelocity.IsNullOrEmpty()) {
                // Errors:
                var errFields = InstantiateErrorFields();
                foreach (var f in errFields) {
                    base.RegisterField(f, IOListOption.Always);
                }
                // Analytical Solutions
                var AnSolFields = InstantiateAnalyticalSolFields();
                foreach (var f in AnSolFields) {
                    base.RegisterField(f, IOListOption.Always);
                }
            }

            if (Control is XNSEC_MF_Control) {
                Console.WriteLine("Instantiating result fields");
                string[] names = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(Control.NumberOfChemicalSpecies));
                var MixtureFraction = this.m_IOFields.Where(f => f.Identification == VariableNames.MixtureFraction).Single();
                base.RegisterField(MixtureFraction, IOListOption.Always);

                foreach (var id in names) {
                    var field = m_IOFields.Where(f => f.Identification == id).SingleOrDefault();
                    if (field == null) {
                        field = MixtureFraction.CloneAs();
                        field.Identification = id;
                        field.Clear();
                        base.RegisterField(field, IOListOption.Always);
                        //base.IOFields.Add(field);
                    }
                }
            }
        }

        /// <summary>
        /// Operator stability analysis
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis(OperatorAnalysisConfig config) {
            return this.Operator.OperatorAnalysis(this.CurrentStateVector.Mapping, config, this.MultigridOperatorConfig);
        }


        //protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
        //    // Cells Numbers
        //    var CellNumbers = this.m_RegisteredFields.Where(s => s.Identification == "CellNumbers").SingleOrDefault();
        //    if (CellNumbers == null) {
        //        CellNumbers = new SinglePhaseField(new Basis(this.GridData, 0), "CellNumbers");
        //        this.RegisterField(CellNumbers);
        //    }
        //    CellNumbers.Clear();
        //    CellNumbers.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
        //        int K = result.GetLength(1); // No nof Nodes
        //        for (int j = 0; j < Len; j++) {
        //            for (int k = 0; k < K; k++) {
        //                result[j, k] = j0 + j;
        //            }
        //        }
        //    }, new CellQuadratureScheme());


        //    //DGField[] RefinedFields = new[] { Refined_u, Refined_TestData, Refined_Grad_u[0], Refined_Grad_u[1], Refined_MagGrad_u };
        //    //string filename2 = "RefinedGrid." + timestepNo;
        //    //Tecplot.PlotFields(RefinedFields, filename2, physTime, superSampling);
        //}

        /// <summary>
        /// User-defined validation of a solver step, e.g. to prevent the solver to iterate out-of-bounds,
        /// e.g. to avoid un-physical 'solutions' (e.g. negative density).
        /// ('safeguard' for the solver)
        /// </summary>
        /// <param name="varIn"></param>
        /// <param name="varOut"></param>
        private void DelValidationCombustion(DGField[] varIn, DGField[] varOut) {
            CellMask AllCells = CellMask.GetFullMask(GridData);
            var Bounds = Control.VariableBounds;
            varOut.Clear();
            varOut = varIn.CloneAs();
            int idx = 0;
            foreach (var f in varIn) { // Loop over each DGField of the solution array
                foreach (var varName in Bounds.Keys) { // Iterate over dictionary with fields to be repaired
                    if (f.Identification == varName) {
                        double MinBound = Bounds[varName].Item1;
                        double MaxBound = Bounds[varName].Item2;

                        double[] mins = new double[AllCells.NoOfItemsLocally];
                        double[] maxs = new double[AllCells.NoOfItemsLocally];
                        f.GetCellwiseExtremalValues(mins, maxs);

                        // now check if values are inside bounds. If not, repair them
                        int lowcount = 0; int topcount = 0;
                        foreach (var cell in AllCells.ItemEnum) {
                            bool BoundedLow = (mins[cell] > MinBound);
                            bool BoundedTop = (maxs[cell] < MaxBound);

                            if (!BoundedLow) { // repair
                                ((XDGField)(varOut[idx])).SetMeanValueAB(cell, MinBound);
                                //varOut[idx].SetMeanValue(cell, MinBound);
                                lowcount++;
                            }
                            if (!BoundedTop) { // repair
                                ((XDGField)(varOut[idx])).SetMeanValueAB(cell, MaxBound);
                                //varOut[idx].SetMeanValue(cell, MaxBound);
                                topcount++;
                            }
                        }
                        if (lowcount > 0 || topcount > 0) {
                            Console.WriteLine(f.Identification + " is out of bounds in {0} cells.", lowcount + topcount);
                        }
                    }
                }
                idx++;
            }

            return;
        }
    }
}