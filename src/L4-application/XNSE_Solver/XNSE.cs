using BoSSS.Application.XNSE_Solver.LoadBalancing;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using CommandLine;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver {


    /// <summary>
    /// Multiphase-XDG-solver, with features:
    /// - incompressible two-phase flows.
    /// - solid immersed boundaries (planned).
    /// - three phase contact lines at the domain boundary
    /// - three phase contact lines at the intersection of the immersed solid boundary 
    /// - the generic control parameter <typeparamref name="T"/> allows derivations of this solver
    /// </summary>
    /// <remarks>
    /// Development history:
    /// - Current (jan2021) Maintainers: Beck, Rieckmann, Kummer
    /// - successor of the old XNSE solver <see cref="XNSE_SolverMain"/>, which was mainly used for SFB 1194 and PhD thesis of M. Smuda.
    /// - Quadrature order: saye algorithm can be regarded as a nonlinear transformation to the [-1,1] reference Element. 
    ///   We transform $` \int f dx $` to the reference Element, $` \int f dx = \int f(T) |det D(T)| d\hat{x} $`
    ///   Suppose f has degree n and suppose the transformation T has degree p, then the integrand in reference space
    ///   has approximately degree <= n * p + (p - 1)
    ///   This is problematic, because we need to find sqrt(n * p + (p - 1)) roots of the level set function, if we want to integrate f exactly.
    ///   This goes unnoticed when verifying the quadrature method via volume/surface integrals with constant f = 1.
    ///   When evaluating a constant function, n = 0, the degree of the integrand immensely simplifies to (p - 1).
    /// - see also: Extended discontinuous Galerkin methods for two-phase flows: the spatial discretization, F. Kummer, IJNME 109 (2), 2017. 
    /// 
    /// Quick start instructions:
    /// - The Fluid-Fluid-Interphase is level-set No. 0 and is named `Phi` 
    ///   in the initial values <see cref="AppControl.InitialValues"/>, see <see cref="DGField.Identification"/>)
    /// - The immersed boundary interface has the name `Phi2`
    ///   `Phi` is always activated and cannot be turned of; `Phi2` can be activated/deactivated via <see cref="XNSE_Control.UseImmersedBoundary"/>.
    /// - One should not need to set an initial value for `Phi` for simulations where no fluid-fluid interface is required:
    ///   If `Phi` it is identical to 0.0 in the entire domain, it will be set to -1, i.e. the entire domain will be assigned to species "A".
    /// - The vectorial velocity of `Phi2` is set via the fields `VelocityX@Phi2`, `VelocityY@Phi2` and `VelocityZ@Phi2` (in 3D). 
    ///   By specifying a vectorial velocity, one can also specify tangential velocities of the surface. 
    ///   The normal component should match the normal velocity which can be computed from `Phi2`,
    ///   yielding following compatibility condition: 
    /// ```math
    ///     \vec{v} \cdot \frac{\nabla \varphi_2}{| \nabla \varphi_2 |} =  \frac{- \partial_t \varphi_2}{| \nabla \varphi_2 |} 
    /// ```
    /// </remarks>
    public class XNSE : XNSE<XNSE_Control> {

        // ===========
        //  Main file
        // ===========
        static void Main(string[] args) {


            //InitMPI();
            //BoSSS.Application.XNSE_Solver.Tests.LevelSetUnitTests.LevelSetAdvectionTest2D(2, 1, LevelSetEvolution.FastMarching, LevelSetHandling.LieSplitting, false);
            //BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ChannelTest(3, 0.1, ViscosityMode.FullySymmetric, 0.0, XQuadFactoryHelper.MomentFittingVariants.Saye, NonLinearSolverCode.Newton);
            //DeleteOldPlotFiles();
            //BoSSS.Application.XNSE_Solver.Tests.LevelSetUnitTests.LevelSetAdvectionTest2D(3, 2, LevelSetEvolution.StokesExtension, LevelSetHandling.LieSplitting, false);
            //BoSSS.Application.XNSE_Solver.Legacy.LegacyTests.UnitTest.BcTest_PressureOutletTest(2, 1, 0.1d, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, SurfaceStressTensor_IsotropicMode.Curvature_Projected, false);
            //Tests.ASUnitTest.CurvedElementsTest(3);
            //Tests.ASUnitTest.IBMChannelTest(1, 0.0d, NonLinearSolverCode.Newton);
            //Tests.ASUnitTest.MovingDropletTest_rel_p3_Saye_FullySymmetric(0.1, true, SurfaceStressTensor_IsotropicMode.Curvature_Projected, 0.70611, true, false);
            //Tests.LevelSetUnitTests.LevelSetAdvectionTest2D(4, 2, LevelSetEvolution.StokesExtension, LevelSetHandling.LieSplitting, false);
            ////Tests.LevelSetUnitTests.LevelSetAdvectionOnWallTest3D(Math.PI / 4, 2, 0, LevelSetEvolution.FastMarching, LevelSetHandling.LieSplitting);
            ////Tests.LevelSetUnitTests.LevelSetShearingTest(2, 3, LevelSetEvolution.FastMarching, LevelSetHandling.LieSplitting);
            //BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.IBMChannelSolverTest(1, 0.0d, LinearSolverCode.exp_gmres_levelpmg);
            //throw new Exception("Remove me");

            bool Evap = false;
            // not sure if this works always, idea is to determine on startup which solver should be run.
            // default is XNSE<XNSE_Control>
            try {
                // peek at control file and select correct solver depending on controlfile type
                // parse arguments
                args = ArgsFromEnvironmentVars(args);
                CommandLineOptions opt = new CommandLineOptions();
                ICommandLineParser parser = new CommandLine.CommandLineParser(new CommandLineParserSettings(Console.Error));
                bool argsParseSuccess;
                argsParseSuccess = parser.ParseArguments(args, opt);

                if (!argsParseSuccess) {
                    System.Environment.Exit(-1);
                }

                if (opt.ControlfilePath != null) {
                    opt.ControlfilePath = opt.ControlfilePath.Trim();
                }

                XNSE_Control ctrlV2 = null;
                XNSE_Control[] ctrlV2_ParameterStudy = null;

                LoadControlFile(opt.ControlfilePath, out ctrlV2, out ctrlV2_ParameterStudy);
                Evap = ctrlV2 is XNSFE_Control | ctrlV2_ParameterStudy is XNSFE_Control[];
            } catch {
                Console.WriteLine("Error while determining control type, using default behavior for 'XNSE_Control'");
            }

            if (Evap) {
                XNSFE<XNSFE_Control>._Main(args, false, delegate () {
                    var p = new XNSFE<XNSFE_Control>(); 
                    return p;
                });
            } else {
                XNSE._Main(args, false, delegate () {
                    var p = new XNSE();
                    return p;
                });
            }
            
        }
    }

    /// <summary>
    /// Generic versions which should be used for derivatives 
    /// </summary>
    public class XNSE<T> : SolverWithLevelSetUpdater<T> where T : XNSE_Control, new() {

        /// <summary>
        /// - 3x the velocity degree if convection is included (quadratic term in convection times test function yields triple order)
        /// - 2x the velocity degree in the Stokes case
        /// </summary>
        /// <remarks>
        /// Note: 
        /// Sayes algorithm can be regarded as a nonlinear transformation to the [-1,1] reference Element. 
        /// We transform $`\int f dx $` to the reference Element, $`\int f dx = \int f(T) |det D(T)| d\hat{x} $`
        /// Suppose $`f$` has degree $`n$` and suppose the transformation $`T$` has degree $`p$`, then the integrand in reference space
        /// has approximately degree $`\leq n * p + (p - 1) $`
        /// This is problematic, because we need to find $`\sqrt(n * p + (p - 1))$` roots of the level set function, if we want to integrate $`f$` exactly.
        /// This goes unnoticed when verifying the quadrature method via volume/surface integrals with constant $`f = 1$`.
        /// When evaluating a constant function, $`n = 0$`, the degree of the integrand immensely simplifies to $`(p - 1)$`.        
        /// </remarks>
        override public int QuadOrder() {
            if(Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye
               && Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes) {
                throw new ArgumentException($"The XNSE solver is only verified for cut-cell quadrature rules " +
                    $"{XQuadFactoryHelper.MomentFittingVariants.Saye} and {XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes}; " +
                    $"you have set {Control.CutCellQuadratureType}, so you are notified that you reach into unknown territory; " +
                    $"If you do not know how to remove this exception, you should better return now!");
            }

            //QuadOrder
            int degU = VelocityDegree();
            int quadOrder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if(this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                //See remarks
                quadOrder *= 2;
                quadOrder += 1;
            }

            return quadOrder;
        }


        /// <summary>
        /// Current velocity
        /// </summary>
        public VectorField<XDGField> Velocity {
            get {
                int D = this.GridData.SpatialDimension;
                return new VectorField<XDGField>(this.CurrentState.Fields.Take(D).Select(f => ((XDGField)f)).ToArray());
            }
        }

        /// <summary>
        /// Current Pressure
        /// </summary>
        public XDGField Pressure {
            get {
                int D = this.GridData.SpatialDimension;
                return (XDGField)this.CurrentState.Fields[D];
            }
        }



        /// <summary>
        /// Usually, the term "DG order of the calculation" means the velocity degree.
        /// </summary>
        protected int VelocityDegree() {
            int pVel;
            if(this.Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v)) {
                pVel = v.Degree;
            } else if(this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1)) {
                pVel = v1.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pVel;
        }

        

        private IncompressibleMultiphaseBoundaryCondMap m_boundaryMap;

        /// <summary>
        /// Relation between 
        /// - edge tags (<see cref="Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>, passed to equation components via <see cref="BoSSS.Foundation.CommonParams.EdgeTag"/>)
        /// - boundary conditions specified in the control object (<see cref="AppControl.BoundaryValues"/>)
        /// </summary>
        protected IncompressibleMultiphaseBoundaryCondMap boundaryMap {
            get {
                if(m_boundaryMap == null)
                    m_boundaryMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, new string[] { "A", "B" });
                return m_boundaryMap;
            }
        }

        /// <summary>
        /// dirty hack...
        /// </summary>
        protected override IncompressibleBoundaryCondMap GetBcMap() {
            return boundaryMap;
        }


        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            int D = GridData.SpatialDimension;
            
            if(iLevSet == 0) {
                // +++++++++++++++++++
                // the fluid interface 
                // +++++++++++++++++++

                // averaging at interface:
                ILevelSetParameter levelSetVelocity = new LevelSetVelocity(VariableNames.LevelSetCG, D, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                return levelSetVelocity;
            } else if(iLevSet == 1) {
                // +++++++++++++++++++++
                // the immersed boundary
                // +++++++++++++++++++++

                string[] VelocityNames = VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(1), VariableNames.VelocityVector(D)).ToArray();
                ScalarFunctionTimeDep[] VelFuncs = new ScalarFunctionTimeDep[D];
                for(int d = 0; d < D; d++) {
                    Control.InitialValues_EvaluatorsVec.TryGetValue(VelocityNames[d], out VelFuncs[d]);
                }


                ILevelSetParameter levelSetVelocity = new ExplicitLevelSetVelocity(VariableNames.LevelSetCGidx(1), VelFuncs);
                return levelSetVelocity;
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        protected override int NoOfLevelSets {
            get {
                if(Control.UseImmersedBoundary)
                    return 2;
                else 
                    return 1;
            }
        }

        /// <summary>
        /// Either fluids A and B; or A, B and solid C.
        /// </summary>
        protected override Array SpeciesTable {
            get {
                if(Control.UseImmersedBoundary) {
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

        /// Cell-performance classes:
        /// cell performance class equals number of species present in that cell
        /// </summary>
        protected override void GetCellPerformanceClasses(out int NoOfClasses, out int[] CellPerfomanceClasses, int TimeStepNo, double physTime) {
            (NoOfClasses,CellPerfomanceClasses)=CellClassifier.ClassifyCells(this,this.Control.DynamicLoadbalancing_ClassifierType);
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            int pVel = VelocityDegree();
            int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
            int D = this.GridData.SpatialDimension;


            if (this.Control.UseSchurBlockPrec) {
                // using a Schur complement for velocity & pressure
                var confMomConti = new MultigridOperator.ChangeOfBasisConfig();
                for (int d = 0; d < D; d++) {
                    d.AddToArray(ref confMomConti.VarIndex);
                    //Math.Max(1, pVel - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                    pVel.AddToArray(ref confMomConti.DegreeS);
                }
                D.AddToArray(ref confMomConti.VarIndex);
                //Math.Max(0, pPrs - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                pPrs.AddToArray(ref confMomConti.DegreeS);

                confMomConti.mode = MultigridOperator.Mode.SchurComplement;

                configsLevel.Add(confMomConti);
            } else {
                // configurations for velocity
                for (int d = 0; d < D; d++) {
                    var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pVel },
                        //DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.VelocityVector(D)[d]) }
                    };
                    configsLevel.Add(configVel_d);
                }
                // configuration for pressure
                var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pPrs },
                    //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                    mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Pressure) }
                };
                configsLevel.Add(configPres);
            }
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {

            OperatorFactory opFactory = new OperatorFactory();

            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            FinalOperatorSettings(XOP);
            XOP.Commit();

            return XOP;
        }

        /// <summary>
        /// Misc adjustments to the spatial operator before calling <see cref="ISpatialOperator.Commit"/>
        /// </summary>
        /// <param name="XOP"></param>
        protected virtual void FinalOperatorSettings(XSpatialOperatorMk2 XOP) {
            XOP.FreeMeanValue[VariableNames.Pressure] = !GetBcMap().DirichletPressureBoundary;
            XOP.IsLinear = !(this.Control.PhysicalParameters.IncludeConvection || Control.NonlinearCouplingSolidFluid);
            XOP.LinearizationHint = XOP.IsLinear == true ? LinearizationHint.AdHoc : this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard ? LinearizationHint.AdHoc : LinearizationHint.GetJacobiOperator;
            XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;
        }

        /// <summary>
        /// Setup of the incompressible two-phase Navier-Stokes equation
        /// </summary>
        protected virtual void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            int quadOrder = QuadOrder();
            GetBcMap();

            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);

            // === momentum equations === //
            for (int d = 0; d < D; ++d) {
                DefineMomentumEquation(opFactory, config, d, D);

                // Add Gravitation
                if(config.isGravity) {
                    var GravA = Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A, Control.GetGravity("A", d));
                    opFactory.AddParameter(GravA);
                    var GravB = Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B, Control.GetGravity("B", d));
                    opFactory.AddParameter(GravB);
                }

                // Add additional volume forces
                if(config.isVolForce) {
                    var VolForceA = VolumeForce.CreateFrom("A", d, D, Control, Control.GetVolumeForce("A", d));
                    opFactory.AddParameter(VolForceA);
                    var VolForceB = VolumeForce.CreateFrom("B", d, D, Control, Control.GetVolumeForce("B", d));
                    opFactory.AddParameter(VolForceB);
                }
            }

            // === continuity equation === //
            if (config.isContinuity) {
                DefineContinuityEquation(opFactory, config, D);
            }

            // === additional parameters === //
            opFactory.AddCoefficient(new SlipLengths(config, VelocityDegree()));
            Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
            if (((config.physParams.IncludeConvection && config.isTransport) | (config.thermParams.IncludeConvection )) & this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                opFactory.AddParameter(new Velocity0(D));
                opFactory.AddParameter(v0Mean);
            }

            // === level set related parameters === //
            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);            

            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, v0Mean);
            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, normalsParameter);
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

            if (Control.UseImmersedBoundary)
                DefineSystemImmersedBoundary(D, opFactory, lsUpdater);
        }

        /// <summary>
        /// Override this method to customize the assembly of the momentum equation
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="d">Momentum component index</param>
        /// <param name="D">Spatial dimension (2 or 3)</param>
        virtual protected void DefineMomentumEquation(OperatorFactory opFactory, XNSFE_OperatorConfiguration config, int d, int D) {

            // === linearized or parameter free variants, difference only in convective term === //
            if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                opFactory.AddEquation(new NavierStokes("A", d, D, boundaryMap, config));
                opFactory.AddEquation(new NavierStokes("B", d, D, boundaryMap, config));
                opFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, config, config.isMovingMesh));
            } else if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                opFactory.AddEquation(new NavierStokes_Newton("A", d, D, boundaryMap, config));
                opFactory.AddEquation(new NavierStokes_Newton("B", d, D, boundaryMap, config));
                opFactory.AddEquation(new NSEInterface_Newton("A", "B", d, D, boundaryMap, config, config.isMovingMesh));
            } else {
                throw new NotSupportedException();
            }
            opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, config));
        }

        /// <summary>
        /// Override this method to customize the assembly of the continuity equation
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="D">Spatial dimension (2 or 3)</param>
        virtual protected void DefineContinuityEquation(OperatorFactory opFactory, XNSFE_OperatorConfiguration config, int D) {
            opFactory.AddEquation(new Continuity("A", config, D, boundaryMap));
            opFactory.AddEquation(new Continuity("B", config, D, boundaryMap));
            opFactory.AddEquation(new InterfaceContinuity("A", "B", config, D, config.isMatInt));
        }

        /// <summary>
        /// Definition of the boundary condition on the immersed boundary, <see cref="XNSE_Control.UseImmersedBoundary"/>;
        /// Override to customize.
        /// </summary>
        protected virtual void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);

            if (this.Control.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature) BoSSS.Foundation.XDG.Quadrature.BruteForceSettingsOverride.doubleCutCellOverride = true;

            for (int d = 0; d < D; ++d) {
                // so far only no slip!
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddEquation(new NSEimmersedBoundary("A", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                    opFactory.AddEquation(new NSEimmersedBoundary("B", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                } else if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    opFactory.AddEquation(new NSEimmersedBoundary_Newton("A", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                    opFactory.AddEquation(new NSEimmersedBoundary_Newton("B", "C", 1, d, D, boundaryMap, config, config.isMovingMesh));
                }

                // surface tension on IBM
                if (config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine && config.physParams.Sigma != 0.0) {
                    opFactory.AddEquation(new NSEimmersedBoundary_SurfaceTension("A", "B", d, D, 1));
                }

                // GNBC
                if (config.dntParams.IBM_BoundaryType != IBM_BoundaryType.NoSlip) {
                    opFactory.AddEquation(new NSEimmersedBoundary_GNBC("A", "B", d, D, config.getPhysParams, 1));
                }
            }

            opFactory.AddEquation(new ImmersedBoundaryContinuity("A", "C", 1, config, D));
            opFactory.AddEquation(new ImmersedBoundaryContinuity("B", "C", 1, config, D));

            //throw new NotImplementedException("todo");
            opFactory.AddParameter((ParameterS)GetLevelSetVelocity(1));

            if (config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine || config.dntParams.IBM_BoundaryType != IBM_BoundaryType.NoSlip) {
                var normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[1]).Basis.Degree, VariableNames.LevelSetCGidx(1));
                opFactory.AddParameter(normalsParameter);
                lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCGidx(1), normalsParameter);
            }
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //// check some properties (only for debugging)
            //double r0 = 1;
            //double volumeRef = (1.0 / 3.0) * Math.PI * Math.Pow(r0, 3); // quarter domain
            ////double rCalc = 0.9074;
            ////double volumeCalc = (1.0 / 3.0) * Math.PI * Math.Pow(rCalc, 3); // quarter domain
            //int quadOrder = QuadOrder();
            //Console.WriteLine("quadOrder = {0}", quadOrder);
            //double volume = XNSEUtils.GetSpeciesArea(LsTrk, LsTrk.GetSpeciesId("A"), quadOrder);
            //Console.WriteLine("droplet volume: volume_A = {0} (ref volume = {1}; {2})", volume, volumeRef, 100*(volume - volumeRef)/volumeRef);
            ////Console.WriteLine("droplet volume: volume_A = {0} (calc volume = {1}; {2})", volume, volumeCalc, 100 * (volume - volumeCalc) / volumeCalc);
            // Update Calls
            dt = GetTimestep();
            Console.WriteLine($"Starting time step {TimestepNo}, dt = {dt} ...");
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            Console.WriteLine($"Done with time step {TimestepNo}.");
            return dt;
        }
    }
}