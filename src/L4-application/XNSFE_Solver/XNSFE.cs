using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
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
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;

namespace BoSSS.Application.XNSFE_Solver {

    /// <summary>
    /// Extension of the <see cref="XNSE"/>-solver for additional heat transfer.
    /// (The 'F' stands for Fourier equation, i.e. Heat equation.)
    /// Changed to Newton Solver 4/2021, Picard might give unexpected results - MR
    /// </summary>
    public class XNSFE : XNSFE<XNSFE_Control> {

        // ===========
        //  Main file
        // ===========
        static void Main(string[] args) {
            //InitMPI(args);
           // ilPSP.Environment.InitThreading(true, 8);
            //BoSSS.Application.XNSFE_Solver.Tests.ASUnitTest.InterfaceSlipTestLin(3, 0.0d, ViscosityMode.FullySymmetric, 0.0d, XQuadFactoryHelper.MomentFittingVariants.Saye, NonLinearSolverCode.Newton, 1.0d, 1.0d, 1.2d);
            //Assert.IsTrue(false, "remove me");

            //InitMPI();
            //DeleteOldPlotFiles();
            //Tests.ParameterizedLevelSetTest_Elemental.Test();
            //Tests.ParameterizedLevelSet_Translation.Test();
            //Tests.ASUnitTest.ParameterizedLevelSetTest_Translation();
            //BoSSS.Application.XNSFE_Solver.Tests.ASUnitTest.TransientEvaporationTest(0.0, 3, 0.1, XQuadFactoryHelper.MomentFittingVariants.Saye, SurfaceStressTensor_IsotropicMode.Curvature_Projected, NonLinearSolverCode.Newton, Solution.XdgTimestepping.LevelSetHandling.LieSplitting);
            //BoSSS.Application.XNSFE_Solver.Tests.ASUnitTest.ParameterizedLevelSetTest(2);
            //System.Environment.Exit(111);



            XNSFE._Main(args, false, delegate () {
                var p = new XNSFE();
                return p;
            });
        }
    }

    /// <summary>
    /// Extension of the <see cref="XNSE"/>-solver for additional heat transfer.
    /// (The 'F' stands for Fourier equation, i.e. Heat equation.)
    /// Changed to Newton Solver 4/2021, Picard might give unexpected results - MR
    /// </summary>
    public class XNSFE<T> : XNSE<T> where T : XNSFE_Control, new() {
       
        private void AddXHeatMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            int D = this.GridData.SpatialDimension;

            int pTemp = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Temperature].Degree;
            // configuration for Temperature
            var confTemp = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pTemp }, //Math.Max(1, pTemp - iLevel) },
                mode = MultigridOperator.Mode.LeftInverse_DiagBlock,//MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Temperature) }
            };
            configsLevel.Add(confTemp);

            // configuration for auxiliary heat flux
            if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                int pFlux;
                if (this.Control.FieldOptions.TryGetValue("HeatFlux*", out FieldOpts f)) {
                    pFlux = f.Degree;
                } else if (this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.HeatFluxX, out FieldOpts f1)) {
                    pFlux = f1.Degree;
                } else {
                    throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of HeatFlux not found");
                }
                for (int d = 0; d < D; d++) {
                    var confHeatFlux = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pFlux }, // Math.Max(1, pFlux - iLevel) },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.HeatFluxVectorComponent(d)) }
                    };
                    configsLevel.Add(confHeatFlux);
                }
            }
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            base.AddMultigridConfigLevel(configsLevel, iLevel);
            AddXHeatMultigridConfigLevel(configsLevel);
        }

        private ThermalMultiphaseBoundaryCondMap m_thermBoundaryMap;

        /// <summary>
        /// Relation between 
        /// - edge tags (<see cref="Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>, passed to equation components via <see cref="BoSSS.Foundation.CommonParams.EdgeTag"/>)
        /// - boundary conditions specified in the control object (<see cref="AppControl.BoundaryValues"/>)
        /// </summary>
        protected ThermalMultiphaseBoundaryCondMap thermBoundaryMap {
            get {
                if (m_thermBoundaryMap == null) {
                    List<string> SpeciesList = new List<string>() { "A", "B" };
                    if (this.Control.UseImmersedBoundary)
                        SpeciesList.Add("C");
                    m_thermBoundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, SpeciesList.ToArray());
                }
                return m_thermBoundaryMap;
            }
        }


        protected override void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystem(D, opFactory, lsUpdater);
            AddXHeat(D, opFactory, lsUpdater);
        }

        /// <summary>
        /// Add Heat and in case of evaporation extensions for momentum and continuity equations
        /// <see cref="DefineHeatEquation(OperatorFactory, XNSFE_OperatorConfiguration, int)"/>
        /// <see cref="DefineMomentumEquation(OperatorFactory, XNSFE_OperatorConfiguration, int, int)"/>
        /// <see cref="DefineContinuityEquation(OperatorFactory, XNSFE_OperatorConfiguration, int)"/>
        /// </summary>
        /// <param name="D"></param>
        /// <param name="opFactory"></param>
        /// <param name="lsUpdater"></param>
        private void AddXHeat(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {            
            int quadOrder = QuadOrder();
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            ThermalMultiphaseBoundaryCondMap heatBoundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            // === heat equation === //
            DefineHeatEquation(opFactory, config, D);

            // === additional parameters === //
            if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) { 
                opFactory.AddParameter(new Velocity0(D)); 
            }

            // Add HeatSource
            if (config.isHeatSource) {
                var HeatA = HeatSource.CreateFrom("A", Control, Control.GetHeatSource("A"));
                opFactory.AddParameter(HeatA);
                var HeatB = HeatSource.CreateFrom("B", Control, Control.GetHeatSource("B"));
                opFactory.AddParameter(HeatB);
            }

            if (config.isEvaporation) {
                //Console.WriteLine("Including mass transfer.");

                if (this.Control.NonLinearSolver.SolverCode != NonLinearSolverCode.Newton) throw new ApplicationException("Evaporation only implemented with use of Newton-solver!");

                var MassFluxExt = new MassFluxExtension_Evaporation(config);
                lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, MassFluxExt);
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddParameter(MassFluxExt);
                }
            }

            opFactory.AddCoefficient(new EvapMicroRegion());
            opFactory.AddCoefficient(new SlipLengths(config));

            if (config.prescribedMassflux != null)
                opFactory.AddCoefficient(new PrescribedMassFlux(config));

            // When using LDG Formulation
            //opFactory.AddParameter(new Temperature0());
            //opFactory.AddParameter(new HeatFlux0(D, lsUpdater.Tracker, config)); // also for test reasons        
        }

        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            if(iLevSet == 0) {
                // Main Difference to base implementation:
                //var levelSetVelocity = new LevelSetVelocityEvaporative("Phi", GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters, new XNSFE_OperatorConfiguration(Control);

                var config = new XNSFE_OperatorConfiguration(Control);
                ILevelSetParameter levelSetVelocity;

                if (config.isEvaporation) {
                    levelSetVelocity = new LevelSetVelocityGeneralNonMaterial(VariableNames.LevelSetCG, GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters, config);
                } else {
                    levelSetVelocity = new LevelSetVelocity(VariableNames.LevelSetCG, GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters); 
                }

                return levelSetVelocity;
            } else {
                return base.GetLevelSetVelocity(iLevSet);
            }
        }

        /// <summary>
        /// override of <see cref="DefineMomentumEquation(OperatorFactory, XNSE_OperatorConfiguration, int, int)"/>
        /// adding evaporation extension for Navier-Stokes equations
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        protected override void DefineMomentumEquation(OperatorFactory opFactory, XNSE_OperatorConfiguration config, int d, int D) {
            base.DefineMomentumEquation(opFactory, config, d, D);

            var extendedConfig = new XNSFE_OperatorConfiguration(this.Control);
            // === evaporation extension === //
            if (extendedConfig.isEvaporation) {
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddEquation(new InterfaceNSE_Evaporation("A", "B", D, d, extendedConfig));                    
                } else if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    opFactory.AddEquation(new InterfaceNSE_Evaporation_Newton("A", "B", D, d, extendedConfig));
                }                    
            }
            if (extendedConfig.isBuoyancy && config.isGravity) {
                opFactory.AddEquation(new NavierStokesBuoyancy("A", D, d, extendedConfig));
                opFactory.AddEquation(new NavierStokesBuoyancy("B", D, d, extendedConfig));
            }
        }


        protected override void DefineContinuityEquation(OperatorFactory opFactory, XNSE_OperatorConfiguration config, int D) {
            base.DefineContinuityEquation(opFactory, config, D);

            var extendedConfig = new XNSFE_OperatorConfiguration(this.Control);
            // === evaporation extension === //
            if (extendedConfig.isEvaporation) {
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddEquation(new InterfaceContinuity_Evaporation("A", "B", D, extendedConfig));
                } else if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    opFactory.AddEquation(new InterfaceContinuity_Evaporation_Newton("A", "B", D, extendedConfig));
                }
            }
        }

        protected override void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystemImmersedBoundary(D, opFactory, lsUpdater);

            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);

            // Add HeatSource
            if (config.isHeatSource) {                
                var HeatC = HeatSource.CreateFrom("C", Control, Control.GetHeatSource("C"));
                opFactory.AddParameter(HeatC);
            }

            opFactory.AddEquation(new SolidHeat("C", D, thermBoundaryMap, config));
            opFactory.AddEquation(new ImmersedBoundaryHeat("A", "C", 1, D, config));
            opFactory.AddEquation(new ImmersedBoundaryHeat("B", "C", 1, D, config));

            // we need these "dummy" equations, otherwise the matrix has zero rows/columns
            // If this is to be used frequently something more sophisticated should be implemented, e.g. strike out rows for unused variables...
            for(int d = 0; d < D; d++) {
                opFactory.AddEquation(new ImmersedBoundaryDummyMomentum("C", d));
            }
            opFactory.AddEquation(new ImmersedBoundaryDummyConti("C"));

        }

        /// <summary>
        /// Override this method to customize the assembly of the heat equation
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="d">Momentum component index</param>
        /// <param name="D">Spatial dimension (2 or 3)</param>
        virtual protected void DefineHeatEquation(OperatorFactory opFactory, XNSFE_OperatorConfiguration config, int D) {

            // === linearized or parameter free variants, difference only in convective term === //
            if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                opFactory.AddEquation(new Heat("A", D, thermBoundaryMap, config));
                opFactory.AddEquation(new Heat("B", D, thermBoundaryMap, config));

                if (config.isEvaporation) {
                    opFactory.AddEquation(new HeatInterface_Evaporation("A", "B", D, thermBoundaryMap, config));
                } else {
                    opFactory.AddEquation(new HeatInterface("A", "B", D, thermBoundaryMap, config));
                }

                if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                    for (int d = 0; d < D; ++d) {
                        throw new ApplicationException("Warning using LDG Formulation for Heat, this is untested. Remove this statement only if you now what you are doing!");
                        opFactory.AddEquation(new HeatFlux("A", d, D, thermBoundaryMap, config));
                        opFactory.AddEquation(new HeatFlux("B", d, D, thermBoundaryMap, config));
                        opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, thermBoundaryMap, config));
                    }
                }                
            } else if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                opFactory.AddEquation(new Heat_Newton("A", D, thermBoundaryMap, config));
                opFactory.AddEquation(new Heat_Newton("B", D, thermBoundaryMap, config));

                if (config.isEvaporation) {
                    opFactory.AddEquation(new HeatInterface_Evaporation_Newton("A", "B", D, thermBoundaryMap, config));
                } else {
                    opFactory.AddEquation(new HeatInterface_Newton("A", "B", D, thermBoundaryMap, config));
                }

                if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                    for (int d = 0; d < D; ++d) {
                        opFactory.AddEquation(new HeatFlux("A", d, D, thermBoundaryMap, config));
                        opFactory.AddEquation(new HeatFlux("B", d, D, thermBoundaryMap, config));
                        opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, thermBoundaryMap, config));
                    }
                    throw new ApplicationException("Warning using LDG Formulation for Heat, this is untested. Remove this statement only if you now what you are doing!");
                }
            } else {
                throw new NotSupportedException();
            }
        }

        protected override XDifferentialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {

            OperatorFactory opFactory = new OperatorFactory();

            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator            
            XDifferentialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());
            var config = new XNSFE_OperatorConfiguration(this.Control);
            //final settings
            XOP.FreeMeanValue[VariableNames.Pressure] = !GetBcMap().DirichletPressureBoundary;
            XOP.IsLinear = !(this.Control.PhysicalParameters.IncludeConvection || this.Control.ThermalParameters.IncludeConvection || (config.isEvaporation & Control.IncludeRecoilPressure) ||Control.NonlinearCouplingSolidFluid); // when using Newton solver the employed coupling between heat and momentum is nonlinear
            XOP.LinearizationHint = XOP.IsLinear == true ? LinearizationHint.AdHoc : this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard ? LinearizationHint.AdHoc : LinearizationHint.GetJacobiOperator;

            XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;

            if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton & this.Control.ThermalParameters.IncludeConvection) {
                //XOP.HomotopyUpdate.Add(delegate (double HomotopyScalar) {
                //    if (HomotopyScalar < 0.0)
                //        throw new ArgumentOutOfRangeException();
                //    if (HomotopyScalar > 1.0)
                //        throw new ArgumentOutOfRangeException();

                //    //Console.WriteLine("Updating Homotopy Scalar aka. scaling for heat convection, Value = " + HomotopyScalar);
                //});
            }
            XOP.FluxesAreNOTMultithreadSafe = true;
            XOP.Commit();

            if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) Console.WriteLine("Warning using Picard iteration, this is not recommended!");
            PrintConfiguration();
            return XOP;
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            //if (Control.InitialValues_EvaluatorsVec.TryGetValue("Temperature#B", out var scalarFunctionTimeDep) && this.Control.SkipSolveAndEvaluateResidual) {
            //    ScalarFunction T_ex = null;
            //    T_ex = scalarFunctionTimeDep.SetTime(phystime);
            //    ((XDGField)this.CurrentState.Fields.Single(s => s.Identification == "Temperature")).GetSpeciesShadowField("B").ProjectField(T_ex);
            //}            

            // Set timestep as minimum of capillary timestep restriction or level set CFL
            //SetTimestep();
            return base.RunSolverOneStep(TimestepNo, phystime, dt);
        }

        protected override List<DGField> GetInterfaceVelocity() {
            int D = this.GridData.SpatialDimension;
            var cm = this.LsTrk.Regions.GetCutCellMask4LevSet(0);
            var ThermParam = this.Control.ThermalParameters;
            var gdat = this.GridData;
            var Tracker = this.LsTrk;

            VectorField<XDGField> Velocity = new VectorField<XDGField>(D.ForLoop(d => (XDGField)m_RegisteredFields.SingleOrDefault(s => s.Identification == VariableNames.Velocity_d(d))));
            XDGField Temperature = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == VariableNames.Temperature).SingleOrDefault();
            var LevelSet = this.m_RegisteredFields.Where(s => s.Identification == VariableNames.LevelSetCG).SingleOrDefault();
            int p = LevelSet.Basis.Degree;

            Basis b = new Basis(gdat, 2 * p);
            XDGBasis xb = new XDGBasis(Tracker, 2 * p);

            VectorField<DGField> Normal = new VectorField<DGField>(D.ForLoop(d => m_RegisteredFields.SingleOrDefault(s => s.Identification == VariableNames.NormalVectorComponent(d)))).CloneAs();
            // this normal is the direct gradient of phi, we normalize first, but only in cut cells!!
            {
                Normal.Clear();
                Normal.Gradient(1.0, LevelSet);
                SinglePhaseField Normalizer = new SinglePhaseField(b);
                for (int d = 0; d < D; d++) {
                    Normalizer.ProjectPow(1.0, Normal[d], 2.0, cm);
                }
                var NormalizerTemp = Normalizer.CloneAs();
                Normalizer.Clear();
                Normalizer.ProjectPow(1.0, NormalizerTemp, -0.5, cm);

                Normal.ScalePointwise(1.0, Normalizer);
            }

            // Heat Flux
            List<XDGField> HeatFlux = new List<XDGField>();
            for (int d = 0; d < D; d++) {
                var xField = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == VariableNames.HeatFluxVectorComponent(d)).SingleOrDefault();
                if (xField == null) {
                    xField = new XDGField(xb, VariableNames.HeatFluxVectorComponent(d));
                    this.RegisterField(xField);
                }
                HeatFlux.Add(xField);
            }

            for (int d = 0; d < D; d++) {
                HeatFlux[d].Clear();
                HeatFlux[d].Derivative(1.0, Temperature, d);
                ((XDGField)HeatFlux[d]).GetSpeciesShadowField("A").Scale(-ThermParam.k_A);
                ((XDGField)HeatFlux[d]).GetSpeciesShadowField("B").Scale(-ThermParam.k_B);
            }

            // Mass Flux
            var MassFlux = this.m_RegisteredFields.Where(s => s.Identification == "MassFlux").SingleOrDefault();
            if (MassFlux == null) {
                MassFlux = new SinglePhaseField(b, "MassFlux");
                this.RegisterField(MassFlux);
            }

            MassFlux.Clear();
            for (int d = 0; d < D; d++) {
                MassFlux.ProjectProduct(1.0 / ThermParam.hVap, Normal[d], HeatFlux[d].GetSpeciesShadowField("A"), cm);
                MassFlux.ProjectProduct(-1.0 / ThermParam.hVap, Normal[d], HeatFlux[d].GetSpeciesShadowField("B"), cm);
            }


            List<XDGField> xInterfaceVelocity = new List<XDGField>(); // = this.m_RegisteredFields.Where(s => s.Identification.Contains("IntefaceVelocityXDG")).ToList();
            List<DGField> InterfaceVelocity = new List<DGField>(); // = this.m_RegisteredFields.Where(s => s.Identification.Contains("IntefaceVelocityDG")).ToList();
            for (int d = 0; d < D; d++) {
                var xField = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == $"IntefaceVelocityXDG_{d}").SingleOrDefault();
                if (xField == null) {
                    xField = new XDGField(xb, $"IntefaceVelocityXDG_{d}");
                    this.RegisterField(xField);
                }
                xInterfaceVelocity.Add(xField);

                var Field = (DGField)this.m_RegisteredFields.Where(s => s.Identification == $"IntefaceVelocityDG_{d}").SingleOrDefault();
                if (Field == null) {
                    Field = new SinglePhaseField(b, $"IntefaceVelocityDG_{d}");
                    this.RegisterField(Field);
                }
                InterfaceVelocity.Add(Field);
            }
            
            for (int d = 0; d < D; d++) {
                var xField = xInterfaceVelocity[d];
                xField.Clear();
                xField.AccLaidBack(1.0, Velocity[d], cm);
                xField.GetSpeciesShadowField("A").ProjectProduct(-1.0 / ThermParam.rho_A, MassFlux, Normal[d], cm);
                xField.GetSpeciesShadowField("B").ProjectProduct(-1.0 / ThermParam.rho_B, MassFlux, Normal[d], cm);

                var Field = InterfaceVelocity[d];
                Field.Clear();
                Field.AccLaidBack(ThermParam.rho_A, xField.GetSpeciesShadowField("A"), cm);
                Field.AccLaidBack(ThermParam.rho_B, xField.GetSpeciesShadowField("B"), cm);
                Field.Scale(1.0 / (ThermParam.rho_A + ThermParam.rho_B), cm);
            }

            // Project normal
            //var NormalInterfaceVelocity = InterfaceVelocity.Select(s => s.CloneAs()).ToArray();
            //NormalInterfaceVelocity.ForEach(s => s.Clear());
            //for (int i = 0; i < D; i++) {
            //    var temp = NormalInterfaceVelocity[i].CloneAs();
            //    temp.Clear();
            //    for (int j = 0; j < D; j++) {                    
            //        temp.ProjectProduct(1.0, Normal[j], InterfaceVelocity[i], cm, true);
            //    }
            //    NormalInterfaceVelocity[i].ProjectProduct(1.0, temp, Normal[i], cm, false);
            //}
            //for (int i = 0; i < D; i++) {
            //    InterfaceVelocity[i].Clear();
            //    InterfaceVelocity[i].Acc(1.0, NormalInterfaceVelocity[i]);
            //}
            return InterfaceVelocity;
        }

        private void SetTimestep() {
            if (this.Control.PhysicalParameters.Sigma != 0.0) {
                double dt_cfl, dt_cap;
                dt_cap = GetCapillaryTimestep();
                //extend by levelset cfl

                if (Math.Abs(dt_cap - this.Control.dtFixed) > 1e-3 * Math.Abs(dt_cap)) {
                    Console.WriteLine($"Setting new capillary timestep restriction, new timestep size: {dt_cap}");
                    this.Control.dtFixed = dt_cap;
                }
            }
        }

        private double GetCapillaryTimestep() {
            var C = this.Control;
            double h = this.GridData.iGeomCells.h_min.Min();
            int p = VelocityDegree();
            double safety = 5;
            return 1 / safety * Math.Sqrt((C.PhysicalParameters.rho_A + C.PhysicalParameters.rho_B) * Math.Pow(h / (p + 1), 3) / (2 * Math.PI * Math.Abs(C.PhysicalParameters.Sigma)));
        }

        private void PlotAdditionalFields(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            #region additional fields

            var ThermParam = this.Control.ThermalParameters;
            var PhysParam = this.Control.PhysicalParameters;
            int D = this.GridData.SpatialDimension;
            var gdat = this.GridData;


            XDGField Temperature = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == VariableNames.Temperature).SingleOrDefault();
            var Tracker = Temperature.Basis.Tracker;
            var cm = Tracker.Regions.GetCutCellMask();
            VectorField<XDGField> Velocity = new VectorField<XDGField>(D.ForLoop(d => (XDGField)m_RegisteredFields.SingleOrDefault(s => s.Identification == VariableNames.Velocity_d(d))));
            VectorField<XDGField>[] VelocityGrad = new VectorField<XDGField>[D];
            for (int d = 0; d < D; d++) {
                VelocityGrad[d] = new VectorField<XDGField>(D.ForLoop(_d => new XDGField(Velocity[d].Basis, VariableNames.Velocity_GradientVector(D)[d, _d])));
                VelocityGrad[d].Gradient(1.0, Velocity[d]);
            }
            XDGField Pressure = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == VariableNames.Pressure).SingleOrDefault();
            var LevelSet = this.m_RegisteredFields.Where(s => s.Identification == VariableNames.LevelSetCG).SingleOrDefault();
            int p = LevelSet.Basis.Degree;
            VectorField<DGField> Normal = new VectorField<DGField>(D.ForLoop(d => m_RegisteredFields.SingleOrDefault(s => s.Identification == VariableNames.NormalVectorComponent(d)))).CloneAs();

            Basis b = new Basis(gdat, 2 * p);
            XDGBasis xb = new XDGBasis(Tracker, 2 * p);
            List<DGField> TrueNormal = new List<DGField>();
            for (int d = 0; d < D; d++) {
                var Field = this.m_RegisteredFields.Where(s => s.Identification == $"TrueNormal_{d}").SingleOrDefault();
                if (Field == null) {
                    Field = new SinglePhaseField(b, $"TrueNormal_{d}");
                    this.RegisterField(Field);
                }
                TrueNormal.Add(Field);
            }
            // this normal is the direct gradient of phi, we normalize first, but only in cut cells!!
            {
                Normal.Clear();
                Normal.Gradient(1.0, LevelSet);
                SinglePhaseField Normalizer = new SinglePhaseField(b);
                for (int d = 0; d < D; d++) {
                    Normalizer.ProjectPow(1.0, Normal[d], 2.0, cm);
                }
                var NormalizerTemp = Normalizer.CloneAs();
                Normalizer.Clear();
                Normalizer.ProjectPow(1.0, NormalizerTemp, -0.5, cm);

                Normal.ScalePointwise(1.0, Normalizer);
                for (int d = 0; d < D; d++) {
                    TrueNormal[d].Clear();
                    TrueNormal[d].AccLaidBack(1.0, Normal[d], cm);
                }
            }
            DGField Curvature = this.m_RegisteredFields.Where(s => s.Identification == VariableNames.Curvature).SingleOrDefault();

            // Heat Flux
            List<XDGField> HeatFlux = new List<XDGField>();
            for (int d = 0; d < D; d++) {
                var xField = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == VariableNames.HeatFluxVectorComponent(d)).SingleOrDefault();
                if (xField == null) {
                    xField = new XDGField(xb, VariableNames.HeatFluxVectorComponent(d));
                    this.RegisterField(xField);
                }
                HeatFlux.Add(xField);
            }

            for (int d = 0; d < D; d++) {
                HeatFlux[d].Clear();
                HeatFlux[d].Derivative(1.0, Temperature, d);
                ((XDGField)HeatFlux[d]).GetSpeciesShadowField("A").Scale(-ThermParam.k_A);
                ((XDGField)HeatFlux[d]).GetSpeciesShadowField("B").Scale(-ThermParam.k_B);
            }

            // Mass Flux
            var MassFlux = this.m_RegisteredFields.Where(s => s.Identification == "MassFlux").SingleOrDefault();
            if (MassFlux == null) {
                MassFlux = new SinglePhaseField(b, "MassFlux");
                this.RegisterField(MassFlux);
            }

            MassFlux.Clear();
            for (int d = 0; d < D; d++) {
                MassFlux.ProjectProduct(1.0 / ThermParam.hVap, Normal[d], HeatFlux[d].GetSpeciesShadowField("A"), cm);
                MassFlux.ProjectProduct(-1.0 / ThermParam.hVap, Normal[d], HeatFlux[d].GetSpeciesShadowField("B"), cm);
            }

            // Interface velocity
            List<XDGField> xInterfaceVelocity = new List<XDGField>(); // = this.m_RegisteredFields.Where(s => s.Identification.Contains("IntefaceVelocityXDG")).ToList();
            List<DGField> InterfaceVelocity = new List<DGField>(); // = this.m_RegisteredFields.Where(s => s.Identification.Contains("IntefaceVelocityDG")).ToList();
            for (int d = 0; d < D; d++) {
                var xField = (XDGField)this.m_RegisteredFields.Where(s => s.Identification == $"IntefaceVelocityXDG_{d}").SingleOrDefault();
                if (xField == null) {
                    xField = new XDGField(xb, $"IntefaceVelocityXDG_{d}");
                    this.RegisterField(xField);
                }
                xInterfaceVelocity.Add(xField);

                var Field = (DGField)this.m_RegisteredFields.Where(s => s.Identification == $"IntefaceVelocityDG_{d}").SingleOrDefault();
                if (Field == null) {
                    Field = new SinglePhaseField(b, $"IntefaceVelocityDG_{d}");
                    this.RegisterField(Field);
                }
                InterfaceVelocity.Add(Field);
            }

            for (int d = 0; d < D; d++) {
                var xField = xInterfaceVelocity[d];
                xField.Clear();
                xField.AccLaidBack(1.0, Velocity[d], cm);
                xField.GetSpeciesShadowField("A").ProjectProduct(-1.0 / ThermParam.rho_A, MassFlux, Normal[d], cm);
                xField.GetSpeciesShadowField("B").ProjectProduct(-1.0 / ThermParam.rho_B, MassFlux, Normal[d], cm);

                var Field = InterfaceVelocity[d];
                Field.Clear();
                Field.AccLaidBack(ThermParam.rho_A, xField.GetSpeciesShadowField("A"), cm);
                Field.AccLaidBack(ThermParam.rho_B, xField.GetSpeciesShadowField("B"), cm);
                Field.Scale(1.0 / (ThermParam.rho_A + ThermParam.rho_B), cm);
            }

            //// Conti
            //var ContiRes = this.m_RegisteredFields.Where(s => s.Identification == "ContinuityResidual").SingleOrDefault();
            //if (ContiRes == null) {
            //    ContiRes = new XDGField(xb, "ContinuityResidual");
            //    this.RegisterField(ContiRes);
            //}

            //ContiRes.Clear();
            //for (int d = 0; d < D; d++) {
            //    ContiRes.Derivative(1.0, Velocity[d], d);
            //}

            //// Momentum

            //// Heat

            // Divergence jump
            var DivJump = this.m_RegisteredFields.Where(s => s.Identification == "DivergenceJump").SingleOrDefault();
            if (DivJump == null) {
                DivJump = new SinglePhaseField(b, "DivergenceJump");
                this.RegisterField(DivJump);
            }

            DivJump.Clear();
            for (int d = 0; d < D; d++) {
                DivJump.ProjectProduct(ThermParam.rho_A, Velocity[d].GetSpeciesShadowField("A"), Normal[d], cm);
                DivJump.ProjectProduct(-ThermParam.rho_A, InterfaceVelocity[d], Normal[d], cm);

                DivJump.ProjectProduct(-ThermParam.rho_B, Velocity[d].GetSpeciesShadowField("B"), Normal[d], cm);
                DivJump.ProjectProduct(ThermParam.rho_B, InterfaceVelocity[d], Normal[d], cm);
            }

            // Tangential velocity jump
            var TanJump = this.m_RegisteredFields.Where(s => s.Identification == "TangentialJump").SingleOrDefault();
            if (TanJump == null) {
                TanJump = new SinglePhaseField(b, "TangentialJump");
                this.RegisterField(TanJump);
            }

            TanJump.Clear();
            // Assuming 2D
            TanJump.ProjectProduct(-1.0, Velocity[0].GetSpeciesShadowField("A"), Normal[1], cm);
            TanJump.ProjectProduct(1.0, Velocity[0].GetSpeciesShadowField("B"), Normal[1], cm);
            TanJump.ProjectProduct(1.0, Velocity[1].GetSpeciesShadowField("A"), Normal[0], cm);
            TanJump.ProjectProduct(-1.0, Velocity[1].GetSpeciesShadowField("B"), Normal[0], cm);

            // Momentum jump
            List<DGField> MomentumJump = new List<DGField>();
            for (int d = 0; d < D; d++) {
                var Field = this.m_RegisteredFields.Where(s => s.Identification == $"MomentumJump_{d}").SingleOrDefault();
                if (Field == null) {
                    Field = new SinglePhaseField(b, $"MomentumJump_{d}");
                    this.RegisterField(Field);
                }
                MomentumJump.Add(Field);
            }

            for (int d = 0; d < D; d++) {
                MomentumJump[d].Clear();

                // recoil pressure
                MomentumJump[d].ProjectProduct(-1.0, MassFlux, Velocity[d].GetSpeciesShadowField("A"), cm);
                MomentumJump[d].ProjectProduct(1.0, MassFlux, Velocity[d].GetSpeciesShadowField("B"), cm);

                // pressure jump
                MomentumJump[d].ProjectProduct(-1.0, Pressure.GetSpeciesShadowField("A"), Normal[d], cm);
                MomentumJump[d].ProjectProduct(1.0, Pressure.GetSpeciesShadowField("B"), Normal[d], cm);

                // stress jump
                for (int _d = 0; _d < D; _d++) {
                    MomentumJump[d].ProjectProduct(PhysParam.mu_A, VelocityGrad[d][_d].GetSpeciesShadowField("A"), Normal[_d], cm);
                    MomentumJump[d].ProjectProduct(PhysParam.mu_A, VelocityGrad[_d][d].GetSpeciesShadowField("A"), Normal[_d], cm);

                    MomentumJump[d].ProjectProduct(-PhysParam.mu_B, VelocityGrad[d][_d].GetSpeciesShadowField("B"), Normal[_d], cm);
                    MomentumJump[d].ProjectProduct(-PhysParam.mu_B, VelocityGrad[_d][d].GetSpeciesShadowField("B"), Normal[_d], cm);
                }

                if (Curvature != null) {
                    MomentumJump[d].ProjectProduct(PhysParam.Sigma, Curvature, Normal[d], cm);
                } else {
                    Console.WriteLine("Warning, no curvature detected jump condition vizualization is incorrect!");
                }
            }

            // Energy jump, should be fullfilled, as we use this to evaluate the mass flux
            var EnergyJump = this.m_RegisteredFields.Where(s => s.Identification == "EnergyJump").SingleOrDefault();
            if (EnergyJump == null) {
                EnergyJump = new SinglePhaseField(b, "EnergyJump");
                this.RegisterField(EnergyJump);
            }

            EnergyJump.Clear();
            EnergyJump.AccLaidBack(-ThermParam.hVap, MassFlux);
            for (int d = 0; d < D; d++) {
                EnergyJump.ProjectProduct(1.0, HeatFlux[d].GetSpeciesShadowField("A"), Normal[d]);

                EnergyJump.ProjectProduct(-1.0, HeatFlux[d].GetSpeciesShadowField("B"), Normal[d]);
            }

            // Plot Points on Interface
            var spaceMetrics = Tracker.GetXDGSpaceMetrics(Tracker.SpeciesIdS.ToArray(), this.QuadOrder());
            var surfRule = spaceMetrics.XQuadSchemeHelper.GetLevelSetquadScheme(0, cm);
            var log = new StreamWriter($"JumpConditions_{timestepNo}.csv");
            log.Write("Cell\tNode");
            for (int d = 0; d < D; d++) {
                log.Write($"\tx{d}");
            }
            for (int d = 0; d < D; d++) {
                log.Write($"\tInterfaceVelocity_{d}");
            }
            log.Write($"\tDivergenceJump");
            for (int d = 0; d < D; d++) {
                log.Write($"\tMomentumJump_{d}");
            }
            log.Write($"\tEnergyJump");
            log.Write($"\tTangentialJump");
            log.Write("\n");
            log.Flush();

            double vMax = 0.0;
            foreach (var factory in surfRule.FactoryChain) {
                foreach (var chunkRulePair in factory.RuleFactory.GetQuadRuleSet(cm.ToGeometicalMask(), this.QuadOrder())) {
                    foreach (int cell in chunkRulePair.Chunk.Elements) {
                        QuadRule rule = chunkRulePair.Rule;

                        MultidimensionalArray globalVertices = MultidimensionalArray.Create(
                            1, rule.NoOfNodes, Grid.SpatialDimension);
                        MultidimensionalArray jumpConditions = MultidimensionalArray.Create(
                            1, rule.NoOfNodes, 2 * D + 3);

                        for (int d = 0; d < D; d++) {
                            InterfaceVelocity[d].Evaluate(cell, 1, rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, d));
                            MomentumJump[d].Evaluate(cell, 1, rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, D + 1 + d));
                        }
                        DivJump.Evaluate(cell, 1, rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, D));
                        EnergyJump.Evaluate(cell, 1, rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, 2 * D + 1));
                        TanJump.Evaluate(cell, 1, rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, 2 * D + 2));


                        MultidimensionalArray metrics = Tracker.DataHistories[0].Current.GetLevelSetNormalReferenceToPhysicalMetrics(
                            rule.Nodes, cell, 1);
                        GridData.TransformLocal2Global(rule.Nodes, cell, 1, globalVertices, 0);

                        for (int k = 0; k < rule.NoOfNodes; k++) {
                            double weight = rule.Weights[k];

                            log.Write("{0}\t{1}", cell, k);
                            for (int d = 0; d < D; d++) {
                                log.Write($"\t{globalVertices[0, k, d]}");
                            }
                            double v = 0.0;
                            for (int d = 0; d < D; d++) {
                                log.Write($"\t{jumpConditions[0, k, d]}");
                                v += jumpConditions[0, k, d].Pow2();
                            }
                            v = v.Sqrt();
                            vMax = v > vMax ? v : vMax;
                            log.Write($"\t{jumpConditions[0, k, D]}");
                            for (int d = 0; d < D; d++) {
                                log.Write($"\t{jumpConditions[0, k, D + 1 + d]}");
                            }
                            log.Write($"\t{jumpConditions[0, k, 2 * D + 1]}");
                            log.Write($"\t{jumpConditions[0, k, 2 * D + 2]}");
                            log.Write("\n");
                            log.Flush();

                        }
                    }
                }
            }
            Console.WriteLine("Maximum Absolute interface Velocity: {0}", vMax);

            //// Plot Points on Interface           
            //log = new StreamWriter($"MicroRegion_{timestepNo}.csv");
            //log.Write("Cell\tNode");
            //for (int d = 0; d < D; d++) {
            //    log.Write($"\tx{d}");
            //}            
            //log.Write($"\tHeatFlux");
            //log.Write($"\tWallTemperature");
            //log.Write("\n");
            //log.Flush();
            //{
            //    cm = Tracker.Regions.GetCutCellMask4LevSet(1).Intersect(Tracker.Regions.GetSpeciesMask("A")).Intersect(Tracker.Regions.GetSpeciesMask("C"));
            //    CellQuadratureScheme SurfIntegrationA = spaceMetrics.XQuadSchemeHelper.GetLevelSetquadScheme(1, Tracker.GetSpeciesId("A"), cm);
            //    var rule = SurfIntegrationA.Compile(GridData, this.QuadOrder());
            //    foreach (var chunkRulePair in rule) {
            //        foreach (int cell in chunkRulePair.Chunk.Elements) {
            //            MultidimensionalArray globalVertices = MultidimensionalArray.Create(
            //                1, chunkRulePair.Rule.NoOfNodes, Grid.SpatialDimension);
            //            MultidimensionalArray jumpConditions = MultidimensionalArray.Create(
            //                1, chunkRulePair.Rule.NoOfNodes, 2);

            //            Temperature.Evaluate(cell, 1, chunkRulePair.Rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, 0));
            //            HeatFlux[1].Evaluate(cell, 1, chunkRulePair.Rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, 1));

            //            GridData.TransformLocal2Global(chunkRulePair.Rule.Nodes, cell, 1, globalVertices, 0);

            //            for (int k = 0; k < chunkRulePair.Rule.NoOfNodes; k++) {

            //                log.Write("{0}\t{1}", cell, k);
            //                for (int d = 0; d < D; d++) {
            //                    log.Write($"\t{globalVertices[0, k, d]}");
            //                }
            //                for (int d = 0; d < 2; d++) {
            //                    log.Write($"\t{jumpConditions[0, k, d]}");
            //                }
            //                log.Write("\n");
            //                log.Flush();

            //            }
            //        }
            //    }
            //    cm = Tracker.Regions.GetCutCellMask4LevSet(1).Intersect(Tracker.Regions.GetSpeciesMask("B")).Intersect(Tracker.Regions.GetSpeciesMask("C"));
            //    CellQuadratureScheme SurfIntegrationB = spaceMetrics.XQuadSchemeHelper.GetLevelSetquadScheme(1, Tracker.GetSpeciesId("B"), cm);
            //    rule = SurfIntegrationB.Compile(GridData, this.QuadOrder());
            //    foreach (var chunkRulePair in rule) {
            //        foreach (int cell in chunkRulePair.Chunk.Elements) {
            //            MultidimensionalArray globalVertices = MultidimensionalArray.Create(
            //                1, chunkRulePair.Rule.NoOfNodes, Grid.SpatialDimension);
            //            MultidimensionalArray jumpConditions = MultidimensionalArray.Create(
            //                1, chunkRulePair.Rule.NoOfNodes, 2);

            //            Temperature.Evaluate(cell, 1, chunkRulePair.Rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, 0));
            //            HeatFlux[1].Evaluate(cell, 1, chunkRulePair.Rule.Nodes, jumpConditions.ExtractSubArrayShallow(-1, -1, 1));

            //            GridData.TransformLocal2Global(chunkRulePair.Rule.Nodes, cell, 1, globalVertices, 0);

            //            for (int k = 0; k < chunkRulePair.Rule.NoOfNodes; k++) {

            //                log.Write("{0}\t{1}", cell, k);
            //                for (int d = 0; d < D; d++) {
            //                    log.Write($"\t{globalVertices[0, k, d]}");
            //                }
            //                for (int d = 0; d < 2; d++) {
            //                    log.Write($"\t{jumpConditions[0, k, d]}");
            //                }
            //                log.Write("\n");
            //                log.Flush();

            //            }
            //        }
            //    }
            //}
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {

            //Basis b = new Basis(this.GridData, 0);
            //// Cut Cells 0
            //var CC0 = this.m_RegisteredFields.Where(s => s.Identification == "CC0").SingleOrDefault();
            //if (CC0 == null) {
            //    CC0 = new SinglePhaseField(b, "CC0");
            //    this.RegisterField(CC0);
            //}
            //var ccm0 = this.LsTrk.Regions.GetCutCellMask4LevSet(0);
            //CC0.Clear();
            //CC0.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //    int K = result.GetLength(1); // No nof Nodes
            //    for (int j = 0; j < Len; j++) {
            //        for (int k = 0; k < K; k++) {
            //            if(ccm0.Contains(j0+j))
            //                result[j, k] = 1;
            //        }
            //    }
            //}, new CellQuadratureScheme());

            ////Cut Cells 1
            //var CC1 = this.m_RegisteredFields.Where(s => s.Identification == "CC1").SingleOrDefault();
            //if (CC1 == null) {
            //    CC1 = new SinglePhaseField(b, "CC1");
            //    this.RegisterField(CC1);
            //}
            //var ccm1 = this.LsTrk.Regions.GetCutCellMask4LevSet(1);
            //CC1.Clear();
            //CC1.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //    int K = result.GetLength(1); // No nof Nodes
            //    for (int j = 0; j < Len; j++) {
            //        for (int k = 0; k < K; k++) {
            //            if (ccm1.Contains(j0 + j))
            //                result[j, k] = 1;
            //        }
            //    }
            //}, new CellQuadratureScheme());

            //// Double Cut Cells
            //var DCC = this.m_RegisteredFields.Where(s => s.Identification == "DCC").SingleOrDefault();
            //if (DCC == null) {
            //    DCC = new SinglePhaseField(b, "DCC");
            //    this.RegisterField(DCC);
            //}
            //var dccm = ccm0.Intersect(ccm1);
            //DCC.Clear();
            //DCC.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //    int K = result.GetLength(1); // No nof Nodes
            //    for (int j = 0; j < Len; j++) {
            //        for (int k = 0; k < K; k++) {
            //            if (dccm.Contains(j0 + j))
            //                result[j, k] = 1;
            //        }
            //    }
            //}, new CellQuadratureScheme());

            //// Cells Numbers
            //var CellNumbers = this.m_RegisteredFields.Where(s => s.Identification == "CellNumbers").SingleOrDefault();
            //if (CellNumbers == null) {
            //    CellNumbers = new SinglePhaseField(b, "CellNumbers");
            //    this.RegisterField(CellNumbers);
            //}
            //CellNumbers.Clear();
            //CellNumbers.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //    int K = result.GetLength(1); // No nof Nodes
            //    for (int j = 0; j < Len; j++) {
            //        for (int k = 0; k < K; k++) {
            //            result[j, k] = j0 + j;
            //        }
            //    }
            //}, new CellQuadratureScheme());

            //PlotAdditionalFields(physTime, timestepNo, superSampling);

            //XDGField GradT_X = new XDGField((XDGBasis)this.CurrentStateVector.Fields.Where(s => s.Identification == "Temperature").First().Basis, "GradT_X");
            //XDGField GradT_Y = new XDGField((XDGBasis)this.CurrentStateVector.Fields.Where(s => s.Identification == "Temperature").First().Basis, "GradT_Y");
            //VectorField<XDGField> GradT = new VectorField<XDGField>(GradT_X, GradT_Y);
            //GradT.Gradient(1.0, this.CurrentStateVector.Fields.Where(s => s.Identification == "Temperature").First());

            //DGField CellNumbers = new SinglePhaseField(new Basis(this.GridData, 0));
            //CellNumbers.ProjectField(1.0, delegate(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //    int K = result.GetLength(1); // No nof Nodes
            //    for (int j = 0; j < Len; j++) {
            //        for (int k = 0; k < K; k++) {
            //            result[j, k] = j0 + j;
            //        }
            //    }
            //}, new CellQuadratureScheme());
            #endregion
            //Tecplot.PlotFields(new List<DGField> { CellNumbers }, "XNSFE_GradT-" + timestepNo, physTime, 3);
            //Tecplot.PlotFields(GradT, "XNSFE_GradT-" + timestepNo, physTime, 3);

            base.PlotCurrentState(physTime, timestepNo, superSampling);
        }

        /// <summary>
        /// automatized analysis of condition number 
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis(OperatorAnalysisConfig config) {

            int D = this.GridData.SpatialDimension;

            //int[] varGroup_convDiff = Enumerable.Range(0, D).ToArray();
            //int[] varGroup_Stokes = Enumerable.Range(0, D + 1).ToArray();
            //int[] varGroup_Temperature = Enumerable.Range(D + 1, 1).ToArray();
            //int[] varGroup_all = Enumerable.Range(0, D + 2).ToArray();
            //var res = this.Timestepping.OperatorAnalysis(new[] { varGroup_convDiff, varGroup_Stokes, varGroup_Temperature, varGroup_all });

            int[] varGroup = new int[] { 0, 1, 2, 3};
            var res = this.Timestepping.OperatorAnalysis(config, new[] { varGroup });

            return res;
        }

        bool PrintOnlyOnce = true;
        private void PrintConfiguration() {
            if (PrintOnlyOnce) {
                PrintOnlyOnce = false;
                XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);

                Console.WriteLine("=============== {0} ===============", "Operator Configuration");
                PropertyInfo[] properties = typeof(XNSFE_OperatorConfiguration).GetProperties();
                foreach (PropertyInfo property in properties) {
                    if (property.PropertyType == typeof(bool)) {
                        bool s = (bool)property.GetValue(config);
                        Console.WriteLine("     {0,-30}:{1,3}", property.Name, "[" + (s == true ? "x" : " ") + "]");
                    }
                }

                Console.WriteLine("=============== {0} ===============", "Linear Solver Configuration");
                Console.WriteLine("     {0,-30}:{1}", "Solvercode", this.Control.LinearSolver.Shortname);

                Console.WriteLine("=============== {0} ===============", "Nonlinear Solver Configuration");
                Console.WriteLine("     {0,-30}:{1}", "Solvercode", this.Control.NonLinearSolver.SolverCode);
                Console.WriteLine("     {0,-30}:{1}", "Convergence Criterion", this.Control.NonLinearSolver.ConvergenceCriterion);
                Console.WriteLine("     {0,-30}:{1}", "Globalization", this.Control.NonLinearSolver.Globalization);
                Console.WriteLine("     {0,-30}:{1}", "Minsolver Iterations", this.Control.NonLinearSolver.MinSolverIterations);
                Console.WriteLine("     {0,-30}:{1}", "Maxsolver Iterations", this.Control.NonLinearSolver.MaxSolverIterations);


            }
        }     
    }
}