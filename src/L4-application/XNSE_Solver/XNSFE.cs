using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
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
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Extension of the <see cref="XNSE"/>-solver for additional heat transfer.
    /// (The 'F' stands for Fourier equation, i.e. Heat equation.)
    /// Changed to Newton Solver 4/2021, Picard might give unexpected results - MR
    /// </summary>
    public class XNSFE : XNSE<XNSFE_Control> {
       
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
        /// override of <see cref="DefineMomentumEquation(OperatorFactory, XNSFE_OperatorConfiguration, int, int)"/>
        /// adding evaporation extension for Navier-Stokes equations
        /// </summary>
        /// <param name="opFactory"></param>
        /// <param name="config"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        protected override void DefineMomentumEquation(OperatorFactory opFactory, XNSFE_OperatorConfiguration config, int d, int D) {
            base.DefineMomentumEquation(opFactory, config, d, D);

            // === evaporation extension === //
            if (config.isEvaporation) {
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddEquation(new InterfaceNSE_Evaporation("A", "B", D, d, config));                    
                } else if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    NonlinearCouplingEvaporation = true; // When using Newton evaporation coupling is always nonlinear
                    opFactory.AddEquation(new InterfaceNSE_Evaporation_Newton("A", "B", D, d, config));
                }                    
            }
            if (config.isBuoyancy && config.isGravity) {
                opFactory.AddEquation(new NavierStokesBuoyancy("A", D, d, config));
                opFactory.AddEquation(new NavierStokesBuoyancy("B", D, d, config));
            }
        }


        protected override void DefineContinuityEquation(OperatorFactory opFactory, XNSFE_OperatorConfiguration config, int D) {
            base.DefineContinuityEquation(opFactory, config, D);

            // === evaporation extension === //
            if (config.isEvaporation) {
                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Picard) {
                    opFactory.AddEquation(new InterfaceContinuity_Evaporation("A", "B", D, config));
                } else if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    opFactory.AddEquation(new InterfaceContinuity_Evaporation_Newton("A", "B", D, config));
                }
            }
        }

        protected override void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystemImmersedBoundary(D, opFactory, lsUpdater);

            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);

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
                        throw new ApplicationException("Warning using LDG Formulation for Heat, this is untested. Remove this statement only if you now what you are doing!");
                        opFactory.AddEquation(new HeatFlux("A", d, D, thermBoundaryMap, config));
                        opFactory.AddEquation(new HeatFlux("B", d, D, thermBoundaryMap, config));
                        opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, thermBoundaryMap, config));
                    }
                }
            } else {
                throw new NotSupportedException();
            }
        }

        bool NonlinearCouplingEvaporation = false;
        protected override XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {

            OperatorFactory opFactory = new OperatorFactory();

            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            XOP.FreeMeanValue[VariableNames.Pressure] = !GetBcMap().DirichletPressureBoundary;
            XOP.IsLinear = !(this.Control.PhysicalParameters.IncludeConvection || this.Control.ThermalParameters.IncludeConvection || Control.NonlinearCouplingSolidFluid || NonlinearCouplingEvaporation); // when using Newton solver the employed coupling between heat and momentum is nonlinear
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

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            base.PlotCurrentState(physTime, timestepNo, superSampling);

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

            //Tecplot.PlotFields(new List<DGField> { CellNumbers }, "XNSFE_GradT-" + timestepNo, physTime, 3);
            //Tecplot.PlotFields(GradT, "XNSFE_GradT-" + timestepNo, physTime, 3);

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
                Console.WriteLine("     {0,-30}:{1}", "Solvercode", this.Control.LinearSolver.SolverCode);
                if (this.Control.LinearSolver.SolverCode != LinearSolverCode.classic_mumps & this.Control.LinearSolver.SolverCode != LinearSolverCode.classic_pardiso) {
                    Console.WriteLine("TODO");
                }

                Console.WriteLine("=============== {0} ===============", "Nonlinear Solver Configuration");
                Console.WriteLine("     {0,-30}:{1}", "Solvercode", this.Control.NonLinearSolver.SolverCode);
                Console.WriteLine("     {0,-30}:{1}", "Convergence Criterion", this.Control.NonLinearSolver.ConvergenceCriterion);
                Console.WriteLine("     {0,-30}:{1}", "Globalization", this.Control.NonLinearSolver.Globalization);
                Console.WriteLine("     {0,-30}:{1}", "Minsolver Iterations", this.Control.NonLinearSolver.MinSolverIterations);
                Console.WriteLine("     {0,-30}:{1}", "Maxsolver Iterations", this.Control.NonLinearSolver.MaxSolverIterations);


            }
        }
        /*

        protected override LevelSetUpdater InstantiateLevelSetUpdater() {
            int levelSetDegree = Control.FieldOptions["Phi"].Degree;
            LevelSetUpdater lsUpdater;

           
            
            switch(Control.Option_LevelSetEvolution) {
                case LevelSetEvolution.Fourier: {
                    if(Control.EnforceLevelSetConservation) {
                        throw new NotSupportedException("mass conservation correction currently not supported");
                    }
                    FourierLevelSet fourrierLevelSet = new FourierLevelSet(Control.FourierLevSetControl, new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                    fourrierLevelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, fourrierLevelSet, VariableNames.LevelSetCG);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, levelSetVelocity);
                    break;
                }
                case LevelSetEvolution.FastMarching: {
                    LevelSet levelSet = new LevelSet(new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                    levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                    var fastMarcher = new FastMarchingEvolver(VariableNames.LevelSetCG, QuadOrder(), levelSet.GridDat.SpatialDimension);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet, VariableNames.LevelSetCG);
                    lsUpdater.AddEvolver(VariableNames.LevelSetCG, fastMarcher);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, levelSetVelocity);
                    break;
                }
                case LevelSetEvolution.StokesExtension: {

                    throw new NotImplementedException("todo");
                }
                case LevelSetEvolution.SplineLS: {
                    int nodeCount = 30;
                    Console.WriteLine("Achtung, Spline node count ist hart gesetzt. Was soll hier hin?");
                    SplineLevelSet SplineLevelSet = new SplineLevelSet(Control.Phi0Initial, new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG, nodeCount);
                    var SplineEvolver = new SplineLevelSetEvolver(VariableNames.LevelSetCG, (GridData)SplineLevelSet.GridDat);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, SplineLevelSet, VariableNames.LevelSetCG);
                    lsUpdater.AddEvolver(VariableNames.LevelSetCG, SplineEvolver);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, levelSetVelocity);
                    break;
                }
                case LevelSetEvolution.None: {
                    LevelSet levelSet1 = new LevelSet(new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                    levelSet1.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet1, VariableNames.LevelSetCG);
                    break;
                }
                default:
                throw new NotImplementedException();
            }
            
            return lsUpdater;
        }

        */



    }
}