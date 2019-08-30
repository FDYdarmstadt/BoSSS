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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using XDGShock;
using XDGShock.Fluxes;
using XDGShock.ShockCapturing;
using XDGShock.TimeStepping;
using XDGShock.Variables;

namespace XDGShockTest {

    public class Program : Application<XDGShockControl> {
        static void Main(string[] args) {
            _Main(args, false, () => new Program());
        }

        /// <summary>
        /// The density $\rho$
        /// </summary>
        public XDGField Density {
            get;
            private set;
        }

        /// <summary>
        /// The momentum field $\rho \vec{u}$
        /// </summary>
        public VectorField<XDGField> Momentum {
            get;
            private set;
        }

        /// <summary>
        /// The total energy per volume $\rho E$
        /// </summary>
        public XDGField Energy {
            get;
            private set;
        }

        /// <summary>
        /// An XDG shock sensor
        /// </summary>
        public IXDGSensor Sensor {
            get;
            private set;
        }

        public SinglePhaseField ArtificialViscosityField {
            get;
            private set;
        }

        public IArtificialViscosityLaw ArtificialViscosityLaw {
            get;
            private set;
        }

        /// <summary>
        /// All conservative fields (density, momentum vector, energy)
        /// </summary>
        public XDGField[] ConservativeFields {
            get {
                int D = CompressibleEnvironment.NumberOfDimensions;
                XDGField[] fields = new XDGField[D + 2];

                fields[0] = Density;
                for (int d = 0; d < D; d++) {
                    fields[d + 1] = Momentum[d];
                }
                fields[D + 1] = Energy;

                return fields;
            }
        }

        public Dictionary<DerivedVariable<XDGField>, XDGField> DerivedVariableToXDGFieldMap {
            get;
            private set;
        }

        public Dictionary<DerivedVariable<SinglePhaseField>, SinglePhaseField> DerivedVariableToSinglePhaseFieldMap {
            get;
            private set;
        }

        private XDGShockTimeStepping TimeStepper {
            get;
            set;
        }

        public LevelSetTracker LevelSetTracker {
            get;
            private set;
        }

        public LevelSet LevelSet {
            get;
            private set;
        }

        public XSpatialOperatorMk2 XSpatialOperator {
            get;
            private set;
        }

        public XDGField[] Residuals {
            get;
            private set;
        }

        public XDGField Check {
            get;
            private set;
        }

        protected override void CreateFields() {

            // Create level set and level set tracker
            this.LevelSet = new LevelSet(new Basis(this.GridData, this.Control.LevelSetDegree), XDGShockVariables.LevelSet);
            this.LevelSetTracker = new LevelSetTracker((GridData)this.GridData, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 1, new string[] { "A", "B" }, this.LevelSet);
            base.LsTrk = this.LevelSetTracker;

            // Create mandatory conservative fields
            int D = CompressibleEnvironment.NumberOfDimensions;

            this.Density = new XDGField(new XDGBasis(this.LevelSetTracker, this.Control.DensityDegree), CompressibleVariables.Density);
            XDGField[] momentumFields = new XDGField[D];
            XDGBasis momentumBasis = new XDGBasis(this.LevelSetTracker, this.Control.MomentumDegree);
            for (int d = 0; d < D; d++) {
                string variableName = CompressibleVariables.Momentum[d];
                momentumFields[d] = new XDGField(momentumBasis, variableName);
            }
            this.Momentum = new VectorField<XDGField>(momentumFields);
            this.Energy = new XDGField(new XDGBasis(this.LevelSetTracker, this.Control.EnergyDegree), CompressibleVariables.Energy);

            // Create fields for residuals
            XDGBasis basis = new XDGBasis(this.LevelSetTracker, this.Control.DensityDegree);

            this.Residuals = new XDGField[D + 2];
            this.Residuals[0] = new XDGField(basis, "residual_" + CompressibleVariables.Density);
            for (int d = 0; d < D; d++) {
                this.Residuals[d + 1] = new XDGField(basis, "residual_" + CompressibleVariables.Momentum[d]);
            }
            this.Residuals[D + 1] = new XDGField(basis, "residual_" + CompressibleVariables.Energy);

            // Create derived fields
            this.DerivedVariableToXDGFieldMap = new Dictionary<DerivedVariable<XDGField>, XDGField>();
            this.DerivedVariableToSinglePhaseFieldMap = new Dictionary<DerivedVariable<SinglePhaseField>, SinglePhaseField>();
            foreach (KeyValuePair<Variable, int> pair in this.Control.VariableToDegreeMap) {
                Variable variable = pair.Key;
                int degree = pair.Value;

                if (variable is DerivedVariable<XDGField> xdgVar) {
                    this.DerivedVariableToXDGFieldMap.Add(xdgVar, new XDGField(new XDGBasis(this.LevelSetTracker, degree), variable.Name));
                } else if (variable is DerivedVariable<SinglePhaseField> singleVar) {
                    this.DerivedVariableToSinglePhaseFieldMap.Add(singleVar, new SinglePhaseField(new Basis(this.GridData, degree), variable.Name));
                }
            }

            // Register all fields for plotting
            // Fields for the conservative variables are plotted in any case
            this.m_IOFields.Add(LevelSet);
            this.m_IOFields.AddRange(ConservativeFields);
            this.m_IOFields.AddRange(DerivedVariableToXDGFieldMap.Values);
            this.m_IOFields.AddRange(DerivedVariableToSinglePhaseFieldMap.Values);

            // Register fields for, e.g. applying the initial conditions;
            this.m_RegisteredFields.Add(LevelSet);
            this.m_RegisteredFields.AddRange(ConservativeFields);
            this.m_RegisteredFields.AddRange(DerivedVariableToXDGFieldMap.Values);
            this.m_RegisteredFields.AddRange(DerivedVariableToSinglePhaseFieldMap.Values);

            // Add sensor
            if (this.Control.SensorVariable != null) {
                switch (this.Control.SensorType) {
                    case SensorTypes.PerssonSensor:
                        this.Sensor = new XDGPerssonSensor(ConservativeFields.First(s => s.Identification.Equals(this.Control.SensorVariable)), this.Control.SensorLimit, this.LevelSetTracker, this.NonlinearQuadratureDegree);
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }

            // Initialize artificial viscosity field, if artificial viscosity has been activated in the control file
            if (this.Control.VariableToDegreeMap.ContainsKey(XDGShockVariables.ArtificialViscosity)) {
                this.ArtificialViscosityField = new SinglePhaseField(new Basis(this.GridData, 2), "artificialViscosity");

                switch (this.Control.ArtificialViscosityLawType) {
                    case ArtificialViscosityLawTypes.SmoothedHeaviside:
                        this.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(this.Sensor, ConservativeFields.First(s => s.Identification.Equals(this.Control.SensorVariable)).Basis.Degree, this.Control.SensorLimit);
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        protected override IGrid CreateOrLoadGrid() {
            IGrid grid = base.CreateOrLoadGrid();
            CompressibleEnvironment.Initialize(grid.SpatialDimension);
            return grid;
        }

        private int NonlinearQuadratureDegree {
            get {
                return Math.Max(2, 3 * this.Density.Basis.Degree);
            }
        }

        private void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {

            IList<DGField> parameterFields = null;
            if (this.ArtificialViscosityField != null) {
                parameterFields = new DGField[] { this.ArtificialViscosityField };
            }

            if (OpMtx != null) {
                var mtxB = this.XSpatialOperator.GetMatrixBuilder(this.LevelSetTracker, Mapping, parameterFields, Mapping, this.LevelSetTracker.SpeciesIdS.ToArray());
                mtxB.ComputeMatrix(OpMtx, OpAffine);
            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin ev = this.XSpatialOperator.GetEvaluatorEx(this.LevelSetTracker, CurrentState, parameterFields, Mapping, this.LevelSetTracker.SpeciesIdS.ToArray());
                ev.Evaluate(1.0, 0.0, OpAffine, null);
            }
        }

        private double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {
            this.LevelSetTracker.UpdateTracker();
            return 0.0;
        }

        private IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                Dictionary<SpeciesId, IEnumerable<double>> massScaleToSpeciesMap = new Dictionary<SpeciesId, IEnumerable<double>>();

                int D = this.Grid.SpatialDimension;
                double[] ones = new double[D + 2];
                ones.SetAll(1.0);

                foreach (SpeciesId speciesId in this.LsTrk.SpeciesIdS) {
                    massScaleToSpeciesMap.Add(speciesId, ones);
                }
                return massScaleToSpeciesMap;
            }
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            // Boundary condition map
            Material material = this.Control.GetMaterial();
            IBoundaryConditionMap boundaryMap = new XDGCompressibleBoundaryCondMap(this.GridData, this.Control, material, new string[] { "A", "B" });

            // Operator
            string[] variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };

            if (ArtificialViscosityField != null) {
                this.XSpatialOperator = new XSpatialOperatorMk2(variables, new string[] { ArtificialViscosityField.Identification }, variables, (int[] A, int[] B, int[] C) => NonlinearQuadratureDegree, LsTrk.SpeciesIdS.ToArray());
            } else {
                this.XSpatialOperator = new XSpatialOperatorMk2(variables, null, variables, (int[] A, int[] B, int[] C) => NonlinearQuadratureDegree, LsTrk.SpeciesIdS.ToArray());
            }

            // Bulk fluxes
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new OptimizedHLLCDensityFlux(boundaryMap, material));
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new OptimizedHLLCMomentumFlux(boundaryMap, 0, material));
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new OptimizedHLLCMomentumFlux(boundaryMap, 1, material));
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new OptimizedHLLCEnergyFlux(boundaryMap, material));

            // Interface fluxes
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new OptimizedHLLCFlux_Interface(this.LevelSetTracker, boundaryMap, material, FluxVariables.Density));
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new OptimizedHLLCFlux_Interface(this.LevelSetTracker, boundaryMap, material, FluxVariables.Momentum, 0));
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new OptimizedHLLCFlux_Interface(this.LevelSetTracker, boundaryMap, material, FluxVariables.Momentum, 1));
            this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new OptimizedHLLCFlux_Interface(this.LevelSetTracker, boundaryMap, material, FluxVariables.Energy));

            // Artificial viscosity bulk fluxes
            if (this.ArtificialViscosityField != null) {
                GridData gridData = (GridData)this.GridData;

                this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(
                    new OptimizedLaplacianArtificialViscosityFlux(
                        gridData,
                        CompressibleVariables.Density,
                        penaltySafetyFactor: 1.0,
                        penaltyFactor: (Density.Basis.Degree + 1) * (Density.Basis.Degree + gridData.SpatialDimension) / gridData.SpatialDimension,
                        cellLengthScales: gridData.Cells.CellLengthScale
                        ));

                for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum[d]].Add(
                        new OptimizedLaplacianArtificialViscosityFlux(
                            gridData,
                            CompressibleVariables.Momentum[d],
                            penaltySafetyFactor: 1.0,
                            penaltyFactor: (Momentum[0].Basis.Degree + 1) * (Momentum[0].Basis.Degree + gridData.SpatialDimension) / gridData.SpatialDimension,
                            cellLengthScales: gridData.Cells.CellLengthScale
                            ));
                }
                this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(
                    new OptimizedLaplacianArtificialViscosityFlux(
                        gridData,
                        CompressibleVariables.Energy,
                        penaltySafetyFactor: 1.0,
                        penaltyFactor: (Energy.Basis.Degree + 1) * (Energy.Basis.Degree + gridData.SpatialDimension) / gridData.SpatialDimension,
                        cellLengthScales: gridData.Cells.CellLengthScale
                        ));
            }

            this.XSpatialOperator.Commit();

            // Timestepper
            this.TimeStepper = new XDGShockTimeStepping(
                this.ConservativeFields,
                this.Residuals,
                this.LevelSetTracker,
                this.DelComputeOperatorMatrix,
                this.DelUpdateLevelset,
                RungeKuttaScheme.ExplicitEuler,
                LevelSetHandling.None,
                MassMatrixShapeandDependence.IsIdentity,
                SpatialOperatorType.Nonlinear,
                this.MassScale,
                this.MultigridOperatorConfig,
                base.MultigridSequence,
                this.LevelSetTracker.SpeciesIdS.ToArray(),
                this.NonlinearQuadratureDegree,
                this.Control.AgglomerationThreshold,
                true,
                this.Control.NonLinearSolver,
                this.Control.LinearSolver);
        }

        MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int p = this.Momentum[0].Basis.Degree;
                int D = this.GridData.SpatialDimension;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 2];

                    // configuration for continuity
                    configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, p - iLevel),
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { 0 }
                    };


                    // configurations for momentum equation
                    for (int d = 0; d < D; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            Degree = Math.Max(0, p - iLevel),
                            mode = MultigridOperator.Mode.Eye,
                            VarIndex = new int[] { d + 1 }
                        };
                    }

                    // configuration for energy
                    configs[iLevel][D + 1] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, p - iLevel),
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { D + 1 }
                    };
                }

                return configs;
            }
        }

        protected override void SetInitial() {
            int D = CompressibleEnvironment.NumberOfDimensions;
            CellQuadratureScheme scheme = new CellQuadratureScheme(true, CellMask.GetFullMask(this.GridData));

            // Level set
            this.LevelSet.ProjectField(1.0, this.Control.LevelSetPos, scheme);
            this.LevelSetTracker.UpdateTracker();

            if (this.Control.uAInitPrimitive != null && this.Control.uBInitPrimitive != null) {
                // Initial conditions are specified in primitive variables

                // Density
                this.Density.Clear();
                this.Density.GetSpeciesShadowField("A").ProjectField(1.0, this.Control.uAInitPrimitive[0]);
                this.Density.GetSpeciesShadowField("B").ProjectField(1.0, this.Control.uBInitPrimitive[0]);

                // Momentum (rho * momentum vector)
                for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                    this.Momentum[d].Clear();
                    this.Momentum[d].GetSpeciesShadowField("A").ProjectField(1.0, X => this.Control.uAInitPrimitive[0](X) * this.Control.uAInitPrimitive[d + 1](X), scheme);
                    this.Momentum[d].GetSpeciesShadowField("B").ProjectField(1.0, X => this.Control.uBInitPrimitive[0](X) * this.Control.uBInitPrimitive[d + 1](X), scheme);
                }

                // Energy
                this.Energy.Clear();
                Energy.GetSpeciesShadowField("A").ProjectField(
                    1.0,
                    delegate (double[] X) {
                        double rho = this.Control.uAInitPrimitive[0](X);
                        double p = this.Control.uAInitPrimitive[D + 1](X);
                        Vector u = new Vector(D);
                        for (int d = 0; d < D; d++) {
                            u[d] = this.Control.uAInitPrimitive[d + 1](X);
                        }

                        StateVector state = StateVector.FromPrimitiveQuantities(this.Control.GetMaterial(), rho, u, p);
                        return state.Energy;
                    },
                    scheme);
                Energy.GetSpeciesShadowField("B").ProjectField(
                    1.0,
                    delegate (double[] X) {
                        double rho = this.Control.uBInitPrimitive[0](X);
                        double p = this.Control.uBInitPrimitive[D + 1](X);
                        Vector u = new Vector(D);
                        for (int d = 0; d < D; d++) {
                            u[d] = this.Control.uBInitPrimitive[d + 1](X);
                        }

                        StateVector state = StateVector.FromPrimitiveQuantities(this.Control.GetMaterial(), rho, u, p);
                        return state.Energy;
                    },
                     scheme);
            } else {
                // Initial conditions are specified in conservative variables

                // Density
                this.Density.Clear();
                this.Density.GetSpeciesShadowField("A").ProjectField(1.0, this.Control.uAInitConservative[0]);
                this.Density.GetSpeciesShadowField("B").ProjectField(1.0, this.Control.uBInitConservative[0]);

                // Momentum
                for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                    this.Momentum[d].Clear();
                    this.Momentum[d].GetSpeciesShadowField("A").ProjectField(1.0, this.Control.uAInitConservative[d + 1]);
                    this.Momentum[d].GetSpeciesShadowField("B").ProjectField(1.0, this.Control.uBInitConservative[d + 1]);
                }

                // Energy
                this.Energy.Clear();
                this.Energy.GetSpeciesShadowField("A").ProjectField(1.0, this.Control.uAInitConservative[D + 1]);
                this.Energy.GetSpeciesShadowField("B").ProjectField(1.0, this.Control.uBInitConservative[D + 1]);
            }

            // Update sensor values
            this.Sensor?.UpdateSensorValues(this.Density);

            // Update artificial viscosity
            if (this.ArtificialViscosityField != null) {
                UpdateArtificialViscosity();
            }
        }

        private PlotDriver plotDriver;

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            if (plotDriver == null) {
                plotDriver = new Tecplot(GridData, true, false, (uint)superSampling);
            }
            UpdateDerivedVariables();
            plotDriver.PlotFields("XDGShockTest-" + timestepNo, physTime, m_IOFields);
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            int printInterval = Control.PrintInterval;
            if (DatabaseDriver.MyRank == 0 && TimestepNo % printInterval == 0) {
#if DEBUG
                Console.WriteLine();
#endif
                Console.Write("Starting time step #" + TimestepNo + "...");
            }

            // Set time step size manually
            if (dt <= 0) {
                dt = this.Control.dtFixed;
            }

            Exception e = null;
            try {
                TimeStepper.Solve(phystime, dt);
            } catch (Exception ee) {
                e = ee;
            }
            e.ExceptionBcast();

            // Update sensor
            this.Sensor?.UpdateSensorValues(this.Density);

            // Update artificial viscosity
            if (this.ArtificialViscosityField != null) {
                UpdateArtificialViscosity();
            }

            if (DatabaseDriver.MyRank == 0 && TimestepNo % printInterval == 0) {
                if (TimestepNo % printInterval == 0) {
                    Console.WriteLine(" done. PhysTime: {0:0.#######E-00}, dt: {1:0.#######E-00}", phystime, dt);
                }
            }

            return dt;
        }

        /// <summary>
        /// Makes sure all derived variables are updated before saving
        /// </summary>
        /// <param name="timestepno"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        protected override ITimestepInfo SaveToDatabase(TimestepNumber timestepno, double t) {
            UpdateDerivedVariables();
            return base.SaveToDatabase(timestepno, t);
        }

        private void UpdateDerivedVariables() {
            //// Update sensor _field_ when used
            //this.Sensor?.UpdateSensorValues(this.Density);

            //// Update artificial viscosit when used
            //if (this.ArtificialViscosityField != null) {
            //    UpdateArtificialViscosity();
            //}

            // Update derived variables (including sensor _variable_)
            foreach (KeyValuePair<DerivedVariable<XDGField>, XDGField> pair in this.DerivedVariableToXDGFieldMap) {
                pair.Key.UpdateFunction(pair.Value, this);
            }

            foreach (KeyValuePair<DerivedVariable<SinglePhaseField>, SinglePhaseField> pair in this.DerivedVariableToSinglePhaseFieldMap) {
                pair.Key.UpdateFunction(pair.Value, this);
            }
        }

        private SpecFemBasis avSpecFEMBasis;

        private void UpdateArtificialViscosity() {
            // Determine piecewise constant viscosity
            this.ArtificialViscosityField.Clear();

            // Cell masks 
            CellMask cutCells = this.LevelSetTracker.Regions.GetCutCellMask();
            CellMask speciesA_nonCut = this.LevelSetTracker.Regions.GetSpeciesMask("A").Except(cutCells);
            CellMask speciesB_nonCut = this.LevelSetTracker.Regions.GetSpeciesMask("B").Except(cutCells);

            // Calculate AV
            SetAVForSpecies(speciesA_nonCut, "A");
            SetAVForSpecies(speciesB_nonCut, "B");
            SetAVForCutCells(cutCells);

            // Set AV manually for testing
            //this.ArtificialViscosityField.Clear();
            //this.ArtificialViscosityField.AccConstant(1.0);

            // Project visocsity onto continuous, multilinear space
            if (CompressibleEnvironment.NumberOfDimensions < 3) {
                // Standard version
                if (avSpecFEMBasis == null || !this.ArtificialViscosityField.Basis.Equals(avSpecFEMBasis.ContainingDGBasis)) {
                    avSpecFEMBasis = new SpecFemBasis((GridData)this.GridData, 2);
                }
                SpecFemField specFemField = new SpecFemField(avSpecFEMBasis);
                specFemField.ProjectDGFieldMaximum(1.0, this.ArtificialViscosityField);
                this.ArtificialViscosityField.Clear();
                specFemField.AccToDGField(1.0, this.ArtificialViscosityField);
            } else {
                throw new NotImplementedException("Artificial viscosity has only been tested for 2D");
            }
        }

        private void SetAVForSpecies(CellMask speciesMask, string speciesName) {
            int D = CompressibleEnvironment.NumberOfDimensions;

            foreach (Chunk chunk in speciesMask) {
                MultidimensionalArray meanValues = MultidimensionalArray.Create(chunk.Len, D + 2);
                this.Density.GetSpeciesShadowField(speciesName).EvaluateMean(chunk.i0, chunk.Len, meanValues.ExtractSubArrayShallow(-1, 0), 0, 0.0);
                for (int d = 0; d < D; d++) {
                    this.Momentum[d].GetSpeciesShadowField(speciesName).EvaluateMean(chunk.i0, chunk.Len, meanValues.ExtractSubArrayShallow(-1, d + 1), 0, 0.0);
                }
                this.Energy.GetSpeciesShadowField(speciesName).EvaluateMean(chunk.i0, chunk.Len, meanValues.ExtractSubArrayShallow(-1, D + 1), 0, 0.0);

                for (int i = 0; i < chunk.Len; i++) {
                    int cell = chunk.i0 + i;
                    StateVector state = new StateVector(meanValues.ExtractSubArrayShallow(i, -1).To1DArray(), this.Control.GetMaterial());

                    Debug.Assert(!double.IsNaN(state.SpeedOfSound), "State vector is NaN! Okay if level set is outside of the domain");

                    double localViscosity = this.ArtificialViscosityLaw.GetViscosity(cell, ((GridData)this.GridData).Cells.h_min[cell], state);
                    Debug.Assert(localViscosity >= 0.0);

                    this.ArtificialViscosityField.SetMeanValue(cell, localViscosity);
                }
            }
        }

        private void SetAVForCutCells(CellMask cutCells) {
            int D = CompressibleEnvironment.NumberOfDimensions;

            foreach (Chunk chunk in cutCells) {
                // Species A
                MultidimensionalArray meanValuesA = MultidimensionalArray.Create(chunk.Len, D + 2);
                this.Density.GetSpeciesShadowField("A").EvaluateMean(chunk.i0, chunk.Len, meanValuesA.ExtractSubArrayShallow(-1, 0), 0, 0.0);
                for (int d = 0; d < D; d++) {
                    this.Momentum[d].GetSpeciesShadowField("A").EvaluateMean(chunk.i0, chunk.Len, meanValuesA.ExtractSubArrayShallow(-1, d + 1), 0, 0.0);
                }
                this.Energy.GetSpeciesShadowField("A").EvaluateMean(chunk.i0, chunk.Len, meanValuesA.ExtractSubArrayShallow(-1, D + 1), 0, 0.0);

                // Species B
                MultidimensionalArray meanValuesB = MultidimensionalArray.Create(chunk.Len, D + 2);
                this.Density.GetSpeciesShadowField("B").EvaluateMean(chunk.i0, chunk.Len, meanValuesB.ExtractSubArrayShallow(-1, 0), 0, 0.0);
                for (int d = 0; d < D; d++) {
                    this.Momentum[d].GetSpeciesShadowField("B").EvaluateMean(chunk.i0, chunk.Len, meanValuesB.ExtractSubArrayShallow(-1, d + 1), 0, 0.0);
                }
                this.Energy.GetSpeciesShadowField("B").EvaluateMean(chunk.i0, chunk.Len, meanValuesB.ExtractSubArrayShallow(-1, D + 1), 0, 0.0);

                // Loop over cells
                for (int i = 0; i < chunk.Len; i++) {
                    int cell = chunk.i0 + i;
                    StateVector stateA = new StateVector(meanValuesA.ExtractSubArrayShallow(i, -1).To1DArray(), this.Control.GetMaterial());
                    StateVector stateB = new StateVector(meanValuesB.ExtractSubArrayShallow(i, -1).To1DArray(), this.Control.GetMaterial());

                    StateVector stateToUse;
                    if (stateA.SpeedOfSound + stateA.Velocity.Abs() >= stateB.SpeedOfSound + stateB.Velocity.Abs()) {
                        stateToUse = stateA;
                    } else {
                        stateToUse = stateB;
                    }
                    Debug.Assert(!double.IsNaN(stateToUse.SpeedOfSound), "Speed of sound is NaN");

                    double localViscosity = this.ArtificialViscosityLaw.GetViscosity(cell, ((GridData)this.GridData).Cells.h_min[cell], stateToUse);
                    Debug.Assert(localViscosity >= 0.0);

                    this.ArtificialViscosityField.SetMeanValue(cell, localViscosity);
                }
            }
        }
    }
}
