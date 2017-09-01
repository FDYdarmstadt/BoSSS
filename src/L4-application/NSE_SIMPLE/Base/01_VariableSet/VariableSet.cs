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

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using NSE_SIMPLE.LowMach;
using NSE_SIMPLE.Multiphase;
using BoSSS.Foundation.Grid.Classic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Physical variables.
    /// Available variables depend on <see cref="PhysicsMode"/>.
    /// </summary>
    public class VariableSet {

        private VariableSetFlowField WorkingSetFlowField;
        private VariableSetVariableDensity WorkingSetVariableDensity = null;
        private VariableSetLowMach WorkingSetLowMach = null;
        private VariableSetMultiphase WorkingSetMultiphase = null;

        private GridData GridDat;
        private SIMPLEControl Control;
        private ICollection<DGField> IOFields;
        private ICollection<DGField> RegisteredFields;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="GridDat"></param>
        /// <param name="Control"></param>
        /// <param name="IOFields"></param>
        /// <param name="RegisteredFields"></param>        
        public VariableSet(GridData GridDat, SIMPLEControl Control, ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields) {

            this.GridDat = GridDat;
            this.Control = Control;
            this.IOFields = IOFields;
            this.RegisteredFields = RegisteredFields;

            WorkingSetFlowField = new VariableSetFlowField(GridDat, Control, IOFields, RegisteredFields);
        }

        #region VariablesFlowField

        /// <summary>
        /// <see cref="VariableSetFlowField.Velocity"/>
        /// </summary>
        public VectorFieldHistory<SinglePhaseField> Velocity {
            get {
                return WorkingSetFlowField.Velocity;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.Velocity_Intrmed"/>
        /// </summary>
        public VectorField<SinglePhaseField> Velocity_Intrmed {
            get {
                return WorkingSetFlowField.Velocity_Intrmed;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.Velocity_Correction"/>
        /// </summary>
        public VectorField<SinglePhaseField> Velocity_Correction {
            get {
                return WorkingSetFlowField.Velocity_Correction;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.VelRes"/>
        /// </summary>
        public VectorField<SinglePhaseField> VelRes {
            get {
                return WorkingSetFlowField.VelRes;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.VelocityMean"/>
        /// </summary>
        public VectorField<SinglePhaseField> VelocityMean {
            get {
                return WorkingSetFlowField.VelocityMean;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.Pressure"/>
        /// </summary>
        public SinglePhaseField Pressure {
            get {
                return WorkingSetFlowField.Pressure;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.Pressure_Correction"/>
        /// </summary>
        public SinglePhaseField Pressure_Correction {
            get {
                return WorkingSetFlowField.Pressure_Correction;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.PressureRes"/>
        /// </summary>
        public SinglePhaseField PressureRes {
            get {
                return WorkingSetFlowField.PressureRes;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.DivB4"/>
        /// </summary>
        public SinglePhaseField DivB4 {
            get {
                return WorkingSetFlowField.DivB4;
            }
        }

        /// <summary>
        /// <see cref="VariableSetFlowField.DivAfter"/>
        /// </summary>
        public SinglePhaseField DivAfter {
            get {
                return WorkingSetFlowField.DivAfter;
            }
        }

        /// <summary>
        /// DG basis for velocity.
        /// </summary>
        public Basis VelBasis {
            get {
                return Velocity.Current[0].Basis;
            }
        }

        /// <summary>
        /// DG basis of <see cref="Pressure"/>.
        /// </summary>
        public Basis PressureBasis {
            get {
                return Pressure.Basis;
            }
        }

        #endregion

        #region VariablesVariableDensity

        /// <summary>
        /// <see cref="VariableSetVariableDensity.Rho"/>
        /// </summary>
        public SinglePhaseField Rho {
            get {
                return WorkingSetVariableDensity.Rho;
            }
        }

        /// <summary>
        /// <see cref="VariableSetVariableDensity.Eta"/>
        /// </summary>
        public SinglePhaseField Eta {
            get {
                return WorkingSetVariableDensity.Eta;
            }
        }

        public VectorField<SinglePhaseField> OperatorTest {
            get {
                return WorkingSetVariableDensity.OperatorTest;
            }
        }

        public VectorField<SinglePhaseField> OperatorAna {
            get {
                return WorkingSetVariableDensity.OperatorAna;
            }
        }

        #endregion

        #region VariablesLowMach

        /// <summary>
        /// <see cref="VariableSetLowMach.Temperature"/>
        /// </summary>
        public ScalarFieldHistory<SinglePhaseField> Temperature {
            get {
                return WorkingSetLowMach.Temperature;
            }
        }

        /// <summary>
        /// <see cref="VariableSetLowMach.TemperatureRes"/>
        /// </summary>
        public SinglePhaseField TemperatureRes {
            get {
                return WorkingSetLowMach.TemperatureRes;
            }
        }

        /// <summary>
        /// <see cref="VariableSetLowMach.TemperatureMean"/>
        /// </summary>
        public SinglePhaseField TemperatureMean {
            get {
                return WorkingSetLowMach.TemperatureMean;
            }
        }

        /// <summary>
        /// DG basis for temperature.
        /// </summary>
        public Basis TemperatureBasis {
            get {
                return Temperature.Current.Basis;
            }
        }

        /// <summary>
        /// <see cref="VariableSetLowMach.ThermodynamicPressure"/>
        /// </summary>
        public ScalarFieldHistory<SinglePhaseField> ThermodynamicPressure {
            get {
                return WorkingSetLowMach.ThermodynamicPressure;
            }
        }

        #endregion

        #region VariablesMultiphase

        /// <summary>
        /// <see cref="VariableSetMultiphase.Phi"/>
        /// </summary>
        public ScalarFieldHistory<SinglePhaseField> Phi {
            get {
                return WorkingSetMultiphase.Phi;
            }
        }

        /// <summary>
        /// <see cref="VariableSetMultiphase.PhiRes"/>
        /// </summary>
        public SinglePhaseField PhiRes {
            get {
                return WorkingSetMultiphase.PhiRes;
            }
        }

        /// <summary>
        /// <see cref="VariableSetMultiphase.PhiMean"/>
        /// </summary>
        public SinglePhaseField PhiMean {
            get {
                return WorkingSetMultiphase.PhiMean;
            }
        }

        /// <summary>
        /// DG basis for level-set.
        /// </summary>
        public Basis LevelSetBasis {
            get {
                return Phi.Current.Basis;
            }
        }

        #endregion

        /// <summary>
        /// Variables for flow field.
        /// </summary>
        private class VariableSetFlowField : BaseVariableSet {

            /// <summary>
            /// Auxiliary field to calculate Lambda in Convective : LinearFlux.
            /// Initialized in ctor.
            /// </summary>        
            public VectorField<SinglePhaseField> VelocityMean;

            /// <summary>
            /// Ctor.
            /// </summary>
            /// <param name="GridDat"></param>
            /// <param name="Control"></param>
            /// <param name="IOFields"></param>
            /// <param name="RegisteredFields"></param>
            /// <param name="SolverConf"></param>
            public VariableSetFlowField(GridData GridDat, SIMPLEControl Control, ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields)
                : base(GridDat, Control, IOFields, RegisteredFields) {

                // Auxiliary field to calculate Lambda in Convective : LinearFlux
                Basis MeanBasis = new Basis(GridDat, 0);
                SinglePhaseField[] _VelocityMean = new SinglePhaseField[GridDat.SpatialDimension];
                _VelocityMean[0] = new SinglePhaseField(MeanBasis, "VelocityX_Mean");
                _VelocityMean[1] = new SinglePhaseField(MeanBasis, "VelocityY_Mean");
                if (GridDat.SpatialDimension == 3)
                    _VelocityMean[2] = new SinglePhaseField(MeanBasis, "VelocityZ_Mean");
                VelocityMean = new VectorField<SinglePhaseField>(_VelocityMean);
            }

#pragma warning disable 649

            /// <summary>
            /// Velocity Field 
            /// </summary>
            /// <remarks>
            /// This field is primary, i.e. it is initialized in every possible configuration of
            /// the application.
            /// </remarks>
            [InstantiateFromControlFile(
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true, IOListOption.Always)]
            public VectorFieldHistory<SinglePhaseField> Velocity;

            /// <summary>
            /// Intermediate velocity field (before pressure correction), i.e. u*.
            /// </summary>
            [InstantiateFromControlFile(
                new string[] { VariableNames.VelocityIntermedX, VariableNames.VelocityIntermedY, VariableNames.VelocityIntermedZ },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true, IOListOption.ControlFileDetermined)]
            public VectorField<SinglePhaseField> Velocity_Intrmed;

            /// <summary>
            /// Velocity correction, i.e. u'.
            /// </summary>
            [InstantiateFromControlFile(
                new string[] { VariableNames.VelocityCorX, VariableNames.VelocityCorY, VariableNames.VelocityCorZ },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true, IOListOption.ControlFileDetermined)]
            public VectorField<SinglePhaseField> Velocity_Correction;

            /// <summary>
            /// Used for calculating Vel(n) - Vel(n-1), where n is the number of the SIMPLE iteration.
            /// </summary>
            [InstantiateFromControlFile(
                new string[] { "u_res", "v_res", "w_res" },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true, IOListOption.ControlFileDetermined)]
            public VectorField<SinglePhaseField> VelRes;

            /// <summary>
            /// Pressure field 
            /// </summary>
            /// <remarks>
            /// This field is primary in the SIMPLE algorithm.
            /// </remarks>
            [InstantiateFromControlFile(VariableNames.Pressure, VariableNames.Pressure, IOListOption.Always)]
            public SinglePhaseField Pressure;

            /// <summary>
            /// Pressure correction
            /// </summary>
            [InstantiateFromControlFile(VariableNames.PresCor, VariableNames.Pressure, IOListOption.ControlFileDetermined)]
            public SinglePhaseField Pressure_Correction;

            /// <summary>
            /// Used for calculating p(n) - p(n-1), where n is the number of the SIMPLE iteration.
            /// </summary>
            [InstantiateFromControlFile("PressureRes", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
            public SinglePhaseField PressureRes;

            /// <summary>
            /// Divergence of velocity field before pressure correction, i.e. after predictor.
            /// </summary>
            [InstantiateFromControlFile("DivB4", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
            public SinglePhaseField DivB4;

            /// <summary>
            /// Divergence of velocity field after pressure correction, if desired.
            /// </summary>
            [InstantiateFromControlFile("DivAfter", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
            public SinglePhaseField DivAfter;

#pragma warning restore 649
        }

        /// <summary>
        /// Common variables for variable density flows,
        /// i.e. LowMach and multiphase flows.
        /// </summary>
        private class VariableSetVariableDensity : BaseVariableSet {

            /// <summary>
            /// Ctor
            /// </summary>
            /// <param name="GridDat"></param>
            /// <param name="Control"></param>
            /// <param name="IOFields"></param>
            /// <param name="RegisteredFields"></param>
            public VariableSetVariableDensity(GridData GridDat, SIMPLEControl Control, ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields)
                : base(GridDat, Control, IOFields, RegisteredFields) { }

#pragma warning disable 649

            /// <summary>
            /// Density field for variable density solver.
            /// </summary>
            [InstantiateFromControlFile(VariableNames.Rho, VariableNames.VelocityX, IOListOption.Always)]
            public SinglePhaseField Rho;

            /// <summary>
            /// Dynamic viscosity for variable density solver.
            /// </summary>
            [InstantiateFromControlFile(VariableNames.ViscosityMolecular, VariableNames.VelocityX, IOListOption.Always)]
            public SinglePhaseField Eta;

            // Variables for testing purposes
            // ToDo[Klein]: can be removed after testing
            [InstantiateFromControlFile(
                new string[] { "OperatorTest_x", "OperatorTest_y", "OperatorTest_z" },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true, IOListOption.ControlFileDetermined)]
            public VectorField<SinglePhaseField> OperatorTest;

            [InstantiateFromControlFile(
                new string[] { "OperatorAna_x", "OperatorAna_y", "OperatorAna_z" },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true, IOListOption.ControlFileDetermined)]
            public VectorField<SinglePhaseField> OperatorAna;

#pragma warning restore 649
        }

        /// <summary>
        /// Variables for low Mach number flows.
        /// </summary>
        private class VariableSetLowMach : BaseVariableSet {

            /// <summary>
            /// Auxiliary field to calculate Lambda in Convective : LinearFlux.
            /// Initialized in ctor.
            /// </summary>   
            public SinglePhaseField TemperatureMean;

            /// <summary>
            /// Ctor
            /// </summary>
            /// <param name="GridDat"></param>
            /// <param name="Control"></param>
            /// <param name="IOFields"></param>
            /// <param name="RegisteredFields"></param>
            public VariableSetLowMach(GridData GridDat, SIMPLEControl Control, ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields)
                : base(GridDat, Control, IOFields, RegisteredFields) {

                // Auxiliary field to calculate Lambda in Convective : LinearFlux
                Basis MeanBasis = new Basis(GridDat, 0);
                TemperatureMean = new SinglePhaseField(MeanBasis, "TemperatureMean");

                // Thermodynamic pressure - the correct value will be calculated in VariableSet.Initialize()
                if (ThermodynamicPressure.Current.Basis.Degree != 0)
                    throw new ApplicationException("The degree for ThermodynamicPressure in control-file has to be set to zero, since the ThermodynamicPressure is constant in space.");
            }

#pragma warning disable 649

            /// <summary>
            /// Temperature field.
            /// </summary>
            [InstantiateFromControlFile(VariableNames.Temperature, VariableNames.Temperature, IOListOption.Always)]
            public ScalarFieldHistory<SinglePhaseField> Temperature;

            /// <summary>
            /// ThermodynamicPressure, i.e. p0.
            /// </summary>
            [InstantiateFromControlFile(VariableNames.ThermodynamicPressure, VariableNames.ThermodynamicPressure, IOListOption.Always)]
            public ScalarFieldHistory<SinglePhaseField> ThermodynamicPressure;

            /// <summary>
            /// Used for calculating Temperature(n) - Temperature(n-1), where n is the number of the SIMPLE iteration.
            /// </summary>
            [InstantiateFromControlFile("T_res", VariableNames.Temperature, IOListOption.Never)]
            public SinglePhaseField TemperatureRes;

#pragma warning restore 649
        }

        /// <summary>
        /// Variables for multiphase flows.
        /// </summary>
        private class VariableSetMultiphase : BaseVariableSet {

            /// <summary>
            /// Auxiliary field to calculate Lambda in Convective : LinearFlux.
            /// Initialized in ctor.
            /// </summary>   
            public SinglePhaseField PhiMean;

            /// <summary>
            /// Ctor.
            /// </summary>
            /// <param name="GridDat"></param>
            /// <param name="Control"></param>
            /// <param name="IOFields"></param>
            /// <param name="RegisteredFields"></param>
            public VariableSetMultiphase(GridData GridDat, SIMPLEControl Control, ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields)
                : base(GridDat, Control, IOFields, RegisteredFields) {

                // Auxiliary field to calculate Lambda in Convective : LinearFlux
                Basis MeanBasis = new Basis(GridDat, 0);
                PhiMean = new SinglePhaseField(MeanBasis, "PhiMean");
            }

#pragma warning disable 649

            /// <summary>
            /// Level-set field.
            /// </summary>
            [InstantiateFromControlFile(VariableNames.LevelSet, VariableNames.LevelSet, IOListOption.Always)]
            public ScalarFieldHistory<SinglePhaseField> Phi;

            /// <summary>
            /// Used for calculating Phi(n) - Phi(n-1), where n is the number of the SIMPLE iteration.
            /// </summary>
            [InstantiateFromControlFile("Phi_res", VariableNames.LevelSet, IOListOption.Never)]
            public SinglePhaseField PhiRes;

#pragma warning restore 649
        }

        /// <summary>
        /// Create extended variables for variable density flows.
        /// </summary>
        /// <param name="control"></param>
        public void CreateExtendedVariables(SIMPLEControl control) {
            switch (control.PhysicsMode) {
                case PhysicsMode.Incompressible:
                    break;
                case PhysicsMode.LowMach:
                    WorkingSetVariableDensity = new VariableSetVariableDensity(GridDat, Control, IOFields, RegisteredFields);
                    WorkingSetLowMach = new VariableSetLowMach(GridDat, Control, IOFields, RegisteredFields);
                    break;
                case PhysicsMode.Multiphase:
                    WorkingSetVariableDensity = new VariableSetVariableDensity(GridDat, Control, IOFields, RegisteredFields);
                    WorkingSetMultiphase = new VariableSetMultiphase(GridDat, Control, IOFields, RegisteredFields);
                    break;
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Initialize all variables, which are not 
        /// initialized by the control-file.
        /// Has to be called after <see cref="Application.SetInitial()"/>
        /// </summary>
        public void Initialize(SIMPLEControl SolverConf) {

            // Set history length according to time order
            if (SolverConf.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE) {
                this.Velocity.IncreaseHistoryLength(SolverConf.TimeOrder);
            }

            this.VelocityMean.Clear();
            this.VelocityMean.AccLaidBack(1.0, this.Velocity.Current);

            switch (SolverConf.PhysicsMode) {
                case PhysicsMode.Incompressible:
                    break;

                case PhysicsMode.LowMach: {
                        LowMachSIMPLEControl lowMachConf = SolverConf as LowMachSIMPLEControl;

                        // Set history length according to time order
                        if (SolverConf.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE) {
                            this.Temperature.IncreaseHistoryLength(SolverConf.TimeOrder);
                        }

                        // Initialize TemperatureMean                
                        this.TemperatureMean.Clear();
                        this.TemperatureMean.AccLaidBack(1.0, this.Temperature.Current);

                        // Initialize ThermodynamicPressure
                        this.ThermodynamicPressure.Current.Clear();
                        switch (lowMachConf.ThermodynamicPressureMode) {
                            case ThermodynamicPressureMode.Constant:
                                this.ThermodynamicPressure.Current.AccConstant(1.0);
                                break;
                            case ThermodynamicPressureMode.MassDetermined:
                                if (SolverConf.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE) {
                                    this.ThermodynamicPressure.IncreaseHistoryLength(SolverConf.TimeOrder);
                                }
                                this.ThermodynamicPressure.Current.AccConstant(
                                    lowMachConf.EoS.GetMassDeterminedThermodynamicPressure(
                                        lowMachConf.InitialMass.Value, this.Temperature.Current));
                                break;
                            default:
                                throw new NotImplementedException();
                        }

                        // Initialize EoS
                        lowMachConf.EoS.Initialize(this.ThermodynamicPressure);

                        // Initialize Density
                        this.Rho.Clear();
                        this.Rho.ProjectFunction(1.0, lowMachConf.EoS.GetDensity, null, this.Temperature.Current);

                        // Initialize Viscosity
                        this.Eta.Clear();
                        this.Eta.ProjectFunction(1.0, lowMachConf.EoS.GetViscosity, null, this.Temperature.Current);
                    }
                    break;

                case PhysicsMode.Multiphase: {
                        MultiphaseSIMPLEControl multiphaseConf = SolverConf as MultiphaseSIMPLEControl;

                        // Set history length according to time order
                        if (SolverConf.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE) {
                            this.Phi.IncreaseHistoryLength(SolverConf.TimeOrder);
                        }

                        // Initialize PhiMean                
                        this.PhiMean.Clear();
                        this.PhiMean.AccLaidBack(1.0, this.Phi.Current);

                        // Initialize Density
                        this.Rho.Clear();
                        this.Rho.ProjectFunction(1.0, multiphaseConf.EoS.GetDensity, null, this.Phi.Current);

                        // Initialize Viscosity
                        this.Eta.Clear();
                        this.Eta.ProjectFunction(1.0, multiphaseConf.EoS.GetViscosity, null, this.Phi.Current);
                    }
                    break;

                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Push WorkingSet to next time step.
        /// </summary>
        public void Push(SIMPLEControl SolverConf) {
            if (SolverConf.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE) {
                this.Velocity.Push();
                switch (SolverConf.PhysicsMode) {
                    case PhysicsMode.Incompressible:
                        break;
                    case PhysicsMode.LowMach:
                        this.Temperature.Push();
                        this.ThermodynamicPressure.Push();
                        break;
                    case PhysicsMode.Multiphase:
                        this.Phi.Push();
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// Checks if any field contains NaN's or Inf's and throws an exception if this is the case.
        /// </summary>
        public bool CheckForNanOrInf(SIMPLEControl SolverConf) {
            bool FoundNanOrInf = false;

            switch (SolverConf.PhysicsMode) {
                case PhysicsMode.Incompressible:
                    if (WorkingSetFlowField.CheckForNanOrInf())
                        FoundNanOrInf = true;
                    break;
                case PhysicsMode.LowMach:
                    if (WorkingSetFlowField.CheckForNanOrInf() || WorkingSetVariableDensity.CheckForNanOrInf() || WorkingSetLowMach.CheckForNanOrInf())
                        FoundNanOrInf = true;
                    break;
                case PhysicsMode.Multiphase:
                    if (WorkingSetFlowField.CheckForNanOrInf() || WorkingSetVariableDensity.CheckForNanOrInf() || WorkingSetMultiphase.CheckForNanOrInf())
                        FoundNanOrInf = true;
                    break;
                default:
                    throw new NotImplementedException();
            }

            return FoundNanOrInf;
        }
    }
}
