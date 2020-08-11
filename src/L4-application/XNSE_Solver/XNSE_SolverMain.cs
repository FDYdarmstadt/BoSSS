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
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Numerics;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.Tracing;
using ilPSP.LinSolvers;

using BoSSS.Platform;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.XDG;

using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.Grid.Aggregation;
using NUnit.Framework;
using MPI.Wrappers;
using System.Collections;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Solver for Incompressible Multiphase flows
    /// Optional: coupled heat equation with evaporation
    /// Optional: kinetic energy equation 
    /// </summary>
    public partial class XNSE_SolverMain : BoSSS.Solution.Application<XNSE_Control> {

        //===========
        // Main file
        //===========

        static void Main(string[] args) {

            //InitMPI();
            //DeleteOldPlotFiles();
            //BoSSS.Application.XNSE_Solver.Tests.UnitTest.BcTest_PressureOutletTest(1, 0.0d, true);
            ////Tests.UnitTest.OneTimeTearDown();
            //return;


            _Main(args, false, delegate () {
                var p = new XNSE_SolverMain();
                return p;
            });
        }

        //=====================================
        // Field declaration and instantiation
        //=====================================
        #region fields

#pragma warning disable 649

        /// <summary>
        /// Bundling of variables which are either DG or XDG (see <see cref="XNSE_Control.UseXDG4Velocity"/>);
        /// </summary>
        class VelocityRelatedVars<TX> where TX : DGField {
            /// <summary>
            /// velocity
            /// </summary>
            [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                null,
                true, true,
                IOListOption.ControlFileDetermined)]
            public VectorField<TX> Velocity;


            /// <summary>
            /// Volume Force, dimension is acceleration, i.e. length per time-square.
            /// </summary>
            [InstantiateFromControlFile(
                new string[] { VariableNames.GravityX, VariableNames.GravityY, VariableNames.GravityZ },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true,
                IOListOption.ControlFileDetermined)]
            public VectorField<TX> Gravity;

            /// <summary>
            /// Residual in the momentum equation.
            /// </summary>
            [InstantiateFromControlFile(new string[] { "ResidualMomentumX", "ResidualMomentumY", "ResidualMomentumZ" },
                new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
                true, true,
                IOListOption.ControlFileDetermined)]
            public VectorField<TX> ResidualMomentum;
        }
                
        /// <summary>
        /// Velocity and related variables for the extended case.
        /// </summary>
        VelocityRelatedVars<XDGField> XDGvelocity;


        /// <summary>
        /// Pressure
        /// </summary>
        //[InstantiateFromControlFile(VariableNames.Pressure, null, IOListOption.ControlFileDetermined)]
        XDGField Pressure;

        /// <summary>
        /// Residual of the continuity equation
        /// </summary>
        //[InstantiateFromControlFile("ResidualConti", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        XDGField ResidualContinuity;


        /// <summary>
        /// Artificial force term at the fluid interface, usually only to support manufactured solutions.
        /// </summary>
        [InstantiateFromControlFile(
            new string[] { VariableNames.SurfaceForceX, VariableNames.SurfaceForceY, VariableNames.SurfaceForceZ },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        VectorField<SinglePhaseField> SurfaceForce;

#pragma warning restore 649


        protected override void CreateFields() {
            using (new FuncTrace()) {
                base.CreateFields();


                this.CreateLevelSetFields();


                this.Pressure = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree), VariableNames.Pressure);
                base.RegisterField(this.Pressure);
                this.ResidualContinuity = new XDGField(this.Pressure.Basis, "ResidualConti");
                base.RegisterField(this.ResidualContinuity);

                this.XDGvelocity = new VelocityRelatedVars<XDGField>();
                InitFromAttributes.CreateFieldsAuto(this.XDGvelocity, this.GridData, base.Control.FieldOptions, base.Control.CutCellQuadratureType, base.IOFields, base.m_RegisteredFields);

                this.CreateUtilityFields();

                this.CreateEnergyFields();

                if (this.Control.solveCoupledHeatEquation) 
                    this.CreateHeatFields();

            }
        }


        #endregion



        //==========================================
        // operator related members
        // (create and update operator/mass matrix)
        //==========================================
        #region operator


        /// <summary>
        /// the spatial operator (momentum and continuity equation)
        /// Optional: kinetic energy / heat equation
        /// </summary>
        XNSFE_OperatorFactory XNSFE_Operator;

        /// <summary>
        /// OperatorConfiguration for the <see cref="XNSFE_Operator"/>
        /// </summary>
        XNSFE_OperatorConfiguration XOpConfig;



        /// <summary>
        /// output of <see cref="AssembleMatrix"/>;
        /// </summary>
        MassMatrixFactory MassFact;
        
        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double rho_A = this.Control.PhysicalParameters.rho_A,
                    rho_B = this.Control.PhysicalParameters.rho_B;

                double c_A = this.Control.ThermalParameters.c_A,
                    c_B = this.Control.ThermalParameters.c_B;

                int D = this.GridData.SpatialDimension;


                double[] scale_A = new double[D + 1];
                double[] scale_B = new double[D + 1];
                int mD = D + 1;
                if (this.Control.solveKineticEnergyEquation) {
                    scale_A = new double[D + 2];
                    scale_B = new double[D + 2];
                    mD = D + 2;
                }
                if (this.Control.solveCoupledHeatEquation) {
                    scale_A = new double[mD + 1];
                    scale_B = new double[mD + 1];
                    if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                        scale_A = new double[mD + 1 + D];
                        scale_B = new double[mD + 1 + D];
                    }
                }

                scale_A.SetAll(rho_A); // mass matrix in momentum equation (kinetic energy equation)
                scale_A[D] = 0; // no  mass matrix for continuity equation
                scale_B.SetAll(rho_B); // mass matrix in momentum equation (kinetic energy equation)
                scale_B[D] = 0; // no  mass matrix for continuity equation

                if (this.Control.solveCoupledHeatEquation) {
                    scale_A[mD + 1] = rho_A * c_A;
                    scale_B[mD + 1] = rho_B * c_B;
                    if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                        scale_A.GetSubVector(mD + 1, D).SetAll(0.0);
                        scale_B.GetSubVector(mD + 1, D).SetAll(0.0);
                    }
                }

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), scale_A);
                R.Add(this.LsTrk.GetSpeciesId("B"), scale_B);


                return R;
            }
        }



        /// <summary>
        /// Boundary conditions.
        /// </summary>
        IncompressibleMultiphaseBoundaryCondMap BcMap {
            get {
                if (m_BcMap == null) {
                    m_BcMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
                }
                return m_BcMap;
            }
        }

        IncompressibleMultiphaseBoundaryCondMap m_BcMap;


        /// <summary>
        /// Current velocity and pressure;
        /// </summary>
        internal CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(this.CurrentVel, this.Pressure));
                    if (this.Control.solveKineticEnergyEquation) {
                        m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.KineticEnergy));
                    }
                    if (this.Control.solveCoupledHeatEquation) {
                        //m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(this.CurrentVel, this.Pressure, this.Temperature));
                        m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.Temperature));
                        if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
                            m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(m_CurrentSolution.Mapping.Fields.ToArray(), this.HeatFlux));
                    }
                } else {
                    for (int d = 0; d < base.GridData.SpatialDimension; d++) {
                        Debug.Assert(object.ReferenceEquals(m_CurrentSolution.Mapping.Fields[d], this.CurrentVel[d]));
                    }
                }

                return m_CurrentSolution;
            }
        }

        CoordinateVector m_CurrentSolution;

        /// <summary>
        /// Current residual for momentum and continuity equation.
        /// </summary>
        internal CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat<DGField>(XDGvelocity.ResidualMomentum, ResidualContinuity));
                    if (this.Control.solveKineticEnergyEquation) {
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.ResidualKineticEnergy));
                    }
                    if (this.Control.solveCoupledHeatEquation) {
                        //m_CurrentResidual = new CoordinateVector(ArrayTools.Cat<DGField>(XDGvelocity.ResidualMomentum, ResidualContinuity, ResidualHeat));
                        m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.ResidualHeat));
                        if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
                            m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(m_CurrentResidual.Mapping.Fields.ToArray(), this.ResidualAuxHeatFlux));
                    }
                }
                return m_CurrentResidual;
            }
        }

        CoordinateVector m_CurrentResidual;


        /// <summary>
        /// Current Velocity
        /// </summary>
        XDGField[] CurrentVel {
            get {
                return this.XDGvelocity.Velocity.ToArray();
            }
        }


        /// <summary>
        /// HMF order/degree which is used globally in this solver.
        /// </summary>
        int m_HMForder;



        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int pVel = this.CurrentVel[0].Basis.Degree;
                int pPrs = this.Pressure.Basis.Degree;
                int D = this.GridData.SpatialDimension;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {

                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 1];
                    int mD = D + 1;
                    if (this.Control.solveKineticEnergyEquation) {
                        configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 2];
                        mD = D + 2;
                    }
                    if (this.Control.solveCoupledHeatEquation) {
                        configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[mD + 1];
                        if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
                            configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[mD + 1 + D];
                    }

                    //if (this.Control.solveCoupledHeatEquation) {
                    //    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 2];
                    //    if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
                    //        configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 2 + D];
                    //} else {
                    //    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 1];
                    //}

                    // configurations for velocity
                    for (int d = 0; d < D; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                            VarIndex = new int[] { d }
                        };
                    }
                    // configuration for pressure
                    configs[iLevel][D] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { D }
                    };
                    
                    if (this.Control.solveKineticEnergyEquation) {
                        int pKinE = this.KineticEnergy.Basis.Degree;
                        // configuration for kinetic energy
                        configs[iLevel][D + 1] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(1, pKinE - iLevel) },
                            mode = this.Control.KineticEnergyeBlockPrecondMode,
                            VarIndex = new int[] { D + 1 }
                        };
                    }

                    if (this.Control.solveCoupledHeatEquation) {

                        int pTemp = this.Temperature.Basis.Degree;
                        // configuration for Temperature
                        configs[iLevel][mD + 1] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(1, pTemp - iLevel) },
                            mode = this.Control.TemperatureBlockPrecondMode,
                            VarIndex = new int[] { mD + 1 }
                        };

                        // configuration for auxiliary heat flux
                        if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                            int pFlux = this.HeatFlux[0].Basis.Degree;
                            for (int d = 0; d < D; d++) {
                                configs[iLevel][mD + 1 + d] = new MultigridOperator.ChangeOfBasisConfig() {
                                    DegreeS = new int[] { Math.Max(1, pFlux - iLevel) },
                                    mode = MultigridOperator.Mode.Eye,
                                    VarIndex = new int[] { mD + 1 + d }
                                };
                            }
                        }

                    }

                }


                return configs;
            }
        }


        /// <summary>
        /// Create XOperator and Timestepper
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            #region Checks
            // CreateEquationsAndSolvers might be called multiple times
            // exit if so, and no LoadBalancing
            if (XNSFE_Operator != null && L == null)
                return;


            if (Control.TimesteppingMode == AppControl._TimesteppingMode.Steady) {
                if (Control.Timestepper_LevelSetHandling != LevelSetHandling.None)
                    throw new ApplicationException(string.Format("Illegal control file: for a steady computation ({0}), the level set handling must be {1}.", AppControl._TimesteppingMode.Steady, LevelSetHandling.None));
            }

            int degU = this.CurrentVel[0].Basis.Degree;

            #endregion


            #region Config and Generate XOperator

            //Quadrature Order
            //----------------

            m_HMForder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if (this.Control.solveKineticEnergyEquation)
                m_HMForder *= 2;


            // Create Spatial Operator
            // ======================= 

            XOpConfig = new XNSFE_OperatorConfiguration(this.Control);

            XNSFE_Operator = new XNSFE_OperatorFactory(XOpConfig, this.LsTrk, this.m_HMForder, this.BcMap, this.thermBcMap, degU);
            updateSolutionParams = new bool[CurrentResidual.Mapping.Fields.Count];

            #endregion


            #region Create Timestepper
            // =======================
            if (L == null) {

                this.CreateTimestepper();

            } else {

                //PlotCurrentState(hack_Phystime, new TimestepNumber(hack_TimestepIndex, 12), 2);

                Debug.Assert(object.ReferenceEquals(this.MultigridSequence[0].ParentGrid, this.GridData));


                DGField[] flds = ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure);
                DGField[] resi = ArrayTools.Cat<DGField>(this.XDGvelocity.ResidualMomentum.ToArray(), this.ResidualContinuity);

                if (this.Control.solveKineticEnergyEquation) {
                    flds = ArrayTools.Cat<DGField>(flds, this.KineticEnergy);
                    resi = ArrayTools.Cat<DGField>(resi, this.ResidualKineticEnergy);
                }

                if (this.Control.solveCoupledHeatEquation) {
                    flds = ArrayTools.Cat<DGField>(flds, this.Temperature);
                    resi = ArrayTools.Cat<DGField>(resi, this.ResidualHeat);
                    if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                        flds = ArrayTools.Cat<DGField>(flds, this.HeatFlux);
                        resi = ArrayTools.Cat<DGField>(resi, this.ResidualAuxHeatFlux);
                    }
                }

                m_BDF_Timestepper.DataRestoreAfterBalancing(L, flds, resi, this.LsTrk, this.MultigridSequence);


                //PlotCurrentState(hack_Phystime, new TimestepNumber(hack_TimestepIndex, 13), 2);

                ContinuityEnforcer = new ContinuityProjection(
                    ContBasis: this.LevSet.Basis,
                    DGBasis: this.DGLevSet.Current.Basis,
                    gridData: GridData,
                    Option: Control.LSContiProjectionMethod);

                if (this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
                    ReInitPDE = new EllipticReInit(this.LsTrk, this.Control.ReInitControl, DGLevSet.Current);
                    FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                    ExtVelMover = new ExtensionVelocityBDFMover(LsTrk, DGLevSet.Current, DGLevSetGradient, new VectorField<DGField>(XDGvelocity.Velocity.ToArray()),
                    Control.EllipticExtVelAlgoControl, BcMap, bdfOrder, ExtensionVelocity.Current, new double[2] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B });
                }

            }
            #endregion

        }


        void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {
            using (var tr = new FuncTrace()) {
                int D = this.GridData.SpatialDimension;

                // ============================
                // treatment of surface tension
                // ============================

                VectorField<SinglePhaseField> filtLevSetGradient;
                switch (this.Control.AdvancedDiscretizationOptions.SST_isotropicMode) {
                    case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                    case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                    case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine: {
                            CurvatureAlgorithms.LaplaceBeltramiDriver(
                                this.Control.AdvancedDiscretizationOptions.SST_isotropicMode,
                                this.Control.AdvancedDiscretizationOptions.FilterConfiguration,
                                out filtLevSetGradient, this.LsTrk,
                                this.DGLevSet.Current);
                            if ((this.Control.solveKineticEnergyEquation && !this.LsTrk.Regions.GetCutCellMask().IsEmptyOnRank) || XOpConfig.isEvaporation) {
                                VectorField<SinglePhaseField> filtLevSetGradient_dummy;
                                CurvatureAlgorithms.CurvatureDriver(
                                    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                                    CurvatureAlgorithms.FilterConfiguration.Default,
                                    this.Curvature, out filtLevSetGradient_dummy, this.LsTrk,
                                    this.m_HMForder,
                                    this.DGLevSet.Current);
                            }
                            break;
                        }
                    case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                    case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                    case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                        CurvatureAlgorithms.CurvatureDriver(
                            this.Control.AdvancedDiscretizationOptions.SST_isotropicMode,
                            this.Control.AdvancedDiscretizationOptions.FilterConfiguration,
                            this.Curvature, out filtLevSetGradient, this.LsTrk,
                            this.m_HMForder,
                            this.DGLevSet.Current);
                        //CurvatureAlgorithms.MakeItConservative(LsTrk, this.Curvature, this.Control.PhysicalParameters.Sigma, this.SurfaceForce, filtLevSetGradient, MomentFittingVariant, this.m_HMForder);
                        break;

                    case SurfaceStressTensor_IsotropicMode.Curvature_Fourier:
                        if (Fourier_LevSet != null) {
                            Fourier_LevSet.ProjectToDGCurvature(this.Curvature, out filtLevSetGradient, this.LsTrk.Regions.GetCutCellMask());
                        } else {
                            throw new NotImplementedException("Curvature_Fourier needs an instance of Fourier_LevSet");
                        }
                        break;

                    default: throw new NotImplementedException("Unknown SurfaceTensionMode");
                }

                if (filtLevSetGradient != null) {
                    if (this.Control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource == CurvatureAlgorithms.LevelSetSource.fromC0) {
                        //this.LevSetGradient.Clear();
                        //this.LevSetGradient.Acc(1.0, filtLevSetGradient);
                        this.LevSetGradient = filtLevSetGradient;
                    } else {
                        //this.DGLevSetGradient.Clear();
                        //this.DGLevSetGradient.Acc(1.0, filtLevSetGradient);
                        this.DGLevSetGradient = filtLevSetGradient;
                    }
                }

                // =================================================
                // Construct evaporative mass flux (extension field) 
                // =================================================

                // heat flux for evaporation
                //DGField[] HeatFluxParam = new DGField[D];
                //if (XOpConfig.solveHeat) {
                //    if (XOpConfig.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                //        HeatFluxParam = new VectorField<XDGField>(D, CurrentState.ToArray()[D + 1].Basis, "HeatFlux0_", XDGField.Factory).ToArray();
                //        Dictionary<string, double> kSpc = new Dictionary<string, double>();
                //        kSpc.Add("A", -this.Control.ThermalParameters.k_A);
                //        kSpc.Add("B", -this.Control.ThermalParameters.k_B);
                //        XNSEUtils.ComputeGradientForParam(CurrentState.ToArray()[D + 1], HeatFluxParam, this.LsTrk, kSpc, this.LsTrk.Regions.GetNearFieldSubgrid(1));
                //    } else {
                //        var HeatFluxMap = new CoordinateMapping(CurrentState.ToArray().GetSubVector(D + 2, D));
                //        HeatFluxParam = HeatFluxMap.Fields.ToArray();
                //    }
                //}
                //ConventionalDGField[] HeatFluxAParam = new VectorField<ConventionalDGField>(D.ForLoop(d => (HeatFluxParam[d] as XDGField).GetSpeciesShadowField("A"))).ToArray();
                //ConventionalDGField[] HeatFluxBParam = new VectorField<ConventionalDGField>(D.ForLoop(d => (HeatFluxParam[d] as XDGField).GetSpeciesShadowField("B"))).ToArray();

                //SinglePhaseField[] HeatFluxAExt = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(new Basis(this.GridData, HeatFluxParam[d].Basis.Degree), "HeatFluxExt" + d))).ToArray();
                //SinglePhaseField[] HeatFluxBExt = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(new Basis(this.GridData, HeatFluxParam[d].Basis.Degree), "HeatFluxExt" + d))).ToArray();
                //double[][] ExtVelMin = new double[HeatFluxAExt.Length][];
                //double[][] ExtVelMax = new double[HeatFluxAExt.Length][];
                //int J = this.LsTrk.GridDat.Cells.NoOfLocalUpdatedCells;
                //for (int i = 0; i < HeatFluxAExt.Length; i++) {
                //    ExtVelMin[i] = new double[J];
                //    ExtVelMax[i] = new double[J];
                //}
                //NarrowMarchingBand.ConstructExtVel_PDE(this.LsTrk, this.LsTrk.Regions.GetCutCellSubGrid(), HeatFluxAExt, HeatFluxAParam, this.LevSet, this.LevSetGradient, ExtVelMin, ExtVelMax, m_HMForder);
                //NarrowMarchingBand.ConstructExtVel_PDE(this.LsTrk, this.LsTrk.Regions.GetCutCellSubGrid(), HeatFluxBExt, HeatFluxBParam, this.LevSet, this.LevSetGradient, ExtVelMin, ExtVelMax, m_HMForder);

                //DGField[] HeatFluxExtParam = new VectorField<XDGField>(D, CurrentState.ToArray()[D + 1].Basis, "HeatFluxExt_", XDGField.Factory).ToArray();
                //for (int d = 0; d < D; d++) {
                //    (HeatFluxExtParam[d] as XDGField).GetSpeciesShadowField("A").Acc(1.0, HeatFluxAExt[d]);
                //    (HeatFluxExtParam[d] as XDGField).GetSpeciesShadowField("B").Acc(1.0, HeatFluxBExt[d]);
                //}

                //Tecplot.PlotFields(HeatFluxParam, "HeatFluxParam" + hack_TimestepIndex, hack_Phystime, 2);
                //Tecplot.PlotFields(HeatFluxExtParam, "HeatFluxExtParam" + hack_TimestepIndex, hack_Phystime, 2);

                // ============================
                // matrix assembly
                // ============================

                var codMap = Mapping;
                var domMap = Mapping;

                using (new BlockTrace("XdgMatrixAssembly",tr)) {
                    this.XNSFE_Operator.AssembleMatrix(
                       OpMtx, OpAffine, codMap, domMap,
                       CurrentState, AgglomeratedCellLengthScales, phystime,
                       this.m_HMForder, SurfaceForce, filtLevSetGradient, Curvature,
                       updateSolutionParams, this.XDGvelocity.Gravity.ToArray()); //, HeatFluxExtParam);
                                                                                  //(this.Control.solveCoupledHeatEquation ? this.Temperature.ToEnumerable() : null),
                                                                                  //(this.Control.solveCoupledHeatEquation ? this.DisjoiningPressure.ToEnumerable() : null));                   
                }


                // ====================================
                // something with surface tension ?????
                // ====================================

                {
                    if (this.Control.PhysicalParameters.useArtificialSurfaceForce == true)
                        throw new NotSupportedException("Not supported for this hack.");


                    var TmpRhs = new CoordinateVector(CurrentState.Select(f => (DGField)f.Clone()).ToArray());
                    TmpRhs.Clear();

                    var VelA = new CoordinateVector(TmpRhs.Mapping.Fields.Take(D).Select(f => (DGField)(((XDGField)f).GetSpeciesShadowField("A"))).ToArray());
                    var VelB = new CoordinateVector(TmpRhs.Mapping.Fields.Take(D).Select(f => (DGField)(((XDGField)f).GetSpeciesShadowField("B"))).ToArray());

                    int N = ((XDGBasis)(CurrentState[0].Basis)).NonX_Basis.Length;

                    //foreach (int jCell in this.LsTrk.Regions.GetCutCellMask4LevSet(0).ItemEnum) {
                    //    for (int d = 0; d < D; d++) {
                    //        for (int n = 0; n < N; n++) {
                    //            ((XDGField)(TmpRhs.Mapping.Fields[d])).GetSpeciesShadowField("A").Coordinates[jCell, n] = 0.5 * SurfaceForce[d].Coordinates[jCell, n];
                    //            ((XDGField)(TmpRhs.Mapping.Fields[d])).GetSpeciesShadowField("B").Coordinates[jCell, n] = 0.5 * SurfaceForce[d].Coordinates[jCell, n];
                    //        }
                    //    }
                    //}

                    OpAffine.AccV(1.0, TmpRhs);
                }

                // so far, 'SaddlePointRHS' is on the left-hand-side, since it is the output of ComputeMatrix
                // multiply by -1 to make it RHS
                OpAffine.ScaleV(-1.0);


                // ============================
                // Generate MassMatrix
                // ============================

                // mass matrix factory
                MassFact = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), m_HMForder, 1).MassMatrixFactory;// new MassMatrixFactory(maxB, CurrentAgg);
                var WholeMassMatrix = MassFact.GetMassMatrix(Mapping, MassScale); // mass matrix scaled with density rho


                // ============================
                //  Add Gravity
                // ============================
                // Dimension: [ rho * G ] = mass / time^2 / len^2 == [ d/dt rho U ]
                var WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity.Gravity.ToArray<DGField>(), new XDGField(this.Pressure.Basis)));
                if (this.Control.solveKineticEnergyEquation) {
                    WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(WholeGravity.Mapping.Fields, new XDGField(this.KineticEnergy.Basis)));
                }
                if (this.Control.solveCoupledHeatEquation) {
                    WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(WholeGravity.Mapping.Fields, new XDGField(this.Temperature.Basis)));
                    if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
                        WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(WholeGravity.Mapping.Fields,
                            new XDGField(this.Temperature.Basis), new VectorField<XDGField>(D, this.HeatFlux[0].Basis, XDGField.Factory)));
                }
                WholeMassMatrix.SpMV(1.0, WholeGravity, 1.0, OpAffine);


                // ============================
                // Set Pressure Reference Point
                // ============================

                if (OpMtx != null) {
                    if (!this.BcMap.DirichletPressureBoundary) {
                        XNSEUtils.SetPressureReferencePoint(
                            Mapping,
                            this.GridData.SpatialDimension,
                            this.LsTrk, OpMtx, OpAffine);
                    }
                } else {
                    if (!this.BcMap.DirichletPressureBoundary) {
                        XNSEUtils.SetPressureReferencePointResidual(
                            new CoordinateVector(CurrentState),
                            this.GridData.SpatialDimension,
                            this.LsTrk, OpAffine);
                    }
                }

                // transform from RHS to Affine
                OpAffine.ScaleV(-1.0);

            }
        }
        
        bool[] updateSolutionParams;

        bool lockUpdate;

        bool m_TransformedResi = true;

        protected void SolutionParamsUpdate(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {

            int NF = this.CurrentResidual.Mapping.Fields.Count;
            int D = this.GridData.SpatialDimension;

            int len = (this.Control.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) ? 1 : 1 + D;
            for (int l = 0; l < len; l++)
                updateSolutionParams[D + 1 + l] = false;

            double[] L2Res = new double[NF];

            if (m_TransformedResi) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // transform current solution and residual back to the DG domain
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                var R = m_BDF_Timestepper.Residuals;
                R.Clear();

                Mgop.TransformRhsFrom(R, currentRes);
                this.LsTrk.GetAgglomerator(this.LsTrk.SpeciesIdS.ToArray(), m_HMForder, this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                    AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true).Extrapolate(R.Mapping);

                for (int i = 0; i < NF; i++) {
                    L2Res[i] = R.Mapping.Fields[i].L2Norm();
                }

            } else {
                // +++++++++++++++++++++++
                // un-transformed residual
                // +++++++++++++++++++++++

                var VarIdx = NF.ForLoop(i => Mgop.Mapping.GetSubvectorIndices(i));

                for (int i = 0; i < VarIdx.Length; i++) {
                    double L2 = 0.0;
                    foreach (int idx in VarIdx[i])
                        L2 += currentRes[idx - Mgop.Mapping.i0].Pow2();
                    L2Res[i] = L2.MPISum().Sqrt();
                }
            }

            double NSE_L2Res = 0.0;
            for (int i = 0; i <= D; i++)
                NSE_L2Res += L2Res[i].Pow2();
            NSE_L2Res.Sqrt();

            //Console.WriteLine("NSE L2 residual = {0}", NSE_L2Res);
            if (!lockUpdate && (NSE_L2Res < this.Control.NonLinearSolver.ConvergenceCriterion || iterIndex == 0)) {
                Console.WriteLine("update solution param");
                for (int l = 0; l < len; l++)
                    updateSolutionParams[D + 1 + l] = true;

                lockUpdate = (iterIndex > 0);
            }


        }


        #endregion



        //=============================
        // timestepper related members
        //=============================
        #region timestepper


        /// <summary>
        /// Implicit timestepping using Backward-Differentiation-Formulas (BDF),
        /// specialized for XDG applications.
        /// </summary>
        XdgBDFTimestepping m_BDF_Timestepper;


        ///// <summary>
        ///// Explicit or implicit timestepping using Runge-Kutta formulas,
        ///// specialized for XDG applications.
        ///// </summary>
        //XdgRKTimestepping m_RK_Timestepper;

        RungeKuttaScheme rksch = null;

        int bdfOrder = -1000;


        private void CreateTimestepper() {

            switch (this.Control.TimeSteppingScheme) {
                case TimeSteppingScheme.RK_ImplicitEuler: {
                        rksch = RungeKuttaScheme.ImplicitEuler;
                        break;
                    }
                case TimeSteppingScheme.RK_CrankNic: {
                        rksch = RungeKuttaScheme.CrankNicolson;
                        break;
                    }
                case TimeSteppingScheme.CrankNicolson: {
                        //do not instantiate rksch, use bdf instead
                        bdfOrder = -1;
                        break;
                    }
                case TimeSteppingScheme.ImplicitEuler: {
                        //do not instantiate rksch, use bdf instead
                        bdfOrder = 1;
                        break;
                    }
                default: {
                        if (this.Control.TimeSteppingScheme.ToString().StartsWith("BDF")) {
                            //do not instantiate rksch, use bdf instead
                            bdfOrder = Convert.ToInt32(this.Control.TimeSteppingScheme.ToString().Substring(3));
                            break;
                        } else
                            throw new NotImplementedException();
                    }

            }


            if (rksch == null) {
                m_BDF_Timestepper = new XdgBDFTimestepping(
                    this.CurrentSolution.Mapping.Fields,
                    this.CurrentResidual.Mapping.Fields,
                    LsTrk,
                    true,
                    DelComputeOperatorMatrix, null, DelUpdateLevelSet,
                    (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? bdfOrder : 1,
                    this.Control.Timestepper_LevelSetHandling,
                    this.XOpConfig.mmsd,
                    (this.Control.PhysicalParameters.IncludeConvection) ? SpatialOperatorType.Nonlinear : SpatialOperatorType.LinearTimeDependent,
                    MassScale,
                    this.MultigridOperatorConfig, base.MultigridSequence,
                    this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder,
                    this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                    true,
                    this.Control.NonLinearSolver,
                    this.Control.LinearSolver
                    );
                m_BDF_Timestepper.m_ResLogger = base.ResLogger;
                m_BDF_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
                m_BDF_Timestepper.Timestepper_Init = (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? this.Control.Timestepper_BDFinit : TimeStepperInit.SingleInit;
                m_BDF_Timestepper.incrementTimesteps = this.Control.incrementTimesteps;
                m_BDF_Timestepper.PushLevelSet = this.PushLevelSetAndRelatedStuff;
                m_BDF_Timestepper.IterUnderrelax = this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? this.Control.LSunderrelax : 1.0;

                m_BDF_Timestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;
                //m_BDF_Timestepper.CustomIterationCallback += this.PlotOnIterationCallback;
                if (this.Control.useSolutionParamUpdate)
                    m_BDF_Timestepper.CustomIterationCallback += this.SolutionParamsUpdate;

                // solver 
                this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations; //m_BDF_Timestepper.config_NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;

                if (this.Control.NonLinearSolver.SolverCode == NonLinearSolverCode.Newton) {
                    m_BDF_Timestepper.XdgSolverFactory.Selfmade_precond =
                                        new Schwarz() {
                                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                                NoOfPartsPerProcess = this.CurrentSolution.Count / 10000,
                                            },
                                            Overlap = 1,
                                            CoarseSolver = new SparseSolver() { WhichSolver = SparseSolver._whichSolver.MUMPS }
                                        };
                } else {
                    //m_BDF_Timestepper.Config_linearSolver = new DirectSolver() { WhichSolver = this.Control.LinearSolver };
                }

                //Console.WriteLine("noofpartsperprocess = {0}", this.CurrentSolution.Count / 10000);   

            } else {

                throw new NotSupportedException();

                //m_RK_Timestepper = new XdgRKTimestepping(
                //    this.CurrentSolution.Mapping.Fields.ToArray(),
                //    this.CurrentResidual.Mapping.Fields.ToArray(),
                //    LsTrk,
                //    DelComputeOperatorMatrix, DelUpdateLevelSet, DelUpdateCutCellMetrics,
                //    rksch,
                //    this.Control.Timestepper_LevelSetHandling,
                //    mmsd,
                //    (this.Control.PhysicalParameters.IncludeConvection) ? SpatialOperatorType.Nonlinear : SpatialOperatorType.LinearTimeDependent,
                //    MassScale,
                //    this.MultigridOperatorConfig, base.MultigridSequence,
                //    this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold, 
                //    true);
                //m_RK_Timestepper.m_ResLogger = base.ResLogger;
                //m_RK_Timestepper.m_ResidualNames = this.CurrentResidual.Mapping.Fields.Select(f => f.Identification).ToArray();
            }

        }


        /// <summary>
        /// delegate for the initialization of previous timesteps from an analytic solution
        /// </summary>
        /// <param name="TimestepIndex"></param>
        /// <param name="Time"></param>
        /// <param name="St"></param>
        private void BDFDelayedInitSetIntial(int TimestepIndex, double Time, DGField[] St) {
            using (new FuncTrace()) {
                Console.WriteLine("Timestep index {0}, time {1} ", TimestepIndex, Time);

                // level-set
                // ---------
                this.DGLevSet.Current.ProjectField(X => this.Control.Phi(X, Time));
                this.LevSet.ProjectField(X => this.Control.Phi(X, Time));

                this.LsTrk.UpdateTracker(incremental: true);

                // solution
                // --------
                int D = this.LsTrk.GridDat.SpatialDimension;

                for (int d = 0; d < D; d++) {
                    XDGField _u = (XDGField)St[d];
                    _u.Clear();
                    _u.GetSpeciesShadowField("A").ProjectField(X => this.Control.ExactSolutionVelocity["A"][d](X, Time));
                    _u.GetSpeciesShadowField("B").ProjectField((X => this.Control.ExactSolutionVelocity["B"][d](X, Time)));
                }
                XDGField _p = (XDGField)St[D];
                _p.Clear();
                _p.GetSpeciesShadowField("A").ProjectField(X => this.Control.ExactSolutionPressure["A"](X, Time));
                _p.GetSpeciesShadowField("B").ProjectField((X => this.Control.ExactSolutionPressure["B"](X, Time)));
            }
        }


        /// <summary>
        /// delegate for the initialization of previous timesteps from restart session
        /// </summary>
        /// <param name="TimestepIndex"></param>
        /// <param name="time"></param>
        /// <param name="St"></param>
        private void BDFDelayedInitLoadRestart(int TimestepIndex, double time, DGField[] St) {

            Console.WriteLine("Timestep index {0}, time {1} ", TimestepIndex, time);

            ITimestepInfo tsi_toLoad;
            if (TimestepIndex < 0) {
                throw new ArgumentOutOfRangeException("Not enough Timesteps to restart with desired Timestepper");
            } else {
                ISessionInfo reloadSession = GetDatabase().Controller.GetSessionInfo(this.CurrentSessionInfo.RestartedFrom);
                tsi_toLoad = reloadSession.Timesteps.Single(t => t.TimeStepNumber.Equals(new TimestepNumber(TimestepIndex)));
            }
            DatabaseDriver.LoadFieldData(tsi_toLoad, this.GridData, this.IOFields);

            // level-set
            // ---------
            this.DGLevSet.Current.Clear();
            this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);

            //this.LsTrk.UpdateTracker(incremental: true);

            // solution
            // --------
            int D = this.LsTrk.GridDat.SpatialDimension;

            for (int d = 0; d < D; d++) {
                St[d] = this.XDGvelocity.Velocity[d].CloneAs();
            }
            St[D] = this.Pressure.CloneAs();
        }


        #endregion



        //===================================================
        // application related members
        // (RunSolverOneStep, SetInitial, LoadRestart, etc.)
        //===================================================
        #region application

        int hack_TimestepIndex;
        double hack_Phystime;

        /// <summary>
        /// Depending on settings <see cref="AppControl.CompMode"/>, computes either one timestep or a steady-state solution.
        /// </summary>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using(var tr = new FuncTrace()) {

                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;
                base.ResLogger.TimeStep = TimestepInt;
                hack_TimestepIndex = TimestepInt;
                hack_Phystime = phystime;

  
                Preprocessing(TimestepInt, phystime, dt, TimestepNo);


                if(Control.SkipSolveAndEvaluateResidual) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    // setup: project exact solution -- for consistency tests
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    foreach(string spc in LsTrk.SpeciesNames) {
                        for(int d = 0; d < this.GridData.SpatialDimension; d++) {
                            ConventionalDGField Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);
                            Vel_d.ProjectField(Control.ExactSolutionVelocity[spc][d].Convert_Xt2X(phystime + dt));
                        }
                        Pressure.GetSpeciesShadowField(spc).ProjectField(Control.ExactSolutionPressure[spc].Convert_Xt2X(phystime + dt));
                    }
                }


                // =====================================================
                // setup stationary 
                // =====================================================


                if(base.Control.TimesteppingMode == AppControl._TimesteppingMode.Steady) {
                    dt = 1.0e100;
                    Console.WriteLine("Steady-state solve ...", TimestepNo, dt);

                    if(this.Control.Option_LevelSetEvolution != LevelSetEvolution.None) {
                        throw new ApplicationException("For steady-state solutions, the only allowed level-set-evolution option is '" + LevelSetEvolution.None + "'.");
                    }



                // =====================================================
                // setup transient 
                // =====================================================
                } else if(base.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {

                    // push stacks
                    // -----------

                    PushLevelSetAndRelatedStuff();


                    // backup old velocity/kinetic energy for energy checks
                    // -----------------------------------------------------
                    if(this.Control.ComputeEnergyProperties && this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
                        for(int d = 0; d < D; d++) {
                            this.prevVel[d].Clear();
                            this.prevVel[d].Acc(1.0, this.CurrentVel[d]);
                        }
                        if (this.Control.ComputeEnergyProperties) {
                            this.prevKineticEnergy.Clear();
                            this.prevKineticEnergy.Acc(1.0, this.DerivedKineticEnergy);
                        }
                    }


                    // fields setup
                    // ------------
                    for (int d = 0; d < D; d++) {
                        // Gravity must be set up like this to avoid regions of zero gravity when updating the level-set
                        this.XDGvelocity.Gravity[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
                    }


                    // +++++++++++++++++++++++++++++++++++++
                    // compute/check time step restrictions
                    // +++++++++++++++++++++++++++++++++++++

                    dt = base.Control.dtFixed;

                    // Level-Set motion-CFL
                    double LevSet_Deg2 = this.DGLevSet.Current.Basis.Degree;
                    LevSet_Deg2 = LevSet_Deg2 * LevSet_Deg2;
                    double dt_LevSetCFL = base.GridData.ComputeCFLTime(this.ExtensionVelocity.Current, dt * LevSet_Deg2);
                    dt_LevSetCFL = dt_LevSetCFL / LevSet_Deg2;
                    if(this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative && this.Control.LSunderrelax == 1.0) {
                        if(dt / dt_LevSetCFL > 1.0) {
                            double underrelax = Math.Round(dt_LevSetCFL / dt, 1);
                            m_BDF_Timestepper.IterUnderrelax = underrelax;
                            Console.WriteLine("Exceeding Level-Set CFL: Setting underrelaxation factor to {0}", underrelax);
                            if (this.Control.solveKineticEnergyEquation && (dt / (dt_LevSetCFL / 4.0)) > 1.0) {
                                Console.WriteLine("Exceeding Level-Set CFL for kinetic energy equation: dt = {0}, dt_sigma = {1}, frac = {2}", dt, dt_LevSetCFL/4.0, (dt / (dt_LevSetCFL / 4.0)));
                            }
                        } else {
                            m_BDF_Timestepper.IterUnderrelax = 1.0;
                        }
                    }

                    if (this.Control.solveKineticEnergyEquation) { 
                        double kinE_Deg2 = this.KineticEnergy.Basis.Degree;
                        kinE_Deg2 = kinE_Deg2 * kinE_Deg2;
                        double dt_kinECFL = base.GridData.ComputeCFLTime(this.KineticEnergy.ToEnumerable(), dt * kinE_Deg2);
                        dt_kinECFL = dt_kinECFL / kinE_Deg2;
                        if (dt / dt_kinECFL > 1.0) {
                            Console.WriteLine("Exceeding CFL-condition for kinetic energy equation: dt = {0}, dt_sigma = {1}, frac = {2}", dt, dt_kinECFL, dt / dt_kinECFL);
                        }
                    }

                    //dt = Math.Min(dt, dt_LevSetCFL);

                    // Capillary Timestep restriction
                    if (this.Control.PhysicalParameters.Sigma != 0.0) {
                        MultidimensionalArray h_mins = ((GridData)this.GridData).Cells.h_min;
                        double h = h_mins.Min();
                        double LevSet_Deg = this.LevSet.Basis.Degree + 1;
                        h /= LevSet_Deg;
                        double dt_sigma = Math.Sqrt((this.Control.PhysicalParameters.rho_A + this.Control.PhysicalParameters.rho_B)
                            * Math.Pow(h, 3) / (2 * Math.PI * this.Control.PhysicalParameters.Sigma));
                        if(dt > dt_sigma)
                            Console.WriteLine("Warning: exceeding Capillary timestep: dt = {0}, dt_sigma = {1}, frac = {2}", dt, dt_sigma, dt / dt_sigma);
                    }


                    // elo
                    // ---

                    Console.WriteLine("Instationary solve, timestep #{0}, dt = {1} ...", TimestepNo, dt);

                } else {
                    throw new NotImplementedException("Option " + base.Control.TimesteppingMode + " not supported yet.");
                }

                // =======================================================================
                // call timestepper
                // =======================================================================

                //if ((m_BDF_Timestepper == null) == (m_RK_Timestepper == null))
                //    throw new ApplicationException();

                //CurvatureAlgorithms.CurvatureDriver(
                //    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                //    CurvatureAlgorithms.FilterConfiguration.NoFilter,
                //    this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                //    this.m_HMForder, this.DGLevSet.Current);

                //double[] momBal_Norm = XNSEUtils.MomentumBalanceNormAtInterface(this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                //    this.Control.PhysicalParameters, this.Control.AdvancedDiscretizationOptions.SurfStressTensor, this.m_HMForder);

                //Console.WriteLine("x-momentum balance norm = {0}", momBal_Norm[0]);
                //Console.WriteLine("y-momentum balance norm = {0}", momBal_Norm[1]);


                // ++++++++++++++++++++++++++++++++++++++++++
                // The actual solution of the System
                // ++++++++++++++++++++++++++++++++++++++++++

                using (new BlockTrace("Solve", tr)) {

                    if(m_BDF_Timestepper != null) {

                        updateSolutionParams.SetAll(true);
                        lockUpdate = false;

                        m_BDF_Timestepper.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
                        
                    } else {
                        //m_RK_Timestepper.Solve(phystime, dt);
                    }

                }


                if (this.Control.solveCoupledHeatEquation && (this.Control.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP))
                    this.ComputeHeatFlux();


                Postprocessing(TimestepInt, phystime, dt, TimestepNo);



                // ================
                // Good bye
                // ================
                if(this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
                    ExtVelMover.FinishTimeStep();
                }

#if DEBUG
            // in case of Debugging Save first Timesteps
            //if(TimestepNo[1] <= 2) {
            //    this.SaveToDatabase(TimestepNo, phystime);
            //}
#endif

                Console.WriteLine("done.");
                return dt;
            }
        }


        protected override void SetInitial() {
            base.SetInitial();

            this.InitLevelSet();

            this.CreateEquationsAndSolvers(null);

            // =========================================
            // XDG BDF Timestepper initialization
            // =========================================

            if (m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DelayedTimestepperInit(0.0, 0, this.Control.GetFixedTimestep(),
                    // delegate for the initialization of previous timesteps from an analytic solution
                    BDFDelayedInitSetIntial);
            }

            After_SetInitialOrLoadRestart(0.0, 0);

        }


        protected override void ResetInitial() {
            base.SetInitial();
            this.InitLevelSet();

            if (this.Control.solveCoupledHeatEquation) {
                if (this.Control.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP)
                    m_BDF_Timestepper.ResetDataAfterBalancing(ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure, this.Temperature));
                else
                    m_BDF_Timestepper.ResetDataAfterBalancing(ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure, this.Temperature, this.HeatFlux));
            } else {
                m_BDF_Timestepper.ResetDataAfterBalancing(ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray(), this.Pressure));
            }

            m_BDF_Timestepper.DelayedTimestepperInit(0.0, 0, this.Control.GetFixedTimestep(),
                // delegate for the initialization of previous timesteps from an analytic solution
                BDFDelayedInitSetIntial);

            //if (this.Control.solveCoupledHeatEquation) {
            //    if (this.Control.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP)
            //        m_BDF_coupledTimestepper.ResetDataAfterBalancing(this.Temperature.ToEnumerable());
            //    else
            //        m_BDF_coupledTimestepper.ResetDataAfterBalancing(ArrayTools.Cat<DGField>(this.Temperature.ToEnumerable(), this.HeatFlux.ToArray()));

            //    m_BDF_coupledTimestepper.SingleInit();
            //}

        }


        private void After_SetInitialOrLoadRestart(double PhysTime, int TimestepNo)
        {

            // =============================================
            // LogFile initialization
            // =============================================  

            if (this.Control.TestMode == true)
            {
                LogQueryValue(PhysTime);
            }
            else
            {
                if (this.Control.LogValues != XNSE_Control.LoggingValues.None && this.CurrentSessionInfo.ID != Guid.Empty && base.MPIRank == 0)
                {
                    InitLogFile(this.CurrentSessionInfo.ID);
                    WriteLogLine(TimestepNo, PhysTime);
                }
            }
        }


        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            base.LoadRestart(out Time, out TimestepNo);

            //this.InitLevelSet();
            if(this.Control.ReInitPeriod > 0) {
                Console.WriteLine("FastMarchReInit performing FirstOrderReInit");
                FastMarchReinitSolver = new FastMarchReinit(DGLevSet.Current.Basis);
                CellMask Accepted = LsTrk.Regions.GetCutCellMask();
                CellMask ActiveField = LsTrk.Regions.GetNearFieldMask(1);
                CellMask NegativeField = LsTrk.Regions.GetSpeciesMask("A");
                FastMarchReinitSolver.FirstOrderReinit(DGLevSet.Current, Accepted, NegativeField, ActiveField);
            }

            ContinuityEnforcer = new ContinuityProjection(
                    ContBasis: this.LevSet.Basis,
                    DGBasis: this.DGLevSet.Current.Basis,
                    gridData: GridData,
                    Option: Control.LSContiProjectionMethod
                    );

            if(this.Control.Option_LevelSetEvolution == LevelSetEvolution.ExtensionVelocity) {
                ExtVelMover = new ExtensionVelocityBDFMover(LsTrk, DGLevSet.Current, DGLevSetGradient, new VectorField<DGField>(XDGvelocity.Velocity.ToArray()),
                    Control.EllipticExtVelAlgoControl, BcMap, bdfOrder, ExtensionVelocity.Current, new double[2] { Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B });
            }

            //this.LsTrk.UpdateTracker();

            this.CreateEquationsAndSolvers(null);

            // =========================================
            // XDG BDF Timestepper initialization
            // =========================================

            if (m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DelayedTimestepperInit(Time, TimestepNo.MajorNumber, this.Control.GetFixedTimestep(),
                    // delegate for the initialization of previous timesteps from restart session
                    BDFDelayedInitLoadRestart );
            }

            After_SetInitialOrLoadRestart(Time, TimestepNo.MajorNumber);

        }


        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

           
            // Load the sample Points for the restart of the Fourier LevelSet
            if (this.Control.FourierLevSetControl != null) {

                Guid sessionGuid = this.Control.RestartInfo.Item1;
                string FourierLog_path = this.Control.DbPath;
                FourierLog_path = FourierLog_path + "/sessions/" + sessionGuid + "/Log_FourierLS.txt";

                IList<Guid> spUids = new List<Guid>();
                try {
                    using (StreamReader samplPLogReader = new StreamReader(FourierLog_path)) {

                        while (!samplPLogReader.EndOfStream) {
                            spUids.Add(Guid.Parse(samplPLogReader.ReadLine()));
                        }
                    }
                } catch (FileNotFoundException) {
                    spUids = new Guid[0];
                }
                Guid[] samplPUids = spUids.ToArray();

                TimestepNumber tsNmbr = this.Control.RestartInfo.Item2;
                Guid samplPUid_restart = Guid.Empty;
                if (tsNmbr == null) {
                    samplPUid_restart = samplPUids[samplPUids.Length - 1];
                } else {
                    samplPUid_restart = samplPUids[tsNmbr.MajorNumber];
                }
                Partitioning p = null;
                double[] samplP = base.DatabaseDriver.LoadVector<double>(samplPUid_restart, ref p).ToArray();

                if (this.Control.FourierLevSetControl.FType == FourierType.Polar) {
                    double[] center = samplP.GetSubVector(0, 2);
                    samplP = samplP.GetSubVector(2, samplP.Length - 2);
                    this.Control.FourierLevSetControl.center = center;
                    this.Control.FourierLevSetControl.samplP = samplP;
                } else {
                    this.Control.FourierLevSetControl.samplP = samplP;
                }
            }

        }


        protected override void Bye() {
            base.Bye();
            if (EnergyLogger != null)
                EnergyLogger.Close();
        }

        #endregion



        //==========================
        // adaptive mesh refinement
        //==========================
        #region AMR

        CellMask NScm;

        CellMask NSbuffer;

        /// <summary>
        /// refinement indicator for a constant near band refinement
        /// </summary>
        int LevelIndicator(int j, int CurrentLevel) {

            if(this.Control.BaseRefinementLevel == 0)
                return 0;

            CellMask ccm = this.LsTrk.Regions.GetCutCellMask();
            CellMask near = this.LsTrk.Regions.GetNearFieldMask(1);
            CellMask nearBnd = near.AllNeighbourCells();
            CellMask buffer = nearBnd.AllNeighbourCells().Union(nearBnd).Except(near);


            int DesiredLevel_j = CurrentLevel;

            if(near.Contains(j)) {
                if(CurrentLevel < this.Control.BaseRefinementLevel) {
                    DesiredLevel_j++;
                } else {
                    // additional refinement
                    switch(this.Control.RefineStrategy) {
                        case XNSE_Control.RefinementStrategy.CurvatureRefined: {
                                double curv_max = 1.0 / (2.0 * ((GridData)this.GridData).Cells.h_min[j]);
                                double mean_curv = Math.Abs(this.Curvature.GetMeanValue(j));
                                double minCurv, maxCurv;
                                this.Curvature.GetExtremalValuesInCell(out minCurv, out maxCurv, j);
                                double max_AbsCurv = Math.Max(Math.Abs(minCurv), Math.Abs(maxCurv));

                                double curv_thrshld = mean_curv;
                                if(curv_thrshld > curv_max && CurrentLevel == this.Control.RefinementLevel) {
                                    DesiredLevel_j++;
                                } else if(curv_thrshld < (curv_max / 2) && CurrentLevel == this.Control.RefinementLevel + 1) {
                                    DesiredLevel_j--;
                                }
                                break;
                            }
                        case XNSE_Control.RefinementStrategy.ContactLineRefined: {
                                CellMask BCells = ((GridData)this.GridData).BoundaryCells.VolumeMask;
                                if(ccm.Contains(j) && BCells.Contains(j) && CurrentLevel < this.Control.RefinementLevel) {
                                    DesiredLevel_j++;
                                } else if(!BCells.Contains(j)) { // && CurrentLevel == this.Control.RefinementLevel + 1) {
                                    DesiredLevel_j--;
                                }
                                break;
                            }
                        case XNSE_Control.RefinementStrategy.constantInterface:
                        default:
                            break;
                    }
                }

            } else if(NScm.Contains(j)) {
                if(CurrentLevel < this.Control.BaseRefinementLevel)
                    DesiredLevel_j++;

            } else if(buffer.Contains(j) || NSbuffer.Contains(j)) {
                if(CurrentLevel < this.Control.BaseRefinementLevel - 1)
                    DesiredLevel_j++;
            } else {
                DesiredLevel_j = 0;
            }

            return DesiredLevel_j;

        }


        //int LevelIndicator(int j, int CurrentLevel) {

        //    CellMask spc = this.LsTrk.Regions.GetSpeciesMask("A");

        //    if(spc.Contains(j)) {
        //        return 2;
        //    } else {
        //        return 0;
        //    }

        //}


        /// <summary>
        /// refinement indicator
        /// </summary>
        //int LevelIndicator(int j, int CurrentLevel) {

        //    int minRefineLevelLS = 1;
        //    int maxRefineLevelLS = 2;

        //    CellMask ccm = this.LsTrk.Regions.GetCutCellMask();
        //    CellMask near = this.LsTrk.Regions.GetNearFieldMask(1);

        //    double curv_max = 1.0 / this.GridData.Cells.h_min[j];

        //    int DesiredLevel_j = CurrentLevel;

        //    if(near.Contains(j)) {

        //        if(DesiredLevel_j < minRefineLevelLS) {
        //            // set minimum refinement level for the interface
        //            DesiredLevel_j = minRefineLevelLS;

        //        } else if (ccm.Contains(j)) {
        //            // further localized refinement

        //            // check for high curvature
        //            //int DesiredLevelj_highCurv = DesiredLevel_j;
        //            //this.Curvature.GetExtremalValuesInCell(out double curv_jMin, out double curv_jMax, j);
        //            //if((curv_jMax >= curv_max || Math.Abs(curv_jMin) >= curv_max) && DesiredLevel_j < maxRefineLevelLS) {
        //            //    DesiredLevelj_highCurv++;
        //            //} else if((curv_jMax < curv_max / 2) || (Math.Abs(curv_jMin) < curv_max / 2)) {
        //            //    DesiredLevelj_highCurv--;
        //            //}

        //            double mean_curv = Math.Abs(this.Curvature.GetMeanValue(j));
        //            if((mean_curv >= curv_max) && CurrentLevel < maxRefineLevelLS) {
        //                DesiredLevel_j++;
        //            } else if(mean_curv < curv_max / 2 && CurrentLevel > minRefineLevelLS) {
        //                DesiredLevel_j--;
        //            }

        //            //// check for small cut cells
        //            //int DesiredLevelj_agglom = DesiredLevel_j;
        //            //double cellVol = this.GridData.Cells.GetCellVolume(j);
        //            //var spcIds = this.LsTrk.SpeciesIdS.ToArray();
        //            //double ratioVolSpcMin = 1.0;
        //            //foreach(SpeciesId spc in this.LsTrk.SpeciesIdS) {
        //            //    double cellVolSpc = this.LsTrk.GetXDGSpaceMetrics(spcIds, m_HMForder, 1).CutCellMetrics.CutCellVolumes[spc][j];
        //            //    double ratioVolSpc = cellVolSpc / cellVol;
        //            //    if(ratioVolSpc < ratioVolSpcMin)
        //            //        ratioVolSpcMin = ratioVolSpc;
        //            //}
        //            //double thrshld = this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold;
        //            //if(ratioVolSpcMin < thrshld && DesiredLevel_j < maxRefineLevelLS) {
        //            //    DesiredLevelj_agglom++;
        //            //} else if(ratioVolSpcMin > 4 * thrshld) {
        //            //    DesiredLevelj_agglom--;
        //            //}

        //            //// check for a change of sign in the curvature
        //            //int DesiredLevelj_inflection = DesiredLevel_j;
        //            ////this.Curvature.GetExtremalValuesInCell(out double curv_jMin, out double curv_jMax, j);
        //            //if(Math.Sign(curv_jMin) != Math.Sign(curv_jMax) && DesiredLevel_j < maxRefineLevelLS)
        //            //    DesiredLevelj_inflection++;

        //            //DesiredLevel_j = (new int[] { DesiredLevelj_highCurv, DesiredLevelj_agglom, DesiredLevelj_inflection }).Max();

        //        }

        //    } else {
        //        // non cut cells don't need to be refined
        //        DesiredLevel_j = 0;
        //    }

        //    return DesiredLevel_j;

        //}

        //CellMask refinedInterfaceCells;

        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            using(new FuncTrace()) {

                if(this.Control.AdaptiveMeshRefinement) {

                    //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 0), 2);

                    // Check grid changes
                    // ==================

                    CellMask BlockedCells;
                    if(this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Once
                        || this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                        int prevInd = LsTrk.PopulatedHistoryLength;
                        CellMask prevNear = LsTrk.RegionsHistory[-prevInd + 1].GetNearFieldMask(1);
                        BlockedCells = (TimestepNo > 1) ? prevNear : null; // CellMask.Union(currNear, prevNear);
                    } else {
                        CellMask currNear = LsTrk.Regions.GetNearFieldMask(1);
                        BlockedCells = currNear;
                    }

                    // compute curvature for levelindicator 
                    if(this.Control.RefineStrategy == XNSE_Control.RefinementStrategy.CurvatureRefined) {
                        CurvatureAlgorithms.CurvatureDriver(
                            SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint,
                            CurvatureAlgorithms.FilterConfiguration.Default,
                            this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                            this.m_HMForder, this.DGLevSet.Current);
                    }


                    // navier slip boundary cells
                    NScm = new CellMask(this.GridData);
                    NSbuffer = new CellMask(this.GridData);
                    if(this.Control.RefineNavierSlipBoundary) {
                        BitArray NSc = new BitArray(((GridData)this.GridData).Cells.Count);
                        CellMask bnd = ((GridData)this.GridData).BoundaryCells.VolumeMask;
                        int[][] c2e = ((GridData)this.GridData).Cells.Cells2Edges;
                        foreach(Chunk cnk in bnd) {
                            for(int i = cnk.i0; i < cnk.JE; i++) {
                                foreach(int e in c2e[i]) {
                                    int eId = (e < 0) ? -e - 1 : e - 1;
                                    byte et = ((GridData)this.GridData).Edges.EdgeTags[eId];
                                    if(this.GridData.EdgeTagNames[et].Contains("navierslip_linear"))
                                        NSc[i] = true;
                                }
                            }
                        }
                        NScm = new CellMask(this.GridData, NSc);
                        CellMask bndNScm = NScm.AllNeighbourCells();
                        int bndLvl = 2;
                        for(int lvl = 1; lvl < bndLvl; lvl++) {
                            NScm = NScm.Union(bndNScm);
                            bndNScm = NScm.AllNeighbourCells();
                            NSbuffer = NScm.Union(bndNScm);
                        }
                        NSbuffer = NSbuffer.AllNeighbourCells();
                    }


                    //PlotCurrentState(hack_Phystime, new TimestepNumber(TimestepNo, 1), 2);


                bool AnyChange = GridRefinementController.ComputeGridChange((BoSSS.Foundation.Grid.Classic.GridData) this.GridData, BlockedCells, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                int NoOfCellsToRefine = 0;
                int NoOfCellsToCoarsen = 0;
                if (AnyChange) {
                    int[] glb = (new int[] {

                    CellsToRefineList.Count,
                    Coarsening.Sum(L => L.Length),
                }).MPISum();
                        NoOfCellsToRefine = glb[0];
                        NoOfCellsToCoarsen = glb[1];
                    }
                    int oldJ = this.GridData.CellPartitioning.TotalLength;

                    // Update Grid
                    // ===========

                    if(AnyChange) {

                        //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 1 }), 2);

                        Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                        Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                        newGrid = ((GridData)this.GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                        //PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, 2 }), 2);


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


        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
            //if(this.Control.solveCoupledHeatEquation)
            //    m_BDF_coupledTimestepper.DataBackupBeforeBalancing(L);
        }


        #endregion



        //===========================
        // I/O (saving and plotting)
        //===========================
        #region IO


        /// <summary>
        /// 
        /// </summary>
        protected override ITimestepInfo SaveToDatabase(TimestepNumber timestepno, double t) {
            var tsi = base.SaveToDatabase(timestepno, t);

            if (tsi != null && m_BDF_Timestepper != null) {
                int S = m_BDF_Timestepper.GetNumberOfStages;

                SinglePhaseField LsBkUp = new SinglePhaseField(this.LevSet.Basis);
                LsBkUp.Acc(1.0, this.LevSet);

                ICollection<DGField>[] restartFields = m_BDF_Timestepper.GetRestartInfos();

                if (S > 1 && this.Control.saveperiod >= S && restartFields != null) {

                    // save additional timesteps/information for restart
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    for (int ti = 1; ti < S; ti++) {

                        //SinglePhaseField LsBkUp = new SinglePhaseField(this.LevSet.Basis);
                        //LsBkUp.Acc(1.0, this.LevSet);

                        ICollection<DGField> restartIOFields = new List<DGField>();
                        foreach (DGField f in this.IOFields) {

                            int rfidx = restartFields[ti - 1].IndexWhere(rf => rf.Identification == f.Identification);
                            if (rfidx > -1) {
                                DGField rf = restartFields[ti - 1].ElementAt(rfidx);
                                if (f.Identification == "Phi") {
                                    this.LevSet.Clear();
                                    this.LevSet.Acc(1.0, rf);
                                    restartIOFields.Add(this.LevSet);
                                } else {
                                    restartIOFields.Add(rf);
                                }
                            } else {
                                DGField rf = f.CloneAs();
                                rf.Clear();
                                restartIOFields.Add(rf);
                            }
                        }

                        //this.LevSet.Clear();
                        //this.LevSet.Acc(1.0, LsBkUp);

                        ITimestepInfo rtsi;
                        TimestepNumber tsn = new TimestepNumber(timestepno.MajorNumber - ti);

                        //Console.WriteLine("saving to Database ...");

                        //Exception e = null;
                        try {
                            rtsi = new TimestepInfo(
                                t - ti * this.Control.GetFixedTimestep(),
                                this.CurrentSessionInfo,
                                tsn,
                                restartIOFields);
                        } catch (Exception ee) {
                            Console.Error.WriteLine(ee.GetType().Name + " on rank " + this.MPIRank + " saving time-step " + tsn + ": " + ee.Message);
                            Console.Error.WriteLine(ee.StackTrace);
                            //tsi = null;
                            //e = ee;

                            if (ContinueOnIOError) {
                                Console.WriteLine("Ignoring IO error: " + DateTime.Now);
                            } else {
                                throw ee;
                            }

                            tsi = null;
                        }
                        // e.ExceptionBcast();
                        csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                    }

                }

                this.LevSet.Clear();
                this.LevSet.Acc(1.0, LsBkUp);
            }

#if DEBUG
            //Debug/Test code for XDG database interaction

            if(tsi != null) {
                // checking some neccessary reference-equalities BEFORE serialisation
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                LevelSet.LevelSetInitializer lsi_1 = (LevelSet.LevelSetInitializer)(tsi.FieldInitializers.Single(fi => fi.Identification == this.LevSet.Identification));
                XDGField.FieldInitializer pri = (XDGField.FieldInitializer)(tsi.FieldInitializers.Single(fi => fi.Identification == this.Pressure.Identification));

                LevelSetTracker.LevelSetTrackerInitializer trki = ((XDGBasis.XDGBasisInitializer)(pri.BasisInfo)).TrackerInitializer;

                LevelSet.LevelSetInitializer lsi_2 = trki.LevelSets[0];

                Debug.Assert(object.ReferenceEquals(lsi_1, lsi_2));

                foreach(XDGField.FieldInitializer fi in tsi.FieldInitializers.Where(fii => fii is XDGField.XDGFieldInitializer)) {
                    LevelSetTracker.LevelSetTrackerInitializer trki_alt = ((XDGBasis.XDGBasisInitializer)(fi.BasisInfo)).TrackerInitializer;
                    Debug.Assert(object.ReferenceEquals(trki, trki_alt));
                }
            }


            if(tsi != null) {
                // checking some neccessary equalities AFTER serialisation
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                

                var tsi_alt = this.DatabaseDriver.LoadTimestepInfo(tsi.ID, base.CurrentSessionInfo, base.GetDatabase());

                Debug.Assert(!object.ReferenceEquals(tsi, tsi_alt));


                LevelSet.LevelSetInitializer lsi_1 = (LevelSet.LevelSetInitializer)(tsi_alt.FieldInitializers.Single(fi => fi.Identification == "Phi"));
                XDGField.FieldInitializer pri = (XDGField.FieldInitializer)(tsi_alt.FieldInitializers.Single(fi => fi.Identification == this.Pressure.Identification));

                LevelSetTracker.LevelSetTrackerInitializer trki = ((XDGBasis.XDGBasisInitializer)(pri.BasisInfo)).TrackerInitializer;

                LevelSet.LevelSetInitializer lsi_2 = trki.LevelSets[0];

                Debug.Assert(lsi_1.Equals(lsi_2));

                foreach(XDGField.FieldInitializer fi in tsi_alt.FieldInitializers.Where(fii => fii is XDGField.XDGFieldInitializer)) {
                    LevelSetTracker.LevelSetTrackerInitializer trki_alt = ((XDGBasis.XDGBasisInitializer)(fi.BasisInfo)).TrackerInitializer;
                    Debug.Assert(trki.Equals(trki_alt));
                }


                var Fields = this.DatabaseDriver.LoadFields(tsi_alt, this.GridData);
                LevelSet Rphi_1 = (LevelSet)(Fields.Single(f => f.Identification.Equals(this.LevSet.Identification)));

                XDGField Rpressure = (XDGField)(Fields.Single(f => f.Identification.Equals(this.Pressure.Identification)));

                LevelSetTracker Rtracker = Rpressure.Basis.Tracker;
                Debug.Assert(!object.ReferenceEquals(this.LsTrk, Rtracker));
                Debug.Assert(object.ReferenceEquals(Rtracker.LevelSets[0], Rphi_1));

                foreach(XDGField xf in Fields.Where(fii => fii is XDGField)) {
                    Debug.Assert(object.ReferenceEquals(xf.Basis.Tracker, Rtracker));
                }
            }

#endif

            return tsi;
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {
            Tecplot.PlotFields(base.m_RegisteredFields, "XNSE_Solver" + timestepNo, physTime, superSampling);
            //Tecplot.PlotFields(new DGField[] { this.LevSet }, "grid" + timestepNo, physTime, 0);
        }


        protected void PlotOnIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            PlotCurrentState(hack_Phystime, new TimestepNumber(new int[] { hack_TimestepIndex, iterIndex }), 2);
        }

        #endregion

        /// <summary>
        /// makes direct use of <see cref="XdgTimesteppingBase.OperatorAnalysis"/>; aids the condition number scaling analysis
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis() {
            return this.m_BDF_Timestepper.OperatorAnalysis();
        }

    }
}
