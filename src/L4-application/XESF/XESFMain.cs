using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Residual;
using BoSSS.Solution.Statistic;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using XESF.Fluxes;
using XESF.Variables;
using ApplicationWithIDT;
using BoSSS.Solution.CompressibleFlowCommon.ShockFinding;
using ApplicationWithIDT.OptiLevelSets;
using BoSSS.Solution.Statistic.QuadRules;


namespace XESF
{
    /// <summary>
    /// Implements shock fitting for stationary XDG Euler (2D) which is solved by the routines defined in <see cref="ApplicationWithIDT"/>
    /// Naming: X(DG) - E(uler) - S(hock) - F(itting)
    /// Concrete configurations of solver (Initial Guess, optimization parameters,...) are set in a <see cref="XESFControl.cs"/> object
    /// Fluxes are implemented in <see cref="XESF.Fluxes"/>. This solver supports:
    /// - Roe Flux (with smoothing factor alpha)
    /// - HLLC Flux
    /// - Godunov Flux
    /// - Central Flux (only for interface)
    /// 
    /// Author: Jakob Vandergrift 
    /// Date of Creation/Maintenance: 08-2022 until at least 08-2024
    /// </summary>
    public class XESFMain : ApplicationWithIDT<XESFControl>
    {

        /// <summary>
        /// </summary>
        /// <param name="args">string pointing to a control file, i.e. 'cs:XESF.XESFHardCodedControl.XDGWedgeFlow_TwoLs_Base()' </param>
        static void Main(string[] args)
        {
            //IDTTestRunner.RunTests();
            //XESF.Tests.XESFTestProgram.XDGBowShockFromOldRun(tsNumber:161);
            //XESF.Tests.XESFTestProgram.XDGBowShockFromDB(5, 16, 1, 0);
            XESFMain._Main(args, false, () => new XESFMain());
        }

        #region XDG Fields
        /// <summary>
        /// The density $\rho$
        /// </summary>
        public XDGField Density
        {
            get;
            private set;
        }

        /// <summary>
        /// The momentum field $\rho \vec{u}$
        /// </summary>
        public VectorField<XDGField> Momentum
        {
            get;
            private set;
        }

        /// <summary>
        /// The total energy per volume $\rho E$
        /// </summary>
        public XDGField Energy
        {
            get;
            private set;
        }
        /// <summary>
        /// All conservative fields (density, momentum vector, energy)
        /// </summary>

        #endregion

        /// <summary>
        /// Initial Guess Level Set loaded from DB (i.e. from Reconstruction) 
        /// </summary>
        public SinglePhaseField ShockLevelSetField;

        /// <summary>
        /// Time Step serving as initial guess
        /// </summary>
        public ITimestepInfo tsiFromDb;

        /// <summary>
        /// residual Logger
        /// </summary>
        public IResidualLogger[] residualLoggers;

        /// <summary>
        /// This is a map between quantities of interest or variables stored as XDGFields ( enthalpy, pressure, ...)
        /// </summary>
        public Dictionary<DerivedVariable<XDGField>, XDGField> DerivedVariableToXDGFieldMap
        {
            get;
            private set;
        }
        /// <summary>
        /// This is a map between quantities of interest or variables stored as Single Phase Fields ( PerssonSensor, Indicators, etc. ...)
        /// </summary>
        public Dictionary<DerivedVariable<SinglePhaseField>, SinglePhaseField> DerivedVariableToSinglePhaseFieldMap
        {
            get;
            private set;
        }
        /// <summary>
        /// This is a map between quantities of interest or variables stored as doubles (error)
        /// </summary>
        public Dictionary<DerivedVariable<double>, double> DerivedVariableToDoubleMap
        {
            get;
            private set;
        }
        /// <summary>
        /// creates or loads the grid
        /// </summary>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        protected override IGrid CreateOrLoadGrid()
        {
            IGrid grid = null;

            // Optional: load grid from CNS calculation from which the shock level set has been reconstructed
            if (Control.ShockLevelSet_UseShockGrid)
            {
                IDatabaseInfo dbi = DatabaseInfo.Open(Control.ShockLevelSet_Db);
                ISessionInfo si = dbi.Controller.GetSessionInfo(Control.ShockLevelSet_Info.Item1);
                ITimestepInfo tsi;

                if (Control.ShockLevelSet_Info.Item2.MajorNumber < 0.0)
                {
                    tsi = si.Timesteps.Last();
                }
                else
                {
                    throw new NotSupportedException("Currently, only the last time step can be used.");
                }

                grid = dbi.Controller.DBDriver.LoadGrid(tsi.GridID, dbi);

                Console.WriteLine("Using grid from LEVEL SET RECONSTRUCTION");
            }
            else
            {
                grid = base.CreateOrLoadGrid();
            }

            CompressibleEnvironment.Initialize(grid.SpatialDimension);

            return grid;
        }
        /// <summary>
        /// Initializes the Optimization Problem and computes and initial residual, also if chosen initial value is computed by P0 projection
        /// 1. Here the objects for the discretized equation are constructed by constructing a Spatial Operator and adding the chosen fluxes defined in XESF.FLUXES
        /// 2. The Functional needed for the objective function is constructed
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L)
        {
            #region Init operator

            GridData gridData = (GridData)this.GridData;
            Material material = this.Control.GetMaterial();
            IBoundaryConditionMap boundaryMap = new XDGCompressibleBoundaryCondMap(this.GridData, this.Control, material, this.Control.SpeciesToEvaluate);

            string[] variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };


            this.XSpatialOperator = new XDifferentialOperatorMk2(variables, null, variables, Control.quadOrderFunc, this.SpeciesToEvaluate);
            this.Op_obj = new XDifferentialOperatorMk2(variables, null, variables, Control.quadOrderFunc, this.SpeciesToEvaluate);
            #endregion



            switch (Control.ConvectiveBulkFlux)
            {
                case ConvectiveBulkFluxes.OptimizedHLLC:
                    switch (Control.FluxVersion)
                    {
                        case FluxVersion.NonOptimized:
                            foreach (SpeciesId id in SpeciesToEvaluate_Ids)
                            {
                                string spcNmn = LsTrk.GetSpeciesName(id);
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.HLLCDensityFlux_XDG(boundaryMap, material, spcNmn, id));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.HLLCMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.HLLCMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.HLLCEnergyFlux_XDG(boundaryMap, material, spcNmn, id));
                            }
                            break;
                        case FluxVersion.Optimized:
                            foreach (SpeciesId id in SpeciesToEvaluate_Ids)
                            {
                                string spcNmn = LsTrk.GetSpeciesName(id);
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG(boundaryMap, material, spcNmn, id));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG(boundaryMap, material, spcNmn, id));
                            }
                            break;
                    }
                    break;
                case ConvectiveBulkFluxes.Roe:
                    foreach (SpeciesId id in SpeciesToEvaluate_Ids)
                    {
                        string spcNmn = LsTrk.GetSpeciesName(id);
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.RoeDensityFlux_XDG(boundaryMap, material, spcNmn, id, Control.flux_s_alpha));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.RoeMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id, Control.flux_s_alpha));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.RoeMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id, Control.flux_s_alpha));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.RoeEnergyFlux_XDG(boundaryMap, material, spcNmn, id, Control.flux_s_alpha));
                    }
                    break;
                case ConvectiveBulkFluxes.Godunov:
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new GodunovFlux(this.Control, boundaryMap, new EulerDensityComponent(), material));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new GodunovFlux(this.Control, boundaryMap, new EulerMomentumComponent(0, material.EquationOfState.HeatCapacityRatio, Control.MachNumber, gridData.SpatialDimension), material));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new GodunovFlux(this.Control, boundaryMap, new EulerMomentumComponent(1, material.EquationOfState.HeatCapacityRatio, Control.MachNumber, gridData.SpatialDimension), material));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new GodunovFlux(this.Control, boundaryMap, new EulerEnergyComponent(), material));
                    break;
                case ConvectiveBulkFluxes.CentralFlux:
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new CentralDensityFlux(boundaryMap, material));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new CentralMomentumFlux(boundaryMap, material,0));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new CentralMomentumFlux(boundaryMap, material, 1));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new CentralEnergyFlux(boundaryMap, material));
                    break;
                default:
                    
                    throw new NotImplementedException("The convectiveBulkFlux you chose is not implemented");
            }

            switch (Control.ConvectiveInterfaceFlux_LsOne)
            {

                case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                    switch (Control.FluxVersion)
                    {
                        case FluxVersion.NonOptimized:
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;
                        case FluxVersion.Optimized:
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;
                    }
                    break;
                case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                case ConvectiveInterfaceFluxes.OptimizedHLLCWall:
                    switch (Control.FluxVersion)
                    {
                        case FluxVersion.NonOptimized:
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;
                        case FluxVersion.Optimized:
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;
                    }
                    break;
                case ConvectiveInterfaceFluxes.CentralFluxInterface:
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    break;

                case ConvectiveInterfaceFluxes.RoeInterface:
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, int.MinValue, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, int.MinValue, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                    break;

                case ConvectiveInterfaceFluxes.RoeWall:
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RoeFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, int.MinValue, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RoeFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RoeFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RoeFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, int.MinValue, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies, s_alpha: Control.flux_s_alpha));
                    break;
                case ConvectiveInterfaceFluxes.GodunovInterface:
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                    break;

                case ConvectiveInterfaceFluxes.None:
                    break;

                default:
                    throw new NotSupportedException("This should never happen");
            }

            if (Control.IsTwoLevelSetRun)
            {
                switch (Control.ConvectiveInterfaceFlux_LsTwo)
                {
                    case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                        switch (Control.FluxVersion)
                        {
                            case FluxVersion.NonOptimized:
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                break;
                            case FluxVersion.Optimized:
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                break;
                        }
                        break;
                    case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                        switch (Control.FluxVersion)
                        {
                            case FluxVersion.NonOptimized:
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                break;
                            case FluxVersion.Optimized:
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                break;
                        }
                        break;

                    case ConvectiveInterfaceFluxes.CentralFluxInterface:
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));

                        break;

                    case ConvectiveInterfaceFluxes.RoeInterface:
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, Int16.MinValue, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, Int16.MinValue, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies, s_alpha: Control.flux_s_alpha, D: 2));
                        break;
                    case ConvectiveInterfaceFluxes.GodunovInterface:
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                        break;


                    case ConvectiveInterfaceFluxes.None:
                        break;

                    default:
                        throw new NotSupportedException("This should never happen");
                }
            }
            if (Control.FluxVersion == FluxVersion.NonOptimized)
            {
                this.XSpatialOperator.LinearizationHint = LinearizationHint.GetJacobiOperator;
                this.XSpatialOperator.AgglomerationThreshold = this.Control.AgglomerationThreshold;
                this.XSpatialOperator.Commit();

                r_JacobiOperator = XSpatialOperator.GetJacobiOperator(SpatialDimension: 2);
                R_JacobiOperator = XSpatialOperator.GetJacobiOperator(SpatialDimension: 2);

            }
            else
            {
                if (Control.ConvectiveBulkFlux != ConvectiveBulkFluxes.None || Control.ConvectiveInterfaceFlux_LsOne != ConvectiveInterfaceFluxes.None)
                {
                    this.XSpatialOperator.LinearizationHint = LinearizationHint.FDJacobi;
                    this.Op_obj.LinearizationHint = LinearizationHint.FDJacobi;
                }
                else
                {
                    this.XSpatialOperator.LinearizationHint = LinearizationHint.AdHoc;
                    this.Op_obj.LinearizationHint = LinearizationHint.AdHoc;
                }
                this.XSpatialOperator.AgglomerationThreshold = this.Control.AgglomerationThreshold;
                this.Op_obj.AgglomerationThreshold = this.Control.AgglomerationThreshold;

                // here the Optimization problem is set.
                // - Rankine Hugoniot sets and objective function that is the weak form of either evaluating the RH condition on the interface only or on all faces.
                // - EnRes uses the objective function defined by the enriched Residual. The Obj Operator is then the same
                switch (Control.optProblemType)
                {
                    case OptProblemType.RankineHugoniotFull:
                    case OptProblemType.RankineHugoniotOnlyInterface:
                        if (Control.optProblemType == OptProblemType.RankineHugoniotFull)
                        { // add a Bulk Flux 
                            foreach (SpeciesId id in SpeciesToEvaluate_Ids)
                            {
                                string spcNmn = LsTrk.GetSpeciesName(id);
                                this.Op_obj.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RHDensityFlux_XDG(boundaryMap, material, spcNmn, id));
                                this.Op_obj.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RHMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id));
                                this.Op_obj.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RHMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id));
                                this.Op_obj.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RHEnergyFlux_XDG(boundaryMap, material, spcNmn, id));
                            }
                        }
                        switch (Control.IsTwoLevelSetRun)
                        {
                            case true:
                                this.Op_obj.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.Op_obj.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.Op_obj.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                this.Op_obj.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                                break;
                            case false:
                                this.Op_obj.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                                this.Op_obj.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                                this.Op_obj.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                                this.Op_obj.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RHFlux_XDG_Interface(boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                                break;
                        }
                        break;
                    default:
                        Op_obj = XSpatialOperator.CloneAs();
                        break;
                }
                this.Op_obj.Commit();
                this.XSpatialOperator.AgglomerationThreshold = this.Control.AgglomerationThreshold;
                this.XSpatialOperator.Commit();
            }
            if (Control.GetInitialValue == GetInitialValue.FromP0Timestepping)
            {
                ComputeP0Solution(); //needs to be done here as operator is now assembled
            }


            #region Residual logging
            // Configure residual handling
            if (L == null && Control.ResidualInterval != 0)
            {
                // Do not change these settings upon repartitioning
                ResLogger.WriteResidualsToTextFile = false;
                ResLogger.WriteResidualsToConsole = false;
            }

            residualLoggers = Control.ResidualLoggerType.Instantiate(
                this,
                Control,
                null,
                ConservativeFields,
                null).ToArray();
            UpdateDerivedVariables();
            //PlotCurrentState("debug",5);
            //Eval_R = XSpatialOperator.GetEvaluatorEx(LsTrk, ConservativeFields, null, obj_f_map);
            //Eval_R.Evaluate(1.0, 0.0, obj_f_vec);
            //if(CurrentAgglo > 0) {
            //    MultiphaseAgglomerator.ManipulateMatrixAndRHS(default(MsrMatrix), obj_f_vec, obj_f_map, ResidualMap);
            //    MultiphaseAgglomerator.ManipulateMatrixAndRHS(default(MsrMatrix), ResidualVector, ResidualMap, ResidualMap);
            //}
            //init OptProblem
            ChooseOptProblem();
            //Initialize empty vectors and matrices
            InitializeMatricesAndVectors();
            //// Cell agglomerator (cell length scales are needed for diffusive AV fluxes)
            UpdateAgglomerator();

            (res_l2, obj_f, res_L2) = ComputeResiduals();

            //obj_f = obj_f_vec.MPI_L2Norm();
            //Eval_r = XSpatialOperator.GetEvaluatorEx(LsTrk, ConservativeFields, null, ResidualMap);
            //Eval_r.Evaluate(1.0, 0.0, ResidualVector);
            //GetResFromEnRes(obj_f_vec, ResidualVector);
            //res_l2 = ResidualVector.MPI_L2Norm();
            InitResNorm = res_l2;
            Init_obj_f = obj_f;
            #endregion
            ResNorms.Add(res_l2);
            obj_f_vals.Add(obj_f);
            UpdateDerivedVariables();

        }

        /// <summary>
        /// Loads all fields in <see cref="m_IOFields"/> from the database
        /// using the given <see cref="AppControl.RestartInfo"/> (<see cref="Control"/>)
        /// </summary>
        /// <param name="time">
        /// On exit, contains the physical time represented by the time-step
        /// </param>
        /// <returns>
        /// Returns the actual time-step number of the loaded time-step
        /// </returns>
        protected override TimestepNumber RestartFromDatabase(out double time)
        {
            using (var tr = new FuncTrace())
            {

                // obtain session timesteps:
                var sessionToLoad = this.Control.RestartInfo.Item1;
                ISessionInfo session = m_Database.Controller.GetSessionInfo(sessionToLoad);

                var all_ts = session.Timesteps;

                // find timestep to load
                Guid tsi_toLoad_ID;
                tsi_toLoad_ID = GetRestartTimestepID();
                ITimestepInfo tsi_toLoad = all_ts.Single(t => t.ID.Equals(tsi_toLoad_ID));

                time = tsi_toLoad.PhysicalTime;


                if (tsi_toLoad is BoSSS.Foundation.IO.TimestepProxy tp)
                {
                    var tsiI = tp.GetInternal() as TimestepInfo;
                    if (tsiI != null)
                    {
                        OnRestartTimestepInfo(tsiI);
                    }
                }

                DatabaseDriver.LoadFieldData(tsi_toLoad, ((GridData)(this.GridData)), this.IOFields);

                if (((TimestepProxy)tsi_toLoad).GetInternal() is IDTTimeStepInfo IDTtsInfo)
                {
                    this.gamma = IDTtsInfo.Gamma;
                    CurrentStepNo = (int)IDTtsInfo.TimeStepNumbers.Last();
                    this.Gammas = IDTtsInfo.GammaHistory;
                    this.Alphas = IDTtsInfo.AlphaHistory.ToList();
                    this.ResNorms = IDTtsInfo.ResHistory.ToList();
                    if (LevelSetOpti is SplineOptiLevelSet splineLs)
                    {                       
                        var LSParams = IDTtsInfo.LevelSetParams;
                        if (LSParams.Length == splineLs.m_AllParams.Length)
                        {
                            for (int i = 0; i < LSParams.Length; i++)
                            {
                                splineLs.m_AllParams[i] = LSParams[i];
                            }
                            splineLs.Interpolate();
                            splineLs.ProjectOntoLevelSet(LsTBO);

                        }
                        else
                        {
                            throw new NotSupportedException($"{GetLevelSet.DirectyFromTimestep} not supported if Grid has not equal y - Cells");
                        }
                    }
                }
                if (Control.IsTwoLevelSetRun)
                {
                    this.LevelSet.ProjectField(1.0, this.Control.LevelSetOneInitialValue);
                    LsTBO = LevelSetTwo;
                }
                else
                {
                    // Level set one
                    LsTBO = LevelSet;
                }
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                LsTrk.UpdateTracker(CurrentStepNo);

                time = CurrentStepNo;
                // return
                return tsi_toLoad.TimeStepNumber;
            }
        }

        /// <summary>
        /// Initializes multi-grid operator configuration object. The MGP is used only in the process of agglomeration so far to transform from solution to agglomeration space.
        /// </summary>
        public override void InitializeMultiGridOpConfig()
        {
            int p = this.Momentum[0].Basis.Degree;
            int D = this.GridData.SpatialDimension;
            // set the MultigridOperator configuration for each level:
            // it is not necessary to have exactly as many configurations as actual multi grid levels:
            // the last configuration entry will be used for all higher level
            MultiGridOperatorConfig = new MultigridOperator.ChangeOfBasisConfig[3][];
            for (int iLevel = 0; iLevel < MultiGridOperatorConfig.Length; iLevel++)
            {
                MultiGridOperatorConfig[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 2];

                // configuration for continuity
                MultiGridOperatorConfig[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig()
                {
                    DegreeS = new int[] { Math.Max(0, p - iLevel) },
                    mode = MultigridOperator.Mode.Eye,
                    VarIndex = new int[] { 0 }
                };


                // configurations for momentum equation
                for (int d = 0; d < D; d++)
                {
                    MultiGridOperatorConfig[iLevel][d + 1] = new MultigridOperator.ChangeOfBasisConfig()
                    {
                        DegreeS = new int[] { Math.Max(0, p - iLevel) },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { d + 1 }
                    };
                }

                // configuration for energy
                MultiGridOperatorConfig[iLevel][D + 1] = new MultigridOperator.ChangeOfBasisConfig()
                {
                    DegreeS = new int[] { Math.Max(0, p - iLevel) },
                    mode = MultigridOperator.Mode.Eye,
                    VarIndex = new int[] { D + 1 }
                };
            }
        }

        /// <summary>
        /// sets initial value
        /// </summary>
        /// <param name="t"></param>
        /// <exception cref="NotSupportedException"></exception>
        /// <exception cref="NotImplementedException"></exception>
        protected override void SetInitial(double t)
        {
            if (LsTrk.GridDat.SpatialDimension > 2)
            {
                throw new NotSupportedException("Only valid for 2D simulation runs.");
            }
            CurrentStepNo = 0;
            CellQuadratureScheme scheme = new CellQuadratureScheme(true, CellMask.GetFullMask(this.GridData));
            #region DB Choice
            // Switch to Load the right DB if a DB is loaded
            switch (Control.GetInitialValue)
            {
                case GetInitialValue.FromDBSinglePhase:
                case GetInitialValue.FromDBXDG:
                    DatabaseInfo dbi = DatabaseInfo.Open(Control.ShockLevelSet_Db);
                    ISessionInfo si = dbi.Controller.GetSessionInfo(Control.ShockLevelSet_Info.Item1);
                    tsiFromDb = Control.IVTimestepNumber >= si.Timesteps.Count ? si.Timesteps.Last() : si.Timesteps.Pick(Control.IVTimestepNumber);
                    break;
                case GetInitialValue.FromAVRun:
                    DatabaseInfo dbi_2 = DatabaseInfo.Open(Control.SeedFromAV_Db);
                    ISessionInfo si_2 = dbi_2.Controller.GetSessionInfo(Control.SeedFromAV_Db_Info.Item1);
                    tsiFromDb = Control.IVTimestepNumber >= si_2.Timesteps.Count ? si_2.Timesteps.Last() : si_2.Timesteps.Pick(Control.IVTimestepNumber);
                    break;
                default:
                    break;

            }
            #endregion

            #region Initialize the LevelSet
            if (Control.IsTwoLevelSetRun)
            {
                this.LevelSet.ProjectField(1.0, this.Control.LevelSetOneInitialValue, scheme);
                LsTBO = LevelSetTwo;
            }
            else
            {
                // Level set one
                LsTBO = LevelSet;
            }

            switch (Control.GetLevelSet)
            {

                case GetLevelSet.FromParams:
                    LevelSetOpti.AssembleTransMat(LsTBO);
                    LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                    break;

                case GetLevelSet.FromFunction:
                    LevelSetOpti.AssembleTransMat(LsTBO);
                    LevelSetOpti.ProjectFromFunction(Control.LevelSetTwoInitialValue);
                    LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                    break;

                case GetLevelSet.FromOldSimulation: //Assuming that we already Loaded the ShockLevelSet
                case GetLevelSet.FromReconstruction:
                    LevelSetOpti.AssembleTransMat(LsTBO);
                    LsTBO.Clear();
                    this.LsTBO.Acc(1.0, this.ShockLevelSetField);
                    LevelSetOpti.ProjectFromLevelSet(LsTBO);
                    LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                    break;

                case GetLevelSet.FromReconstructionFromPoints: //here we need to have a density field and points to work with

                    if (Control.PointPath == null)
                    {
                        throw new NotSupportedException("must specify a PointPath");
                    }
                    if (Control.GetInitialValue == GetInitialValue.FromFunctionPerSpecies)
                    {
                        throw new NotSupportedException("GetLevelSet.FromReconstructionFromPoints if Initial value is loaded from DB");
                    }

                    switch (Control.OptiLevelSetType)
                    {
                        // in case of a spline we can directly get the spline from the points
                        case OptiLevelSetType.SplineLevelSet:
                            //get Points
                            if (LevelSetOpti is SplineOptiLevelSet splineLS)
                            {
                                var points2 = IMatrixExtensions.LoadFromTextFile(Control.PointPath);
                                var points3 = points2.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { points2.GetLength(0) - 1, 1 });
                                splineLS.GetSplineOverDetermined(points3);
                                splineLS.Interpolate();
                                splineLS.ProjectOntoLevelSet(LsTBO);
                            }
                            break;
                        default:
                            LevelSetOpti.AssembleTransMat(LsTBO);

                            //load the density of last time step
                            var densityField = tsiFromDb.Fields.Where(f => f.Identification == "rho").SingleOrDefault() as SinglePhaseField;
                            var points = IMatrixExtensions.LoadFromTextFile(Control.PointPath);
                            var tmpLS = ShockFindingExtensions.ReconstructLevelSetField(densityField, points);
                            LevelSetOpti.ProjectFromForeignLevelSet(tmpLS);
                            LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                            break;
                    }
                    break;
                case GetLevelSet.DirectyFromTimestep:
                    if (LevelSetOpti is SplineOptiLevelSet splineLs)
                    {
                        if (((TimestepProxy)tsiFromDb).GetInternal() is IDTTimeStepInfo IDTtsInfo)
                        {
                            var LSParams = IDTtsInfo.LevelSetParams;
                            if (LSParams.Length == splineLs.m_AllParams.Length)
                            {
                                for (int i = 0; i < LSParams.Length; i++)
                                {
                                    splineLs.m_AllParams[i] = LSParams[i];
                                }
                                splineLs.Interpolate();
                                splineLs.ProjectOntoLevelSet(LsTBO);
                                this.gamma = IDTtsInfo.Gamma;
                            }
                            else
                            {
                                throw new NotSupportedException($"{GetLevelSet.DirectyFromTimestep} not supported if Grid has not equal y - Cells");
                            }

                        }

                        //splineLs.x = points3.To1DArray();

                    }
                    else
                    {
                        throw new NotSupportedException($"{Control.GetLevelSet} only supported for SplineOptiLevelSet");
                    }
                    break;
            }

            LsTrk.UpdateTracker(CurrentStepNo);

            #endregion


            #region Initial Values
            // Clear fields (just to be sure)
            this.Density.Clear();
            this.Momentum[0].Clear();
            this.Momentum[1].Clear();
            this.Energy.Clear();

            // Seed from level set reconstruction run or apply initial conditions from the control file
            switch (Control.GetInitialValue)
            {
                #region FromFunctionPerSpecies
                case GetInitialValue.FromFunctionPerSpecies:
                case GetInitialValue.FromP0Timestepping:
                    // Check if initial conditions have been specified in primitive or conservative variables
                    bool primitiveVariables = true;
                    foreach (string key in Control.InitialValues_Evaluators.Keys)
                    {
                        primitiveVariables = key.Contains(CompressibleVariables.Density)
                          || key.Contains(XESFVariables.Velocity.xComponent)
                          || key.Contains(XESFVariables.Velocity.yComponent)
                          || key.Contains(XESFVariables.Pressure);

                        if (!primitiveVariables)
                        {
                            break;
                        }
                    }
                    if (primitiveVariables)
                        // Initial conditions are specified in primitive variables
                        foreach (SpeciesId id in this.SpeciesToEvaluate_Ids)
                        {
                            SeedSpecies_PrimVars(id, scheme);
                        }
                    else
                    {
                        // Initial conditions are specified in conservative variables
                        foreach (SpeciesId id in this.SpeciesToEvaluate_Ids)
                        {
                            // Helpers
                            Func<double[], double> density = Control.InitialValues_Evaluators[CompressibleVariables.Density + "#" + LsTrk.GetSpeciesName(id)];
                            Func<double[], double> momentumX = Control.InitialValues_Evaluators[CompressibleVariables.Momentum.xComponent + "#" + LsTrk.GetSpeciesName(id)];
                            Func<double[], double> momentumY = Control.InitialValues_Evaluators[CompressibleVariables.Momentum.yComponent + "#" + LsTrk.GetSpeciesName(id)];
                            Func<double[], double> energy = Control.InitialValues_Evaluators[CompressibleVariables.Energy + "#" + LsTrk.GetSpeciesName(id)];

                            // Density
                            this.Density.GetSpeciesShadowField(id).ProjectField(1.0, density);

                            // Momentum
                            this.Momentum[0].GetSpeciesShadowField(id).ProjectField(1.0, momentumX);
                            this.Momentum[1].GetSpeciesShadowField(id).ProjectField(1.0, momentumY);

                            // Energy
                            this.Energy.GetSpeciesShadowField(id).ProjectField(1.0, energy);
                        }
                    }
                    break;
                #endregion
                #region FromDB SinglePhase
                case GetInitialValue.FromDBSinglePhase:

                    ConventionalDGField rho = tsiFromDb.Fields.Where(f => f.Identification == "rho").SingleOrDefault() as ConventionalDGField;
                    ConventionalDGField m0 = tsiFromDb.Fields.Where(f => f.Identification == "m0").SingleOrDefault() as ConventionalDGField;
                    ConventionalDGField m1 = tsiFromDb.Fields.Where(f => f.Identification == "m1").SingleOrDefault() as ConventionalDGField;
                    ConventionalDGField rhoE = tsiFromDb.Fields.Where(f => f.Identification == "rhoE").SingleOrDefault() as ConventionalDGField;

                    //var rho = tsiFromDb.Fields.Where(f => f.Identification == "rho").SingleOrDefault();
                    //var m0 = tsiFromDb.Fields.Where(f => f.Identification == "m0").SingleOrDefault();
                    //var m1 = tsiFromDb.Fields.Where(f => f.Identification == "m1").SingleOrDefault();
                    //var rhoE = tsiFromDb.Fields.Where(f => f.Identification == "rhoE").SingleOrDefault();

                    foreach (string s in Control.SpeciesToEvaluate)
                    {
                        Density.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rho);
                        Momentum[0].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m0);
                        Momentum[1].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m1);
                        Energy.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rhoE);
                    }
                    //ReinitializeAllCutCells();

                    //ReinitializeNearBand();
                    break;
                #endregion
                #region FromDB XDG
                case GetInitialValue.FromDBXDG:

                    XDGField rhoXDG = (XDGField)tsiFromDb.Fields.Where(f => f.Identification == "rho").SingleOrDefault();
                    XDGField m0XDG = (XDGField)tsiFromDb.Fields.Where(f => f.Identification == "m0").SingleOrDefault();
                    XDGField m1XDG = (XDGField)tsiFromDb.Fields.Where(f => f.Identification == "m1").SingleOrDefault();
                    XDGField rhoEXDG = (XDGField)tsiFromDb.Fields.Where(f => f.Identification == "rhoE").SingleOrDefault();

                    foreach (string s in Control.SpeciesToEvaluate)
                    {
                        Density.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rhoXDG.GetSpeciesShadowField(s));
                        Momentum[0].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m0XDG.GetSpeciesShadowField(s));
                        Momentum[1].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m1XDG.GetSpeciesShadowField(s));
                        Energy.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rhoEXDG.GetSpeciesShadowField(s));
                    }
                    break;
                #endregion
                #region FromDB Partial SeedInflowExactly
                case GetInitialValue.FromDB_Partial_SeedInflowExactly:
                    if (Control.ShockLevelSet_SeedFromDb)
                    {
                        ConventionalDGField rho_inflow = tsiFromDb.Fields.Where(f => f.Identification == "rho").SingleOrDefault() as ConventionalDGField;
                        ConventionalDGField m0_inflow = tsiFromDb.Fields.Where(f => f.Identification == "m0").SingleOrDefault() as ConventionalDGField;
                        ConventionalDGField m1_inflow = tsiFromDb.Fields.Where(f => f.Identification == "m1").SingleOrDefault() as ConventionalDGField;
                        ConventionalDGField rhoE_inflow = tsiFromDb.Fields.Where(f => f.Identification == "rhoE").SingleOrDefault() as ConventionalDGField;

                        // On the left of the shock: Seed exact inflow values
                        Console.WriteLine("Using exact inflow conditions as pre-shock state (\"L\")...");
                        SeedSpecies_PrimVars(LsTrk.GetSpeciesId("L"), scheme);

                        // On the right of the shock: Seed values from previous simulation run
                        Density.GetSpeciesShadowField("R").ProjectFromForeignGrid(1.0, rho_inflow);
                        Momentum[0].GetSpeciesShadowField("R").ProjectFromForeignGrid(1.0, m0_inflow);
                        Momentum[1].GetSpeciesShadowField("R").ProjectFromForeignGrid(1.0, m1_inflow);
                        Energy.GetSpeciesShadowField("R").ProjectFromForeignGrid(1.0, rhoE_inflow);
                    }
                    break;
                #endregion
                #region FromDB AVRun
                case GetInitialValue.FromAVRun:

                    XDGField rho_2 = tsiFromDb.Fields.Where(f => f.Identification == "rho").SingleOrDefault() as XDGField;
                    XDGField m0_2 = tsiFromDb.Fields.Where(f => f.Identification == "m0").SingleOrDefault() as XDGField;
                    XDGField m1_2 = tsiFromDb.Fields.Where(f => f.Identification == "m1").SingleOrDefault() as XDGField;
                    XDGField rhoE_2 = tsiFromDb.Fields.Where(f => f.Identification == "rhoE").SingleOrDefault() as XDGField;

                    foreach (string s in Control.SpeciesToEvaluate)
                    {
                        Density.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rho_2.GetSpeciesShadowField("A"));
                        Momentum[0].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m0_2.GetSpeciesShadowField("A"));
                        Momentum[1].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m1_2.GetSpeciesShadowField("A"));
                        Energy.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rhoE_2.GetSpeciesShadowField("A"));
                    }
                    break;
                #endregion
                default:
                    throw new NotImplementedException();
            }
            #endregion

        }
        /// <summary>
        /// sets initial value prescribed in primitive variables
        /// </summary>
        /// <param name="id"></param>
        /// <param name="scheme"></param>
        private void SeedSpecies_PrimVars(SpeciesId id, CellQuadratureScheme scheme)
        {
            // Helpers
            Func<double[], double> density = Control.InitialValues_Evaluators[CompressibleVariables.Density + "#" + LsTrk.GetSpeciesName(id)];
            Func<double[], double> velocityX = Control.InitialValues_Evaluators[XESFVariables.Velocity.xComponent + "#" + LsTrk.GetSpeciesName(id)];
            Func<double[], double> velocityY = Control.InitialValues_Evaluators[XESFVariables.Velocity.yComponent + "#" + LsTrk.GetSpeciesName(id)];
            Func<double[], double> pressure = Control.InitialValues_Evaluators[XESFVariables.Pressure + "#" + LsTrk.GetSpeciesName(id)];

            // Density
            this.Density.GetSpeciesShadowField(id).ProjectField(1.0, density);

            // Momentum (rho * momentum vector)
            this.Momentum[0].GetSpeciesShadowField(id).ProjectField(1.0, X => density(X) * velocityX(X));
            this.Momentum[1].GetSpeciesShadowField(id).ProjectField(1.0, X => density(X) * velocityY(X));

            // Energy
            Energy.GetSpeciesShadowField(id).ProjectField(
                1.0,
                delegate (double[] X)
                {
                    Vector velocityVec = new Vector(LsTrk.GridDat.SpatialDimension);
                    velocityVec[0] = velocityX(X);
                    velocityVec[1] = velocityY(X);

                    StateVector state = StateVector.FromPrimitiveQuantities(this.Control.GetMaterial(), density(X), velocityVec, pressure(X));
                    return state.Energy;
                },
                    scheme);
        }
        ///
        #region XDG auxiliary functions
        private IDictionary<SpeciesId, IEnumerable<double>> MassScale
        {
            get
            {
                Dictionary<SpeciesId, IEnumerable<double>> massScaleToSpeciesMap = new Dictionary<SpeciesId, IEnumerable<double>>();

                int D = this.Grid.SpatialDimension;
                double[] ones = new double[D + 2];
                ones.SetAll(1.0);

                foreach (SpeciesId speciesId in this.SpeciesToEvaluate_Ids)
                {
                    massScaleToSpeciesMap.Add(speciesId, ones);
                }
                return massScaleToSpeciesMap;
            }
        }
        #endregion
        /// <summary>
        /// initializes Conservative Fields for the Euler Equations
        /// </summary>
        /// <param name="in_LsTrk"></param>
        /// <param name="in_DgDegree"></param>
        /// <exception cref="NotSupportedException"></exception>
        /// <exception cref="NotImplementedException"></exception>
        public override void CreateConservativeFields(LevelSetTracker in_LsTrk, int in_DgDegree)
        {
            #region Optional: load shock level set DG field from a CNS simulation stored in a database
            switch (Control.GetLevelSet)
            {
                case GetLevelSet.FromOldSimulation:
                case GetLevelSet.FromReconstruction:

                    DatabaseInfo dbi = DatabaseInfo.Open(Control.ShockLevelSet_Db);
                    //DatabaseInfo dbi = DatabaseInfo.Open(@"C:\BoSSS-experimental\internal\src\private-mag\XESF\Tests\bosss_db_levelSets.zip");
                    ISessionInfo si = dbi.Controller.GetSessionInfo(Control.ShockLevelSet_Info.Item1);
                    if (Control.ShockLevelSet_Info.Item2.MajorNumber < 0.0)
                    {
                        tsiFromDb = si.Timesteps.Last();
                    }
                    else
                    {
                        throw new NotSupportedException("Currently, only the last time step can be used.");
                    }

                    // Check if grid from shock level set reconstruction equals current grid
                    if (tsiFromDb.Grid.ID != this.GridData.GridID)
                    {
                        Console.WriteLine("Different grids used for LEVEL SET RECONSTRUCTION and XESF");
                        IGrid grid = dbi.Controller.DBDriver.LoadGrid(tsiFromDb.GridID, dbi);

                        // Load grid and shock level set field
                        var fieldsFromDb = dbi.Controller.DBDriver.LoadFields(tsiFromDb, grid.iGridData);

                        // Project onto current grid
                        if (Control.IsTwoLevelSetRun)
                        {
                            ShockLevelSetField = new SinglePhaseField(new Basis(this.GridData, Control.LevelSetTwoDegree), "shockLevelSetField");
                        }
                        else
                        {
                            ShockLevelSetField = new SinglePhaseField(new Basis(this.GridData, Control.LevelSetDegree), "shockLevelSetField");
                        }

                        ShockLevelSetField.ProjectFromForeignGrid(1.0, (ConventionalDGField)fieldsFromDb.Pick(1));

                    }
                    else
                    {
                        ShockLevelSetField = (SinglePhaseField)dbi.Controller.DBDriver.LoadFields(tsiFromDb, this.GridData, new string[] { Control.ShockLevelSet_FieldName }).Single();
                    }
                    break;
                case GetLevelSet.FromReconstructionFromPoints:
                case GetLevelSet.DirectyFromTimestep:
                case GetLevelSet.FromParams:
                case GetLevelSet.FromFunction:
                    break;
                default:
                    throw new NotImplementedException();
            }
            #endregion
            #region Create mandatory conservative fields
            int D = in_LsTrk.GridDat.SpatialDimension;

            this.Density = new XDGField(new XDGBasis(in_LsTrk, in_DgDegree), CompressibleVariables.Density);
            this.Density.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            XDGField[] momentumFields = new XDGField[D];
            XDGBasis momentumBasis = new XDGBasis(in_LsTrk, in_DgDegree);
            for (int d = 0; d < D; d++)
            {
                string variableName = CompressibleVariables.Momentum[d];
                momentumFields[d] = new XDGField(momentumBasis, variableName);
                momentumFields[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            }
            this.Momentum = new VectorField<XDGField>(momentumFields);
            this.Energy = new XDGField(new XDGBasis(in_LsTrk, in_DgDegree), CompressibleVariables.Energy);
            this.Energy.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;

            ConservativeFields = new XDGField[D + 2];



            ConservativeFields[0] = Density;
            for (int d = 0; d < D; d++)
            {
                ConservativeFields[d + 1] = Momentum[d];
            }
            ConservativeFields[D + 1] = Energy;

            #endregion

            #region Create derived fields
            this.DerivedVariableToXDGFieldMap = new Dictionary<DerivedVariable<XDGField>, XDGField>();
            this.DerivedVariableToSinglePhaseFieldMap = new Dictionary<DerivedVariable<SinglePhaseField>, SinglePhaseField>();
            this.DerivedVariableToDoubleMap = new Dictionary<DerivedVariable<double>, double>();
            foreach (KeyValuePair<Variable, int> pair in this.Control.VariableToDegreeMap)
            {
                Variable variable = pair.Key;
                int degree = pair.Value;

                if (variable is DerivedVariable<XDGField> xdgVar)
                {
                    this.DerivedVariableToXDGFieldMap.Add(xdgVar, new XDGField(new XDGBasis(in_LsTrk, degree), variable.Name));
                }
                else if (variable is DerivedVariable<SinglePhaseField> singleVar)
                {
                    this.DerivedVariableToSinglePhaseFieldMap.Add(singleVar, new SinglePhaseField(new Basis(in_LsTrk.GridDat, degree), variable.Name));
                }
                else if (variable is DerivedVariable<double> doubleVar)
                {
                    this.DerivedVariableToDoubleMap.Add(doubleVar, new double());
                }
            }
            #endregion
            #region register fields
            #endregion

            #region register the derived Fields
            this.m_IOFields.AddRange(DerivedVariableToXDGFieldMap.Values);
            this.m_IOFields.AddRange(DerivedVariableToSinglePhaseFieldMap.Values);

            this.m_RegisteredFields.AddRange(DerivedVariableToXDGFieldMap.Values);
            this.m_RegisteredFields.AddRange(DerivedVariableToSinglePhaseFieldMap.Values);
            #endregion 

        }
        /// <summary>
        /// updates derived variables such as pressure, velocity, enthalpy, ...
        /// used in every time step
        /// </summary>
        public override void UpdateDerivedVariables()
        {
            // Update derived variables (including sensor _variable_)
            foreach (KeyValuePair<DerivedVariable<XDGField>, XDGField> pair in this.DerivedVariableToXDGFieldMap)
            {
                pair.Key.UpdateFunction(pair.Value, this);
            }

            foreach (KeyValuePair<DerivedVariable<SinglePhaseField>, SinglePhaseField> pair in this.DerivedVariableToSinglePhaseFieldMap)
            {
                pair.Key.UpdateFunction(pair.Value, this);
            }

            foreach (KeyValuePair<DerivedVariable<double>, double> pair in this.DerivedVariableToDoubleMap)
            {
                using (var tr = new FuncTrace())
                {
                    pair.Key.UpdateFunction(pair.Value, this);
                    tr.Info(pair.Key.ToString() + "=" + pair.Value.ToString());
                }
            }
        }
        /// <summary>
        /// Creates the spatial operator from a Control and a LevelSetTracker
        /// Helper function for debugging
        /// </summary>
        /// <param name="Control"></param>
        /// <param name="LsTrk"></param>
        /// <returns></returns>
        public static XDifferentialOperatorMk2 BuildOperatorFrom_Control_LsTrk(XESFControl Control, LevelSetTracker LsTrk)
        {
            var grid = Control.GridFunc();
            return BuildOperatorFrom_Control_LsTrk_Grid(Control, LsTrk, grid);
        }
        /// <summary>
        /// Creates the spatial operator from a Control, a LevelSetTracker and a grid
        /// Helper function for debugging
        /// DOes not have all the changes done in CreateEquationsAndSolvers, will throw errors for new fluxes.
        /// TODO: define a common function that is used in both methods remove duplicate code
        /// </summary>
        /// <param name="Control">Main control</param>
        /// <param name="LsTrk"> LevelSetTracker</param>
        /// <param name="grid">grid </param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"> gets thrown if for a FluxType there is no implementation</exception>
        public static XDifferentialOperatorMk2 BuildOperatorFrom_Control_LsTrk_Grid(XESFControl Control, LevelSetTracker LsTrk, IGrid grid)
        {
            GridData GridData = (GridData)grid.iGridData;
            Material material = Control.GetMaterial();
            IBoundaryConditionMap boundaryMap = new XDGCompressibleBoundaryCondMap(GridData, Control, material, Control.SpeciesToEvaluate);
            string[] variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };
            var XSpatialOperator = new XDifferentialOperatorMk2(variables, null, variables, Control.quadOrderFunc, Control.SpeciesToEvaluate);
            CompressibleEnvironment.Initialize(2);
            switch (Control.FluxVersion)
            {
                case FluxVersion.NonOptimized:
                    foreach (string spcNmn in Control.SpeciesToEvaluate)
                    {
                        var id = LsTrk.GetSpeciesId(spcNmn);
                        XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.HLLCDensityFlux_XDG(boundaryMap, material, spcNmn, id));
                        XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.HLLCMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id));
                        XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.HLLCMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id));
                        XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.HLLCEnergyFlux_XDG(boundaryMap, material, spcNmn, id));
                    }
                    switch (Control.ConvectiveInterfaceFlux_LsOne)
                    {
                        case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;
                        case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;

                        case ConvectiveInterfaceFluxes.None:
                            break;

                        default:
                            throw new NotSupportedException("This should never happen");
                    }
                    switch (Control.ConvectiveInterfaceFlux_LsTwo)
                    {
                        case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            break;

                        case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.HLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.HLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.HLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            break;
                        case ConvectiveInterfaceFluxes.None:
                            break;
                        default:
                            throw new NotSupportedException("This should never happen");
                    }
                    break;
                #region Markus FLuxes

                case FluxVersion.Optimized:
                    #region Convective bulk fluxes
                    switch (Control.ConvectiveBulkFlux)
                    {
                        case ConvectiveBulkFluxes.OptimizedHLLC:
                            foreach (string spcNmn in Control.SpeciesToEvaluate)
                            {
                                var id = LsTrk.GetSpeciesId(spcNmn);
                                XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG(boundaryMap, material, spcNmn, id));
                                XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id));
                                XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id));
                                XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG(boundaryMap, material, spcNmn, id));
                            }
                            break;

                        case ConvectiveBulkFluxes.Godunov:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new GodunovFlux(Control, boundaryMap, new EulerDensityComponent(), material));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new GodunovFlux(Control, boundaryMap, new EulerMomentumComponent(0, material.EquationOfState.HeatCapacityRatio, Control.MachNumber, GridData.SpatialDimension), material));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new GodunovFlux(Control, boundaryMap, new EulerMomentumComponent(1, material.EquationOfState.HeatCapacityRatio, Control.MachNumber, GridData.SpatialDimension), material));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new GodunovFlux(Control, boundaryMap, new EulerEnergyComponent(), material));
                            break;

                        case ConvectiveBulkFluxes.None:
                            break;

                        default:
                            throw new NotSupportedException("This should never happen");
                    }
                    #endregion
                    #endregion


                    #region Convective interface flux (level set one)
                    switch (Control.ConvectiveInterfaceFlux_LsOne)
                    {
                        case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;

                        case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            break;

                        //case ConvectiveInterfaceFluxes.OptimizedHLLCWall:
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Density));
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Momentum, 0));
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Momentum, 1));
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Energy));
                        //    break;

                        case ConvectiveInterfaceFluxes.GodunovInterface:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: Control.LsOne_NegSpecies, posSpecies: Control.LsOne_PosSpecies));

                            // Interface as wall
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Density));
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Momentum, 0));
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Momentum, 1));
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Energy));
                            break;

                        case ConvectiveInterfaceFluxes.None:
                            break;

                        default:
                            throw new NotSupportedException("This should never happen");
                    }
                    #endregion

                    #region Convective interface flux (level set two)
                    switch (Control.ConvectiveInterfaceFlux_LsTwo)
                    {
                        case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            break;

                        case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            break;

                        //case ConvectiveInterfaceFluxes.OptimizedHLLCWall:
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Density));
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Momentum, 0));
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Momentum, 1));
                        //     XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Energy));
                        //    break;

                        case ConvectiveInterfaceFluxes.GodunovInterface:
                            XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));
                            XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.GodunovFlux_Interface(LsTrk, Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: Control.LsTwo_NegSpecies, posSpecies: Control.LsTwo_PosSpecies));

                            // Interface as wall
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Density));
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Momentum, 0));
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Momentum, 1));
                            // XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new GodunovFlux_Wall(LsTrk,  Control, boundaryMap, material, FluxVariables.Energy));
                            break;

                        case ConvectiveInterfaceFluxes.None:
                            break;

                        default:
                            throw new NotSupportedException("This should never happen");
                    }
                    #endregion
                    break;
                default:
                    throw new NotSupportedException("This should never happen");
            }
            return XSpatialOperator;
        }



    }

    /// <summary>
    /// Runs test cases from all IDTSolvers shown in a publication (except bow shock, because takes to long)
    /// </summary>
    public static class IDTTestRunner
    {
        public static void RunTests()
        {
            BoSSS.Solution.Application.InitMPI(num_threads:1);
            XESFMain.DeleteOldPlotFiles();
            List<string> fails = new List<string>();

            try { SAIDT.Tests.SAIDTTestProgram.CurvedShock_Eccomas22(); } catch { fails.Add("SAIDTTestProgram.CurvedShock_Eccomas22"); }
            try { BUIDT.Tests.BUIDTTestProgram.StraightShockCurvedStart_Eccomas22(); } catch { fails.Add("BUIDTTestProgram.StraightShockCurvedStart_Eccomas22()"); }
            try { XESF.Tests.XESFTestProgram.XDG_SWF_OneLs_Cart(); } catch { fails.Add("XESFTestProgram.XDG_SWF_OneLs_Cart"); }
            try { BUIDT.Tests.BUIDTTestProgram.AcceleratingShock(); } catch { fails.Add("BUIDTTestProgram.AcceleratingShock"); }
            try { XESF.Tests.XESFTestProgram.XDG_SWF_TwoLs(); } catch { fails.Add("XESFTestProgram.XDG_SWF_TwoLs"); }
            if (fails.Count > 0)
            {
                throw new Exception(fails.ToConcatString("Fails: (", ", ", ")"));
            }
        }
    }
}



