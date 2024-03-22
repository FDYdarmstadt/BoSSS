using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Residual;
using BoSSS.Solution.Statistic;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using ilPSP.Tracing;
using XESF.Fluxes;
using ApplicationWithIDT;
using BoSSS.Solution.CompressibleFlowCommon.ShockFinding;
using ApplicationWithIDT.OptiLevelSets;
using ilPSP.LinSolvers;
using XESTSF.Variables;
using XESTSF.SPFluxes;
using XESTSF.Fluxes;
using static BoSSS.Solution.GridImport.NASTRAN.NastranFile;

namespace XESTSF {
    public class XESTSFMain : ApplicationWithIDT<XESTSFControl> {
        /// <summary>
        /// Implements shock fitting for XDG Space Time Euler which is solved by the routines defined in ApplicationWithIDT 
        /// Naming: X(DG) - E(uler) - S(pace) - T(ime) - S(hock) - F(itting)
        /// Concrete configurations of solver (Initial Guess, optimization parameters,...) are set in a XESFControl.cs object
        /// Fluxes are implemented in XESTSF.Fluxes. This sovler supports:
        /// - Space Time Roe Flux (with smoothing factor alpha)
        /// 
        /// Author: Jakob Vandergrift 
        /// Date of Creation/Maintanance: 07-2023 until at least 08-2024
        /// </summary>
        /// <param name="args">string pointing to a control file, i.e. 'cs:XESF.XESFHardCodedControl.XDGWedgeFlow_TwoLs_Base()' </param>
        static void Main(string[] args) {

            //Test00();
            //TestSpaceTimeSlabs();
            //BoSSS.Solution.Application.InitMPI();
            //DeleteOldPlotFiles();

            XESTSFMain._Main(args, false, () => new XESTSFMain());
        }
        public override int GetGlobalQuadOrder() {
            //return 15;
            return base.GetGlobalQuadOrder();
        }
        static void Test00() {
            BoSSS.Solution.Application.InitMPI();
            DeleteOldPlotFiles();
            var p = new XESTSFMain();
            //var c = XESTSFHardCodedControl.XDG_1D_Acoustic_Wave_Interaction_TwoLs_Base();
            //var c = XESTSFHardCodedControl.MovingShockWave2(shkspd:0.4);
            //var c = XESTSFHardCodedControl.AcousticWave1D(withLevelSet:false,acoustic_amplitude_neg:0.001);
            //var c = XESTSFHardCodedControl.AcousticWave1D(withShock: true, withLevelSet: true, p_amp_neg: 1e-5, scaling:100);
            //var c = XESTSFHardCodedControl.StationaryShockWave();
            var c = XESTSFHardCodedControl.ShuOsher1D(numOfCellsX: 70, numOfCellsT: 10,PlotInterval:1);
            p.Init(c);
            p.RunSolverMode();
            p.ResNorms.SaveToTextFile("ResNorms.txt");
            p.GetResPlot().SaveTabular("ResNormsTab.txt",false);
            p.GetResPlot().SavePgfplotsFile("ResNormsTab.txt");
        }
        static void TestSpaceTimeSlabs() {
            BoSSS.Solution.Application.InitMPI();
            DeleteOldPlotFiles();
            var p = new XESTSFMain();
            //var c = XESTSFHardCodedControl.XDG_1D_Acoustic_Wave_Interaction_TwoLs_Base();
            //var c = XESTSFHardCodedControl.MovingShockWave2(shkspd:0.4);
            string waveform= "bump";
            var c = XESTSFHardCodedControl.AcousticWave1D(waveform: waveform, withLevelSet: false, p_amp_neg: 0.001,MaxIterations:1,numOfCellsX:20, numOfCellsT:20);
            //var c = XESTSFHardCodedControl.StationaryShockWave()
            p.Init(c);
            p.RunSolverMode();

            

            //// artificially reset u_previous to Initial Value
            DeleteOldPlotFiles();
            ////plot converged state
            p.PlotCurrentState(0);
            for(int i = 0; i < p.ConservativeFields.Length; i++) {
                p.ConservativeFields[i].Clear();
            }
            // reset to initial value
            p.SetInitial(0.0);
            DGField[] u_previous = new DGField[p.ConservativeFields.Length];
            for(int i = 0; i < p.ConservativeFields.Length; i++) {
                u_previous[i] = p.ConservativeFields[i].CloneAs();
            }
            p.PlotCurrentState(1);

            var c2 = XESTSFHardCodedControl.AcousticWave1D(waveform: waveform, withLevelSet: false, p_amp_neg: 0.001, MaxIterations: 20, numOfCellsX: 20, numOfCellsT: 20);
            XESTSFControl.UpdateControl(c2, 1.0, u_previous, _tNref: 10);
            p = new XESTSFMain();
            p.Init(c2);
            p.RunSolverMode();

            p.ResNorms.SaveToTextFile("ResNorms.txt");
            p.GetResPlot().SaveTabular("ResNormsTab.txt", false);
            p.GetResPlot().SavePgfplotsFile("ResNormsTab.txt");
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            if(Control.DoSpaceTimeSlabs) {
                if(TimestepNo == 1) { //first timestep has no previous U
                    while(base.TerminationKey == false && base.CurrentStepNo < Control.MaxIterations) {
                        base.RunSolverOneStep(base.CurrentStepNo, phystime, dt);
                    }
                } else { // do space time slabs

                    //Copy the actual ConsVariables
                    DGField[] u_previous = new DGField[ConservativeFields.Length];
                    for(int i = 0; i < ConservativeFields.Length; i++) {
                        u_previous[i] = ConservativeFields[i].CloneAs();
                    }

                    //Do the timestep running an XESTSF solver
                    XESTSFControl.UpdateControl(Control, dt, u_previous, _tNref: 1, Control.MaxIterations);
                    var p = new XESTSFMain();
                    p.Init(Control);
                    p.RunSolverMode();

                    // Merge Grids
                    var grid_prev = (GridCommons) ConservativeFields[0].GridDat.Grid;
                    var grid_add = (GridCommons) p.ConservativeFields[0].GridDat.Grid;
                    var grd_updated = GridCommons.MergeLogically(grid_prev, grid_add);
                    LevelSetTracker lstrk_updated;
                    if(Control.IsTwoLevelSetRun) {
                        lstrk_updated = new LevelSetTracker(grd_updated.GridData, LsTrk.CutCellQuadratureType, LsTrk.NearRegionWidth, Control.SpeciesTable, this.LevelSet, this.LevelSetTwo);
                    } else {
                        lstrk_updated = new LevelSetTracker(grd_updated.GridData, LsTrk.CutCellQuadratureType, LsTrk.NearRegionWidth, Control.SpeciesTable, this.LevelSet);
                    }

                    CreateConservativeFields(lstrk_updated, Control.SolDegree);
                    //Initialize all dependent Fields
                    InitDependentFields();
                    InitObjFields();
                    //Register them (also registers other fields)
                    RegisterFields();


                }
                return 0;
            } else {
                return base.RunSolverOneStep(TimestepNo, phystime, dt);
            } 
        }

        #region XDG Fields
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
        /// All conservative fields (density, momentum vector, energy)
        /// </summary>

        #endregion

        CoordinateVector m_UnknownsVector;
        CoordinateVector UnknownsVector {
            get {
                if(m_UnknownsVector == null)
                    m_UnknownsVector = new CoordinateVector(ConservativeFields);
                return m_UnknownsVector;
            }
        }


        public SinglePhaseField ShockLevelSetField;
        public ITimestepInfo tsiFromDb;

        public IResidualLogger[] residualLoggers;

        public Dictionary<DerivedVariable<XDGField>, XDGField> DerivedVariableToXDGFieldMap {
            get;
            private set;
        }

        public Dictionary<DerivedVariable<SinglePhaseField>, SinglePhaseField> DerivedVariableToSinglePhaseFieldMap {
            get;
            private set;
        }

        public Dictionary<DerivedVariable<double>, double> DerivedVariableToDoubleMap {
            get;
            private set;
        }

        protected override IGrid CreateOrLoadGrid() {
            IGrid grid = null;

            // Optional: load grid from CNS calculation from which the shock level set has been reconstructed
            if(Control.ShockLevelSet_UseShockGrid) {
                IDatabaseInfo dbi = DatabaseInfo.Open(Control.ShockLevelSet_Db);
                ISessionInfo si = dbi.Controller.GetSessionInfo(Control.ShockLevelSet_Info.Item1);
                ITimestepInfo tsi;

                if(Control.ShockLevelSet_Info.Item2.MajorNumber < 0.0) {
                    tsi = si.Timesteps.Last();
                } else {
                    throw new NotSupportedException("Currently, only the last time step can be used.");
                }

                grid = dbi.Controller.DBDriver.LoadGrid(tsi.GridID, dbi);

                Console.WriteLine("Using grid from LEVEL SET RECONSTRUCTION");
            } else {
                grid = base.CreateOrLoadGrid();
            }

            CompressibleEnvironment.Initialize(grid.SpatialDimension);

            return grid;
        }
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            #region Init operator
            //// Cell agglomerator (cell length scales are needed for diffusive AV fluxes
            var D = Grid.SpatialDimension;

            //if(this.Control.SaveAgglomerationPairs) {
            //    MultiphaseAgglomerator.PlotAgglomerationPairs("agglomerationPairs");
            //}

            GridData gridData = (GridData)this.GridData;
            Material material = this.Control.GetMaterial();
            IBoundaryConditionMap boundaryMap = new XDGCompressibleBoundaryCondMap(this.GridData, this.Control, material, this.Control.SpeciesToEvaluate);
            string[] variables;
            if(Grid.SpatialDimension == 2) {
                variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Energy };
            } else if(Grid.SpatialDimension == 3) {
                variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };
            } else {
                variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent,CompressibleVariables.Momentum.zComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };
            }


            this.XSpatialOperator = new XDifferentialOperatorMk2(variables, null, variables, Control.quadOrderFunc, this.SpeciesToEvaluate);
            this.XSpatialOperator.FluxesAreNOTMultithreadSafe = true;
            #endregion
            DGField[] prevU;
            #region Convective bulk fluxes
            switch(Control.ConvectiveBulkFlux) {
                case ConvectiveBulkFluxes.OptimizedHLLC:
                throw new NotImplementedException();
                case ConvectiveBulkFluxes.Roe:

                //foreach(SpeciesId id in SpeciesToEvaluate_Ids) {
                //    string spcNmn = LsTrk.GetSpeciesName(id);
                //    var Sdensityflux = new XESF.Fluxes.RoeDensityFlux_XDG(boundaryMap, material, spcNmn, id,10, D: D - 1);
                //    var SmomXflux = new XESF.Fluxes.RoeMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id,10, m_D: D - 1);
                //    var SmomYflux = new XESF.Fluxes.RoeMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id,10, m_D: D - 1);
                //    var Senergyflux = new XESF.Fluxes.RoeEnergyFlux_XDG(boundaryMap, material, spcNmn, id,10, D: D - 1);

                //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPBulkFlux(boundaryMap,(IEdgeForm)Sdensityflux, (IVolumeForm)Sdensityflux, 0, id, spcNmn, material));
                //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPBulkFlux(boundaryMap,(IEdgeForm)SmomXflux, (IVolumeForm)SmomXflux, 1, id, spcNmn, material));
                //    int next_comp = 2;
                //    if(Grid.SpatialDimension == 3) {
                //        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPBulkFlux(boundaryMap,(IEdgeForm)SmomYflux, (IVolumeForm)SmomYflux, next_comp, id, spcNmn, material));
                //        next_comp++;
                //    }
                //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPBulkFlux(boundaryMap,(IEdgeForm)Senergyflux, (IVolumeForm)Senergyflux, next_comp, id, spcNmn, material));
                //}
                foreach(SpeciesId id in SpeciesToEvaluate_Ids) {
                    string spcNmn = LsTrk.GetSpeciesName(id);
                    if(Control.previous_u is XDGField[] uXDG) {
                        prevU = uXDG.Select(u => (DGField)u.GetSpeciesShadowField(spcNmn)).ToArray();
                    } else {
                        prevU = Control.previous_u;
                    }
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new RoeSTDensityFlux_XDG(boundaryMap, material, spcNmn, id, Control.flux_s_alpha, _previous_u: prevU, m_D:D-1));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new RoeSTMomentumFlux_XDG(boundaryMap, 0, material, spcNmn, id, Control.flux_s_alpha, _previous_u: prevU, D: D - 1));
                    if(Grid.SpatialDimension >= 3) {
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new RoeSTMomentumFlux_XDG(boundaryMap, 1, material, spcNmn, id, Control.flux_s_alpha, prevU, D: D - 1));
                    }
                    if(Grid.SpatialDimension >= 4) {
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.zComponent].Add(new RoeSTMomentumFlux_XDG(boundaryMap, 2, material, spcNmn, id, Control.flux_s_alpha, prevU, D: D - 1));
                    }
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new RoeSTEnergyFlux_XDG(boundaryMap, material, spcNmn, id, Control.flux_s_alpha, prevU, m_D: D - 1));
                }
                break;

                case ConvectiveBulkFluxes.Godunov:
                throw new NotImplementedException();
                //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new GodunovFlux(this.Control, boundaryMap, new EulerDensityComponent(), material));
                //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new GodunovFlux(this.Control, boundaryMap, new EulerMomentumComponent(0, material.EquationOfState.HeatCapacityRatio, Control.MachNumber, gridData.SpatialDimension), material));
                //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new GodunovFlux(this.Control, boundaryMap, new EulerMomentumComponent(1, material.EquationOfState.HeatCapacityRatio, Control.MachNumber, gridData.SpatialDimension), material));
                //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new GodunovFlux(this.Control, boundaryMap, new EulerEnergyComponent(), material));
                break;

                case ConvectiveBulkFluxes.None:
                break;

                default:
                throw new NotSupportedException("This should never happen");
            }
            #endregion
            #region Convective interface flux (level set one)
            if(Control.LsOne_SpeciesPairs == null) {
                throw new ArgumentException("need to provide species pairs for level set one");
            }
            for(int i = 0; i < Control.LsOne_SpeciesPairs.GetLength(0); i++) {
                string negSpec = Control.LsOne_SpeciesPairs[i, 0];
                string posSpec = Control.LsOne_SpeciesPairs[i, 1];
                switch(Control.ConvectiveInterfaceFlux_LsOne) {
                    case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                    throw new NotImplementedException();
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec));
                    break;

                    case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                    throw new NotImplementedException();
                    //var Sdensityflux = new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var SmomXflux = new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var SmomYflux = new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var Senergyflux = new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPInterfaceFlux((ILevelSetForm)Sdensityflux, 0));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomXflux, 1));
                    //int next_comp = 2;
                    //if(Grid.SpatialDimension == 3) {
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux, next_comp));
                    //    next_comp++;
                    //}
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPInterfaceFlux((ILevelSetForm)Senergyflux, next_comp));
                    break;

                    case ConvectiveInterfaceFluxes.GodunovInterface:
                    throw new NotImplementedException();
                    //var Sdensityflux2 = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var SmomXflux2 = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var SmomYflux2 = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var Senergyflux2 = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPInterfaceFlux((ILevelSetForm)Sdensityflux2, 0));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomXflux2, 1));
                    //int next_comp2 = 2;
                    //if(Grid.SpatialDimension == 3) {
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux2, next_comp2));
                    //    next_comp2++;
                    //}
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPInterfaceFlux((ILevelSetForm)Senergyflux2, next_comp2));
                    break;
                    case ConvectiveInterfaceFluxes.CentralFluxInterface:
                    throw new NotImplementedException();
                    //var Sdensityflux3 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var SmomXflux3 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var SmomYflux3 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //var Senergyflux3 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec);
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPInterfaceFlux((ILevelSetForm)Sdensityflux3, 0));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomXflux3, 1));
                    //int next_comp3 = 2;
                    //if(Grid.SpatialDimension == 3) {
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux3, next_comp3));
                    //    next_comp3++;
                    //}
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPInterfaceFlux((ILevelSetForm)Senergyflux3, next_comp3));
                    break;
                    case ConvectiveInterfaceFluxes.RoeInterface:
                    if(Control.previous_u is XDGField[] uXDG) {
                        prevU = uXDG.Select(u => (DGField)u.GetSpeciesShadowField(posSpec)).ToArray();
                    } else {
                        prevU = Control.previous_u;
                    }
                    //var Sdensityflux4 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, component: 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D - 1);
                    //var SmomXflux4 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, component: 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D - 1);
                    //var SmomYflux4 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, component: 1, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D - 1);
                    //var SmomZflux4 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, component: 2, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D - 1);
                    //var Senergyflux4 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, component: 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D-1);
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPInterfaceFlux((ILevelSetForm)Sdensityflux4, 0));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomXflux4, 1));
                    //int next_comp4 = 2;
                    //if(Grid.SpatialDimension == 3) {
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux4, next_comp4));
                    //    next_comp4++;
                    //}
                    //if(Grid.SpatialDimension == 4) {
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux4, next_comp4));
                    //    next_comp4++;
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomZflux4, next_comp4));
                    //    next_comp4++;
                    //}
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPInterfaceFlux((ILevelSetForm)Senergyflux4, next_comp4));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Density, int.MinValue, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1, hasDirichletBoundary: Control.hasDirichletBoundary, DirichletBoundaryFunc: Control.DirichletBoundaryFunc));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Momentum, 0, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1, hasDirichletBoundary: Control.hasDirichletBoundary, DirichletBoundaryFunc: Control.DirichletBoundaryFunc));
                    if(Grid.SpatialDimension >= 3) {
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Momentum, 1, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1, hasDirichletBoundary:Control.hasDirichletBoundary, DirichletBoundaryFunc:Control.DirichletBoundaryFunc));
                    }
                    if(Grid.SpatialDimension >= 4) {
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.zComponent].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Momentum, 2, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1, hasDirichletBoundary: Control.hasDirichletBoundary, DirichletBoundaryFunc: Control.DirichletBoundaryFunc));
                    }
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Energy, int.MinValue, levelSetIndex: 0, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1, hasDirichletBoundary: Control.hasDirichletBoundary, DirichletBoundaryFunc: Control.DirichletBoundaryFunc));
                    //// Interface as wall
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Density));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Momentum, 0));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Momentum, 1));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Energy));
                    break;

                    case ConvectiveInterfaceFluxes.None:
                    break;

                    default:
                    throw new NotSupportedException("This should never happen");
                }
            }
            #endregion
            #region Convective interface flux (level set two)
            if(Control.IsTwoLevelSetRun) {
                if(Control.LsTwo_SpeciesPairs == null) {
                    throw new ArgumentException("need to provide species pairs for level set two");
                }
                for(int i = 0; i < Control.LsTwo_SpeciesPairs.GetLength(0); i++) {
                    string negSpec = Control.LsTwo_SpeciesPairs[i, 0];
                    string posSpec = Control.LsTwo_SpeciesPairs[i, 1];
                    switch(Control.ConvectiveInterfaceFlux_LsTwo) {
                    case ConvectiveInterfaceFluxes.OptimizedHLLCInterface:
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCFlux_XDG_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //break;

                    case ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var:
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new XESF.Fluxes.OptimizedHLLCDensityFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new XESF.Fluxes.OptimizedHLLCMomentumFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new XESF.Fluxes.OptimizedHLLCEnergyFlux_XDG_Wall(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec));
                    //break;

                    //case ConvectiveInterfaceFluxes.OptimizedHLLCWall:
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Density));
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Momentum, 0));
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Momentum, 1));
                    //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new OptimizedHLLCFlux_XDG_Wall(LsTrk, boundaryMap, material, FluxVariables.Energy));
                    //    break;
                    throw new NotImplementedException();
                    case ConvectiveInterfaceFluxes.GodunovInterface:
                    var Sdensityflux = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    var SmomXflux = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    var SmomYflux = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    var Senergyflux = new XESF.Fluxes.GodunovFlux_Interface(LsTrk, this.Control, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPInterfaceFlux((ILevelSetForm)Sdensityflux, 0));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomXflux, 1));
                    int next_comp = 2;
                    if(Grid.SpatialDimension == 3) {
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux, next_comp));
                        next_comp++;
                    }
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPInterfaceFlux((ILevelSetForm)Senergyflux, next_comp));
                    break;
                    case ConvectiveInterfaceFluxes.CentralFluxInterface:
                    var Sdensityflux2 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Density, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    var SmomXflux2 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    var SmomYflux2 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    var Senergyflux2 = new XESF.Fluxes.CentralFlux_XDG_Interface(boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec);
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPInterfaceFlux((ILevelSetForm)Sdensityflux2, 0));
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomXflux2, 1));
                    int next_comp2 = 2;
                    if(Grid.SpatialDimension == 3) {
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux2, next_comp2));
                        next_comp2++;
                    }
                    this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPInterfaceFlux((ILevelSetForm)Senergyflux2, next_comp2));
                    break;
                    case ConvectiveInterfaceFluxes.RoeInterface:
                        //var Sdensityflux3 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Density, component:0,levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec,s_alpha:10,D:D-1);
                        //var SmomXflux3 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, component: 0, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D-1);
                        //var SmomYflux3 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, component: 1, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D-1);
                        //var SmomZflux3 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Momentum, component: 2, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D-1);
                        //var Senergyflux3 = new XESF.Fluxes.RoeFlux_Interface(LsTrk, boundaryMap, material, XESF.Fluxes.FluxVariables.Energy, component: 0, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: 10, D: D - 1);
                        //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new SPInterfaceFlux((ILevelSetForm)Sdensityflux3, 0));
                        //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomXflux3, 1));
                        //int next_comp3 = 2;
                        //if (Grid.SpatialDimension==3) {
                        //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux3, next_comp3));
                        //    next_comp3++;
                        //}
                        //if(Grid.SpatialDimension == 4) {
                        //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomYflux3, next_comp3));
                        //    next_comp3++;
                        //    this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new SPInterfaceFlux((ILevelSetForm)SmomZflux3, next_comp3));
                        //    next_comp3++;
                        //}
                        //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new SPInterfaceFlux((ILevelSetForm)Senergyflux3, next_comp3));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Density, int.MinValue, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1));
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Momentum, 0, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1));
                        if(Grid.SpatialDimension >= 3) {
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Momentum, 1, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1));
                        }
                        if(Grid.SpatialDimension >= 4) {
                            this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.zComponent].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Momentum, 2, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1));
                        }
                        this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new Fluxes.RoeSTFlux_Interface(LsTrk, boundaryMap, material, FluxVariables.Energy, int.MinValue, levelSetIndex: 1, negSpecies: negSpec, posSpecies: posSpec, s_alpha: Control.flux_s_alpha, D: D - 1));


                        // Interface as wall
                        //this.XSpatialOperator.EquationComponents[CompressibleVariables.Density].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Density));
                        //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Momentum, 0));
                        //this.XSpatialOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Momentum, 1));
                        //this.XSpatialOperator.EquationComponents[CompressibleVariables.Energy].Add(new GodunovFlux_Wall(LsTrk, this.Control, boundaryMap, material, FluxVariables.Energy));
                        break;

                    case ConvectiveInterfaceFluxes.None:
                    break;

                    default:
                    if(Control.IsTwoLevelSetRun) {
                        throw new NotSupportedException("This should never happen");
                    }
                    break;

                }
            }
        }
            #endregion

            this.XSpatialOperator.LinearizationHint = LinearizationHint.FDJacobi;
            this.XSpatialOperator.AgglomerationThreshold = this.Control.AgglomerationThreshold;
            CurrentAgglo = this.Control.AgglomerationThreshold;
            this.XSpatialOperator.Commit();

            if(Control.GetInitialValue == GetInitialValue.FromP0Timestepping) {
                ComputeP0Solution(); //needs to be done here as operator is now assembled
            }

            //untested_ComputeP0Solution();
            #region Residual logging
            // Configure residual handling
            if(L == null && Control.ResidualInterval != 0) {
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
   
            ChooseOptProblem();
            //Initialize empty vectors and matrices
            InitializeMatricesAndVectors();
            //// Cell agglomerator (cell length scales are needed for diffusive AV fluxes)
            UpdateAgglomerator();

            ComputeResiduals();

            InitResNorm = res_l2;
            Init_obj_f = obj_f;
            #endregion
            ResNorms.Add(res_l2);
            obj_f_vals.Add(obj_f);
            UpdateDerivedVariables();
        }

        public override void InitializeMultiGridOpConfig() {
            int p = this.Momentum[0].Basis.Degree;
            int D = this.GridData.SpatialDimension;
            // set the MultigridOperator configuration for each level:
            // it is not necessary to have exactly as many configurations as actual multigrid levels:
            // the last configuration enty will be used for all higher level
            MultiGridOperatorConfig = new MultigridOperator.ChangeOfBasisConfig[3][];
            for(int iLevel = 0; iLevel < MultiGridOperatorConfig.Length; iLevel++) {
                MultiGridOperatorConfig[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D +1];

                // configuration for continuity
                MultiGridOperatorConfig[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { Math.Max(0, p - iLevel) },
                    mode = MultigridOperator.Mode.Eye,
                    VarIndex = new int[] { 0 }
                };


                // configurations for momentum equation
                for(int d = 0; d < D-1; d++) {
                    MultiGridOperatorConfig[iLevel][d + 1] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Math.Max(0, p - iLevel) },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { d + 1 }
                    };
                }

                // configuration for energy
                MultiGridOperatorConfig[iLevel][D ] = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { Math.Max(0, p - iLevel) },
                    mode = MultigridOperator.Mode.Eye,
                    VarIndex = new int[] { D }
                };
            }
        }


        protected override void SetInitial(double t) {
            if(LsTrk.GridDat.SpatialDimension > 2) {
                throw new NotSupportedException("Only valid for 2D simulation runs.");
            }
            CurrentStepNo = 0;
            CellQuadratureScheme scheme = new CellQuadratureScheme(true, CellMask.GetFullMask(this.GridData));
            #region DB Choice
            // Switch to Load the right DB if a DB is loaded
            switch(Control.GetInitialValue) {
                case GetInitialValue.FromDBSinglePhase:
                case GetInitialValue.FromDBXDG:
                DatabaseInfo dbi = DatabaseInfo.Open(Control.ShockLevelSet_Db);
                ISessionInfo si = dbi.Controller.GetSessionInfo(Control.ShockLevelSet_Info.Item1);
                tsiFromDb = Control.IVTimestepNumber >= si.Timesteps.Count ? si.Timesteps.Last() : si.Timesteps.Pick(Control.IVTimestepNumber);//if(Control.ShockLevelSet_Info.Item2.MajorNumber < 0.0) {
                //    tsiFromDb = si.Timesteps.Last();
                //} else {
                //    throw new NotSupportedException("Currently, only the last time step can be used.");
                //}
                break;
                case GetInitialValue.FromAVRun:
                DatabaseInfo dbi_2 = DatabaseInfo.Open(Control.SeedFromAV_Db);
                ISessionInfo si_2 = dbi_2.Controller.GetSessionInfo(Control.SeedFromAV_Db_Info.Item1);
                tsiFromDb = Control.IVTimestepNumber >= si_2.Timesteps.Count ? si_2.Timesteps.Last() : si_2.Timesteps.Pick(Control.IVTimestepNumber);
                //if(Control.SeedFromAV_Db_Info.Item2.MajorNumber < 0.0) {
                //    tsiFromDb = si_2.Timesteps.Last();
                //} else {
                //    throw new NotSupportedException("Currently, only the last time step can be loaded.");
                //}
                break;
                default:
                break;

            }
            #endregion

            #region Initialiize the LevelSet
            if(Control.IsTwoLevelSetRun) {
                this.LevelSet.ProjectField(1.0, this.Control.LevelSetOneInitialValue, scheme);
                LsTBO = LevelSetTwo;
            } else {
                // Level set one
                LsTBO = LevelSet;
            }

            switch(Control.GetLevelSet) {

                case GetLevelSet.FromParams:
                LevelSetOpti.AssembleTransMat(LsTBO);
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                break;

                case GetLevelSet.FromFunction:
                LevelSetOpti.AssembleTransMat(LsTBO);
                if(Control.IsTwoLevelSetRun) {
                    LevelSetOpti.ProjectFromFunction(Control.LevelSetTwoInitialValue);
                } else {
                    LevelSetOpti.ProjectFromFunction(Control.LevelSetOneInitialValue);
                }
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                break;

                case GetLevelSet.FromOldSimulation: //Assuming that we allready Loaded the ShockLevelSet
                case GetLevelSet.FromReconstruction:
                LevelSetOpti.AssembleTransMat(LsTBO);
                LsTBO.Clear();
                this.LsTBO.Acc(1.0, this.ShockLevelSetField);
                LevelSetOpti.ProjectFromLevelSet(LsTBO);
                LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                break;

                case GetLevelSet.FromReconstructionFromPoints: //here we need to have a density field and points to work with

                if(Control.PointPath == null) {
                    throw new NotSupportedException("must specify a PointPath");
                }
                if(Control.GetInitialValue == GetInitialValue.FromFunctionPerSpecies) {
                    throw new NotSupportedException("GetLevelSet.FromReconstructionFromPoints if Initial value is loaded from DB");
                }

                switch(Control.OptiLevelSetType) {
                    // in case of a spline we can directly get the spline from the points (Hopefully)
                    case OptiLevelSetType.SplineLevelSet:
                    //get Points

                    if(LevelSetOpti is SplineOptiLevelSet splineLS) {
                        var points2 = IMatrixExtensions.LoadFromTextFile(Control.PointPath);
                        var points3 = points2.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { points2.GetLength(0) - 1, 1 });
                        splineLS.GetSplineOverDetermined(points3);
                        splineLS.Interpolate();
                        splineLS.ProjectOntoLevelSet(LsTBO);
                    }
                    break;
                    default:
                    LevelSetOpti.AssembleTransMat(LsTBO);

                    //load the density of last timestep
                    var densityField = tsiFromDb.Fields.Where(f => f.Identification == "rho").SingleOrDefault() as SinglePhaseField;
                    var points = IMatrixExtensions.LoadFromTextFile(Control.PointPath);
                    var tmpLS = ShockFindingExtensions.ReconstructLevelSetField(densityField, points);

                    // Continuity Projection does not work
                    //var jCells = new double[points.Lengths[0]];
                    //for(int i=0; i< points.Lengths[0]; i++) {
                    //    IGridData_Extensions.LocatePoint(LsTBO.GridDat, new double[] { points[i, 0], points[i, 1] }, out long GlobalId, out long GlobalIndex, out bool inSide, out bool onthisProcess);
                    //    jCells[i] = GlobalIndex;
                    //}
                    //tmpLS = ShockFindingExtensions.ContinuousLevelSet(tmpLS, jCells);
                    //LsTBO.ProjectFromForeignGrid(1.0, tmpLS);
                    //LevelSetOpti.ProjectFromLevelSet(LsTBO);
                    LevelSetOpti.ProjectFromForeignLevelSet(tmpLS);
                    LevelSetOpti.ProjectOntoLevelSet(LsTBO);
                    //PlotCurrentState(0.0, Grid.NumberOfCells + LsTBO.Basis.Degree, 3);
                    break;
                }
                break;
                case GetLevelSet.DirectyFromTimestep:
                if(LevelSetOpti is SplineOptiLevelSet splineLs) {
                    if(((TimestepProxy)tsiFromDb).GetInternal() is IDTTimeStepInfo IDTtsInfo) {
                        var LSParams = IDTtsInfo.LevelSetParams;
                        if(LSParams.Length == splineLs.m_AllParams.Length) {
                            for(int i = 0; i < LSParams.Length; i++) {
                                splineLs.m_AllParams[i] = LSParams[i];
                            }
                            splineLs.Interpolate();
                            splineLs.ProjectOntoLevelSet(LsTBO);
                            //this.Mus = IDTtsInfo.MuHistory.ToList();
                            //this.Alphas = IDTtsInfo.AlphaHistory.ToList();
                            //this.Gammas = IDTtsInfo.GammaHistory.ToList();
                            //this.ResNorms = IDTtsInfo.ResHistory.ToList();
                            //this.EnResNorms = IDTtsInfo.EnResHistory.ToList();
                            //this.CurrentStepNo = Control.IVTimestepNumber + 1;
                            this.gamma = IDTtsInfo.Gamma;
                        } else {
                            throw new NotSupportedException($"{GetLevelSet.DirectyFromTimestep} not supported if Grid has not equal y - Cells");
                        }

                    }

                    //splineLs.x = points3.To1DArray();

                } else {
                    throw new NotSupportedException($"{Control.GetLevelSet} only supported for SplineOptiLevelSet");
                }
                break;
            }
            //PlotCurrentState("00");
            LsTrk.UpdateTracker(CurrentStepNo);

            #endregion


            #region Initial Values
            // Clear fields (just to be sure)
            this.Density.Clear();
            for(int d=0; d < Grid.SpatialDimension - 2; d++) {
                this.Momentum[d].Clear();
            }
            this.Energy.Clear();

            // Seed from level set reconstruction run or apply initial conditions from the control file
            switch(Control.GetInitialValue) {
                #region FromFunctionPerSpecies
                case GetInitialValue.FromFunctionPerSpecies:
                case GetInitialValue.FromP0Timestepping:
                // Check if initial conditions have been specified in primitive or conservative variables
                bool primitiveVariables = true;
                foreach(string key in Control.InitialValues_Evaluators.Keys) {
                    primitiveVariables = key.Contains(CompressibleVariables.Density)
                      || key.Contains(XESTSFVariables.Velocity.xComponent)
                      || key.Contains(XESTSFVariables.Velocity.yComponent)
                      || key.Contains(XESTSFVariables.Pressure);

                    if(!primitiveVariables) {
                        break;
                    }
                }
                if(primitiveVariables)
                    // Initial conditions are specified in primitive variables
                    foreach(SpeciesId id in this.SpeciesToEvaluate_Ids) {
                        SeedSpecies_PrimVars(id, scheme);
                    }
                else {
                    // Initial conditions are specified in conservative variables
                    foreach(SpeciesId id in this.SpeciesToEvaluate_Ids) {
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

                foreach(string s in Control.SpeciesToEvaluate) {
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

                foreach(string s in Control.SpeciesToEvaluate) {
                    Density.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rhoXDG.GetSpeciesShadowField(s));
                    Momentum[0].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m0XDG.GetSpeciesShadowField(s));
                    Momentum[1].GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, m1XDG.GetSpeciesShadowField(s));
                    Energy.GetSpeciesShadowField(s).ProjectFromForeignGrid(1.0, rhoEXDG.GetSpeciesShadowField(s));
                }
                break;
                #endregion
                #region FromDB Partial SeedInflowExactly
                case GetInitialValue.FromDB_Partial_SeedInflowExactly:
                if(Control.ShockLevelSet_SeedFromDb) {
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

                foreach(string s in Control.SpeciesToEvaluate) {
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

        private void SeedSpecies_PrimVars(SpeciesId id, CellQuadratureScheme scheme) {
            // Helpers
            Func<double[], double> density = Control.InitialValues_Evaluators[CompressibleVariables.Density + "#" + LsTrk.GetSpeciesName(id)];
            Func<double[], double> velocityX = Control.InitialValues_Evaluators[XESTSFVariables.Velocity.xComponent + "#" + LsTrk.GetSpeciesName(id)];
            //Func<double[], double> velocityY = Control.InitialValues_Evaluators[XESTSFVariables.Velocity.yComponent + "#" + LsTrk.GetSpeciesName(id)];
            Func<double[], double> pressure = Control.InitialValues_Evaluators[XESTSFVariables.Pressure + "#" + LsTrk.GetSpeciesName(id)];

            // Density
            this.Density.GetSpeciesShadowField(id).ProjectField(1.0, density);

            // Momentum (rho * momentum vector)
            this.Momentum[0].GetSpeciesShadowField(id).ProjectField(1.0, X => density(X) * velocityX(X));
            //this.Momentum[1].GetSpeciesShadowField(id).ProjectField(1.0, X => density(X) * velocityY(X));

            // Energy
            Energy.GetSpeciesShadowField(id).ProjectField(
                1.0,
                delegate (double[] X) {
                    Vector velocityVec = new Vector(LsTrk.GridDat.SpatialDimension);
                    velocityVec[0] = velocityX(X);
                    //velocityVec[1] = velocityY(X);

                    StateVector state = StateVector.FromPrimitiveQuantities(this.Control.GetMaterial(), density(X), velocityVec, pressure(X));
                    return state.Energy;
                },
                    scheme);
        }

        #region XDG auxiliary functions
        private IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                Dictionary<SpeciesId, IEnumerable<double>> massScaleToSpeciesMap = new Dictionary<SpeciesId, IEnumerable<double>>();

                int D = this.Grid.SpatialDimension;
                double[] ones = new double[D + 2];
                ones.SetAll(1.0);

                foreach(SpeciesId speciesId in this.SpeciesToEvaluate_Ids) {
                    massScaleToSpeciesMap.Add(speciesId, ones);
                }
                return massScaleToSpeciesMap;
            }
        }
        #endregion
        public override void CreateConservativeFields(LevelSetTracker in_LsTrk, int in_DgDegree) {
            #region Optional: load shock level set DG field from a CNS simulation stored in a database
            switch(Control.GetLevelSet) {
                case GetLevelSet.FromOldSimulation:
                case GetLevelSet.FromReconstruction:

                DatabaseInfo dbi = DatabaseInfo.Open(Control.ShockLevelSet_Db);
                //DatabaseInfo dbi = DatabaseInfo.Open(@"C:\BoSSS-experimental\internal\src\private-mag\XESF\Tests\bosss_db_levelSets.zip");
                ISessionInfo si = dbi.Controller.GetSessionInfo(Control.ShockLevelSet_Info.Item1);
                if(Control.ShockLevelSet_Info.Item2.MajorNumber < 0.0) {
                    tsiFromDb = si.Timesteps.Last();
                } else {
                    throw new NotSupportedException("Currently, only the last time step can be used.");
                }

                // Check if grid from shock level set reconstruction equals current grid
                if(tsiFromDb.Grid.ID != this.GridData.GridID) {
                    Console.WriteLine("Different grids used for LEVEL SET RECONSTRUCTION and XESF");
                    IGrid grid = dbi.Controller.DBDriver.LoadGrid(tsiFromDb.GridID, dbi);

                    // Load grid and shock level set field
                    var fieldsFromDb = dbi.Controller.DBDriver.LoadFields(tsiFromDb, grid.iGridData);

                    // Project onto current grid
                    if(Control.IsTwoLevelSetRun) {
                        ShockLevelSetField = new SinglePhaseField(new Basis(this.GridData, Control.LevelSetTwoDegree), "shockLevelSetField");
                    } else {
                        ShockLevelSetField = new SinglePhaseField(new Basis(this.GridData, Control.LevelSetDegree), "shockLevelSetField");
                    }

                    ShockLevelSetField.ProjectFromForeignGrid(1.0, (ConventionalDGField)fieldsFromDb.Pick(1));

                } else {
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
            int D = in_LsTrk.GridDat.SpatialDimension - 1;

            this.Density = new XDGField(new XDGBasis(in_LsTrk, in_DgDegree), CompressibleVariables.Density);
            this.Density.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            XDGField[] momentumFields = new XDGField[D];
            XDGBasis momentumBasis = new XDGBasis(in_LsTrk, in_DgDegree);
            for(int d = 0; d < D; d++) {
                string variableName = CompressibleVariables.Momentum[d];
                momentumFields[d] = new XDGField(momentumBasis, variableName);
                momentumFields[d].UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            }
            this.Momentum = new VectorField<XDGField>(momentumFields);
            this.Energy = new XDGField(new XDGBasis(in_LsTrk, in_DgDegree), CompressibleVariables.Energy);
            this.Energy.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;

            ConservativeFields = new XDGField[D + 2];



            ConservativeFields[0] = Density;
            for(int d = 0; d < D; d++) {
                ConservativeFields[d + 1] = Momentum[d];
            }
            ConservativeFields[D + 1] = Energy;

            #endregion

            #region Create derived fields
            this.DerivedVariableToXDGFieldMap = new Dictionary<DerivedVariable<XDGField>, XDGField>();
            this.DerivedVariableToSinglePhaseFieldMap = new Dictionary<DerivedVariable<SinglePhaseField>, SinglePhaseField>();
            this.DerivedVariableToDoubleMap = new Dictionary<DerivedVariable<double>, double>();
            foreach(KeyValuePair<Variable, int> pair in this.Control.VariableToDegreeMap) {
                Variable variable = pair.Key;
                int degree = pair.Value;

                if(variable is DerivedVariable<XDGField> xdgVar) {
                    this.DerivedVariableToXDGFieldMap.Add(xdgVar, new XDGField(new XDGBasis(in_LsTrk, degree), variable.Name));
                } else if(variable is DerivedVariable<SinglePhaseField> singleVar) {
                    this.DerivedVariableToSinglePhaseFieldMap.Add(singleVar, new SinglePhaseField(new Basis(in_LsTrk.GridDat, degree), variable.Name));
                } else if(variable is DerivedVariable<double> doubleVar) {
                    this.DerivedVariableToDoubleMap.Add(doubleVar, new double());
                }
            }
            #endregion

            #region register the derived Fields
            this.m_IOFields.AddRange(DerivedVariableToXDGFieldMap.Values);
            this.m_IOFields.AddRange(DerivedVariableToSinglePhaseFieldMap.Values);

            this.m_RegisteredFields.AddRange(DerivedVariableToXDGFieldMap.Values);
            this.m_RegisteredFields.AddRange(DerivedVariableToSinglePhaseFieldMap.Values);
            #endregion 

        }
        

        public override void UpdateDerivedVariables() {
            // Update derived variables (including sensor _variable_)
            foreach(KeyValuePair<DerivedVariable<XDGField>, XDGField> pair in this.DerivedVariableToXDGFieldMap) {
                pair.Key.UpdateFunction(pair.Value, this);
            }

            foreach(KeyValuePair<DerivedVariable<SinglePhaseField>, SinglePhaseField> pair in this.DerivedVariableToSinglePhaseFieldMap) {
                pair.Key.UpdateFunction(pair.Value, this);
            }

            foreach(KeyValuePair<DerivedVariable<double>, double> pair in this.DerivedVariableToDoubleMap) {
                using(var tr = new FuncTrace()) {
                    pair.Key.UpdateFunction(pair.Value, this);
                    tr.Info(pair.Key.ToString() + "=" + pair.Value.ToString());
                }
            }
        }

    }
}




