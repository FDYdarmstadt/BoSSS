// See https://aka.ms/new-console-template for more information
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.EllipticExtension;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XNSECommon;
using FreeXNSE.Test;
using ilPSP.Tracing;
using MPI.Wrappers;
using ilPSP.Utils;
using Newtonsoft.Json;
using System;
using System.Reflection;
using System.Runtime.Serialization;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.LoadBalancing;
using BoSSS.Solution.Tecplot;
using ilPSP;
using System.Diagnostics;
using BoSSS.Solution.Statistic.QuadRules;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.LevelSetTools;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using System.Reflection.Metadata;

namespace FreeXNSE {

    public class FreeXNSE : SolverWithLevelSetUpdater<FreeXNSE_Control> {
        static void Main(string[] args) {

            //InitMPI();
            //DeleteOldPlotFiles();
            //FreeXNSE_test.FreeXNSE_ChannelTest();
            //FreeXNSE_test.FreeXNSE_CircleTest();
            //FreeXNSE_test.FreeXNSE_EllipseTest();
            //FreeXNSE_test.FreeXNSE_EllipseTestParameterized();
            //FreeXNSE_test.FreeXNSE_TestParameterizedAdvection();
            //FreeXNSE_test.FreeXNSE_EllipsoidTest();
            //FreeXNSE_Contactline_test.FreeXNSE_SlugInChannel_Equilibrium();
            //FreeXNSE_Contactline_test.FreeXNSE_SlugInChannel_Couette();
            //FreeXNSE_Contactline_test.FreeXNSE_StaticDroplet_FixedInterface();
            //FreeXNSE_Contactline_test.FreeXNSE_StaticDroplet_DynamicContactAngle();
            //FreeXNSE_Contactline_test.FreeXNSE_PlanarInterface_DynamicContactAngle();
            //FreeXNSE_Contactline_test.FreeXNSE_SlidingDroplet_TiltedPlane(Math.PI/180 * 30);
            //FreeXNSE_Contactline_test.FreeXNSE_SlidingDroplet_ContactAngleHysteresis(105.0/180.0*Math.PI, 0.0);
            //FreeXNSE_Contactline_test.FreeXNSE_SlidingDroplet_ContactAngleHysteresis(Math.PI, 89.0 / 180.0 * Math.PI);
            //FreeXNSE_Contactline_test.AggregationFail();

            //FinalizeMPI();
            //System.Environment.Exit(0);

            {
                FreeXNSE._Main(args, false, delegate () {
                    var p = new FreeXNSE();
                    return p;
                });
            }
        }
        protected override int NoOfLevelSets => Control.UseImmersedBoundary ? 2 : 1;

        /// <summary>
        /// Either fluids A and B; or A, B and solid C.
        /// </summary>
        protected override Array SpeciesTable {
            get {
                if(Control.UseImmersedBoundary) {
                    var r = new string[2, 2];
                    r[0, 0] = SpeciesNames[0];
                    r[0, 1] = SpeciesNames[2]; // solid
                    r[1, 0] = SpeciesNames[1];
                    r[1, 1] = SpeciesNames[2]; // also solid
                    return r;
                } else {
                    return SpeciesNames;
                }
            }
        }

        public virtual string[] SpeciesNames {
            get {
                if(Control.UseImmersedBoundary) {
                    return new[] { "A", "B", "C" };
                } else {
                    return new[] { "A", "B" };
                }
            }
        }

        public override int QuadOrder() {
            if(Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye) {
                throw new ArgumentException($"Please use Saye-Rules!");
            }

            //QuadOrder
            int degU = Control.Degree;
            int quadOrder = degU * (Control.ActiveTerms.Convective != Convective.Off ? 3 : 2);
            if(this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                //See remarks
                quadOrder *= 2;
                quadOrder += 1;
            }

            return quadOrder;
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            using(var tr = new FuncTrace()) {
                int pVel = Control.Degree;
                int pPrs = Control.EqualOrder ? Control.Degree : Control.Degree - 1;
                int D = this.GridData.SpatialDimension;

                tr.Info($"pre-precond, level {iLevel}: using {MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite} and {MultigridOperator.Mode.IdMass_DropIndefinite}");

                // configurations for velocity
                for(int d = 0; d < D; d++) {
                    var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pVel },
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.VelocityVector(D)[d]) }
                    };
                    configsLevel.Add(configVel_d);
                }
                // configuration for pressure
                var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pPrs },
                    //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) }, // p-multigrid reduction
                    mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Pressure) }
                };
                configsLevel.Add(configPres);

            }
        }

        private FreeXNSE_BoundaryCondMap m_boundaryMap;

        /// <summary>
        /// Relation between 
        /// - edge tags (<see cref="Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>, passed to equation components via <see cref="BoSSS.Foundation.CommonParams.EdgeTag"/>)
        /// - boundary conditions specified in the control object (<see cref="AppControl.BoundaryValues"/>)
        /// </summary>
        protected FreeXNSE_BoundaryCondMap boundaryMap {
            get {
                if(m_boundaryMap == null)
                    m_boundaryMap = new FreeXNSE_BoundaryCondMap(this.GridData, this.Control.BoundaryValues, new string[] { "A", "B" });
                return m_boundaryMap;
            }
        }

        /// <summary>
        /// dirty hack...
        /// </summary>
        protected override IncompressibleBoundaryCondMap GetBcMap() {
            return boundaryMap.TranslateToLegacy(this.GridData);
        }        

        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            int D = GridData.SpatialDimension;

            if(iLevSet == 0) {
                // +++++++++++++++++++
                // the fluid interface 
                // +++++++++++++++++++

                // averaging at interface:
                ILevelSetParameter levelSetVelocity = new FreeSurfaceVelocity(VariableNames.LevelSetCG, D, Control.Degree);
                return levelSetVelocity;
            } else if(iLevSet == 1) {
                // +++++++++++++++++++++
                // the immersed boundary
                // +++++++++++++++++++++

                throw new NotImplementedException();
                //string[] VelocityNames = VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(1), VariableNames.VelocityVector(D)).ToArray();
                //ScalarFunctionTimeDep[] VelFuncs = new ScalarFunctionTimeDep[D];
                //for(int d = 0; d < D; d++) {
                //    Control.InitialValues_EvaluatorsVec.TryGetValue(VelocityNames[d], out VelFuncs[d]);
                //}

                //ILevelSetParameter levelSetVelocity = new PrescribedLevelSetVelocity(VariableNames.LevelSetCGidx(1), VelFuncs);
                //return levelSetVelocity;
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        protected override XDifferentialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {

            OperatorFactory opFactory = new OperatorFactory();
            GetBcMap();

            #region Momentum     

            for(int d = 0; d < D; d++) {
                opFactory.AddEquation(new BulkMomentum(d, D, SpeciesNames, boundaryMap, Control));
                opFactory.AddEquation(new InterfaceMomentum(d, D, SpeciesNames, boundaryMap, Control));
                //opFactory.AddEquation(new BulkDummy(SpeciesNames[1],EquationNames.MomentumEquationComponent(d), VariableNames.Velocity_d(d)));
            }

            #endregion

            #region Conti

            opFactory.AddEquation(new BulkContinuity(D, SpeciesNames, boundaryMap, Control));
            opFactory.AddEquation(new InterfaceContinuity(D, SpeciesNames, boundaryMap, Control));
            //opFactory.AddEquation(new BulkDummy(SpeciesNames[1], EquationNames.ContinuityEquation, VariableNames.Pressure));


            #endregion

            #region additional parameters

            if(Control.ActiveTerms.Convective == Convective.LaxFriedrich) {
                Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, QuadOrder());
                opFactory.AddParameter(v0Mean);
                levelSetUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, v0Mean);
            }

            if(Control.ActiveTerms.Viscous != Viscous.Off) {
                opFactory.AddCoefficient(new Bulkfriction(D, Control.DimensionlessNumbers.beta, Control.SlipScaling, levelSetUpdater));
            }

            if(Control.ActiveTerms.Viscous != Viscous.Off || Control.EqualOrder) {
                opFactory.AddCoefficient(new Bulkviscosityfield(Control.ViscosityScaling, levelSetUpdater));
            }

            Normals normalsParameter = new Normals(D, ((LevelSet)levelSetUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);
            levelSetUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, normalsParameter);

            if(Control.ActiveTerms.SurfaceTension != SurfaceTension.Off) {
                opFactory.AddCoefficient(new Contactlinefriction(D, Control.DimensionlessNumbers.alpha, Control.ContactAngleScaling));
                opFactory.AddCoefficient(new Contactangle(D, Control.DimensionlessNumbers.Theta));
                opFactory.AddCoefficient(new Contactangle(D, Control.DimensionlessNumbers.ThetaAdv, 1));
                opFactory.AddCoefficient(new Contactangle(D, Control.DimensionlessNumbers.ThetaRec, -1));
                opFactory.AddCoefficient(new Surfacetensionfield());
                if(Control.ActiveTerms.SurfaceTension == SurfaceTension.LaplaceBeltrami_BoussinesqScriven) {
                    opFactory.AddCoefficient(new Surfaceviscosityfields());
                }
            }

            #endregion

            XDifferentialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            #region Final Settings
            if(this.Control.FixedInterface) {
                Console.WriteLine("Solving with fixed interface position");
                InterfaceComponent_FreeXNSE.FixedInterface = this.Control.FixedInterface;
            }

            /// <summary>
            /// Misc adjustments to the spatial operator before calling <see cref="ISpatialOperator.Commit"/>
            /// </summary>
            XOP.FreeMeanValue[VariableNames.Pressure] = levelSetUpdater.Tracker.Regions.GetCutCellMask().NoOfItemsLocally.MPISum() == 0 & !GetBcMap().DirichletPressureBoundary; // no Neumann boundary present.
            XOP.IsLinear = (this.Control.ActiveTerms.Convective == Convective.Off && (this.Control.ActiveTerms.SurfaceTension == SurfaceTension.Off || this.Control.ContactAngleScaling == 1.0));
            XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;

            if(XOP.IsLinear == true) {
                XOP.LinearizationHint = LinearizationHint.AdHoc;
            } else {
                XOP.LinearizationHint = LinearizationHint.GetJacobiOperator;
            }

            // elementary checks on operator
            if(XOP.CodomainVar.IndexOf(EquationNames.ContinuityEquation) != D)
                throw new ApplicationException("Operator configuration messed up.");
            if(XOP.DomainVar.IndexOf(VariableNames.Pressure) != D)
                throw new ApplicationException("Operator configuration messed up.");
            for(int d = 0; d < D; d++) {
                if(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationComponent(d)) != d)
                    throw new ApplicationException("Operator configuration messed up.");
                if(XOP.DomainVar.IndexOf(VariableNames.Velocity_d(d)) != d)
                    throw new ApplicationException("Operator configuration messed up.");
            }
            PrintConfiguration();
            
            void PrintConfiguration() {
                if(!ConfigurationPrinted) {
                    ConfigurationPrinted = true;
                    switch(this.Control.DimensionlessNumbers.Oh) {
                        case double Oh when Oh > 1: {
                            Console.WriteLine("Value of Ohnesorge number : {0} - using viscous-inertial timescale", Oh);
                            double hmin = this.GridData.iGeomCells.h_min.Min();
                            double tmin = this.Control.dtFixed;
                            if(tmin > this.Control.DimensionlessNumbers.Oh * Math.Sqrt(Math.Pow(hmin, 3)))
                                Console.WriteLine("Possible breach of capillary timestep detected");
                            break;
                        }
                        case double Oh when Oh <= 1: {
                            Console.WriteLine("Value of Ohnesorge number : {0} - using capillary-inertial timescale", Oh);
                            double hmin = this.GridData.iGeomCells.h_min.Min();
                            double tmin = this.Control.dtFixed;
                            if(tmin > Math.Sqrt(Math.Pow(hmin, 3)))
                                Console.WriteLine("Possible breach of capillary timestep detected");
                            break;
                        }
                    }
                    Console.WriteLine("=============== {0} ===============", "TODO : Control Configuration");
                }
            }
            #endregion

            XOP.Commit();
                return XOP;
        }
        static bool ConfigurationPrinted = false;

        /// <summary>
        /// Plot using Tecplot
        /// </summary>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {

            //// Cells Numbers - Local
            //var CellNumbers = this.m_RegisteredFields.Where(s => s.Identification == "CellNumbers").SingleOrDefault();
            //if(CellNumbers == null) {
            //    CellNumbers = new SinglePhaseField(new Basis(this.GridData, 0), "CellNumbers");
            //    this.RegisterField(CellNumbers);
            //}
            //CellNumbers.Clear();
            //CellNumbers.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //    int K = result.GetLength(1); // No nof Nodes
            //    for(int j = 0; j < Len; j++) {
            //        for(int k = 0; k < K; k++) {
            //            result[j, k] = j0 + j;
            //        }
            //    }
            //}, new CellQuadratureScheme());

            List<DGField> Fields2Plot = new List<DGField>();
            foreach(var field in this.m_RegisteredFields) {
                if(field is XDGField xField) {
                    Fields2Plot.Add(xField.GetSpeciesShadowField(SpeciesNames[0]));                    
                } else {
                    Fields2Plot.Add(field);
                }
            }
            Tecplot.PlotFields(Fields2Plot, this.GetType().Name.Split('`').First() + "-" + timestepNo, physTime, superSampling);
        }


        Queue<double> Etot = new Queue<double>();
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            if(Etot.Count == 0) {
                Etot.Enqueue(this.GetKineticEnergy() + this.GetSurfaceEnergy());
                Console.WriteLine("Etot: {0}", Etot.Peek());
            }

            double ret = base.RunSolverOneStep(TimestepNo, phystime, dt);         


            foreach(XDGField xfield in this.CurrentState.Fields) {
                for(int i = 1; i < this.SpeciesNames.Length; i++) {
                    string spc = this.SpeciesNames[i];
                    double norm = xfield.GetSpeciesShadowField(spc).L2Norm();
                    if(norm >= BLAS.MachineEps) { throw new ApplicationException("Inactive Phases should stay zero at all times!"); }
                }
            }

            // instantaneous energy changerate

            // kinetic energy
            double Ekin = this.GetKineticEnergy();
            // surface energy
            double Esurf = this.GetSurfaceEnergy();
            // total energy
            Etot.Enqueue(Ekin + Esurf);

            // effective energy changerate
            double dEtot = Etot.Last() - Etot.Dequeue();

            Console.WriteLine("Ekin: {0}, Esurf: {1}, Etot: {2}, dEtot: {3}", Ekin, Esurf, Etot.Single(), dEtot);

            return ret;
        }        

        protected override ILevelSetEvolver CustomEvolver(int iLevSet) {
            var LevelSetCG = this.LevelSetNames[iLevSet].ContLs;

            if(!Control.ParameterLevelSet.IsNullOrEmpty()) {
                var parameter = new ParameterizedLevelSetEvolver(
                                LevelSetCG,
                                QuadOrder(),
                                this.Grid.SpatialDimension,
                                this.GridData);

                return parameter;
            } else if(Control.DualSplinePhi0Initial != null) {
                int nodeCount = Control.NoOfNodes;
                var SplineEvolver = new DualSplineLevelSetEvolver(LevelSetCG, (GridData)(this.GridData));
                return SplineEvolver;
            } else if(Control.SemiCircleSplinePhi0Initial != null) {
                int nodeCount = Control.NoOfNodes;
                var SplineEvolver = new SemiCircleSplineLevelSetEvolver(LevelSetCG, (GridData)(this.GridData));
                return SplineEvolver;
            } else {
                throw new NotImplementedException();
            }
        }

        protected override LevelSet CustomLevelSet(int iLevSet) {
            var LevelSetCG = this.LevelSetNames[iLevSet].ContLs;
            var LevelSetDG = this.LevelSetNames[iLevSet].DgLs;
            
            if(!Control.ParameterLevelSet.IsNullOrEmpty()) {
                ParameterizedLevelSet parameterLevelSetDG = new ParameterizedLevelSet(this.Control.ParameterLevelSet, new Basis(this.Grid, Control.FieldOptions[LevelSetCG].Degree), LevelSetDG);
                return parameterLevelSetDG;
            } else if(Control.DualSplinePhi0Initial != null) {
                DualSplineLevelSet levelSetDG = new DualSplineLevelSet(Control, new Basis(this.Grid, Control.FieldOptions[LevelSetCG].Degree), LevelSetDG, Control.NoOfNodes);
                return levelSetDG;
            } else if(Control.SemiCircleSplinePhi0Initial != null) {
                SemiCircleSplineLevelSet levelSetDG = new SemiCircleSplineLevelSet(Control, new Basis(this.Grid, Control.FieldOptions[LevelSetCG].Degree), LevelSetDG, Control.NoOfNodes);
                return levelSetDG;
            } else {
                throw new NotImplementedException();
            }
        }

        protected override void CustomInitializeLevelSet(DualLevelSet pair) {
            if(!Control.ParameterLevelSet.IsNullOrEmpty()) {
                return;
            } else if(Control.DualSplinePhi0Initial != null) {
                DualSplineLevelSet SplineLevelSet = new DualSplineLevelSet(Control, new Basis(this.Grid, Control.FieldOptions[pair.CGLevelSet.Identification].Degree), VariableNames.LevelSetDG, Control.NoOfNodes);
                pair.DGLevelSet = SplineLevelSet;
            } else if(Control.SemiCircleSplinePhi0Initial != null) {
                Console.WriteLine("Careful! very experimental!");
                SemiCircleSplineLevelSet SplineLevelSet = new SemiCircleSplineLevelSet(Control, new Basis(this.Grid, Control.FieldOptions[pair.CGLevelSet.Identification].Degree), VariableNames.LevelSetDG, Control.NoOfNodes);
                pair.DGLevelSet = SplineLevelSet;
            } else {
                throw new NotImplementedException();
            }
        }

        /* Override of the level set related methods, unnecessary, see implementation of LevelSetEvolution.CustomLevelSet

        /// <summary>
        /// Instantiate the level-set-system (fields for storing, evolution operators, ...) 
        /// Before creating XDG-fields one need to
        /// initialize the level-set fields <see cref="InitializeLevelSets"/>
        /// </summary>
        protected override LevelSetUpdater InstantiateLevelSetUpdater() {
            if(!this.GridData.IsAlive())
                throw new ApplicationException("invalid grid -- most likely something went wrong during mesh adaptation/redistribution");
            int D = this.Grid.SpatialDimension;
            var lsNames = this.LevelSetNames;
            //ISpatialOperator test = this.Operator; hier ist noch kein OP
            int NoOfLevelSets = lsNames.Length;
            if(NoOfLevelSets != this.NoOfLevelSets)
                throw new ApplicationException();
            //bool isRestart = Control.RestartInfo != null;


            // phase 1: create DG level-sets
            // ======================================
            LevelSet[] DGlevelSets = new LevelSet[NoOfLevelSets];
            for(int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                var LevelSetCG = lsNames[iLevSet].ContLs;
                var LevelSetDG = lsNames[iLevSet].DgLs;

                int levelSetDegree = Control.FieldOptions[LevelSetCG].Degree;    // need to change naming convention of old XNSE_Solver

                switch(Control.Get_Option_LevelSetEvolution(iLevSet)) {
                    case LevelSetEvolution.CustomLevelSet:
                        ParameterizedLevelSet parameterLevelSetDG = new ParameterizedLevelSet(this.Control.ParameterLevelSet, new Basis(GridData, levelSetDegree), LevelSetDG);
                        DGlevelSets[iLevSet] = parameterLevelSetDG;
                    break;
                    case LevelSetEvolution.Prescribed:
                    case LevelSetEvolution.StokesExtension:
                    case LevelSetEvolution.Phasefield: {
                        LevelSet levelSetDG = new LevelSet(new Basis(GridData, levelSetDegree), LevelSetDG);
                        DGlevelSets[iLevSet] = levelSetDG;
                        break;
                    }
                    default:
                        throw new NotImplementedException($"Unknown option for level-set evolution: {Control.Option_LevelSetEvolution}");
                }

            }


            // phase 2: create updater
            // =======================
            LevelSetUpdater lsUpdater;
            switch(NoOfLevelSets) {
                case 1:
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, Control.LS_TrackerWidth,
                    (string[])this.SpeciesTable,
                    base.GetLsUpdaterInputFields,
                    DGlevelSets[0], lsNames[0].ContLs, Control.LSContiProjectionMethod);
                break;

                case 2:
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, Control.LS_TrackerWidth,
                    (string[,])this.SpeciesTable,
                    base.GetLsUpdaterInputFields,
                    DGlevelSets[0], lsNames[0].ContLs, DGlevelSets[1], lsNames[1].ContLs, Control.LSContiProjectionMethod);
                break;

                default:
                throw new NotImplementedException("Unsupported number of level-sets: " + NoOfLevelSets);
            }


            // phase 3: instantiate evolvers
            // ============================
            for(int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                var LevelSetCG = lsNames[iLevSet].ContLs;
                var LevelSetDG = lsNames[iLevSet].DgLs;

                // create evolver:
                switch(Control.Get_Option_LevelSetEvolution(iLevSet)) {
                    case LevelSetEvolution.CustomLevelSet: {
                        var parameter = new ParameterizedLevelSetEvolver(
                            VariableNames.LevelSetCG,
                            QuadOrder(),
                            D,
                            this.GridData);

                        lsUpdater.AddEvolver(LevelSetCG, parameter);
                        break;
                    }
                    case LevelSetEvolution.StokesExtension: {
                        ILevelSetEvolver stokesExtEvo;
                        if(LevelSetHandling == BoSSS.Solution.XdgTimestepping.LevelSetHandling.Coupled_Iterative) {
                            stokesExtEvo = new ImplicitStokesExtensionEvolver(LevelSetCG, QuadOrder(), D,
                            GetBcMap(),
                            this.Control.AgglomerationThreshold, this.GridData);
                        } else {
                            stokesExtEvo = new StokesExtensionEvolver(LevelSetCG, QuadOrder(), D,
                            GetBcMap(),
                            this.Control.AgglomerationThreshold, this.GridData,
                            ReInitPeriod: Control.ReInitPeriod);
                        }
                        lsUpdater.AddEvolver(LevelSetCG, stokesExtEvo);
                        break;
                    }
                    case LevelSetEvolution.Phasefield: {
                        var PhasefieldEvolver = new PhasefieldEvolver(LevelSetCG, QuadOrder(), D,
                            GetBcMap(), this.Control,
                            this.Control.AgglomerationThreshold, this.GridData);

                        lsUpdater.AddEvolver(LevelSetCG, PhasefieldEvolver);
                        break;
                    }
                    case LevelSetEvolution.Prescribed: {
                        var prescrEvo = new PrescribedEvolver(this.Control.InitialValues_EvaluatorsVec[LevelSetCG]);
                        lsUpdater.AddEvolver(LevelSetCG, prescrEvo);
                        break;
                    }
                    case LevelSetEvolution.None: {
                        break;
                    }
                    default:
                    throw new NotImplementedException($"Unknown option for level-set evolution: {Control.Option_LevelSetEvolution}");
                }

                // add velocity parameter:
                var levelSetVelocity = GetLevelSetVelocity(iLevSet);
                if(levelSetVelocity != null) {
                    if(!ArrayTools.ListEquals(levelSetVelocity.ParameterNames,
                        BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(LevelSetCG, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)))) {
                        throw new ApplicationException($"Parameter names for the level-set velocity provider for level-set #{iLevSet} ({LevelSetCG}) does not comply with convention.");
                    }
                    lsUpdater.AddLevelSetParameter(LevelSetCG, levelSetVelocity);
                }

            }

            // return
            // ======
            return lsUpdater;
        }

        /// <summary>
        /// Corresponding to <see cref="LevelSetEvolution"/> initialization of LevelSetDG
        /// and projection on continuous LevelSetCG
        /// calls <see cref="LevelSetTracker.UpdateTracker(double, int, bool, int[])">
        /// </summary>
        protected override void InitializeLevelSets(LevelSetUpdater lsUpdater, double time) {

            var lsNames = this.LevelSetNames;
            int NoOfLevelSets = lsNames.Length;
            if(NoOfLevelSets != lsUpdater.LevelSets.Count)
                throw new ApplicationException();


            for(int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                string LevelSetDG = lsNames[iLevSet].DgLs;
                string LevelSetCG = lsNames[iLevSet].ContLs;
                DualLevelSet pair = LsUpdater.LevelSets[LevelSetCG];

                int levelSetDegree = Control.FieldOptions[LevelSetCG].Degree;    // need to change naming convention of old XNSE_Solver
                if(levelSetDegree != pair.DGLevelSet.Basis.Degree)
                    throw new ApplicationException();


                ScalarFunction Phi_InitialValue = null;
                if(Control.InitialValues_EvaluatorsVec.TryGetValue(LevelSetCG, out var scalarFunctionTimeDep)) {
                    Phi_InitialValue = scalarFunctionTimeDep.SetTime(0.0);
                }

                switch(Control.Get_Option_LevelSetEvolution(iLevSet)) {
                    case LevelSetEvolution.CustomLevelSet: {
                        // already done when instantiating the updater!
                        break;
                    }
                    case LevelSetEvolution.Prescribed:
                    case LevelSetEvolution.StokesExtension:
                    case LevelSetEvolution.None: {
                        pair.DGLevelSet.Clear();
                        if(Phi_InitialValue != null)
                            pair.DGLevelSet.ProjectField(Control.InitialValues_EvaluatorsVec[LevelSetCG].SetTime(time));
                        break;
                    }
                    case LevelSetEvolution.Phasefield: {
                        pair.DGLevelSet.Clear();
                        if(Phi_InitialValue != null)
                            pair.DGLevelSet.ProjectField(Control.InitialValues_EvaluatorsVec[LevelSetCG].SetTime(time));

                        break;
                    }
                    default:
                    throw new NotImplementedException($"Unknown option for level-set evolution: {Control.Option_LevelSetEvolution}");
                }

                if(pair.DGLevelSet.L2Norm() == 0.0) {
                    Console.WriteLine($"Level-Set field {LevelSetCG} is **exactly** zero: setting entire field to -1.");
                    pair.DGLevelSet.AccConstant(-1.0);
                }

                pair.CGLevelSet.Clear();
                pair.CGLevelSet.AccLaidBack(1.0, pair.DGLevelSet);

            }

            LsUpdater.Tracker.UpdateTracker(time); // update the tracker **before** pushing

            LsUpdater.Tracker.PushStacks();

            LsUpdater.Tracker.UpdateTracker(time);

            //LsUpdater.Tracker.PushStacks();

        }

        */
    }

    internal static class FreeXNSE_utils {

        internal static double[,] Projection<T>(this T vec) where T : IList<double> {
            int D = vec.Count;
            double[,] P = new double[D, D];

            for(int d = 0; d < D; d++) {
                for(int dd = 0; dd < D; dd++) {
                    if(dd == d)
                        P[d, dd] = (1 - vec[d] * vec[dd]);
                    else
                        P[d, dd] = (0 - vec[d] * vec[dd]);
                }
            }

            return P;
        }

        internal static double GetKineticEnergy(this FreeXNSE solver) {
            using(new FuncTrace()) {
                var lsTrk = solver.LsTrk;
                var SchemeHelper = lsTrk.GetXDGSpaceMetrics(lsTrk.SpeciesIdS.ToArray(), solver.QuadOrder()).XQuadSchemeHelper;

                double kinE = 0.0;

                string spc = solver.SpeciesNames[0];
                var scheme = SchemeHelper.GetVolumeQuadScheme(lsTrk.GetSpeciesId(spc));

                int D = solver.GridData.SpatialDimension;
                for(int d = 0; d < D; d++) {
                    DGField U = solver.CurrentState.Fields.Where(f => f.Identification == VariableNames.VelocityVector(D)[d]).Single();
                    if(U is XDGField) {
                        if(!object.ReferenceEquals((U as XDGField).Basis.Tracker, lsTrk))
                            throw new ArgumentException();
                        U = (U as XDGField).GetSpeciesShadowField(spc);
                    }
                    kinE += U.L2Error(null, solver.QuadOrder(), scheme).Pow2() * 0.5;
                }

                Debug.Assert(kinE >= 0.0);

                return kinE;
            }
        }
        internal static double GetSurfaceEnergy(this FreeXNSE solver) {
            using(new FuncTrace()) {
                var lsTrk = solver.LsTrk;
                var SchemeHelper = lsTrk.GetXDGSpaceMetrics(lsTrk.SpeciesIdS.ToArray(), solver.QuadOrder()).XQuadSchemeHelper;

                double surfE = 0.0;

                var scheme = SchemeHelper.GetLevelSetquadScheme(0, lsTrk.Regions.GetCutCellMask4LevSet(0));
                CellQuadrature.GetQuadrature(new int[] { 1 }, lsTrk.GridDat,
                    scheme.Compile(lsTrk.GridDat, solver.QuadOrder()),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        EvalResult.SetAll(1.0);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++)
                            surfE += ResultsOfIntegration[i, 0];
                    }
                ).Execute();

                Debug.Assert(surfE >= 0.0);

                return 1.0 / solver.Control.DimensionlessNumbers.We * surfE;
            }
        }
    }



   
}

