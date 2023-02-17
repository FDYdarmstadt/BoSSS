using BoSSS.Application.XNSFE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// Solver using the mass fraction approach for calculation of combustion
    /// </summary>
    public class XNSEC_MixtureFraction : XNSEC {

        /// <summary>
        /// dirty hack...
        /// </summary>
        protected override IncompressibleBoundaryCondMap GetBcMap() {
            if (boundaryMap == null)
                base.boundaryMap = new LowMachMixtureFractionMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, new string[] { "A", "B" });
            return boundaryMap;
        }

        public override void DefineScalarEquations(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater) {
            Console.WriteLine("==================");
            Console.WriteLine("Solving the mixture fraction equations!");
            Console.WriteLine("==================");

            //================================
            // Mixture fraction
            //================================
            opFactory.AddEquation(new BulkMixtureFraction_MF("A", D, boundaryMap, config, EoS_A));
            opFactory.AddEquation(new BulkMixtureFraction_MF("B", D, boundaryMap, config, EoS_B));
            opFactory.AddEquation(new MixtureFractionInterface_MF(config, D, EoS_A, EoS_B));

            opFactory.AddParameter(new ThermodynamicPressure(Control.InitialMass, Control.ThermodynamicPressureMode, EoS_A));
        }

        protected override void DefineMomentumEquations(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int d, int D, LevelSetUpdater lsUpdater) {
            opFactory.AddEquation(new BulkNavierStokes_MF("A", d, D, boundaryMap, config, EoS_A));
            opFactory.AddEquation(new BulkNavierStokes_MF("B", d, D, boundaryMap, config, EoS_B));
            opFactory.AddEquation(new NSEInterface_MF("A", "B", d, D, boundaryMap, config, EoS_A, EoS_B, NoOfChemComp: -1));
            if (config.isEvaporation) {
                opFactory.AddEquation(new InterfaceNSE_Evaporation_MF("A", "B", D, d, config));
            }
        }

        protected override void DefineContinuityEquation(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater) {
            opFactory.AddEquation(new BulkContinuity_MF(D, "A", config, boundaryMap, EoS_A, Control.dtFixed));
            opFactory.AddEquation(new BulkContinuity_MF(D, "B", config, boundaryMap, EoS_B, Control.dtFixed));
            //=== evaporation extension === //
            if (config.isEvaporation) {
                opFactory.AddEquation(new InterfaceContinuity_Evaporation_Newton_LowMach_MF("A", "B", D, config));
            } else {
                opFactory.AddEquation(new InterfaceContinuityLowMach(config, D, LsTrk, config.isMatInt));

            }
        }

        /// <summary>
        /// Definition of the boundary condition on the immersed boundary (fluid-solid boundary, level-set 1),
        /// <see cref="XNSE_Control.UseImmersedBoundary"/>;
        /// Override to customize.
        /// </summary>
        protected override void DefineSystemImmersedBoundary(int D, OperatorFactory opFactory, XNSEC_OperatorConfiguration config, LevelSetUpdater lsUpdater) {
            //////////////////////////////////////////////////
            /////     Momentum
            //////////////////////////////////////////////////
            for (int d = 0; d < D; ++d) {
                // so far only no slip!
                opFactory.AddEquation(new NSEimmersedBoundary_Newton_LowMach("A", "C", 1, d, D, boundaryMap, config, EoS_A, config.isMovingMesh, Control.physicsMode));
                opFactory.AddEquation(new NSEimmersedBoundary_Newton_LowMach("B", "C", 1, d, D, boundaryMap, config, EoS_B, config.isMovingMesh, Control.physicsMode));
            }
            ////////////////////////////////////////////////
            ///     Conti
            ////////////////////////////////////////////////
            opFactory.AddEquation(new ImmersedBoundaryContinuity_LowMach("A", "C", 1, config, D));
            opFactory.AddEquation(new ImmersedBoundaryContinuity_LowMach("B", "C", 1, config, D));
            ////////////////////////////////////////////////
            ///     Mixture Fraction
            ////////////////////////////////////////////////
            opFactory.AddEquation(new ImmersedBoundaryMixtureFraction_LowMach("A", "C", 1, D, config, EoS_A));
            opFactory.AddEquation(new ImmersedBoundaryMixtureFraction_LowMach("B", "C", 1, D, config, EoS_B));

            opFactory.AddParameter((ParameterS)GetLevelSetVelocity(1));
        }

        protected override void DefineAditionalParameters(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater, int quadOrder) {
            // ============================== //
            // === additional parameters === //
            // ============================= //
            if (config.PlotAdditionalParameters) {
                opFactory.AddParameter(new DensityMF(EoS_A, EoS_B, config.NoOfChemicalSpecies));
                //opFactory.AddParameter(new Viscosity(EoS_A, EoS_B));
                //opFactory.AddParameter(new HeatCapacity(EoS_A, EoS_B));
            }

            var v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
            if (config.physParams.IncludeConvection && config.isTransport) {                
                opFactory.AddParameter(v0Mean);
            }
            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, v0Mean);

            if (config.isEvaporation) {
                var MassFluxExt = new MassFluxExtension_Evaporation_MixtureFraction(config, Control.HeatRelease, Control.zSt);
                lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, MassFluxExt);
            }
            opFactory.AddCoefficient(new EvapMicroRegion());
            if (config.prescribedMassflux != null)
                opFactory.AddCoefficient(new PrescribedMassFlux(config));
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            int D = this.GridData.SpatialDimension;
            int pVel = VelocityDegree();
            int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
            int pMxtFrctn = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.MixtureFraction].Degree;

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

            // configuration for mixture fraction
            var confTemp = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pMxtFrctn },
                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                //mode = MultigridOperator.Mode.Eye,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.MixtureFraction) }
            };
            configsLevel.Add(confTemp);
        }

        /// <summary>
        /// Update of Massflux Parameter
        /// Massflux is defined to be positive when pointing in LS normal direction, i.e.
        /// mass flows from - phase to + phase.
        /// Keep in mind maybe it is more stable to only update massflux once per timestep, not in every nonlinear iteration.
        /// </summary>
        public class MassFluxExtension_Evaporation_MixtureFraction : ParameterS, ILevelSetParameter {
            private XNSFE_OperatorConfiguration config;
            private DualLevelSet levelSet;
            private double time;
            private double Q;
            private double zSt;
            public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension };

            public override DelParameterFactory Factory => ParameterFactory;

            public override DelPartialParameterUpdate Update {
                get {
                    //return MassFluxExtension_Evaporation_Update;  // seems more stable to update once per timestep
                    return null;
                }
            }

            public MassFluxExtension_Evaporation_MixtureFraction(XNSFE_OperatorConfiguration config, double Q, double zST) {
                this.config = config;
                this.Q = Q;
                this.zSt = zST;
            }

            public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
                var massfluxext = new (string, DGField)[1];
                string paramName = BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension;
                Basis b = new Basis(DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
                DGField MassFluxExtension = new SinglePhaseField(b, paramName);
                massfluxext[0] = (paramName, MassFluxExtension);
                return massfluxext;
            }

            public void MassFluxExtension_Evaporation_Update(double phystime, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
                MassFluxExtension_Evaporation_Update(DomainVarFields, ParameterVarFields);
            }

            public void MassFluxExtension_Evaporation_Update(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
                var thermalParams = config.getThermParams;
                double kA = 0, kB = 0;
                foreach (var spc in levelSet.Tracker.SpeciesNames) {
                    switch (spc) {
                        case "A": { kA = thermalParams.k_A; break; }
                        case "B": { kB = thermalParams.k_B; break; }
                        case "C": { break; }
                        default: { throw new ArgumentException("unknown species"); }
                    }
                }
                var paramName = BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension;

                SinglePhaseField MassFluxField = new SinglePhaseField(ParameterVarFields[paramName].Basis);
                int order = MassFluxField.Basis.Degree * MassFluxField.Basis.Degree + 2;

                XDGField MixtureFraction = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.MixtureFraction];

                int D = levelSet.Tracker.GridDat.SpatialDimension;
                MassFluxField.Clear();
                MassFluxField.ProjectField(1.0,
                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                        int K = result.GetLength(1); // No nof Nodes

                        MultidimensionalArray GradMixtFrac_A_Res = MultidimensionalArray.Create(Len, K, D);
                        MultidimensionalArray GradMixtFrac_B_Res = MultidimensionalArray.Create(Len, K, D);

                        MixtureFraction.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradMixtFrac_A_Res);
                        MixtureFraction.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradMixtFrac_B_Res);

                        var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);
                        double fact = Q * zSt / (1 - zSt);// TODO
                        for (int j = 0; j < Len; j++) {
                            MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                            levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);

                            for (int k = 0; k < K; k++) {
                                double qEvap = 0.0;
                                //macro region
                                for (int dd = 0; dd < D; dd++) {
                                    qEvap += ((-kA) * fact * GradMixtFrac_A_Res[j, k, dd] - (-kB) * GradMixtFrac_B_Res[j, k, dd]) * Normals[j, k, dd];
                                }

                                //Console.WriteLine("qEvap delUpdateLevelSet = {0}", qEvap);
                                double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };
                                double mEvap = (config.prescribedMassflux != null) ? config.prescribedMassflux(globX, time) : qEvap / thermalParams.hVap; // mass flux

                                //Console.WriteLine("mEvap - delUpdateLevelSet = {0}", mEvap);
                                result[j, k] = mEvap;
                            }
                        }
                    }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));

                // no extension
                ParameterVarFields[paramName].Clear();
                ParameterVarFields[paramName].Acc(1.0, MassFluxField);

                ParameterVarFields[paramName].CheckForNanOrInf(true, true, true);
            }

            public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
                this.levelSet = levelSet;
                this.time = time;

                MassFluxExtension_Evaporation_Update(DomainVarFields, ParameterVarFields);
            }
        }
    }
}