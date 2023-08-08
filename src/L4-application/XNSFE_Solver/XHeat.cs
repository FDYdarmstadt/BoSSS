using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XNSFE_Solver {
    public class XHeat : SolverWithLevelSetUpdater<XNSFE_Control> {
        ThermalMultiphaseBoundaryCondMap boundaryMap;

        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {

            int D = this.GridData.SpatialDimension;

            int pTemp = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Temperature].Degree;
            // configuration for Temperature
            var confTemp = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pTemp }, //Math.Max(1, pTemp - iLevel) },
                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Temperature) }
            };
            configsLevel.Add(confTemp);

            // configuration for auxiliary heat flux
            if(this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                int pFlux;
                if(this.Control.FieldOptions.TryGetValue("HeatFlux*", out FieldOpts f)) {
                    pFlux = f.Degree;
                } else if(this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.HeatFluxX, out FieldOpts f1)) {
                    pFlux = f1.Degree;
                } else {
                    throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of HeatFlux not found");
                }
                for(int d = 0; d < D; d++) {
                    var confHeatFlux = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pFlux }, // Math.Max(1, pFlux - iLevel) },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.HeatFluxVectorComponent(d)) }
                    };
                    configsLevel.Add(confHeatFlux);
                }
            }
        }

        public override int QuadOrder() {
            //QuadOrder
            int degT;
            if(Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.Temperature, out FieldOpts field)) {
                degT = field.Degree;
            } else {
                throw new Exception("Temperature not found!");
            }
            int quadOrder = degT * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if(this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                quadOrder *= 2;
                quadOrder += 1;
            }
            return quadOrder;
        }

        int VelocityDegree() {
            int pVel;
            if(this.Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v)) {
                pVel = v.Degree;
            } else if(this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1)) {
                pVel = v1.Degree;
            } else if(this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.Temperature, out FieldOpts t1)) {
                Console.WriteLine("Degree of Velocity not found, using Temperature Degree");
                pVel = t1.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pVel;
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {
            OperatorFactory opFactory = new OperatorFactory();
            boundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;
            XOP.IsLinear = true;
            XOP.Commit();

            return XOP;
        }

        public void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);

            int quadOrder = QuadOrder();
            // add Heat equation components
            // ============================
            opFactory.AddEquation(new Heat("A", D, boundaryMap, config));
            opFactory.AddEquation(new Heat("B", D, boundaryMap, config));

            if(config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                for(int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new HeatFlux("A", d, D, boundaryMap, config));
                    opFactory.AddEquation(new HeatFlux("B", d, D, boundaryMap, config));
                    opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, boundaryMap, config));
                }
            }
            opFactory.AddEquation(new HeatInterface("A", "B", D, boundaryMap, config));
            opFactory.AddCoefficient(new SlipLengths(config));
            opFactory.AddCoefficient(new EvapMicroRegion());


            // Add Velocity parameters as prescribed variables (needed for level set)
            for (int d = 0; d < D; d++)
                opFactory.AddParameter(Velocity0Prescribed.CreateFrom(lsUpdater.Tracker, d, D, Control));

            if (this.Control.ThermalParameters.IncludeConvection) {
                Velocity0MeanPrescribed v0Mean = new Velocity0MeanPrescribed(D, lsUpdater.Tracker, quadOrder);
                opFactory.AddParameter(v0Mean);
                lsUpdater.AddLevelSetParameter("Phi", v0Mean);
            }

            // Level set evolver need the gradient
            BeltramiGradient lsGradient = FromControl.BeltramiGradient(Control, "Phi", D);
            lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, lsGradient);            
        }        

        protected override IncompressibleBoundaryCondMap GetBcMap() {
            // required only for the stokes extension
            throw new NotImplementedException();
        }

        protected override int NoOfLevelSets {
            get {
                return 1;
            }
        }

        protected override Array SpeciesTable {
            get {
                return new[] { "A", "B" };
            }
        }

        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            return new LevelSetVelocity(VariableNames.LevelSetCG, GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters); 
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetTimestep();
            Console.WriteLine($"Starting timesetp {TimestepNo}, dt={dt}");
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            Console.WriteLine($"done with timestep {TimestepNo}");
            return dt;
        }
    }
}
