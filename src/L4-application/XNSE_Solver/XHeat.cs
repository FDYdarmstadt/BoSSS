using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.SpatialOperator;

namespace BoSSS.Application.XNSE_Solver {
    class XHeat : XCommon<XNSE_Control> {

        public override void MultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {

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

        protected override int QuadOrder() {
            //QuadOrder
            int degT;
            if (Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.Temperature, out FieldOpts field)) {
                degT = field.Degree;
            } else {
                throw new Exception("Temperature not found!");
            }
            int quadOrder = degT * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if (this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                quadOrder *= 2;
                quadOrder += 1;
            }
            return quadOrder;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {
            Tecplot.PlotFields(this.m_RegisteredFields, "XHEAT_Solver" + timestepNo, physTime, superSampling);
            if (Timestepping?.Parameters != null) {
                Tecplot.PlotFields(Timestepping.Parameters, "XHEAT_Solver_Params" + timestepNo, physTime, superSampling);
            }
        }

        public override void SetOperatorEquations(int D, OperatorFactory opFactory) {

            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            ThermalMultiphaseBoundaryCondMap boundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            XHeatOperatorProvider.SetOperatorEquations(D, opFactory, QuadOrder(), boundaryMap, this.lsUpdater, this.Control, config);
        }

        public override void SetOperatorParameter(int D, OperatorFactory opFactory) {

            opFactory.SetCoefficient(XHeatCoefficients);

            XHeatOperatorProvider.SetOperatorParameter(D, opFactory, QuadOrder(), null, this.lsUpdater, this.Control, null);

        }
        private CoefficientSet XHeatCoefficients(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time){

            var r = new CoefficientSet() {
                GrdDat = lstrk.GridDat
            };
            var g = lstrk.GridDat;
            if (g is Foundation.Grid.Classic.GridData cgdat) {
                r.CellLengthScales = cgdat.Cells.CellLengthScale;
                r.EdgeLengthScales = cgdat.Edges.h_min_Edge;

            } else {
                Console.Error.WriteLine("Rem: still missing cell length scales for grid type " + g.GetType().FullName);
            }

            //Console.WriteLine("Heat configured without Evaporation and no fixed Interface Temperature");
            BitArray EvapMicroRegion = lstrk.GridDat.GetBoundaryCells().GetBitMask();
            EvapMicroRegion.SetAll(true);
            r.UserDefinedValues["EvapMicroRegion"] = EvapMicroRegion;

            return r;
        }

        public override void SetSpatialOperator(out XSpatialOperatorMk2 XOP, int D, OperatorFactory opFactory) {
            XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            XOP.LinearizationHint = LinearizationHint.AdHoc;
            XOP.Commit();
        }        

    }
}
