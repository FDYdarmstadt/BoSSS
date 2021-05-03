using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.NSECommon.Operator.Viscosity;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {


    /// <summary>
    /// Level-set evolution using a div-free velocity extension computed from Stokes equation.
    /// Due to the div-free velocity extension, the level-set evolution can be implemented as scalar convection,
    /// which has proven to be very stable in DG.
    /// </summary>
    public class XStokesExtension {

        /// <summary>
        /// ctor
        /// </summary>
        public XStokesExtension(int D, IncompressibleBoundaryCondMap map, int cutCellQuadOrder, double AgglomerationThrshold) {
            this.D = D;
            this.map = map;
            this.m_CutCellQuadOrder = cutCellQuadOrder;
            this.AgglomerationThreshold = AgglomerationThrshold;
        }

        int D;
        int m_CutCellQuadOrder;
        IncompressibleBoundaryCondMap map;
        double AgglomerationThreshold;

        const double penalty_safety = 4.0;

        const double viscosity = 1.0;

        CoefficientSet BulkOperatorCoefficientsProvider(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {

            var r = new CoefficientSet() {
                GrdDat = lstrk.GridDat
            };

            if (r.GrdDat is Foundation.Grid.Classic.GridData cgdat) {
                r.CellLengthScales = LatestLengthScales;
                r.EdgeLengthScales = null;
            } else {
                Console.Error.WriteLine("Rem: still missing cell length scales for grid type " + r.GrdDat.GetType().FullName);
            }


            //foreach(var kv in UserDefinedValues) {
            //    r.UserDefinedValues[kv.Key] = kv.Value;
            //}


            r.HomotopyValue = 1.0;

            return r;
        }

        /// <summary>
        /// Interface part of the Stokes extension; 
        /// this provides the coupling of the artificial Stokes equation for the extension
        /// to the physical equation, <see cref="InteriorVelocityBoundary"/>.
        /// </summary>
        XSpatialOperatorMk2 GetOperator(LevelSetTracker LsTrk, DGField[] InterfaceVelocity) {

            var Op = new XSpatialOperatorMk2(
                VariableNames.VelocityVector(D).Cat(VariableNames.Pressure),
                VariableNames.AsLevelSetVariable("Interface", VariableNames.VelocityVector(D)),
                EquationNames.MomentumEquations(D).Cat(EquationNames.ContinuityEquation),
                (int[] a, int[] b, int[] c) => m_CutCellQuadOrder,
                new[] { "A", "B" });

            // Momentum, Viscous:
            for (int d = 0; d < D; d++) {
                var visc = new SipViscosity_GradU(penalty_safety, d, D, map, ViscosityOption.ConstantViscosity, constantViscosityValue: viscosity);
                Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(visc);
            }
            // Momentum, Pressure gradient:
            for (int d = 0; d < D; d++) {
                var PresDeriv = new PressureGradientLin_d(d, map);
                Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(PresDeriv);
            }

            // Continuity:
            for (int d = 0; d < D; d++) {
                var divVol = new Divergence_DerivativeSource(d, D);
                var divEdg = new Divergence_DerivativeSource_Flux(d, map);
                Op.EquationComponents[EquationNames.ContinuityEquation].Add(divVol);
                Op.EquationComponents[EquationNames.ContinuityEquation].Add(divEdg);
            }
            //Interior Surface
            for (int d = 0; d < D; d++) {
                Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(
                    new InteriorVelocityBoundary( d, InterfaceVelocity[d])
                    );
            }
            //Immersed Boundary
            for (int d = 0; d < D; d++) {
                Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(
                    new ViscosityAtSlipIB(d, D, penalty_safety, viscosity, 1, "A", "C", false)
                    );
                Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(
                    new ViscosityAtSlipIB(d, D, penalty_safety, viscosity, 1, "B", "C", false)
                    );
                Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(
                    new BoSSS.Solution.NSECommon.Operator.Pressure.PressureFormAtIB(d, D, 1, "A", "C")
                    );
                Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(
                    new BoSSS.Solution.NSECommon.Operator.Pressure.PressureFormAtIB(d, D, 1, "B", "C")
                    );
            }
            Op.EquationComponents[EquationNames.ContinuityEquation].Add(
                new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, 1, "A", "C", false));
            Op.EquationComponents[EquationNames.ContinuityEquation].Add(
                new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, 1, "B", "C", false));

            Op.AgglomerationThreshold = 0.0;
            Op.TemporalOperator = null;
            Op.IsLinear = true;
            Op.FreeMeanValue[Op.DomainVar.Last()] = !this.map.DirichletPressureBoundary;
            Op.OperatorCoefficientsProvider = BulkOperatorCoefficientsProvider;
            Op.Commit();

            return Op;
        }

        MultiphaseCellAgglomerator m_LatestAgglom;


        MultidimensionalArray m_LatestLengthScales;

        /// <summary>
        /// Cell Length scales 
        /// </summary>
        /// <remarks>
        /// The Stokes-Extension is some kind of hybrid DG-XDG method:
        /// the main bulk equation is certainly DG but the source/coupling term -- which is enforced by penalties (!) -- is
        /// defined on cut cells, resp. cut cell boundaries.
        /// Therefore, we chose the length scales pessimistic (minimum over cut cells) and hope for the best.
        /// </remarks>
        MultidimensionalArray LatestLengthScales {
            get {
                if (m_LatestLengthScales == null) {
                    var cgdat = m_LatestAgglom.Tracker.GridDat;

                    MultidimensionalArray[] allLenghtScales = new[] { cgdat.Cells.CellLengthScale };
                    foreach (var s in m_LatestAgglom.SpeciesList) {
                        m_LatestAgglom.CellLengthScales[s].AddToArray(ref allLenghtScales);
                    }

                    MultidimensionalArray lenScale = allLenghtScales[0].CloneAs();
                    int JE = lenScale.GetLength(0);
                    for (int k = 1; k < allLenghtScales.Length; k++) {
                        var len_k = allLenghtScales[k];

                        for (int j = 0; j < JE; j++) {
                            double len_kj = len_k[j];
                            if (!len_kj.IsNaNorInf()) {
                                lenScale[j] = Math.Min(lenScale[j], len_kj);
                            }
                        }
                    }

                    m_LatestLengthScales = lenScale;
                    m_LatestLengthScales.CheckForNanOrInf();
                }
                return m_LatestLengthScales;
            }
        }

        /// <summary>
        /// actually computing the extension velocity
        /// </summary>
        /// <param name="lsTrk"></param>
        /// <param name="VelocityAtInterface">
        /// input: must contain the velocity field to extend, only relevant in cut cells.
        /// </param>
        /// <param name="ExtensionVelocity">
        /// output
        /// </param>
        public void SolveExtension(LevelSetTracker lsTrk, DGField[] VelocityAtInterface, SinglePhaseField[] ExtensionVelocity) {
            var gDat = lsTrk.GridDat;
            int deg = ExtensionVelocity[0].Basis.Degree;


            m_LatestAgglom = lsTrk.GetAgglomerator(lsTrk.SpeciesIdS.ToArray(), this.m_CutCellQuadOrder, this.AgglomerationThreshold);

            SinglePhaseField dummyPressure = new SinglePhaseField(new Basis(gDat, deg - 1), "DummyPressure");

            CoordinateVector ExtenstionSolVec = new CoordinateVector(ExtensionVelocity.Cat(dummyPressure));
            var Op = GetOperator(lsTrk, VelocityAtInterface);

            UniSolver.NonXDG_LsTrk = lsTrk; //A so called freak case!
            Op.Solve(ExtenstionSolVec.Mapping);
            UniSolver.NonXDG_LsTrk = null;
        }
    }
}

