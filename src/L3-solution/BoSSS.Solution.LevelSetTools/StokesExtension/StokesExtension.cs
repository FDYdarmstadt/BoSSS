﻿using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
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
    public class StokesExtension {

        /// <summary>
        /// ctor
        /// </summary>
        public StokesExtension(int D, IncompressibleBoundaryCondMap map, int cutCellQuadOrder, double AgglomerationThrshold) {
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

        /// <summary>
        /// Bulk of the artificial Stokes extension equation;
        /// 
        /// Note that it would also be possible to integrate this operator with <see cref="GetInterfaceOperator(int)"/>:
        /// We could simply use an <see cref="XSpatialOperatorMk2"/> **on single phase fields** for the bulk matrix assembly;
        /// Then, in any cut background cell, the contribution of all cut cells would be summed up (because
        /// the single-phase field has only a single set of DOFs), in sum mimicking a standard cell integration.
        /// Obviously, this is overkill/unnecessary expense.
        /// It could be resolved by tinkering with the quadrature scheme providers 
        /// (<see cref="XSpatialOperatorMk2.VolumeQuadraturSchemeProvider"/> and <see cref="XSpatialOperatorMk2.EdgeQuadraturSchemeProvider"/>).
        /// However, since the X-Navier-Stokes terms are implemented in the XNSECommon-package, 
        /// **we cannot use them here anyway, because we would run into circular reference**.
        /// </summary>
        SpatialOperator GetBulkOperator() {
            
            SpatialOperator Op = new SpatialOperator();
            Op.QuadOrderFunction = QuadOrderFunc.Linear();

           
            {
                // Momentum, Viscous:
                for(int d = 0; d < D; d++) {
                    var visc = new ExtensionSIP(penalty_safety, d, D, map, ViscosityOption.ConstantViscosity, constantViscosityValue: viscosity);
                    Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(visc);
                }
                // Momentum, Pressure gradient:
                for(int d = 0; d < D; d++) {
                    var PresDeriv = new ExtensionPressureGradient(d, map);
                    Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(PresDeriv);
                }
                
                // Continuity:
                for(int d = 0; d < D; d++) {
                    var divVol = new Divergence_DerivativeSource(d, D);
                    var divEdg = new ExtensionDivergenceFlux(d, map);
                    Op.EquationComponents[EquationNames.ContinuityEquation].Add(divVol);
                    Op.EquationComponents[EquationNames.ContinuityEquation].Add(divEdg);
                }
            }


            Op.OperatorCoefficientsProvider = BulkOperatorCoefficientsProvider;

            Op.TemporalOperator = null;
            Op.IsLinear = true;
            Op.Commit();
            return Op;
        }

        CoefficientSet BulkOperatorCoefficientsProvider(Foundation.Grid.IGridData g, double time) {

            var r = new CoefficientSet() {
                GrdDat = g
            };

            if(g is Foundation.Grid.Classic.GridData cgdat) {
                r.CellLengthScales = LatestLengthScales;
                r.EdgeLengthScales = null;
            } else {
                Console.Error.WriteLine("Rem: still missing cell length scales for grid type " + g.GetType().FullName);
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
        XSpatialOperatorMk2 GetInterfaceOperator(int levelSetIndex, LevelSetTracker LsTrk, DGField[] InterfaceVelocity) {
            var BulkOp = GetBulkOperator();

            IEnumerable<string> requiredSpecies = LsTrk.GetSpeciesSeparatedByLevSet(levelSetIndex);
            var Op = new XSpatialOperatorMk2(
                BulkOp.DomainVar, 
                VariableNames.AsLevelSetVariable("Interface", VariableNames.VelocityVector(D)),
                BulkOp.CodomainVar,
                (int[] a, int[] b, int[] c) => m_CutCellQuadOrder,
                requiredSpecies);
            
            for(int d = 0; d < D; d++) {
                foreach ((string, string) speciesPair in LsTrk.GetSpeciesPairsSeparatedByLevSet(levelSetIndex)) {
                    string negativeSpecies = speciesPair.Item1;
                    string positiveSpecies = speciesPair.Item2;
                    Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(
                        new InteriorVelocityBoundary(positiveSpecies, negativeSpecies, levelSetIndex, d, D, InterfaceVelocity[d])
                        );
                }
            }

            Op.AgglomerationThreshold = 0.0;
            Op.TemporalOperator = null;
            Op.IsLinear = true;
            Op.FreeMeanValue[Op.DomainVar.Last()] = false;
            Op.Commit();

            return Op;
        }

        (BlockMsrMatrix OpMtx, double[] RHS) ComputeMatrix(int levelSetIndex, LevelSetTracker lsTrk, UnsetteledCoordinateMapping mapping, DGField[] Velocity) {
            BlockMsrMatrix opmtx = new BlockMsrMatrix(mapping, mapping);
            double[] RHS = new double[mapping.LocalLength];

            {
                var BulkOp = GetBulkOperator();
                BulkOp.GetMatrixBuilder(mapping, null, mapping).ComputeMatrix(opmtx, RHS);
            }

            {
                var IntfOp = GetInterfaceOperator(levelSetIndex, lsTrk, Velocity);
                var builder = IntfOp.GetMatrixBuilder(lsTrk, mapping, Velocity, mapping);

                foreach(var s in m_LatestAgglom.SpeciesList)
                    builder.CellLengthScales[s] = m_LatestAgglom.CellLengthScales[s];

               
                builder.ComputeMatrix(opmtx, RHS);
            }

    
            RHS.ScaleV(-1.0); // change from affine (lhs) to RHS
            return (opmtx, RHS);
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
                if(m_LatestLengthScales == null) {
                    var cgdat = m_LatestAgglom.Tracker.GridDat;

                    MultidimensionalArray[] allLenghtScales = new[] { cgdat.Cells.CellLengthScale };
                    foreach(var s in m_LatestAgglom.SpeciesList) {
                        m_LatestAgglom.CellLengthScales[s].AddToArray(ref allLenghtScales);
                    }

                    MultidimensionalArray lenScale = allLenghtScales[0].CloneAs();
                    int JE = lenScale.GetLength(0);
                    for(int k = 1; k < allLenghtScales.Length; k++) {
                        var len_k = allLenghtScales[k];

                        for(int j = 0; j < JE; j++) {
                            double len_kj = len_k[j];
                            if(!len_kj.IsNaNorInf()) {
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
        public void SolveExtension(int levelSetIndex, LevelSetTracker lsTrk, DGField[] VelocityAtInterface, SinglePhaseField[] ExtensionVelocity) {
            var gDat = lsTrk.GridDat;
            int deg = ExtensionVelocity[0].Basis.Degree;

            
            m_LatestAgglom = lsTrk.GetAgglomerator(lsTrk.SpeciesIdS.ToArray(), this.m_CutCellQuadOrder, this.AgglomerationThreshold);

            SinglePhaseField dummyPressure = new SinglePhaseField(new Basis(gDat, deg - 1), "DummyPressure");

            CoordinateVector ExtenstionSolVec = new CoordinateVector(ExtensionVelocity.Cat(dummyPressure));

            (BlockMsrMatrix OpMtx, double[] RHS) = ComputeMatrix(levelSetIndex, lsTrk, ExtenstionSolVec.Mapping, VelocityAtInterface);

            // should be replaced by something more sophisticated
            OpMtx.Solve_Direct(ExtenstionSolVec, RHS);
        }
    }
}
