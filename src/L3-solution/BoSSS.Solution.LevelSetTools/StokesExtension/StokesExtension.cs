using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
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
        public StokesExtension(int D, IncompressibleBoundaryCondMap map, int cutCellQuadOrder, double AgglomerationThrshold, bool fullStokes, StokesExtentionBoundaryOption useBCMap = StokesExtentionBoundaryOption.FreeSlipAtWall) {
            this.D = D;
            this.map = map;
            this.m_CutCellQuadOrder = cutCellQuadOrder;
            this.AgglomerationThreshold = AgglomerationThrshold;
            this.fullStokes = fullStokes;
            this.m_useBCMap = useBCMap; // if false, Neumann boundary conditions are applied everywhere. If true, the boundary conditions from map are applied.
        }

        int D;
        int m_CutCellQuadOrder;
        IncompressibleBoundaryCondMap map;
        double AgglomerationThreshold;
        bool fullStokes;
        StokesExtentionBoundaryOption m_useBCMap;

        const double penalty_safety = 4.0;

        const double viscosity = 1000.0;

        /// <summary>
        /// If set, the extension's linear system is solved through PETSc
        /// (<see cref="BoSSS.Solution.AdvancedSolvers.PetscGMRESConfig"/>) using this options string,
        /// instead of the default direct-solver fallback chain (<see cref="ilPSP.LinSolvers.SimpleSolversInterface.Solve_Direct"/>).
        /// If the PETSc solve throws, the existing direct-solver fallback in <see cref="SolveExtension"/> still applies.
        /// Null/empty (the default) preserves the original behavior.
        /// </summary>
        public string PetscOptions = null;

        /// <summary>
        /// If true, <see cref="SolveExtension"/> uses whatever is already in
        /// <c>ExtensionVelocity</c> on entry (typically the previous time step's converged
        /// extension velocity, since callers no longer clear it beforehand) as the initial guess
        /// for the iterative solve, transformed into the solver's equilibrated basis via
        /// <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.TransformSolInto{T1, T2}"/>.
        /// If false (the default), the solve always starts from zero, matching the original behavior.
        /// </summary>
        public bool WarmStart = false;

        /// <summary>
        /// Only relevant if <see cref="fullStokes"/> is set (via the constructor). If true, the
        /// velocity and pressure <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig"/>
        /// are combined into a single block using
        /// <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.SchurComplement"/>, mirroring
        /// <c>XNSE_Control.UseSchurBlockPrec</c> (see <c>XNSE.cs</c>'s <c>AddMultigridConfigLevel</c>).
        /// If false (the default), velocity and pressure keep their separate
        /// <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite"/>/
        /// <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.IdMass_DropIndefinite"/> configs, matching
        /// the original behavior. Best paired with a field-split-aware <see cref="PetscOptions"/> string; combining it
        /// with the internal <see cref="BoSSS.Solution.AdvancedSolvers.SchurComplementForStokes"/> fallback solver
        /// applies Schur-complement structure twice.
        /// </summary>
        public bool UseSchurBlockPrec = false;

        /// <summary>
        /// If true, skips BoSSS's own local (per-cell) preconditioning transform entirely -- neither
        /// <see cref="UseSchurBlockPrec"/>'s combined block nor the "standard"
        /// <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite"/>/
        /// <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.IdMass_DropIndefinite"/> configs are
        /// used -- instead every field gets <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.Eye"/>
        /// (identity: no change of basis at all), handing PETSc/the fallback solver the raw <c>OpMtx</c>
        /// unmodified. Takes precedence over <see cref="UseSchurBlockPrec"/> if both are set. Note: without
        /// <see cref="UsePseudoTimestepping"/> also regularizing the raw matrix, the pressure-pressure
        /// diagonal block stays structurally exactly zero here (see remarks on <see cref="UseSchurBlockPrec"/>),
        /// which a plain ILU factorization cannot handle. False (the default) preserves the original behavior.
        /// </summary>
        public bool DisableLocalPreconditioning = false;

        /// <summary>
        /// Only relevant if <see cref="PetscOptions"/> is set: forwarded to
        /// <see cref="BoSSS.Solution.AdvancedSolvers.PetscGMRESConfig.UseGeometricMultigrid"/>.
        /// False (the default) preserves the original flat-matrix behavior.
        /// </summary>
        public bool UseGeometricMultigrid = false;

        /// <summary>
        /// Only relevant if <see cref="UseGeometricMultigrid"/> is true: forwarded to
        /// <see cref="BoSSS.Solution.AdvancedSolvers.PetscGMRESConfig.CoarseKickIn"/>.
        /// </summary>
        public int CoarseKickIn = 90000;

        /// <summary>
        /// If true, <see cref="SolveExtension"/> solves via pseudo-transient continuation instead of a
        /// single solve: repeatedly solves <c>(OpMtx + I/dtau)*u = RHS + u/dtau</c> for a growing
        /// <c>dtau</c> (switched-evolution-relaxation between <see cref="PseudoDtau0"/> and
        /// <see cref="PseudoDtauMax"/>), until the TRUE (un-regularized) residual <c>RHS - OpMtx*u</c>,
        /// relative to <c>‖RHS‖</c>, falls below <see cref="PseudoTimeConvergenceCriterion"/> or
        /// <see cref="MaxPseudoTimeSteps"/> is reached. False (the default) preserves the original
        /// single-solve behavior.
        /// </summary>
        public bool UsePseudoTimestepping = false;

        /// <summary>
        /// Initial pseudo-time step size. Only relevant if <see cref="UsePseudoTimestepping"/> is set.
        /// </summary>
        public double PseudoDtau0 = 1.0;

        /// <summary>
        /// Upper bound for the pseudo-time step size's growth (switched-evolution-relaxation).
        /// Only relevant if <see cref="UsePseudoTimestepping"/> is set.
        /// </summary>
        public double PseudoDtauMax = 1.0e8;

        /// <summary>
        /// Maximum number of pseudo-time steps. Only relevant if <see cref="UsePseudoTimestepping"/> is set.
        /// </summary>
        public int MaxPseudoTimeSteps = 20;

        /// <summary>
        /// Convergence criterion on the true (un-regularized) residual norm, RELATIVE to <c>‖RHS‖</c>
        /// (i.e. stops once <c>‖RHS - OpMtx*u‖ / ‖RHS‖</c> falls below this value) -- an absolute
        /// criterion would be meaningless here, since ‖RHS‖ itself can vary by many orders of magnitude
        /// depending on the problem's physical scale. Only relevant if <see cref="UsePseudoTimestepping"/>
        /// is set.
        /// </summary>
        public double PseudoTimeConvergenceCriterion = 1.0e-10;

        /// <summary>
        /// Bulk of the artificial Stokes extension equation;
        /// 
        /// Note that it would also be possible to integrate this operator with <see cref="GetInterfaceOperator(int)"/>:
        /// We could simply use an <see cref="XDifferentialOperatorMk2"/> **on single phase fields** for the bulk matrix assembly;
        /// Then, in any cut background cell, the contribution of all cut cells would be summed up (because
        /// the single-phase field has only a single set of DOFs), in sum mimicking a standard cell integration.
        /// Obviously, this is overkill/unnecessary expense.
        /// It could be resolved by tinkering with the quadrature scheme providers 
        /// (<see cref="XDifferentialOperatorMk2.VolumeQuadraturSchemeProvider"/> and <see cref="XDifferentialOperatorMk2.EdgeQuadraturSchemeProvider"/>).
        /// However, since the X-Navier-Stokes terms are implemented in the XNSECommon-package, 
        /// **we cannot use them here anyway, because we would run into circular reference**.
        /// </summary>
        DifferentialOperator GetBulkOperator() {
            var DomVar = VariableNames.VelocityVector(D);
            var CoDomVar = EquationNames.MomentumEquations(D);
            if (fullStokes) {
                DomVar = DomVar.Cat(VariableNames.Pressure);
                CoDomVar = CoDomVar.Cat(EquationNames.ContinuityEquation);
            }
            DifferentialOperator Op = new DifferentialOperator(DomVar, CoDomVar, QuadOrderFunc.Linear() );
            //Op.QuadOrderFunction = QuadOrderFunc.Linear();

            {
                // Momentum, Viscous:
                for(int d = 0; d < D; d++) {
                    var visc = new ExtensionSIP(penalty_safety, d, D, map, ViscosityOption.ConstantViscosity, constantViscosityValue: viscosity, useBCMap: m_useBCMap);
                    Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(visc);
                }

                if (fullStokes) {
                    // Momentum, Pressure gradient:
                    for (int d = 0; d < D; d++) {
                        var PresDeriv = new ExtensionPressureGradient(d, map, useBCMap: m_useBCMap);
                        Op.EquationComponents[EquationNames.MomentumEquationComponent(d)].Add(PresDeriv);
                    }

                    // Continuity:
                    for (int d = 0; d < D; d++) {
                        var divVol = new Divergence_DerivativeSource(d, D);
                        var divEdg = new ExtensionDivergenceFlux(d, map, useBCMap: m_useBCMap);
                        Op.EquationComponents[EquationNames.ContinuityEquation].Add(divVol);
                        Op.EquationComponents[EquationNames.ContinuityEquation].Add(divEdg);
                    }
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

            MultidimensionalArray SlipLengths = MultidimensionalArray.Create(g.iGeomCells.NoOfLocalUpdatedCells);
            SlipLengths.AccConstant(-1.0); // freeslip on all slipboundaries
            r.UserDefinedValues["SlipLengths"] = SlipLengths;

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
        XDifferentialOperatorMk2 GetInterfaceOperator(int levelSetIndex, LevelSetTracker LsTrk, DGField[] InterfaceVelocity) {
            var BulkOp = GetBulkOperator();

            IEnumerable<string> requiredSpecies = LsTrk.GetSpeciesSeparatedByLevSet(levelSetIndex);
            var Op = new XDifferentialOperatorMk2(
                BulkOp.DomainVar, 
                VariableNames.AsLevelSetVariable("Interface", VariableNames.VelocityVector(D)),
                BulkOp.CodomainVar,
                (int[] a, int[] b, int[] c) => m_CutCellQuadOrder,
                requiredSpecies);

            //var cc = LsTrk.Regions.GetCutCellMask();
            //Console.WriteLine("Stokes Extension No of cut cells " + cc.NoOfItemsLocally.MPISum());
            //Console.WriteLine(InterfaceVelocity.Select(vel => vel.L2Norm(cc)).ToConcatString("[", ",", "]"));
            //Console.WriteLine(InterfaceVelocity.Select(vel => vel.L2Norm()).ToConcatString("[", ",", "]"));
            
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

        static int timestepNo = 0;
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
            using(FuncTrace tr = new FuncTrace()) {
                tr.InfoToConsole = false;
                var gDat = lsTrk.GridDat;
                int deg = ExtensionVelocity[0].Basis.Degree;

                double[] InputVelL2 = VelocityAtInterface.Select(vel => vel.L2Norm()).ToArray();


                m_LatestAgglom = lsTrk.GetAgglomerator(lsTrk.SpeciesIdS.ToArray(), this.m_CutCellQuadOrder, this.AgglomerationThreshold);

                int J = gDat.iLogicalCells.NoOfLocalUpdatedCells;
                var Gradients = lsTrk.DataHistories[0].Current.GetLevelSetGradients(gDat.Grid.RefElements[0].Center, 0, J);
                var yValues = Gradients.ExtractSubArrayShallow(-1, -1, 1);

                DGField[] DummySolFields = ExtensionVelocity;
                if (fullStokes) {
                    SinglePhaseField dummyPressure = new SinglePhaseField(new Basis(gDat, deg - 1), "DummyPressure");
                    DummySolFields = DummySolFields.Cat(dummyPressure);
                }
                CoordinateVector ExtenstionSolVec = new CoordinateVector(DummySolFields);

                (BlockMsrMatrix OpMtx, double[] RHS) = ComputeMatrix(levelSetIndex, lsTrk, ExtenstionSolVec.Mapping, VelocityAtInterface);

                // Level-set-restricted (surface) quadrature for the interface-velocity-mismatch
                // diagnostic below: a plain CellMask-based L2Error would integrate over the FULL
                // volume of each cut cell, not just the interface itself. GetLevelSetQuadScheme
                // gives the same kind of surface quadrature that InteriorVelocityBoundary (an
                // ILevelSetForm) is actually integrated with when the interface condition is
                // assembled into OpMtx/RHS, so this measures compliance with exactly what's enforced.
                var ifaceSchemeHelper = lsTrk.GetXDGSpaceMetrics(lsTrk.SpeciesIdS.ToArray(), this.m_CutCellQuadOrder, 1).XQuadSchemeHelper;
                var ifaceQuadRule = ifaceSchemeHelper.GetLevelSetQuadScheme(levelSetIndex, lsTrk.Regions.GetCutCellMask4LevSet(levelSetIndex))
                    .Compile(gDat, this.m_CutCellQuadOrder);

                double InterfaceVelocityMismatchNorm() {
                    double acc = 0.0;
                    for(int d = 0; d < D; d++) {
                        double err_d = ExtensionVelocity[d].LxError(VelocityAtInterface[d].Evaluate, null, ifaceQuadRule).Sqrt();
                        acc += err_d * err_d;
                    }
                    return Math.Sqrt(acc);
                }

                //var RHSdg = new CoordinateVector(DummySolFields.Select(dd => dd.CloneAs()));
                //RHSdg.SetV(RHS, 1.0);
                //double[] RHSL2norm = RHSdg.Fields.Select(f => f.L2Norm()).ToArray();

                if (this.UseSchurBlockPrec && !this.fullStokes)
                    throw new InvalidOperationException($"{nameof(UseSchurBlockPrec)} requires {nameof(fullStokes)} (no pressure to combine with otherwise).");

                var AggGrid = CoarseningAlgorithms.CreateSequence(ExtenstionSolVec.Mapping.GridDat, -1);
                var AggBasis = AggregationGridBasis.CreateSequence(AggGrid, ExtenstionSolVec.Mapping.BasisS);
                var configsLevel = new List<BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig>();
                if (this.DisableLocalPreconditioning) {
                    // no local preconditioning transform at all -- every field gets Mode.Eye (identity),
                    // handing PETSc/the fallback solver the raw OpMtx unmodified.
                    for (int i = 0; i < ExtenstionSolVec.Mapping.Fields.Count; i++) {
                        configsLevel.Add(new BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig() {
                            mode = BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.Eye,
                            VarIndex = new int[] { i },
                            DegreeS = new int[] { ExtenstionSolVec.Mapping.Fields[i].Basis.Degree }
                        });
                    }
                } else if (this.fullStokes && this.UseSchurBlockPrec) {
                    // single combined block for velocity + pressure, mirroring XNSE.cs's UseSchurBlockPrec
                    configsLevel.Add(new BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig() {
                        mode = BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.SchurComplement,
                        VarIndex = Enumerable.Range(0, ExtenstionSolVec.Mapping.Fields.Count).ToArray(),
                        DegreeS = ExtenstionSolVec.Mapping.Fields.Select(f => f.Basis.Degree).ToArray()
                    });
                } else if (this.fullStokes) {
                    for (int i = 0; i < ExtenstionSolVec.Mapping.Fields.Count; i++) {
                        bool isPressure = i == ExtenstionSolVec.Mapping.Fields.Count - 1;
                        configsLevel.Add(new BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig() {
                            mode = isPressure
                                ? BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.IdMass_DropIndefinite
                                : BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                            VarIndex = new int[] { i },
                            DegreeS = new int[] { ExtenstionSolVec.Mapping.Fields[i].Basis.Degree }
                        });
                    }
                } else {
                    for (int i = 0; i < ExtenstionSolVec.Mapping.Fields.Count; i++) {
                        configsLevel.Add(new BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig() {
                            mode = BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                            VarIndex = new int[] { i },
                            DegreeS = new int[] { ExtenstionSolVec.Mapping.Fields[i].Basis.Degree }
                        });
                    }
                }
                var MultigridOperatorConfig = new BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig[1][] { configsLevel.ToArray() };

                var bkup_RHS = RHS.CloneAs();
                Exception e;
                try {
                    // Solves one (possibly pseudo-time-regularized) linearization, dispatching to PETSc
                    // or the internal Schur-complement fallback exactly as before.
                    void SolveOnce(MultigridOperator mgOp_, double[] Pre_Sol_, double[] Pre_RHS_) {
                        if(!string.IsNullOrWhiteSpace(PetscOptions)) {
                            var petscSolver = new PetscGMRESConfig() {
                                PetscOptions = PetscOptions,
                                UseGeometricMultigrid = UseGeometricMultigrid,
                                CoarseKickIn = CoarseKickIn
                            }.CreateInstance(mgOp_);
                            petscSolver.Solve(Pre_Sol_, Pre_RHS_);
                            petscSolver.Dispose();
                        } else {
                            // HARDCODED FOR TESTING: GMRESChebyshevMGConfig instead of the usual
                            // SchurComplementForStokes fallback -- a plain outer GMRES wrapping a
                            // non-accelerating classic V-/W-cycle multigrid (NOT
                            // OrthonormalizationMultigrid's self-accelerating K-cycle), which is
                            // the literal structural analog of the PETSc setup it mirrors as
                            // closely as possible:
                            //   "-ksp_type gmres -ksp_atol 1e-6 -ksp_rtol 1e-10
                            //    -ksp_gmres_restart 500 -ksp_max_it 1000
                            //    -pc_type mg -pc_mg_cycle_type w -mg_levels_ksp_type chebyshev
                            //    -mg_levels_ksp_max_it 5 -mg_levels_pc_type jacobi
                            //    -mg_coarse_ksp_type preonly -mg_coarse_pc_type lu
                            //    -mg_coarse_pc_factor_mat_solver_type mumps -mat_mumps_icntl_24 1"
                            //
                            // Known gaps (no BoSSS-native equivalent currently exposed):
                            //  - "-mat_mumps_icntl_24 1" (null-pivot detection): ilPSP's
                            //    MUMPSSolver wrapper doesn't expose MUMPS's ICNTL array at all.
                            //  - UsePointJacobi = true is a LITERAL match to "-mg_levels_pc_type
                            //    jacobi", but this is a saddle-point system (velocity+pressure)
                            //    with exactly-zero-diagonal continuity rows -- ChebyshevSmoother.Init
                            //    will throw immediately if it hits one, the same failure mode that
                            //    made plain SOR/Jacobi unusable for XNSFE earlier. If this throws,
                            //    switch back to UsePointJacobi = false (BlockJacobi).
                            var mgSolver = new GMRESChebyshevMGConfig() {
                                GMRESRestart = 500,                                    // -ksp_gmres_restart 500
                                MaxSolverIterations = 1000,                            // -ksp_max_it 1000
                                ConvergenceCriterion = 1e-10,                          // -ksp_rtol 1e-10
                                AbsoluteTolerance = 1e-6,                              // -ksp_atol 1e-6
                                TargetBlockSize = 10000,
                                UseWCycle = true,                                      // -pc_mg_cycle_type w
                                ChebyshevDegree = 5,                                   // -mg_levels_ksp_max_it 5
                                UsePointJacobi = true,                                 // -mg_levels_pc_type jacobi
                                CoarseDirectSolver = DirectSolver._whichSolver.MUMPS,  // -mg_coarse_pc_factor_mat_solver_type mumps
                            }.CreateInstance(mgOp_);

                            mgSolver.Solve(Pre_Sol_, Pre_RHS_);
                            mgSolver.Dispose();

                            //var schurSolver = new SchurComplementForStokes.Config() {
                            //    ConvergenceCriterion = 1e-10,
                            //    MaxSolverIterations = 50,
                            //    WriteDebugInfoLevel = 0
                            //}.CreateInstance(mgOp_);
                            //schurSolver.Solve(Pre_Sol_, Pre_RHS_);
                            //schurSolver.Dispose();
                        }
                    }

                    if(UsePseudoTimestepping) {
                        // Pseudo-transient continuation: repeatedly solve (OpMtx + I/dtau)*u = RHS + u/dtau
                        // for a growing dtau (switched-evolution-relaxation), until the TRUE (un-regularized)
                        // residual RHS - OpMtx*u is small, or MaxPseudoTimeSteps is reached. dtau/I are added
                        // to the whole system (velocity and pressure alike), not just the momentum rows.
                        tr.Info("StokesExtension.SolveExtension trying with pseudo-timestepping");
                        double dtau = PseudoDtau0;
                        double[] u = WarmStart ? ExtenstionSolVec.ToArray() : new double[ExtenstionSolVec.Mapping.LocalLength];
                        double resNorm_prev = RHS.MPI_L2Norm();

                        for(int k = 0; k < MaxPseudoTimeSteps; k++) {
                            var stw_k = Stopwatch.StartNew();
                            var OpMtx_k = OpMtx.CloneAs();
                            OpMtx_k.AccEyeSp(1.0 / dtau);

                            var RHS_k = RHS.CloneAs();
                            RHS_k.AccV(1.0 / dtau, u);

                            var mgOp_k = new MultigridOperator(AggBasis, ExtenstionSolVec.Mapping, OpMtx_k, null, MultigridOperatorConfig, GetBulkOperator());

                            var Pre_RHS_k = new double[mgOp_k.Mapping.LocalLength];
                            mgOp_k.TransformRhsInto(RHS_k, Pre_RHS_k, true);
                            var Pre_Sol_k = new double[mgOp_k.Mapping.LocalLength];
                            mgOp_k.TransformSolInto(u, Pre_Sol_k);

                            SolveOnce(mgOp_k, Pre_Sol_k, Pre_RHS_k);

                            mgOp_k.TransformSolFrom(u, Pre_Sol_k);

                            var Res = RHS.CloneAs();
                            OpMtx.SpMV(-1.0, u, 1.0, Res);
                            double resNorm = Res.MPI_L2Norm();
                            double rhsNorm = RHS.MPI_L2Norm();
                            double relResNorm = resNorm / Math.Max(rhsNorm, 1.0e-300);

                            ExtenstionSolVec.SetV(u, 1.0);
                            double ifaceDiffNorm = InterfaceVelocityMismatchNorm();

                            stw_k.Stop();
                            Console.WriteLine($"pseudo-time step {k}: dtau={dtau:0.###E-00}, |RHS|={rhsNorm:0.###E-00}, "
                                + $"true residual={resNorm:0.###E-00} (relative: {relResNorm:0.###E-00}), |u|={u.MPI_L2Norm():0.###E-00}, "
                                + $"|div(u_ext)|={ExtensionVelocity.Divergence().L2Norm():0.###E-00}, |u_ext - u_iface| on interface={ifaceDiffNorm:0.###E-00}, "
                                + $"time={stw_k.Elapsed.TotalSeconds:0.###}sec");

                            if(relResNorm < PseudoTimeConvergenceCriterion)
                                break;

                            dtau = Math.Min(dtau * resNorm_prev / resNorm, PseudoDtauMax);
                            resNorm_prev = resNorm;
                        }

                        ExtenstionSolVec.SetV(u, 1.0);
                    } else {
                        var stw = Stopwatch.StartNew();
                        MultigridOperator mgOp = new MultigridOperator(AggBasis, ExtenstionSolVec.Mapping, OpMtx, null, MultigridOperatorConfig, GetBulkOperator());

                        var Pre_RHS = RHS.CloneAs();
                        double[] Pre_Sol = new double[mgOp.Mapping.LocalLength];
                        tr.Info("StokesExtension.SolveExtension trying with Multigridoperator");
                        mgOp.TransformRhsInto(RHS, Pre_RHS, true);
                        // Warm start: 'ExtenstionSolVec' may already carry the previous time step's
                        // extension velocity (see callers of SolveExtension); it must be moved into
                        // the same equilibrated basis as 'Pre_RHS'/'OperatorMatrix', otherwise a nonzero
                        // initial guess would be handed to the solver in the wrong scaling. 'Pre_Sol' is
                        // zero-initialized above, so leaving this out when WarmStart is false reproduces
                        // the original cold-start behavior.
                        if(WarmStart)
                            mgOp.TransformSolInto(ExtenstionSolVec, Pre_Sol);

                        SolveOnce(mgOp, Pre_Sol, Pre_RHS);

                        mgOp.TransformSolFrom(ExtenstionSolVec, Pre_Sol);
                        double ifaceDiffNorm = InterfaceVelocityMismatchNorm();
                        stw.Stop();
                        Console.WriteLine($"|div(u_ext)|={ExtensionVelocity.Divergence().L2Norm():0.###E-00}, |u_ext - u_iface| on interface={ifaceDiffNorm:0.###E-00}, "
                            + $"time={stw.Elapsed.TotalSeconds:0.###}sec");
                    }
                    e = null;
                } catch(Exception ee) {
                    e = ee;
                }
                try {
                    e.ExceptionBcast(); // MPI-parallel exception handling 
                    
                } catch {
                    // should be replaced by something more sophisticated
                    //var Residual = RHS.CloneAs();
                    //var cExtenstionSolVec = ExtenstionSolVec.CloneAs();
                    tr.Warning($"StokesExtension.SolveExtension with Multigridoperator failed, trying direct method; got exception");
                    OpMtx.Solve_Direct(ExtenstionSolVec, bkup_RHS);
                }

                //{
                //    DGField[] stokesShit = VelocityAtInterface.Cat(DummySolFields);
                //    Tecplot.Tecplot.PlotFields(stokesShit, "StokesExt.ls" + levelSetIndex, 0.0, 3);
                //}



                //Console.WriteLine("Difference: {0}", ExtenstionSolVec.L2Dist(cExtenstionSolVec));

                // Console.WriteLine("StokesExt Vel L2: " + RHS.L2Norm());

                double[] OutputVelL2 = ExtensionVelocity.Select(vel => vel.L2Norm()).ToArray();

                // Tecplot.Tecplot.PlotFields(ExtensionVelocity.Cat(lsTrk.LevelSets[0] as LevelSet), this.GetType().ToString().Split('.').Last() + "-" + timestepNo, (double)timestepNo, 2);
                /*
                {
                    double RhsNorm = Residual.L2NormPow2().MPISum().Sqrt();
                    double MatrixInfNorm = OpMtx.InfNorm();
                    OpMtx.SpMV(-1.0, ExtenstionSolVec, 1.0, Residual);

                    double ResidualNorm = Residual.L2NormPow2().MPISum().Sqrt();
                    double SolutionNorm = ExtenstionSolVec.L2NormPow2().MPISum().Sqrt();
                    double Denom = Math.Max(MatrixInfNorm, Math.Max(RhsNorm, Math.Max(SolutionNorm, Math.Sqrt(BLAS.MachineEps))));
                    double RelResidualNorm = ResidualNorm / Denom;

                    //Console.WriteLine("done: Abs.: {0}, Rel.: {1}", ResidualNorm, RelResidualNorm);

                    if(RelResidualNorm > 1.0e-10) {
                        string ErrMsg;
                        using(var stw = new System.IO.StringWriter()) {
                            stw.WriteLine("Stokes Extension: High residual from direct solver.");
                            stw.WriteLine("    L2 Norm of RHS:         " + RhsNorm);
                            stw.WriteLine("    L2 Norm of Solution:    " + SolutionNorm);
                            stw.WriteLine("    L2 Norm of Residual:    " + ResidualNorm);
                            stw.WriteLine("    Relative Residual norm: " + RelResidualNorm);
                            stw.WriteLine("    Matrix Inf norm:        " + MatrixInfNorm);

                            ErrMsg = stw.ToString();
                        }
                        Console.Error.WriteLine(ErrMsg);

                        string curDir = System.IO.Directory.GetCurrentDirectory();
                        string failVault = System.IO.Path.Combine(curDir, "failVault_" + (DateTime.Now.Ticks));
                        Directory.CreateDirectory(failVault);
                        foreach(var plt in System.IO.Directory.GetFiles(curDir, "*.plt")) {
                            System.IO.File.Copy(plt, Path.Combine(failVault, Path.GetFileName(plt)));
                        }

                        OpMtx.SaveToTextFileSparse(Path.Combine(failVault, "StokesExtMtx.txt"));
                        RHS.SaveToTextFile(Path.Combine(failVault, "StokesExtRHS.txt"));
                        ExtenstionSolVec.SaveToTextFile(Path.Combine(failVault, "StokesExtSOL.txt"));

                        throw new ArithmeticException(ErrMsg);

                    }
                }*/


                // plotting for debug reasons
                //Console.WriteLine("StokesExtension Plotting enabled!");
                //Tecplot.Tecplot.PlotFields(ExtenstionSolVec.Fields.Cat(VelocityAtInterface), this.GetType().ToString().Split('.').Last() + "-" + timestepNo, (double)timestepNo, 2);
                //timestepNo++;
            }
        }
    }
}
