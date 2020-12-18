using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {
    
    /// <summary>
    /// Driver routines for the solver framework
    /// </summary>
    public static class Extensions {


        /// <summary>
        /// General-Purpose driver Routine for performing a steady-state solution of a PDE defined through some operator
        /// </summary>
        /// <param name="op">
        /// The spatial operator for the PDE to be solved (we solve `$ \mathrm(op)(U) = 0 `$, where `$ U `$ stands for <paramref name="Solution"/>).
        /// In order to obtain a non-zero solution, this must typically contain some source terms or inhomogeneous boundary conditions.
        /// </param>
        /// <param name="Solution">
        /// DG fields to store the solution.
        /// </param>
        /// <param name="nsc">
        /// - configuration of the nonlinear solver (see also <see cref="BoSSS.Solution.Control.AppControlSolver.NonLinearSolver"/>)
        /// - if null, an default solver configuration is used
        /// </param>
        /// <param name="lsc">
        /// - configuration of the linear solver (see also <see cref="BoSSS.Solution.Control.AppControlSolver.LinearSolver"/>)
        /// - if null, an default solver configuration is used 
        /// </param>
        /// <param name="MultigridSequence">
        /// Multigrid sequence/hierarchy on which a multigrid solver should operate;
        /// Providing this does not guarantee that a multigrid solver is used, this also depends on the other solver settings.
        /// </param>
        /// <param name="verbose">
        /// - If true, Writes a lot of logging information
        /// - If false, console output should be relatively less
        /// </param>
        /// <param name="MgConfig">
        /// Provisional: will be integrated into the <see cref="ISpatialOperator"/> at some point.
        /// </param>
        static public void Solve(this ISpatialOperator op, CoordinateMapping Solution, MultigridOperator.ChangeOfBasisConfig[][] MgConfig, NonLinearSolverConfig nsc = null, LinearSolverConfig lsc = null, AggregationGridData[] MultigridSequence = null, bool verbose = false) {
            using(var tr = new FuncTrace()) {

                // init
                // ====

                Stopwatch stw = new Stopwatch();
                stw.Start();
                if(verbose) {
                    Console.WriteLine($"Trying to solve equation (starting at {DateTime.Now})...");
                }

                int NoOfVar = Solution.BasisS.Count;

                if(op.DomainVar.Count != op.CodomainVar.Count)
                    throw new ArgumentException("Operator is not square: Number of Domain variables differs form number of Codomain variables; This is not supported.");

                if(op.DomainVar.Count != NoOfVar)
                    throw new ArgumentException("Mismatch between number of Domain variables in operator and number of DG fields in solution mapping");

                DGField[] SolutionFields = Solution.Fields.ToArray();
                DGField[] ParamFields = op.InvokeParameterFactory(SolutionFields);

                double[] SolutionVec = new double[Solution.LocalLength];

                IGridData gdat = Solution.GridDat;

                if(op.TemporalOperator != null)
                    throw new ArgumentException("Operator contains a temporal component - cannot solve this with a steady-state solver.");

                LevelSetTracker LsTrk = null;
                foreach(var b in Solution.BasisS) {
                    if(b is XDGBasis xb)
                        LsTrk = xb.Tracker;
                }

                if(verbose) {
                    if(LsTrk == null)
                        Console.WriteLine("Solver running in DG mode (no XDG).");
                    else
                        Console.WriteLine("Solver running in XDG mode.");
                }

                XSpatialOperatorMk2 xop = op as XSpatialOperatorMk2;

                string[] speciesNames;
                if(xop != null) {
                    speciesNames = xop.Species.ToArray();
                } else {
                    if(LsTrk != null)
                        speciesNames = LsTrk.SpeciesNames.ToArray();
                    else
                        speciesNames = null;
                }

                if(verbose && speciesNames != null) {
                    Console.Write("Solving for species: ");
                    for(int i = 0; i < speciesNames.Length; i++) {
                        Console.Write(speciesNames[i]);
                        if(i < speciesNames.Length - 1)
                            Console.Write(", ");
                    }
                    Console.WriteLine(";");
                }

                SpeciesId[] spcIDs;
                if(speciesNames != null && LsTrk != null) {
                    spcIDs = speciesNames.Select(spcN => LsTrk.GetSpeciesId(spcN)).ToArray();
                } else {
                    spcIDs = null;
                }

                int quadOrder = op.GetOrderFromQuadOrderFunction(Solution.BasisS, ParamFields.Select(p => p != null ? p.Basis : null), Solution.BasisS);
                if(verbose) {
                    Console.WriteLine($"Using quadrature order {quadOrder}.");
                }



                // Verify or Create Multigrid sequence
                // ===================================

                if(MultigridSequence != null) {
                    if(!object.ReferenceEquals(MultigridSequence[0].ParentGrid, gdat))
                        throw new ArgumentException("Multigrid parent must be the same grid on which the solution is defined.");
                    Console.WriteLine($"{MultigridSequence.Length} multigrid levels available.");
                } else {
                    Console.WriteLine("No multigrid sequence available.");
                    using(new BlockTrace("Aggregation_grid_construction", tr)) {
                        MultigridSequence = new AggregationGridData[] { CoarseningAlgorithms.ZeroAggregation(gdat) };
                    }
                }

                // Create Agglomeration
                // ====================

                MultiphaseCellAgglomerator Op_Agglomeration;
                if(LsTrk != null && xop != null) {
                    Op_Agglomeration = LsTrk.GetAgglomerator(spcIDs, quadOrder, 0.1);
                } else {
                    Op_Agglomeration = null;
                }

                // create aggregation basis
                // ========================

                AggregationGridBasis[][] AggBasisS;
                using (new BlockTrace("Aggregation_basis_init", tr)) {
                    AggBasisS = AggregationGridBasis.CreateSequence(MultigridSequence, Solution.BasisS);
                    AggBasisS.UpdateXdgAggregationBasis(Op_Agglomeration);
                }


                // solver
                // ======

                bool Converged = false;
                if(op.IsLinear) {
                    // +++++++++++++++++++++
                    // solve a linear system
                    // +++++++++++++++++++++

                    if(verbose)
                        Console.WriteLine("Operator/PDE is linear.");


                    // assemble the linear system
                    // --------------------------
                    Console.WriteLine($"creating sparse system for {Solution.GlobalCount} DOF's ...");
                    Stopwatch stwAssi = new Stopwatch();
                    stwAssi.Start();

                    BlockMsrMatrix opMtx;
                    double[] opAff;
                    using(new BlockTrace("MatrixAssembly", tr)){
                        

                        var MtxBuilder = op.GetMatrixBuilder(Solution, ParamFields, Solution);

                        op.InvokeParameterUpdate(SolutionFields, ParamFields);

                        opMtx = new BlockMsrMatrix(Solution, Solution);
                        opAff = new double[Solution.LocalLength];
                        MtxBuilder.ComputeMatrix(opMtx, opAff);


                        // agglomeration wahnsinn
                        if(Op_Agglomeration != null) {

                            Op_Agglomeration.ManipulateMatrixAndRHS(opMtx, opAff, Solution, Solution);

                            if(verbose) {
                                foreach(var S in speciesNames) {
                                    Console.WriteLine("  Species {0}: no of agglomerated cells: {1}",
                                        S, Op_Agglomeration.GetAgglomerator(LsTrk.GetSpeciesId(S)).AggInfo.SourceCells.NoOfItemsLocally);
                                }
                            }
                        }
                    }

                    BlockMsrMatrix MassMatrix;
                    using(new BlockTrace("Mass_Matrix_comp", tr)) {
                        if(LsTrk != null) {
                            

                            var massFact = LsTrk.GetXDGSpaceMetrics(spcIDs, quadOrder).MassMatrixFactory;
                            MassMatrix = massFact.GetMassMatrix(Solution, NoOfVar.ForLoop(iVar => 1.0), false, spcIDs);

                            Op_Agglomeration.ManipulateMatrixAndRHS(MassMatrix, default(double[]), Solution, Solution);
                        } else {
                            MassMatrix = null;
                        }
                    }

                    

                    stwAssi.Stop();
                    if(verbose)
                        Console.WriteLine("done {0} sec.", stwAssi.Elapsed.TotalSeconds);

                    if(verbose) {
                        int minBlkSize = Solution.MaxTotalNoOfCoordinatesPerCell;
                        int maxBlkSize = Solution.MinTotalNoOfCoordinatesPerCell;
                        int NoOfMtxBlocks = 0;
                        foreach(int[] Neigs in gdat.iLogicalCells.CellNeighbours) {
                            NoOfMtxBlocks++; //               diagonal block
                            NoOfMtxBlocks += Neigs.Length; // off-diagonal block
                        }
                        NoOfMtxBlocks = NoOfMtxBlocks.MPISum();

                        opMtx.GetMemoryInfo(out long allc, out long used);

                        long MtxSize = opMtx.GetTotalNoOfNonZeros();

                        double MtxStorage = allc / (1024 * 1024); // 12 bytes (double+int) per entry

                        Console.WriteLine("   System size:                 {0}", Solution.TotalLength);
                        Console.WriteLine("   No of blocks:                {0}", Solution.TotalNoOfBlocks);
                        Console.WriteLine("   No of blocks in matrix:      {0}", NoOfMtxBlocks);
                        Console.WriteLine("   min DG coordinates per cell: {0}", minBlkSize);
                        Console.WriteLine("   max DG coordinates per cell: {0}", maxBlkSize);
                        Console.WriteLine("   Total non-zeros in matrix:   {0}", MtxSize);
                        Console.WriteLine("   Approx. matrix storage (MB): {0}", MtxStorage);


                       // base.QueryHandler.ValueQuery("MtxBlkSz", MtxBlockSize, true);
                       // base.QueryHandler.ValueQuery("NNZMtx", MtxSize, true);
                       // base.QueryHandler.ValueQuery("NNZblk", NoOfMtxBlocks, true);
                       // base.QueryHandler.ValueQuery("MtxMB", MtxStorage, true);
                    }




                    // setup of multigrid operator
                    // ---------------------------

                    var solverSetup = new Stopwatch();
                    solverSetup.Start();
                    ISolverSmootherTemplate solver;

                    

                    if(verbose)
                        Console.WriteLine("Setting up multigrid operator...");

                    var MultigridOp = new MultigridOperator(AggBasisS, Solution, 
                        opMtx, MassMatrix, MgConfig, 
                        op.DomainVar.Select(varName => op.FreeMeanValue[varName]).ToArray());


                    SolverFactory SF = new SolverFactory(nsc, lsc);
                    SF.GenerateLinear(out solver, AggBasisS, MgConfig);
                                        
                    using(new BlockTrace("Solver_Init", tr)) {
                        solver.Init(MultigridOp);
                    }

                    solverSetup.Stop();
                    if(verbose)
                        Console.WriteLine("done. (" + solverSetup.Elapsed.TotalSeconds + " sec)");

                    if(verbose) {
                        MultigridOp.GetMemoryInfo(out long AllocMem, out long UsedMem);
                        Console.WriteLine("  Memory reserved|used by multi-grid operator {0:F2} | {1:F2} MB", (double)AllocMem / (1024.0 * 1024.0), (double)UsedMem / (1024.0 * 1024.0));
                    }

                    // call the linear solver
                    // ----------------------

                    if(verbose)
                        Console.WriteLine("Running solver...");
                    var solverIteration = new Stopwatch();
                    solverIteration.Start();
                    using(new BlockTrace("Solver_Run", tr)) {
                        solver.ResetStat();

                        var RHSvec = opAff.CloneAs();
                        RHSvec.ScaleV(-1);
                        
                        MultigridOp.UseSolver(solver, SolutionVec, RHSvec);
                    }
                    solverIteration.Stop();
                    if(verbose) {
                        Console.WriteLine("done. (" + solverIteration.Elapsed.TotalSeconds + " sec)");
                    }


                    Converged = solver.Converged;


                } else {
                    // ++++++++++++++++++++
                    // use nonlinear solver
                    // ++++++++++++++++++++

                    if(verbose)
                        Console.WriteLine("Operator/PDE is nonlinear.");



                    throw new NotImplementedException("Nonlinear solver: todo");

                }

                // Finalize
                // ========

                if(verbose) {
                    Console.WriteLine("  Pardiso phase 11: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_11.Elapsed.TotalSeconds);
                    Console.WriteLine("  Pardiso phase 22: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_22.Elapsed.TotalSeconds);
                    Console.WriteLine("  Pardiso phase 33: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_33.Elapsed.TotalSeconds);
                    Console.WriteLine("  spmm total " + BlockMsrMatrix.multiply.Elapsed.TotalSeconds);
                    Console.WriteLine("  spmm core " + BlockMsrMatrix.multiply_core.Elapsed.TotalSeconds);
                    Console.WriteLine("  spmv total " + BlockMsrMatrix.SPMV_tot.Elapsed.TotalSeconds);
                    Console.WriteLine("  spmv inner " + BlockMsrMatrix.SPMV_inner.Elapsed.TotalSeconds);
                    Console.WriteLine("  spmv outer " + BlockMsrMatrix.SpMV_local.Elapsed.TotalSeconds);
                }

                // store solution
                {
                    (new CoordinateVector(Solution)).SetV(SolutionVec);
                    if(Op_Agglomeration != null)
                        Op_Agglomeration.Extrapolate(Solution);
                }
                // time measurement, statistics
                stw.Stop();
                if(verbose) {
                    Console.WriteLine($"Total runtime: {stw.Elapsed}");

                    if(Converged)
                        Console.WriteLine("Solver SUCCESSFUL!");
                    else
                        Console.WriteLine("Solver FAILED!");
                }

            }
        }
    
    
    
    
        /// <summary>
        /// Easy-to-use driver routine for operator analysis
        /// </summary>
        /// <param name="op"></param>
        /// <param name="Mapping">
        /// 
        /// </param>
        /// <returns></returns>
        /// <seealso cref="BoSSS.Solution.Application{T}.OperatorAnalysis"/>
        static public IDictionary<string,double> OperatorAnalysis(this ISpatialOperator op, CoordinateMapping Mapping) {

            /*
             var ana = new BoSSS.Solution.AdvancedSolvers.Testing.OpAnalysisBase(this.LsTrk, 
                this.Op_Matrix, this.Op_Affine, 
                Mapping, Op_Agglomeration, 
                this.Op_mass.GetMassMatrix(this.u.Mapping, new double[] { 1.0 }, false, this.LsTrk.SpeciesIdS.ToArray()), 
                this.OpConfig, this.Op);

            return ana.GetNamedProperties();
            */
            throw new NotImplementedException("fucking todo");
        }
    }
}
