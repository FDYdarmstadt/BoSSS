using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
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
    public static class UniSolver {


        class MatrixAssembler {

            public MatrixAssembler(ISpatialOperator __op, CoordinateMapping Solution, AggregationGridData[] __MultigridSequence = null, QueryHandler __queryHandler = null, MultigridOperator.ChangeOfBasisConfig[][] __MgConfig = null) {
                using(var tr = new FuncTrace()) {
                    // init
                    // ====
                    this.op = __op;
                    this.queryHandler = __queryHandler;
                    NoOfVar = Solution.BasisS.Count;
                    SolMapping = Solution;

                    if(op.DomainVar.Count != op.CodomainVar.Count)
                        throw new ArgumentException("Operator is not square: Number of Domain variables differs form number of Codomain variables; This is not supported.");

                    if(op.DomainVar.Count != NoOfVar)
                        throw new ArgumentException("Mismatch between number of Domain variables in operator and number of DG fields in solution mapping");

                    SolutionFields = Solution.Fields.ToArray();

                    SolutionVec = new double[Solution.LocalLength];

                    ParamFields = op.InvokeParameterFactory(Solution.Fields);

                    gdat = Solution.GridDat;

                    //if(op.TemporalOperator != null)
                    //    throw new ArgumentException("Operator contains a temporal component - cannot solve this with a steady-state solver.");

                    LsTrk = null;
                    foreach(var b in Solution.BasisS) {
                        if(b is XDGBasis xb)
                            LsTrk = xb.Tracker;
                    }
                    if(NonXDG_LsTrk != null) {
                        if(LsTrk == null) {
                            LsTrk = NonXDG_LsTrk;
                        } else {
                            if(!object.ReferenceEquals(LsTrk, NonXDG_LsTrk))
                                throw new ArgumentException("Level Set Tracker mismatch: Someone has used the freak API and did not cleaned it up properly.");
                        }
                    }


                    xop = op as XSpatialOperatorMk2;


                    if(xop != null) {
                        speciesNames = xop.Species.ToArray();
                    } else {
                        if(LsTrk != null)
                            speciesNames = LsTrk.SpeciesNames.ToArray();
                        else
                            speciesNames = null;
                    }

                    if(speciesNames != null && LsTrk != null) {
                        spcIDs = speciesNames.Select(spcN => LsTrk.GetSpeciesId(spcN)).ToArray();
                    } else {
                        spcIDs = null;
                    }

                    quadOrder = op.GetOrderFromQuadOrderFunction(Solution.BasisS, ParamFields.Select(p => p != null ? p.Basis : null), Solution.BasisS);

                    if(__MgConfig == null) {
                        MgConfig = new MultigridOperator.ChangeOfBasisConfig[1][];
                        MgConfig[0] = new MultigridOperator.ChangeOfBasisConfig[NoOfVar];

                        for(int iVar = 0; iVar < NoOfVar; iVar++) {
                            MgConfig[0][iVar] = new MultigridOperator.ChangeOfBasisConfig() {
                                VarIndex = new int[] { iVar },
                                DegreeS = new int[] { SolutionFields[iVar].Basis.Degree },
                                mode = (SolutionFields[iVar].Basis is XDGBasis) ? MultigridOperator.Mode.IdMass_DropIndefinite : MultigridOperator.Mode.Eye
                            };
                        }
                    } else {
                        MgConfig = __MgConfig;
                    }



                    // Verify or Create Multigrid sequence
                    // ===================================
                    MultigridSequence = __MultigridSequence;
                    if(MultigridSequence != null) {
                        if(!object.ReferenceEquals(MultigridSequence[0].ParentGrid, gdat))
                            throw new ArgumentException("Multigrid parent must be the same grid on which the solution is defined.");
                    } else {
                        using(new BlockTrace("Aggregation_grid_construction", tr)) {
                            MultigridSequence = new AggregationGridData[] { CoarseningAlgorithms.ZeroAggregation(gdat) };
                        }
                    }

                    // Create Agglomeration
                    // ====================


                    if(LsTrk != null && xop != null) {
                        Op_Agglomeration = LsTrk.GetAgglomerator(spcIDs, quadOrder, xop.AgglomerationThreshold);

                        foreach(SpeciesId s in Op_Agglomeration.SpeciesList) {
                            var agg = Op_Agglomeration.GetAgglomerator(s);
                            Console.WriteLine($"Agglom. for species {LsTrk.GetSpeciesName(s)}: {agg.TotalNumberOfAgglomerations}");
                                                       

                            Console.WriteLine(agg.AggInfo.ToString());
                        }


                    } else {
                        Op_Agglomeration = null;
                    }

                    // create aggregation basis
                    // ========================

                    using(new BlockTrace("Aggregation_basis_init", tr)) {
                        AggBasisS = AggregationGridBasis.CreateSequence(MultigridSequence, Solution.BasisS);
                        AggBasisS.UpdateXdgAggregationBasis(Op_Agglomeration);
                    }

                    if(queryHandler != null) {
                        var BS00 = AggBasisS[0][0];
                        int J = BS00.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                        int cntCutCellBlocks = 0;
                        int cutCellDiagBlocks = 0;
                        for(int j = 0; j < J; j++) {
                            int[] Neigs = BS00.AggGrid.iLogicalCells.CellNeighbours[j];

                            if(BS00.GetNoOfSpecies(j) > 1) {
                                cutCellDiagBlocks++;

                                foreach(int jN in Neigs) {
                                    if(BS00.GetNoOfSpecies(jN) > 1)
                                        cntCutCellBlocks++;
                                }
                            }
                        }

                        queryHandler.ValueQuery("NoOfCutCellBlocks", cntCutCellBlocks, true);
                        queryHandler.ValueQuery("NoOfCutCellDiagBlocks", cutCellDiagBlocks, true);
                    }
                }
            }

            public QueryHandler queryHandler;

            public int NoOfVar;

            public IGridData gdat;

            public DGField[] ParamFields;

            DGField[] JacobiParameterVars;

            public double[] SolutionVec;

            public MultiphaseCellAgglomerator Op_Agglomeration;

            public AggregationGridBasis[][] AggBasisS;

            public DGField[] SolutionFields;

            public string[] speciesNames;

            public SpeciesId[] spcIDs;

            public LevelSetTracker LsTrk;

            public XSpatialOperatorMk2 xop;

            public int quadOrder;

            public ISpatialOperator op;

            public bool verbose = false;

            AggregationGridData[] MultigridSequence;

            public CoordinateMapping SolMapping;

            public MultigridOperator.ChangeOfBasisConfig[][] MgConfig;

            public double time = 1.0e30;

            /// <summary>
            /// Implementation of <see cref="OperatorEvalOrLin"/>
            /// </summary>
            public void AssembleMatrix(out BlockMsrMatrix opMtx, out double[] opAff, out BlockMsrMatrix MassMatrix, DGField[] CurrentState, bool Linearization, out ISpatialOperator OberFrickelHack) {
                using(var tr = new FuncTrace()) {

                    var Solution = new CoordinateMapping(CurrentState);
                    OberFrickelHack = op;
                    if(!Solution.EqualsPartition(this.SolMapping))
                        throw new ApplicationException("something is weird.");

                    if(op.IsLinear && op.LinearizationHint != LinearizationHint.AdHoc)
                        throw new NotSupportedException("Configuration Error: for a supposedly linear operator, the linearization hint must be " + LinearizationHint.AdHoc);


                    // assemble the linear system
                    // --------------------------
                    Console.WriteLine($"creating sparse system for {Solution.GlobalCount} DOF's ...");
                    Stopwatch stwAssi = new Stopwatch();
                    stwAssi.Start();

                    if(Linearization) {
                        using(new BlockTrace("MatrixAssembly", tr)) {
                            opMtx = new BlockMsrMatrix(Solution, Solution);
                            opAff = new double[Solution.LocalLength];
                            if(xop != null) {
                                this.XDGMatrixAssembly(opMtx, opAff);
                            } else {
                                this.DGMatrixAssembly(opMtx, opAff);
                            }


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

                        stwAssi.Stop();
                        if(verbose)
                            Console.WriteLine("done {0} sec.", stwAssi.Elapsed.TotalSeconds);

                        if(verbose || queryHandler != null) {
                            int minBlkSize = Solution.MinTotalNoOfCoordinatesPerCell;
                            int maxBlkSize = Solution.MaxTotalNoOfCoordinatesPerCell;
                            int NoOfMtxBlocks = 0;
                            int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
                            for(int j = 0; j < J; j++) {
                                int[] Neigs = gdat.iLogicalCells.CellNeighbours[j];

                                NoOfMtxBlocks++; //               diagonal block
                                NoOfMtxBlocks += Neigs.Length; // off-diagonal block


                            }
                            NoOfMtxBlocks = NoOfMtxBlocks.MPISum();



                            opMtx.GetMemoryInfo(out long allc, out long used);

                            long MtxSize = opMtx.GetTotalNoOfNonZeros();

                            double MtxStorage = allc / (1024 * 1024); // 12 bytes (double+int) per entry

                            if(verbose) {
                                Console.WriteLine("   System size:                 {0}", Solution.TotalLength);
                                Console.WriteLine("   No of blocks:                {0}", Solution.TotalNoOfBlocks);
                                Console.WriteLine("   No of blocks in matrix:      {0}", NoOfMtxBlocks);
                                Console.WriteLine("   min DG coordinates per cell: {0}", minBlkSize);
                                Console.WriteLine("   max DG coordinates per cell: {0}", maxBlkSize);
                                Console.WriteLine("   Total non-zeros in matrix:   {0}", MtxSize);
                                Console.WriteLine("   Approx. matrix storage (MB): {0}", MtxStorage);
                            }

                            if(queryHandler != null) {
                                queryHandler.ValueQuery("NNZMtx", MtxSize, true);
                                queryHandler.ValueQuery("NNZblk", NoOfMtxBlocks, true);
                                queryHandler.ValueQuery("MtxMB", MtxStorage, true);
                                queryHandler.ValueQuery("NumberOfMatrixBlox", NoOfMtxBlocks, true); // the same data under two names, not to crash any worksheet.
                            }
                        }
                    } else {
                        // +++++++++++++++
                        // just evaluation 
                        // +++++++++++++++


                        opAff = new double[Solution.LocalLength];
                        if(xop != null) {
                            this.XDGevaluation(opAff);
                        } else {
                            this.DGevaluation(opAff);
                        }

                        opMtx = null;

                    }

                    // Assembly of mass matrix
                    // -----------------------

                    using(new BlockTrace("Mass_Matrix_comp", tr)) {
                        if(LsTrk != null) {


                            var massFact = LsTrk.GetXDGSpaceMetrics(spcIDs, quadOrder).MassMatrixFactory;
                            MassMatrix = massFact.GetMassMatrix(Solution, NoOfVar.ForLoop(iVar => 1.0), false, spcIDs);

                            Op_Agglomeration.ManipulateMatrixAndRHS(MassMatrix, default(double[]), Solution, Solution);
                        } else {
                            MassMatrix = null;
                        }
                    }
                }
            }

            void XDGMatrixAssembly(BlockMsrMatrix OpMtx, double[] OpAffine) {

                var AgglomeratedCellLengthScales = this.Op_Agglomeration.CellLengthScales;

                switch(xop.LinearizationHint) {

                    case LinearizationHint.AdHoc: {
                        xop.InvokeParameterUpdate(this.time, this.SolutionFields, this.ParamFields.ToArray());

                        var mtxBuilder = xop.GetMatrixBuilder(LsTrk, this.SolMapping, this.ParamFields, this.SolMapping);
                        mtxBuilder.time = 0.0;
                        mtxBuilder.MPITtransceive = true;
                        foreach(var kv in AgglomeratedCellLengthScales) {
                            mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                        }
                        mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                        return;
                    }

                    case LinearizationHint.FDJacobi: {
                        var mtxBuilder = xop.GetFDJacobianBuilder(LsTrk, this.SolMapping, this.ParamFields, this.SolMapping);
                        mtxBuilder.time = 0.0;
                        mtxBuilder.MPITtransceive = true;
                        if(mtxBuilder.Eval is XSpatialOperatorMk2.XEvaluatorNonlin evn) {
                            foreach(var kv in AgglomeratedCellLengthScales) {
                                evn.CellLengthScales[kv.Key] = kv.Value;
                            }
                        }
                        mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                        return;
                    }

                    case LinearizationHint.GetJacobiOperator: {
                        var JacXop = xop.GetJacobiOperator(gdat.SpatialDimension) as XSpatialOperatorMk2;

                        if(JacobiParameterVars == null)
                            JacobiParameterVars = JacXop.InvokeParameterFactory(this.SolutionFields);

                        JacXop.InvokeParameterUpdate(this.time, this.SolutionFields, JacobiParameterVars);

                        var mtxBuilder = JacXop.GetMatrixBuilder(LsTrk, this.SolMapping, this.JacobiParameterVars, this.SolMapping);
                        mtxBuilder.time = 0.0;
                        mtxBuilder.MPITtransceive = true;
                        foreach(var kv in AgglomeratedCellLengthScales) {
                            mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                        }
                        mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                        return;
                    }
                }
            }

            void DGMatrixAssembly(BlockMsrMatrix OpMtx, double[] OpAffine) {
                switch(op.LinearizationHint) {

                    case LinearizationHint.AdHoc: {
                        this.op.InvokeParameterUpdate(this.time, this.SolutionFields, this.ParamFields);

                        var mtxBuilder = op.GetMatrixBuilder(this.SolMapping, this.ParamFields, this.SolMapping);
                        mtxBuilder.time = 0.0;
                        mtxBuilder.MPITtransceive = true;
                        mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                        return;
                    }

                    case LinearizationHint.FDJacobi: {
                        var mtxBuilder = op.GetFDJacobianBuilder(this.SolMapping, this.ParamFields, this.SolMapping);
                        mtxBuilder.time = 0.0;
                        mtxBuilder.MPITtransceive = true;
                        mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                        return;
                    }

                    case LinearizationHint.GetJacobiOperator: {
                        var JacXop = op.GetJacobiOperator(gdat.SpatialDimension);

                        if(JacobiParameterVars == null)
                            JacobiParameterVars = JacXop.InvokeParameterFactory(this.SolutionFields);
                        JacXop.InvokeParameterUpdate(this.time, this.SolutionFields, JacobiParameterVars);

                        var mtxBuilder = JacXop.GetMatrixBuilder(this.SolMapping, this.JacobiParameterVars, this.SolMapping);
                        mtxBuilder.time = 0.0;
                        mtxBuilder.MPITtransceive = true;
                        mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                        return;
                    }
                }
            }

            void XDGevaluation(double[] OpAffine) {
                this.xop.InvokeParameterUpdate(this.time, this.SolutionFields, this.ParamFields);
                var AgglomeratedCellLengthScales = this.Op_Agglomeration.CellLengthScales;

                var eval = xop.GetEvaluatorEx(this.LsTrk, this.SolMapping, this.ParamFields, this.SolMapping);
                eval.time = 0.0;

                eval.MPITtransceive = true;
                foreach(var kv in AgglomeratedCellLengthScales) {
                    eval.CellLengthScales[kv.Key] = kv.Value;
                }
                eval.Evaluate(1.0, 0.0, OpAffine);
            }

            void DGevaluation(double[] OpAffine) {
                this.op.InvokeParameterUpdate(this.time, this.SolutionFields, this.ParamFields);

                var eval = op.GetEvaluatorEx(this.SolMapping, this.ParamFields, this.SolMapping);
                eval.time = 0.0;

                eval.MPITtransceive = true;
                eval.Evaluate(1.0, 0.0, OpAffine);
            }
        }


        /// <summary>
        /// Freak API:
        /// For the very special case, that an XDG-computation should be performed on DG fields (there are some examples),
        /// this can be used to tell the <see cref="Solve"/>-routine which Tracker to use.
        /// This is an **intentionally bad, but good design, because freak cases should never be allowed to spoil good APIs.** 
        /// </summary>
        public static LevelSetTracker NonXDG_LsTrk;


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
        /// - configuration of the nonlinear solver 
        /// - if null, an default solver configuration is used
        /// </param>
        /// <param name="lsc">
        /// - configuration of the linear solver 
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
        /// <param name="queryHandler">
        /// if provided, logging of solver statistics to the unified query table used with the BoSSS database
        /// </param>
        /// <param name="optRHS">
        /// An optional right-hand-side for the system; per Definition, it is independent of the solution.
        /// If provided, the mapping must be partition-equal (<see cref="Partitioning_Extensions.EqualsPartition(IPartitioning, IPartitioning)"/>) to <paramref name="Solution"/>
        /// </param>
        /// <returns>
        /// true on solver success
        /// </returns>
        static public bool Solve(this ISpatialOperator op, CoordinateMapping Solution, CoordinateMapping optRHS = null,
            MultigridOperator.ChangeOfBasisConfig[][] MgConfig = null,
            NonLinearSolverConfig nsc = null, ISolverFactory lsc = null,
            AggregationGridData[] MultigridSequence = null,
            bool verbose = false, QueryHandler queryHandler = null) {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = verbose;

                // init
                // ====

                Stopwatch stw = new Stopwatch();
                stw.Start();
                tr.Info($"Trying to solve equation (starting at {DateTime.Now})...");
                

                if(optRHS != null)
                    optRHS.EqualsPartition(Solution);

                if(nsc == null) {
                    nsc = new NonLinearSolverConfig() {
                        SolverCode = NonLinearSolverCode.Newton,
                        Globalization = Newton.GlobalizationOption.Dogleg,
                    };
                }

                if(lsc == null) {
                    lsc = new DirectSolver.Config() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO
                    };
                }

                var G = new MatrixAssembler(op, Solution, MultigridSequence, queryHandler, __MgConfig: MgConfig);
                G.verbose = verbose;

                if(verbose) {
                    if(G.LsTrk == null)
                        tr.Info("Solver running in DG mode (no XDG).");
                    else
                        tr.Info("Solver running in XDG mode.");


                    if(G.speciesNames != null) {
                        Console.Write("Solving for species: ");
                        for(int i = 0; i < G.speciesNames.Length; i++) {
                            Console.Write(G.speciesNames[i]);
                            if(i < G.speciesNames.Length - 1)
                                Console.Write(", ");
                        }
                        Console.WriteLine(";");
                    }

                    tr.Info($"Using quadrature order {G.quadOrder}.");
                    

                    if(MultigridSequence != null) {
                        tr.Info($"{MultigridSequence.Length} multigrid levels available.");
                    } else {
                        tr.Info("No multigrid sequence available - using only one mesh level.");
                    }
                }

                // solver
                // ======

                bool Converged = false;
                int NoOfIterations = 0;
                long TotalDOF = 0;
                if(op.IsLinear) {
                    // +++++++++++++++++++++
                    // solve a linear system
                    // +++++++++++++++++++++

                    tr.Info("Operator/PDE is linear.");

                    G.verbose = verbose;
                    G.AssembleMatrix(out var opMtx, out double[] opAff, out var MassMatrix, G.SolutionFields, true, out _);



                    // setup of multigrid operator
                    // ---------------------------

                    var solverSetup = new Stopwatch();
                    solverSetup.Start();

                    tr.Info("Setting up multigrid operator...");

                    var MultigridOp = new MultigridOperator(G.AggBasisS, Solution,
                        opMtx, MassMatrix, G.MgConfig,
                        op);

                    //SolverFactory SF = new SolverFactory(nsc, lsc);
                    //SF.GenerateLinear(out solver, G.AggBasisS, G.MgConfig);
                    var Solver_Init = new BlockTrace("Solver_Init", tr);
                    using (ISolverSmootherTemplate solver = lsc.CreateInstance(MultigridOp)) {
                        Solver_Init.Dispose();

                        if(solver is ISolverWithCallback cl) {
                            cl.IterationCallback = delegate (int iIter, double[] X, double[] Res, MultigridOperator _) {
                                double ResNorm = Res.MPI_L2Norm();
                                Console.WriteLine($"{iIter} {ResNorm:0.##e-00}");
                            };
                        }


                        solverSetup.Stop();
                        tr.Info("done. (" + solverSetup.Elapsed.TotalSeconds + " sec)");

                        if(verbose) {
                            MultigridOp.GetMemoryInfo(out long AllocMem, out long UsedMem);
                            tr.Info(string.Format("  Memory reserved|used by multi-grid operator {0:F2} | {1:F2} MB", (double)AllocMem / (1024.0 * 1024.0), (double)UsedMem / (1024.0 * 1024.0)));
                        }
                        //LastMtx = MultigridOp.OperatorMatrix.CloneAs();

                        // call the linear solver
                        // ----------------------

                        tr.Info("Running solver...");
                        var solverIteration = new Stopwatch();
                        solverIteration.Start();
                        using(new BlockTrace("Solver_Run", tr)) {
                            solver.ResetStat();

                            var RHSvec = opAff.CloneAs();
                            RHSvec.ScaleV(-1);

                            if(optRHS != null)
                                RHSvec.AccV(1.0, new CoordinateVector(optRHS));

                            MultigridOp.UseSolver(solver, G.SolutionVec, RHSvec);

                            NoOfIterations = solver.ThisLevelIterations;
                            TotalDOF = MultigridOp.Mapping.TotalLength;
                        }
                        solverIteration.Stop();
                        tr.Info("done. (" + solverIteration.Elapsed.TotalSeconds + " sec, " + NoOfIterations + " iter)");
                        


                        Converged = solver.Converged;

                        
                    }
                } else {
                    // ++++++++++++++++++++
                    // use nonlinear solver
                    // ++++++++++++++++++++

                    tr.Info("Operator/PDE is nonlinear.");



                    throw new NotImplementedException("Nonlinear solver: todo");

                }

                // Finalize
                // ========

                if(verbose) {
                    tr.Info("  Pardiso phase 11: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_11.Elapsed.TotalSeconds);
                    tr.Info("  Pardiso phase 22: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_22.Elapsed.TotalSeconds);
                    tr.Info("  Pardiso phase 33: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_33.Elapsed.TotalSeconds);
                    tr.Info("  spmm total " + BlockMsrMatrix.multiply.Elapsed.TotalSeconds);
                    tr.Info("  spmm core " + BlockMsrMatrix.multiply_core.Elapsed.TotalSeconds);
                    tr.Info("  spmv total " + BlockMsrMatrix.SPMV_tot.Elapsed.TotalSeconds);
                    tr.Info("  spmv inner " + BlockMsrMatrix.SPMV_inner.Elapsed.TotalSeconds);
                    tr.Info("  spmv outer " + BlockMsrMatrix.SpMV_local.Elapsed.TotalSeconds);
                }

                // store solution
                {
                    (new CoordinateVector(Solution)).SetV(G.SolutionVec);
                    if(G.Op_Agglomeration != null)
                        G.Op_Agglomeration.Extrapolate(Solution);
                }
                // time measurement, statistics
                stw.Stop();
                if(verbose) {
                    tr.Info($"Total runtime: {stw.Elapsed}");

                    if(Converged)
                        tr.Info("Solver SUCCESSFUL!");
                    else
                        tr.Info("Solver FAILED!");
                }


                if(G.queryHandler != null) {
                    G.queryHandler.ValueQuery("minSolRunT", stw.Elapsed.TotalSeconds, true);
                    G.queryHandler.ValueQuery("maxSolRunT", stw.Elapsed.TotalSeconds, true);
                    G.queryHandler.ValueQuery(QueryHandler.Conv, Converged ? 1.0 : 0.0, true);
                    G.queryHandler.ValueQuery(QueryHandler.NoIter, NoOfIterations, true);
                    G.queryHandler.ValueQuery(QueryHandler.NoOfCells, G.gdat.CellPartitioning.TotalLength, true);
                    G.queryHandler.ValueQuery(QueryHandler.DOFs, TotalDOF, true); // 'essential' DOF, in the XDG case less than cordinate mapping length 

                    if(Solution.MaxTotalNoOfCoordinatesPerCell == Solution.MinTotalNoOfCoordinatesPerCell)
                        G.queryHandler.ValueQuery("BlockSize", Solution.MaxTotalNoOfCoordinatesPerCell, true);
                    else
                        G.queryHandler.ValueQuery("BlockSize", -99, true);
                    G.queryHandler.ValueQuery("maxBlkSize", Solution.MaxTotalNoOfCoordinatesPerCell, true);
                    G.queryHandler.ValueQuery("minBlkSize", Solution.MinTotalNoOfCoordinatesPerCell, true);
                }

                return Converged;
            }
        }

        //public static BlockMsrMatrix LastMtx;

        /// <summary>
        /// Easy-to-use driver routine for operator analysis
        /// </summary>
        /// <param name="op"></param>
        /// <param name="Mapping">
        /// Current Solution resp. solution approximation;
        /// Any matrix analysis is performed at this linearization point.
        /// </param>
        /// <param name="MgConfig">
        /// provisional; not that the block preconditioning configured here has huge effects on the condition number.
        /// </param>
        /// <returns></returns>
        /// <seealso cref="BoSSS.Solution.Application{T}.OperatorAnalysis"/>
        static public IDictionary<string, double> OperatorAnalysis(this ISpatialOperator op, CoordinateMapping Mapping, MultigridOperator.ChangeOfBasisConfig[][] MgConfig) {

            var G = new MatrixAssembler(op, Mapping, null, null, MgConfig);

            G.AssembleMatrix(out var Op_Matrix, out var Op_Affine, out var MassMatrix, G.SolutionFields, true, out _);

            OpAnalysisBase ana;
            if(G.LsTrk != null) {
                // ++++++++++
                // XDG branch
                // ++++++++++

                ana = new OpAnalysisBase(G.LsTrk,
                   Op_Matrix, Op_Affine,
                   Mapping, G.Op_Agglomeration,
                   MassMatrix,
                   MgConfig, op);

            } else {
                // ++++++++++
                // DG branch
                // ++++++++++

                ana = new OpAnalysisBase(
                   Op_Matrix, Op_Affine,
                   Mapping,
                   MgConfig, op);
            }

            //long J = G.gdat.CellPartitioning.TotalLength;
            //ana.PrecondOpMatrix.SaveToTextFileSparse("OpMatrix-J" + J + ".txt");

            return ana.GetNamedProperties();
        }


        /// <summary>
        /// Easy-to-use driver routine to obtain the multigrid-operator object
        /// </summary>
        /// <param name="op"></param>
        /// <param name="Mapping">
        /// Current Solution resp. solution approximation;
        /// Any matrix analysis is performed at this linearization point.
        /// </param>
        /// <param name="MgConfig">
        /// provisional; note that the block preconditioning configured here has huge effects on the condition number.
        /// </param>
        /// <returns></returns>
        static public MultigridOperator GetMultigridOperator(this ISpatialOperator op, CoordinateMapping Mapping, MultigridOperator.ChangeOfBasisConfig[][] MgConfig) {

            var G = new MatrixAssembler(op, Mapping, null, null, MgConfig);

            G.AssembleMatrix(out var Op_Matrix, out var _, out var MassMatrix, G.SolutionFields, true, out _);

            var MultigridOp = new MultigridOperator(G.AggBasisS, Mapping,
                        Op_Matrix, MassMatrix, G.MgConfig,
                        op);

            return MultigridOp;

        }

        /// <summary>
        /// Easy-to-use driver routine to obtain the block-preconditioned operator matrix.
        /// </summary>
        /// <param name="op"></param>
        /// <param name="Mapping">
        /// Current Solution resp. solution approximation;
        /// Any matrix analysis is performed at this linearization point.
        /// </param>
        /// <param name="MgConfig">
        /// provisional; note that the block preconditioning configured here has huge effects on the condition number.
        /// </param>
        /// <returns></returns>
        static public BlockMsrMatrix GetMatrix(this ISpatialOperator op, CoordinateMapping Mapping, MultigridOperator.ChangeOfBasisConfig[][] MgConfig) {

            return GetMultigridOperator(op, Mapping, MgConfig).OperatorMatrix;

        }

    }
}
