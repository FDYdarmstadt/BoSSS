// using BoSSS.Foundation;
// using BoSSS.Foundation.Grid;
// using BoSSS.Foundation.Grid.Aggregation;
// using BoSSS.Foundation.Grid.Classic;
// using BoSSS.Foundation.XDG;
// using BoSSS.Solution.AdvancedSolvers;
// using BoSSS.Solution.Control;
// using BoSSS.Solution.Queries;
// using BoSSS.Solution.Timestepping;
// using BoSSS.Solution.Gnuplot;
// using ilPSP;
// using ilPSP.LinSolvers;
// using ilPSP.Tracing;
// using ilPSP.Utils;
// using MPI.Wrappers;
// using NUnit.Framework;
// using System;
// using System.Collections.Generic;
// using System.ComponentModel;
// using System.Diagnostics;
// using System.Linq;
// using System.Runtime.Serialization;
// using System.Text;
// using System.Threading.Tasks;

// namespace BoSSS.Solution.XdgTimestepping {

//     public class XdgBDFSubGridTimestepping : XdgBDFTimestepping {

//         /// <summary>
//         /// Constructor for an XDG operator only evaluating on a subgrid (see <see cref="XdgOperator"/>)
//         /// </summary>
//         public XdgBDFSubGridTimestepping(
//             IEnumerable<DGField> Fields,
//             IEnumerable<DGField> __Parameters,
//             IEnumerable<DGField> IterationResiduals,
//             LevelSetTracker LsTrk,
//             bool DelayInit,
//             DelComputeOperatorMatrix _ComputeOperatorMatrix,
//             ISpatialOperator abstractOperator,
//             Func<ISlaveTimeIntegrator> _UpdateLevelset,
//             int BDForder,
//             LevelSetHandling _LevelSetHandling,
//             MassMatrixShapeandDependence _MassMatrixShapeandDependence,
//             SpatialOperatorType _SpatialOperatorType,
//             MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
//             AggregationGridData[] _MultigridSequence,
//             SpeciesId[] _SpId,
//             int _CutCellQuadOrder,
//             double _AgglomerationThreshold, bool _useX, Control.NonLinearSolverConfig nonlinconfig,
//             ISolverFactory linearconfig)
//          : base(
//             Fields,
//             __Parameters,
//             IterationResiduals,
//             LsTrk,
//             DelayInit,
//             _ComputeOperatorMatrix,
//             abstractOperator,
//             _UpdateLevelset,
//             BDForder,
//             _LevelSetHandling,
//             _MassMatrixShapeandDependence,
//             _SpatialOperatorType,
//             _MultigridOperatorConfig,
//             _MultigridSequence,
//             _SpId,
//             _CutCellQuadOrder,
//             _AgglomerationThreshold, _useX, nonlinconfig,
//             linearconfig
// ) //
//         {
//         }

//         /// <summary>
//         /// driver for solver calls
//         /// </summary>
//         /// <returns>
//         /// - true: solver algorithm successfully converged
//         /// - false: something went wrong
//         /// </returns>
//         public override bool Solve(double phystime, double dt) {

//             // bool SkipSolveAndEvaluateResidual = false;
//             if(!initialized)
//                 SingleInit();


//             if (dt <= 0)
//                 throw new ArgumentOutOfRangeException();
//             //if (m_CurrentDt_Timestep > 0 && Math.Abs(dt / m_CurrentDt_Timestep - 1.0) > 1.0e-14)
//             //    throw new ArgumentOutOfRangeException();
//             if (m_CurrentDt_Timestep > 0 && Math.Abs(m_CurrentDt_Timestep - dt) > 1e-14)
//                 AdaptToNewTimestep(dt, m_CurrentDt_Timestep);

//             m_CurrentDt_Timestep = dt;

//             bool success = true;
//             //for (int i = 1; i <= incrementTimesteps; i++) {
//             {
//                 // push levelsets for every incremental timestep
//                 //if (i > 1)
//                 //    PushLevelSet();

//                 //// solve timestep with incremental timestep size
//                 //double incTimestepSize = dt / (double)incrementTimesteps;

//                 //success = success && Solve_Increment(i, phystime, incTimestepSize, ComputeOnlyResidual);
//                 success = success && Solve_Increment(phystime, dt, false);

//                 //phystime += incTimestepSize;
//                 phystime += dt;
//             }

//             return success;
//         }

// //         protected override bool Solve_Increment(double phystime, double dt, bool ComputeOnlyResidual = false) {
// //             using(var tr = new FuncTrace()) {
// //                 if(dt <= 0)
// //                     throw new ArgumentOutOfRangeException();
// //                 //if (m_CurrentDt > 0 && Math.Abs(dt / m_CurrentDt - 1.0) > 1.0e-14)
// //                 //    throw new ArgumentOutOfRangeException();

// //                 m_CurrentPhystime = phystime;
// //                 m_CurrentDt = dt;
// //                 m_IterationCounter = 0;
// //                 m_CoupledIterations = 0;
// //                 m_InnerCoupledIterations = 0;
// //                 PushStack();
// //                 fsiOldPhystime = phystime;
// //                 //if (incrementTimesteps == 1)
// //                 //    dt = m_CurrentDt;
// //                 //else
// //                 //    Console.WriteLine("Increment solve, timestep #{0}, dt = {1} ...", increment, dt);
// //                 dt = m_CurrentDt;


// //                 // ===========================================
// //                 // update level-set (in the case of splitting)
// //                 // ===========================================
// //                 if(this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
// //                     || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
// //                     || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {

// //                     Debug.Assert(m_CurrentAgglomeration == null);

// //                     double ls_dt = dt;
// //                     if(this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting)
// //                         ls_dt *= 0.5;

// //                     // remember which old cells had values
// //                     //var oldCCM = this.UpdateCutCellMetrics();

// //                     // evolve the level set
// //                     if(Config_LevelSetHandling != LevelSetHandling.FSILieSplittingFullyCoupled && !this.coupledOperator) {
// //                         m_LsTrk.IncreaseHistoryLength(1);
// //                         m_LsTrk.PushStacks();
// //                     }

// //                     int oldPushCount = m_LsTrk.PushCount;
// //                     int oldVersion = m_LsTrk.VersionCnt;
// //                     this.MoveLevelSetAndRelatedStuff(m_Stack_u[0].Mapping.Fields.ToArray(), phystime, ls_dt, 1.0);

// //                     int newPushCount = m_LsTrk.PushCount;
// //                     int newVersion = m_LsTrk.VersionCnt;
// //                     if((newPushCount - oldPushCount) != 0 && !coupledOperator)
// //                         throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");
// //                     if((newVersion - oldVersion) != 1 && !coupledOperator)
// //                     {
// //                         throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");
// //                     }

// //                     // in the case of splitting, the fields must be extrapolated
// //                     //var newCCM = this.UpdateCutCellMetrics();
// //                     //var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
// //                     Debug.Assert(m_LsTrk.HistoryLength >= 1);
// //                     var SplittingAgg = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
// //                         __AgglomerationTreshold: 0.0, AgglomerateNewborn: true, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true,
// //                         oldTs__AgglomerationTreshold: new double[] { 0.0 });
// //                     for(int i = 0; i < this.m_Stack_u.Length; i++)
// //                         SplittingAgg.Extrapolate(this.m_Stack_u[i].Mapping);

// //                     // delete new agglomeration; in case of splitting, the agglomeration for the **bulk operator timestep** does not depend on previous time-steps
// //                     m_CurrentAgglomeration = null;
// //                 }

// //                 // ==============================================
// //                 // solve main system
// //                 // ==============================================
// //                 bool success;

// //                 int oldLsTrkPushCount = m_LsTrk.PushCount;

// //                 {
// //                     int[] Jtot =
// //                         (new int[] { base.m_LsTrk.Regions.GetCutCellMask().NoOfItemsLocally, base.m_LsTrk.GridDat.Cells.NoOfLocalUpdatedCells })
// //                         .MPISum();
// //                     //Console.WriteLine("No of cells {0}, No of cut cells {1}.", Jtot[1], Jtot[0]);
// //                     // Console.WriteLine("WARNING: Comment this back in!");
// //                     if(Jtot[0] == Jtot[1])
// //                         throw new ArithmeticException("All cells are cut cells - check your settings!");
// //                 }

// //                 //{
// //                 //    var rnd = new Random(1);
// //                 //    var vec = new CoordinateVector(CurrentStateMapping);
// //                 //    int L = CurrentStateMapping.LocalLength;
// //                 //    for(int i = 0; i < L; i++)
// //                 //        vec[i] = rnd.NextDouble();


// //                 //    double[] Affine;
// //                 //    BoSSS.Foundation.Quadrature.NonLin.Arsch.ShutTheFuckUp = true;
// //                 //    this.AssembleMatrixCallback(out BlockMsrMatrix System, out Affine, out BlockMsrMatrix MaMa, CurrentStateMapping.Fields.ToArray(), false, out var dummy);
// //                 //    Debug.Assert(System == null);

// //                 //    base.Residuals.Clear();
// //                 //    base.Residuals.SetV(Affine, -1.0);


// //                 //}




// //                 if(!ComputeOnlyResidual) {

// //                     // ++++++++++++++++++
// //                     // normal solver run
// //                     // ++++++++++++++++++


// //                     if(RequiresNonlinearSolver) {
// //                         Console.WriteLine("Hello from Requires NLSolver");
// //                         // ++++++++++++++++++++++++++++++++
// //                         // Nonlinear Solver (Navier-Stokes)
// //                         // ++++++++++++++++++++++++++++++++

// //                         // use solver
// //                         var nonlinSolver = GetNonlinSolver();
// //                         Debugger.Launch();
// //                         Console.ReadLine();
// //                         success = nonlinSolver.SolverDriver(m_Stack_u[0], default(double[])); // Note: the RHS is passed as the affine part via 'this.SolverCallback'

// //                         // 'revert' agglomeration
// //                         Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
// //                         m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

// //                     } else {
// //                         Console.WriteLine("Hello from LinSolver");
// //                         // ++++++++++++++++++++++++++++++++
// //                         // Linear Solver (Stokes)
// //                         // ++++++++++++++++++++++++++++++++
// //                         tr.Info("Using linear solver.");

// //                         // build the saddle-point matrix
// //                         BlockMsrMatrix System, MaMa;
// //                         double[] RHS;
// //                         this.AssembleMatrixCallback(out System, out RHS, out MaMa, CurrentStateMapping.Fields.ToArray(), true, out var dummy);
// //                         RHS.ScaleV(-1);

// //                         // update the multigrid operator
// //                         csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
// //                         MultigridOperator mgOperator;
// //                         using(new BlockTrace("MultigridOperator setup", tr)) {
// //                             mgOperator = new MultigridOperator(this.MultigridBasis, CurrentStateMapping,
// //                                 System, MaMa,
// //                                 this.Config_MultigridOperator,
// //                                 dummy);
// //                         }

// //                         using (var linearSolver = GetLinearSolver(mgOperator)) {



// //                             // try to solve the saddle-point system.
// //                             TimeSpan duration;
// //                             using(new BlockTrace("Solver_Run", tr)) {
// //                                 var st = DateTime.Now;
// //                                 mgOperator.UseSolver(linearSolver, m_Stack_u[0], RHS);
// //                                 //mgOperator.ComputeResidual(this.Residuals, m_Stack_u[0], RHS);
// //                                 duration = DateTime.Now - st;
// //                             }
// //                             Console.WriteLine("solver success: " + linearSolver.Converged + "; runtime: " + duration.TotalSeconds + " sec.");
// //                             success = linearSolver.Converged;



// //                             // 'revert' agglomeration
// //                             Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
// //                             m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);


// //                             if(base.QueryHandler != null) {
// //                                 base.QueryHandler.ValueQuery(QueryHandler.Conv, linearSolver.Converged ? 1.0 : 0.0, true);
// //                                 base.QueryHandler.ValueQuery(QueryHandler.NoIter, linearSolver.ThisLevelIterations, true);
// //                                 base.QueryHandler.ValueQuery(QueryHandler.NoOfCells, this.m_LsTrk.GridDat.CellPartitioning.TotalLength, true);
// //                                 base.QueryHandler.ValueQuery(QueryHandler.DOFs, mgOperator.Mapping.TotalLength, true); // 'essential' DOF, in the XDG case less than cordinate mapping length

// //                             }

// //                         }

// //                     //ExtractSomeSamplepoints("samples");
// //                 }

// //             } else {
// //                 // ++++++++++++++++++++++++++++++++++++
// //                 // compute residual of actual solution
// //                 // ++++++++++++++++++++++++++++++++++++


// //                     double[] Affine;
// //                     this.AssembleMatrixCallback(out BlockMsrMatrix System, out Affine, out BlockMsrMatrix MaMa, CurrentStateMapping.Fields.ToArray(), false, out var dummy);
// //                     Debug.Assert(System == null);

// //                     base.Residuals.Clear();
// //                     base.Residuals.SetV(Affine, -1.0);

// //                     success = true;

// // #if DEBUG
// //                     {

// //                         this.AssembleMatrixCallback(out BlockMsrMatrix checkSystem, out double[] checkAffine, out BlockMsrMatrix MaMa1, CurrentStateMapping.Fields.ToArray(), true, out var dummy2);

// //                         double[] checkResidual = new double[checkAffine.Length];
// //                         checkResidual.SetV(checkAffine, -1.0);
// //                         checkSystem.SpMV(-1.0, m_Stack_u[0], +1.0, checkResidual);

// //                         Console.WriteLine("Norm of evaluated residual: " + base.Residuals.MPI_L2Norm());
// //                         Console.WriteLine("Norm of reference residual: " + checkResidual.MPI_L2Norm());


// //                         double distL2 = GenericBlas.L2DistPow2(checkResidual, base.Residuals).MPISum().Sqrt();
// //                         double refL2 = (new double[] { GenericBlas.L2NormPow2(m_Stack_u[0]), GenericBlas.L2NormPow2(checkResidual), GenericBlas.L2NormPow2(base.Residuals) }).MPISum().Max().Sqrt();

// //                         if(distL2 >= refL2 * 1.0e-5) {
// //                             double __distL2 = GenericBlas.L2DistPow2(checkAffine, base.Residuals).MPISum().Sqrt();
// //                         }

// //                         Tecplot.Tecplot.PlotFields(base.Residuals.Fields, "resi", 0.0, 2);

// //                         Assert.LessOrEqual(distL2, refL2 * 1.0e-5, "Significant difference between linearized and non-linear evaluation.");

// //                     }
// // #endif


// //                     var ResidualFields = base.Residuals.Mapping.Fields.ToArray();


// //                     for(int i = 0; i < ResidualFields.Length; i++) {
// //                         double L2Res = ResidualFields[i].L2Norm();
// //                         if(m_ResLogger != null) {
// //                             m_ResLogger.CustomValue(L2Res, m_ResidualNames[i]);
// //                         } else {
// //                             Console.WriteLine("Residual {0}: {1}", m_ResidualNames != null ? m_ResidualNames[i] : ResidualFields[i].Identification, L2Res);
// //                         }
// //                     }

// //                     if(Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && m_ResLogger != null) {
// //                         m_ResLogger.CustomValue(0.0, "LevelSet");
// //                     }

// //                     if(m_ResLogger != null)
// //                         m_ResLogger.NextIteration(true);
// //                 }

// //                 int newLsTrkPushCount = m_LsTrk.PushCount;
// //                 if(newLsTrkPushCount != oldLsTrkPushCount)
// //                     throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");

// //                 if(Config_LevelSetHandling != LevelSetHandling.None)
// //                     m_CurrentAgglomeration = null;

// //                 // ===========================================
// //                 // update level-set (in the case of splitting)
// //                 // ===========================================

// //                 if(this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting) {

// //                     Debug.Assert(m_CurrentAgglomeration == null);

// //                     // evolve the level set
// //                     m_LsTrk.IncreaseHistoryLength(1);
// //                     m_LsTrk.PushStacks();

// //                     int oldPushCount = m_LsTrk.PushCount;
// //                     int oldVersion = m_LsTrk.VersionCnt;

// //                     this.MoveLevelSetAndRelatedStuff(m_Stack_u[0].Mapping.Fields.ToArray(), phystime + dt * 0.5, dt * 0.5, 1.0);

// //                     int newPushCount = m_LsTrk.PushCount;
// //                     int newVersion = m_LsTrk.VersionCnt;
// //                     if((newPushCount - oldPushCount) != 0)
// //                         throw new ApplicationException("Calling 'LevelSetTracker.PushStacks()' is not allowed. Level-set-tracker stacks must be controlled by time-stepper.");
// //                     if((newVersion - oldVersion) != 1)
// //                         throw new ApplicationException("Expecting exactly one call to 'UpdateTracker(...)' in 'UpdateLevelset(...)'.");

// //                     // in the case of splitting, the fields must be extrapolated
// //                     //var newCCM = this.UpdateCutCellMetrics();
// //                     //var SplittingAgg = new MultiphaseCellAgglomerator(newCCM, 0.0, true, false, true, new CutCellMetrics[] { oldCCM }, new double[] { 0.0 });
// //                     //SplittingAgg.Extrapolate(this.CurrentStateMapping);
// //                     Debug.Assert(m_LsTrk.HistoryLength >= 1);
// //                     var SplittingAgg = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
// //                         __AgglomerationTreshold: 0.0, AgglomerateNewborn: true, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true,
// //                         oldTs__AgglomerationTreshold: new double[] { 0.0 });
// //                     for(int i = 0; i < this.m_Stack_u.Length; i++)
// //                         SplittingAgg.Extrapolate(this.m_Stack_u[i].Mapping);

// //                     Debug.Assert(m_CurrentAgglomeration == null);
// //                 }

// //                 // ======
// //                 // return
// //                 // ======
// //                 //string path = Directory.GetCurrentDirectory();
// //                 //var dinfo = Directory.CreateDirectory(path+@"\plots");
// //                 //ExecuteWaterfallAnalysis(dinfo.FullName);
// //                 //CreateFAMatrices(dinfo.FullName);


// //                 m_CurrentPhystime = phystime + dt;

// //                 if(Config_LevelSetHandling == LevelSetHandling.None) {
// //                     m_LsTrk.UpdateTracker(m_CurrentPhystime); // call is required to bring the internal time-stamp up-to-date;
// //                 }

// //                 return success;
// //             }
// //         }
//     }
// }
