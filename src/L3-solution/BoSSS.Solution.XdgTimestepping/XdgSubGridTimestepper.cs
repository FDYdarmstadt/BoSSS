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

//     public class XdgSubGridTimestepping : XdgTimestepping {

//         SubGrid m_SubGrid;
//         SubGridBoundaryModes m_SubGridBnd;

//         XdgBDFSubGridTimestepping m_Timestepper;

//         /// <summary>
//         /// Constructor for an XDG operator only evaluating on a subgrid (see <see cref="XdgOperator"/>)
//         /// </summary>
//         public XdgSubGridTimestepping(
//             XSpatialOperatorMk2 op,
//             IEnumerable<DGField> Fields,
//             IEnumerable<DGField> IterationResiduals,
//             TimeSteppingScheme __Scheme,
//             SubGrid subGrid,
//             SubGridBoundaryModes subGridBnd,
//             Func<ISlaveTimeIntegrator> _UpdateLevelset = null,
//             LevelSetHandling _LevelSetHandling = LevelSetHandling.None,
//             MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig = null,
//             AggregationGridData[] _MultigridSequence = null,
//             double _AgglomerationThreshold = 0.1,
//             AdvancedSolvers.ISolverFactory LinearSolver = null, NonLinearSolverConfig NonLinearSolver = null,
//             LevelSetTracker _optTracker = null,
//             IList<DGField> _Parameters = null,
//             QueryHandler queryHandler = null) //
//          : base(
//             op,
//             Fields,
//             IterationResiduals,
//             __Scheme,
//             _UpdateLevelset,
//             _LevelSetHandling,
//             _MultigridOperatorConfig,
//             _MultigridSequence,
//             _AgglomerationThreshold,
//             LinearSolver, NonLinearSolver,
//             _optTracker,
//             _Parameters,
//             queryHandler) //
//         {
//             m_SubGrid = subGrid;
//             m_SubGridBnd = subGridBnd;
//             throw new NotImplementedException("This feature has not been properly implemented, please expect to fix some bugs before it works");
//         }

//         /// <summary>
//         /// Constructor for conventional (non-X, but potentially multiphase) DG
//         /// </summary>
//         public XdgSubGridTimestepping(
//             SpatialOperator op,
//             IEnumerable<DGField> Fields,
//             IEnumerable<DGField> IterationResiduals,
//             TimeSteppingScheme __Scheme,
//             SubGrid subGrid,
//             SubGridBoundaryModes subGridBnd,
//             LevelSetTracker _optTracker = null,
//             Func<ISlaveTimeIntegrator> _UpdateLevelset = null,
//             LevelSetHandling _LevelSetHandling = LevelSetHandling.None,
//             MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig = null,
//             double _AgglomerationThreshold = 0.1,
//             AggregationGridData[] _MultigridSequence = null,
//             ISolverFactory LinearSolver = null, NonLinearSolverConfig NonLinearSolver = null,
//             IList<DGField> _Parameters = null,
//             QueryHandler queryHandler = null) //
//             : base(
//                 op,
//                 Fields,
//                 IterationResiduals,
//                 __Scheme,
//                 _optTracker,
//                 _UpdateLevelset,
//                 _LevelSetHandling,
//                 _MultigridOperatorConfig,
//                 _MultigridSequence,
//                 LinearSolver, NonLinearSolver,
//                 _Parameters,
//                 _AgglomerationThreshold,
//                 queryHandler)
//         {
//             m_SubGrid = subGrid;
//             m_SubGridBnd = subGridBnd;

//             int bdfOrder;
//             RungeKuttaScheme rksch;
//             DecodeScheme(this.Scheme, out rksch, out bdfOrder);
//             SpatialOperatorType _SpatialOperatorType = op.IsLinear ? SpatialOperatorType.LinearTimeDependent : SpatialOperatorType.Nonlinear;

//             SpeciesId[] spcToCompute = _optTracker.SpeciesNames.Select(spcName => LsTrk.GetSpeciesId(spcName)).ToArray();

//             int quadOrder = op.QuadOrderFunction(
//                 Fields.Select(f => f.Basis.Degree).ToArray(),
//                 Parameters.Select(f => f != null ? f.Basis.Degree : 0).ToArray(),
//                 IterationResiduals.Select(f => f.Basis.Degree).ToArray());

//             if(_MultigridOperatorConfig == null) {
//                 int NoOfVar = Fields.Count();
//                 _MultigridOperatorConfig = new MultigridOperator.ChangeOfBasisConfig[1][];
//                 _MultigridOperatorConfig[0] = new MultigridOperator.ChangeOfBasisConfig[NoOfVar];
//                 for(int iVar = 0; iVar < NoOfVar; iVar++) {
//                     _MultigridOperatorConfig[0][iVar] = new MultigridOperator.ChangeOfBasisConfig() {
//                         DegreeS = new int[] { Fields.ElementAt(iVar).Basis.Degree },
//                         mode = MultigridOperator.Mode.Eye,
//                         VarIndex = new int[] { iVar }
//                     };
//                 }

//             }

//             if (_MultigridSequence == null) {
//                 _MultigridSequence = new[] { CoarseningAlgorithms.ZeroAggregation(this.GridDat) };
//             }

//             m_Timestepper = new XdgBDFSubGridTimestepping(Fields, _Parameters, IterationResiduals,
//                 LsTrk, true,
//                 this.ComputeOperatorMatrix, op, _UpdateLevelset,
//                 bdfOrder,
//                 _LevelSetHandling,
//                 MassMatrixShapeandDependence.IsTimeDependent,
//                 _SpatialOperatorType,
//                 _MultigridOperatorConfig, _MultigridSequence,
//                 spcToCompute, quadOrder,
//                 _AgglomerationThreshold, false,
//                 NonLinearSolver,
//                 LinearSolver);

//             m_BDF_Timestepper.Config_AgglomerationThreshold = _AgglomerationThreshold;
//         }

//         /// <summary>
//         /// Operator Evaluation and Linearization, <see cref="DelComputeOperatorMatrix"/>:
//         /// - either update operator linearization matrix
//         /// - or evaluate the operator in the current linearization point
//         /// In both cases, only the spatial component (i.e. no temporal derivatives) are linearized/evaluated.
//         /// /// </summary>
//         public override void ComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] __CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time, int LsTrkHistoryIndex) {
//             using(var ft = new FuncTrace()) {
//                 // compute operator
//                 Debug.Assert(OpAffine.L2Norm() == 0.0);

//                 // all kinds of checks
//                 if(!this.CurrentState.EqualsPartition(Mapping))
//                     throw new ApplicationException("something is weird");
//                 if(OpMtx != null) {
//                     if(!OpMtx.RowPartitioning.EqualsPartition(Mapping))
//                         throw new ArgumentException("Codomain/Matrix Row mapping mismatch.");
//                     if(!OpMtx.ColPartition.EqualsPartition(Mapping))
//                         throw new ArgumentException("Domain/Matrix column mapping mismatch.");
//                 }

//                 if(XdgOperator != null) {
//                     throw new NotImplementedException();
//                     // // +++++++++++++++++++++++++++++++++++++++++++++++
//                     // // XDG Branch: still requires length-scale-hack
//                     // // (should be cleaned some-when in the future)
//                     // // +++++++++++++++++++++++++++++++++++++++++++++++

//                     // //if(XdgOperator.AgglomerationThreshold <= 0)
//                     // //    throw new ArgumentException("Mismatch between agglomeration threshold provided ");

//                     // // throw new NotImplementedException();
//                     // if(OpMtx != null) {
//                     //     // +++++++++++++++++++++++++++++
//                     //     // Solver requires linearization
//                     //     // +++++++++++++++++++++++++++++

//                     //     Debug.Assert(OpMtx.InfNorm() == 0.0);
//                     //     switch(XdgOperator.LinearizationHint) {

//                     //         case LinearizationHint.AdHoc: using(new BlockTrace("XDG-LinearizationHint.AdHoc", ft, false)) {
//                     //             throw new NotImplementedException();
//                     //             // this.XdgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());



//                     //             // var mtxBuilder = XdgOperator.GetMatrixBuilder(LsTrk, Mapping, this.Parameters, Mapping, LsTrkHistoryIndex);
//                     //             // mtxBuilder.time = time;
//                     //             // mtxBuilder.MPITtransceive = true;
//                     //             // foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
//                     //             //     mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
//                     //             // }
//                     //             // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

//                     //             // return;
//                     //         }

//                     //         case LinearizationHint.FDJacobi: using(new BlockTrace("XDG-LinearizationHint.FDJacobi", ft, false)){
//                     //             throw new NotImplementedException();
//                     //             // var mtxBuilder = XdgOperator.GetFDJacobianBuilder(LsTrk, __CurrentState, this.Parameters, Mapping, LsTrkHistoryIndex);
//                     //             // mtxBuilder.time = time;
//                     //             // mtxBuilder.MPITtransceive = true;
//                     //             // if(mtxBuilder.Eval is XSpatialOperatorMk2.XEvaluatorNonlin evn) { // length-scale hack
//                     //             //     foreach(var kv in AgglomeratedCellLengthScales) {
//                     //             //         evn.CellLengthScales[kv.Key] = kv.Value;
//                     //             //     }
//                     //             // }
//                     //             // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
//                     //             // return;
//                     //         }

//                     //         case LinearizationHint.GetJacobiOperator: using(new BlockTrace("XDG-LinearizationHint.GetJacobiOperator", ft, false)){
//                     //             var op = GetJacobiXdgOperator();

//                     //             if(JacobiParameterVars == null)
//                     //                 JacobiParameterVars = op.InvokeParameterFactory(this.CurrentState);

//                     //             op.InvokeParameterUpdate(time, __CurrentState, JacobiParameterVars);

//                     //             var mtxBuilder = op.GetMatrixBuilder(LsTrk, Mapping, this.JacobiParameterVars, Mapping, LsTrkHistoryIndex);
//                     //             mtxBuilder.time = time;
//                     //             mtxBuilder.MPITtransceive = true;
//                     //             foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
//                     //                 mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
//                     //             }
//                     //             mtxBuilder.ActivateSubgridBoundary(this.m_SubGrid.VolumeMask, m_SubGridBnd);
//                     //             mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
//                     //             return;
//                     //         }
//                         // }
//                     // } else {
//                         // ++++++++++++++++++++++++
//                         // only operator evaluation
//                         // ++++++++++++++++++++++++


//                         // throw new NotImplementedException();

//                         // using(new BlockTrace("XDG-Evaluate", ft, false)) {
//                         //     this.XdgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

//                         //     var eval = XdgOperator.GetEvaluatorEx(this.LsTrk, __CurrentState, this.Parameters, Mapping, LsTrkHistoryIndex);
//                         //     eval.time = time;

//                         //     eval.MPITtransceive = true;
//                         //     foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
//                         //         eval.CellLengthScales[kv.Key] = kv.Value;
//                         //     }
//                         //     eval.Evaluate(1.0, 0.0, OpAffine);
//                 } else if(DgOperator != null) {
//                     // +++++++++++++++++++++++++++++++++++++++++++++++
//                     // DG Branch
//                     // +++++++++++++++++++++++++++++++++++++++++++++++

//                     // throw new NotImplementedException();
//                     if(OpMtx != null) {
//                         // +++++++++++++++++++++++++++++
//                         // Solver requires linearization
//                         // +++++++++++++++++++++++++++++

//                         Debug.Assert(OpMtx.InfNorm() == 0.0);
//                         switch(DgOperator.LinearizationHint) {

//                             case LinearizationHint.AdHoc: using(new BlockTrace("DG-LinearizationHint.AdHoc", ft)) {
//                                 throw new NotImplementedException();
//                                 // this.DgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());
//                                 // var mtxBuilder = DgOperator.GetMatrixBuilder(Mapping, this.Parameters, Mapping);
//                                 // mtxBuilder.time = time;
//                                 // mtxBuilder.MPITtransceive = true;
//                                 // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
//                                 // return;
//                             }

//                             case LinearizationHint.FDJacobi: using(new BlockTrace("DG-LinearizationHint.FDJacobi", ft)){
//                                 throw new NotImplementedException();
//                                 // var mtxBuilder = DgOperator.GetFDJacobianBuilder(__CurrentState, this.Parameters, Mapping);
//                                 // mtxBuilder.time = time;
//                                 // mtxBuilder.MPITtransceive = true;
//                                 // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
//                                 // return;
//                             }

//                             case LinearizationHint.GetJacobiOperator: using(new BlockTrace("DG-LinearizationHint.GetJacobiOperator", ft)){
//                                 var op = GetJacobiDgOperator();

//                                 if(JacobiParameterVars == null)
//                                     JacobiParameterVars = op.InvokeParameterFactory(__CurrentState);
//                                 op.InvokeParameterUpdate(time, __CurrentState, JacobiParameterVars);

//                                 var mtxBuilder = op.GetMatrixBuilder(Mapping, this.JacobiParameterVars, Mapping);
//                                 if (this.m_SubGrid != null) {
//                                     mtxBuilder.ActivateSubgridBoundary(this.m_SubGrid.VolumeMask, m_SubGridBnd);
//                                 }
//                                 mtxBuilder.time = time;
//                                 mtxBuilder.MPITtransceive = true;
//                                 mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

//                                 // for (int i = 0; i < OpMtx.NoOfRows; i++){
//                                 //     for (int j = 0; j < OpMtx.NoOfRows; j++){
//                                 //         Console.Write(OpMtx[i,j] == 0.0 ? "o" : "x");
//                                 //     }
//                                 //     Console.WriteLine();
//                                 // }

//                                 // PlotFormat format = new PlotFormat(
//                                 //         Style: Styles.Points,
//                                 //         pointType: PointTypes.Asterisk,
//                                 //         pointSize: 0.5);

//                                 // var gp = new Gnuplot.Gnuplot(baseLineFormat: format);
//                                 // gp.PlotMatrixStructure(OpMtx.ToFullMatrixOnProc0());
//                                 // gp.PlotDataFile("/home/klingenberg/Downloads/matrix.png", deferred: false);
//                                 // // gp.PlotDataFile("matrix.png", deferred: false);
//                                 // gp.Execute();
//                                 // Console.ReadLine();

//                                 for (int i = 0; i < OpMtx.NoOfRows; i++)
//                                 {
//                                     if (!(OpMtx.GetNoOfNonZerosPerRow(i) > 0))
//                                     {
//                                         Console.WriteLine("changing matrix in line " + i);
//                                         OpMtx.SetDiagonalElement(i, 1.0);
//                                     }
//                                 }
//                                 return;
//                             }
//                         }
//                     } else {
//                         // ++++++++++++++++++++++++
//                         // only operator evaluation
//                         // ++++++++++++++++++++++++

//                         // throw new NotImplementedException();
//                         using(new BlockTrace("DG-Evaluate", ft)) {
//                             this.DgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

//                             var eval = DgOperator.GetEvaluatorEx(__CurrentState, this.Parameters, Mapping);
//                             eval.time = time;
//                             eval.MPITtransceive = true;
//                             eval.Evaluate(1.0, 0.0, OpAffine);
//                         }
//                     }
//                 } else {
//                     throw new NotImplementedException();
//                 }
//             }
//         }

//         /// <summary>
//         /// driver for solver calls
//         /// </summary>
//         /// <returns>
//         /// - true: solver algorithm successfully converged
//         /// - false: something went wrong
//         /// </returns>
//         public override bool Solve(double phystime, double dt, bool SkipSolveAndEvaluateResidual = false) {
//             bool success = false;
//             if(m_Timestepper == null)
//                 throw new ApplicationException();
//             if(!ilPSP.DoubleExtensions.ApproxEqual(this.LsTrk.Regions.Time, phystime))
//                 throw new ApplicationException($"Before timestep, mismatch in time between tracker (Regions.Time = {LsTrk.Regions.Time}) and physical time ({phystime})");



//             double[] AvailTimesBefore;
//             if(TimesteppingBase.Config_LevelSetHandling != LevelSetHandling.None) {
//                 AvailTimesBefore = LsTrk.TimeLevelsInStack;
//                 Assert.IsTrue((AvailTimesBefore[0] - phystime).Abs() < dt*1e-7, "Error in Level-Set tracker time");
//             }

//             if(UseAdaptiveTimestepping()) {

//                 TimeLevel TL = new TimeLevel(this, dt, StateAtTime.Obtain(this, phystime));
//                 TL.Compute();
//                 success = true;

//             } else {
//                     success = m_Timestepper.Solve(phystime, dt, SkipSolveAndEvaluateResidual);
//             }

//             double[] AvailTimesAfter;
//             if(TimesteppingBase.Config_LevelSetHandling != LevelSetHandling.None) {
//                 AvailTimesAfter = LsTrk.RegionsHistory.AvailableIndices.Select((int iHist) => LsTrk.RegionsHistory[iHist].Time).ToArray();
//                 if(!AvailTimesAfter[0].ApproxEqual(phystime + dt))
//                     throw new ApplicationException($"Internal algorithm inconsistency: Error in Level-Set tracker time; expecting time {phystime + dt}, but most recent tracker time is {AvailTimesAfter[0]}");
//             }

//             if(!ilPSP.DoubleExtensions.ApproxEqual(this.LsTrk.Regions.Time, phystime + dt))
//                 throw new ApplicationException($"After timestep, mismatch in time between tracker (Regions.Time = {LsTrk.Regions.Time}) and physical time ({phystime + dt}), level-set handling is '{TimesteppingBase.Config_LevelSetHandling}'.");


//             JacobiParameterVars = null;

//             return success;
//         }

//     }
// }
