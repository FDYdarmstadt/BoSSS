using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XdgTimestepping {

    public class XdgSubGridTimestepping : XdgTimestepping {

        SubGrid m_SubGrid;
        SubGridBoundaryModes m_SubGridBnd;

        /// <summary>
        /// Constructor for an XDG operator only evaluating on a subgrid (see <see cref="XdgOperator"/>)
        /// </summary>
        public XdgSubGridTimestepping(
            XSpatialOperatorMk2 op,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            TimeSteppingScheme __Scheme,
            SubGrid subGrid,
            SubGridBoundaryModes subGridBnd,
            Func<ISlaveTimeIntegrator> _UpdateLevelset = null,
            LevelSetHandling _LevelSetHandling = LevelSetHandling.None,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig = null,
            AggregationGridData[] _MultigridSequence = null,
            double _AgglomerationThreshold = 0.1,
            AdvancedSolvers.ISolverFactory LinearSolver = null, NonLinearSolverConfig NonLinearSolver = null,
            LevelSetTracker _optTracker = null,
            IList<DGField> _Parameters = null,
            QueryHandler queryHandler = null) //
         : base(
            op,
            Fields,
            IterationResiduals,
            __Scheme,
            _UpdateLevelset,
            _LevelSetHandling,
            _MultigridOperatorConfig,
            _MultigridSequence,
            _AgglomerationThreshold,
            LinearSolver, NonLinearSolver,
            _optTracker,
            _Parameters,
            queryHandler) //
        {
            m_SubGrid = subGrid;
            m_SubGridBnd = subGridBnd;
            throw new NotImplementedException("This feature has not been properly implemented, please expect to fix some bugs before it works");
        }

        /// <summary>
        /// Constructor for conventional (non-X, but potentially multiphase) DG
        /// </summary>
        public XdgSubGridTimestepping(
            SpatialOperator op,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            TimeSteppingScheme __Scheme,
            SubGrid subGrid,
            SubGridBoundaryModes subGridBnd,
            LevelSetTracker _optTracker = null,
            Func<ISlaveTimeIntegrator> _UpdateLevelset = null,
            LevelSetHandling _LevelSetHandling = LevelSetHandling.None,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig = null,
            double _AgglomerationThreshold = 0.1,
            AggregationGridData[] _MultigridSequence = null,
            ISolverFactory LinearSolver = null, NonLinearSolverConfig NonLinearSolver = null,
            IList<DGField> _Parameters = null,
            QueryHandler queryHandler = null) //
            : base(
                op,
                Fields,
                IterationResiduals,
                __Scheme,
                _optTracker,
                _UpdateLevelset,
                _LevelSetHandling,
                _MultigridOperatorConfig,
                _MultigridSequence,
                LinearSolver, NonLinearSolver,
                _Parameters,
                _AgglomerationThreshold,
                queryHandler)
        {
            m_SubGrid = subGrid;
            m_SubGridBnd = subGridBnd;
        }

        /// <summary>
        /// Operator Evaluation and Linearization, <see cref="DelComputeOperatorMatrix"/>:
        /// - either update operator linearization matrix
        /// - or evaluate the operator in the current linearization point
        /// In both cases, only the spatial component (i.e. no temporal derivatives) are linearized/evaluated.
        /// /// </summary>
        public override void ComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] __CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time, int LsTrkHistoryIndex) {
            using(var ft = new FuncTrace()) {
                // compute operator
                Debug.Assert(OpAffine.L2Norm() == 0.0);

                // all kinds of checks
                if(!this.CurrentState.EqualsPartition(Mapping))
                    throw new ApplicationException("something is weird");
                if(OpMtx != null) {
                    if(!OpMtx.RowPartitioning.EqualsPartition(Mapping))
                        throw new ArgumentException("Codomain/Matrix Row mapping mismatch.");
                    if(!OpMtx.ColPartition.EqualsPartition(Mapping))
                        throw new ArgumentException("Domain/Matrix column mapping mismatch.");
                }

                if(XdgOperator != null) {
                    throw new NotImplementedException();
                    // // +++++++++++++++++++++++++++++++++++++++++++++++
                    // // XDG Branch: still requires length-scale-hack
                    // // (should be cleaned some-when in the future)
                    // // +++++++++++++++++++++++++++++++++++++++++++++++

                    // //if(XdgOperator.AgglomerationThreshold <= 0)
                    // //    throw new ArgumentException("Mismatch between agglomeration threshold provided ");

                    // // throw new NotImplementedException();
                    // if(OpMtx != null) {
                    //     // +++++++++++++++++++++++++++++
                    //     // Solver requires linearization
                    //     // +++++++++++++++++++++++++++++

                    //     Debug.Assert(OpMtx.InfNorm() == 0.0);
                    //     switch(XdgOperator.LinearizationHint) {

                    //         case LinearizationHint.AdHoc: using(new BlockTrace("XDG-LinearizationHint.AdHoc", ft, false)) {
                    //             throw new NotImplementedException();
                    //             // this.XdgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());



                    //             // var mtxBuilder = XdgOperator.GetMatrixBuilder(LsTrk, Mapping, this.Parameters, Mapping, LsTrkHistoryIndex);
                    //             // mtxBuilder.time = time;
                    //             // mtxBuilder.MPITtransceive = true;
                    //             // foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
                    //             //     mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                    //             // }
                    //             // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

                    //             // return;
                    //         }

                    //         case LinearizationHint.FDJacobi: using(new BlockTrace("XDG-LinearizationHint.FDJacobi", ft, false)){
                    //             throw new NotImplementedException();
                    //             // var mtxBuilder = XdgOperator.GetFDJacobianBuilder(LsTrk, __CurrentState, this.Parameters, Mapping, LsTrkHistoryIndex);
                    //             // mtxBuilder.time = time;
                    //             // mtxBuilder.MPITtransceive = true;
                    //             // if(mtxBuilder.Eval is XSpatialOperatorMk2.XEvaluatorNonlin evn) { // length-scale hack
                    //             //     foreach(var kv in AgglomeratedCellLengthScales) {
                    //             //         evn.CellLengthScales[kv.Key] = kv.Value;
                    //             //     }
                    //             // }
                    //             // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                    //             // return;
                    //         }

                    //         case LinearizationHint.GetJacobiOperator: using(new BlockTrace("XDG-LinearizationHint.GetJacobiOperator", ft, false)){
                    //             var op = GetJacobiXdgOperator();

                    //             if(JacobiParameterVars == null)
                    //                 JacobiParameterVars = op.InvokeParameterFactory(this.CurrentState);

                    //             op.InvokeParameterUpdate(time, __CurrentState, JacobiParameterVars);

                    //             var mtxBuilder = op.GetMatrixBuilder(LsTrk, Mapping, this.JacobiParameterVars, Mapping, LsTrkHistoryIndex);
                    //             mtxBuilder.time = time;
                    //             mtxBuilder.MPITtransceive = true;
                    //             foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
                    //                 mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                    //             }
                    //             mtxBuilder.ActivateSubgridBoundary(this.m_SubGrid.VolumeMask, m_SubGridBnd);
                    //             mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                    //             return;
                    //         }
                        // }
                    // } else {
                        // ++++++++++++++++++++++++
                        // only operator evaluation
                        // ++++++++++++++++++++++++


                        // throw new NotImplementedException();

                        // using(new BlockTrace("XDG-Evaluate", ft, false)) {
                        //     this.XdgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

                        //     var eval = XdgOperator.GetEvaluatorEx(this.LsTrk, __CurrentState, this.Parameters, Mapping, LsTrkHistoryIndex);
                        //     eval.time = time;

                        //     eval.MPITtransceive = true;
                        //     foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
                        //         eval.CellLengthScales[kv.Key] = kv.Value;
                        //     }
                        //     eval.Evaluate(1.0, 0.0, OpAffine);
                } else if(DgOperator != null) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++
                    // DG Branch
                    // +++++++++++++++++++++++++++++++++++++++++++++++

                    // throw new NotImplementedException();
                    if(OpMtx != null) {
                        // +++++++++++++++++++++++++++++
                        // Solver requires linearization
                        // +++++++++++++++++++++++++++++

                        Debug.Assert(OpMtx.InfNorm() == 0.0);
                        switch(DgOperator.LinearizationHint) {

                            case LinearizationHint.AdHoc: using(new BlockTrace("DG-LinearizationHint.AdHoc", ft)) {
                                throw new NotImplementedException();
                                // this.DgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());
                                // var mtxBuilder = DgOperator.GetMatrixBuilder(Mapping, this.Parameters, Mapping);
                                // mtxBuilder.time = time;
                                // mtxBuilder.MPITtransceive = true;
                                // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                                // return;
                            }

                            case LinearizationHint.FDJacobi: using(new BlockTrace("DG-LinearizationHint.FDJacobi", ft)){
                                throw new NotImplementedException();
                                // var mtxBuilder = DgOperator.GetFDJacobianBuilder(__CurrentState, this.Parameters, Mapping);
                                // mtxBuilder.time = time;
                                // mtxBuilder.MPITtransceive = true;
                                // mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                                // return;
                            }

                            case LinearizationHint.GetJacobiOperator: using(new BlockTrace("DG-LinearizationHint.GetJacobiOperator", ft)){
                                var op = GetJacobiDgOperator();

                                if(JacobiParameterVars == null)
                                    JacobiParameterVars = op.InvokeParameterFactory(__CurrentState);
                                op.InvokeParameterUpdate(time, __CurrentState, JacobiParameterVars);

                                var mtxBuilder = op.GetMatrixBuilder(Mapping, this.JacobiParameterVars, Mapping);
                                // mtxBuilder.ActivateSubgridBoundary(this.m_SubGrid.VolumeMask, m_SubGridBnd);
                                mtxBuilder.time = time;
                                mtxBuilder.MPITtransceive = true;
                                mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

                                // for (int i = 0; i < OpMtx.NoOfRows; i++){
                                //     for (int j = 0; j < OpMtx.NoOfRows; j++){
                                //         Console.Write(OpMtx[i,j] == 0.0 ? "o" : "x");
                                //     }
                                //     Console.WriteLine();
                                // }

                                PlotFormat format = new PlotFormat(
                                        Style: Styles.Points,
                                        pointType: PointTypes.Asterisk,
                                        pointSize: 0.5);

                                var gp = new Gnuplot.Gnuplot(baseLineFormat: format);
                                gp.PlotMatrixStructure(OpMtx.ToFullMatrixOnProc0());
                                gp.PlotDataFile("/home/klingenberg/Downloads/matrix.png", deferred: false);
                                // gp.PlotDataFile("matrix.png", deferred: false);
                                gp.Execute();
                                Console.ReadLine();
                                for (int i = 0; i < OpMtx.NoOfRows; i++)
                                {
                                    if (OpMtx.GetNoOfNonZerosPerRow(i) > 0)
                                    {
                                        OpMtx.SetDiagonalElement(i, 1.0);
                                    }
                                }
                                return;
                            }
                        }
                    } else {
                        // ++++++++++++++++++++++++
                        // only operator evaluation
                        // ++++++++++++++++++++++++

                        // throw new NotImplementedException();
                        using(new BlockTrace("DG-Evaluate", ft)) {
                            this.DgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

                            var eval = DgOperator.GetEvaluatorEx(__CurrentState, this.Parameters, Mapping);
                            eval.time = time;
                            eval.MPITtransceive = true;
                            eval.Evaluate(1.0, 0.0, OpAffine);
                        }
                    }
                } else {
                    throw new NotImplementedException();
                }
            }
        }
    }
}
