using BoSSS.Foundation.Grid;
using ilPSP.LinSolvers;
using System.Collections.Generic;
using static BoSSS.Foundation.SpatialOperator;

namespace BoSSS.Foundation {

    /// <summary>
    /// Evaluation and matrix assembly of Spatial operators
    /// </summary>
    public interface IEvaluator {

        /// <summary>
        /// Grid, on which this evaluator operates on.
        /// </summary>
        IGridData GridData { get; }

        /// <summary>
        /// coordinate mapping for the co-domain variables (test functions, row mapping for matrix assembly, output of operator).
        /// </summary>
        UnsetteledCoordinateMapping CodomainMapping { get; }

        /// <summary>
        /// coordinate mapping for the domain variables (trial functions, column mapping for matrix assembly, input of operator).
        /// </summary>
        UnsetteledCoordinateMapping DomainMapping { get; }

        /// <summary>
        /// List of parameter fields, correlates with <see cref="SpatialOperator.ParameterVar"/>.
        /// </summary>
        IList<DGField> Parameters { get; }
        
        /// <summary>
        /// Turns on  the sub-grid thing uses for local-time-stepping.
        /// </summary>
        void ActivateSubgridBoundary(BoSSS.Foundation.Grid.CellMask sgrd, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge);

        /// <summary>
        /// Actual value set by <see cref="ActivateSubgridBoundary(CellMask, SubGridBoundaryModes)"/>
        /// </summary>
        SubGridBoundaryModes SubGridBoundaryTreatment { get; }
        
        /// <summary>
        /// Time passed e.g. to <see cref="CommonParams.time"/>, <see cref="CommonParamsBnd.time"/> and <see cref="CommonParamsVol.time"/>.
        /// </summary>
        double time { get; set; }

        /// <summary>
        /// The pointer to a owner object, which totally contradicts the original idea of object-orientation. Hehe.
        /// </summary>
        SpatialOperator Owner { get; }

        /// <summary>
        /// Turn MPI sending/receiving of parameters and domain fields on/off.
        /// </summary>
        bool MPITtransceive { get; set; }

        /// <summary>
        /// Stuff passed to equation components which implement <see cref="IEquationComponentCoefficient"/>.
        /// </summary>
        CoefficientSet OperatorCoefficients { get; set; }
    }

    public interface IEvaluatorNonLin : IEvaluator {

        CoordinateMapping DomainFields { get; }


        void Evaluate<Tout>(double alpha, double beta, Tout output, double[] outputBndEdge = null) where Tout : IList<double>;
        
    }

    /// <summary>
    /// Evaluation of linear/linearized spatial operators, i.e. matrix assembly.
    /// </summary>
    public interface IEvaluatorLinear : IEvaluator {

        void ComputeMatrix<M, V>(M Matrix, V AffineOffset) where M : IMutableMatrixEx where V : IList<double>;

        void ComputeAffine<V>(V AffineOffset) where V : IList<double>;
    }


}