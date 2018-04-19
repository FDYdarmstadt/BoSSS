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

        IList<DGField> Parameters { get; }



        void ActivateSubgridBoundary(BoSSS.Foundation.Grid.CellMask sgrd, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge);

        SubGridBoundaryModes SubGridBoundaryTreatment { get; }


        /// <summary>
        /// Time passed e.g. to <see cref="CommonParams.time"/>, <see cref="CommonParamsBnd.time"/> and <see cref="CommonParamsVol.time"/>.
        /// </summary>
        double time { get; set; }



        SpatialOperator Owner { get; }

        /// <summary>
        /// Turn MPI sending/receiving of parameters and domain fields on/off.
        /// </summary>
        bool MPITtransceive { get; set; }

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