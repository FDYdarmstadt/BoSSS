/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation.Grid;
using ilPSP.LinSolvers;
using System.Collections.Generic;
using static BoSSS.Foundation.DifferentialOperator;

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
        /// List of parameter fields, correlates with <see cref="DifferentialOperator.ParameterVar"/>.
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
        /// Turn MPI sending/receiving of parameters and domain fields on/off.
        /// </summary>
        bool MPITtransceive { get; set; }


        /// <summary>
        /// The pointer to a owner object, which totally contradicts the original idea of object-orientation. Hehe.
        /// </summary>
        IDifferentialOperator Owner { get; }
    }

    /// <summary>
    /// Explicit operator evaluation
    /// </summary>
    public interface IEvaluatorNonLin : IEvaluator {

        /// <summary>
        /// Values for trial variables
        /// </summary>
        CoordinateMapping DomainFields { get; }
        
        /// <summary>
        /// Evaluation call
        /// </summary>
        /// <typeparam name="Tout"></typeparam>
        /// <param name="alpha">result scaling</param>
        /// <param name="beta">scaling of accumulator</param>
        /// <param name="output">accumulation buffer</param>
        /// <param name="outputBndEdge">temporary hack, introduced for Compressible Local-Timestepping solver </param>
        void Evaluate<Tout>(double alpha, double beta, Tout output, double[] outputBndEdge = null) where Tout : IList<double>;
    }

    /// <summary>
    /// Evaluation of linear/linearized spatial operators, i.e. matrix assembly.
    /// </summary>
    public interface IEvaluatorLinear : IEvaluator {

        /// <summary>
        /// Computation of matrix \f$ \mathcal{M} \f$ *and* affine offset \f$ \tilde{b} \f$; the
        /// DG operator <see cref="IEvaluator.Owner"/> is discretized as
        /// \f[
        ///    \mathcal{M} \cdot \tilde{U} + \tilde{b}
        /// \f]
        /// where \f$ \tilde{U} \f$ are the DG coordinates of the trial function.
        /// </summary>
        /// <param name="Matrix">
        /// On entry, some pre-allocated matrix; on exit, the operator matrix, scaled by <paramref name="alpha"/>, will be accumulated.
        /// - <see cref="ISparseMatrix.RowPartitioning"/> correlates with <see cref="IEvaluator.CodomainMapping"/>.
        /// - <see cref="ISparseMatrix.ColPartition"/> correlates with <see cref="IEvaluator.DomainMapping"/>
        /// </param>
        /// <param name="AffineOffset">
        /// On entry, some pre-allocated vector; on exit, the affine vector, scaled by <paramref name="alpha"/>, will be accumulated.
        /// Length correlates with <see cref="IEvaluator.CodomainMapping"/>.
        /// </param>
        /// <param name="alpha">
        /// scaling factor 
        /// </param>
        void ComputeMatrix<M, V>(M Matrix, V AffineOffset, double alpha = 1.0) where M : IMutableMatrixEx where V : IList<double>;

        /// <summary>
        /// only the affine part of <see cref="ComputeMatrix{M, V}(M, V, double)"/>
        /// </summary>
        void ComputeAffine<V>(V AffineOffset) where V : IList<double>;
    }

}