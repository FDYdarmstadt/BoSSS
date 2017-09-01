using System;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using BoSSS.Solution.Timestepping;
using BoSSS.Platform;
using BoSSS.Foundation.Comm;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;



namespace BoSSS.Solution {

#pragma warning disable 1572
#pragma warning disable 1573
#pragma warning disable 1574
#pragma warning disable 1591

    public class ImplicitEulerSubgrid : ImplicitTimeStepperSubgrid {
       
        /// <summary>
        /// Constructor for an implicit Euler scheme that is operating on a subgrid.
        /// No ISparseExt solver is needed!
        /// </summary>
        /// <param name="ctx"></param>
        /// <param name="solver">The ISparse solver to be used</param>
        /// <param name="spatialOpMtx"> Matrix of the operator on the full domain</param>
        /// <param name="spatialOpAffine">Affine part of the operator defined on the full domain</param>
        /// <param name="subgrid">Subgrid where the computation is assumed to be performed</param>
        /// <param name="fields">Subgrid Mapping of all fields</param>
        public ImplicitEulerSubgrid(Context ctx, ISparseSolver solver, MsrMatrix spatialOpMtx,
            IList<double> spatialOpAffine, SubGrid subgrid, SubgridCoordinateMapping fields, double Initialdt)
            : this(ctx, solver, AllTrue(fields.Fields.Count), spatialOpMtx, spatialOpAffine, subgrid, fields, Initialdt) { }

        /// <summary>
        /// Constructor for an implicit Euler scheme that is operating on a subgrid
        /// </summary>
        /// <param name="ctx"></param>
        /// <param name="solver">The ISparse solver to be used</param>
        /// <param name="temporalOp">
        /// Indicates, for each each equation whether it is
        /// <list type="bullet">
        ///   <item>(false) a side condition  i.e. a variable where no time derivative occurs in the equation, or</item>
        ///   <item>(true) a differential equation</item>
        /// </list>
        /// At least one equation must be time-dependent;
        /// Otherwise, the <see cref="BoSSS.Solution.Solvers.LinearSolver"/> should be used;
        /// </param>
        /// <param name="spatialOpMtx"> Matrix of the operator on the full domain</param>
        /// <param name="spatialOpAffine">Affine part of the operator defined on the full domain</param>
        /// <param name="subgrid">Subgrid where the computation is assumed to be performed</param>
        /// <param name="fields">Subgrid mapping of all fields</param>
        public ImplicitEulerSubgrid(Context ctx, ISparseSolver solver, bool[] temporalOp, MsrMatrix spatialOpMtx,
            IList<double> spatialOpAffine, SubGrid subgrid, SubgridCoordinateMapping fields, double Initialdt)
            : base(ctx, solver, temporalOp, spatialOpMtx, spatialOpAffine, subgrid, fields, Initialdt) { }
      
       
        /// <summary>
        /// Create an empty scheme
        /// </summary>
        public ImplicitEulerSubgrid() { }

        /// <summary>
        /// contains the sparse solver statistics acquired in the last call to <see cref="ImplicitTimeStepper.Perform"/>
        /// </summary>
        public SolverResult LastSolverResult;
        /*
        /// <summary>
        /// Solves the linear system (diag(1/<paramref name="dt"/>) +
        /// <em>M</em>) * x =
        /// <see cref=" SubgridMapping.subgridCoordinates"/> /
        /// <paramref name="dt"/> -
        /// <see cref="ImplicitTimeStepper.m_AffineOffset1"/> and writes the
        /// result to <see cref=" SubgridMapping.subgridCoordinates"/>.
        /// </summary>
        /// <param name="dt">The lenght of the timestep</param>
        protected override void PerformTimeStep(double dt) {
            using (var tr = new ilPSP.Tracing.FuncTrace()) {

                int n = SubgridMapping.SubgridNUpdate * SubgridMapping.MaxTotalNoOfCoordinatesPerCell;

                double[] rhs = (double[])m_CompressedAffine.Clone();
                BLAS.dscal(rhs.Length, -1.0, rhs, 1);

                for (int i = 0; i < n; i++) {
                    rhs[i] += 1.0 / dt * SubgridDGCoordinates[i];
                }

                tr.Info("Calling solver");

                LastSolverResult = m_Solver.Solve<double[], double[]>(SubgridDGCoordinates, rhs);
                //Console.WriteLine(LastSolverResult.NoOfIterations);
                SubgridMapping.subgridCoordinates = SubgridDGCoordinates;

            }
        }*/
        protected override void PerformTimeStep(double dt) {
            using (var tr = new ilPSP.Tracing.FuncTrace()) {
                int[] cells = m_Subgrid.SubgridIndex2LocalCellIndex;
                // int n = SubgridMapping.SubgridNUpdate * SubgridMapping.MaxTotalNoOfCoordinatesPerCell;
                int CoordinateOffset = 0;
                double[] rhs = (double[])m_CompressedAffine.Clone();
                //     BLAS.dscal(rhs.Length, -1.0, rhs, 1);
                /*  for (int j = 0; j < cells.Length; j++) {
                      CoordinateOffset = 0;
                      for (int g = 0; g < temporalOp.Length; g++) {

                          if (temporalOp[g]) {*/


                for (int g = 0; g < temporalOp.Length; g++) {

                    if (temporalOp[g]) {

                        for (int j = 0; j < cells.Length; j++) {
                            int loc = j * m_SubgridMapping.MaxTotalNoOfCoordinatesPerCell + CoordinateOffset;
                            for (int i = 0; i < m_SubgridMapping.Fields[g].Basis.MaximalLength; i++) {

                                rhs[loc] += 1.0 / dt * SubgridDGCoordinates[loc];
                                loc++;
                            }
                        }

                    }
                    CoordinateOffset += m_SubgridMapping.Fields[g].Basis.MaximalLength;
                }

                tr.Info("Calling solver");

                LastSolverResult = m_Solver.Solve<double[], double[]>(SubgridDGCoordinates, rhs);
                SubgridMapping.subgridCoordinates = SubgridDGCoordinates;
                Console.WriteLine("Solver" + LastSolverResult.NoOfIterations);
                if(LastSolverResult.Converged)
                Console.WriteLine("SolverConverged?:Yes");
                else
                    Console.WriteLine("SolverConverged?:No");
              
            }
        }
        /// <summary>
        /// Method that is called by the base class constructor. In this implicit euler scheme, it simply
        /// adds the diagonal entries to the operator matrix. 
        /// Currently only working for fields that should be included in the time stepping. 
        /// </summary>
        /// <param name="OperatorMatrix">The operator matrix that should be adjusted and passed onto the solver</param>
        /// <param name="dt">time step size</param>
        protected override void DefineMatrix(MsrMatrix OperatorMatrix, double dt) {

            int[] cells = m_Subgrid.SubgridIndex2LocalCellIndex;
            IList<Field> fields = m_SubgridMapping.Fields;
            int CoordinateOffset = 0;
           
            for (int g = 0; g < temporalOp.Length; g++) {
              
                if (temporalOp[g]) {
                    Field gamma = fields[g];
                 
                    for (int j = 0; j < cells.Length; j++) {
                        int iglob = j * m_SubgridMapping.MaxTotalNoOfCoordinatesPerCell + CoordinateOffset+(int)OperatorMatrix.RowPartition.i0;
                        for (int n = 0; n < gamma.Basis.MaximalLength; n++) {

                            double vadd = 1.0 / dt;
                            double v = OperatorMatrix.GetDiagElement(iglob);
                            v += vadd;
                            OperatorMatrix.SetDiagElement(iglob, v);
                            iglob++;
                        }
                    }
                   
                }
                CoordinateOffset += fields[g].Basis.MaximalLength;
            }
            m_Solver.DefineMatrix(OperatorMatrix);
            OperatorMatrix.SaveToTextFileSparse("Matrix.txt");
        }

    }

#pragma warning restore 1572
#pragma warning restore 1573
#pragma warning restore 1574
#pragma warning restore 1591

}
