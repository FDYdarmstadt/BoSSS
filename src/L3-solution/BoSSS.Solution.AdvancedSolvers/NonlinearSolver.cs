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

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using ilPSP.LinSolvers;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;
using System.CodeDom;
using BoSSS.Solution.Tecplot;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Evaluation or linearization/matrix assembly of the operator;
    /// this delegate is, so-to-say, the connection between the used <see cref="IDifferentialOperator"/> and its evaluation/linearization,
    /// which can be used to build a <see cref="MultigridOperator"/>
    /// </summary>
    /// <param name="Matrix"></param>
    /// <param name="Affine"></param>
    /// <param name="Linearization">
    /// - false: operator evaluation
    /// - true: linearization
    /// </param>
    /// <param name="MassMatrix"></param>
    /// <param name="CurrentState">
    /// Current state of the solution
    /// </param>
    /// <param name="OberFrickelHack">
    /// the original operator that somehow produced the matrix; yes, this API is convoluted piece-of-shit
    /// </param>
    /// <remarks>
    /// As a recipe, this function does:
    /// 1. compute linearization/evaluation at the current state , i.e. 
    /// 1.1 in the case of `Linearization == true`,
    ///     compute the operator matrix (<see cref="IDifferentialOperator.LinearizationHint"/>, <see cref="IDifferentialOperator.GetMatrixBuilder"/>, <see cref="IDifferentialOperator.GetJacobiOperator"/>, <see cref="IDifferentialOperator.GetFDJacobianBuilder"/>)
    /// 1.2 in the case of `Linearization == false`,
    ///     evaluate the operator (<see cref="IDifferentialOperator.GetEvaluatorEx"/>)
    /// 2. add timestepping terms (<see cref="IDifferentialOperator.TemporalOperator"/>), also depending on the actual timestepping scheme (e.g. BDF or Runge-Kutta)
    /// 3. Compute the mass matrix
    /// 4. perform the agglomeration (<see cref="LevelSetTracker.GetAgglomerator"/>
    /// </remarks>
    public delegate void OperatorEvalOrLin(out BlockMsrMatrix Matrix, out double[] Affine, out BlockMsrMatrix MassMatrix, DGField[] CurrentState, bool Linearization, out IDifferentialOperator OberFrickelHack);


    /// <summary>
    /// base-class for non-linear solvers
    /// </summary>
    public abstract class NonlinearSolver {

        /// <summary>
        /// ctor
        /// </summary>
        public NonlinearSolver(OperatorEvalOrLin __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq, MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig) {
            m_AssembleMatrix = __AssembleMatrix;
            m_AggBasisSeq = __AggBasisSeq.ToArray();
            m_MultigridOperatorConfig = __MultigridOperatorConfig;
        }

        /// <summary>
        /// Evaluation and linearization of PDE to solve
        /// </summary>
        protected OperatorEvalOrLin m_AssembleMatrix;


        IDifferentialOperator m_AbstractOperator;

        /// <summary>
        /// spatial operator provided by <see cref="m_AssembleMatrix"/>;
        /// </summary>
        protected IDifferentialOperator AbstractOperator {
            get {
                return m_AbstractOperator;
            }
            set {
                if (m_AbstractOperator == null)
                    m_AbstractOperator = value;
                if (!object.ReferenceEquals(m_AbstractOperator, value)) {
                    throw new ApplicationException("Hey buddy, you should not change the operator during the solver life-time");
                }
            }
        }


        /// <summary>
        /// Multigrid basis
        /// - 1st index: Multigrid level
        /// - 2nd index: variable index
        /// </summary>
        protected AggregationGridBasis[][] m_AggBasisSeq;

        /// <summary>
        /// required for construction of <see cref="CurrentLin"/>
        /// </summary>
        protected MultigridOperator.ChangeOfBasisConfig[][] m_MultigridOperatorConfig;

        /// <summary>
        /// Called at every iteration; the arguments are 
        ///  - iteration index 
        ///  - current solution 
        ///  - current residual 
        ///  - current multigrid operator
        /// </summary>
        public event Action<int, double[], double[], MultigridOperator> IterationCallback;

        /// <summary>
        /// Triggers <see cref="IterationCallback"/>.
        /// </summary>
        protected void OnIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            if (IterationCallback != null)
                IterationCallback(iterIndex, currentSol, currentRes, Mgop);
        }



        /// <summary>
        /// Helper routine for the initial phase of <see cref="SolverDriver{S}"/>
        /// </summary>
        /// <param name="X">initial guess</param>
        /// <param name="RHS"></param>
        /// <param name="Sol1">
        /// The initial solution, transformed to the aggregation multigrid basis, see <see cref="m_AggBasisSeq"/>.
        /// </param>
        /// <param name="Res1">
        /// The residual of the initial solution
        /// </param>
        protected void Init<S>(CoordinateVector X, S RHS, out double[] Sol1, out double[] Res1)
            where S : IList<double> //
        {
            this.ProblemMapping = X.Mapping;

            int Lraw = X.Mapping.LocalLength;  // length of Solution/RHS in original space
            if (RHS != null) {
                if (RHS.Count != Lraw)
                    throw new ArgumentException();
                this.RHSRaw = RHS.ToArray();
            } else {
                this.RHSRaw = null;
            }
            this.UpdateLinearization(X.Mapping.Fields, 1.0);

            int Ltrf = this.CurrentLin.Mapping.LocalLength;


            // set initial guess (input) as first approximation to the solution ...
            Sol1 = new double[Ltrf];
            this.CurrentLin.TransformSolInto(X, Sol1);

            // ... and evaluate its residual
            Res1 = new double[Ltrf];
            this.EvalLinearizedResidual(Sol1, ref Res1);
        }


        /// <summary>
        /// Template for implementation of the solver routine.
        /// </summary>
        /// <typeparam name="S"></typeparam>
        /// <param name="X">
        /// On entry, an initial guess to the linear system.
        /// On exit, hopefully the solution to the nonlinear system.
        /// </param>
        /// <param name="RHS">
        /// If not equal null, must be passed to <see cref="RHSRaw"/>.
        /// </param>
        /// <returns>
        /// - true: solver algorithm successfully converged
        /// - false: something went wrong
        /// </returns>
        abstract public bool SolverDriver<S>(CoordinateVector X, S RHS)
           where S : IList<double>;

        /// <summary>
        /// Number of nonlinear iterations in last call to <see cref="SolverDriver{S}(CoordinateVector, S)"/>
        /// </summary>
        public int NoOfNonlinearIter {
            get;
            protected set;
        }


        /// <summary>
        /// Preconditioner/solver configuration for the linearized problem
        /// </summary>
        public ISolverFactory PrecondConfig;





        /// <summary>
        /// Current linearization of the nonlinear operator: the linearized
        /// system is given as 
        /// <see cref="CurrentLin"/>*X = <see cref="LinearizationRHS"/>.
        /// </summary>
        protected MultigridOperator CurrentLin;

        /// <summary>
        /// number of DOFs in linearization of last iteration
        /// </summary>
        public long EssentialDOFs {
            get {
                return CurrentLin?.Mapping.TotalLength ?? 0;
            }
        }


        /// <summary>
        /// Optional RHS to the nonlinear system, 
        /// </summary>
        protected double[] RHSRaw;

        /// <summary>
        /// Current affine part of the nonlinear operator.
        /// </summary>
        protected double[] LinearizationRHS;



        /// <summary>
        /// Evaluation of the nonlinear operator.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="CurrentState">
        /// Current state of DG fields
        /// </param>
        /// <param name="Output"></param>
        /// <param name="HomotopyValue">
        /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        /// <param name="ApplyRef">
        /// apply additional modification due to free-mean-value fixing (aka. pressure reference point), <see cref="MultigridOperator.FreeMeanValue"/>
        /// </param>        
        protected void EvaluateOperator(double alpha, IEnumerable<DGField> CurrentState, double[] Output, double HomotopyValue, bool ApplyRef = false) {
            if (alpha != 1.0)
                throw new NotSupportedException("some moron has removed this");

            SetHomotopyValue(HomotopyValue);

            // the real call:
            this.m_AssembleMatrix(out BlockMsrMatrix DummyMtx, out double[] OpEvalRaw, out BlockMsrMatrix MassMtxRaw, CurrentState.ToArray(), false, out var abstractOp);
            if (DummyMtx != null)
                // only evaluation ==> OpMatrix must be null
                throw new ApplicationException($"The provided {typeof(OperatorEvalOrLin).Name} is not correctly implemented.");
            this.AbstractOperator = abstractOp;

#if DEBUG
            const int TEST_INTERVALL = 10;
#else
            const int TEST_INTERVALL = 1000;
#endif
            if (EvaluationCounter % TEST_INTERVALL == 0) { // do the following, expensive check only for every TEST_INTERVALL-th evaluation.
                // Comparison of linearization and evaluation:
                // -------------------------------------------
                //
                // Note that, in BoSSS, currently the Linearization of f(u) around u0 is defines as
                //     f(u) ≈ M(u0)*u + b(u0),
                // instead of the typical Taylor series representation f(u) ≈ f(u0) + ∂f(u0)*(u-u0).
                // The relation between the BoSSS-representation and the Taylor series is
                //     M(u0) = ∂f(u0),
                //     b(u0) = f(u0) - ∂f(u0)*u0.
                // Therefore, we check that
                //     M(u0)*u0 + b(u0) = f(u0).
                //


                this.m_AssembleMatrix(out BlockMsrMatrix LinMtx, out double[] OpAffine, out _, CurrentState.ToArray(), true, out _);
                var Check = OpAffine.CloneAs();
                LinMtx.SpMV(1.0, new CoordinateVector(CurrentState), 1.0, Check);

                var err = Check.CloneAs();
                err.AccV(-1.0, OpEvalRaw);
                double l2_err = err.MPI_L2Norm();
                double comp = Math.Sqrt(Math.Max(OpEvalRaw.MPI_L2Norm(), Check.MPI_L2Norm()) * BLAS.MachineEps + BLAS.MachineEps);
                
                if (l2_err > comp) {
                    Console.Error.WriteLine($"Mismatch between operator linearization and evaluation: Operator matrix-Jacobian distance: {l2_err}, relative: {l2_err/comp} (comparison value {comp})");
                    //throw new ArithmeticException($"Mismatch between operator linearization and evaluation: Operator matrix-Jacobian distance: { l2_err }, relative: { l2_err/comp} (comparison value { comp})");
				}
            }
            EvaluationCounter++;

            CurrentLin.TransformRhsInto(OpEvalRaw, Output, ApplyRef);
        }

        int EvaluationCounter = 0;

        protected void SetHomotopyValue(double HomotopyValue) {
            if (HomotopyValue < 0)
                throw new ArgumentOutOfRangeException();
            if (HomotopyValue > 1)
                throw new ArgumentOutOfRangeException();


            if (AbstractOperator == null && HomotopyValue != 1) {
                // do one call just to attain the spatial operator;
                // not very efficient, but ok for the moment

                var dummyState = this.CurrentLin.BaseGridProblemMapping.BasisS.Select(delegate (Basis b) {
                    DGField r;
                    if (b is XDGBasis xb)
                        r = new XDGField(xb);
                    else
                        r = new SinglePhaseField(b);
                    return r;
                }).ToArray();

                this.m_AssembleMatrix(out var OpMtxRaw, out _, out _, dummyState, false, out var _Dummy);
                Debug.Assert(OpMtxRaw == null); // only evaluation ==> OpMatrix must be null
                this.AbstractOperator = _Dummy;
            }

            if (AbstractOperator == null && HomotopyValue != 1.0) {
                throw new NotSupportedException("unable to attain operator for homotopy update");
            }

            if (AbstractOperator != null && HomotopyValue != AbstractOperator.CurrentHomotopyValue)
                // set homotopy value
                AbstractOperator.CurrentHomotopyValue = HomotopyValue;
        }



        /// <summary>
        /// Updating the <see cref="CurrentLin"/> -- operator;
        /// </summary>
        /// <param name="CurrentState">input, linearization point</param>
        /// <param name="U0">output, linearization point, after external update, transformed back</param>
        /// <param name="HomotopyValue">
        /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void Update(IEnumerable<DGField> CurrentState, double[] U0, double HomotopyValue) {

            this.UpdateLinearization(CurrentState, HomotopyValue);
            CurrentLin.TransformSolInto(new CoordinateVector(CurrentState), U0);
        }


        /// <summary>
        /// Updating the <see cref="CurrentLin"/> -- operator;
        /// </summary>
        /// <param name="CurrentState">linearization point</param>
        /// <param name="HomotopyValue">
        /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void UpdateLinearization(IEnumerable<DGField> CurrentState, double HomotopyValue) {
            if (!(this.ProblemMapping.BasisS.Count == CurrentState.Count()))
                throw new ArgumentException("mismatch in number of fields.");

            SetHomotopyValue(HomotopyValue);

            // the real call:
            this.m_AssembleMatrix(out BlockMsrMatrix OpMtxRaw, out double[] OpAffineRaw, out BlockMsrMatrix MassMtxRaw, CurrentState.ToArray(), true, out IDifferentialOperator abstractOperator);
            AbstractOperator = abstractOperator;

            // blabla:
            CurrentLin = new MultigridOperator(this.m_AggBasisSeq, this.ProblemMapping,
                OpMtxRaw.CloneAs(), MassMtxRaw,
                this.m_MultigridOperatorConfig,
                AbstractOperator);



            OpAffineRaw = OpAffineRaw.CloneAs();
            if (this.RHSRaw != null)
                OpAffineRaw.AccV(-1.0, this.RHSRaw);
            if (LinearizationRHS == null || LinearizationRHS.Length != this.CurrentLin.Mapping.LocalLength)
                LinearizationRHS = new double[this.CurrentLin.Mapping.LocalLength];
            else
                LinearizationRHS.ClearEntries();
            CurrentLin.TransformRhsInto(OpAffineRaw, this.LinearizationRHS, true);
            this.LinearizationRHS.ScaleV(-1.0);
        }


        /// <summary>
        /// This method tests that, 
        /// if any entry in <see cref="IDifferentialOperator.FreeMeanValue"/> is true, 
        /// a change in the mean/average value of the respective variable must **not** have any effect on the residual.
        /// </summary>
        public void TestFreeMeanValue(CoordinateVector SolutionVec, double HomotopyValue) {

            int L = this.CurrentLin.Mapping.LocalLength;

            if (CurrentLin.FreeMeanValue.Any()) {

                double[] ResidualBeforMeanCor = new double[L];
                double[] ResidualAfterMeanCor = new double[L];

                EvaluateOperator(1, SolutionVec.Mapping.Fields, ResidualBeforMeanCor, HomotopyValue);
                double SolNormA = SolutionVec.MPI_L2Norm();


                DGField[] flds = SolutionVec.Mapping.Fields.ToArray();
                bool[] FreeMeanValue = CurrentLin.FreeMeanValue;
                if (flds.Length != FreeMeanValue.Length)
                    throw new ApplicationException();

                const double arbitrary_distortion_value = 2000.1234;

                for (int iFld = 0; iFld < flds.Length; iFld++) {
                    if (FreeMeanValue[iFld]) {
                        flds[iFld].AccConstant(arbitrary_distortion_value);
                    }
                }

                EvaluateOperator(1, SolutionVec.Mapping.Fields, ResidualAfterMeanCor, HomotopyValue);
                double SolNormB = SolutionVec.L2Norm();

                for (int iFld = 0; iFld < flds.Length; iFld++) {
                    if (FreeMeanValue[iFld]) {
                        flds[iFld].AccConstant(-arbitrary_distortion_value);
                    }
                }
                

                //double[] ResidualDifference = ResidualAfterMeanCor.CloneAs();
                //ResidualDifference.AccV(-1.0, ResidualBeforMeanCor);
                //DGField[] ResidualDifferenceDg = this.CurrentLin.ProlongateRhsToDg(ResidualDifference, "residualDifference");
                //Tecplot.Tecplot.PlotFields(ResidualDifferenceDg, "ResidualDifference", 0, 2);

                double Dist = ResidualBeforMeanCor.MPI_L2Dist(ResidualAfterMeanCor);
                double ResNormA = ResidualBeforMeanCor.MPI_L2Norm();
                double ResNormB = ResidualAfterMeanCor.MPI_L2Norm();



                double RefVal = Math.Max(Math.Max(BLAS.MachineEps.Sqrt(), Math.Max(ResNormA, ResNormB) * 1e-7), Math.Abs(SolNormA - SolNormB) * 1e-7);
                if (Dist > RefVal) {
                    Console.Error.WriteLine($"Something seems wrong with `FreeMeanValue`: drastic change of operator residual; Original residual: {ResNormA}; residual after distortion {ResNormB}; distance is {Dist}, reference value {RefVal}");
                    //throw new ArithmeticException($"Something seems wrong with `FreeMeanValue`: drastic change of operator residual; Original residual: {ResNormA}; residual after distortion {ResNormB}; distance is {Dist}, reference value {RefVal}");
                }
            }

        }


        /// <summary>
        /// DG Coordinate mapping for the base DG/XDG space, without any multigrid fuzz
        /// </summary>
        protected UnsetteledCoordinateMapping ProblemMapping {
            get;
            set;
        }

        /// <summary>
        /// Residual of linearized system, i.e.
        /// <paramref name="out_Resi"/> := RHS - M*<paramref name="in_U"/>
        /// </summary>
        protected void EvalLinearizedResidual(double[] in_U, ref double[] out_Resi) {
            if (this.LinearizationRHS.Length != in_U.Length)
                throw new ApplicationException("internal error");
            if (out_Resi.Length != in_U.Length)
                out_Resi = new double[in_U.Length];
            out_Resi.SetV(this.LinearizationRHS, 1.0);
            CurrentLin.OperatorMatrix.SpMV(-1.0, in_U, 1.0, out_Resi);
        }

        /// <summary>
        /// Inner product with respect to the current mass matrix.
        /// 
        /// Note: this is the canonical inner product of the underlying DG-space, since 
        /// for a DG/XDG field represented in an arbitrary basis $` \phi_{j} $` one verifies that
        /// ```math
        ///         (u, v) = 
        ///     \int_\Omega 
        ///         \left( \sum_{j} \phi_{j} \tilde{u}_{j} \right) 
        ///         \left( \sum_{l} \phi_{l} \tilde{v}_{l} \right) 
        ///     dV = 
        ///     \sum_{j l} \tilde{u}_{j} \tilde{v}_{l} ( \phi_{j}, \phi_{l} )
        ///     =
        ///       \tilde{u}^T M \tilde{v},
        /// ```
        /// where $`M `$ denotes the mass matrix ($` M_{j l} = ( \phi_ { j}, \phi_ { l} )  `$).
        /// </summary>
        protected double InnerProduct<T1, T2>(T1 vecA, T1 vecB)
            where T1 : IList<double>
            where T2 : IList<double> //
        {
            int L = this.CurrentLin.Mapping.LocalLength;

            if (vecA.Count != L)
                throw new ArgumentException("wrong length.");
            if (vecB.Count != L)
                throw new ArgumentException("wrong length.");


            var MaMa = this.CurrentLin.MassMatrix;

            if (MaMa == null) {
                return vecA.MPI_InnerProd(vecB);
            } else {
                double[] tmp = new double[L];
                MaMa.SpMV(1.0, vecB, 0.0, tmp);
                return vecA.MPI_InnerProd(tmp);
            }
        }

        /// <summary>
        /// Norm induced by <see cref="InnerProduct{T1, T2}(T1, T1)"/>
        /// </summary>
        protected double Norm<T1>(T1 vec)
            where T1 : IList<double> //
        {
            return InnerProduct<T1,T1>(vec, vec).Sqrt();

            /*
            int L0 = CurrentLin.BaseGridProblemMapping.LocalLength;
            double[] tmp = new double[L0];
            //CurrentLin.TransformSolFrom(tmp, vec);
            CurrentLin.TransformRhsFrom(tmp, vec);

            var XFields = CurrentLin.BaseGridProblemMapping.BasisS.Select(bs => new XDGField(bs as XDGBasis)).ToArray();
            var Coordv = new CoordinateVector(XFields);
            Coordv.SetV(tmp);

            double norm_acc = 0;
            for (int i = 0; i < XFields.Length; i++) {
                norm_acc += XFields[i].L2NormSpecies("A").Pow2();
                norm_acc += XFields[i].L2NormSpecies("B").Pow2();
            }

            //Console.WriteLine("New art of norm: " + norm_acc.Sqrt());
            return norm_acc.Sqrt();
            */
        }
    }
}
