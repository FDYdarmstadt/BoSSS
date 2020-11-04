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

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Evaluation or linearization/matrix assembly of the operator;
    /// this delegate is, so-to-say, the connection between the used <see cref="ISpatialOperator"/> and its evaluation
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
    public delegate void OperatorEvalOrLin(out BlockMsrMatrix Matrix, out double[] Affine, out BlockMsrMatrix MassMatrix, DGField[] CurrentState, bool Linearization, out ISpatialOperator OberFrickelHack);
              

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


        ISpatialOperator m_AbstractOperator;

        /// <summary>
        /// spatial operator provided by <see cref="m_AssembleMatrix"/>;
        /// </summary>
        protected ISpatialOperator AbstractOperator {
            get {
                return m_AbstractOperator;
            }
            set {
                if(m_AbstractOperator == null)
                    m_AbstractOperator = value;
                if(!object.ReferenceEquals(m_AbstractOperator, value)) {
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
            if(RHS != null) {
                if(RHS.Count != Lraw)
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
        /// Preconditioner/solver for the linearized problem
        /// </summary>
        public ISolverSmootherTemplate Precond;
        


        /// <summary>
        /// Current linearization of the nonlinear operator: the linearized
        /// system is given as 
        /// <see cref="CurrentLin"/>*X = <see cref="LinearizationRHS"/>.
        /// </summary>
        protected MultigridOperator CurrentLin;

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
        /// <see cref="ISpatialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void EvaluateOperator(double alpha, IEnumerable<DGField> CurrentState, double[] Output, double HomotopyValue) {
            if(alpha != 1.0)
                throw new NotSupportedException("some moron has removed this");

            SetHomotopyValue(HomotopyValue);

            // the real call:
            this.m_AssembleMatrix(out BlockMsrMatrix OpMtxRaw, out double[] OpAffineRaw, out BlockMsrMatrix MassMtxRaw, CurrentState.ToArray(), false, out var Dummy);
            if(OpMtxRaw != null)
                // only evaluation ==> OpMatrix must be null
                throw new ApplicationException($"The provided {typeof(OperatorEvalOrLin).Name} is not correctly implemented.");
            this.AbstractOperator = Dummy;

            CurrentLin.TransformRhsInto(OpAffineRaw, Output, false);
        }

        protected void SetHomotopyValue(double HomotopyValue) {
            if(HomotopyValue < 0)
                throw new ArgumentOutOfRangeException();
            if(HomotopyValue > 1)
                throw new ArgumentOutOfRangeException();


            if(AbstractOperator == null && HomotopyValue != 1) {
                // do one call just to attain the spatial operator;
                // not very efficient, but ok for the moment

                var dummyState = this.CurrentLin.BaseGridProblemMapping.BasisS.Select(delegate (Basis b) {
                    DGField r;
                    if(b is XDGBasis xb)
                        r = new XDGField(xb);
                    else
                        r = new SinglePhaseField(b);
                    return r;
                }).ToArray();
          
                this.m_AssembleMatrix(out var OpMtxRaw, out _, out _, dummyState, false, out var _Dummy);
                Debug.Assert(OpMtxRaw == null); // only evaluation ==> OpMatrix must be null
                this.AbstractOperator = _Dummy;
            }

            if(AbstractOperator == null && HomotopyValue != 1.0) {
                throw new NotSupportedException("unable to attain operator for homotopy update");
            }

            if(AbstractOperator != null && HomotopyValue != AbstractOperator.CurrentHomotopyValue)
                // set homotopy value
                AbstractOperator.CurrentHomotopyValue = HomotopyValue;
        }



        /// <summary>
        /// Updating the <see cref="CurrentLin"/> -- operator;
        /// </summary>
        /// <param name="CurrentState">input, linearization point</param>
        /// <param name="U0">output, linearization point, after external update, transformed back</param>
        /// <param name="HomotopyValue">
        /// <see cref="ISpatialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void Update(IEnumerable<DGField> CurrentState, double[] U0, double HomotopyValue) {
            /*
            DGField[] U0fields = this.ProblemMapping.BasisS.Select(
                delegate(Basis b) {
                    DGField ret;
                    if (b is XDGBasis) {
                        XDGField xf = new XDGField(b as XDGBasis);
                        xf.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
                        ret = xf;
                    } else {
                        ret = new SinglePhaseField(b);
                    }
                    return ret;
                }).ToArray();

            CoordinateVector u0Raw = new CoordinateVector(U0fields);

            CurrentLin.TransformSolFrom(u0Raw, U0);
            */
           
            this.UpdateLinearization(CurrentState, HomotopyValue);
            CurrentLin.TransformSolInto(new CoordinateVector(CurrentState), U0);
        }


        /// <summary>
        /// Updating the <see cref="CurrentLin"/> -- operator;
        /// </summary>
        /// <param name="CurrentState">linearization point</param>
        /// <param name="HomotopyValue">
        /// <see cref="ISpatialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void UpdateLinearization(IEnumerable<DGField> CurrentState, double HomotopyValue) {
            if(!(this.ProblemMapping.BasisS.Count == CurrentState.Count()))
                throw new ArgumentException("mismatch in number of fields.");

            SetHomotopyValue(HomotopyValue);

            // the real call:
            this.m_AssembleMatrix(out BlockMsrMatrix OpMtxRaw, out double[] OpAffineRaw, out BlockMsrMatrix MassMtxRaw, CurrentState.ToArray(), true, out ISpatialOperator abstractOperator);
            AbstractOperator = abstractOperator;

            // blabla:
            CurrentLin = new MultigridOperator(this.m_AggBasisSeq, this.ProblemMapping,
                OpMtxRaw.CloneAs(), MassMtxRaw,
                this.m_MultigridOperatorConfig,
                AbstractOperator.DomainVar.Select(varName => AbstractOperator.FreeMeanValue[varName]).ToArray()); 

            OpAffineRaw = OpAffineRaw.CloneAs();
            if (this.RHSRaw != null)
                OpAffineRaw.AccV(-1.0, this.RHSRaw);
            if(LinearizationRHS == null || LinearizationRHS.Length != this.CurrentLin.Mapping.LocalLength)
                LinearizationRHS = new double[this.CurrentLin.Mapping.LocalLength];
            else
                LinearizationRHS.ClearEntries();
            CurrentLin.TransformRhsInto(OpAffineRaw, this.LinearizationRHS, true);
            this.LinearizationRHS.ScaleV(-1.0);
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



    }
}
