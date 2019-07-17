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
using ilPSP.LinSolvers;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Memory-intensive ortho-normalization scheme; should converge with similar rate as GMRES (at twice the memory consumption), 
    /// but is able to use multiple pre-conditioners, see <see cref="PrecondS"/>.
    /// </summary>
    public class OrthonormalizationScheme : ISolverSmootherTemplate, ISolverWithCallback {
        
        MultigridOperator m_mgop;

        public void Init(MultigridOperator op) {
            var M = op.OperatorMatrix;
            var MgMap = op.Mapping;
            this.m_mgop = op;

            if(!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if(!M.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");

            foreach(var pc in PrecondS) {
                pc.Init(m_mgop);
            }
        }

        
        public ISolverSmootherTemplate[] PrecondS;



        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (new FuncTrace()) {
                // init 
                // ====

                int L = this.m_mgop.Mapping.LocalLength;
                if (X.Count != L)
                    throw new ArgumentException();
                if (B.Count != L)
                    throw new ArgumentException();

                var Mtx = this.m_mgop.OperatorMatrix;


                // residual of initial guess
                // =========================

                // history of solutions and residuals (max vector length 'MaxKrylovDim')
                List<double[]> SolHistory = new List<double[]>();
                List<double[]> MxxHistory = new List<double[]>();

                double[] Correction = new double[L];
                double[] Mxx = new double[L];
                double[] CurrentSol = new double[L];
                double[] CurrentRes = new double[L];

                CurrentSol.SetV(X, 1.0);
                CurrentRes.SetV(B, 1.0);
                Mtx.SpMV(-1.0, CurrentSol, 1.0, CurrentRes);
                int KrylovDim = 0;

                double[] Residual0 = CurrentRes.CloneAs();
                double[] Solution0 = CurrentSol.CloneAs();

                List<double> _R = new List<double>();

                // diagnostic output
                if (this.IterationCallback != null)
                    this.IterationCallback(0, CurrentSol.CloneAs(), CurrentRes.CloneAs(), this.m_mgop);

                // iterations...
                // =============
                double[] PreviousRes = new double[L];

                //MultidimensionalArray raw_Mxx = MultidimensionalArray.Create(L, MaxIter + 1);
                //MultidimensionalArray ortho_Mxx = MultidimensionalArray.Create(L, MaxIter + 1);

                MultidimensionalArray MassMatrix = MultidimensionalArray.Create(MaxKrylovDim, MaxKrylovDim);

                int PCcounter = 0;
                double[] prevAlpha = null;
                for (int iIter = 0; iIter < MaxIter; iIter++) {
                    Debug.Assert(SolHistory.Count == MxxHistory.Count);
                    Debug.Assert(SolHistory.Count == KrylovDim);

                    // select preconditioner
                    var Precond = PrecondS[PCcounter];
                    PCcounter++;
                    if (PCcounter >= PrecondS.Length) {
                        PCcounter = 0;
                        m_ThisLevelIterations++; // because we abuse the Orthonormalization to do some multi-grid stuff, 
                        //                          we only count every full cycle of preconditiones.
                    }

                    // solve the residual equation: M*Correction = prev. Residual
                    PreviousRes.SetV(CurrentRes);
                    Correction.ClearEntries();
                    Precond.Solve(Correction, PreviousRes);

                    // compute M*Correction
                    Mtx.SpMV(1.0, Correction, 0.0, Mxx);

                    // orthonormalize the Mxx -- vector with respect to the previous ones.
                    Debug.Assert(KrylovDim == MxxHistory.Count);
                    Debug.Assert(KrylovDim == SolHistory.Count);

                    //raw_Mxx.SetColumn(KrylovDim, Mxx);

                    for (int i = 0; i < KrylovDim; i++) {
                        Debug.Assert(!object.ReferenceEquals(Mxx, MxxHistory[i]));
                        double beta = GenericBlas.InnerProd(Mxx, MxxHistory[i]).MPISum();
                        Mxx.AccV(-beta, MxxHistory[i]);
                        Correction.AccV(-beta, SolHistory[i]);
                    }
                    {
                        double gamma = 1.0 / GenericBlas.L2NormPow2(Mxx).MPISum().Sqrt();
                        Mxx.ScaleV(gamma);
                        Correction.ScaleV(gamma);
                    }

                    // the following lines should produce the identity matrix
                    for (int i = 0; i < KrylovDim; i++) {
                        MassMatrix[i, KrylovDim] = GenericBlas.InnerProd(Mxx, MxxHistory[i]).MPISum();
                    }
                    MassMatrix[KrylovDim, KrylovDim] = GenericBlas.L2NormPow2(Mxx).MPISum();
                    //


                    //ortho_Mxx.SetColumn(KrylovDim, Mxx);

                    MxxHistory.Add(Mxx.CloneAs());
                    SolHistory.Add(Correction.CloneAs());
                    KrylovDim++;

                    

                    bool updateEveryIteration = false;


                    // RHS of the minimization problem (LHS is identity matrix)
                    if (!updateEveryIteration) {
                        _R.Add(GenericBlas.InnerProd(MxxHistory.Last(), Residual0).MPISum());
                    } else {
                        _R.Clear();
                        for (int i = 0; i < KrylovDim; i++) {
                            _R.Add(GenericBlas.InnerProd(MxxHistory[i], Residual0).MPISum());
                        }
                    }

                    // compute accelerated solution
                    //double[] alpha = _R.ToArray(); // factors for re-combining solutions
                    double[] alpha;
                    {
                        double[] minimi_rhs = _R.ToArray();
                        var minimi_lhs = MassMatrix.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { KrylovDim - 1, KrylovDim - 1 }).CloneAs();
                        alpha = new double[KrylovDim];
                        minimi_lhs.Solve(alpha, minimi_rhs);
                    }
                    if(prevAlpha != null) {
                        var del = alpha.GetSubVector(0, prevAlpha.Length);
                        del.AccV(-1.0, prevAlpha);
                    }
                    prevAlpha = alpha;

                    Console.WriteLine("Correction factor: " + alpha.Last() + ", solution " + SolHistory.Last().L2Norm() + " Resi " + Mxx.L2Norm());


                    Debug.Assert(alpha.Length == SolHistory.Count);
                    Debug.Assert(alpha.Length == MxxHistory.Count);
                    Debug.Assert(alpha.Length == KrylovDim);
                    CurrentSol.SetV(Solution0, 1.0);
                    for (int i = 0; i < KrylovDim; i++)
                        CurrentSol.AccV(alpha[i], SolHistory[i]);

                    // compute new Residual
                    CurrentRes.SetV(B);
                    Mtx.SpMV(-1.0, CurrentSol, 1.0, CurrentRes);
                    double crL2 = CurrentRes.L2Norm();

                    // diagnostic output
                    if (this.IterationCallback != null)
                        this.IterationCallback(iIter + 1, CurrentSol.CloneAs(), CurrentRes.CloneAs(), this.m_mgop);

                    //{
                    //    var gdat = m_mgop.BaseGridProblemMapping.GridDat;
                    //    var basis = m_mgop.BaseGridProblemMapping.BasisS[0];

                    //    var dgCurrentSol = new Foundation.SinglePhaseField(basis, "Solution");
                    //    var dgResidual = new Foundation.SinglePhaseField(basis, "Residual");
                    //    var dgCorrection = new Foundation.SinglePhaseField(basis, "Correction");

                    //    m_mgop.TransformRhsFrom(dgResidual.CoordinateVector, CurrentRes);
                    //    m_mgop.TransformSolFrom(dgCurrentSol.CoordinateVector, CurrentSol);
                    //    m_mgop.TransformSolFrom(dgCorrection.CoordinateVector, SolHistory.Last());
                    //    dgCorrection.Scale(alpha.Last());

                    //    Tecplot.Tecplot.PlotFields(new Foundation.DGField[] { dgCurrentSol, dgResidual, dgCorrection},  "OrthoScheme-" + iIter, iIter, 2);

                    //}


                    if (crL2 < Tolerance) {
                        //Console.WriteLine("    Kcy converged:");
                        //for (int iii = 0; iii < KrylovDim; iii++) {
                        //    Console.WriteLine("       fac #" + iii + "  :  " + alpha[iii]);
                        //}
                        if (PCcounter > 0)
                            m_ThisLevelIterations += 1;

                        m_Converged = true;
                        break;
                    }

                    if (updateEveryIteration) {
                        Solution0.SetV(CurrentSol);
                        Residual0.SetV(CurrentRes);
                    }

                    if (KrylovDim >= MaxKrylovDim) {
                        if (this.Restarted) {
                            // restarted version of the algorithm
                            // ++++++++++++++++++++++++++++++++++

                            MxxHistory.Clear();
                            SolHistory.Clear();
                            _R.Clear();
                            KrylovDim = 0;
                            Residual0.SetV(CurrentRes);
                            Solution0.SetV(CurrentSol);
                        } else {
                            // throw-away version of the algorithm
                            // +++++++++++++++++++++++++++++++++++

                            int i_leastSig = alpha.IndexOfMin(x => x.Abs());
                            MxxHistory.RemoveAt(i_leastSig);
                            SolHistory.RemoveAt(i_leastSig);
                            KrylovDim--;

                            for (int i = i_leastSig; i < KrylovDim; i++) {
                                for (int j = 0; j <= KrylovDim; j++) {
                                    MassMatrix[i, j] = MassMatrix[i + 1, j];
                                }
                            }
                            for (int i = i_leastSig; i < KrylovDim; i++) {
                                for (int j = 0; j <= KrylovDim; j++) {
                                    MassMatrix[j, i] = MassMatrix[j, i + 1];
                                }
                            }

                            Residual0.SetV(CurrentRes);
                            Solution0.SetV(CurrentSol);

                            _R.Clear();
                            foreach (double[] mxx in MxxHistory) {
                                _R.Add(GenericBlas.InnerProd(mxx, Residual0).MPISum());
                            }
                        }
                    }
                }


                X.SetV(CurrentSol, 1.0);
                //raw_Mxx.SaveToTextFile("C:\\temp\\raw_Mxx.txt");
                //ortho_Mxx.SaveToTextFile("C:\\temp\\ortho_Mxx.txt");
            }
        }

        /// <summary>
        /// If true, the orthonormalization is restarted when the maximum Krylov dimension is reached;
        /// Otherwise, Krylov dimension prevented to grow by just removing the least significant solution vector.
        /// </summary>
        public bool Restarted {
            get;
            set;
        }
        
        /// <summary>
        /// Maximum number of FlexGMRES iterations.
        /// </summary>
        public int MaxIter = 100;

        /// <summary>
        /// Threshold for convergence detection
        /// </summary>
        public double Tolerance = 1E-10;

        public int MaxKrylovDim = 80;



        bool m_Converged = false;
        int m_ThisLevelIterations = 0;

        public int IterationsInNested {
            get {
                if(this.PrecondS != null)
                    return this.PrecondS.Sum(pc => pc.IterationsInNested + pc.ThisLevelIterations);
                else
                    return 0;
            }
        }

        public int ThisLevelIterations {
            get { return this.m_ThisLevelIterations; }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public void ResetStat() {
            m_Converged = false;
            m_ThisLevelIterations = 0;
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }
        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }
    }
}
