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
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using MPI.Wrappers;
using System.Numerics;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;
using System.IO;
using System.Runtime.Serialization;
using BoSSS.Solution.Control;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// 
    /// </summary>
    /// <remarks>
    /// First described in: 
    ///  Vinsome, P.K.W.: 
	///  Orthomin, an Iterative Method for Solving Sparse Sets of Simultaneous Linear Equations,
    ///  SPE Symposium on Numerical Simulation of Reservoir Performance, 
    ///  doi: 10.2118/5729-MS, Los Angeles, California, 1976
    /// </remarks>
    class CoreOrthonormalizationProcedure {


        public CoreOrthonormalizationProcedure(ISparseMatrix __OpMatrix) {
            MxxHistory = new List<double[]>();
            SolHistory = new List<double[]>();
            Alphas = new List<(double, double, string)>();
            OpMatrix = __OpMatrix;
        }

        public void Clear() {
            MxxHistory.Clear();
            SolHistory.Clear();
            Alphas.Clear();
        }


        public long GetMemUsage() {
            int SizeSol = SolHistory != null && SolHistory.Count() > 0 ? this.SolHistory.Count() * this.SolHistory[0].Length * sizeof(double) : 0;
            int SizeMxx = MxxHistory != null && MxxHistory.Count() > 0 ? this.MxxHistory.Count() * this.MxxHistory[0].Length * sizeof(double) : 0;
            int SizeAlpha = Alphas != null && Alphas.Count() > 0 ? this.Alphas.Count() * (sizeof(double) * 2 + sizeof(int)) : 0;
            return SizeAlpha + SizeMxx + SizeSol;
        }


        /// <summary>
        /// 
        /// </summary>
        public ISparseMatrix OpMatrix {
            get;
        }


        /// <summary>
        /// solution guesses from smoothers
        /// </summary>
        public List<double[]> SolHistory = new List<double[]>();

        /// <summary>
        /// - orthonormal system of matrix-vector products;
        /// - the i-th entry is equal to  <see cref="OpMatrix"/>*<see cref="SolHistory"/>[i]
        /// </summary>
        public List<double[]> MxxHistory = new List<double[]>();


        /// <summary>
        /// scaling factors which were applied to <see cref="SolHistory"/> to approximate the solution
        /// - 1st item: the scaling applied to the respective solutions provided through <see cref="AddSol"/>
        /// - 2nd: relative residual reduction (typically greater or equal to 1)
        /// - 3rd: user id data
        /// </summary>
        public List<(double alpha_i, double RelResReduction, string id)> Alphas = new List<(double, double, string)>();


        /// <summary>
        /// 
        /// </summary>
        public double FullMinimization(double[] outX, double[] Sol0, double[] Res0, double[] outRes) {
            //var ids = Alphas.Select(ttt => ttt.id);
            Alphas.Clear();
            double LastResi = -1;
            while(Alphas.Count < SolHistory.Count) {
                LastResi = MinimizeResidual(outX, Sol0, Res0, outRes, "lulu");
                //ids = ids.Skip(1);
            }
            return LastResi;
        }


        /// <summary>
        /// Residual minimization in the space spanned by <see cref="SolHistory"/>/<see cref="MxxHistory"/>
        /// </summary>
        /// <param name="outX">
        /// - input: the solution build upon all vectors in the Krylov space **except the last one**
        /// - output: on exit, updated by the contribution from the most recent Krylov basis
        /// </param>
        /// <param name="Res0">
        /// residual with respect to the initial guess
        /// </param>
        /// <param name="Sol0">
        /// initial guess; required for checking <paramref name="outRes"/>,
        /// resp. for re-computing <paramref name="outRes"/> in the case of severe round-of errors.
        /// </param>
        /// <param name="outRes">
        /// - input: the residual with respect to the input value of <paramref name="outX"/>
        /// - output: on exit, the residual with respect to the output value of <paramref name="outX"/>
        /// </param>
        /// <param name="id">
        /// name/information for debugging and analysis
        /// </param>
        private double MinimizeResidual(double[] outX, double[] Sol0, double[] Res0, double[] outRes, string id) {
            using (var ft = new FuncTrace()) {
                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(outX.Length == OpMatrix.ColPartition.LocalLength);
                //Debug.Assert(Sol0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(Res0.Length == OpMatrix.RowPartitioning.LocalLength);
                Debug.Assert(outRes.Length == OpMatrix.RowPartitioning.LocalLength);



                int KrylovDim = SolHistory.Count;
                int L = outX.Length;

                //if (Alphas.Count != KrylovDim - 1) {
                //    throw new ApplicationException();
                //}

                // note: this implementation is based on the assumption, that the 
                // factors alpha_0 ... alpha_(i-1) are equal to the previous call;
                // therefore, it is only necessary to apply the contributions from the most recent Krylov vector

                var oldResiNorm = Norm(outRes);

                int i = KrylovDim - 1;
                double alpha_i = InnerProduct(MxxHistory[i], Res0);
                BLAS.daxpy(L, alpha_i, SolHistory[i], 1, outX, 1); // accumulate solution correction...
                BLAS.daxpy(L, -alpha_i, MxxHistory[i], 1, outRes, 1); // ...and its effect on the residual

                double ResNorm = Norm(outRes);
                double tol = Math.Max(ResNorm, oldResiNorm) * 0.2;
                tol = Math.Max(tol, 1.0e-7);
                if (ResNorm > oldResiNorm + tol)
                    throw new ArithmeticException($"residual increase (L = {L}): old norm: {oldResiNorm}, new norm {ResNorm}, reduction factor: {ResNorm / oldResiNorm}");


                if (m_UsedMoreThanOneOrthoCycle || KrylovDim % 40 == 0) { // 40 should be a decent compromise between performance and stability
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // once in a while, re-compute the residual to get rid of round-off errors which might accumulate
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                    outRes.SetV(Res0);
                    Sol0.AccV(-1.0, outX);
                    this.OpMatrix.SpMV(1.0, Sol0, 1.0, outRes);
                    Sol0.AccV(+1.0, outX);
                }



                /*
                if (m_MgOperator is MultigridOperator mgop && mgop.LevelIndex == 0) {
                    // prevent to much info dropping from lower levels
                    double xNorm = Norm(outX);
                    ft.Info("|x|: " + xNorm + ", last alpha was " + alpha_i);

                    //var iterSol = mgop.ProlongateRhsToDg(outX, "Sol");
                    //Tecplot.Tecplot.PlotFields(iterSol, "itersol-" + cnt, 0, 2);
                    Console.WriteLine("residual: " + cnt + "    " + ResNorm);
                    cnt++;
                }
                */


                /*
                // -------------------
                var rjCheck = Res0.CloneAs();
                OpMatrix.SpMV(-1.0, outX, 1.0, rjCheck);
                double dist = rjCheck.L2Dist(outRes);

                
                var OrthoCheck = MultidimensionalArray.Create(i + 1, i + 1);
                for(int k = 0; k <= i; k++) {
                    for(int w = 0; w <= i; w++) {
                        OrthoCheck[k, w] = MxxHistory[k].InnerProd(MxxHistory[w]).MPISum();
                    }
                }
                OrthoCheck.AccEye(-1.0);
                var orthoErr = OrthoCheck.InfNorm();

                double newResiNorm = outRes.L2Norm();
                Debug.Assert(dist <= Math.Max(outRes.L2Norm(), rjCheck.L2Norm()) * 1e-10);
                Debug.Assert(newResiNorm < oldResiNorm * 1.0001, "residual should shrink");
                Debug.Assert(orthoErr < 1.0e-7);
                // ------------------------
                */




                Alphas.Add((alpha_i, oldResiNorm / ResNorm, id));
                return ResNorm;
            }
        }

 

        bool m_UsedMoreThanOneOrthoCycle = false;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="X"></param>
        /// <param name="name"></param>
        private double AddSol(ref double[] X, string name) {
            using (var ft = new FuncTrace()) {

                m_UsedMoreThanOneOrthoCycle = false;

                double FillXwithRandom(double[] __X, double[] __Mxx) {
                    double __NormMxx;
                    ft.Error("Solution norm is exactly 0.0 -- OR -- re-orthonormalization failed; trying with a random vector instead to recover.");
                    __X.FillRandom();

                    __Mxx.ClearEntries();
                    OpMatrix.SpMV(1.0, __X, 0.0, __Mxx);
                    __NormMxx = Norm(__Mxx);

                    if (__NormMxx == 0)
                        throw new ArithmeticException("Numerical breakdown: norm of matrix * solution after using RANDOM NUMBERS (!) is " + __NormMxx);
                    if (__NormMxx.IsNaNorInf())
                        throw new ArithmeticException("Numerical breakdown: norm of matrix * solution after using RANDOM NUMBERS (!) is " + __NormMxx);
                    return __NormMxx;
                }

                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(X.Length == OpMatrix.RowPartitioning.LocalLength);
                int L = X.Length;
                int KrylovDim = SolHistory.Count;


                double[] Mxx = new double[L];
                OpMatrix.SpMV(1.0, X, 0.0, Mxx);
                double NormMxx = Norm(Mxx);
                if (NormMxx.IsNaNorInf())
                    throw new ArithmeticException("Numerical breakdown: norm of matrix * solution is " + NormMxx);

                if (NormMxx == 0) {
                    // solution is really strange: try with a random vector.
                    NormMxx = FillXwithRandom(X, Mxx);
                }

                // scale Mxx to norm 1, should allow better stability when compared to the Krylov vectors who have all norm 1;
                BLAS.dscal(L, 1.0 / NormMxx, Mxx, 1);
                BLAS.dscal(L, 1.0 / NormMxx, X, 1);

                const int MaxOrtho = 10;
                double MxxRemainderAfterOrthoNorm = 0;
                for (int jj = 0; jj <= MaxOrtho; jj++) { // re-orthogonalisation, loop-limit to 10; See also book of Saad, p 156, section 6.3.2

                    m_UsedMoreThanOneOrthoCycle = jj > 0;

                    for (int i = 0; i < KrylovDim; i++) {
                        Debug.Assert(!object.ReferenceEquals(Mxx, MxxHistory[i]));
                        double beta = InnerProduct(Mxx, MxxHistory[i]);
                        BLAS.daxpy(L, -beta, SolHistory[i], 1, X, 1);
                        BLAS.daxpy(L, -beta, MxxHistory[i], 1, Mxx, 1);
                    }

                    //double NormAfter = Norm(Mxx);
                    //Console.WriteLine("   orthonormalization norm reduction: " + (NormAfter/NormInitial));

                    double NewMxxNorm = Norm(Mxx);
                    if (jj == 0)
                        MxxRemainderAfterOrthoNorm = NewMxxNorm;

                    if (NewMxxNorm <= 1E-5) {
                        using (new BlockTrace("re-orthonormalization", ft)) {

                            if (jj < MaxOrtho) {

                                // a lot of canceling out has occurred; |Mxx| dropped by several magnitudes.
                                // do another loop, to ensure ortho-normality:
                                // Mxx and X lost more than 5 Magnitudes, so we might have lost a couple of digits of accuracy

                                // Note: although |Mxx| might be small, |X| can be quite large (and probably random).
                                // To ensure stability, we must start over with a re-scaled X!
                                // We have to re-scale what is remaining of X:
                                double Xnorm = Norm(X);
                                Console.Error.WriteLine("Orthonormalization: Severe cancellation may have occurred after " + name + ", norm after orthogonalization is " + NewMxxNorm + "; norm of X: " + Xnorm + " Doing Re-orthonormalization (" + (jj + 1) + "); L = " + X.Length);
                                ft.Info("Orthonormalization: Severe cancellation may have occurred after " + name + ", norm after orthogonalization is " + NewMxxNorm + "; norm of X: " + Xnorm + " Doing Re-orthonormalization (" + (jj + 1) + "); L = " + X.Length);

                                if ((Xnorm < 1e-200) || (jj == MaxOrtho - 1)) {
                                    // prohibits div by 0, if we got zero solution  
                                    NormMxx = FillXwithRandom(X, Mxx);
                                } else {
                                    X.ScaleV(1.0 / Xnorm);
                                    Mxx.ClearEntries();
                                    OpMatrix.SpMV(1.0, X, 0.0, Mxx);
                                    NormMxx = Norm(Mxx);

                                    if (NormMxx == 0)
                                        // solution is really strange: try with a random vector.
                                        NormMxx = FillXwithRandom(X, Mxx);
                                }

                                double gamma = 1 / NewMxxNorm;
                                BLAS.dscal(L, gamma, Mxx, 1);
                                BLAS.dscal(L, gamma, X, 1);
                            } else {
                                throw new ArithmeticException("numerical breakdown");
                            }
                        }
                    } else {
                        double gamma = 1 / NewMxxNorm;
                        BLAS.dscal(L, gamma, Mxx, 1);
                        BLAS.dscal(L, gamma, X, 1);
                        break; // normally, we should terminate after 1 orthonormalization cycle
                    }
                }


                SolHistory.Add(X);
                MxxHistory.Add(Mxx);
                X = null;

                return MxxRemainderAfterOrthoNorm;
            }
        }

        MPI_Comm comm {
            get {
                return OpMatrix.MPI_Comm;
            }
        }

        public double[] Diagonal {
            get;
            set;
        }


        /// <summary>
        /// Inner Product use for orthonormalization (not necessarily the standard l2-product)
        /// </summary>
        public double InnerProduct(double[] A, double[] B) {
            if (Diagonal == null) {
                Debug.Assert(A.Length == B.Length);
                return BLAS.ddot(A.Length, A, 1, B, 1).MPISum(comm);
            } else {
                int L = Diagonal.Length; 
                Debug.Assert(A.Length == L);
                Debug.Assert(B.Length == L);

                double acc = 0;
                for(int i = 0; i < L; i++) {
                    if (Diagonal[i] <= 0)
                        throw new ArithmeticException("fuced up diagonal");
                    acc += A[i] * B[i] * Diagonal[i];
                }
                return acc.MPISum(comm);
            }
        }

        /// <summary>
        /// norm induced by <see cref="InnerProduct(double[], double[])"/>
        /// </summary>
        public double Norm(double[] A) {
            if (Diagonal == null) {
                return A.MPI_L2Norm(comm);
            } else {
                return InnerProduct(A, A).Sqrt();
            }
        }


        /// <summary>
        /// Residual minimization in the space spanned by <see cref="SolHistory"/>/<see cref="MxxHistory"/>
        /// </summary>
        /// <param name="outX">
        /// - input: the solution build upon all vectors in the Krylov space **except the last one**
        /// - output: on exit, updated by the contribution from the most recent Krylov basis
        /// </param>
        /// <param name="Res0">
        /// residual with respect to the initial guess
        /// </param>
        /// <param name="Sol0">
        /// initial guess; required for checking <paramref name="outRes"/>,
        /// resp. for re-computing <paramref name="outRes"/> in the case of severe round-of errors.
        /// </param>
        /// <param name="outRes">
        /// - input: the residual with respect to the input value of <paramref name="outX"/>
        /// - output: on exit, the residual with respect to the output value of <paramref name="outX"/>
        /// </param>
        /// <param name="name">
        /// identifier 
        /// </param>
        /// <param name="xGuess">
        /// guess for the solution (typically, preconditioner output); will be modified during the orthonormalization
        /// </param>
        public double AddSolAndMinimizeResidual(ref double[] xGuess, double[] outX, double[] Sol0, double[] Res0, double[] outRes, string name) {            
            double UsablePart = AddSol(ref xGuess, name);
            double newResNorm = MinimizeResidual(outX, Sol0, Res0, outRes, name);

            //if(name != null)
            //    Console.WriteLine($"   {name} \t\t\t{UsablePart:0.####e-00}\t{Alphas.Last().RelResReduction:0.####e-00}\t{newResNorm:0.####e-00}");


            return newResNorm;
        }
    }
    /*
    class DecompositionAccelerator {

        CoreOrthonormalizationProcedure m_kSpace;

        public DecompositionAccelerator(MultigridOperator op, CoreOrthonormalizationProcedure kSpace) {
            m_kSpace = kSpace;
            
            int pLowVel = 4;
            int pLowPres = pLowVel - 1;

            SubBlockSelector VelLo = new SubBlockSelector(op.DgMapping);
            VelLo.VariableSelector(0, 1);
            VelLo.SetModeSelector(k => k <= pLowVel);
            BlockMask maskVelLo = new BlockMask(VelLo);

            SubBlockSelector VelHi = new SubBlockSelector(op.DgMapping);
            VelHi.VariableSelector(0, 1);
            VelHi.SetModeSelector(k => k > pLowVel);
            BlockMask maskVelHi = new BlockMask(VelHi);

            SubBlockSelector PresLo = new SubBlockSelector(op.DgMapping);
            PresLo.VariableSelector(2);
            PresLo.SetModeSelector(k => k <= pLowPres);
            BlockMask maskPresLo = new BlockMask(PresLo);

            SubBlockSelector PresHi = new SubBlockSelector(op.DgMapping);
            PresHi.VariableSelector(2);
            PresHi.SetModeSelector(k => k > pLowPres);
            BlockMask maskPresHi = new BlockMask(PresHi);

            int L = op.Mapping.LocalLength;


            bool[] check = new bool[L];
            foreach (var msk in new BlockMask[] { maskVelLo, maskVelHi, maskPresLo, maskPresHi }) {
                foreach (int i in msk.LocalIndices) {
                    if (check[i] == true)
                        throw new Exception("double hit " + i);
                    check[i] = true;
                }
            }
            for (int i = 0; i < L; i++) {
                if (check[i] == false)
                    throw new Exception("miss " + i);

            }
        
            m_Decompositions.Add(maskVelLo.LocalIndices);
            m_Decompositions.Add(maskVelHi.LocalIndices);
            m_Decompositions.Add(maskPresLo.LocalIndices);
            m_Decompositions.Add(maskPresHi.LocalIndices);
        }

        List<int[]> m_Decompositions = new List<int[]>();                                                   


        public void Solve(double[] X, double[] RHS) {

            var X0 = X.CloneAs();
            var Res = RHS.CloneAs();

            int L = X.Length;
            X.ClearEntries();

            var myOrtho = new CoreOrthonormalizationProcedure(m_kSpace.OpMatrix);

            //myOrtho.AddSolAndMinimizeResidual()

            foreach(var Y in m_kSpace.SolHistory) {
                int i = -1;
                foreach(int[] Decomp in m_Decompositions) {
                    i++;

                    double[] Yd = new double[L];
                    Yd.AccV(1.0, Y, Decomp, Decomp);
                    if (Yd.MPI_L2Norm(m_kSpace.OpMatrix.MPI_Comm) <= 0)
                        continue;
                    double resNorm = myOrtho.AddSolAndMinimizeResidual(ref Yd, X, X0, RHS, Res, "decomp" + i);
                    
                    //Console.WriteLine("decomp" + i + ": " + resNorm);
                }
            }

        }


    }
    */


    /// <summary>
    /// A recursive multigrid method, where the convergence (i.e. non-divergence) of each mesh level is guaranteed by an
    /// orthonormalization approach, similar to flexible GMRES, 
    /// as described in:
    ///   BoSSS: A package for multigrid extended discontinuous Galerkin methods; Kummer, Florian; Weber, Jens; Smuda, Martin; Computers &amp; Mathematics with Applications 81
    ///   see https://www.sciencedirect.com/science/article/abs/pii/S0898122120301917?via%3Dihub
    ///   see also https://tubiblio.ulb.tu-darmstadt.de/121465/
    ///   
    /// One instance of this class represents only one level of the multigrid method; coarser 
    /// levels must be added by configuring the <see cref="OrthonormalizationMultigrid.CoarserLevelSolver"/> member.
    /// </summary>
    /// <remarks>
    /// Residual minimization through orthonormalization (ORTHOMIN) is first described in:
    ///  Vinsome, P.K.W.: 
    ///  Orthomin, an Iterative Method for Solving Sparse Sets of Simultaneous Linear Equations,
    ///  SPE Symposium on Numerical Simulation of Reservoir Performance, 
    ///  doi: 10.2118/5729-MS, Los Angeles, California, 1976
    /// </remarks>    
    public class OrthonormalizationMultigrid : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination, ISubsystemSolver {


        /// <summary>
        /// Individual configuration of <see cref="OrthonormalizationMultigrid"/>
        /// </summary>
        [DataContract]
        [Serializable]
        public class Config : IterativeSolverConfig {


            /// <summary>
            /// Ctor
            /// </summary>
            public Config() {
            }



            /// <summary>
            /// config cctor of <see cref="OrthonormalizationMultigrid"/>
            /// </summary>
            public int MaxKrylovDim = int.MaxValue;

            /// <summary>
            /// - True: the default value: <see cref="OrthonormalizationMultigrid.CoarserLevelSolver"/> is initialized and solved on coarser level
            /// - false: <see cref="OrthonormalizationMultigrid.CoarserLevelSolver"/> is initialized on the same level, but it may perform tis own restriction
            /// </summary>
            [DataMember]
            public bool CoarseOnLovwerLevel = true;

            /// <summary>
            /// - if set to 1, a this performs a V-cycle
            /// - if set to 2 or higher, a W-cycle, resp. WW-cycle, WWW-cycle, etc.
            /// </summary>
            [DataMember]
            public int m_omega = 1;

            /// <summary>
            /// 
            /// </summary>
            override public string Name {
                get { return "OrthonormalizationMultigrid"; }
            }

            /// <summary>
            /// 
            /// </summary>
            override public string Shortname {
                get { return "OrthMG"; }
            }

            /// <summary>
            /// factory
            /// </summary>
            override public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var instance = new OrthonormalizationMultigrid();
                instance.myConfig = this;
                instance.Init(level);
                instance.TerminationCriterion = base.DefaultTermination;
                return instance;
            }

            /// <summary>
            /// 
            /// </summary>
            public int NoOfPostSmootherSweeps = 2;
        }


        Config myConfig;

        /// <summary>
        /// Solver configuration
        /// </summary>
        public Config config {
            get {
                return myConfig;
            }
        }

        static int counter_objectId = 1;

        public int objectId;

        /// <summary>
        /// ctor
        /// </summary>
        public OrthonormalizationMultigrid() {
            myConfig = new Config();
            TerminationCriterion = myConfig.DefaultTermination;
            objectId = counter_objectId;
            counter_objectId++;
        }

        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
        }


        /// <summary>
        /// The matrix at this level.
        /// </summary>
        public BlockMsrMatrix OpMatrix => m_MgOperator.OperatorMatrix;

        /// <summary>
        /// passed in <see cref="InitImpl"/>
        /// </summary>
        IOperatorMappingPair m_MgOperator;


        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(MultigridOperator op) {
            InitImpl(op);

        }

        /*
        public void ScalingKakke(MultigridOperator op) { 
            int pLowVel = 4;
            int pLowPres = pLowVel - 1;

            SubBlockSelector VelLo = new SubBlockSelector(op.DgMapping);
            VelLo.SetVariableSelector(0, 1);
            VelLo.SetModeSelector(k => k <= pLowVel);
            BlockMask maskVelLo = new BlockMask(VelLo);

            SubBlockSelector VelHi = new SubBlockSelector(op.DgMapping);
            VelHi.SetVariableSelector(0, 1);
            VelHi.SetModeSelector(k => k > pLowVel);
            BlockMask maskVelHi = new BlockMask(VelHi);
            
            SubBlockSelector PresLo = new SubBlockSelector(op.DgMapping);
            PresLo.SetVariableSelector(2);
            PresLo.SetModeSelector(k => k <= pLowPres);
            BlockMask maskPresLo = new BlockMask(PresLo);

            SubBlockSelector PresHi = new SubBlockSelector(op.DgMapping);
            PresHi.SetVariableSelector(2);
            PresHi.SetModeSelector(k => k > pLowPres);
            BlockMask maskPresHi = new BlockMask(PresHi);

            int L = op.Mapping.LocalLength;
            ortho.Diagonal = new double[L];
            ortho.Diagonal.SetAll(1.0);


            bool[] check = new bool[L];
            foreach (var msk in new BlockMask[] { maskVelLo, maskVelHi, maskPresLo, maskPresHi }) {
                foreach (int i in msk.LocalIndices) {
                    if (check[i] == true)
                        throw new Exception("double hit " + i);
                    check[i] = true;
                }
            }
            for(int i = 0; i < L; i++) {
                if (check[i] == false)
                    throw new Exception("miss " + i);

            }

            Console.WriteLine($"    Inner Product diagonal factors: {velScaleLo}  {velScaleHi}  {presScaleLo}  {presScaleHi}");
            
            ortho.Diagonal.ScaleV(velScaleLo, maskVelLo.LocalIndices);
            ortho.Diagonal.ScaleV(velScaleHi, maskVelHi.LocalIndices);
            ortho.Diagonal.ScaleV(presScaleLo, maskPresLo.LocalIndices);
            ortho.Diagonal.ScaleV(presScaleHi, maskPresHi.LocalIndices);

        

            //Console.WriteLine("     vL: " + maskVelLo.GetSubVec(vec).L2Norm());
            //Console.WriteLine("     vH: " + maskVelHi.GetSubVec(vec).L2Norm());
            //Console.WriteLine("     pL: " + maskPresLo.GetSubVec(vec).L2Norm());
            //Console.WriteLine("     pH: " + maskPresHi.GetSubVec(vec).L2Norm());


            Console.WriteLine("check passed");
        }

        
        public static double velScaleLo = 1;
        public static double presScaleLo = 1;
        public static double velScaleHi = 1;
        public static double presScaleHi = 1;
        //*/


        CoreOrthonormalizationProcedure ortho;

        void InitImpl(IOperatorMappingPair op) {
            using (var tr = new FuncTrace()) {
                if (object.ReferenceEquals(op, m_MgOperator))
                    return; // already initialized
                else
                    this.Dispose();


                this.m_MgOperator = op;
                var Mtx = op.OperatorMatrix;
                var MgMap = op.DgMapping;
                
                if (!Mtx.RowPartitioning.EqualsPartition(MgMap))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mtx.ColPartition.EqualsPartition(MgMap))
                    throw new ArgumentException("Column partitioning mismatch.");



                ortho = new CoreOrthonormalizationProcedure(Mtx);


                // set operator
                // ============

                // initiate coarser level
                // ======================
                if (this.CoarserLevelSolver == null) {
                    //throw new NotSupportedException("Missing coarse level solver.");
                    Console.WriteLine("OrthonormalizationMultigrid: running without coarse solver.");
                } else {
                    if (op is MultigridOperator mgOp) {
                        if (myConfig.CoarseOnLovwerLevel && mgOp.CoarserLevel != null) {
                            this.CoarserLevelSolver.Init(mgOp.CoarserLevel);
                        } else {
                            Console.WriteLine("OrthonormalizationMultigrid: running coarse solver on same level.");
                            this.CoarserLevelSolver.Init(mgOp);
                        }
                    } else {
                        if(myConfig.CoarseOnLovwerLevel == false && this.CoarserLevelSolver is ISubsystemSolver ssCoarse) {
                            ssCoarse.Init(op);
                        } else {
                            throw new NotSupportedException($"Unable to initialize coarse-level-solver if operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }

                // init smoother
                // =============
                if (PreSmoother != null) {
                    if (PreSmoother is ISubsystemSolver ssPreSmother) {
                        ssPreSmother.Init(op);
                    } else {
                        if (op is MultigridOperator mgOp) {
                            PreSmoother.Init(mgOp);
                        } else {
                            throw new NotSupportedException($"Unable to initialize pre-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }
                if (PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother)) {
                    //PostSmoother.Init(op);
                    if (PostSmoother is ISubsystemSolver ssPostSmother) {
                        ssPostSmother.Init(op);
                    } else {
                        if (op is MultigridOperator mgOp) {
                            PostSmoother.Init(mgOp);
                        } else {
                            throw new NotSupportedException($"Unable to initialize post-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }
                if(AdditionalPostSmoothers != null) {
                    foreach (var ps in AdditionalPostSmoothers) {
                        if (ps is ISubsystemSolver ssPostSmother) {
                            ssPostSmother.Init(op);
                        } else {
                            if (op is MultigridOperator mgOp) {
                                ps.Init(mgOp);
                            } else {
                                throw new NotSupportedException($"Unable to initialize post-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// coarse-level correction; can be defined either
        /// - on this level (then the coarse solver may perform its of prolongation/restriction), or
        /// - on coarser level, then prolongation/restriction is handled by this solver.
        /// </summary>
        public ISolverSmootherTemplate CoarserLevelSolver;

        /// <summary>
        /// high frequency solver before coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate PreSmoother;

        /// <summary>
        /// high frequency solver after coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate PostSmoother;

        /// <summary>
        /// high frequency solver after coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate[] AdditionalPostSmoothers;



        /// <summary>
        /// 
        /// </summary>
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        /// <summary>
        /// computes the residual on this level.
        /// </summary>
        /// <param name="B">input: RHS of the system</param>
        /// <param name="X">input: solution approximation</param>
        /// <param name="Res">output: on exit <paramref name="B"/> - <see cref="OpMatrix"/>*<paramref name="X"/></param>
        public void Residual(double[] Res, double[] X, double[] B) {
            using (new FuncTrace()) {
                Debug.Assert(Res.Length == m_MgOperator.DgMapping.LocalLength);
                Debug.Assert(X.Length == m_MgOperator.DgMapping.LocalLength);
                Debug.Assert(B.Length == m_MgOperator.DgMapping.LocalLength);
                int L = Res.Length;
                //Res.SetV(B);
                Array.Copy(B, Res, L);
                OpMatrix.SpMV(-1.0, X, 1.0, Res);
                //Res.AccV(1.0, B);


            }
        }

        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="_xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="_B">the right-hand-side of the problem</param>
        public void Solve<U, V>(U _xl, V _B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (var f = new FuncTrace()) {
                double[] B, X;
                if (_B is double[])
                    B = _B as double[];
                else
                    B = _B.ToArray();
                if (_xl is double[])
                    X = _xl as double[];
                else
                    X = _xl.ToArray();

                int L = X.Length;
                int Lc;
                if (this.CoarserLevelSolver != null && myConfig.CoarseOnLovwerLevel)
                    Lc = ((MultigridOperator)m_MgOperator).CoarserLevel.Mapping.LocalLength;
                else
                    Lc = -1;

                double[] Res = new double[L];
                double[] ResCoarse = new double[Lc];


                double[] Sol0 = X.CloneAs();
                double[] Res0 = new double[L];
                Residual(Res0, Sol0, B);
                Array.Copy(Res0, Res, L);

                bool isTopLevel = ((m_MgOperator as MultigridOperator)?.LevelIndex ?? -1) == 0;

                double iter0_resNorm = ortho.Norm(Res0);
                double resNorm = iter0_resNorm;
                this.IterationCallback?.Invoke(0, Sol0, Res0, this.m_MgOperator as MultigridOperator);
                                
                // clear history of coarse solvers
                ortho.Clear();
                bool bIterate = true;

                for (int iIter = 1; bIterate; iIter++) {
                    var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if (!termState.bNotTerminate) {
                        Converged = termState.bSuccess;
                        break;
                    } else {

                    }

                    // pre-smoother
                    // ------------

                    {
                        VerivyCurrentResidual(X, B, Res, iIter);


                        // compute correction
                        //var oldRl = rl.CloneAs();

                        if (PreSmoother != null) {
                            double[] PreCorr = new double[L];
                            PreSmoother.Solve(PreCorr, Res); // Vorglättung

                            // orthonormalization and residual minimization
                            resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorr, X, Sol0, Res0, Res, isTopLevel ? "presmooth" : null);

                        }


                        //SpecAnalysisSample(iIter, X, "ortho1");
                        var termState2 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                        if (!termState2.bNotTerminate) {
                            Converged = termState2.bSuccess;
                            break;
                        }
                    }

                    // coarse grid correction
                    // ----------------------
                    // Test: Residual on this level / already computed by 'MinimizeResidual' above
                    VerivyCurrentResidual(X, B, Res, iIter);

                    double resNorm_b4Coarse = resNorm;
                    for (int i = 0; i < myConfig.m_omega; i++) {
                        if (this.CoarserLevelSolver != null) {

                            double[] vl = new double[L];
                            if (myConfig.CoarseOnLovwerLevel) {

                                var _MgOperator = m_MgOperator as MultigridOperator;

                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // coarse grid solver defined on COARSER MESH LEVEL:
                                // this solver must perform restriction and prolongation
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                using (new BlockTrace("Restriction", f)) {
                                    // restriction of residual
                                    _MgOperator.CoarserLevel.Restrict(Res, ResCoarse);
                                }
                                // Berechnung der Grobgitterkorrektur
                                double[] vlc = new double[Lc];
                                this.CoarserLevelSolver.Solve(vlc, ResCoarse);
                                using (new BlockTrace("Prolongation", f)) {
                                    // Prolongation der Grobgitterkorrektur
                                    _MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);
                                }

                            } else {
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // coarse grid solver defined on the SAME MESH LEVEL:
                                // performs (probably) its own restriction/prolongation
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                // Berechnung der Grobgitterkorrektur
                                this.CoarserLevelSolver.Solve(vl, Res);
                            }



                            // orthonormalization and residual minimization
                            //ortho.AddSol(ref vl, "coarsecor");
                            //resNorm = MinimizeResidual(X, Sol0, Res0, Res, 2);
                            resNorm = ortho.AddSolAndMinimizeResidual(ref vl, X, Sol0, Res0, Res, isTopLevel ? "coarsecor" : null);



                        }
                    } // end of coarse-solver loop

                    var termState3 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if (!termState3.bNotTerminate) {
                        Converged = termState3.bSuccess;
                        break;
                    }



                    // post-smoother
                    // -------------
                    if (PostSmoother != null || AdditionalPostSmoothers != null) {
                        bool termPost = false;

                        ISolverSmootherTemplate[] allSmooters;
                        {
                            allSmooters = new ISolverSmootherTemplate[(PostSmoother != null ? 1 : 0) + (AdditionalPostSmoothers?.Length ?? 0)];
                            if (PostSmoother != null)
                                allSmooters[allSmooters.Length - 1] = PostSmoother;
                            if (AdditionalPostSmoothers != null)
                                Array.Copy(AdditionalPostSmoothers, allSmooters, AdditionalPostSmoothers.Length);
                        }

                        for (int g = 0; g < config.NoOfPostSmootherSweeps; g++) { // Test: Residual on this level / already computed by 'MinimizeResidual' above

                            ISolverSmootherTemplate _PostSmoother = allSmooters[g % allSmooters.Length];

                            VerivyCurrentResidual(X, B, Res, iIter); // 

                            string id = "";
                            if (_PostSmoother is CellILU cilu)
                                id = cilu.id;

                            double[] PostCorr = new double[L];
                            //if (g % 2 == 1)
                            //    PreSmoother.Solve(PostCorr, Res);
                            //else

                            _PostSmoother.Solve(PostCorr, Res); // compute correction (Nachglättung)

                            if (PostCorr.ContainsForNanOrInfV()) {
                                Console.Error.WriteLine("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g);

                                if (Res.ContainsForNanOrInfV())
                                    Console.WriteLine("... so does RHS");
                                else
                                    Console.WriteLine("... although RHS is regular");

                                var viz = new MGViz(m_MgOperator as MultigridOperator);
                                var __RHS = viz.ProlongateRhsToDg(Res, "RES");
                                var __SOL = viz.ProlongateRhsToDg(PostCorr, "SOL");
                                Tecplot.Tecplot.PlotFields(__SOL.Cat(__RHS), "iilufail", 0.0, 0);
                                throw new ArithmeticException();


                            } else {
                                resNorm = ortho.AddSolAndMinimizeResidual(ref PostCorr, X, Sol0, Res0, Res, isTopLevel ? ("pstsmth" + id + "-" + g) : null);



                                var termState4 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                                if (!termState4.bNotTerminate) {
                                    Converged = termState4.bSuccess;
                                    termPost = true;
                                    break;
                                }
                            }
                        }
                        

                        if (termPost)
                            break;

                    } // end of post-smoother loop


                    
                    // iteration callback
                    // ------------------
                    this.ThisLevelIterations++;
                    IterationCallback?.Invoke(iIter, X, Res, this.m_MgOperator as MultigridOperator);


      

                    //ScalingKakke(this.m_MgOperator as MultigridOperator, Res);

                    //SpecAnalysisSample(iIter, X, "_");

                    

                } // end of solver iterations

                // solution copy
                // =============
                if (!ReferenceEquals(_xl, X)) {
                    _xl.SetV(X);
                }

            } // end of functrace
        }


      
        /// <summary>
        /// For performance optimization, the <see cref="OrthonormalizationMultigrid"/>
        /// assumes that <see cref="PreSmoother"/> and <see cref="PostSmoother"/>
        /// update the residual on exit.
        /// </summary>
        /// <param name="X"></param>
        /// <param name="B"></param>
        /// <param name="Res"></param>
        /// <param name="iter"></param>
        private void VerivyCurrentResidual(double[] X, double[] B, double[] Res, int iter) {
            using (var tr = new FuncTrace()) {
#if DEBUG
            {
#else
                if (iter % 20 == 0 && iter > 1) {
#endif
                    double[] rTest = new double[Res.Length];
                    Residual(rTest, X, B); // Residual on this level; 
                                           // Test also fails if convergence criterium is to strict because then machine accuracy is reached
                                           //Console.WriteLine("verified Residual: " + resDist);
                    var ss = (new double[] { rTest.L2DistPow2(Res), Res.L2NormPow2(), X.L2NormPow2() }).MPISum(m_MgOperator.DgMapping.MPI_Comm);
                    double resDist = ss[0].Sqrt();
                    double resNormTst = ss[1].Sqrt();
                    double XnormTest = ss[2].Sqrt();
                    tr.Info($"Residual vector check iter {iter}: distance is {resDist}, reference value {resNormTst}");
                    if (resDist > resNormTst * 10e-5 + XnormTest * 1e-5)
                        throw new ArithmeticException($"Residual vector (after pre-smoother/before coarse-correction) is not up-to-date: distance is {resDist}, reference value {resNormTst}");
                    //Debug.Assert(resDist <= resNormTst * 10e-5, $"Residual vector is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                }
            }
        }

        /*
        private double[] cloneofX;

        /// <summary>
        /// RealX is left unchanged, no worries
        /// </summary>
        /// <param name="iter"></param>
        /// <param name="RealX"></param>
        /// <param name="name"></param>
        private void SpecAnalysisSample(int iter, double[] RealX, string name) {

            var _mgOp = this.m_MgOperator as MultigridOperator;

            if (ExtractSamples == null || _mgOp.LevelIndex > 0)
                return; // zero overhead

            if (cloneofX == null)
                cloneofX = new double[RealX.Length];
            cloneofX.SetV(RealX);
            ExtractSamples.Invoke(iter, cloneofX, name);
        }

        /// <summary>
        /// gefrickel, could be integrated into IterationCallback
        /// without individual string of course
        /// </summary>
        public Action<int, double[], string> ExtractSamples {
            get;
            set;
        }
        */

        /// <summary>
        /// ~
        /// </summary>
        public int IterationsInNested {
            get {
                int iter = 0;

                iter += (this.PreSmoother?.IterationsInNested ?? 0) + (this.PreSmoother?.ThisLevelIterations ?? 0);

                if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                    iter += this.PostSmoother.IterationsInNested + this.PostSmoother.ThisLevelIterations;

                iter += (this.CoarserLevelSolver?.IterationsInNested ?? 0) + (this.CoarserLevelSolver?.ThisLevelIterations ?? 0);

                return iter;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public int ThisLevelIterations {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public bool Converged {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public void ResetStat() {
            this.Converged = false;
            this.ThisLevelIterations = 0;
            if (this.PreSmoother != null)
                this.PreSmoother.ResetStat();
            if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                this.PostSmoother.ResetStat();
            if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.ResetStat();
        }

        /// <summary>
        /// %
        /// </summary>
        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

        bool m_verbose = true;

        public long UsedMemory() {
            long Memory = 0;
            Memory += MemoryOfMultigrid();
            Memory += MemoryOfSmoother();
            return Memory;
        }

        public long MemoryOfSmoother() {
            long Memory = 0;
            if (this.CoarserLevelSolver is OrthonormalizationMultigrid)
                Memory += (this.CoarserLevelSolver as OrthonormalizationMultigrid).MemoryOfSmoother();
            if (PreSmoother != null) Memory += PreSmoother.UsedMemory();
            if (PostSmoother != null) Memory += PostSmoother.UsedMemory();
            return Memory;
        }

        public long MemoryOfMultigrid() {
            long Memory = 0;
            if (this.CoarserLevelSolver is OrthonormalizationMultigrid)
                Memory += (this.CoarserLevelSolver as OrthonormalizationMultigrid).MemoryOfMultigrid();

            Memory += ortho?.GetMemUsage() ?? 0;
            return Memory;
        }

        public void Dispose() {
            //if (m_MTracker != null) m_MTracker.Dispose();
            if (m_verbose && m_MgOperator != null && m_MgOperator is MultigridOperator _mgop && _mgop.LevelIndex == 0) {
                Console.WriteLine($"OrthoMG - total memory: {UsedMemory() / (1024 * 1024)} MB");
                Console.WriteLine($"OrthoMG - internal memory: {MemoryOfMultigrid() / (1024 * 1024)} MB");
                Console.WriteLine($"OrthoMG - smoother memory: {MemoryOfSmoother() / (1024 * 1024)} MB");
            }
            if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.Dispose();
            //this.CoarserLevelSolver = null; // don't delete - we need this again for the next init

            ortho?.Clear();
            ortho = null;
            this.m_MgOperator = null;

            if (PreSmoother != null)
                this.PreSmoother.Dispose();
            if (PostSmoother != null)
                this.PostSmoother.Dispose();
            //this.PreSmoother = null; // don't delete - we need this again for the next init
            //this.PostSmoother = null;  // don't delete - we need this again for the next init
        }


    }
}
