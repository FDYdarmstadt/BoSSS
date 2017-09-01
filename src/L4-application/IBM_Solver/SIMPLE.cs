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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP.Tracing;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.SpecFEM;
using MPI.Wrappers;
using BoSSS.Foundation.Quadrature;
using ilPSP.LinSolvers.PARDISO;
using System.Diagnostics;
using NUnit.Framework;
using System.IO;
using System.Globalization;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.Multigrid;
using ilPSP;
using ilPSP.Connectors.Matlab;

namespace BoSSS.Application.IBM_Solver {
    class SIMPLE : ISolverSmootherTemplate, ISolverWithCallback {


        MultigridOperator m_MgOp;

        public void Init(MultigridOperator op) {
            this.m_MgOp = op;
            int D = this.LsTrk.GridDat.SpatialDimension;
 
            int[] VelVarIdx = D.ForLoop(d => d);

            this.USubMatrixIdx_Row = this.m_MgOp.Mapping.GetSubvectorIndices(VelVarIdx);
            this.PSubMatrixIdx_Row = this.m_MgOp.Mapping.GetSubvectorIndices(new int[] { D });

            this.USpcSubMatrix_Row = new int[this.LsTrk.SpeciesIdS.Count][];
            this.PSpcSubMatrix_Row = new int[this.LsTrk.SpeciesIdS.Count][];
            for(int iSpc = 0; iSpc < this.LsTrk.SpeciesIdS.Count; iSpc++) {
                this.USpcSubMatrix_Row[iSpc] = this.m_MgOp.Mapping.GetSubvectorIndices(this.LsTrk.SpeciesIdS[iSpc], VelVarIdx);
                this.PSpcSubMatrix_Row[iSpc] = this.m_MgOp.Mapping.GetSubvectorIndices(this.LsTrk.SpeciesIdS[iSpc], D);
            }
                        
            this.ExtractMatrices();
            this.ApproximationMatrix();
        }

        /// <summary>
        /// maximum bumber of SIMPLE iteraions
        /// </summary>
        public int MaxIterations = 1;

        void DoCallBack<V>(double[] Velocity, double[] Pressure, V RHS)
            where V : IList<double> //
        {
            if(this.IterationCallback != null) {
                int LL = this.m_MgOp.Mapping.LocalLength;
                double[] sol = new double[LL];
                double[] res = new double[LL];

                res.SetV(RHS);
                sol.AccV(1.0, Velocity, USubMatrixIdx_Row, default(int[]));
                sol.AccV(1.0, Pressure, PSubMatrixIdx_Row, default(int[]));

                this.m_MgOp.OperatorMatrix.SpMV(-1.0, sol, 1.0, res);

                this.IterationCallback(this.NoOfIterations, sol, res, this.m_MgOp);
            }
        }

        /// <summary>
        /// performs multiple SIMPLE iterations, at max <see cref="MaxIterations"/>.
        /// </summary>
        /// <param name="X">input/output: solution guess</param>
        /// <param name="RHS">RHS of the saddle point problem</param>
        public void Solve<U, V>(U X, V RHS)
            where U : IList<double>
            where V : IList<double> //
        {
            // check some settings
            if(m_SIMPLEOptions.relax_p < 0.0)
                throw new ArithmeticException("Illegal pressure correction relaxation parameter: " + m_SIMPLEOptions.relax_p);
             
            // memalloc/split X and RHS in velocity/pressure resp. momentum/continuity components
            // ----------------------------------------------------------------------------------
            double[] RHSMomentum = new double[USubMatrixIdx_Row.Length];
            double[] RHSContinuity = new double[PSubMatrixIdx_Row.Length];
            RHS.GetSubVector(RHSMomentum, USubMatrixIdx_Row, default(int[]));
            RHS.GetSubVector(RHSContinuity, PSubMatrixIdx_Row, default(int[]));

            double[] Velocity = X.GetSubVector(USubMatrixIdx_Row, default(int[]));
            double[] Pressure = X.GetSubVector(PSubMatrixIdx_Row, default(int[]));

            double[] IntermediateVelocity = new double[Velocity.Length];
            double[] PressureCorrection = new double[Pressure.Length];

            this.DoCallBack(Velocity, Pressure, RHS);

            // SIMPLE(R) iterations
            // --------------------
            for(int iIter = 0; iIter < this.MaxIterations; iIter++) {
                // clone current velocity and pressure
                var VelocityIncrement = Velocity.CloneAs();
                var PressureIncrement = Pressure.CloneAs();

                // 
                double[] StaticPressure = null;
                if(this.m_SIMPLEOptions.RunSIMPLER == true) {
                    StaticPressure = new double[Pressure.Length];
                    
                    this.PressureSolver(StaticPressure, Velocity, RHSMomentum);

                    if(this.m_SIMPLEOptions.CorrrectionMode == SimpleCorrectionMode.classic) {
                        Pressure.ScaleV(1.0 - m_SIMPLEOptions.relax_p);
                        Pressure.AccV(m_SIMPLEOptions.relax_p, StaticPressure);
                    }

                }

                // predictor
                IntermediateVelocity.ClearEntries(); // sicherheitshalber
                VelocityPredictor(Pressure, Velocity, IntermediateVelocity, RHSMomentum);
                
                // pressure correction
                PCorrL2Norm = PressureCorrector(Pressure, IntermediateVelocity, PressureCorrection, RHSContinuity);

                if(this.m_SIMPLEOptions.CorrrectionMode == SimpleCorrectionMode.classic) {
                    VelocityUpdate(PressureCorrection, IntermediateVelocity, Velocity);
                    
                    if(this.m_SIMPLEOptions.RunSIMPLER == false) {
                        Pressure.AccV(m_SIMPLEOptions.relax_p, PressureCorrection);
                    }
                } else if(this.m_SIMPLEOptions.CorrrectionMode == SimpleCorrectionMode.ResidualMinimization || this.m_SIMPLEOptions.CorrrectionMode == SimpleCorrectionMode.ResidualMinimization_perSpecies) {
                    Minimizer(Velocity, Pressure, StaticPressure, IntermediateVelocity, PressureCorrection, RHS.ToArray());
                } else {
                    throw new NotSupportedException();
                }

                // evaluate 'Delta'
                VelocityIncrement.AccV(-1.0, Velocity);
                PressureIncrement.AccV(-1.0, Pressure);
                double VelocityIncrementNorm = VelocityIncrement.L2Norm();
                double PressureIncrementNorm = PressureIncrement.L2Norm();

                this.NoOfIterations++;
                this.DoCallBack(Velocity, Pressure, RHS);
            }

            X.ClearEntries();
            X.AccV(1.0, Velocity, USubMatrixIdx_Row, default(int[]));
            X.AccV(1.0, Pressure, PSubMatrixIdx_Row, default(int[]));
        }

        
        

        void Minimizer(
            double[] Velocity, double[] Pressure, // input: old solution; output: solution with 'optimally' reduced residual
            double[] StaticPressure, // optional input: pressure from additional SIMLER pressure computation.
            double[] VelocityPred, // velocity from predictor
            double[] PressureCorrection, 
            double[] b
            ) //
        {
            
            // collect all vectors from which the new solution will be constructed
            // ===================================================================

            List<double[]> Z = new List<double[]>(); // basis for the new solution
            List<double[]> Q = new List<double[]>(); // basis for the residual


            // add old solution
            recomb(Velocity, null, Z, Q);
            recomb(null, Pressure, Z, Q);

            // predictor velocity
            recomb(VelocityPred, null, Z, Q);

            // SIMPLER pressure
            if(StaticPressure != null)
                recomb(null, StaticPressure, Z, Q);

            // correction
            {
                double[] CorrectionVelocity = new double[Velocity.Length];
                double[] temp = new double[Velocity.Length];
                PressureGrad.SpMVpara(-1.0, PressureCorrection, 0.0, temp);
                AapproxInverse.SpMVpara(1.0, temp, 1.0, CorrectionVelocity);

                recomb(CorrectionVelocity, null, Z, Q);
                recomb(null, PressureCorrection, Z, Q);
            }

            // solve mimimization problem
            // ==========================

            Debug.Assert(Z.Count == Q.Count);
            int K = Z.Count;
            MultidimensionalArray lhsMinimi = MultidimensionalArray.Create(K, K);
            double[] rhsMinimi = new double[K];

            for(int k = 0; k < K; k++) {
                for(int l = 0; l < k; l++) {
                    lhsMinimi[k, l] = GenericBlas.InnerProd(Q[k], Q[l]).MPISum();
                    lhsMinimi[l, k] = lhsMinimi[k, l];
                }
                lhsMinimi[k, k] = GenericBlas.L2NormPow2(Q[k]).MPISum();
                rhsMinimi[k] = GenericBlas.InnerProd(Q[k], b).MPISum();
            }

            double[] alpha = new double[K];
            lhsMinimi.LeastSquareSolve(alpha, rhsMinimi); // use least-square solve in the case some 'Q'-vectors are linear dependent.

            // construct new solution
            // ======================

            Velocity.ScaleV(alpha[0]);
            Pressure.ScaleV(alpha[0]);
            for(int k = 1; k < K; k++) {
                Velocity.AccV(alpha[k], Z[k], default(int[]), this.USubMatrixIdx_Row);
                Pressure.AccV(alpha[k], Z[k], default(int[]), this.PSubMatrixIdx_Row);
            }
            
            
            /*

            int L = this.USubMatrixIdx_Row.Length + this.PSubMatrixIdx_Row.Length;

            double[] Res = b.CloneAs();
            {
                double[] X0 = new double[L];
                X0.AccV(1.0, Velocity, this.USubMatrixIdx_Row, default(int[]));
                X0.AccV(1.0, Pressure, this.PSubMatrixIdx_Row, default(int[]));

                this.m_MgOp.OperatorMatrix.SpMVpara(-1.0, X0, 1.0, Res);
            }



            // collect all vectors from which the new solution will be constructed
            // ===================================================================

            List<double[]> Z = new List<double[]>(); // basis for the new solution
            List<double[]> Q = new List<double[]>(); // basis for the residual


            // predictor velocity
            recomb(VelocityPred, null, Z, Q);

            // SIMPLER pressure
            if(StaticPressure != null)
                recomb(null, StaticPressure, Z, Q);


            // correction
            {
                double[] CorrectionVelocity = new double[Velocity.Length];
                double[] temp = new double[Velocity.Length];
                PressureGrad.SpMVpara(-1.0, PressureCorrection, 0.0, temp);
                AapproxInverse.SpMVpara(1.0, temp, 1.0, CorrectionVelocity);

                recomb(CorrectionVelocity, null, Z, Q);
                recomb(null, PressureCorrection, Z, Q);
                recomb(CorrectionVelocity, PressureCorrection, Z, Q);
            }

            // solve mimimization problem
            // ==========================

            Debug.Assert(Z.Count == Q.Count);
            int K = Z.Count;
            MultidimensionalArray lhsMinimi = MultidimensionalArray.Create(K, K);
            double[] rhsMinimi = new double[K];

            for(int k = 0; k < K; k++) {
                for(int l = 0; l < k; l++) {
                    lhsMinimi[k, l] = GenericBlas.InnerProd(Q[k], Q[l]).MPISum();
                    lhsMinimi[l, k] = lhsMinimi[k, l];
                }
                lhsMinimi[k, k] = GenericBlas.L2NormPow2(Q[k]).MPISum();
                rhsMinimi[k] = GenericBlas.InnerProd(Q[k], Res);
            }

            double[] alpha = new double[K];
            lhsMinimi.LeastSquareSolve(alpha, rhsMinimi);

            // construct new solution
            // ======================

            for(int k = 0; k < K; k++) {
                Velocity.AccV(alpha[k], Z[k], default(int[]), this.USubMatrixIdx_Row);
                Pressure.AccV(alpha[k], Z[k], default(int[]), this.PSubMatrixIdx_Row);
            }
         
             */
        }

        private void recomb(double[] Velocity, double[] Pressure, List<double[]> Z, List<double[]> Q) {
            int L = this.USubMatrixIdx_Row.Length + this.PSubMatrixIdx_Row.Length;

            if(this.m_SIMPLEOptions.CorrrectionMode == SimpleCorrectionMode.ResidualMinimization) {
                double[] Z1 = new double[L];
                if(Velocity != null)
                    Z1.AccV(1.0, Velocity, this.USubMatrixIdx_Row, default(int[]));
                if(Pressure != null)
                    Z1.AccV(1.0, Pressure, this.PSubMatrixIdx_Row, default(int[]));


                double[] Q1 = new double[L];
                this.m_MgOp.OperatorMatrix.SpMV(1.0, Z1, 0.0, Q1);

                Z.Add(Z1);
                Q.Add(Q1);
            } else if(this.m_SIMPLEOptions.CorrrectionMode == SimpleCorrectionMode.ResidualMinimization_perSpecies) {
                if(this.PSpcSubMatrix_Row.Length != this.USpcSubMatrix_Row.Length)
                    throw new ApplicationException();

                int NoOfSpc = PSpcSubMatrix_Row.Length;


                for(int iSpc = 0; iSpc < NoOfSpc; iSpc++) {
                    double[] Z4 = new double[L];

                    if(Velocity != null)
                        Z4.AccV(1.0, Velocity, this.USpcSubMatrix_Row[iSpc], default(int[]));
                    if(Pressure != null)
                        Z4.AccV(1.0, Pressure, this.PSpcSubMatrix_Row[iSpc], default(int[]));

                    double[] Q4 = new double[L];
                    this.m_MgOp.OperatorMatrix.SpMV(1.0, Z4, 0.0, Q4);

                    Z.Add(Z4);
                    Q.Add(Q4);
                }



            } else {
                throw new NotImplementedException();
            }
            

        }


        /// <summary>
        /// SIMPLE iteration counter
        /// </summary>
        private int NoOfIterations = 0;

        public int IterationsInNested {
            get {
                return 0;
            }
        }

        public int ThisLevelIterations {
            get {
                return this.NoOfIterations;
            }
        }

        public bool Converged {
            get {
                return m_Converged;
            }
        }

        bool m_Converged = false;

        public void ResetStat() {
            this.m_Converged = false;
            this.NoOfIterations = 0;
        }

        
        
        
        LevelSetTracker LsTrk;


        public SIMPLEOptions m_SIMPLEOptions = null;

        
        /// <summary>
        /// The inverse of <see cref="AapproxInverse"/>.
        /// </summary>
        MsrMatrix Aapprox;

        /// <summary>
        /// Approximation to <br/>
        /// \f[ 
        ///   \left( \frac{1}{\Delta t} M_{\textrm{mass}}  + M_{\textrm{C/D}} \right)^{-1}.
        /// \f]
        /// In the steady case, we assume $\Delta t = \infty$, therfor $1/\Delta t = 0$.
        /// </summary>
        MsrMatrix AapproxInverse;

        //private int[] UWholeMatrixIdx;
        //private int[] PWholeMatrixIdx;
        
        MsrMatrix VelocityDiv;
        MsrMatrix PressureGrad;
        MsrMatrix ConvDiff;
        MsrMatrix Stab;
        
        

        double PCorrL2Norm;
        private int[] PSubMatrixIdx_Row;
        private int[] USubMatrixIdx_Row;
        
        private int[][] USpcSubMatrix_Row;
        private int[][] PSpcSubMatrix_Row;

        //IDictionary<SpeciesId, IEnumerable<double>> Rho {
        //    get {
        //        int D = this.LsTrk.GridDat.SpatialDimension;
        //        double rho_A = 0, rho_B = 0;
        //        this.m_SIMPLEOptions.Control.confSolver.GetExtProperty("rho_A", true, ref rho_A);
        //        this.m_SIMPLEOptions.Control.confSolver.GetExtProperty("rho_B", true, ref rho_B);

        //        double[] _rho_A = new double[D];
        //        _rho_A.SetAll(rho_A);
        //        double[] _rho_B = new double[D];
        //        _rho_B.SetAll(rho_B);

        //        Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
        //        R.Add(this.LsTrk.GetSpeciesId("A"), _rho_A);
        //        R.Add(this.LsTrk.GetSpeciesId("B"), _rho_B);

        //        return R;
        //    }
        //}


        public SIMPLE(LevelSetTracker _LsTrk){
            this.LsTrk = _LsTrk;
        }

        void ExtractMatrices() {

            if(this.m_PressureSolver != null) {
                this.m_PressureSolver.Dispose();
                this.m_PressureSolver = null;
            }


            // sub-matrices for the SIMPLE
            // ===========================
            int VelocityLength = this.USubMatrixIdx_Row.Length;
            int PressureLength = this.PSubMatrixIdx_Row.Length;


            this.PressureGrad = new MsrMatrix(VelocityLength, PressureLength, 1, 1);
            this.VelocityDiv = new MsrMatrix(PressureLength, VelocityLength, 1, 1);
            this.ConvDiff = new MsrMatrix(VelocityLength, VelocityLength, 1, 1);
            this.Stab = new MsrMatrix(PressureLength, PressureLength, 1, 1);

            var WholeSystemMatrix = this.m_MgOp.OperatorMatrix;

            WholeSystemMatrix.WriteSubMatrixTo(PressureGrad, USubMatrixIdx_Row, default(int[]), PSubMatrixIdx_Row, default(int[]));
            WholeSystemMatrix.WriteSubMatrixTo(VelocityDiv, PSubMatrixIdx_Row, default(int[]), USubMatrixIdx_Row, default(int[]));
            WholeSystemMatrix.WriteSubMatrixTo(ConvDiff, USubMatrixIdx_Row, default(int[]), USubMatrixIdx_Row, default(int[]));
            WholeSystemMatrix.WriteSubMatrixTo(Stab, PSubMatrixIdx_Row, default(int[]), PSubMatrixIdx_Row, default(int[]));

            double condNo = ConvDiff.condest();
            Console.WriteLine("Convection/Diffusion Condition number: {0:0.####E-00}", condNo);


            //double[] ones = new double[D + 1];
            //ones.SetAll(1.0);
            //MassMatrix = MassFact.GetMassMatrix(Velocity.Current.Mapping, Rho, false)._ToMsrMatrix();
            //MassMatrixInv = MassFact.GetMassMatrix(Velocity.Current.Mapping, Rho, true)._ToMsrMatrix();
        }

        
        private void ApproximationMatrix() {

            this.Aapprox = null;
            this.AapproxInverse = null;

            //MsrMatrix AapproxComp = null;

            /*
            switch(m_SIMPLEOptions.Option_Approximation_Predictor) {
                case ApproxPredictor.MassMatrix:
                case ApproxPredictor.LocalizedOperator: { break; }
                default: {
                    Aapprox = ConvDiff.CloneAs();
                    //switch (m_SIMPLEOptions.Option_Timestepper) {
                    //    case Timestepper.Steady: break;
                    //    case Timestepper.ImplicitEuler: {
                    //            Aapprox.Acc(1.0 / dt, MassMatrix);
                    //            break;
                    //        }
                    //    default: {
                    //            throw new NotImplementedException("Unknown Timestepper");
                    //        }
                    //}
                    AapproxComp = new MsrMatrix(USubMatrixIdx.Length, USubMatrixIdx.Length, MassMatrix.RowPartitioning.BlockSize / (2 * D), MassMatrix.ColPartition.BlockSize / (2 * D));
                    Aapprox.WriteSubMatrixTo(AapproxComp, USubMatrixIdx, default(int[]), USubMatrixIdx, default(int[]));
                    break;
                }
            }
             */
            switch(m_SIMPLEOptions.Option_Approximation_Predictor) {
                case ApproxPredictor.MassMatrix: {

                    BlockMsrMatrix MM;
                    if(!double.IsPositiveInfinity(this.m_SIMPLEOptions.dt)) {
                        // instationary SIMPLE

                        //MM = this.m_MgOp.MassMatrix.CloneAs();

                        
                        // hier muss ich mir nochmal was überlegen --
                        // für einige Präkond.-Optionen 
                        // (genau jene, welche die XDG-Basen für beide Phasen in Cut-Zellen mischen),
                        // wie etwa 
                        //   MultigridOperator.Mode.SymPart_DiagBlockEquilib
                        // ist eine Block-Skalierung mit rho_A und rho_B
                        // inkonsistent!
                        // 

                        throw new NotImplementedException("todo");
                    } else {
                        MM = this.m_MgOp.MassMatrix;
                    }
                    Aapprox = new MsrMatrix(this.ConvDiff.RowPartitioning);
                    this.m_MgOp.MassMatrix.WriteSubMatrixTo(Aapprox, this.USubMatrixIdx_Row, default(int[]), this.USubMatrixIdx_Row, default(int[]));
                    
                    //AapproxInverse = MassMatrixInv._ToMsrMatrix();// Aapprox.Invert();
                    //switch(m_SIMPLEOptions.Option_Timestepper) {
                    //    case Timestepper.Steady: break;
                    //    case Timestepper.ImplicitEuler: {
                    //        Aapprox.Scale(1 + 1.0 / dt);
                    //        AapproxInverse.Scale(1 / (1 + 1.0 / dt));
                    //        /*#if DEBUG
                    //                                            var CheckMX = AapproxInverse * Aapprox;
                    //                                            foreach (int i in USubMatrixIdx) {
                    //                                                if (Math.Abs(CheckMX.GetDiagonalElement(i) - 1.0) > 1e-10) throw new ArithmeticException("AapproxInverse is not the Inverse of the Aapprox-Matrix");
                    //                                            }

                    //        #endif*/
                    //        break;
                    //    }
                    //    default: {
                    //        throw new NotImplementedException("Unknown Timestepper");
                    //    }
                    //}
                    break;
                }
                case ApproxPredictor.Exact: {
                    if (this.LsTrk.GridDat.CellPartitioning.MpiSize > 1)
                        throw new NotSupportedException("Not implemented for MPI-parallel runs.");

                    //RowIdx = VelocityMapping.GetSubvectorIndices(this.LsTrk, D.ForLoop(d => d), _SpcIds: this.LsTrk.SpeciesIdS, drk: this.TransportAgglomerator);
                    //ColIdx = RowIdx;

                    if(USubMatrixIdx_Row.Length > 4500)
                        Console.WriteLine(string.Format("WARNING: you don't really want to invert a {0}x{0} matrix.", USubMatrixIdx_Row.Length));

                    this.Aapprox = this.ConvDiff;

                    MultidimensionalArray AapproxFull = Aapprox.ToFullMatrixOnProc0();
                    MultidimensionalArray AapproxInverseFull = AapproxFull.GetInverse();
                    this.AapproxInverse = new MsrMatrix(new Partitioning(USubMatrixIdx_Row.Length));
                    this.AapproxInverse.AccDenseMatrix(1.0, AapproxInverseFull);
                    break;
                }
                case ApproxPredictor.Diagonal: {
                    /*
                    Aapprox = new MsrMatrix(ConvDiff.RowPartitioning, ConvDiff.ColPartition);
                    AapproxInverse = new MsrMatrix(ConvDiff.RowPartitioning, ConvDiff.ColPartition);
                    foreach(int i in USubMatrixIdx) {
                        int[] j = new int[] { i };
                        double Value = ConvDiff.GetValues(i, j)[0];
                        //Value += MassMatrix.GetValues(i, j)[0];
                        Aapprox.SetDiagonalElement(i, Value);
                        if(Value == 0) {
                            AapproxInverse.SetDiagonalElement(i, 0);
                        } else {
                            AapproxInverse.SetDiagonalElement(i, 1 / Value);
                        }
                    }
#if DEBUG
                    Aapprox.VerifyDataStructure();
                    AapproxInverse.VerifyDataStructure();
                    var CheckMX = AapproxInverse * Aapprox;
                    foreach(int i in USubMatrixIdx) {
                        if(Math.Abs(CheckMX.GetDiagonalElement(i) - 1.0) > 1e-13) throw new ArithmeticException("AapproxInverse is not the Inverse of the Operator Matrix");
                    }
#endif
                    break;
                    */
                    throw new NotImplementedException("todo");
                }
                case ApproxPredictor.BlockDiagonal: {
                    /*
                    AapproxInverse = new MsrMatrix(Aapprox.RowPartitioning, Aapprox.ColPartition);
                    var AapproxCompBD = new BlockDiagonalMatrix(AapproxComp);
                    var AapproxCompInvBD = AapproxCompBD.Invert();
                    AapproxCompBD._ToMsrMatrix().WriteSubMatrixTo(Aapprox, default(int[]), USubMatrixIdx, default(int[]), USubMatrixIdx);
                    AapproxCompInvBD._ToMsrMatrix().WriteSubMatrixTo(AapproxInverse, default(int[]), USubMatrixIdx, default(int[]), USubMatrixIdx);
                    //#if DEBUG
                    Aapprox.VerifyDataStructure();
                    AapproxInverse.VerifyDataStructure();
                    var CheckMX = AapproxInverse * Aapprox;
                    foreach(int i in USubMatrixIdx) {
                        if(Math.Abs(CheckMX.GetDiagonalElement(i) - 1.0) > 1e-12) throw new ArithmeticException("AapproxInverse is not the Inverse of the Operator Matrix");
                    }
                    //#endif
                    
                    break;
                     */
                    throw new NotImplementedException("todo");
                }
                case ApproxPredictor.BlockSum: {
                    /*
                    Console.WriteLine("BlockSum is not properly tested yet and did not work in previous tests");
                    int AccdBlockSize = ConvDiff.RowPartitioning.BlockSize / (2 * D);
                    var AapproxBD = new BlockDiagonalMatrix(ConvDiff.RowPartitioning);
                    var Aapprox = new MsrMatrix(ConvDiff.RowPartitioning);
                    int[] indexer = new int[AccdBlockSize];
                    for(int i = 0; i < indexer.Length; i++) {
                        indexer[i] = i;
                    }
                    int[] rowindexer = indexer.CloneAs();


                    for(int i = 0; i < ConvDiff.NoOfRows / AccdBlockSize; i++) {
                        int[] colindexer = indexer.CloneAs();
                        for(int j = 0; j < ConvDiff.NoOfCols / AccdBlockSize; j++) {
                            ConvDiff.AccSubMatrixTo(1.0, Aapprox, rowindexer, rowindexer, colindexer, rowindexer);
                            for(int r = 0; r < indexer.Length; r++) {
                                colindexer[r] += AccdBlockSize;
                            }
                        }
                        for(int r = 0; r < indexer.Length; r++) {
                            rowindexer[r] += AccdBlockSize;
                        }
                    }
                    //if (m_SIMPLEOptions.Option_Timestepper == Timestepper.ImplicitEuler) {
                    //    Aapprox.Acc(1 / dt, MassMatrix);
                    //}
                    AapproxComp = new MsrMatrix(USubMatrixIdx.Length, USubMatrixIdx.Length, AccdBlockSize, AccdBlockSize);
                    Aapprox.WriteSubMatrixTo(AapproxComp, USubMatrixIdx, default(int[]), USubMatrixIdx, default(int[]));
                    AapproxInverse = new MsrMatrix(Aapprox.RowPartitioning, Aapprox.ColPartition);
                    var AapproxCompBD = new BlockDiagonalMatrix(AapproxComp);
                    var AapproxCompInvBD = AapproxCompBD.Invert();
                    AapproxCompInvBD._ToMsrMatrix().WriteSubMatrixTo(AapproxInverse, default(int[]), USubMatrixIdx, default(int[]), USubMatrixIdx);
#if DEBUG
                    Aapprox.VerifyDataStructure();
                    AapproxInverse.VerifyDataStructure();

                    // Check copying back and forth
                    var CheckAapproxComp = new MsrMatrix(Aapprox.RowPartitioning, Aapprox.ColPartition);
                    AapproxComp.WriteSubMatrixTo(CheckAapproxComp, default(int[]), USubMatrixIdx, default(int[]), USubMatrixIdx);
                    CheckAapproxComp.Acc(-1.0, Aapprox);
                    if(CheckAapproxComp.InfNorm() > 1e-14) throw new ArithmeticException("Something went wrong while copying the Aapprox Matrix");

                    //Check Transformation to BlockdiagonalMatrix
                    var CheckAapproxCompBD = new MsrMatrix(Aapprox.RowPartitioning, Aapprox.ColPartition);
                    AapproxCompBD._ToMsrMatrix().WriteSubMatrixTo(CheckAapproxCompBD, default(int[]), USubMatrixIdx, default(int[]), USubMatrixIdx);
                    CheckAapproxCompBD.Acc(-1.0, Aapprox);
                    if(CheckAapproxCompBD.InfNorm() > 1e-14) throw new ArithmeticException("Something went wrong while copying the Aapprox Matrix");

                    //Check Matrix Inversion
                    var CheckMX = AapproxInverse * Aapprox;
                    foreach(int i in USubMatrixIdx) {
                        if(Math.Abs(CheckMX.GetDiagonalElement(i) - 1.0) > 1e-12) throw new ArithmeticException("AapproxInverse is not the Inverse of the Operator Matrix");
                    }
#endif
                    break;
                     */
                    throw new NotImplementedException("todo");
                }
                case ApproxPredictor.Neumann: {
                    /*
                    Console.WriteLine("Neumann did not work in previous Tests, Series does typically not converge");
                    int serieslength = 10;

                    Aapprox = ConvDiff.CloneAs();
                    //if (m_SIMPLEOptions.Option_Timestepper != Timestepper.Steady) {
                    //    Aapprox.Acc(1.0 / dt, MassMatrix);
                    //}
                    var B = Aapprox.CloneAs();
                    double ScalingFactor = Aapprox.InfNorm();
                    B.Scale((-1.0 / ScalingFactor.Pow(0))); //Scaling Power up to 4 tried
                    B.AccEyeSp(1.0);
                    AapproxInverse = B;
                    AapproxInverse.AccEyeSp(1.0);
                    MsrMatrix OldMoment = B;
                    for(int i = 0; i <= serieslength; i++) {
                        var NewMoment = OldMoment * B;
                        AapproxInverse.Acc(1.0, NewMoment);
                        //Debug
                        Console.WriteLine("MomentNumber #{0}, InfNormOf Inverse #{1}", i, NewMoment.InfNorm());
                        OldMoment = NewMoment;
                    }
                    AapproxInverse.Scale((1.0 / ScalingFactor.Pow(0))); //Scaling Power up to 4 tried
                    break;
                     */
                    throw new NotImplementedException("todo");
                }
                case ApproxPredictor.LocalizedOperator: {
                    /*
                    Console.WriteLine("Localized Operator did not work in previous Tests");
                    double[] LocalizedOpAffine;
                    MultiphaseCellAgglomerator LocalizedAgglomerator;
                    Aapprox = new MsrMatrix(MassMatrix.RowPartitioning, MassMatrix.ColPartition);
                    TransportOpLocalized.AssembleMatrix(
                        out Aapprox, out LocalizedOpAffine,
                        out LocalizedAgglomerator, out TransportMassFact,
                        this.Velocity.Current, null,
                        this.LevSet, null, Curv,
                        VelocityMapping, VelocityMapping);
                    if(Option_Timestepper == Timestepper.ImplicitEuler) {
                        Aapprox.Acc(1 / dt, MassMatrix);
                    }

                    var DiagAverage = Aapprox.GetDiagVector().Average();
                    foreach(int i in RowIdx) {
                        if(Aapprox.GetDiagonalElement(i) == 0.0) {
                            //Aapprox.SetDiagonalElement(i, MassMatrix.GetDiagonalElement(i));
                            Aapprox.SetDiagonalElement(i, DiagAverage);
                        }
                    }
                    AapproxComp = new MsrMatrix(RowIdx.Length, RowIdx.Length, MassMatrix.RowPartitioning.BlockSize / (2 * D), MassMatrix.ColPartition.BlockSize / (2 * D));
                    Aapprox.WriteSubMatrixTo(AapproxComp, RowIdx, default(int[]), ColIdx, default(int[]));
                    AapproxInverse = new MsrMatrix(Aapprox.RowPartitioning, Aapprox.ColPartition);
                    var AapproxCompBD = new BlockDiagonalMatrix(AapproxComp);
                    var AapproxCompInvBD = AapproxCompBD.Invert();
                    AapproxCompInvBD._ToMsrMatrix().WriteSubMatrixTo(AapproxInverse, default(int[]), RowIdx, default(int[]), ColIdx);
                    //AapproxComp._ToMsrMatrix().WriteSubMatrixTo(Aapprox, default(int[]), RowIdx, default(int[]), ColIdx);


#if DEBUG
                    Aapprox.VerifyDataStructure();
                    AapproxInverse.VerifyDataStructure();
                    var CheckAapproxCompBD = new MsrMatrix(Aapprox.RowPartitioning, Aapprox.ColPartition);
                    AapproxCompBD._ToMsrMatrix().WriteSubMatrixTo(CheckAapproxCompBD, default(int[]), RowIdx, default(int[]), ColIdx);
                    CheckAapproxCompBD.Acc(-1.0, Aapprox);
                    if(CheckAapproxCompBD.InfNorm() > 1e-14) throw new ArithmeticException("Something went wrong while copying the Aapprox Matrix");
                    var CheckMX = AapproxInverse * Aapprox;
                    foreach(int i in RowIdx) {
                        if(Math.Abs(CheckMX.GetDiagonalElement(i) - 1.0) > 1e-12) throw new ArithmeticException("AapproxInverse is not the Inverse of the Operator Matrix");
                    }
#endif
                    break;
                    */
                    throw new NotImplementedException("todo");
                }
                    
                default:
                throw new NotImplementedException("todo");
            }


            if(this.AapproxInverse == null) {
                // block-inversion is required.
                // ++++++++++++++++++++++++++++

                int D = this.LsTrk.GridDat.SpatialDimension;
                this.AapproxInverse = new MsrMatrix(this.Aapprox.RowPartitioning, this.Aapprox.ColPartition);

                //int N = this.m_MgOp.Mapping.AggBasis.GetMinimalLength(this.m_MgOp.Mapping.DgDegree[0]);
                //Debug.Assert(this.Aapprox.RowPartitioning.LocalLength % N == 0);
                MultidimensionalArray Block = new MultidimensionalArray(2);
                MultidimensionalArray InvBlock = new MultidimensionalArray(2);


                int iRow0 =  this.Aapprox.RowPartitioning.i0;
                int JAGG = this.m_MgOp.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                int[] DegreeS = this.m_MgOp.Mapping.DgDegree;
                for(int jagg = 0; jagg < JAGG; jagg++) { // loop over aggregate cells...
                    for(int d = 0; d < D; d++) { // loop over velocity components...
                        int N = this.m_MgOp.Mapping.AggBasis[d].GetLength(jagg, DegreeS[d]);
                        if(Block.GetLength(0) != N) {
                            Block.Allocate(N, N);
                            InvBlock.Allocate(N, N);
                        }

                        for(int n = 0; n < N; n++) {
#if DEBUG
                            int iRow = iRow0 + n;
                            {
                                int[] Cols = null;
                                double[] Vals = null;
                                int LR = Aapprox.GetRow(iRow, ref Cols, ref Vals);
                                int cMin = int.MaxValue;
                                int cMax = int.MinValue;
                                for (int lr = 0; lr < LR; lr++) {
                                    if (Vals[lr] != 0.0) {
                                        cMin = Math.Min(cMin, Cols[lr]);
                                        cMax = Math.Max(cMax, Cols[lr]);
                                    }
                                }
                                Debug.Assert(cMin >= iRow0);
                                Debug.Assert(cMax < iRow0 + N);
                            }
                            
#endif
                            for (int m = 0; m < N; m++) {
                                Block[n, m] = this.Aapprox[iRow0 + n, iRow0 + m];
                            }
                        }

                        Block.InvertTo(InvBlock);

                        for(int n = 0; n < N; n++) {
                            for(int m = 0; m < N; m++) {
                                this.AapproxInverse[iRow0 + n, iRow0 + m] = InvBlock[n, m];
                            }
                        }

                        iRow0 += N;
                    }
                }

                Debug.Assert(iRow0 == this.Aapprox.RowPartitioning.iE);
            }
            
#if DEBUG
            var CheckMX = AapproxInverse * Aapprox;
            double TRESH = Math.Max(AapproxInverse.InfNorm(), Aapprox.InfNorm()) * 1.0e-10;
            for(int iRow = CheckMX.RowPartitioning.i0; iRow < CheckMX.RowPartitioning.iE; iRow++) {
                if(Math.Abs(CheckMX.GetDiagonalElement(iRow) - 1.0) > TRESH)
                    throw new ArithmeticException("AapproxInverse is not the Inverse of the Aapprox-Matrix");
            }

#endif
        }


        void PressureSolver(double[] Pressure, double[] Velocity, double[] RHSMomentum) {
            using(new FuncTrace()) {
                // solve Poisson equation
                // ======================

                // Matrix: Div*A^-1*Grad + S
                InitPressureSolver();
                
                // RHS
                double[] RHS = new double[Pressure.Length];
                {
                    double[] tmp1 = RHSMomentum.CloneAs();
                    this.ConvDiff.SpMVpara(-1.0, Velocity, 1.0, tmp1);
                    double[] tmp2 = new double[RHSMomentum.Length];

                    this.AapproxInverse.SpMVpara(1.0, tmp1, 0.0, tmp2);

                    this.VelocityDiv.SpMVpara(1.0, tmp2, 0.0, RHS);
                }

                // solve
                m_PressureSolver.Solve(Pressure, RHS);

            }
        }


        /// <summary>
        /// computes the pressure correction.
        /// </summary>
        /// <param name="Pressure">input: current pressure</param>
        /// <param name="VelocityPCorrIn">input: intermediate velocity</param>
        /// <param name="PressurePCorr">output: pressure correction</param>
        /// <param name="RHSContinuity">input: RHS of the continuity equation</param>
        /// <returns></returns>
        double PressureCorrector(double[] Pressure, double[] VelocityPCorrIn, double[] PressurePCorr, double[] RHSContinuity) {
            using(new FuncTrace()) {

                // solve Poisson equation
                // ======================

                // Matrix: Div*A^-1*Grad + S
                InitPressureSolver();

                // RHS
                var RHS = new double[Pressure.Length];

                {
                    RHS.AccV(-1.0, RHSContinuity);
                    VelocityDiv.SpMVpara(+1.0, VelocityPCorrIn, 1.0, RHS);
                    Stab.SpMVpara(+1.0, Pressure, 1.0, RHS);
                }

                // solve
                m_PressureSolver.Solve(PressurePCorr, RHS);
                
                
                // return
                return PressurePCorr.L2Norm();

            }
        }

        private void InitPressureSolver() {
            if(m_PressureSolver == null) {
                var MXCorrector = VelocityDiv * (AapproxInverse * PressureGrad);
                MXCorrector.Acc(-1.0, Stab);
                m_PressureSolver = this.m_SIMPLEOptions.PressureSolver;
                if(m_PressureSolver is PARDISOSolver) {
                    ((PARDISOSolver)m_PressureSolver).CacheFactorization = true;
                }
                m_PressureSolver.DefineMatrix(MXCorrector);
            }
        }

        ISparseSolver m_PressureSolver;


        public void VelocityPredictor(double[] PressureEstimateIn, double[] VelocityEstimateIn, double[] VelocityPrediction, double[] RHSMomentum) {
            using (new FuncTrace()) {
                double m_relax_vel = m_SIMPLEOptions.relax_v;

                //build Matrix
                var PredictorMX = ConvDiff.CloneAs();

                
                //underrelaxation LHS
                if(m_relax_vel <= 0.0)
                    throw new ArithmeticException("Illegal velocity underrelaxation parameter: " + m_relax_vel);
                
                PredictorMX.Acc((1 - m_relax_vel) / m_relax_vel, Aapprox);


#if DEBUG
                PredictorMX.CheckForNanOrInfM();
#endif


                // build RHS: b1-grad*p (+ oldVelocities/dt)
                double[] RHS = new double [RHSMomentum.Length];
                RHS.AccV(1.0, RHSMomentum);
                PressureGrad.SpMVpara(-1.0, PressureEstimateIn, 1.0, RHS);

                //underrelaxation RHS
                Aapprox.SpMVpara((1 - m_relax_vel) / m_relax_vel, VelocityEstimateIn, 1.0, RHS);

                // solve
                using(ISparseSolver solver = m_SIMPLEOptions.ViscousSolver) {
                    //Agglomerator.ClearAgglomerated(RHS, VelocityMapping);
                    //double SolverResidual = PC.SolveDirect(VelocityPrediction.CoordinateVector, RHS, solver, false);
                    //solver.Dispose();
                    solver.DefineMatrix(PredictorMX);
                    solver.Solve(VelocityPrediction, RHS);
                }
            }
        }

        void VelocityUpdate(double[] PressureCorrection, double[] VelocityPrediction, double[] FinalVelocity) {
            FinalVelocity.SetV(VelocityPrediction);

            double[] temp = new double[FinalVelocity.Length];
            PressureGrad.SpMVpara(-1.0, PressureCorrection, 0.0, temp);
            AapproxInverse.SpMVpara(1.0, temp, 1.0, FinalVelocity);
        }

        //void PressureUpdate(double[] PressureCorrection, double[] FinalPressure) {
        //    if(m_SIMPLEOptions.relax_p < 0.0)
        //        throw new ArithmeticException("Illegal pressure correction relaxation parameter: " + m_SIMPLEOptions.relax_p);
        //    FinalPressure.AccV(m_SIMPLEOptions.relax_p, PressureCorrection);
        //}

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }
    }
}
