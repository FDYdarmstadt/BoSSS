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
using BoSSS.Foundation;
using ilPSP.Utils;
using BoSSS.Solution.Multigrid;
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP;
using MPI.Wrappers;
using BoSSS.Foundation.XDG;

namespace BoSSS.Solution.Multigrid {


   
    public class NonlinGMRES : NonlinearSolver {

        public NonlinGMRES(OperatorEvalOrLin __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq, MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig) 
            : base(__AssembleMatrix, __AggBasisSeq, __MultigridOperatorConfig) //
        {
        }

        public int MaxIter = 400;
        public double ConvCrit = -1.0;


        public int MaxKrylovDim = 80;


#pragma warning disable 414
        // from the Oosterlee & Whasio paper.
        double gamma_A = 2.0;
        double epsilon_B = 0.1;
        double delta_B = 0.9;
#pragma warning restore 414

        override public void SolverDriver<S>(CoordinateVector VelocityAndPressure, S RHS) {
           
            // history of solutions and residuals (max vector length 'MaxKrylovDim')
            List<double[]> SolHistory = new List<double[]>();
            List<double[]> ResHistory = new List<double[]>();

            // initial guess and its residual:
            double[] Sol1, Res1;
            base.Init(VelocityAndPressure, RHS, out Sol1, out Res1);
            SolHistory.Add(Sol1.CloneAs());
            Sol1.ClearEntries();
            ResHistory.Add(Res1.CloneAs());
            Res1.ClearEntries();

            // norm of initial residual
            MultidimensionalArray InnerProds = MultidimensionalArray.Create(200, 200); // memory for saving the inner products of the residuals
            InnerProds[0, 0] = GenericBlas.L2NormPow2(ResHistory[0]).MPISum();

            // diagnostic output
            OnIterationCallback(0, SolHistory.Last().CloneAs(), ResHistory.Last().CloneAs(), this.CurrentLin);
            throw new NotImplementedException("todo: missing trafo of solution history under linearization update");
    /*
            // let's iterate ...
            for(int iIter = 0; iIter < MaxIter; iIter++) {
                Debug.Assert(SolHistory.Count == ResHistory.Count);
                Debug.Assert(SolHistory.Count >= 1);

                // (approximately) solve the linearized equation:
                Precond.Init(this.CurrentLin);
                Sol1.SetV(SolHistory.Last(), 1.0);
                Precond.Solve(Sol1, this.LinearizationRHS);
                SolHistory.Add(Sol1.CloneAs());
                Sol1.ClearEntries();

                // update linearization of nonlinear equation
                this.CurrentLin.TransformSolFrom(VelocityAndPressure, SolHistory.Last());
                this.Update(VelocityAndPressure.Mapping.Fields);
                throw new NotImplementedException("todo: missing trafo of solution history under linearization update");
                
                
                // compute the new residual
                this.EvalResidual(SolHistory.Last(), ref Res1);
                ResHistory.Add(Res1.CloneAs());
                Res1.ClearEntries();

                // compute the new inner products
                Debug.Assert(SolHistory.Count == ResHistory.Count);
                int m = ResHistory.Count - 1;
                for(int iInnerIter = 0; iInnerIter <= m; iInnerIter++) {
                    if(m != iInnerIter) {
                        InnerProds[iInnerIter, m] = GenericBlas.InnerProd(ResHistory[iInnerIter], ResHistory[m]).MPISum();
                        InnerProds[m, iInnerIter] = InnerProds[iInnerIter, m];
                    } else {
                        InnerProds[m, m] = GenericBlas.L2NormPow2(ResHistory[m]).MPISum();
                    }
                }
                

                // compute 'accelerated solution'
                {
                    // try to find an accelerated solution
                    // +++++++++++++++++++++++++++++++++++

                    double[] alpha = Minimi(InnerProds.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { m, m }));
                    //alpha.ClearEntries();
                    Debug.Assert(alpha.Length == m);

                    double[] SolM = SolHistory.Last();
                    double[] ResM = ResHistory.Last();
                    double Norm_ResM = GenericBlas.L2NormPow2(ResM).MPISum().Sqrt();
                    

                    double[] SolA = Sol1; // re-use the mem.
                    double[] ResA = Res1;
                    Debug.Assert(object.ReferenceEquals(SolA, SolM) == false);
                    Debug.Assert(object.ReferenceEquals(ResA, ResM) == false);
                    Debug.Assert(SolA.L2Norm() == 0.0);
                    SolA.AccV(1.0, SolM);

                    
                    for(int i = 0; i < m; i++) {
                        SolA.AccV(alpha[i], SolHistory[i]);
                        SolA.AccV(-alpha[i], SolM);
                    }

                    //this.Update(SolA);
                    this.EvalResidual(SolA, ref ResA);

                    SolHistory.Last().SetV(SolA, 1.0);
                    ResHistory.Last().SetV(ResA, 1.0);

                    for(int iInnerIter = 0; iInnerIter <= m; iInnerIter++) {
                        if(m != iInnerIter) {
                            InnerProds[iInnerIter, m] = GenericBlas.InnerProd(ResHistory[iInnerIter], ResHistory[m]).MPISum();
                            InnerProds[m, iInnerIter] = InnerProds[iInnerIter, m];
                        } else {
                            InnerProds[m, m] = GenericBlas.L2NormPow2(ResHistory[m]).MPISum();
                        }
                    }



                    //double Norm_ResA = GenericBlas.L2NormPow2(ResA).MPISum().Sqrt();
                    //Console.WriteLine("                 {0}", Norm_ResM/Norm_ResA);

                    double minNorm_Res = double.MaxValue;
                    for(int i = 0; i <= m; i++) {
                        minNorm_Res = Math.Min(minNorm_Res, InnerProds[i, i].Sqrt());
                    }


                    //if(Norm_ResA < this.gamma_A * minNorm_Res  // first condition for accepting accelerated solution
                    //    &&
                    //    (Norm_ResA < this.delta_B * minNorm_Res
                    //        || MinDistCond(SolA, SolM, SolHistory.Take(m).ToArray(), this.epsilon_B))) //
                    if(SolHistory.Count > this.MaxKrylovDim) {
                    //if(true) {
                        // also the second criterion is fulfilled

                        // accept accelerated solution and Restart;
                        SolHistory.Clear();
                        SolHistory.Add(SolA.CloneAs());
                        ResHistory.Clear();
                        ResHistory.Add(ResA.CloneAs());

                        Console.WriteLine("restart.", m);

                    } 
                    //else {
                    //    // we do not accept the accelerated solution, but we go on...

                    //    this.Update(SolA);

                    //    //Console.WriteLine("Inner GMRES iteration {0}: Accelertated solution NOT accepted.", m);

                    //}


                }

                

                OnIterationCallback(iIter + 1, SolHistory.Last().CloneAs(), ResHistory.Last().CloneAs(), this.CurrentLin);
                
            }

            this.CurrentLin.TransformSolFrom(VelocityAndPressure, SolHistory.Last());
            */
        }

        static bool MinDistCond(double[] SolA, double[] SolM, double[][] PrevSol, double epsilon_B) {
            Debug.Assert(object.ReferenceEquals(SolA, SolM) == false);

            double dist_solAsolM = GenericBlas.L2DistPow2(SolA, SolM).MPISum().Sqrt();


            int m = PrevSol.Length;

            double minDist_SolA_Soli = double.MaxValue;
            for(int i = 0; i < m; i++) {
                Debug.Assert(object.ReferenceEquals(PrevSol[i], SolA) == false);
                Debug.Assert(object.ReferenceEquals(PrevSol[i], SolM) == false);

                double Dist_SolA_Soli = GenericBlas.L2DistPow2(SolA, PrevSol[i]);
                minDist_SolA_Soli = Math.Min(minDist_SolA_Soli, Dist_SolA_Soli);
            }

            return epsilon_B * dist_solAsolM < minDist_SolA_Soli;
        }


        /// <summary>
        /// residual minimization
        /// </summary>
        static double[] Minimi(MultidimensionalArray InnerProds) {
            int m = InnerProds.GetLength(0) - 1;
            Debug.Assert(InnerProds.Dimension == 2);
            Debug.Assert(InnerProds.GetLength(0) == InnerProds.GetLength(1));

            MultidimensionalArray EqSys = MultidimensionalArray.Create(m, m);
            double[] RHS = new double[m];
            double[] alpha = new double[m];

//#if DEBUG
//            // orginal-formulierung im paper
//            MultidimensionalArray OrgEqSys = MultidimensionalArray.Create(m, m);
            
//#endif

            for(int i = 0; i < m; i++) {
                for(int k = 0; k < m; k++) {
                    EqSys[i, k] = InnerProds[i, k] - InnerProds[m, i] - InnerProds[m, k] + InnerProds[m, m];
                    //OrgEqSys[i, k] = InnerProds[m - i, m - k] - InnerProds[m, m - i] - InnerProds[m, m - k] + InnerProds[m, m];
                }

                RHS[i] = InnerProds[m, m] - InnerProds[m, i];
                //RHS[i] = InnerProds[m, m] - InnerProds[m, m - i];
            }

            EqSys.Solve(alpha, RHS);
            return alpha;
        }
    }

}
