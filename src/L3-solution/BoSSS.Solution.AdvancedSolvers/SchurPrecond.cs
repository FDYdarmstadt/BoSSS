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
using MPI.Wrappers;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using BoSSS.Foundation;
using ilPSP.Connectors.Matlab;
using BoSSS.Solution.NSECommon;
using System.Diagnostics;
using ilPSP.Tracing;
using ilPSP.LinSolvers.PARDISO;
using BoSSS.Solution.Control;

namespace BoSSS.Solution.AdvancedSolvers {

	[Serializable]
	public class SchurPrecondConfig : IterativeSolverConfig {
		public override string Name => "Schur complement with Uzawa algorithm";

		public override string Shortname => "Uzawa";

		public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {

			var templinearSolve = new SchurPrecond(this);


			templinearSolve.Init(level);
			return templinearSolve;
		}
	}


    public class SchurPrecond : ISolverSmootherTemplate, ISolverWithCallback {
        
        /// <summary>
        /// ctor
        /// </summary>
        public SchurPrecond(SchurPrecondConfig config) {
            m_config = config;
        }

		readonly SchurPrecondConfig m_config;

		public int IterationsInNested {
            get {
                return m_IterationsInNested;
                //throw new NotImplementedException();
            }
        }

        public int ThisLevelIterations {
            get {
                return m_ThisLevelIterations;
                //throw new NotImplementedException();
            }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get {
                return null;
                //throw new NotImplementedException();
            }
            set {

                //throw new NotImplementedException();
            }
        }

        MultigridOperator m_mgop;
        MultigridMapping MgMap {
            get {
                return m_mgop.Mapping;	}
            }
            

		BlockMsrMatrix Mtx;
		//BlockMsrMatrix operatorM;

        MsrMatrix P;
        MsrMatrix ConvDiff, pGrad, divVel, SchurMtx, SchurRHSMtx, PoissonMtx_T, PoissonMtx_H, SchurConvMtx, invVelMassMatrix, invVelMassMatrixSqrt, simpleSchur, velMassMatrix, pMassMatrix;
        long[] Uidx, Pidx;
		int[] UidxInt, PidxInt;

        int D;

        int m {
            get { return Uidx.Length; }
        }

		int n {
			get { return Pidx.Length; }
		}

		public enum SchurOptions { Uzawa = 0, exact = 1, decoupledApprox = 2, SIMPLE = 3, exact_matlab = 4, }

        public bool ApproxScaling = false;
        
        public SchurOptions SchurOpt = SchurOptions.Uzawa;

        public void Init(MultigridOperator op)
        {
			this.m_mgop = op;
			D = op.Mapping.GridData.SpatialDimension;
            var M = op.OperatorMatrix;

           
			if (!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!M.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");

            Uidx = MgMap.GetSubvectorIndices( D.ForLoop(i => i));
            Pidx = MgMap.GetSubvectorIndices( D);
            UidxInt = Array.ConvertAll(Uidx, item => (int)item);
			PidxInt = Array.ConvertAll(Pidx, item => (int)item);

			int Upart = Uidx.Length;
            int Ppart = Pidx.Length;
            ConvDiff = new MsrMatrix(Upart, Upart, 1, 1);
            pGrad = new MsrMatrix(Upart, Ppart, 1, 1);
            divVel = new MsrMatrix(Ppart, Upart, 1, 1);
            var PxP = new MsrMatrix(Ppart, Ppart, 1, 1);


            M.AccSubMatrixTo(1.0, ConvDiff, Uidx, default(long[]), Uidx, default(long[]));//, default(int[]), default(int[]));
            M.AccSubMatrixTo(1.0, pGrad, Uidx, default(long[]), Pidx, default(long[]));//, default(int[]), default(int[]));
            M.AccSubMatrixTo(1.0, divVel, Pidx, default(long[]), Uidx, default(long[]));//, default(int[]), default(int[]));
            M.AccSubMatrixTo(1.0, PxP, Pidx, default(long[]), Pidx, default(long[]));//, default(int[]), default(int[]));

			Mtx = M;

            int L = M.RowPartitioning.LocalLength;

            long i0 = Mtx.RowPartitioning.i0;

            P = new MsrMatrix(Mtx);
            P.Clear();

            velMassMatrix = new MsrMatrix(Upart, Upart, 1, 1);
            op.MassMatrix.AccSubMatrixTo(1.0, velMassMatrix, Uidx, default(long[]), Uidx, default(long[]), default(long[]), default(long[]));
			pMassMatrix = new MsrMatrix(Ppart, Ppart, 1, 1);
			op.MassMatrix.AccSubMatrixTo(1.0, pMassMatrix, Pidx, default(long[]), Pidx, default(long[]), default(long[]), default(long[]));

			//ConvDiff.SaveToTextFileSparse("ConvDiff");
			//pGrad.SaveToTextFileSparse("pGrad");
			//divVel.SaveToTextFileSparse("divVel");
			//PxP.SaveToTextFileSparse("PxP");
			//velMassMatrix.SaveToTextFileSparse("velMassMatrix");
			//pMassMatrix.SaveToTextFileSparse("pMassMatrix");

			switch (SchurOpt)
            {
                case SchurOptions.exact_matlab:
                    {
                        // Building complete Schur and Approximate Schur
                        MultidimensionalArray Poisson = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);
                        MultidimensionalArray SchurConvPart = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);
                        MultidimensionalArray Schur = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);
                        using (BatchmodeConnector bmc = new BatchmodeConnector())
                        {
                            bmc.PutSparseMatrix(ConvDiff, "ConvDiff");
                            bmc.PutSparseMatrix(velMassMatrix, "MassMatrix");
                            bmc.PutSparseMatrix(divVel, "divVel");
                            bmc.PutSparseMatrix(pGrad, "pGrad");
                            bmc.Cmd("Qdiag = diag(diag(MassMatrix))");
                            bmc.Cmd("invT= inv(Qdiag)");
                            bmc.Cmd("Poisson = full(invT)*pGrad");
                            bmc.Cmd("ConvPart = ConvDiff*Poisson");
                            bmc.Cmd("ConvPart= full(invT)*ConvPart");
                            bmc.Cmd("ConvPart= divVel*ConvPart");
                            bmc.Cmd("Poisson = divVel*Poisson");
                            bmc.Cmd("ConvDiffInv = inv(full(ConvDiff))");
                            bmc.Cmd("Schur = divVel*ConvDiffInv");
                            bmc.Cmd("Schur = Schur*pGrad");
                            bmc.GetMatrix(Poisson, "Poisson");
                            bmc.GetMatrix(SchurConvPart, "ConvPart");
                            bmc.GetMatrix(Schur, "-Schur");
                            bmc.Execute(false);
                        }
                        PoissonMtx_T = Poisson.ToMsrMatrix();
                        PoissonMtx_H = Poisson.ToMsrMatrix();
                        SchurConvMtx = SchurConvPart.ToMsrMatrix();
                        SchurMtx = Schur.ToMsrMatrix();
                        SchurMtx.Acc(PxP, 1);

                        ConvDiff.AccSubMatrixTo(1.0, P, default(long[]), Uidx, default(long[]), Uidx);
                        pGrad.AccSubMatrixTo(1.0, P, default(long[]), Uidx, default(long[]), Pidx);
                        SchurMtx.AccSubMatrixTo(1.0, P, default(long[]), Pidx, default(long[]), Pidx);
                        return;
                    }
                case SchurOptions.decoupledApprox:
                    {
                        // Do assembly for approximate Schur inverse
                        invVelMassMatrix = velMassMatrix.CloneAs();
                        invVelMassMatrix.Clear();
                        invVelMassMatrixSqrt = invVelMassMatrix.CloneAs();
                        for (long i = velMassMatrix.RowPartitioning.i0; i < velMassMatrix.RowPartitioning.iE; i++)
                        {
                            if (ApproxScaling)
                            {
                                invVelMassMatrix.SetDiagonalElement(i, 1 / (velMassMatrix[i, i]));
                                invVelMassMatrixSqrt.SetDiagonalElement(i, 1 / (Math.Sqrt(velMassMatrix[i, i])));
                            }
                            else
                            {
                                invVelMassMatrix.SetDiagonalElement(i, 1);
                                invVelMassMatrixSqrt.SetDiagonalElement(i, 1);
                            }
                        }

                        //invVelMassMatrix.SaveToTextFileSparse("invVelMassMatrix");
                        //velMassMatrix.SaveToTextFileSparse("velMassMatrix");


                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(ConvDiff, pGrad);
                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(divVel, ConvDiffPoissonMtx);

                        // Inverse of mass matrix in Matlab
                        //MultidimensionalArray temp = MultidimensionalArray.Create(Uidx.Length, Uidx.Length);
                        //using (BatchmodeConnector bmc = new BatchmodeConnector())
                        //{
                        //    bmc.PutSparseMatrix(velMassMatrix, "velMassMatrix");
                        //    bmc.Cmd("invVelMassMatrix = inv(full(velMassMatrix))");
                        //    bmc.GetMatrix(temp, "invVelMassMatrix");
                        //    bmc.Execute(false);
                        //}
                        //invVelMassMatrix = temp.ToMsrMatrix();

                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(ConvDiffPoissonMtx, PoissonMtx);
                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(PoissonMtx, ConvDiffPoissonMtx);

                        //ConvDiff.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Uidx);
                        //pGrad.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Pidx);
                        //ConvDiffPoissonMtx.AccSubMatrixTo(1.0, P, default(int[]), Pidx, default(int[]), Pidx);

                        //op.MassMatrix.SaveToTextFileSparse("MassMatrix");
                        //velMassMatrix.SaveToTextFileSparse("velMassMatrix2");


                        // Possion scaled by inverse of the velocity mass matrix 
                        PoissonMtx_T = MsrMatrix.Multiply(invVelMassMatrix, pGrad); 
                        PoissonMtx_T = MsrMatrix.Multiply(divVel, PoissonMtx_T);
                        //PoissonMtx_T.Acc(PxP, 1); // p.379

                        // Poisson scaled by sqrt of inverse of velocity mass matrix
                        PoissonMtx_H = MsrMatrix.Multiply(invVelMassMatrixSqrt, pGrad);
                        PoissonMtx_H = MsrMatrix.Multiply(divVel, PoissonMtx_H);
                        //PoissonMtx_H.Acc(PxP, 1); // p.379
                        return;
                    }
                case SchurOptions.SIMPLE:
                    {

                        var invdiag_ConvDiff = ConvDiff.CloneAs();
                        invdiag_ConvDiff.Clear();
                        for (long i = ConvDiff.RowPartitioning.i0; i < ConvDiff.RowPartitioning.iE; i++)
                        {
                            invdiag_ConvDiff[i, i] = 1 / ConvDiff[i, i];
                        }

                        simpleSchur = MsrMatrix.Multiply(invdiag_ConvDiff, pGrad);
                        simpleSchur = MsrMatrix.Multiply(divVel, simpleSchur);

                        return;
                    }
				case SchurOptions.Uzawa: {
						Console.WriteLine("Uzawa is set");
						return;
					}
				case SchurOptions.exact: {
                        // Building complete Schur and Approximate Schur

                        // USING MATLAB
                        //using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                        //	bmc.PutSparseMatrix(ConvDiff, "A");
                        //	bmc.PutSparseMatrix(pGrad, "B");
                        //	bmc.PutSparseMatrix(divVel, "C");
                        //	bmc.Cmd("invA = inv(full(A));");
                        //	bmc.Cmd("Schur = -C *full(A \\ B);");
                        //	bmc.Cmd("SchurRHS = C*invA;"); 
                        //	bmc.GetMatrix(Schur, "Schur");
                        //	bmc.GetMatrix(SchurRHS, "SchurRHS");

                        //	bmc.Execute(false);
                        //}
                        //SchurMtx = Schur.ToMsrMatrix();
                        //SchurMtx.Acc(PxP, 1);

                        //SchurMtx.SaveToTextFileSparse("returnedSchur");

                        //SchurRHSMtx = SchurRHS.ToMsrMatrix();
                        //SchurRHSMtx.SaveToTextFileSparse("SchurRHSMtx");

                        //Using direct solver
                        var Ainv = new MsrMatrix(ConvDiff.RowPartitioning, ConvDiff.ColPartition);
                        int rank;
						csMPI.Raw.Comm_Rank(Ainv.MPI_Comm, out rank);

						ilPSP.Environment.StdoutOnlyOnRank0 = false;
						Console.WriteLine($"proc-{rank} - {Ainv.RowPartitioning.i0} {Ainv.RowPartitioning.iE}");
                        double inc = 0.0;
                        var solvPar = new PARDISOSolver();
						solvPar.DefineMatrix(ConvDiff);
						using (var tr = new FuncTrace()) { 
                            for (int i = (int)Ainv.RowPartitioning.i0; i < (int)Ainv.RowPartitioning.iE; i++) {
                                var b = new double[m];
								var x = new double[m];
								b[i - (int)Ainv.RowPartitioning.i0] = 1;
								solvPar.Solve(x, b);
                                Ainv.SetValues(i, Enumerable.Range(0, m).Select(r => (long)r).ToArray(), x);
                                if (i >= i0 + inc * Ainv.RowPartitioning.LocalLength) {
                                    Console.WriteLine($"{inc *100}%");
                                    inc += 0.1;
                                }
							}
						}
						csMPI.Raw.Barrier(Ainv.MPI_Comm);
						//Ainv = Ainv.Transpose(); (for stokes we do not need this at this point)

						var SchurSelfMSR = MsrMatrix.Multiply(Ainv, pGrad);
						SchurSelfMSR = MsrMatrix.Multiply(divVel, SchurSelfMSR);
						SchurSelfMSR.Scale(-1.0);

						//SchurSelfMSR.Acc(PxP, 1); //this is already zero
						//SchurSelfMSR.SaveToTextFileSparse("SchurSelfMSR");

						var SchurRHSselfMSR = MsrMatrix.Multiply(divVel, Ainv);
						//SchurRHSselfMSR.SaveToTextFileSparse("SchurRHSselfMSR");


						// Using multidimensional array
						//MultidimensionalArray Schur = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);
						//MultidimensionalArray SchurRHS = MultidimensionalArray.Create(Pidx.Length, Uidx.Length);
						//var A = ConvDiff.ToFullMatrixOnProc0();
						//                  var B = pGrad.ToFullMatrixOnProc0();
						//                  var C = divVel.ToFullMatrixOnProc0();

						//                  A.InvertInPlace();

						//var SchurSelf = A.MatMatMul(B);
						//var SchurSelf2 = C.MatMatMul(SchurSelf);
						//SchurSelf2.Scale(-1.0);
						//var SchurSelfMtx = SchurSelf2.ToMsrMatrix();
						//SchurSelfMtx.Acc(PxP, 1);
						//SchurSelfMtx.SaveToTextFileSparse("SchurSelfMtx");


						//var SchurRHSself = C.MatMatMul(A);
						//                  var SchurRHSselfMtx = SchurRHSself.ToMsrMatrix();
						//SchurRHSselfMtx.SaveToTextFileSparse("SchurRHSselfMtx");

						SchurMtx = SchurSelfMSR;
                        //SchurMtx.SaveToTextFileSparse("returnedSchur");

						SchurRHSMtx = SchurRHSselfMSR;
						//SchurRHSMtx.SaveToTextFileSparse("SchurRHSMtx");



                        //var configs = Enumerable.Repeat(op.Config, op.NoOfLevels).ToArray();

                        //                  op.B

                        //MultigridOperator mgOp = new MultigridOperator(MgBasis, this.CurrentSolution.Mapping,
                        //                   Mtx, null, configs, null);



						Console.WriteLine("Uzawa is set");
						return;
					}
				default:
                    throw new NotImplementedException("SchurOption");
            }


            //var ConvDiffInvMtx = ConvDiffInv.ToMsrMatrix();


            //// x= inv(P)*b !!!!! To be done with approximate Inverse
            // P.SpMVpara(1, B, 0, X);
        }

        public void ResetStat()
        {
            m_Converged = false;
            m_ThisLevelIterations = 0;
        }

        bool m_Converged = false;
        int m_ThisLevelIterations = 0;
		int m_IterationsInNested = 0;

		public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double>
        {

            //using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
            //    solver.DefineMatrix(P);
            //    solver.Solve(X, B);
            //}

            m_ThisLevelIterations++;

            switch (SchurOpt)
            {
                case SchurOptions.exact_matlab:
                    {
                        // Directly invert Preconditioning Matrix
                        using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
                        {
                            solver.DefineMatrix(P);
                            solver.Solve(X, B);
                        }
                        return;
                    }

                case SchurOptions.decoupledApprox:
                    {
                        SolveSubproblems(X, B);
                        return;
                    }

                case SchurOptions.SIMPLE:
                    {
                        SolveSIMPLE(X, B);
                        return;
                    }

				case SchurOptions.exact: {
                        Console.WriteLine("starting uzawa");

                        var b1 = Uidx.Select(ind => B[(int)ind]);
						var b2 = Pidx.Select(ind => B[(int)ind]);
                        var vecb1 = b1.ToArray();
						var vecb2 = b2.ToArray();

						//b1.SaveToTextFile("b1");
						//b2.SaveToTextFile("b2");

						SchurRHSMtx.SpMVpara(-1.0, vecb1, 1.0, vecb2);
                        var P = new double[Pidx.Length];
						var Usol = new double[Uidx.Length];

						//vecb2.SaveToTextFile("schurb2");

                        var Setled = (CoordinateMapping)m_mgop.BaseGridProblemMapping;

                        int D = m_mgop.Mapping.GridData.SpatialDimension;
                        var op = m_mgop;

                        var presureField = Setled.Fields[D];

                        var dummy = new DifferentialOperator(
                                       new string[] { "pressure" },
                                       new string[] { "div" },
                                       QuadOrderFunc.Linear());
                        dummy.Commit();

                        List<AggregationGridBasis[]> leveledBases = new List<AggregationGridBasis[]>();
                        List<MultigridOperator.ChangeOfBasisConfig[]> leveledConfigs = new List<MultigridOperator.ChangeOfBasisConfig[]>();



						for (var mo = op; mo != null; mo = mo.CoarserLevel) {
							Debug.Assert(mo.Mapping.AggBasis.Length == D + 1);

							AggregationGridBasis[] bases = new AggregationGridBasis[1] { mo.Mapping.AggBasis[D] };
							leveledBases.Add(bases);

							var conf = mo.Config[D];
							conf.VarIndex = new int[] { 0 };

							MultigridOperator.ChangeOfBasisConfig[] configs = new MultigridOperator.ChangeOfBasisConfig[1] { conf };
							leveledConfigs.Add(configs);
						}

						var co = presureField.Mapping;
                        var pressureMGmapping = new MultigridMapping(co, new[] { op.Mapping.AggBasis[D] }, new[] { op.Mapping.DgDegree[D] });



						List<AggregationGridBasis[]> aggBasisSeq = new List<AggregationGridBasis[]>();
						for (var mo = op; mo != null; mo = mo.CoarserLevel) {
							aggBasisSeq.Add(mo.Mapping.AggBasis);
						}

						int[] Degrees = op.BaseGridProblemMapping.BasisS.Select(b => b.Degree).ToArray();

						MultigridOperator.ChangeOfBasisConfig[][] config = new MultigridOperator.ChangeOfBasisConfig[1][];
						config[0] = new MultigridOperator.ChangeOfBasisConfig[op.BaseGridProblemMapping.BasisS.Count];
						for (int iVar = 0; iVar < config[0].Length; iVar++) {
							config[0][iVar] = new MultigridOperator.ChangeOfBasisConfig() {
								DegreeS = new int[] { Degrees[iVar] },
								mode = MultigridOperator.Mode.IdMass,
								VarIndex = new int[] { iVar }
							};
						}


                        var MsrOp =new MsrMatrix(Pidx.Length*2, Pidx.Length * 2, 1, 1);

						var rawOp = new BlockMsrMatrix(presureField.Mapping, presureField.Mapping);
						rawOp.Clear();
						var rawMaMa = new BlockMsrMatrix(presureField.Mapping, presureField.Mapping);
						rawMaMa.AccEyeSp(1);



						var MultigridOp = new MultigridOperator(leveledBases, presureField.Mapping,
						rawOp, rawMaMa, leveledConfigs,
                        dummy);

                        MultigridOp.m_RawOperatorMatrix = SchurMtx.ToBlockMsrMatrix(pressureMGmapping, pressureMGmapping);

                        var OrthoMgConfig = new OrthoMGSchwarzConfig() {
                            TargetBlockSize = 100,
                            CoarseKickIn = 200,
                            MinSolverIterations =0,
						};

                        OrthoMgConfig.ConvergenceCriterion = m_config.ConvergenceCriterion;

                        var solver = OrthoMgConfig.CreateInstance(MultigridOp);
						solver.Solve(P, vecb2);


						pGrad.SpMVpara(-1.0, P, 1.0, vecb1);

						using (var Usolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver()) {
							Usolver.DefineMatrix(ConvDiff);
							Usolver.Solve(Usol, vecb1);
						}
						
						//var MultigridOp = UniSolver.GetMultigridOperator(dummy, presureField.Mapping, leveledConfigs.ToArray());

						//                  var itP = SchurMtx.Solve_CG(P, vecb2);

						//pGrad.SpMVpara(-1.0, P, 1.0, vecb1);

						//                  var itU = ConvDiff.Solve_CG(Usol, vecb1);

						//                  Console.WriteLine($"Pressure it: {itP} Velocity it: {itU}");


						for (int i = 0; i < Uidx.Length; i++)
                            X[(int)Uidx[i]] = Usol[i];


                        for (int i = 0; i < Pidx.Length; i++)
                            X[(int)Pidx[i]] = P[i];

						//Usol.SaveToTextFile("CalculatedU");
						//P.SaveToTextFile("CalculatedP");
						//X.SaveToTextFile("CalculatedX");
                        this.m_ThisLevelIterations = solver.ThisLevelIterations;
						this.m_IterationsInNested = solver.IterationsInNested;
						this.m_Converged = solver.Converged;
						return;
					}

				case SchurOptions.Uzawa: {
						Console.WriteLine("starting uzawa");
						var b1 = Uidx.Select(ind => B[MgMap.Global2Local(ind)]);
						var b2 = Pidx.Select(ind => B[MgMap.Global2Local(ind)]);
						var vecb1 = b1.ToArray();
						var vecb2 = b2.ToArray();

						//b1.SaveToTextFile("b1");
						//b2.SaveToTextFile("b2");

						//vecb2 = schur rhs2 = b2-C*A^-1*b1
						var Ainvb1 = new double[m];
                        SolveWithMatrix(ConvDiff, Ainvb1, vecb1);
                        divVel.SpMVpara(-1.0, Ainvb1, 1.0, vecb2);
						//vecb2.SaveToTextFile("schurb2");

                        //schur res2 = schur rhs2 if the initial guess is zero vector (i.e., )
                        var res2 = vecb2.CloneAs();

						var Psol = new double[n];
						var Usol = new double[m];

                        Action<double[], double[]> multipWithPgrad = (input,output) => {
                            pGrad.SpMVpara(1.0,input,0.0, output);
                        };

						Action<double[], double[]> multipWithDivVel = (input, output) => {
							divVel.SpMVpara(-1.0, input, 0.0, output); //attention to minus sign
						};

						
						var InnerSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
                        InnerSolver.DefineMatrix(ConvDiff);

						//preconditioner matrix
						var PrecondM = new MsrMatrix(n, n, 1, 1);

                        for (int i = 0; i < n; i++) {
                            var row = pGrad.GetRowShallow((long)i);
                            double rowDotColumn = 0;

                            foreach (var e in row)
                                rowDotColumn += e.Value * e.Value;

							rowDotColumn = rowDotColumn != 0 ? rowDotColumn : 1;

							PrecondM.SetDiagonalElement(i, rowDotColumn);


						}

						using (var Psolver = new BoSSS.Solution.AdvancedSolvers.SoftPCG(false)) {
                            Psolver.ConvergenceCriterion = m_config.ConvergenceCriterion;
							Psolver.MaxIterations = m_config.MaxSolverIterations;

                            var PardisoConfig = new AdvancedSolvers.DirectSolver.Config() { WhichSolver = AdvancedSolvers.DirectSolver._whichSolver.PARDISO };

                            if (!m_mgop.MassMatrix.CheckIfUnitMatrix()) { 
                                Psolver.Precond = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
                                Psolver.Precond.DefineMatrix(pMassMatrix);
							}

							Psolver.InnerIterBefore = multipWithPgrad;
							Psolver.InnerIterAfter = multipWithDivVel;

							Psolver.InnerCycle = InnerSolver;
							Psolver.Solve(Psol, res2);

							this.m_ThisLevelIterations = Psolver.ThisLevelIterations;
							this.m_IterationsInNested = Psolver.IterationsInNested;
							this.m_Converged = Psolver.Converged;
						}

                        //Update the rhs1
						pGrad.SpMVpara(-1.0, Psol, 1.0, vecb1);

                        //Solver velocity
						using (var Usolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver()) {
							Usolver.DefineMatrix(ConvDiff);
							Usolver.Solve(Usol, vecb1);
						}

                        //Re-assign variables
						for (int i = 0; i < Uidx.Length; i++)
							X[MgMap.Global2Local(Uidx[i])] = Usol[i];


						for (int i = 0; i < Pidx.Length; i++)
							X[MgMap.Global2Local(Pidx[i])] = Psol[i];

                        //Usol.SaveToTextFile("CalculatedU");
                        //Psol.SaveToTextFile("CalculatedP");
                        //X.SaveToTextFile("CalculatedX");
                        Console.WriteLine($"total iteration number: {this.m_ThisLevelIterations}");
                        return;
					}

			}
        }



        BlockMsrMatrix GetPreconditioningMatrix() {
            var invDiagConvDiff = new MsrMatrix(m, m, 1, 1);
			var P = new MsrMatrix(n, n, 1, 1);

			for (int i = 0; i < m; i++)
                invDiagConvDiff.SetDiagonalElement(i, 1/ConvDiff[i, i]);

			var Setled = (CoordinateMapping)m_mgop.BaseGridProblemMapping;
            var VelocityFields = Setled.Fields.Take(D).ToArray();
			var PressureField = Setled.Fields.Last();

			var coP = PressureField.Mapping;
			var pressureMGmapping = new MultigridMapping(coP, new[] { MgMap.AggBasis[D] }, new[] { MgMap.DgDegree[D] });

			var coU = new CoordinateMapping(VelocityFields);
			var velocityMGmapping = new MultigridMapping(coU,  MgMap.AggBasis.Take(D).ToArray() , MgMap.DgDegree.Take(D).ToArray());

            var invDiagConvDiffBlock = ConvDiff.ToBlockMsrMatrix(velocityMGmapping, velocityMGmapping);
			//invDiagConvDiffBlock.InvertBlocks();
			var pGradBlock = pGrad.ToBlockMsrMatrix(velocityMGmapping, pressureMGmapping);
			var intermediate = pGrad.ToBlockMsrMatrix(velocityMGmapping, pressureMGmapping);
			var divVelBlock = divVel.ToBlockMsrMatrix(pressureMGmapping, velocityMGmapping);
            var PBlock = P.ToBlockMsrMatrix(pressureMGmapping, pressureMGmapping);


            BlockMsrMatrix.Multiply(intermediate, invDiagConvDiffBlock, pGradBlock);
			BlockMsrMatrix.Multiply(PBlock, divVelBlock, intermediate);
            return PBlock;
		}

		public void SolveWithMatrix<U, V>(IMutableMatrixEx M, U X, V B)
	        where U : IList<double>
	        where V : IList<double> {

			using (var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver()) {
				solver.DefineMatrix(M);
				solver.Solve(X, B);
			}
		}

		/// <summary>
		/// Solve Preconditioning Matrix in Subsystems with ConvDiff, pGrad and Schur
		/// </summary>
		/// <typeparam name="U"></typeparam>
		/// <typeparam name="V"></typeparam>
		/// <param name="X"></param>
		/// <param name="B"></param>
		public void SolveSubproblems<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double>
        {
            // For MPI
            //var idxU = Uidx[0];
            //for (int i = 0; i < Uidx.Length; i++)
            //    Uidx[i] -= idxU;
            //var idxP = Pidx[0];
            //for (int i = 0; i < Pidx.Length; i++)
            //    Pidx[i] -= idxP;

            var Bu = new double[Uidx.Length];
            var Xu = Bu.CloneAs();
            Bu = B.GetSubVector(UidxInt);
            var Bp = new double[Pidx.Length];
            var Xp = Bp.CloneAs();
            Bp = B.GetSubVector(PidxInt);

            Xu = X.GetSubVector(UidxInt);
            Xp = X.GetSubVector(PidxInt);

            ApproxAndSolveSchur(Xp, Bp);

            // Solve ConvDiff*w=v-q*pGrad
            pGrad.SpMVpara(-1, Xp, 1, Bu);

            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(ConvDiff);
                solver.Solve(Xu, Bu);
            }

            var temp = new double[Uidx.Length + Pidx.Length];

            for (int i = 0; i < Uidx.Length; i++)
            {
                temp[Uidx[i]] = Xu[i];
            }

            for (int i = 0; i < Pidx.Length; i++)
            {
                temp[Pidx[i]] = Xp[i];
            }

            X.SetV(temp);
        }


        /// <summary>
        /// Approximate the inverse of the Schur matrix and perform two Poisson solves and Matrix-Vector products. Finite elements and Fast Iterative Solvers p.383
        /// </summary>
        public void ApproxAndSolveSchur<U, V>(U Xp, V Bp)
           where U : IList<double>
           where V : IList<double>
        {

            var temp = new double[Xp.Count];
            var sol = new double[pGrad.RowPartitioning.LocalLength];

            // Poisson solve
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(PoissonMtx_T);
                solver.Solve(temp, Bp);
            }

            // Schur Convective part with scaling
            pGrad.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            sol.Clear();
            invVelMassMatrix.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            sol.Clear();
            ConvDiff.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            sol.Clear();
            invVelMassMatrixSqrt.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            divVel.SpMVpara(1, temp, 0, Xp);

            // Poisson solve
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(PoissonMtx_H);
                solver.Solve(Xp, Xp);
            }
        }

        /// <summary>
        /// Solve with a SIMPLE like approcation of the Schur complement
        /// </summary>
        /// <typeparam name="U"></typeparam>
        /// <typeparam name="V"></typeparam>
        /// <param name="X"></param>
        /// <param name="B"></param>
        public void SolveSIMPLE<U, V>(U X, V B)
           where U : IList<double>
           where V : IList<double>
        {

            var Bu = new double[Uidx.Length];
            var Xu = Bu.CloneAs();

            // For MPI
            //var idxU = Uidx[0];
            //for (int i = 0; i < Uidx.Length; i++)
            //    Uidx[i] -= idxU;
            //var idxP = Pidx[0];
            //for (int i = 0; i < Pidx.Length; i++)
            //    Pidx[i] -= idxP;

            Bu = B.GetSubVector(UidxInt);
            var Bp = new double[Pidx.Length];
            var Xp = Bp.CloneAs();
            Bp = B.GetSubVector(PidxInt);

            Xu = X.GetSubVector(UidxInt);
            Xp = X.GetSubVector(PidxInt);

            // Solve SIMPLE Schur
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(simpleSchur);
                solver.Solve(Xp, Bp);
            }


            // Solve ConvDiff*w=v-q*pGrad
            pGrad.SpMVpara(-1, Xp, 1, Bu);
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(ConvDiff);
                solver.Solve(Xu, Bu);
            }

            var temp = new double[Uidx.Length + Pidx.Length];

            for (int i = 0; i < Uidx.Length; i++)
            {
                temp[Uidx[i]] = Xu[i];
            }

            for (int i = 0; i < Pidx.Length; i++)
            {
                temp[Pidx[i]] = Xp[i];
            }

            X.SetV(temp);
        }
        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

		public long UsedMemory() {
			throw new NotImplementedException();
		}

		public void Dispose() {
            Console.WriteLine("Returned");

			// Clear all the matrices
			P?.Clear();
			ConvDiff?.Clear();
			pGrad?.Clear();
			divVel?.Clear();
			SchurMtx?.Clear();
			PoissonMtx_T?.Clear();
			PoissonMtx_H?.Clear();
			SchurConvMtx?.Clear();
			invVelMassMatrix?.Clear();
			invVelMassMatrixSqrt?.Clear();
			simpleSchur?.Clear();
			velMassMatrix?.Clear();
		}
	}
    
}
