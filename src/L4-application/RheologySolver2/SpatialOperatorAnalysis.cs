using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using System;
using BoSSS.Foundation.Quadrature;
using ilPSP.LinSolvers;
using System.Linq;
using ilPSP.Utils;
using NUnit.Framework;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.Connectors.Matlab;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Spatial operator matrix analysis
    /// </summary>
    public class SpatialOperatorAnalysis {

        /// <summary>
        /// Spatial operator matrix analysis method
        /// </summary>
        public void SpatialOperatorMatrixAnalysis(bool CheckAssertions, int AnalysisLevel) {
            using (var solver = new Rheology()) {
                int D = solver.Grid.SpatialDimension;

                if (AnalysisLevel < 0 || AnalysisLevel > 2)
                    throw new ArgumentException();


                BlockMsrMatrix OpMatrix;
                double[] OpAffine;

                solver.AssembleMatrix(out OpMatrix, out OpAffine, solver.CurrentSolution.Mapping.ToArray(), true);


                // =============================
                // AnalysisLevel 0
                // =============================
                {
                    var OpMatrixT = OpMatrix.Transpose();

                    CoordinateVector TestVec = new CoordinateVector(solver.CurrentSolution.Mapping.Fields.Select(f => f.CloneAs()).ToArray());

                    double testsumPos = 0.0;
                    double testsumNeg = 0.0;
                    for (int rnd_seed = 0; rnd_seed < 20; rnd_seed++) {

                        // fill the pressure components of the test vector
                        TestVec.Clear();
                        Random rnd = new Random(rnd_seed);
                        DGField Pressack = TestVec.Mapping.Fields[D] as DGField;
                        int J = solver.gridData.iLogicalCells.NoOfLocalUpdatedCells;
                        for (int j = 0; j < J; j++) {
                            int N = Pressack.Basis.GetLength(j);

                            for (int n = 0; n < N; n++)
                                Pressack.Coordinates[j, n] = rnd.NextDouble();
                        }

                        // Gradient times P:
                        double[] R1 = new double[TestVec.Count];
                        OpMatrix.SpMV(1.0, TestVec, 0.0, R1);       // R1 = Grad * P
                        //Console.WriteLine("L2 of 'Grad * P': " + R1.L2Norm());

                        // transpose of Divergence times P: 
                        double[] R2 = new double[TestVec.Count];
                        OpMatrix.SpMV(1.0, TestVec, 0.0, R2);      // R2 = divT * P
                        //Console.WriteLine("L2 of 'divT * P': " + R2.L2Norm());

                        TestVec.Clear();
                        TestVec.Acc(1.0, R1);
                        TestVec.Acc(1.0, R2);


                        // analyze!
                        testsumNeg += GenericBlas.L2Dist(R1, R2);

                        R2.ScaleV(-1.0);
                        testsumPos += GenericBlas.L2Dist(R1, R2);

                    }

                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (+): " + testsumPos);
                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (-): " + testsumNeg);

                    if (CheckAssertions)
                        Assert.LessOrEqual(Math.Abs(testsumNeg), testsumPos * 1.0e-13);
                }


                // =============================
                // AnalysisLevel 1 and 2
                // =============================

                if (AnalysisLevel > 0) {
                    AggregationGridBasis[][] MgBasis = AggregationGridBasis.CreateSequence(solver.MultigridSequence, solver.CurrentSolution.Mapping.BasisS);

                    MultigridOperator mgOp = new MultigridOperator(MgBasis, solver.CurrentSolution.Mapping, OpMatrix, null, solver.MultigridOperatorConfig);

                    // extract
                    ////////////

                    MsrMatrix FullMatrix = mgOp.OperatorMatrix.ToMsrMatrix();

                    MsrMatrix DiffMatrix;
                    {
                        int[] VelVarIdx = D.ForLoop(d => d);

                        int[] USubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int[] USubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int L = USubMatrixIdx_Row.Length;

                        DiffMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(int[]), USubMatrixIdx_Col, default(int[]));

                        double DiffMatrix_sd = DiffMatrix.SymmetryDeviation();
                        Console.WriteLine("Diffusion assymetry:" + DiffMatrix_sd);
                    }

                    MsrMatrix SaddlePointMatrix;
                    {
                        int[] VelPVarIdx = new int[] { 0, 1, 2 };

                        int[] VelPSubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelPVarIdx);
                        int[] VelPSubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelPVarIdx);
                        int L = VelPSubMatrixIdx_Row.Length;

                        SaddlePointMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(SaddlePointMatrix, VelPSubMatrixIdx_Row, default(int[]), VelPSubMatrixIdx_Col, default(int[]));
                    }
                    //SaddlePointMatrix.SaveToTextFileSparse("C:\\Users\\kikker\\Documents\\MATLAB\\spm.txt");

                    MsrMatrix ConstitutiveMatrix;
                    {
                        int[] StressVarIdx = new int[] { 3, 4, 5 };

                        int[] StressSubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(StressVarIdx);
                        int[] StressSubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(StressVarIdx);
                        int L = StressSubMatrixIdx_Row.Length;

                        ConstitutiveMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(ConstitutiveMatrix, StressSubMatrixIdx_Row, default(int[]), StressSubMatrixIdx_Col, default(int[]));
                    }

                    // operator analysis
                    //////////////////////

                    bool posDef;
                    if (AnalysisLevel > 1) {
                        // +++++++++++++++++++++++++++++++
                        // check condition number, etc
                        // +++++++++++++++++++++++++++++++

                        MultidimensionalArray ret = MultidimensionalArray.Create(1, 5);
                        Console.WriteLine("Calling MATLAB/Octave...");
                        using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                            bmc.PutSparseMatrix(FullMatrix, "FullMatrix");
                            bmc.PutSparseMatrix(SaddlePointMatrix, "SaddlePointMatrix");
                            bmc.PutSparseMatrix(ConstitutiveMatrix, "ConstitutiveMatrix");
                            bmc.PutSparseMatrix(DiffMatrix, "DiffMatrix");

                            bmc.Cmd("DiffMatrix = 0.5*(DiffMatrix + DiffMatrix');");

                            bmc.Cmd("condNoFullMatrix = condest(FullMatrix);");
                            bmc.Cmd("condNoSaddlePointMatrix = condest(SaddlePointMatrix);");
                            bmc.Cmd("condNoConstitutiveMatrix = condest(ConstitutiveMatrix);");
                            bmc.Cmd("condNoDiffMatrix = condest(DiffMatrix);");

                            //bmc.Cmd("eigiMaxiSaddle = 1.0; % eigs(SaddlePointMatrix,1,'lm')");
                            //bmc.Cmd("eigiMiniSaddle = 1.0; % eigs(SaddlePointMatrix,1,'sm')");
                            //bmc.Cmd("eigiMaxiConst = 1.0; % eigs(ConstitutiveMatrix,1,'lm')");
                            //bmc.Cmd("eigiMiniConst = 1.0; % eigs(ConstitutiveMatrix,1,'sm')");
                            //bmc.Cmd("eigiMaxiDiff = 1.0; % eigs(DiffMatrix,1,'lm')");
                            //bmc.Cmd("eigiMiniDiff = 1.0; % eigs(DiffMatrix,1,'sm')");

                            bmc.Cmd("lasterr");
                            bmc.Cmd("[V,r]=chol(SaddlePointMatrix);");
                            bmc.Cmd("[V,r]=chol(ConstitutiveMatrix);");
                            bmc.Cmd("ret = [condNoFullMatrix, condNoSaddlePointMatrix, condNoConstitutiveMatrix, condNoDiffMatrix, r]"); //eigiMaxiSaddle, eigiMiniSaddle, eigiMaxiConst, eigiMiniConst, eigiMaxiDiff, eigiMiniDiff,
                            bmc.GetMatrix(ret, "ret");

                            bmc.Execute(false);
                        }

                        double condNoFullMatrix = ret[0, 0];
                        double condNoSaddlePMatrix = ret[0, 1];
                        double condNoConstitutiveMatrix = ret[0, 2];
                        double condNoDiffMatrix = ret[0, 3];
                        //double eigiMaxiSaddle = ret[0, 4];
                        //double eigiMiniSaddle = ret[0, 5];
                        //double eigiMaxiConst = ret[0, 6];
                        //double eigiMiniConst = ret[0, 7];
                        //double eigiMaxiDiff = ret[0, 8];
                        //double eigiMiniDiff = ret[0, 9];
                        posDef = ret[0, 4] == 0;

                        //Console.WriteLine("Eigenvalue range of saddle point matrix: {0} to {1}", eigiMiniSaddle, eigiMaxiSaddle);
                        //Console.WriteLine("Eigenvalue range of constitutive matrix: {0} to {1}", eigiMiniConst, eigiMaxiConst);
                        //Console.WriteLine("Eigenvalue range of diffusion matrix: {0} to {1}", eigiMiniDiff, eigiMaxiDiff);

                        Console.WriteLine("Condition number full operator: {0:0.####E-00}", condNoFullMatrix);
                        Console.WriteLine("Condition number saddle point operator: {0:0.####E-00}", condNoSaddlePMatrix);
                        Console.WriteLine("Condition number constitutive operator: {0:0.####E-00}", condNoConstitutiveMatrix);
                        Console.WriteLine("Condition number diffusion operator: {0:0.####E-00}", condNoDiffMatrix);

                        //base.QueryHandler.ValueQuery("ConditionNumber", condNoFullMatrix);

                    } else {
                        // +++++++++++++++++++++++++++++++++++++++
                        // test only for positive definiteness
                        // +++++++++++++++++++++++++++++++++++++++

                        var SaddlePMatrixFull = SaddlePointMatrix.ToFullMatrixOnProc0();
                        var ConstMatrixFull = ConstitutiveMatrix.ToFullMatrixOnProc0();


                        posDef = true;
                        try {
                            SaddlePMatrixFull.Cholesky();
                        } catch (ArithmeticException) {
                            posDef = false;
                        }

                        posDef = true;
                        try {
                            ConstMatrixFull.Cholesky();
                        } catch (ArithmeticException) {
                            posDef = false;
                        }
                    }


                    double SaddlePSymm = SaddlePointMatrix.SymmetryDeviation();
                    Console.WriteLine("Symmetry deviation of saddle point matrix: " + SaddlePSymm);

                    if (posDef)
                        Console.WriteLine("Good news: Saddle point operator matrix seems to be positive definite.");
                    else
                        Console.WriteLine("WARNING: Saddle point operator matrix is not positive definite.");


                    double ConstSymm = ConstitutiveMatrix.SymmetryDeviation();
                    Console.WriteLine("Symmetry deviation of constitutive matrix: " + ConstSymm);

                    if (posDef)
                        Console.WriteLine("Good news: constitutive operator matrix seems to be positive definite.");
                    else
                        Console.WriteLine("WARNING: constitutive operator matrix is not positive definite.");

                    //if (CheckAssertions) {
                    //    if (Control.AdvancedDiscretizationOptions.ViscosityMode == ViscosityMode.FullySymmetric && Control.PhysicalParameters.IncludeConvection == false) {
                    //        Assert.IsTrue(posDef, "Positive definiteness test failed.");
                    //        double compVal = DiffMatrix.InfNorm() * 1e-13;
                    //        Assert.LessOrEqual(DiffSymm, compVal, "Diffusion matrix seems to be non-symmetric.");
                    //    }
                    //}
                }
            }
        }
    }
}