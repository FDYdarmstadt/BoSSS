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
using System.Threading.Tasks;

namespace BoSSS.Application.Rheology {
    class OldVersionMain {

        //// Operator test
        //static public RheologyControl OperatorTest() {
        //    RheologyControl C = new RheologyControl();

        //    // Solver Options
        //    C.ViscousPenaltyScaling = 1;
        //    C.NoOfTimesteps = 1;
        //    C.savetodb = false;
        //    C.DbPath = @"C:\AnnesBoSSSdb\OperatorTest";
        //    C.ProjectName = "OperatorTest";
        //    C.Stokes = true;
        //    C.FixedStreamwisePeriodicBC = true;

        //    // Create Fields
        //    C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        //    //C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        //    C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        //    C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        //    C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        //    C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

        //    // Create Grid
        //    C.GridFunc = delegate {
        //        var _xNodes = GenericBlas.Linspace(-1, 1, 201);

        //        var grd = Grid1D.LineGrid(_xNodes, C.FixedStreamwisePeriodicBC);

        //        return grd;
        //    };

        //    // Set Initial Conditions
        //    C.InitialValues_Evaluators.Add("ExternalStressXX", X => (X[0] * X[0]));//Math.Cos(X[0]));
        //    C.InitialValues_Evaluators.Add("ExternalStressXY", X => 0); //Math.Sin(X[0]));
        //    C.InitialValues_Evaluators.Add("ExternalStressYY", X => 0);// Math.Sin(X[0]));

        //    return C;
        //}


        //______________________________________________________________________________________________________________
        // Extra Stresses Parameter (external, 2D)
        //[InstantiateFromControlFile("ExternalStressXX", "StressXX", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField ExternalStressXX;

        //[InstantiateFromControlFile("ExternalStressXY", "StressXY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField ExternalStressXY;

        //[InstantiateFromControlFile("ExternalStressYY", "StressYY", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField ExternalStressYY;


        // OLD VERSION
        //______________________________________________________________________________________________________________

        //dt = 0;

        //MsrMatrix OpMatrix; // Sparse Matrix M
        //double[] OpAffine;  // RHS b

        //SpatialOperatorMatrixAnalysis(false, this.Control.AnalysisLevel);

        ////STOKES
        ////===================================================================

        //if (this.Control.Stokes == true) {

        //    AssembleMatrix(out OpMatrix, out OpAffine, this.Velocity);

        //    //Solving with MATLAB
        //    //_________________________________________________________________________________________________
        //    //double[] RHS = OpAffine.CloneAs();              //build RHS b
        //    //RHS.ScaleV(-1.0);                               // multiply with (-1), since M x u + b = 0
        //    //OpMatrix.SolveMATLAB(CurrentSolution, RHS);
        //    //_________________________________________________________________________________________________

        //    // Solving with FullMatrix
        //    //_________________________________________________________________________________________________
        //    //var fullOpMatrix = OpMatrix.ToFullMatrixOnProc0();
        //    //double[] RHS = OpAffine.CloneAs();            //build RHS b
        //    //RHS.ScaleV(-1.0);                             // multiply with (-1), since M x u + b = 0
        //    //var test = CurrentSolution.ToArray();
        //    //fullOpMatrix.Solve(test, RHS);

        //    //CurrentSolution.ClearEntries();
        //    //CurrentSolution.AccV(1.0, test);
        //    //_________________________________________________________________________________________________

        //    // Solving with MUMPS
        //    //_________________________________________________________________________________________________
        //    using (var slv = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
        //        double[] RHS = OpAffine.CloneAs();        //build RHS b
        //        RHS.ScaleV(-1.0);                         // multiply with (-1), since M x u + b = 0
        //        slv.DefineMatrix(OpMatrix);               // build sparse matrix M
        //        slv.Solve(CurrentSolution, RHS);          // Solve System M x u + b = 0
        //    }

        //    //_________________________________________________________________________________________________


        //    // FOR OPERATOR TEST ONLY
        //    //_________________________________________________________________________________________________

        //var TAU_ext = new CoordinateVector(Velocity[0], Velocity[1], Pressure, ExternalStressXX, ExternalStressXY, ExternalStressYY).ToArray();
        //DGField result = new SinglePhaseField(new Basis(GridData, 2));
        //var test = new CoordinateVector(result, ExternalStressXY, ExternalStressYY);


        ////OpMatrix.SaveToTextFile("FullMatrix");
        //OpMatrix.SpMV(1, m_CurrentSolution, 0, test);

        //double test_norm = test.L2Norm();


        //Gnuplot Result = new Gnuplot();
        //Result.PlotField(result, 10);
        //Result.Execute();

        //Console.WriteLine(test_norm);


        //    //_________________________________________________________________________________________________

        //} else {

        //    //NAVIER_STOKES (STEADY)
        //    //=================================================================================

        //    int NoOfIterations = 0;
        //    bool LastIteration_Converged = false;
        //    double ResidualNorm;

        //    // initial guess and its residual
        //    AssembleMatrix(out OpMatrix, out OpAffine, this.Velocity);

        //    double[] RHS = OpAffine.CloneAs();              //build RHS b
        //    RHS.ScaleV(-1.0);                               // multiply with (-1), since M x u + b = 0

        //    CurrentResidual.Clear();
        //    CurrentResidual.Acc(-1.0, OpAffine);


        //    //OpMatrix.SaveToTextFile("FullMatrix");


        //    OpMatrix.SpMVpara(-1.0, CurrentSolution, 1.0, CurrentResidual);


        //    double Cond = OpMatrix.condest();
        //    //base.QueryHandler.ValueQuery("ConditionNumber", Cond);
        //    //Console.WriteLine("ConditionNumber:" + Cond);
        //    //double Rank = OpMatrix.rank();
        //    //base.QueryHandler.ValueQuery("Rank", Rank);
        //    //OpMatrix.SaveToTextFile("SparseMatrix");


        //    ResidualNorm = CurrentResidual.L2NormPow2().MPISum().Sqrt();
        //    Console.WriteLine("{0}:\t{1:0.####E-00}", NoOfIterations, ResidualNorm);



        //    // iterate...
        //    while ((!LastIteration_Converged && NoOfIterations < this.Control.MaxIter) || (NoOfIterations < this.Control.MinIter)) {
        //        NoOfIterations++;

        //        //Solving with MUMPS
        //        //_________________________________________________________________________________________________
        //        //using (var slv = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
        //        //    slv.DefineMatrix(OpMatrix);               // build sparse matrix M
        //        //    slv.Solve(CurrentSolution, RHS);
        //        //}
        //        //_________________________________________________________________________________________________


        //        //Solving with MATLAB
        //        //_________________________________________________________________________________________________
        //        OpMatrix.SolveMATLAB(CurrentSolution, RHS);
        //        //_________________________________________________________________________________________________


        //        // update linearization
        //        this.AssembleMatrix(out OpMatrix, out OpAffine, this.Velocity);
        //        RHS = OpAffine.CloneAs();
        //        RHS.ScaleV(-1.0);

        //        // residual evaluation
        //        CurrentResidual.Clear();
        //        CurrentResidual.Acc(-1.0, OpAffine);
        //        OpMatrix.SpMVpara(-1.0, CurrentSolution, 1.0, CurrentResidual);

        //        ResidualNorm = CurrentResidual.L2NormPow2().MPISum().Sqrt();
        //        //Console.WriteLine("{0}:\t{1:0.####E-00}", NoOfIterations, ResidualNorm);
        //        base.ResLogger.CustomValue(ResidualNorm, "Overall");
        //        base.ResLogger.NextIteration(true);

        //        PlotCurrentState(0, NoOfIterations, 0);

        //        if (ResidualNorm < Control.ConvCrit)
        //            LastIteration_Converged = true;
        //    }
        //}


        //ComputeResidual();

        //base.ResLogger.NextTimestep(false);

        //return dt;
        //______________________________________________________________________________________________________________

        //private void ComputeResidual() {

        //    // matrix assembly
        //    MsrMatrix OpMatrix;
        //    double[] OpAffine;
        //    AssembleMatrix(out OpMatrix, out OpAffine, this.Velocity.Current);

        //    // evaluate the operator

        //    CurrentResidual.Clear();
        //    CurrentResidual.Acc(-1.0, OpAffine);
        //    OpMatrix.SpMVpara(-1.0, CurrentSolution, 1.0, CurrentResidual);

        //    // print 
        //    foreach (var F in CurrentResidual.Mapping.Fields) {
        //        Console.WriteLine("L2 norm of residual for '" + F.Identification + "' = " + F.L2Norm());
        //    }
        //}
    }
}
