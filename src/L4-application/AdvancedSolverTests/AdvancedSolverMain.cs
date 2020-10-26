using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using BoSSS.Solution;
using ilPSP.Connectors.Matlab;
using MPI.Wrappers;
using NUnit.Framework;
using AdvancedSolverTests.SubBlocking;

namespace AdvancedSolverTests {
    
    public class AdvancedSolverMain {

      

        public static void Main() {
            BoSSS.Solution.Application.InitMPI();
            Test();
            BoSSS.Solution.Application.FinalizeMPI();
        }

        public static void Test() {
           //SubBlockTests.LocalIndexTest(XDGusage.all,2);
            //SubBlockTests.ExternalIndexTest(XDGusage.all, 2);
            //SubBlockTests.MapConsistencyTest(XDGusage.all, 2);
            //SubBlockTests.WriteOutTestMatrices();
            //SubBlockTests.SubMatrixExtractionWithCoupling(XDGusage.all, 2,MatrixShape.diagonal_var_spec);
            //SubBlockTests.SplitVectorOperations(XDGusage.none, 2, MatrixShape.diagonal_var);
            //LocalTests.SubSelection(SelectionType.species);
            //LocalTests.CellwiseSubSelection(SelectionType.species);
            //ExternalTests.ExternalIndexTest(XDGusage.all,2);
            //ExternalTests.ExternalIndexTest(XDGusage.all, 2,4);
            //ExternalTests.GetExternalRowsTest(XDGusage.all, 2,4);
            //ExternalTests.FastSubMatrixExtraction(XDGusage.all, 2,MatrixShape.laplace,4);
            //ExternalTests.SubMatrixIgnoreCoupling(XDGusage.all, 2, MatrixShape.diagonal_var_spec,4);
            //ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var_spec,4);
            //ExternalTests.SubBlockExtraction(XDGusage.all, 2, MatrixShape.laplace, 4);
            //ExternalTests.VectorCellwiseOperation(XDGusage.all, 2, MatrixShape.diagonal_var_spec, 4);
            //ExternalTests.SubSelection(XDGusage.all, 2, MatrixShape.full_var_spec, 4);
            //ExternalTests.VectorSplitOperation(XDGusage.all, 2, MatrixShape.diagonal_var_spec, 4);
            //AdvancedSolverTests.SubBlocking.LocalTests.CellBlockVectorOperations(XDGusage.all, 2, MatrixShape.diagonal_var_spec);
            //AdvancedSolverTests.SubBlocking.LocalTests.CellwiseSubSelection(SelectionType.all_combined);
            //AdvancedSolverTests.SubBlocking.LocalTests.LocalIndexTest(XDGusage.all, 2);
            //AdvancedSolverTests.SubBlocking.LocalTests.SplitVectorOperations(XDGusage.none, 2, MatrixShape.full_var);
            //AdvancedSolverTests.SubBlocking.LocalTests.SubMatrixExtractionWithCoupling(XDGusage.all, 2, MatrixShape.full);
            //AdvancedSolverTests.SubBlocking.ExternalTests.VectorCellwiseOperation(XDGusage.none, 2, MatrixShape.diagonal_var_spec, 4);
            AdvancedSolverTests.SolverChooser.ConfigTest.TestLinearSolverConfigurations();
            
            //AdvancedSolverTests.SolverChooser.ConfigTest.TestNonLinearSolverConfigurations();
        }
    }
}
