using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Tests {
    struct SolutionTestField {
        public string Name;
        public double AcceptableL2Error;
        public DGField ExactSolution;

        public SolutionTestField(string name, DGField exactSolution) {
            Name = name;
            AcceptableL2Error = -1;
            ExactSolution = exactSolution;
        }

        public SolutionTestField(string name, double acceptableL2Error, DGField exactSolution) {
            Name = name;
            AcceptableL2Error = acceptableL2Error;
            ExactSolution = exactSolution;
        }
    }

    struct ResidualTestField {
        public string Name;
        public double AcceptableResidual;

        public ResidualTestField(string name) {
            Name = name;
            AcceptableResidual = -1;
        }

        public ResidualTestField(string name, double acceptableResidual) {
            Name = name;
            AcceptableResidual = acceptableResidual;
        }
    }

    class Test<T> where T : AppControl {

        public IList<SolutionTestField> SolutionTestFields { get; private set; }
        public IList<ResidualTestField> ResidualTestFields { get; private set; }

        public T Control { get; private set; }

        public Test(T control) {
            this.Control = control;
            SolutionTestFields = new List<SolutionTestField>();
            ResidualTestFields = new List<ResidualTestField>();
        }

        public void AddSolutionTestField(SolutionTestField field) {
            SolutionTestFields.Add(field);
        }

        public void AddResidualTestField(ResidualTestField field) {
            ResidualTestFields.Add(field);
        }
    }

    class ConvergenceTest<T> where T : AppControl {
        public SortedList<double, Test<T>> Tests { get; private set; }

        public bool UseExactSolution { get; private set; }

        public IList<double> ExpectedSlopes { get; private set; }

        public ConvergenceTest(SortedList<double, Test<T>> tests, bool useExactSolution, IList<double> expectedSlopes) {

            Tests = tests;
            UseExactSolution = useExactSolution;
            ExpectedSlopes = expectedSlopes;
        }

        public void AddTest(double characteristicLength, Test<T> test) {
            Tests.Add(characteristicLength, test);
        }
    }

    static class TestRunner<S, T> where S : ApplicationWithSolver<T>, new() where T : AppControlSolver, new() {

        public static void SolverTest(Test<T> test) {
            using(var solver = new S()) {
                solver.Init(test.Control);
                solver.RunSolverMode();
                IDictionary<string, DGField> solutionFields = GetSolutionFields(solver);
                bool success = EvaluateL2Error(test.SolutionTestFields, solutionFields);

                IDictionary<string, DGField> residualFields = GetResidualFields(solver);
                success &= EvaluateResidual(test.ResidualTestFields, residualFields);
                Assert.IsTrue(success);
            }
        }

        static Dictionary<string, DGField> GetSolutionFields(S solver) {
            IList<DGField> variableFields = solver.CurrentState.Fields;
            IList<DGField> parameterFields = solver.Parameters;
            Dictionary<string, DGField> solutionFields = new Dictionary<string, DGField>(variableFields.Count + parameterFields.Count);
            //Which one is the correct name?
            foreach(DGField dgF in variableFields) {
                solutionFields.Add(dgF.Identification, dgF);
            }
            foreach(DGField dgF in parameterFields) {
                solutionFields.Add(dgF.Identification, dgF);
            }
            return solutionFields;
        }

        static Dictionary<string, DGField> GetResidualFields(S solver) {
            IList<DGField> variableFields = solver.CurrentResidual.Fields;
            Dictionary<string, DGField> solutionFields = new Dictionary<string, DGField>(variableFields.Count);
            //Which one is the correct name?
            foreach(DGField dgF in variableFields) {
                solutionFields.Add(dgF.Identification, dgF);
            }
            return solutionFields;
        }

        static bool EvaluateL2Error(IList<SolutionTestField> testFields, IDictionary<string, DGField> solutionFields) {
            bool success = true;
            foreach(SolutionTestField field in testFields) {
                DGField solutionField = solutionFields[field.Name];
                double error = solutionField.L2Error(field.ExactSolution);
                Console.Write("L2 error, '{0}': \t{1}", field.Name, error);
                if(field.AcceptableL2Error >= 0 && error > field.AcceptableL2Error) {
                    Console.WriteLine("   Above Threshold (" + field.AcceptableL2Error + ")");
                    success = false;
                }
            }
            return success;
        }

        static bool EvaluateResidual(IList<ResidualTestField> testFields, IDictionary<string, DGField> residualFields) {
            bool success = true;
            foreach(ResidualTestField field in testFields) {
                DGField residualField = residualFields[field.Name];
                double residual = residualField.L2Norm();
                Console.Write("L2 norm, '{0}': \t{1}", field.Name, residual);
                if(field.AcceptableResidual >= 0 && residual > field.AcceptableResidual) {
                    Console.WriteLine("   Above Threshold (" + field.AcceptableResidual + ")");
                    success = false;
                }
            }
            return success;
        }

        static void ConvergenceTest<T>(ConvergenceTest<T> test) where T : AppControl {

        }
    }
}
