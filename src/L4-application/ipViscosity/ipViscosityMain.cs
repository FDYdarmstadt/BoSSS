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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Solvers;
using ilPSP.LinSolvers;
using System.Globalization;
using ilPSP.Tracing;
using NUnit.Framework;
using System.IO;
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using System.Linq;
using BoSSS.Foundation.IO;
using MPI.Wrappers;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.ipViscosity {

    /// <summary>
    /// What should the test program do for you?
    /// </summary>
    public enum TestMode {
        /// <summary>
        /// only check the residual of the known solution
        /// </summary>
        CheckResidual,

        /// <summary>
        /// run a solver
        /// </summary>
        Solve
    }

    /// <summary>
    /// which terms should be invoked?
    /// </summary>
    [Flags]
    public enum Terms {

        /// <summary>
        /// \f[ 
        ///   - \operatorname{div} \left( \mu \nabla u_d \right)
        /// \f]
        /// </summary>
        T1 = 0x1,

        /// <summary>
        /// \f[ 
        ///   - \operatorname{div} \left( \mu (\partial_d \vec{u}) \right)
        /// \f]
        /// </summary>
        T2 = 0x2,

        /// <summary>
        /// \f[ 
        ///    \frac{2}{3} \operatorname{div} \left( \mu \myMatrix{I} \operatorname{div} ( \vec{u} )  \right)
        /// \f]
        /// </summary>
        T3 = 0x4
    }


    class ipViscosityMain : BoSSS.Solution.Application {

       
        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            
            BoSSS.Solution.Application._Main(args, true, delegate() {
                ipViscosityMain p = new ipViscosityMain(); 
                return p;
            });
            
            //BoSSS.Solution.Application.InitMPI(new string[0]);
            //_Test.ConsistencyTest(Terms.T1, ViscosityImplementation.H, 1);
            //csMPI.Raw.mpiFinalize();
        }

        internal TestMode mode = TestMode.CheckResidual;

        internal TestSolution solution = new Polynomial2D_VariableVisc();

        internal int PolynomialDegree = 2;

        internal ITestGrid grid = new MixedBcGrid();

        internal Terms whichTerms = Terms.T1;

        /// <summary>
        /// dependent variable
        /// </summary>
        VectorField<SinglePhaseField> U;

        /// <summary>
        /// Viscosity terms for <see cref="U"/>, by linear evaluation.
        /// </summary>
        VectorField<SinglePhaseField> ViscU_linear;

        /// <summary>
        /// Viscosity terms for <see cref="U"/>, by nonlinear evaluation.
        /// </summary>
        VectorField<SinglePhaseField> ViscU_nonlinear;

        /// <summary>
        /// RHS
        /// </summary>
        VectorField<SinglePhaseField> RHS;

        /// <summary>
        /// Residual
        /// </summary>
        VectorField<SinglePhaseField> Residual;

        /// <summary>
        /// inhomogeneous boundary conditions
        /// </summary>
        VectorField<SinglePhaseField> bnd;

        /// <summary>
        /// error
        /// </summary>
        VectorField<SinglePhaseField> Err;


        /// <summary>
        /// molecular viscosity
        /// </summary>
        SinglePhaseField mu;


        protected override void CreateFields() {
            var b = new Basis(this.GridData, this.PolynomialDegree);
            int D = this.GridData.SpatialDimension;

            var _U = new SinglePhaseField[D] ;
            for(int d = 0; d < D; d++) 
                _U[d] = new SinglePhaseField(b,VariableNames.Velocity_d(d));
            this.U = new VectorField<SinglePhaseField>(_U);
            this.ViscU_nonlinear = new VectorField<SinglePhaseField>(D, b, "ViscU_nonlinear", SinglePhaseField.Factory);
            this.ViscU_linear = new VectorField<SinglePhaseField>(D, b, "ViscU_linear", SinglePhaseField.Factory);

            this.RHS = new VectorField<SinglePhaseField>(D, b, "RHS", SinglePhaseField.Factory);
            this.bnd = new VectorField<SinglePhaseField>(D, b, "bnd", (c, s) => (new SinglePhaseField(c, s)));
            this.Residual = new VectorField<SinglePhaseField>(D, b, "Residual", (c, s) => (new SinglePhaseField(c, s)));
            this.Err = new VectorField<SinglePhaseField>(D, b, "Error", (c, s) => (new SinglePhaseField(c, s)));
            this.mu = new SinglePhaseField(b, VariableNames.ViscosityMolecular);
        }
        
        IncompressibleBoundaryCondMap BcMap;

        protected override GridCommons CreateOrLoadGrid() {
            GridCommons grd = grid.GetGrid();
            return grd;
        }

        SpatialOperator Operator;
        MsrMatrix OperatorMtx;
        
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            using(FuncTrace tr = new FuncTrace()) {

                this.BcMap = new IncompressibleBoundaryCondMap(this.GridData, grid.GetBoundaryConfig(), PhysicsMode.Incompressible);


                // assemble system, create matrix
                // ------------------------------

                var volQrSch = new CellQuadratureScheme(true, CellMask.GetFullMask(this.GridData));
                var edgQrSch = new EdgeQuadratureScheme(true, EdgeMask.GetFullMask(this.GridData));
                
                //var volQrSch = new CellQuadratureScheme(true, CellMask.GetEmptyMask(this.GridData));
                //var edgQrSch = new EdgeQuadratureScheme(true, EdgeMask.GetEmptyMask(this.GridData));
                //var edgQrSch = new EdgeQuadratureScheme(true, this.GridDat.BoundaryEdges);
                //var edgQrSch = new EdgeQuadratureScheme(true, this.GridDat.BoundaryEdges.Complement());

                int D = GridData.SpatialDimension;
                double penalty_base = ((double)((U[0].Basis.Degree + 1) * (U[0].Basis.Degree + D))) / ((double)D);
                double penalty_factor = 1.2;



                // equation assembly
                // -----------------
                string[] CodNames = D.ForLoop(i => "C" + i);
                Operator = new SpatialOperator(VariableNames.VelocityVector(D), new string[] { VariableNames.ViscosityMolecular }, CodNames,QuadOrderFunc.Linear());

                for(int d = 0; d < D; d++) {
                    if((this.whichTerms & Terms.T1) != 0) {
                        var flx1 = new swipViscosity_Term1(penalty_base * penalty_factor, this.GridData.Cells.cj, d, D, BcMap,ViscosityOption.VariableViscosity);

                        flx1.g_Diri_Override = this.solution.U;
                        flx1.g_Neu_Override = this.solution.dU;
                        Operator.EquationComponents[CodNames[d]].Add(flx1);
                    }
                    if((this.whichTerms & Terms.T2) != 0) {
                        var flx2 = new swipViscosity_Term2(penalty_base * penalty_factor, this.GridData.Cells.cj, d, D, BcMap, ViscosityOption.VariableViscosity);

                        flx2.g_Diri_Override = this.solution.U;
                        flx2.g_Neu_Override = this.solution.dU;
                        Operator.EquationComponents[CodNames[d]].Add(flx2);
                    }
                    if((this.whichTerms & Terms.T3) != 0) {
                        var flx3 = new swipViscosity_Term3(penalty_base * penalty_factor, this.GridData.Cells.cj, d, D, BcMap, ViscosityOption.VariableViscosity);

                        flx3.g_Diri_Override = this.solution.U;
                        flx3.g_Neu_Override = this.solution.dU;
                        Operator.EquationComponents[CodNames[d]].Add(flx3);
                    }
                } // */
                Operator.Commit();


                var map = this.U.Mapping;
                OperatorMtx = new MsrMatrix(map, map);
                Operator.ComputeMatrixEx(map, new DGField[] { this.mu }, map,
                                         OperatorMtx, this.bnd.CoordinateVector,
                                         volQuadScheme: volQrSch, edgeQuadScheme: edgQrSch);

                // test for matrix symmetry
                // ========================

                if(base.MPISize == 1) {
                    double MatrixAssymmetry = OperatorMtx.SymmetryDeviation();
                    Console.WriteLine("Matrix asymmetry: " + MatrixAssymmetry);
                    Assert.LessOrEqual(Math.Abs(MatrixAssymmetry), 1.0e-10);
                }
            }
        }

        /// <summary>
        /// checks whether the linear and nonlinear implementation of operator evaluation are mathematically equal
        /// </summary>
        void LinearNonlinComparisonTest() {
            int L = this.bnd.CoordinateVector.Count();

            // need to assure to use the same quadrature oder on both evaluation variants
            var volQrSch = (new CellQuadratureScheme(false, CellMask.GetFullMask(this.GridData)))
                                    .AddFixedOrderRules(this.GridData, this.PolynomialDegree*3);
            var edgQrSch = new EdgeQuadratureScheme(true, EdgeMask.GetFullMask(this.GridData))
                                    .AddFixedOrderRules(this.GridData, this.PolynomialDegree*3);

            //var volQrSch = new CellQuadratureScheme(true, CellMask.GetEmptyMask(this.GridData));
            //var edgQrSch = new EdgeQuadratureScheme(true, EdgeMask.GetEmptyMask(this.GridData));
            //var edgQrSch = new EdgeQuadratureScheme(true, this.GridData.BoundaryEdges)
            //                        .AddFixedOrderRules(this.GridData, this.PolynomialDegree * 3);
            //var edgQrSch = new EdgeQuadratureScheme(true, this.GridData.BoundaryEdges.Complement())
            //                          .AddFixedOrderRules(this.GridData, this.PolynomialDegree * 3);


            
            for(int run = 0; run < 1; run++) {

                // setup a random test vector
                Random rnd = new Random();
                var TestArgument = this.bnd.CloneAs().CoordinateVector;
                for(int i = 0; i < L; i++) {
                    TestArgument[i] = rnd.NextDouble();
                }
                
                Stopwatch lin = new Stopwatch();
                Stopwatch nol = new Stopwatch();
                
                // linear evaluation
                CoordinateVector LinResult = this.ViscU_linear.CoordinateVector;
                LinResult.Clear();
                lin.Start();
                {
                    var map = this.U.Mapping;
                    var tempOperatorMtx = new MsrMatrix(map, map);
                    var tempAffine = new double[L];
                    Operator.ComputeMatrixEx(map, new DGField[] { this.mu }, map,
                                             tempOperatorMtx, tempAffine,
                                             volQuadScheme: volQrSch, edgeQuadScheme: edgQrSch);


                    tempOperatorMtx.SpMVpara(1.0, TestArgument, 0.0, LinResult);
                    LinResult.AccV(1.0, tempAffine);
                }
                lin.Stop();

                // nonliner evaluation
                CoordinateVector NolResult = this.ViscU_nonlinear.CoordinateVector;
                NolResult.Clear();
                nol.Start();
                {
                    var evaluator = Operator.GetEvaluatorEx(TestArgument.Mapping, new DGField[] { this.mu }, this.Residual.Mapping,
                        volQrCtx: volQrSch, edgeQrCtx: edgQrSch);
                    evaluator.Evaluate(1.0, 0.0, NolResult);
                }
                nol.Stop();

                double L2Dist = GenericBlas.L2DistPow2(LinResult, NolResult).MPISum().Sqrt();
                Console.WriteLine("L2 dist of linear/Nonlinear evaluation comparison: {0}", L2Dist);
                
                LinResult.Acc(-1.0, NolResult);

                foreach(SinglePhaseField DGfield in LinResult.Mapping.Fields) {
                    for(int p = 0; p <= DGfield.Basis.Degree; p++) {
                        double L2err_p = DGfield.L2NormPerMode(p);
                        Console.WriteLine("   ERR{2} {1} \t{0}", L2err_p, DGfield.Identification, p);
                    }
                }

                Console.WriteLine("Time linear {0}, time nonlinear: {1}", lin.Elapsed, nol.Elapsed);

                Assert.LessOrEqual(L2Dist, 1.0e-4, "L2 distance between linear and nonlinear evaluation of the same flux.");
                
            }
        }



        protected override void SetInitial() {
            int D = this.GridData.SpatialDimension;

            var sol = this.solution;

            for (int d = 0; d < D; d++) {
                this.U[d].ProjectField((Func<double[], double>)(X => (sol.U(d, X))));

                if ((this.whichTerms & Terms.T1) != 0)
                    this.RHS[d].ProjectField(1.0, ((Func<double[], double>)(X => (sol.Term1(d, X)))).Vectorize());
                if ((this.whichTerms & Terms.T2) != 0)
                    this.RHS[d].ProjectField(1.0, ((Func<double[], double>)(X => (sol.Term2(d, X)))).Vectorize());
                if ((this.whichTerms & Terms.T3) != 0)
                    this.RHS[d].ProjectField(1.0, ((Func<double[], double>)(X => (sol.Term3(d, X)))).Vectorize());
            }

            this.mu.ProjectField(sol.mu);
        }




        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                base.NoOfTimesteps = -1;
                base.TerminationKey = true;
                int D = this.GridData.SpatialDimension;

                this.Err.Clear();
                this.Err.Acc(-1.0, this.U);

                // compare linear vs. nonlinear evaluation
                // =======================================

                LinearNonlinComparisonTest();
                
                // call solver
                // ===========
                if (this.mode == TestMode.Solve) {

                    var solver = new ilPSP.LinSolvers.monkey.CG();
                    solver.MaxIterations = 20000;
                    solver.Tolerance = 1.0e-12;
                    solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.CPU;

                    solver.DefineMatrix(this.OperatorMtx);

                    double[] _rhs = this.RHS.CoordinateVector.ToArray();
                    _rhs.AccV(-1.0, this.bnd.CoordinateVector);

                    var solRes = solver.Solve(this.U.CoordinateVector, _rhs);

                    Console.WriteLine("sparse solver result: " + solRes.ToString());
                    Assert.IsTrue(solRes.Converged);

                } else if (this.mode == TestMode.CheckResidual) {
                    // residual computation is done anyway...
                } else {
                    throw new NotImplementedException();
                }

                // Residual computation
                // ====================
                {
                    this.Residual.Clear();
                    this.Residual.Acc(1.0, this.bnd);
                    this.OperatorMtx.SpMVpara(1.0, this.U.CoordinateVector, 1.0, this.Residual.CoordinateVector);

                    this.Residual.Acc(-1.0, this.RHS);

                    L2ResidualNorm = D.ForLoop(delegate(int i) {
                        var f = this.Residual[i];
                        double L2Norm = f.L2Norm();
                        Console.WriteLine("L2 " + f.Identification + ": " + L2Norm);
                        return L2Norm;
                    });
                }

                // Error computation
                // =================
                {
                    this.Err.Acc(1.0, this.U);

                    L2ErrorNorm = D.ForLoop(delegate(int i) {
                        var f = this.Err[i];
                        double L2Norm = f.L2Norm();
                        Console.WriteLine("L2 " + f.Identification + ": " + L2Norm);
                        return L2Norm;
                    });
                }

                return 0.0;
            }
        }

        internal double[] L2ErrorNorm;

        internal double[] L2ResidualNorm;


        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0) {
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(
                ArrayTools.Cat<DGField>(this.U, this.RHS, this.bnd, this.Residual, this.Err, this.mu, this.ViscU_linear, this.ViscU_nonlinear), "ipViscosity", phystime, superSampling);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(
                ArrayTools.Cat<DGField>(this.U, this.RHS, this.bnd, this.Residual, this.Err, this.mu, this.ViscU_linear, this.ViscU_nonlinear), "grid", phystime, 0);
        }
    }
}