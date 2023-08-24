using ApplicationWithIDT;
using ApplicationWithIDT.OptiLevelSets;
using ilPSP.Utils;
using NUnit.Framework;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using System;
using BUIDT.Fluxes;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using System.Linq;

namespace BUIDT.Tests
{
    [TestFixture]
    public static class BUIDTTestProgram
    {
        #region NUnit stuff
        [OneTimeSetUp]
        static public void Init() {
            BoSSS.Solution.Application.InitMPI();
        }

        [OneTimeTearDown]
        static public void Cleanup() {
        }
        #endregion

        [Test]
        /// <summary>
        /// Test for the StraightShock presented at Eccomas22
        /// </summary>
        public static void StraightShockCurvedStart_Eccomas22()
        {
            using (var p = new BUIDTMain())
            {
                var C = BUIDTHardCodedControl.StraightShockCurvedStart_Eccomas22(
                    //dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\BUIDT\AcceleratingShock\BUIDT_db",
                    dbPath: null,
                    MaxIterations: 50,
                    dgDegree: 0,
                    numOfCellsX: 10,
                    numOfCellsY: 10,
                    linearization: Linearization.FD,
                    agg:0.1,
                    ImmediatePlotPeriod: -1
                    );
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-09 && p.ResidualVector.MPI_L2Norm() < 1e-09), System.String.Format("the L2 Error is greater than 1e-10 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(),p.obj_f_vec.MPI_L2Norm()));

            }
        }

        [Test]
        /// <summary>
        /// Test for the AcceleratingShock
        /// </summary>
        public static void AcceleratingShock() {
            using(var p = new BUIDTMain()) {
                var C = BUIDTHardCodedControl.AccShock(
                    //dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\BUIDT\AcceleratingShock\BUIDT_db",
                    dbPath:null,
                    MaxIterations: 100,
                    dgDegree: 3,
                    numOfCellsX: 10,
                    numOfCellsY: 10,
                    OptiNumOfCellsX: 10,
                    OptiNumOfCellsY: 10,
                    linearization: Linearization.FD,
                    agg: 0.4,
                    ImmediatePlotPeriod: -1,
                    optiLevelSetType: OptiLevelSetType.SplineLevelSet,
                    getLevelSet: GetLevelSet.FromFunction,
                    applyReInit: true,
                    solverRunType: SolverRunType.Staggerd
                    ) ;
                p.Init(C);
                p.RunSolverMode();
                
                
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-01 && p.ResidualVector.MPI_L2Norm() < 1e-02), System.String.Format("the L2 Error is greater than 1e-10 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));

            }
        }
        /// <summary>
        /// Utility function to save the Explicit surface described by a spline LevelSet
        /// </summary>
        /// <param name="p"></param>
        /// <param name="filename"></param>
        public static void SaveIsoContourToTextFile(BUIDTMain p, string filename) {
            if(p.LevelSetOpti is SplineOptiLevelSet spliny) {
                spliny.GetSpline();
                if(spliny.Spline is CubicSpline cSpliny) {

                    var yPoints = GenericBlas.Linspace(0, 1.2, 100);
                    var xPoints = new double[100];
                    var cPoints = new double[100];
                    var allPoints = MultidimensionalArray.Create(100, 3);
                    for(int i = 0; i < 100; i++) {
                        allPoints[i, 0] = yPoints[i];

                        xPoints[i] = cSpliny.Interpolate(yPoints[i]);
                        allPoints[i, 1] = xPoints[i];

                        cPoints[i] = cSpliny.Differentiate(yPoints[i]);
                        allPoints[i, 2] = cPoints[i];


                    }
                    allPoints.SaveToTextFile(filename);

                }
            }
        }
        /// <summary>
        /// does a convergence study for the l2 projection for the Accelerating shock
        /// </summary>
        public static void AccShockXDGConvergenceStudy() {
            double mu1 = 4.0;
            double mu2 = 3.0;
            double ShockSpeed(double[] x) {
                return x[0] - ((mu1 / mu2 + 1) * (1 - Math.Sqrt(1 + mu2 * x[1])) + mu1 * x[1]);
            }
            double ExactSolUnsmooth(double[] X,string spc) {
                if(spc =="L") {
                    return mu1;
                } else {
                    return mu2 * (X[0] - 1) / (mu2 * X[1] + 1);
                }
            }
            var results = DoConverGenceStudy(LSDegree:6, ExactSolUnsmooth, ShockSpeed,
                                                        nPol:6,nGrid:5,is_nf_smth:false,NonLinQuadDegree:12);
        }
        /// <summary>
        /// does a convergence study for the l2 projection for some Exact solution
        /// </summary>
        /// <param name="LSDegree">level Set degree tested</param>
        /// <param name="ExactSol">an Exact solution for ST Burgers equation</param>
        /// <param name="ShockSpeed">a function that maps the graph of the discontinuity </param>
        /// <param name="nPol">polynomial Degree of highest degree solution</param>
        /// <param name="nGrid">number of reinements</param>
        /// <param name="is_nf_smth">ture if uwind smooth numerical flux is smooth</param>
        /// <param name="NonLinQuadDegree">quadrature degree</param>
        /// <returns></returns>
        public static MultidimensionalArray DoConverGenceStudy(int LSDegree, Func<double[],string, double> ExactSol, Func<double[], double> ShockSpeed, int nPol=6,int nGrid=6, bool is_nf_smth=true, int NonLinQuadDegree=6){

            //create the spatial Operator
            var XSpatialOperator = new XSpatialOperatorMk2(new string[] { "c" }, null, new string[] { "c" }, (int[] A, int[] B, int[] C) => NonLinQuadDegree, new string[] { "L","R"});

            //create BndValue Map from EcatSol
            Func<double[], double> BndVals = delegate (double[] x) {
                if(0 > ShockSpeed(x)) {
                    return ExactSol(x, "L");
                } else {
                    return ExactSol(x, "R");
                }
            };
            //add equations/fluxes
            XSpatialOperator.EquationComponents["c"].Add(new BUIDT.Fluxes.STBurgersUpwindFlux("L", is_nf_smth, 10, BndVals));
            XSpatialOperator.EquationComponents["c"].Add(new BUIDT.Fluxes.STBurgersUpwindFlux("R", is_nf_smth, 10, BndVals));
            XSpatialOperator.EquationComponents["c"].Add(new BurgersUpwindFlux_Interface( is_nf_smth, 10));
            XSpatialOperator.Commit();

            //create some grids
            GridCommons[] Grids = new GridCommons[nGrid];
            for(int i = 1; i < nGrid+1; i++) {
                int res = (int) Math.Pow(2, i);
                double[] xNodes = GenericBlas.Linspace(-0.2, 1.0, res + 1);
                double[] tNodes = GenericBlas.Linspace(0.0, 1.2, res + 1);
                Grids[i-1] = Grid2D.Cartesian2DGrid(xNodes, tNodes);
            }
            //Choose Polynomial Degrees
            int[] pDegrees = new int[nPol+1];
            for(int i = 0; i < nPol + 1; i++) { 
                pDegrees[i] = i;
            }
            var results = MultidimensionalArray.Create(Grids.Length,pDegrees.Length );
            for(int iGrid =0; iGrid <Grids.Length;iGrid++) {
                var grid = Grids[iGrid];    
                for(int iDeg = 0; iDeg < pDegrees.Length; iDeg++) {
                    var p=pDegrees[iDeg];
                    Basis basis = new Basis(grid.GridData, LSDegree);
                    LevelSet levelSet = new LevelSet(basis, "LevelSet");
                    levelSet.ProjectField(x => ShockSpeed(x));

                    LevelSetTracker LsTrK = new LevelSetTracker(grid.GridData, XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "L", "R" }, new LevelSet[] { levelSet });
                    LsTrK.UpdateTracker(0.0);
                    XDGBasis pBasis = new XDGBasis(LsTrK, p);
                    XDGField pField= new XDGField(pBasis,"c");
                    XDGField pField2 = new XDGField(pBasis, "c_fullProj");
                    XDGField Residual = new XDGField(pBasis,"c_res");
                    pField.GetSpeciesShadowField("L").ProjectField(x => ExactSol(x, "L"));
                    pField.GetSpeciesShadowField("R").ProjectField(x => ExactSol(x, "R"));
                    pField2.ProjectField(BndVals);
                    var p_eval = XSpatialOperator.GetEvaluatorEx(new DGField[] { pField }, null, pField.Mapping);
                    p_eval.Evaluate(1.0, 0.0, Residual.CoordinateVector);
                    results[iGrid, iDeg] = Residual.L2NormAllSpecies();

                    //var tp = new Tecplot(grid.iGridData, 3);
                    //tp.PlotFields("AccShockXDGConvergenceStudy_" + iDeg + iGrid, 0.0, new DGField[] { pField, pField2,levelSet, Residual });
                }
            }
            results.SaveToTextFile("AccShockXDGConvergenceStudyResults.txt");
            return results;

        }
        
    }

}
