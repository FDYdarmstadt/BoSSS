using ApplicationWithIDT;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace BUIDT
{
    public class BUIDTHardCodedControl
    {
        /// <summary>
        /// a base control setting some values 
        /// TODO: Move as many as possible to BUIDTControl / IDTControl
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="ImmediatePlotPeriod"></param>
        /// <returns></returns>
        public static BUIDTControl BaseControl(string dbPath, int MaxIterations, int dgDegree, int numOfCellsX,
        int numOfCellsY, int ImmediatePlotPeriod = -1)
        {
            var c = new BUIDTControl();
            c.DbPath = dbPath;
            c.savetodb = c.DbPath != null;
            c.NoOfTimesteps = MaxIterations;
            
            c.ImmediatePlotPeriod = ImmediatePlotPeriod;


            // Adaptive Regularization
            c.Gamma_Start = 1;
            c.Gamma_Min = 1e-06;
            c.Gamma_Max = 1;
            c.tauGamma = 1.5;
            c.L = 1;
            c.sigma_1 = 1e-2;
            c.sigma_2 = 1e-1;
            c.Kappa_Xi = 10;

            //Globalization Parameters
            c.Alpha_Start = 1;
            c.Alpha_Min = 1e-08;
            c.Mu_Start = 1;

            // ### Quadrature ###
            c.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            c.SolDegree = dgDegree;
            c.LsOne_NegSpecies = "L";
            c.LsOne_PosSpecies = "R";


            // ### Grid ###
            double xMin = 0;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 1.0;


            c.GridFunc = delegate
            {
                double[] xNodes;
                //if(levelSetPos == 0.5) {
                //    xNodes = new double[] { 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0 };
                //} else {
                xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };

            c.IsTwoLevelSetRun = false;
            c.LevelSetDegree = 1;
            c.NonlinearQuadratureDegree = 10;

            c.OptiLevelSet_ParamNames = new List<string>();
            c.OptiLevelSet_ParamValues = new List<double>();
            c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>();

            c.OptiLevelSet_ParamNames.Add("L");
            c.OptiLevelSet_ParamValues.Add(-0.2);
            c.OptiLevelSet_Param_Functions.Add((x, a) => a);

            c.OptiLevelSet_ParamNames.Add("R");
            c.OptiLevelSet_ParamValues.Add(1);
            c.OptiLevelSet_Param_Functions.Add((x, b) => x[0] * b);

            c.OptiLevelSet_ParamNames.Add("c");
            c.OptiLevelSet_ParamValues.Add(-0.6);
            c.OptiLevelSet_Param_Functions.Add((x, c) => x[1] * c);

            c.DirichletBoundaryMap = x => x[0] >= 0.25 ? 0.25 : 0.75; 

            c.SpeciesToEvaluate = new string[] { "L", "R" };

            // ### Project and sessions name ###
            c.ProjectName = "BurgersIDT";
            c.SessionName = c.ProjectName;

            return c;

        }
        /// <summary>
        /// This control corresponds to a straight shock, presented in Eccomas Paper
        /// - initial value is piecewise constant
        /// - initially a curved level Set is initialized
        /// TODO: Add section from paper, after publication
        public static BUIDTControl StraightShockCurvedStart_Eccomas22(string dbPath=null, int MaxIterations=50, int dgDegree=0, int numOfCellsX=10,
        int numOfCellsY=10, Linearization linearization=Linearization.FD, int lsDegree = 2,  int ImmediatePlotPeriod = -1, double agg = 0.4, GetLevelSet getLevelSet=GetLevelSet.FromFunction, OptiLevelSetType optiLevelSetType=OptiLevelSetType.SplineLevelSet)
        {
            var c = BUIDTHardCodedControl.BaseControl(
                    dbPath: dbPath,
                    MaxIterations: MaxIterations,
                    dgDegree: dgDegree,
                    numOfCellsX: numOfCellsX,
                    numOfCellsY: numOfCellsY,
                    ImmediatePlotPeriod: ImmediatePlotPeriod
                    );
            c.Linearization = linearization;
            c.LevelSetDegree = lsDegree;
            c.OptiLevelSetDegree = lsDegree;
            c.OptiLevelSetType= optiLevelSetType;
            c.GetLevelSet= getLevelSet;
            switch (optiLevelSetType)
            {
                case OptiLevelSetType.SplineLevelSet:
                    c.InitialShockPostion = y => 0.4 + y[1] * 0.6 - 0.2 * y[1] * y[1];
                    break;
                case OptiLevelSetType.GlobalLevelSet:
                    switch (getLevelSet)
                    {
                        case GetLevelSet.FromFunction:
                            c.InitialShockPostion = x => x[0] - (0.4 + x[1] * 0.6 - 0.2 * x[1] * x[1]);
                            break;
                        case GetLevelSet.FromParams:
                            c.OptiLevelSet_ParamNames = new List<string>();
                            c.OptiLevelSet_ParamValues = new List<double>();
                            c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>();

                            c.OptiLevelSet_ParamNames.Add("L");
                            c.OptiLevelSet_ParamValues.Add(1.0);
                            c.OptiLevelSet_Param_Functions.Add((x, a) => x[0] * a);

                            c.OptiLevelSet_ParamNames.Add("R");
                            c.OptiLevelSet_ParamValues.Add(-0.4);
                            c.OptiLevelSet_Param_Functions.Add((x, b) => b);

                            c.OptiLevelSet_ParamNames.Add("c");
                            c.OptiLevelSet_ParamValues.Add(-0.6);
                            c.OptiLevelSet_Param_Functions.Add((x, c) => x[1] * c);

                            c.OptiLevelSet_ParamNames.Add("d");
                            c.OptiLevelSet_ParamValues.Add(0.2);
                            c.OptiLevelSet_Param_Functions.Add((x, d) => x[1] * x[1] * d);
                            break;
                        default:
                            throw new Exception($"{getLevelSet.ToString()} not supported for this control");
                    }
                    break;
                default:
                    c.InitialShockPostion = x => x[0] - (0.4 + x[1] * 0.6 - 0.2 * x[1] * x[1]);
                break;
            }
            c.InitialValueFunctionsPerSpecies.Clear();
            c.InitialValueFunctionsPerSpecies.Add("L", x => x[0] < 0.5 * x[1] ? 0.75 : 0.25);
            c.InitialValueFunctionsPerSpecies.Add("R", x => x[0] < 0.5 * x[1] ? 0.75 : 0.25);
            c.UseP0ProjectionAsInitialGuess = true;
            c.AgglomerationThreshold = agg;
            return c;
        }
        /// <summary>
        /// This control corresponds to a straight shock and a PWC solution
        /// - initial value is piecewise constant
        /// - as a first guess a curved Level Set is initialized 
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="linearization"></param>
        /// <param name="ImmediatePlotPeriod"></param>
        /// <param name="agg"></param>
        /// <param name="OptiNumOfCellsX"></param>
        /// <param name="OptiNumOfCellsY"></param>
        /// <param name="optiLevelSetType"></param>
        /// <returns></returns>
        public static BUIDTControl StraightShock(string dbPath, int MaxIterations, int dgDegree, int numOfCellsX,
        int numOfCellsY, Linearization linearization, int ImmediatePlotPeriod = -1, double agg = 0.0, int OptiNumOfCellsX=1, int OptiNumOfCellsY=1, OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet) {
            var c = BUIDTHardCodedControl.BaseControl(
                    dbPath: dbPath,
                    MaxIterations: MaxIterations,
                    dgDegree: dgDegree,
                    numOfCellsX: numOfCellsX,
                    numOfCellsY: numOfCellsY,                    
                    ImmediatePlotPeriod: ImmediatePlotPeriod
                    );

            // ### Grid ###
            double xMin = 0;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 1.0;

            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = 1;
            c.LevelSetDegree = 1;
            c.LevelSetGridFunc = delegate {
                double[] xNodes;
                //if(levelSetPos == 0.5) {
                //    xNodes = new double[] { 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0 };
                //} else {
                xNodes = GenericBlas.Linspace(xMin, xMax, OptiNumOfCellsX + 1);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, OptiNumOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };

            c.InitialShockPostion = x => x[0] - 0.44;
            c.UseP0ProjectionAsInitialGuess = true;
            c.NonlinearQuadratureDegree = 10;
            c.SuperSampling = 4;
            c.AgglomerationThreshold = agg;
            c.Linearization = linearization;
            return c;
        }
        /// <summary>
        /// This control corresponds to the Accelerating Shock Case, (taken from Huang et al.-  A robust, high-order implicit shock tracking method for simulation of complex, high-speed flows) and presented in 2023 paper
        /// - Shock is curved and solution is non-polynomial
        /// - Staggered solver is used to obtain solution
        /// 
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="linearization"></param>
        /// <param name="applyReInit"></param>
        /// <param name="ImmediatePlotPeriod"></param>
        /// <param name="agg"></param>
        /// <param name="OptiNumOfCellsX"></param>
        /// <param name="OptiNumOfCellsY"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="LSDregree"></param>
        /// <param name="getLevelSet"></param>
        /// <param name="solverRunType"></param>
        /// <param name="staggeredTS"></param>
        /// <returns></returns>
        public static BUIDTControl AccShock(string dbPath, int MaxIterations, int dgDegree, int numOfCellsX,
        int numOfCellsY, Linearization linearization, bool applyReInit=false, int ImmediatePlotPeriod = -1, double agg = 0.0, int OptiNumOfCellsX = 1, int OptiNumOfCellsY = 1, 
        OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet, int LSDregree=3,GetLevelSet getLevelSet = GetLevelSet.FromFunction,
            SolverRunType solverRunType = SolverRunType.Standard, int[] staggeredTS=null) {
            var c = BUIDTHardCodedControl.BaseControl(
                    dbPath: dbPath,
                    MaxIterations: MaxIterations,
                    dgDegree: dgDegree,
                    numOfCellsX: numOfCellsX,
                    numOfCellsY: numOfCellsY,
                    ImmediatePlotPeriod: ImmediatePlotPeriod
                    );

            // ### Grid ###
            double xMin = -0.2;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 1.2;

            c.GridFunc = delegate {
                double[] xNodes;
                //if(levelSetPos == 0.5) {
                //    xNodes = new double[] { 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0 };
                //} else {
                xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };

            double s_alpha = 10;
            c.s_alpha = s_alpha;
            double mu1 = 4.0;
            double mu2 = 3.0;
            double x0 = 0.2501;
            c.fphiType = FphiType.None; // does not work with any fphi, bc. shock is "too" curved
            c.GetLevelSet = getLevelSet;
            double SmoothedHeaviSide(double x) {
                return 1 / (1 + Math.Exp(-2 * s_alpha * x));
            }
            double ShockSpeed(double t) {
                return (mu1 / mu2 + 1) * (1 - Math.Sqrt(1 + mu2 * t)) + mu1 * t;
            }
            double ExactSolUnsmooth(double[] X) {
                if(X[0] < ShockSpeed(X[1])) {
                    return mu1;
                } else {
                    return mu2 * (X[0] - 1) / (mu2 * X[1] + 1);
                }
            }
            c.InitialValueFunctionsPerSpecies = new Dictionary<string, Func<double[], double>>();
            c.InitialValueFunctionsPerSpecies.Add("L", delegate ( double[] X) {
                return ExactSolUnsmooth(X);
                //return mu1;
            });
            c.InitialValueFunctionsPerSpecies.Add("R", delegate (double[] X) {
                return ExactSolUnsmooth(X);
                //return mu2 * (X[0] - 1)  / (mu2 * X[1] + 1);
            });
            c.UseP0ProjectionAsInitialGuess = false;
            //double ExactSolSmooth(double[] X) {
            //    return mu1 * SmoothedHeaviSide(ShockSpeed(X[1]) - X[0]) + mu2 * (X[0] - 1) * (1 - SmoothedHeaviSide(X[0] - ShockSpeed(X[1]))) / (mu2 * X[1] + 1);
            //}

            c.MeritFunctionType = MeritFunctionType.FullyL2Merit;
            c.DirichletBoundaryMap = x => ExactSolUnsmooth(x);
            c.ApplyReiInit = applyReInit;
            c.OptiLevelSetType = optiLevelSetType;
            c.is_nf_smth = false;
            c.OptiLevelSetDegree = LSDregree ;
            c.LevelSetDegree = LSDregree;
            c.LevelSetGridFunc = delegate {
                double[] xNodes;
                //if(levelSetPos == 0.5) {
                //    xNodes = new double[] { 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0 };
                //} else {
                xNodes = GenericBlas.Linspace(xMin, xMax, OptiNumOfCellsX + 1);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, OptiNumOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };

            //Exact LevelSetPosition
            if(optiLevelSetType == OptiLevelSetType.SplineLevelSet) {
                //c.InitialShockPostion = x =>ShockSpeed(x[1]);
                c.InitialShockPostion = x =>  x[1]*1.2;
                //Incase of Point Reconstruction
                c.SLSPointPath = @"../../../SplinePointsCurvedBurgers.txt";
            } else {
                c.InitialShockPostion = x => x[0] - ShockSpeed(x[1]);
            }

            

#region     Initial Vlaue for Global LevelSet
            c.OptiLevelSet_ParamNames = new List<string>();
            c.OptiLevelSet_ParamValues = new List<double>();
            c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>();

            c.OptiLevelSet_ParamNames.Add("L");
            c.OptiLevelSet_ParamValues.Add(1.0);
            c.OptiLevelSet_Param_Functions.Add((x, a) => x[0] * a);

            c.OptiLevelSet_ParamNames.Add("R");
            c.OptiLevelSet_ParamValues.Add(-0.0);
            c.OptiLevelSet_Param_Functions.Add((x, b) => b);

            c.OptiLevelSet_ParamNames.Add("c");
            c.OptiLevelSet_ParamValues.Add(-0.6);
            c.OptiLevelSet_Param_Functions.Add((x, c) => x[1] * c);

            c.OptiLevelSet_ParamNames.Add("d");
            c.OptiLevelSet_ParamValues.Add(-0.6);
            c.OptiLevelSet_Param_Functions.Add((x, d) => x[1] * x[1] * d);

            c.OptiLevelSet_ParamNames.Add("e");
            c.OptiLevelSet_ParamValues.Add(-0.6);
            c.OptiLevelSet_Param_Functions.Add((x, d) => x[1] * x[1] * x[1] * d);

            c.OptiLevelSet_ParamNames.Add("f");
            c.OptiLevelSet_ParamValues.Add(-0.05);
            c.OptiLevelSet_Param_Functions.Add((x, d) => x[1] * x[1] * x[1] * x[1] * d);
#endregion
            c.solRunType = solverRunType;
            if(staggeredTS != null)
                c.staggeredTimeSteps = staggeredTS;
            c.Gamma_Max = 1;
            c.MinPIter = new int[] { 20, 15, 10, 10, 10 };
            c.SuperSampling = 4;
            c.AgglomerationThreshold = agg;
            c.Linearization = linearization;
            c.Queries.Add("L2ErrorExactSol", BUIDTQueries.L2Error("c", referenceSolution: (X, t) => ExactSolUnsmooth(X)));
            c.Queries.Add("L2ErrorLevelSet", BUIDTQueries.L2Error1D(referenceSolution: y => ShockSpeed(y)));
            return c;
        }
      
    }
}
