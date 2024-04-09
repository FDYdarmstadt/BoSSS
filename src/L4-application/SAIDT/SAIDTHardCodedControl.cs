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
using System.Data;

namespace SAIDT {
    /// <summary>
    /// This class defines static functions creating Control Objects for the SAIDT solver.
    /// </summary>
    public class SAIDTHardCodedControl {

        /// <summary>
        /// Base Control used as a starting point for other Controls
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="ImmediatePlotPeriod"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="ApplyReInit"></param>
        /// <param name="meritFunctionType"></param>
        /// <returns></returns>
        public static SAIDTControl BaseControl(string dbPath, int MaxIterations, int dgDegree, int numOfCellsX,
        int numOfCellsY, int ImmediatePlotPeriod = -1, OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet, bool ApplyReInit = true, MeritFunctionType meritFunctionType = MeritFunctionType.L1Merit) {


            var c = new SAIDTControl();
            c.DbPath = dbPath;
            c.savetodb = c.DbPath != null;
            c.NoOfTimesteps = MaxIterations;
            c.ImmediatePlotPeriod = ImmediatePlotPeriod;

            #region Optimization Parameters
            c.MeritFunctionType = meritFunctionType;
            c.OptiLevelSetType = optiLevelSetType;
            c.fphiType = FphiType.None;
            // Reinitialization
            c.ApplyReiInit = ApplyReInit;
            #endregion
            #region  Grid
            double xMin = 0;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 1.0;
            c.GridFunc = delegate {
                double[] xNodes;
                xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };
            #endregion
            #region Problem Definition 
            c.DirichletBoundaryMap = (x => x[0] <= 0.25 ? 1 : 0);
            c.FlowFunc = (x => 0.5);
            c.SolDegree = dgDegree;
            #endregion
            #region Initial Guess for state 
            c.InitialValueFunctionsPerSpecies.Clear();
            c.InitialValueFunctionsPerSpecies.Add("L", x => x[0] < 0.5 * x[1] ?1 : 0);
            c.InitialValueFunctionsPerSpecies.Add("R", x => x[0] < 0.5 * x[1] ? 1 : 0);
            #endregion
            #region  LevelSet 
            // names for sub-domains
            c.LsOne_NegSpecies = "L";
            c.LsOne_PosSpecies = "R";
            c.SpeciesToEvaluate = new string[] { "L", "R" };
            c.IsTwoLevelSetRun = false;
            //polynomial degree
            c.LevelSetDegree = 1;
            c.OptiLevelSetDegree = 1;
            //set up the initial guess for the Global level Set (only used if chosen)
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
            #endregion
            #region  Project and sessions name ###
            c.ProjectName = "ScalarAdvection";
            c.SessionName = c.ProjectName;
            #endregion
            return c;

        }
        /// <summary>
        /// Curved Shock example featured in ECCOMAS22 paper.
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="ImmediatePlotPeriod"></param>
        /// <param name="agg"></param>
        /// <param name="optiLevelSetType"></param>
        /// <returns></returns>
        public static SAIDTControl CurvedShock_Eccomas22(string dbPath = null, int MaxIterations = 50, int dgDegree = 0, int numOfCellsX = 10,
        int numOfCellsY = 10, int ImmediatePlotPeriod = -1, double agg = 0.2, OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet) {
            var c = SAIDTHardCodedControl.BaseControl(
                    dbPath: dbPath,
                    MaxIterations: MaxIterations,
                    dgDegree: dgDegree,
                    numOfCellsX: numOfCellsX,
                    numOfCellsY: numOfCellsY, ImmediatePlotPeriod: ImmediatePlotPeriod,
                    optiLevelSetType: optiLevelSetType
                    );
            // Problem Def
            c.FlowFunc = (t => 3 * t * t - 3 * t + 0.5);
            c.fphiType = FphiType.None;
            #region  LevelSet 
            //set up the initial guess for the Global level Set (if chosen)
            c.LevelSetDegree = 3;
            c.OptiLevelSetDegree = 3;
            c.OptiLevelSet_ParamNames = new List<string>();
            c.OptiLevelSet_ParamValues = new List<double>();
            c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>();
            c.OptiLevelSet_ParamNames.Add("L");
            c.OptiLevelSet_ParamValues.Add(1.0);
            c.OptiLevelSet_Param_Functions.Add((x, a) => x[0] * a);

            c.OptiLevelSet_ParamNames.Add("R");
            c.OptiLevelSet_ParamValues.Add(-0.1);
            c.OptiLevelSet_Param_Functions.Add((x, b) => b);

            c.OptiLevelSet_ParamNames.Add("c");
            c.OptiLevelSet_ParamValues.Add(-0.7);
            c.OptiLevelSet_Param_Functions.Add((x, c) => x[1] * c);

            c.OptiLevelSet_ParamNames.Add("d");
            c.OptiLevelSet_ParamValues.Add(1.0);
            c.OptiLevelSet_Param_Functions.Add((x, d) => x[1] * x[1] * d);

            c.OptiLevelSet_ParamNames.Add("e");
            c.OptiLevelSet_ParamValues.Add(-0.7);
            c.OptiLevelSet_Param_Functions.Add((x, e) => x[1] * x[1] * x[1] * e);

            //set up the initial guess for the Spline level Set (if chosen)
            if(c.OptiLevelSetType == OptiLevelSetType.SplineLevelSet) {
                c.LevelSetTwoInitialValue = x => 0.1 + 0.7 * x[1] - x[1] * x[1] + 0.7 * x[1] * x[1] * x[1];
            } else {
                c.LevelSetTwoInitialValue = x => x[0] - 0.1 - 0.7 * x[1] + x[1] * x[1] - 0.7 * x[1] * x[1] * x[1];
            }
            #endregion
            c.GetInitialValue = GetInitialValue.OneFullNewtonStep;
            Func<double, double> exactShock = t => t*t*t-3.0/2.0*t*t+0.5*t+0.25;
            c.InitialValueFunctionsPerSpecies.Clear();
            c.InitialValueFunctionsPerSpecies.Add("L", x => exactShock(x[1])-x[0]  > 0 ? 1 : 0);
            c.InitialValueFunctionsPerSpecies.Add("R", x => exactShock(x[1]) - x[0] > 0 ? 1 : 0);
            c.UseP0ProjectionAsInitialValue = true;
            c.SuperSampling = 4;
            c.AgglomerationThreshold = agg;
            return c;
        }
        /// <summary>
        /// Example with straight shock
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="OptiNumOfCellsX"></param>
        /// <param name="OptiNumOfCellsY"></param>
        /// <param name="ImmediatePlotPeriod"></param>
        /// <param name="agg"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="LSDegree"></param>
        /// <param name="withReInit"></param>
        /// <param name="meritFunctionType"></param>
        /// <param name="isFarConfig"></param>
        /// <returns></returns>
        public static SAIDTControl StraightShock(string dbPath = null, int MaxIterations = 50, int dgDegree = 0, int numOfCellsX = 10,
        int numOfCellsY = 10,
        int OptiNumOfCellsX = 1,
        int OptiNumOfCellsY = 1, int ImmediatePlotPeriod = -1, double agg = 0.0, OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet, int LSDegree = 1, bool withReInit = true, MeritFunctionType meritFunctionType=MeritFunctionType.L1Merit,bool isFarConfig =false) {
            var c = SAIDTHardCodedControl.BaseControl(
                    dbPath: dbPath,
                    MaxIterations: MaxIterations,
                    dgDegree: dgDegree,
                    numOfCellsX: numOfCellsX,
                    numOfCellsY: numOfCellsY, ImmediatePlotPeriod: ImmediatePlotPeriod,
                    optiLevelSetType: optiLevelSetType,
                    ApplyReInit: withReInit,
                    meritFunctionType: meritFunctionType
                    ) ;
            // ### Grid ###
            double xMin = 0;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 1.0;
            

            c.LevelSetGridFunc = delegate {
                double[] xNodes;
                xNodes = GenericBlas.Linspace(xMin, xMax, OptiNumOfCellsX + 1);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, OptiNumOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };
            c.isFarConfig = isFarConfig;
            if(isFarConfig) {
                
                if(optiLevelSetType == OptiLevelSetType.SplineLevelSet) {
                    c.LevelSetTwoInitialValue = x => 0.2 + x[1] * 0.2 * x[1] * x[1] * 0.5;
                } else {
                    c.LevelSetTwoInitialValue = x => x[0] -(0.2 + x[1] * 0.2 * x[1] * x[1] * 0.5);
                }
                
            } else {
                if(optiLevelSetType == OptiLevelSetType.SplineLevelSet) {
                    c.LevelSetTwoInitialValue = x => 0.2 + x[1] * 0.6;
                } else {
                    c.LevelSetTwoInitialValue = x => x[0] - 0.2 - x[1] * 0.6;
                }
            }
            c.UseP0ProjectionAsInitialValue = true;
            c.NonlinearQuadratureDegree = 10;

            c.SuperSampling = 4;
            
            
            c.AgglomerationThreshold = agg;
            c.OptiLevelSetDegree = LSDegree;
            c.LevelSetDegree = LSDegree;
            if(optiLevelSetType == OptiLevelSetType.SpecFemField) {
                c.LevelSetDegree = LSDegree*2;
            }
            c.SessionName = string.Format("SAIDT-J{0}_p{1}_q{2}_{3}_{4}_isFarConfig{5}", numOfCellsX*numOfCellsY, dgDegree, LSDegree, optiLevelSetType, meritFunctionType,isFarConfig);
            return c;
        }
        /// <summary>
        /// Example with sinusoidal shock
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="OptiNumOfCellsX"></param>
        /// <param name="OptiNumOfCellsY"></param>
        /// <param name="ImmediatePlotPeriod"></param>
        /// <param name="agg"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="LSDegree"></param>
        /// <returns></returns>
        public static SAIDTControl SinusShock(string dbPath = null, int MaxIterations = 50, int dgDegree = 0,
        int numOfCellsX = 10,
        int numOfCellsY = 10,
        int OptiNumOfCellsX = 1,
        int OptiNumOfCellsY = 1,
        int ImmediatePlotPeriod = -1,
        double agg = 0.0,
        OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet,
        int LSDegree = 2) {
            var c = SAIDTHardCodedControl.BaseControl(
                    dbPath: dbPath,
                    MaxIterations: MaxIterations,
                    dgDegree: dgDegree,
                    numOfCellsX: numOfCellsX,
                    numOfCellsY: numOfCellsY, ImmediatePlotPeriod: ImmediatePlotPeriod,
                    optiLevelSetType: optiLevelSetType
                    );

            c.FlowFunc = (x => 0.5 * Math.Cos(2 * Math.PI * x));
            // ### Grid ###
            double xMin = 0;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 1.0;

            c.Alpha_Min = 0.1;
            c.OptiLevelSetDegree = LSDegree;
            c.LevelSetDegree = LSDegree;
            c.LevelSetGridFunc = delegate {
                double[] xNodes;
                xNodes = GenericBlas.Linspace(xMin, xMax, OptiNumOfCellsX + 1);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, OptiNumOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };

            if(c.OptiLevelSetType == OptiLevelSetType.SplineLevelSet) {
                c.LevelSetTwoInitialValue = x => 0.5 / (2 * Math.PI) * Math.Sin(2 * Math.PI * x[1]) +0.24 +0.1 * x[1];
            } else {
                c.LevelSetTwoInitialValue = x => x[0] - 0.5 / (2 * Math.PI) * Math.Sin(2 * Math.PI * x[1]) - 0.24 - 0.1 * x[1];
            }
            c.UseP0ProjectionAsInitialValue = true;
            c.NonlinearQuadratureDegree = LSDegree * 4 > 20 ? 20 : LSDegree * 4;
            c.SuperSampling = 4;
            c.AgglomerationThreshold = agg;
            return c;
        }
    }
}
