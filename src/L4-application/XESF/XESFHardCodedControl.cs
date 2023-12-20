using ApplicationWithIDT;
using ApplicationWithIDT.OptiLevelSets;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using XESF.Variables;
using XESF.Fluxes;


namespace XESF {
    public class XESFHardCodedControl {

        /// <summary>
        /// Base Control for supersonic flow over an inclined plane (aka. WedgeFlow) discretized using 2 level sets (TwoLs).
        /// - one LS is used to represent the immersed boundary (i.e. the inclined Plane)
        /// - the second level Set represents the straight shock that is obtained
        /// -  as long as the shock does not detach or vanish the solution is piecewise constant
        /// </summary>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="wedge_angle"></param>
        /// <param name="initialAngle_shockLS"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsDegree"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="interfaceFluxLS2"></param>
        /// <param name="FluxVersion"></param>
        /// <param name="shocksetup"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="optiLSDegree"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        public static XESFControl XDGWedgeFlow_TwoLs_Base(int MaxIterations=20, int dgDegree=0, int numOfCellsX=10,
        int numOfCellsY=10, double wedge_angle = 10,double initialAngle_shockLS = 32, int PlotInterval=-1,
        string dbPath = null, int lsDegree = 1, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.OptimizedHLLC,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var,
        ConvectiveInterfaceFluxes interfaceFluxLS2 = ConvectiveInterfaceFluxes.GodunovInterface,
        XESF.Fluxes.FluxVersion FluxVersion = Fluxes.FluxVersion.Optimized, GetLevelSet shocksetup = GetLevelSet.FromParams, OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet, int optiLSDegree = 1,
        GetInitialValue initialValue = GetInitialValue.FromFunctionPerSpecies, bool iVFromShockRelations = false, double agg = 0.2, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch, MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit) {


            XESFControl c = new XESFControl();
            #region DBStuff and Plotting
            c.DbPath = dbPath;
            c.savetodb = true;
            if(dbPath == null)
                c.savetodb = false;
            c.NoOfTimesteps = MaxIterations;
            c.ImmediatePlotPeriod = PlotInterval;
            if(dgDegree == 0 && numOfCellsX < 10 && numOfCellsY < 10) {
                c.SuperSampling = 4;
            } else {
                c.SuperSampling = 3;
            }

            #endregion

            #region SQP constants/parameters
            // ### Time config ###
            c.NoOfTimesteps = MaxIterations;

            //Globalization Parameters
            c.GlobalizationStrategy = globalization;
            c.MeritFunctionType = meritFunctionType;
            c.MinPIter = new int[] { 25, 25, 25, 0, 0 };

            #endregion

            #region LevelSetStuff

            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            c.IsTwoLevelSetRun = true;

            c.LevelSetDegree = 1;
            c.AddVariable(XESFVariables.LevelSet, 1); //this is the levelSet used for the Immersed Boundary (the Wedge)
            c.AddVariable(XESFVariables.LevelSetTwo, lsDegree); // this is the levelSet used for the shock

            c.SpeciesTable[0, 0] = "X"; // 'Forbidden' species
            c.SpeciesTable[0, 1] = "V"; // Void area (inside the blunt body)
            c.SpeciesTable[1, 0] = "L"; // Pre-shock region (on the left of the shock)
            c.SpeciesTable[1, 1] = "R"; // Post-shock region (on the right of the shock)
            c.LsOne_NegSpecies = "V";
            c.LsOne_PosSpecies = "R";
            c.LsTwo_NegSpecies = "L";
            c.LsTwo_PosSpecies = "R";
            c.SpeciesToEvaluate = new string[] { "L", "R" };

            c.LevelSetTwoDegree = lsDegree;
            c.OptiLevelSetType = optiLevelSetType;

            //shockLevelSet_Db = @"C:\BoSSS-experimental\internal\src\private-mag\XDGShock\Tests\bosss_db_levelSets.zip";

            #endregion 

            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.ConvectiveInterfaceFlux_LsTwo = interfaceFluxLS2;
            c.FluxVersion = FluxVersion;
            c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;
            #endregion

            #region Grid
            double xMin = 0;
            double xMax = 1.5;
            double yMin = 0;
            double yMax = 1.0;
            double wedge_angle_radial = wedge_angle * Math.PI / 180.0; // Angle of Wedge

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                if(numOfCellsX > 1 && xNodes[1] >= 0.5) {
                    xNodes[1] = 0.4;
                }
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");
                //grid.EdgeTagNames.Add(4, "SubsonicOutlet");


                grid.DefineEdgeTags(delegate (double[] X) {
                    if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                                                           ////if((-X[0] + 0.5 + (X[1] / Math.Tan(wedge_angle)))>=0){ 
                                                           //if(X[1] >= Math.Tan(39.3139318 * Math.PI / 180.0)) { //-0.01 doesn't work
                                                           //    return 2;
                                                           //} else { // (part of void area)
                        return 2;
                        //return 3; // Wedge Part of Outlet

                    } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                        return 1;
                    } else if(Math.Abs(X[1] - yMax) < 1e-14) { // top boundary
                        return 2;
                        //if(X[0] < 0.5) {
                        //    return 1;
                        //} else {
                        //    return 2;
                        //}

                    } else { //if(Math.Abs(X[1] - yMin) < 1e-14) { //bottom boundary
                        return 3;
                        //if(X[0] <= 1) {
                        //    return 2;
                        //} else {
                        //    return 2;
                        //}
                    }
                });
                return grid;
            };
            #endregion

            #region Initial Position of LevelSets
            double Mach = 2.0;
            double shock_angle_exact = 39.3139318;
            double shock_angle_radial;
            if(iVFromShockRelations) {
                shock_angle_radial = initialAngle_shockLS * Math.PI / 180.0; // this angle is the shock angle assumed in the beginning and used to calculate the IV
            } else {
                shock_angle_radial = shock_angle_exact * Math.PI / 180.0; // when we want to project the initial solution we need to use the exact angle in the calculations
            }
            double LevelSet1Prime = Math.Tan(wedge_angle_radial);
            double LevelSet2Prime = Math.Tan(initialAngle_shockLS * Math.PI / 180.0);

            
            // ### Wedge Level set function ###
            c.LevelSetOneInitialValue = delegate (double[] X) { return -X[0] + 0.5 + (X[1] / 0.17632698070846498); };

            //// Shock level set
            //c.LevelSetOneInitialValue = delegate (double[] X) { return X[0] - 0.5 - (X[1] / LevelSet2Prime); };
            c.GetLevelSet = shocksetup;
            switch(shocksetup) {
                case GetLevelSet.FromParams:
                if(optiLSDegree == 0) {
                    c.OptiLevelSet_ParamNames = new List<string> { "(x-0.5) ", "t" };
                    c.OptiLevelSet_ParamValues = new List<double> { 1.0, -1.0 / LevelSet2Prime };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a * (x[0] - 0.5), (x, a) => x[1] * a };
                } else if(optiLSDegree == 1) {
                    c.OptiLevelSet_ParamNames = new List<string> { " ", "x", "t" };
                    c.OptiLevelSet_ParamValues = new List<double> { -0.5, 1.0, -1.0 / LevelSet2Prime };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a, (x, a) => a * x[0], (x, a) => x[1] * a };
                } else if(optiLSDegree == 2) {
                    c.OptiLevelSet_ParamNames = new List<string> { " ", "x", "t", "t^2" };
                    c.OptiLevelSet_ParamValues = new List<double> { -0.5, 1.0, -1.0 / LevelSet2Prime, 0.2 };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a, (x, a) => a * x[0], (x, a) => x[1] * a, (x, a) => x[1] * x[1] * a };

                } else if((optiLSDegree == 3)) {
                    c.OptiLevelSet_ParamNames = new List<string> { " ", "x", "t", "t^2", "t^3" };
                    c.OptiLevelSet_ParamValues = new List<double> { -0.5, 1.0, -1.0 / LevelSet2Prime, 0.2, -0.2 };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a, (x, a) => a * x[0], (x, a) => x[1] * a, (x, a) => x[1] * x[1] * a, (x, a) => x[1] * x[1] * x[1] * a };
                } else {
                    throw new NotSupportedException("optiLSDegree:" + optiLSDegree + " not supported!");
                }

                break;
                case GetLevelSet.FromFunction:

                GlobalOptiLevelSet tmp_LS = GlobalOptiLevelSet.CreateONBLevelSet(xMin, xMax, yMin, yMax, optiLSDegree);
                c.OptiLevelSet_ParamNames = new List<string>(tmp_LS.m_ParamNames);
                c.OptiLevelSet_ParamValues = new List<double>(tmp_LS.m_ParamValues);
                c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>(tmp_LS.m_phi);
                c.OptiLSIsOrthonormal = true;
                
                
                break;
                case GetLevelSet.FromReconstruction:
                case GetLevelSet.FromReconstructionFromPoints:
                GlobalOptiLevelSet tmp_LS2 = GlobalOptiLevelSet.CreateONBLevelSet(xMin, xMax, yMin, yMax, optiLSDegree);
                c.OptiLevelSet_ParamNames = new List<string>(tmp_LS2.m_ParamNames);
                c.OptiLevelSet_ParamValues = new List<double>(tmp_LS2.m_ParamValues);
                c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>(tmp_LS2.m_phi);
                c.OptiLSIsOrthonormal = true;
                break;
                default:
                throw new NotSupportedException("Not supported Shock Setup");

            }
            //Initial Level Set as function
            switch (optiLevelSetType)
            {
                case OptiLevelSetType.SplineLevelSet:
                    c.LevelSetTwoInitialValue = Y => 0.5 + Y[1] / LevelSet2Prime;
                    break;
                default:
                    c.LevelSetTwoInitialValue = X => X[0] - 0.5 - (X[1] / LevelSet2Prime);
                    break;
            }
            #endregion

            #region Initial Values for Solution + Boundary Values
            c.SolDegree = dgDegree;
            c.GetInitialValue = initialValue;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.ExactEnthalpy = 6.3;
            c.AddVariable(XESFVariables.Enthalpy_Error, dgDegree);
            c.AddVariable(XESFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESFVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(XESFVariables.Pressure, dgDegree);
            //c.AddVariable(XESFVariables.PressureStep, dgDegree); here we have some bug

            c.AddVariable(XESFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;


            double gamma = IdealGas.Air.HeatCapacityRatio;
            double speedOfSound = Math.Sqrt(gamma * pressure / density);

            double velocityX = Mach * speedOfSound;
            double velocityY = 0.0;

            double momentumX = density * velocityX;


            double innerEnergy = pressure / (gamma - 1);
            double kineticEnergy = 0.5 * density * (velocityX * velocityX + velocityY * velocityY);
            double totalEnergy = innerEnergy + kineticEnergy;

            // Boundary conditions in PRIMITIVE variables
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#L", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#R", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => pressure);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => pressure);


            double density_Right = (gamma + 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) / ((gamma - 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) + 2.0);
            double pressure_Right = 1.0 + (2.0 * gamma * (Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) - 1.0)) / (gamma + 1.0);

            // calculate Velocity via Transformation matrix
            MultidimensionalArray TransMat = MultidimensionalArray.Create(2, 2);
            TransMat[0, 0] = 1.0 / Math.Sqrt(1.0 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[0, 1] = Math.Tan(shock_angle_radial) / Math.Sqrt(1 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[1, 0] = Math.Tan(shock_angle_radial) / Math.Sqrt(1 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[1, 1] = -1.0 / Math.Sqrt(1.0 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));


            double[] velocity_left_normal = new double[2];
            velocity_left_normal[0] = velocityX * Math.Cos(shock_angle_radial);
            velocity_left_normal[1] = velocityX * Math.Sin(shock_angle_radial);


            double[] velocity_right_normal = new double[2];
            double[] velocity_right_cartesian = new double[2];


            velocity_right_normal[0] = velocity_left_normal[0];
            velocity_right_normal[1] = velocity_left_normal[1] / density_Right;

            TransMat.MatVecMul(1.0, velocity_right_normal, 0.0, velocity_right_cartesian);

            double velocityX_Right = velocity_right_cartesian[0];
            double velocityY_Right = velocity_right_cartesian[1];

            // calculate Velocity via Right-side Mach number

            //double Mach_Right = 1.6405;
            //double speedOfSoundRight = Math.Sqrt(gamma * pressure_Right / density_Right);
            //double velocity_norm = Mach_Right * speedOfSoundRight;
            //double velocityX_Right = velocity_norm;
            //double velocityY_Right = velocity_norm*Math.Tan(wedge_angle);

            double innerEnergy_Right = pressure_Right / (gamma - 1);
            double kineticEnergy_Right = 0.5 * density_Right * (velocityX_Right * velocityX_Right + velocityY_Right * velocityY_Right);
            double totalEnergy_Right = innerEnergy_Right + kineticEnergy_Right;

            double exact_sol(double[] X, int component) {
                //if(X[0] > 0.5 + (X[1] / LevelSet1Prime)) { //wedge -X[0] + 0.5 + (X[1] / LevelSet1Prime)
                //    return 0;
                if(X[0] >= 0.5 + (X[1] / Math.Tan(39.3139318 * Math.PI / 180.0))) {/*&& (X[0] <= 0.5 + (X[1] / LevelSet1Prime)))*/ //post Shock
                    if(component == 0) {
                        return density_Right;
                    } else if(component == 1) {
                        return velocityX_Right;
                    } else if(component == 2) {
                        return velocityY_Right;
                    } else if(component == 3) {
                        return pressure_Right;
                    }
                } else { //pre Shock
                    if(component == 0) {
                        return density;
                    } else if(component == 1) {
                        return velocityX;
                    } else if(component == 2) {
                        return 0;
                    } else if(component == 3) {
                        return pressure;
                    }
                }
                return 0;
            }

            if(iVFromShockRelations) {

                //// Initial conditions in PRIMITIVE variables -- Pre Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => density);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => velocityX);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => velocityY);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => pressure);

                //// Initial conditions in PRIMITIVE variables -- Post Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => velocityX_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => velocityY_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => pressure_Right);

                //// Initial conditions in Conservative variables -- Post Shock
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.xComponent + "#R", X => density_Right * velocityX_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.yComponent + "#R", X => density_Right * velocityY_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Energy + "#R", X => totalEnergy_Right);
            } else {
                //// Initial conditions in PRIMITIVE variables -- Pre Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => exact_sol(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => exact_sol(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => exact_sol(X, 2));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => exact_sol(X, 3));

                //// Initial conditions in PRIMITIVE variables -- Post Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => exact_sol(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => exact_sol(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => exact_sol(X, 2));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => exact_sol(X, 3));
            }
            #endregion

            #region Queries
            //var enthalpy_inflow = (gamma) / (gamma - 1) + 0.5 * velocityX * velocityX;
            //c.Queries.Add("L2ErrorEntropyL", XESFQueries.L2Error("h", "L", referenceSolution: (X, t) => enthalpy_inflow));
            //c.Queries.Add("L2ErrorEntropyR", XESFQueries.L2Error("h", "R", referenceSolution: (X, t) => enthalpy_inflow));
            //c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density));
            //c.Queries.Add("L2NormMomentumX", QueryLibrary.L2Norm(CompressibleVariables.Momentum[0]));
            //c.Queries.Add("L2NormEnergy", QueryLibrary.L2Norm(CompressibleVariables.Energy));
            #endregion

            #region Project Information
            c.ProjectName = "XESFWedgeFlow_TwoLs";
            c.SessionName = string.Format("WedgeFlow_TwoLs_p{0}_xCells{1}_yCells{2}_{3}_iterations{4}_angle{5}_lSDegree{6}_Glob{7}_agg{8}_",
                dgDegree, numOfCellsX, numOfCellsY, interfaceFluxLS2.ToString(), c.NoOfTimesteps, initialAngle_shockLS, optiLSDegree, globalization.ToString(), agg);

            #endregion

            return c;
        }
        /// <summary>
        /// Control for supersonic flow over an inclined plane (aka. WedgeFlow) discretized using 1 level set.
        /// To obtain the inclined plane a Cartesian grid with bottom boundary is rotated with prescribed wedge_angle.
        /// - the level Set represents the straight shock that is obtained
        /// - as long as the shock does not detach or vanish the solution is piecewise constant
        /// </summary>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="initialAngle_shockLS"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsDegree"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="FluxVersion"></param>
        /// <param name="shocksetup"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="optiLSDegree"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        /// <exception cref="NotImplementedException"></exception>
        public static XESFControl XDGWedgeFlow_OneLs_Rotation(int MaxIterations=100, int dgDegree=0, int numOfCellsX=10,
        int numOfCellsY=10, double initialAngle_shockLS=32, int PlotInterval=-1,
        string dbPath = null, int lsDegree = 1, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.OptimizedHLLC,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.GodunovInterface,
        XESF.Fluxes.FluxVersion FluxVersion = Fluxes.FluxVersion.Optimized, GetLevelSet shocksetup = GetLevelSet.FromParams, OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet, int optiLSDegree = 1,
        GetInitialValue initialValue = GetInitialValue.FromFunctionPerSpecies, bool iVFromShockRelations = true, double agg = 0.2, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch, MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit
        ) {


            XESFControl c = new XESFControl();
            #region DBStuff and Plotting
            c.DbPath = dbPath;
            c.savetodb = true;
            if(dbPath == null)
                c.savetodb = false;
            c.NoOfTimesteps = MaxIterations;
            c.ImmediatePlotPeriod = PlotInterval;
            if(dgDegree == 0 && numOfCellsX < 10 && numOfCellsY < 10) {
                c.SuperSampling = 4;
            } else {
                c.SuperSampling = 3;
            }

            #endregion

            #region SQP constants/parameters
            // ### Time config ###
            c.NoOfTimesteps = MaxIterations;

            //Globalization Parameters
            c.GlobalizationStrategy = globalization;
            c.MeritFunctionType = meritFunctionType;
            c.MinPIter = new int[] { 20, 20, 20 };
            #endregion

            #region LevelSetStuff

            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            c.IsTwoLevelSetRun = false;

            c.LevelSetDegree = lsDegree;
            c.AddVariable(XESFVariables.LevelSet, 1); //this is the levelSet used for the Immersed Boundary (the Wedge)

            c.SpeciesTable[0, 0] = "L"; // 'Forbidden' species
            c.SpeciesTable[1, 0] = "R"; // Void area (inside the blunt body)
            c.LsOne_NegSpecies = "L";
            c.LsOne_PosSpecies = "R";
            c.SpeciesToEvaluate = new string[] { "L", "R" };
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = optiLSDegree;

            #endregion 

            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.FluxVersion = FluxVersion;
            c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;
            #endregion

            #region Grid
            double xMin = 0;
            double xMax = 1.5;
            double yMin = 0;
            double yMax = 1.0;
            double wedge_angle_radial = 10.0 * Math.PI / 180.0; // Angle of Wedge

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                if(numOfCellsX > 1 && xNodes[1] >= 0.5) {
                    xNodes[1] = 0.4;
                }
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                var rotation = AffineTrafo.Some2DRotation(wedge_angle_radial);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");
                //grid.EdgeTagNames.Add(4, "SubsonicOutlet");

                grid = grid.Transform(rotation); //rotate the grid by the prescribed angle
                grid.DefineEdgeTags(delegate (double[] X) {
                    double[] rotX = rotation.Matrix.GetInverse().MatVecMul(1.0, X); //accounts for the rotation
                    //double[] rotX = X;
                    if(Math.Abs(rotX[0] - xMax) < 1e-14) {    // Right boundary
                                                              ////if((-X[0] + 0.5 + (X[1] / Math.Tan(wedge_angle)))>=0){ 
                                                              //if(X[1] >= Math.Tan(39.3139318 * Math.PI / 180.0)) { //-0.01 doesn't work
                                                              //    return 2;
                                                              //} else { // (part of void area)
                        return 2;
                        //return 3; // Wedge Part of Outlet

                    } else if(Math.Abs(rotX[0] - xMin) < 1e-14) { // Left boundary
                        return 1;
                    } else if(Math.Abs(rotX[1] - yMax) < 1e-14) { // top boundary
                        return 1;
                        //if(X[0] < 0.5) {
                        //    return 1;
                        //} else {
                        //    return 2;
                        //}

                    } else { //if(Math.Abs(rotX[1] - yMin) < 1e-14) { //bottom boundary
                        return 3;
                        //if(X[0] <= 1) {
                        //    return 2;
                        //} else {
                        //    return 2;
                        //}
                    }
                });
                return grid;
            };
            #endregion

            #region Initial Position of LevelSets
            double Mach = 2.0;
            double shock_angle_exact = 39.3139318;
            double shock_angle_radial;
            if(iVFromShockRelations) {
                shock_angle_radial = initialAngle_shockLS * Math.PI / 180.0; // this angle is the shock angle assumed in the beginning and used to calculate the IV
            } else {
                shock_angle_radial = shock_angle_exact * Math.PI / 180.0; // when we want to project the initial solution we need to use the exact angle in the calculations
            }
            double LevelSet1Prime = Math.Tan(wedge_angle_radial);
            double LevelSet2Prime = Math.Tan(initialAngle_shockLS * Math.PI / 180.0);


            // ### Wedge Level set function ###
            c.LevelSetOneInitialValue = delegate (double[] X) { return -X[0] + 0.5 + (X[1] / LevelSet2Prime); };

            //1=X[1]


            //// Shock level set
            //c.LevelSetOneInitialValue = delegate (double[] X) { return X[0] - 0.5 - (X[1] / LevelSet2Prime); };
            c.GetLevelSet = shocksetup;
            switch(shocksetup) {
                case GetLevelSet.FromParams:
                if(optiLSDegree == 0) {
                    c.OptiLevelSet_ParamNames = new List<string> { "(x-0.5) ", "t" };
                    c.OptiLevelSet_ParamValues = new List<double> { 1.0, -1.0 / LevelSet2Prime };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a * (x[0] - 0.5), (x, a) => x[1] * a };
                } else if(optiLSDegree == 1) {
                    c.OptiLevelSet_ParamNames = new List<string> { " ", "x", "t" };
                    c.OptiLevelSet_ParamValues = new List<double> { -0.0, 1.0, -1.0 / LevelSet2Prime };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a, (x, a) => a * x[0], (x, a) => x[1] * a };
                } else if(optiLSDegree == 2) {
                    c.OptiLevelSet_ParamNames = new List<string> { " ", "x", "t", "t^2" };
                    c.OptiLevelSet_ParamValues = new List<double> { -0.0, 1.0, -1.0 / LevelSet2Prime, 0.2 };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a, (x, a) => a * x[0], (x, a) => x[1] * a, (x, a) => x[1] * x[1] * a };

                } else if((optiLSDegree == 3)) {
                    c.OptiLevelSet_ParamNames = new List<string> { " ", "x", "t", "t^2", "t^3" };
                    c.OptiLevelSet_ParamValues = new List<double> { -0.0, 1.0, -1.0 / LevelSet2Prime, 0.2, -0.2 };
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a, (x, a) => a * x[0], (x, a) => x[1] * a, (x, a) => x[1] * x[1] * a, (x, a) => x[1] * x[1] * x[1] * a };
                } else {
                    throw new NotSupportedException("optiLSDegree:" + optiLSDegree + " not supported!");
                }

                break;
                case GetLevelSet.FromFunction:
                switch(optiLevelSetType) {
                    case OptiLevelSetType.SplineLevelSet:
                    //c.LevelSetTwoInitialValue = X => 0.0 - (X[1] / LevelSet2Prime);
                    throw new NotImplementedException("Spline Level Set is not supported with rotated Coordinate System");
                    break;
                    default:
                    if(optiLSDegree == 0) {
                        optiLSDegree++;
                        c.LevelSetTwoInitialValue = X => X[0] - 0.0 - (X[1] / LevelSet2Prime);
                    } else if(optiLSDegree == 1) {
                        c.LevelSetTwoInitialValue = X => X[0] - 0.0 - (X[1] / LevelSet2Prime);
                    } else if(optiLSDegree == 2) {
                        c.LevelSetTwoInitialValue = X => X[0] - 0.0 - (X[1] / LevelSet2Prime); //+ X[1] * X[1] * 0.2;
                    } else if((optiLSDegree == 3)) {
                        c.LevelSetTwoInitialValue = X => X[0] - 0.0 - (X[1] / LevelSet2Prime);// + X[1] * X[1] * 0.2 - X[1] * X[1] * X[1] * 0.2;
                    } else {
                        throw new NotSupportedException("optiLSDegree:" + optiLSDegree + " not supported!");
                    }
                    GlobalOptiLevelSet tmp_LS = GlobalOptiLevelSet.CreateONBLevelSet(xMin, xMax, yMin, yMax, optiLSDegree);
                    c.OptiLevelSet_ParamNames = new List<string>(tmp_LS.m_ParamNames);
                    c.OptiLevelSet_ParamValues = new List<double>(tmp_LS.m_ParamValues);
                    c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>(tmp_LS.m_phi);
                    c.OptiLSIsOrthonormal = true;
                    break;
                }

                break;
                case GetLevelSet.FromReconstruction:
                case GetLevelSet.FromReconstructionFromPoints:
                GlobalOptiLevelSet tmp_LS2 = GlobalOptiLevelSet.CreateONBLevelSet(xMin, xMax, yMin, yMax, optiLSDegree);
                c.OptiLevelSet_ParamNames = new List<string>(tmp_LS2.m_ParamNames);
                c.OptiLevelSet_ParamValues = new List<double>(tmp_LS2.m_ParamValues);
                c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>(tmp_LS2.m_phi);
                c.OptiLSIsOrthonormal = true;
                break;
                default:
                throw new NotSupportedException("Not supported Shock Setup");

            }
            #endregion

            #region Initial Values for Solution + Boundary Values
            c.SolDegree = dgDegree;
            c.GetInitialValue = initialValue;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;
            c.ExactEnthalpy=6.3;
            c.AddVariable(XESFVariables.Enthalpy_Error, dgDegree);
            c.AddVariable(XESFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESFVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(XESFVariables.Pressure, dgDegree);
            c.AddVariable(XESFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;


            double gamma = IdealGas.Air.HeatCapacityRatio;
            double speedOfSound = Math.Sqrt(gamma * pressure / density);

            double velocityX = Mach * speedOfSound;
            double velocityY = 0.0;

            double momentumX = density * velocityX;


            double innerEnergy = pressure / (gamma - 1);
            double kineticEnergy = 0.5 * density * (velocityX * velocityX + velocityY * velocityY);
            double totalEnergy = innerEnergy + kineticEnergy;

            // Boundary conditions in PRIMITIVE variables
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#L", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#R", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => pressure);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => pressure);


            double density_Right = (gamma + 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) / ((gamma - 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) + 2.0);
            double pressure_Right = 1.0 + (2.0 * gamma * (Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) - 1.0)) / (gamma + 1.0);

            // calculate Velocity via Transformation matrix
            MultidimensionalArray TransMat = MultidimensionalArray.Create(2, 2);
            TransMat[0, 0] = 1.0 / Math.Sqrt(1.0 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[0, 1] = Math.Tan(shock_angle_radial) / Math.Sqrt(1 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[1, 0] = Math.Tan(shock_angle_radial) / Math.Sqrt(1 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[1, 1] = -1.0 / Math.Sqrt(1.0 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));


            double[] velocity_left_normal = new double[2];
            velocity_left_normal[0] = velocityX * Math.Cos(shock_angle_radial);
            velocity_left_normal[1] = velocityX * Math.Sin(shock_angle_radial);


            double[] velocity_right_normal = new double[2];
            double[] velocity_right_cartesian = new double[2];


            velocity_right_normal[0] = velocity_left_normal[0];
            velocity_right_normal[1] = velocity_left_normal[1] / density_Right;

            TransMat.MatVecMul(1.0, velocity_right_normal, 0.0, velocity_right_cartesian);

            double velocityX_Right = velocity_right_cartesian[0];
            double velocityY_Right = velocity_right_cartesian[1];

            // calculate Velocity via Right-side Mach number

            //double Mach_Right = 1.6405;
            //double speedOfSoundRight = Math.Sqrt(gamma * pressure_Right / density_Right);
            //double velocity_norm = Mach_Right * speedOfSoundRight;
            //double velocityX_Right = velocity_norm;
            //double velocityY_Right = velocity_norm*Math.Tan(wedge_angle);

            double innerEnergy_Right = pressure_Right / (gamma - 1);
            double kineticEnergy_Right = 0.5 * density_Right * (velocityX_Right * velocityX_Right + velocityY_Right * velocityY_Right);
            double totalEnergy_Right = innerEnergy_Right + kineticEnergy_Right;

            double exact_sol(double[] X, int component) {
                //if(X[0] > 0.5 + (X[1] / LevelSet1Prime)) { //wedge -X[0] + 0.5 + (X[1] / LevelSet1Prime)
                //    return 0;
                if(X[0] >= (X[1] / Math.Tan(39.3139318 * Math.PI / 180.0))) {/*&& (X[0] <= 0.5 + (X[1] / LevelSet1Prime)))*/ //post Shock
                    if(component == 0) {
                        return density_Right;
                    } else if(component == 1) {
                        return velocityX_Right;
                    } else if(component == 2) {
                        return velocityY_Right;
                    } else if(component == 3) {
                        return pressure_Right;
                    }
                } else { //pre Shock
                    if(component == 0) {
                        return density;
                    } else if(component == 1) {
                        return velocityX;
                    } else if(component == 2) {
                        return 0;
                    } else if(component == 3) {
                        return pressure;
                    }
                }
                return 0;
            }

            if(iVFromShockRelations) {

                //// Initial conditions in PRIMITIVE variables -- Pre Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => density);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => velocityX);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => velocityY);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => pressure);

                //// Initial conditions in PRIMITIVE variables -- Post Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => velocityX_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => velocityY_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => pressure_Right);

                //// Initial conditions in Conservative variables -- Post Shock
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.xComponent + "#R", X => density_Right * velocityX_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.yComponent + "#R", X => density_Right * velocityY_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Energy + "#R", X => totalEnergy_Right);
            } else {
                //// Initial conditions in PRIMITIVE variables -- Pre Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => exact_sol(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => exact_sol(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => exact_sol(X, 2));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => exact_sol(X, 3));

                //// Initial conditions in PRIMITIVE variables -- Post Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => exact_sol(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => exact_sol(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => exact_sol(X, 2));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => exact_sol(X, 3));
            }
            #endregion

            #region Queries
            //var enthalpy_inflow = (gamma) / (gamma - 1) + 0.5 * velocityX * velocityX;
            //c.Queries.Add("L2ErrorEntropyL", XESFQueries.L2Error("h", "L", referenceSolution: (X, t) => enthalpy_inflow));
            //c.Queries.Add("L2ErrorEntropyR", XESFQueries.L2Error("h", "R", referenceSolution: (X, t) => enthalpy_inflow));
            //c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density));
            //c.Queries.Add("L2NormMomentumX", QueryLibrary.L2Norm(CompressibleVariables.Momentum[0]));
            //c.Queries.Add("L2NormEnergy", QueryLibrary.L2Norm(CompressibleVariables.Energy));
            #endregion

            #region Project Information
            c.ProjectName = "XESFWedgeFlow_OneLs";
            c.SessionName = string.Format("WedgeFlow_OneLs_p{0}_xCells{1}_yCells{2}_{3}_iterations{4}_angle{5}_lSDegree{6}_Flux{7}_agg{8}_",
                dgDegree, numOfCellsX, numOfCellsY, interfaceFluxLS1.ToString(), c.NoOfTimesteps, initialAngle_shockLS, optiLSDegree, interfaceFluxLS1.ToString(), agg);

            #endregion

            return c;
        }
        /// <summary>
        /// Control for supersonic flow over an inclined plane (aka. WedgeFlow) discretized using 1 level set.
        /// To obtain the inclined plane a grid/mesh must be provided via Gmsh (.msh)
        /// - the level Set represents the straight shock that is obtained
        /// - as long as the shock does not detach or vanish the solution is piecewise constant
        /// - may not work ATM, --> problems with non-Cartesian cut-cells
        /// TODO: Look into this issue
        /// </summary>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="initialAngle_shockLS"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsDegree"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="FluxVersion"></param>
        /// <param name="shocksetup"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <param name="meshPath"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        public static XESFControl XDGWedgeFlow_OneLs_Base_GridFromDB(int MaxIterations, int dgDegree, double initialAngle_shockLS, int PlotInterval,
        string dbPath = null, int lsDegree = 1, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.OptimizedHLLC,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.GodunovInterface,
        XESF.Fluxes.FluxVersion FluxVersion = Fluxes.FluxVersion.Optimized, GetLevelSet shocksetup = GetLevelSet.FromFunction, OptiLevelSetType optiLevelSetType = OptiLevelSetType.GlobalLevelSet,
        GetInitialValue initialValue = GetInitialValue.FromFunctionPerSpecies, bool iVFromShockRelations = false, double agg = 0.1, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch, MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit,
        string meshPath = @"..\..\..\Meshes\WFExactMesh4Elements.msh") {

            XESFControl c = new XESFControl();
            #region DBStuff and Plotting
            c.DbPath = dbPath;
            c.savetodb = true;
            if(dbPath == null)
                c.savetodb = false;
            c.NoOfTimesteps = MaxIterations;
            c.ImmediatePlotPeriod = PlotInterval;
            c.SuperSampling = 5;


            #endregion
            #region SQP constants/parameters
            // ### Time config ###
            c.NoOfTimesteps = MaxIterations;

            // Adaptive Regularization
            c.Gamma_Start = 1;
            c.Gamma_Max = 1;
            c.Gamma_Min = 1e-06;
            c.tauGamma = 1.5;
            c.L = 1;
            c.sigma_1 = 1e-2;
            c.sigma_2 = 1e-1;

            //Globalization Parameters
            c.GlobalizationStrategy = globalization;
            c.MeritFunctionType = meritFunctionType;
            c.Alpha_Start = 1;
            c.Alpha_Min = 1e-08;
            #endregion
            #region LevelSetStuff

            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.IsTwoLevelSetRun = false;
            c.LevelSetDegree = lsDegree;
            c.AddVariable(XESFVariables.LevelSet, lsDegree); //this is the levelSet used for the Immersed Boundary (the Wedge)
            c.SpeciesTable = new string[1, 2];
            c.SpeciesTable[0, 0] = "L"; // Left to the Shock
            c.SpeciesTable[0, 1] = "R"; // Right to the Shock
            c.LsOne_NegSpecies = "L";
            c.LsOne_PosSpecies = "R";
            c.SpeciesToEvaluate = new string[] { "L", "R" };
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = lsDegree;
            //shockLevelSet_Db = @"C:\BoSSS-experimental\internal\src\private-mag\XDGShock\Tests\bosss_db_levelSets.zip";
            //c.quadOrderFunc = (int[] A, int[] B, int[] C) => 4;
            #endregion 

            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.FluxVersion = FluxVersion;
            c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;
            #endregion

            #region Grid
            double xMin = 0;
            double xMax = 1.5;
            double yMin = 0;
            double yMax = 1.0;
            double wedge_angle_radial = 10.0 * Math.PI / 180.0; // Angle of Wedge

            c.getGridFrom = GetGridFrom.DB;
            c.MeshPath = meshPath;
            c.EdgTagNames = new string[] { "SupersonicInlet", "SupersonicOutlet", "AdiabaticSlipWall" };
            c.EdgTagFunc = delegate (double[] X) {
                if(Math.Abs(X[0] - xMax) < 1e-14 && Math.Abs(X[1] - 0.18) > 1e-14) {    // Right boundary
                    return 2;
                } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                    return 1;
                } else if(Math.Abs(X[1] - yMax) < 1e-14) { // top boundary
                    return 3;
                } else {  //bottom boundary
                    return 3;
                }
            };
            #endregion

            #region Initial Position of LevelSets
            #region Computation of angles and Primes
            double Mach = 2.0;
            double shock_angle_exact = 39.3139318;
            double shock_angle_radial;
            if(iVFromShockRelations) {
                shock_angle_radial = initialAngle_shockLS * Math.PI / 180.0; // this angle is the shock angle assumed in the beginning and used to calculate the IV
            } else {
                shock_angle_radial = shock_angle_exact * Math.PI / 180.0; // when we want to project the initial solution we need to use the exact angle in the calculations
            }
            double LevelSet1Prime = Math.Tan(wedge_angle_radial);
            double LevelSet2Prime = Math.Tan(initialAngle_shockLS * Math.PI / 180.0);
            #endregion

            // ### Shock Level set function ###
            c.LevelSetOneInitialValue = delegate (double[] X) { return X[0] - 0.5 - (X[1] / LevelSet2Prime); };

            c.GetLevelSet = shocksetup;
            switch(shocksetup) {
                case GetLevelSet.FromParams:
                throw new NotSupportedException("not supported");
                case GetLevelSet.FromFunction:
                GlobalOptiLevelSet tmp_LS = GlobalOptiLevelSet.CreateONBLevelSet(xMin, xMax, yMin, yMax, lsDegree);
                c.OptiLevelSet_ParamNames = new List<string>(tmp_LS.m_ParamNames);
                c.OptiLevelSet_ParamValues = new List<double>(tmp_LS.m_ParamValues);
                c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>(tmp_LS.m_phi);
                c.OptiLSIsOrthonormal = true;
                c.LevelSetTwoInitialValue = X => X[0] - 0.5 - (X[1] / LevelSet2Prime);
                break;
                default:
                     throw new NotSupportedException("Not supported Shock Setup");

            }
            #endregion

            #region Initial Values for Solution + Boundary Values
            c.SolDegree = dgDegree;
            c.GetInitialValue = initialValue;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;
            #region add variables
            c.AddVariable(XESFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESFVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(XESFVariables.Pressure, dgDegree);
            //c.AddVariable(XESFVariables.PressureStep, dgDegree); crashes - and unneeded
            c.AddVariable(XESFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            #endregion
            #region Values to left /at inlet
            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double speedOfSound = Math.Sqrt(gamma * pressure / density);
            double velocityX = Mach * speedOfSound;
            double velocityY = 0.0;
            double momentumX = density * velocityX;
            double innerEnergy = pressure / (gamma - 1);
            double kineticEnergy = 0.5 * density * (velocityX * velocityX + velocityY * velocityY);
            double totalEnergy = innerEnergy + kineticEnergy;

            // Boundary conditions in PRIMITIVE variables for the Inlet
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#L", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#R", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => pressure);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => pressure);
            #endregion
            #region Computation of values to the right of the shock
            double density_Right = (gamma + 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) / ((gamma - 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) + 2.0);
            double pressure_Right = 1.0 + (2.0 * gamma * (Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) - 1.0)) / (gamma + 1.0);

            // calculate Velocity via Transformation matrix
            MultidimensionalArray TransMat = MultidimensionalArray.Create(2, 2);
            TransMat[0, 0] = 1.0 / Math.Sqrt(1.0 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[0, 1] = Math.Tan(shock_angle_radial) / Math.Sqrt(1 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[1, 0] = Math.Tan(shock_angle_radial) / Math.Sqrt(1 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            TransMat[1, 1] = -1.0 / Math.Sqrt(1.0 + Math.Tan(shock_angle_radial) * Math.Tan(shock_angle_radial));
            double[] velocity_left_normal = new double[2];
            velocity_left_normal[0] = velocityX * Math.Cos(shock_angle_radial);
            velocity_left_normal[1] = velocityX * Math.Sin(shock_angle_radial);


            double[] velocity_right_normal = new double[2];
            double[] velocity_right_cartesian = new double[2];


            velocity_right_normal[0] = velocity_left_normal[0];
            velocity_right_normal[1] = velocity_left_normal[1] / density_Right;

            TransMat.MatVecMul(1.0, velocity_right_normal, 0.0, velocity_right_cartesian);

            double velocityX_Right = velocity_right_cartesian[0];
            double velocityY_Right = velocity_right_cartesian[1];

            // calculate Velocity via Right-side Mach number

            //double Mach_Right = 1.6405;
            //double speedOfSoundRight = Math.Sqrt(gamma * pressure_Right / density_Right);
            //double velocity_norm = Mach_Right * speedOfSoundRight;
            //double velocityX_Right = velocity_norm;
            //double velocityY_Right = velocity_norm*Math.Tan(wedge_angle);

            double innerEnergy_Right = pressure_Right / (gamma - 1);
            double kineticEnergy_Right = 0.5 * density_Right * (velocityX_Right * velocityX_Right + velocityY_Right * velocityY_Right);
            double totalEnergy_Right = innerEnergy_Right + kineticEnergy_Right;

            double exact_sol(double[] X, int component) {
                //if(X[0] > 0.5 + (X[1] / LevelSet1Prime)) { //wedge -X[0] + 0.5 + (X[1] / LevelSet1Prime)
                //    return 0;
                if(X[0] >= 0.5 + (X[1] / Math.Tan(39.3139318 * Math.PI / 180.0))) {/*&& (X[0] <= 0.5 + (X[1] / LevelSet1Prime)))*/ //post Shock
                    if(component == 0) {
                        return density_Right;
                    } else if(component == 1) {
                        return velocityX_Right;
                    } else if(component == 2) {
                        return velocityY_Right;
                    } else if(component == 3) {
                        return pressure_Right;
                    }
                } else { //pre Shock
                    if(component == 0) {
                        return density;
                    } else if(component == 1) {
                        return velocityX;
                    } else if(component == 2) {
                        return 0;
                    } else if(component == 3) {
                        return pressure;
                    }
                }
                return 0;
            }
            #endregion
            #region Set initial values
            if(iVFromShockRelations) {

                //// Initial conditions in PRIMITIVE variables -- Pre Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => density);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => velocityX);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => velocityY);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => pressure);

                //// Initial conditions in PRIMITIVE variables -- Post Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => velocityX_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => velocityY_Right);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => pressure_Right);

                //// Initial conditions in Conservative variables -- Post Shock
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.xComponent + "#R", X => density_Right * velocityX_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.yComponent + "#R", X => density_Right * velocityY_Right);
                //c.InitialValues_Evaluators.Add(CompressibleVariables.Energy + "#R", X => totalEnergy_Right);
            } else {
                //// Initial conditions in PRIMITIVE variables -- Pre Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => exact_sol(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => exact_sol(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => exact_sol(X, 2));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => exact_sol(X, 3));

                //// Initial conditions in PRIMITIVE variables -- Post Shock
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => exact_sol(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => exact_sol(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => exact_sol(X, 2));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => exact_sol(X, 3));
            }
            #endregion
            #endregion

            #region Queries
            var enthalpy_inflow = (gamma) / (gamma - 1) + 0.5 * velocityX * velocityX;
            c.Queries.Add("L2ErrorEntropyL", XESFQueries.L2Error("h", "L", referenceSolution: (X, t) => enthalpy_inflow));
            c.Queries.Add("L2ErrorEntropyR", XESFQueries.L2Error("h", "R", referenceSolution: (X, t) => enthalpy_inflow));
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density));
            c.Queries.Add("L2NormMomentumX", QueryLibrary.L2Norm(CompressibleVariables.Momentum[0]));
            c.Queries.Add("L2NormEnergy", QueryLibrary.L2Norm(CompressibleVariables.Energy));
            #endregion

            #region Project Information
            c.ProjectName = "XESFWedgeFlow_OneLs";
            c.SessionName = string.Format("WedgeFlow_OneLs_p{0}_iterations{1}_angle{2}_lSDegree{3}_Glob{4}_agg{5}_",
                dgDegree, c.NoOfTimesteps, initialAngle_shockLS, lsDegree, globalization.ToString(), agg);

            #endregion

            return c;


        }
        /// <summary>
        /// Control for supersonic flow over an inclined plane (aka. WedgeFlow) discretized using 2 Level Sets derived from Base Control.
        /// Here only the selected problem parameters may be controlled
        /// </summary>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="initialAngle_shockLS"></param>
        /// <param name="PlotInterval"></param>
        /// <returns></returns>
        public static XESFControl XDGWedgeFlow_TwoLs(int MaxIterations, int dgDegree, int numOfCellsX, int numOfCellsY, double initialAngle_shockLS, int PlotInterval) {

            XESFControl c = XDGWedgeFlow_TwoLs_Base(MaxIterations, dgDegree, numOfCellsX, numOfCellsY, initialAngle_shockLS, PlotInterval);

            double shock_angle_radial = initialAngle_shockLS * Math.PI / 180.0;
            double LevelSet2Prime = Math.Tan(shock_angle_radial);

            c.OptiLevelSet_ParamNames = new List<string> { "(x-0.5) ", "t" };
            c.OptiLevelSet_ParamValues = new List<double> { 1.0, -1.0 / LevelSet2Prime };
            c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>> { (x, a) => a * (x[0] - 0.5), (x, a) => x[1] * a };

            return c;
        }
        /// <summary>
        /// helper function to obtain controls for various testing and convergence study for the FDY cluster
        /// </summary>
        /// <param name="gammaMin"></param>
        /// <param name="agg"></param>
        /// <param name="tALNR"></param>
        /// <param name="TermN"></param>
        /// <param name="numY"></param>
        /// <param name="numX"></param>
        /// <param name="DegE"></param>
        /// <param name="DegS"></param>
        /// <param name="plotInterval"></param>
        /// <param name="iflux"></param>
        /// <param name="terStrat"></param>
        /// <param name="iRI"></param>
        /// <param name="itALNR"></param>
        /// <param name="itMinPIter"></param>
        /// <param name="iTermN"></param>
        /// <returns></returns>
        public static XESFControl XDGBS_Cluster(bool aRI = false, double gammaMin = 1e-04, double agg = 0.4, double tALNR = 1.001, int TermN = 8, int numY = 22, string dbPath = null,
            int numX = 10, int DegE = 3, int DegS = 0, int plotInterval = -1, int iProb=0,int iflux = 0, int terStrat = 0, int iRI = 0, int itALNR =0, int itMinPIter =0,int iTermN =0, int iRIT=0, int iFphi=0) 
            
            {
            var fphiTypes = new FphiType[] { FphiType.None, FphiType.CurvatureAll, FphiType.CurvatureCut, FphiType.PerssonSensorCut, FphiType.PerssonSensorAll };
            var OProblems = new OptProblemType[] { OptProblemType.FullEnRes, OptProblemType.EnResOnlyNearBand, OptProblemType.RankineHugoniotFull, OptProblemType.RankineHugoniotOnlyInterface };
            var IFluxes = new ConvectiveInterfaceFluxes[] { ConvectiveInterfaceFluxes.GodunovInterface,ConvectiveInterfaceFluxes.RoeInterface, ConvectiveInterfaceFluxes.CentralFluxInterface, ConvectiveInterfaceFluxes.OptimizedHLLCInterface };
            var ITerStrats = new TerminationStrategy[] { TerminationStrategy.Skyline, TerminationStrategy.SkylineWithDifferntTermNs }; // not needed right now but stays in case....
            var TermNArrays = new int[][] { new int[] { TermN, TermN, TermN, TermN, TermN, TermN }, new int[] { TermN, TermN, 5, 5, 5, 5 } };
            var MinPIterArrays = new int[][] { new int[] { 30, 30, 10, 10, 10, 10 }, new int[] { 30, 20, 5, 5, 5, 5 } };
            var tALNRArrays = new double[][] { new double[] { tALNR, tALNR, tALNR, tALNR, tALNR, tALNR }, new double[] { tALNR, 1.1, 1.1, 1.1, 1.1, 1.1 } };
            var MaxReinitArrays = new int[][] { new int[] { 30, 20, 5, 5, 5, 5 }, new int[] { 30, 15, 0, 0, 0, 0 } };
            var ReInitTolsArray = new double[][] { new double[] { -2e-1, -2e-1, -2e-1, -3e-1, -4e-1, 0 },new double[] { -2e-1, -2e-1, -2e-1, -2e-1, -2e-1, -2e-1 }, new double[] { -1, -0.2, -0.3, -0.4, -0.5 }, new double[] { -1, -0.2, -0.4, -0.6, -0.8} };

            var c = XDGBowShock_TwoLs_LSFromDB(
                agg: agg,
                numOfCellsX: numX, numOfCellsY: numY,
                dgDegreeStart: DegS, dgDegreeEnd: DegE,
                dbPath: dbPath,
                shockLevelSet_SessionId: @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2",
                pointPath: @".\BowShockPoints.txt",
                shockLevelSet_Db: @".\bosss_db_levelSets.zip",
                PlotInterval: plotInterval,
                optprob: OProblems[iProb],
                interfaceFluxLS2: IFluxes[iflux],
                terStrat: ITerStrats[terStrat],
                gammaMin: gammaMin,
                tALNRs: tALNRArrays[itALNR],
                TermNs: TermNArrays[iTermN],
                MinPIter: MinPIterArrays[itMinPIter],
                applyReInit:aRI,
                MaxReInits: MaxReinitArrays[iRI],
                ReInitTols: ReInitTolsArray[iRIT],
                fphitype:fphiTypes[iFphi]
                ) ;
                
                c.SessionName = string.Format($"XDGBS-p{DegE}-{numX}x{numY}-agg{agg}-iProb{iProb}-iFlux{iflux}-FphiType{iFphi}-aRI_{aRI}");
            return c;
        }
        /// <summary>
        /// helper function to obtain controls for various testing and convergence study
        /// </summary>
        /// <param name="gammaMin"></param>
        /// <param name="agg"></param>
        /// <param name="tALNR"></param>
        /// <param name="TermN"></param>
        /// <param name="numY"></param>
        /// <param name="numX"></param>
        /// <param name="DegE"></param>
        /// <param name="DegS"></param>
        /// <param name="plotInterval"></param>
        /// <param name="iflux"></param>
        /// <param name="terStrat"></param>
        /// <param name="iRI"></param>
        /// <param name="itALNR"></param>
        /// <param name="itMinPIter"></param>
        /// <param name="iTermN"></param>
        /// <returns></returns>
        public static XESFControl XDGBS_Local(double gammaMin = 1e-04, double agg = 0.4, double tALNR = 1.001, int TermN = 10, int cflux = 0,
            int numY = 22, int numX = 10, int DegE = 3, int DegS = 0, int plotInterval = -1, int iflux = 0, int terStrat = 0, int iRI = 0, int itALNR = 0, int itMinPIter = 0, int iTermN = 0) { 


           
            var IFluxes = new ConvectiveInterfaceFluxes[] { ConvectiveInterfaceFluxes.GodunovInterface, ConvectiveInterfaceFluxes.RoeInterface, ConvectiveInterfaceFluxes.CentralFluxInterface, ConvectiveInterfaceFluxes.OptimizedHLLCInterface };
            var CFluxes = new ConvectiveBulkFluxes[] { ConvectiveBulkFluxes.OptimizedHLLC,ConvectiveBulkFluxes.Roe };
            var ITerStrats = new TerminationStrategy[] { TerminationStrategy.Skyline, TerminationStrategy.MaxVsMin, TerminationStrategy.Experimental };
            var TermNArrays = new int[][] { new int[] { TermN, TermN, TermN, TermN, TermN, TermN }, new int[] { TermN, TermN, 5, 5, 5, 5 } };
            var MinPIterArrays = new int[][] { new int[] { 30, 30, 10, 10, 10, 10 }, new int[] { TermN, TermN, 5, 5, 5, 5 } };
            var tALNRArrays = new double[][] { new double[] { tALNR, tALNR, tALNR, tALNR, tALNR, tALNR }, new double[] { tALNR, 1.1, 1.1, 1.1, 1.1, 1.1 } };
            var MaxReinitArrays = new int[][] { new int[] { 30, 20, 5, 5, 5, 5 }, new int[] { 30, 15, 0, 0, 0, 0 } };
            var c = XDGBowShock_TwoLs_LSFromDB(agg: agg,
                numOfCellsX: numX, numOfCellsY: numY,
                dgDegreeStart: DegS, dgDegreeEnd: DegE,
                //MArkus AV Run *************
                ///Uni PC
                //shockLevelSet_Db: @"C:\experimental\internal\src\private-mag\XDGShock\Tests\bosss_db_levelSets.zip",
                //shockLevelSet_SessionId: @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2",
                //pointPath: @"C:\experimental\internal\src\private-mag\XDGShock\Tests\BowShockPoints.txt",
                initialValue: GetInitialValue.FromDBSinglePhase,
                ///Home PC
                shockLevelSet_Db: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\Databases\bosss_db_levelSets\bosss_db_levelSets",
                shockLevelSet_SessionId: @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2",
                pointPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-mag\XDGShock\Tests\BowShockPoints.txt",
                //initialValue: GetInitialValue.FromDBSinglePhase,PlotInterval: plotInterval,
                interfaceFluxLS2: IFluxes[iflux],
                bulkFlux: CFluxes[cflux],
                terStrat: ITerStrats[terStrat],
                gammaMin: gammaMin,
                tALNRs: tALNRArrays[itALNR],
                TermNs: TermNArrays[iTermN],
                MinPIter: MinPIterArrays[itMinPIter],
                MaxReInits: MaxReinitArrays[iRI],
                PlotInterval:plotInterval
                
                );

        c.SessionName = string.Format($"XDGBS-p{DegE}-{numX}x{numY}-agg{agg}-tALNR{itALNR}-TermN{iTermN}-MinPI{itMinPIter}-maxRI{iRI}"); 
            return c;
        }

        /// <summary>
        /// Control for the BowShock Problem (setup from HiOCFD5 workshop 2017 https://how5.cenaero.be/content/ci1-inviscid-bow-shock). Initial guess is loaded from a database from a previous AV run. 
        /// Level Set is loaded from (x,y) points describing the Interface position. 
        /// Originally also obtained from Level Set Reconstruction procedure (see. Markus Geisenhofer Dissertation (2021)
        /// </summary>
        /// <param name="optiLSDegree"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegreeEnd"></param>
        /// <param name="dgDegreeStart"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="IVtsNumber"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsTwoDegree"></param>
        /// <param name="lsOneDegree"></param>
        /// <param name="optprob"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="fphitype"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="interfaceFluxLS2"></param>
        /// <param name="FluxVersion"></param>
        /// <param name="shockLevelSet_Db"></param>
        /// <param name="shockLevelSet_SessionId"></param>
        /// <param name="pointPath"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="tALNRs"></param>
        /// <param name="TermNs"></param>
        /// <param name="MaxReInits"></param>
        /// <param name="ReInitTols"></param>
        /// <param name="gammaMax"></param>
        /// <param name="gammaMin"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <param name="solverRunType"></param>
        /// <param name="MinPIter"></param>
        /// <param name="terStrat"></param>
        /// <param name="staggeredTS"></param>
        /// <param name="applyReInit"></param>
        /// <param name="getLevelSet"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        public static XESFControl XDGBowShock_TwoLs_LSFromDB(int optiLSDegree=3, int MaxIterations=500, int dgDegreeEnd=3, int dgDegreeStart=0, int numOfCellsX=10,
        int numOfCellsY=22, int IVtsNumber = 39, int PlotInterval = -1,
        string dbPath = null, int lsTwoDegree = 3, int lsOneDegree = 4, OptProblemType optprob= OptProblemType.FullEnRes, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.OptimizedHLLC, FphiType fphitype=FphiType.None,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var,
        ConvectiveInterfaceFluxes interfaceFluxLS2 = ConvectiveInterfaceFluxes.GodunovInterface,
        Fluxes.FluxVersion FluxVersion = Fluxes.FluxVersion.Optimized, string shockLevelSet_Db = null, string shockLevelSet_SessionId = @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2", string pointPath = null,
         OptiLevelSetType optiLevelSetType = OptiLevelSetType.SplineLevelSet, double[] tALNRs = null, int[] TermNs = null, int[] MaxReInits = null, double[] ReInitTols=null, double gammaMax=1, double gammaMin = 1e-2,
        GetInitialValue initialValue = GetInitialValue.FromDBSinglePhase, bool iVFromShockRelations = false,
        double agg = 0.2, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch,
        MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit, SolverRunType solverRunType = SolverRunType.PContinuation, int[] MinPIter= null, TerminationStrategy terStrat =TerminationStrategy.Skyline,
        int[] staggeredTS = null, bool applyReInit=false, GetLevelSet getLevelSet = GetLevelSet.FromReconstructionFromPoints) {
            XESFControl c = new XESFControl();

            int dgDegree = dgDegreeEnd;
            c.SolDegree = dgDegreeEnd;
            #region database stuff
            c.DbPath = dbPath;
            c.savetodb = true;
            if(dbPath == null)
                c.savetodb = false;
            c.saveperiod = 1;
            c.ImmediatePlotPeriod = PlotInterval;
            c.SuperSampling = 2;
            c.NoOfTimesteps = MaxIterations;
            bool restart = false;
            c.IVTimestepNumber = IVtsNumber;
            #endregion
            #region Optimization variables
            c.solRunType = solverRunType;
            if(staggeredTS != null)
                c.staggeredTimeSteps = staggeredTS;
            c.MinPIter = MinPIter;
            // Adaptive Regularization
            c.Gamma_Max = gammaMax;
            c.Gamma_Min = gammaMin;
            c.optProblemType = optprob;
            c.GlobalizationStrategy = globalization;
            c.GetInitialValue = initialValue;
            c.MeritFunctionType = meritFunctionType;
            c.fphiType = fphitype;

            
            // ### Time config ###
            c.NoOfTimesteps = MaxIterations;

            // Reinitialization
            c.ApplyReiInit = applyReInit;
            c.reInitTols = ReInitTols != null? ReInitTols: new double[] { -2e-1, -2e-1, 0, 0, 0, 0 };
            //Termination
            c.terStrat = terStrat;
            c.TerminationMinNs = TermNs != null ? TermNs:new int[] { 8, 8, 8, 8, 8, 8 }; ;
            c.tALNRs = tALNRs != null ? tALNRs:new double[] { 1.005, 1.005, 1.001, 1.01, 1.01, 1.01 };
            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            //c.NonlinearQuadratureDegree = 6;
            #endregion
            #region Grid
            double xMin = -2.0;
            double xMax = 0.0;
            double yMin = -4.0;
            double yMax = 4.0;
            // Shift grid (the same way the stored simulation is shifted
            xMin = xMin - 0.025;
            xMax = xMax - 0.025;

            if(restart == true) {
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("3b7ec178-3859-48c6-91c6-13c566ee6246"), new TimestepNumber(1));
                c.GridGuid = new Guid("650eb82f-3642-4150-9cc6-f6925f46ef63");
            } else {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                            if(Math.Abs(X[1]) - 0.1 < 1e-14) { // Right boundary (part of void area)
                                return 2;
                            } else {
                                return 2;
                            }
                            //return 2;
                        } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                            return 1;
                        } else {    // Top and bottom boundary
                            return 3;
                        }
                    });

                    //var gDat = new GridData(grid);
                    //var em1 = gDat.GetBoundaryEdges();
                    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                    return grid;
                };
            }
            #endregion
            #region LevelSets
            /// Species Setting
            c.SpeciesTable[0, 0] = "X"; // 'Forbidden' species
            c.SpeciesTable[0, 1] = "V"; // Void area (inside the blunt body)
            c.SpeciesTable[1, 0] = "L"; // Pre-shock region (on the left of the shock)
            c.SpeciesTable[1, 1] = "R"; // Post-shock region (on the right of the shock)

            c.LsOne_NegSpecies = "V";
            c.LsOne_PosSpecies = "R";
            c.LsTwo_NegSpecies = "L";
            c.LsTwo_PosSpecies = "R";
            c.SpeciesToEvaluate = new string[] { "L", "R" };

            //degree Setting
            c.AddVariable(XESFVariables.LevelSet, lsOneDegree);
            c.AddVariable(XESFVariables.LevelSetTwo, lsTwoDegree);
            c.LevelSetDegree = lsOneDegree;
            c.LevelSetTwoDegree = lsTwoDegree;
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = optiLSDegree;

            //Initial Value for LevelSet
            c.GetLevelSet = getLevelSet;
            //string shockLevelSet_SessionId = @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2";
            //string shockLevelSet_SessionId = @"3a743670-0220-4f90-81b6-b28a027944d5";
            c.ShockLevelSet_Db = shockLevelSet_Db ?? throw new NotSupportedException("Shock level set DB is null.");
            c.ShockLevelSet_Info = new Tuple<Guid, TimestepNumber>(new Guid(shockLevelSet_SessionId), -1);
            c.ShockLevelSet_FieldName = "levelSet_recon_prc_cont";
            c.ShockLevelSet_SeedFromDb = true;

            /// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var;
            c.ConvectiveInterfaceFlux_LsTwo = interfaceFluxLS2;

            //c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;

            c.PointPath = pointPath;
            if(pointPath == null){
                c.PointPath = @"C:\experimental\internal\src\private-mag\XDGShock\Tests\BowShockPoints.txt";
                //c.PointPath = @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_P0_10x20.txt";
                c.PointPath = @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_P0_10x20.txt";
            }
            
            

            if(optiLevelSetType == OptiLevelSetType.GlobalLevelSet) {
                GlobalOptiLevelSet tmp_LS = GlobalOptiLevelSet.CreateONBLevelSet(xMin, xMax, yMin, yMax, optiLSDegree);
                c.OptiLevelSet_ParamNames = new List<string>(tmp_LS.m_ParamNames);
                c.OptiLevelSet_ParamValues = new List<double>(tmp_LS.m_ParamValues);
                c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>(tmp_LS.m_phi);
                c.OptiLSIsOrthonormal = true;
                c.OptiLevelSetDegree = optiLSDegree;
            }

            // ### Level set function of the blunt body###
            c.LevelSetOneInitialValue = delegate (double[] X) {
                // Circle 1
                double x0 = 0.0;
                double y0 = 0.5;
                double r0 = 0.5;

                // Circle 2
                double x1 = 0.0;
                double y1 = -0.5;
                double r1 = 0.5;

                // Signed distance formulation
                //if (X[1] >= 0.5) {
                //    return Math.Sqrt((X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0)) - r0;
                //} else if (X[1] <= -0.5) {
                //    return Math.Sqrt((X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1)) - r1;
                //} else {
                //    return -(X[0] + 0.5);
                //}

                // Quadratic formulation
                if(X[1] >= 0.5) {
                    return (X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0) - r0 * r0;
                } else if(X[1] <= -0.5) {
                    return (X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1) - r1 * r1;
                } else {
                    return X[0] * X[0] - 0.5 * 0.5;
                }
            };
            #endregion
            #region some general stuff (mostly unused)


            
 
            #endregion
            #region     Boundary Value stuff

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            //c.AddVariable(XESFVariables.Velocity.xComponent, dgDegree);
            //c.AddVariable(XESFVariables.Velocity.yComponent, dgDegree);
            //c.AddVariable(XESFVariables.Pressure, dgDegree);

            //c.AddVariable(XESFVariables.Enthalpy, dgDegree);
            //c.AddVariable(XESFVariables.LocalMachNumber, dgDegree);

            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;
            double Mach = 4.0;
            double velocityX = Mach * Math.Sqrt(c.EquationOfState.HeatCapacityRatio * pressure / density);
            double velocityY = 0.0;

            double momentumX = density * velocityX;


            double gamma = IdealGas.Air.HeatCapacityRatio;
            double innerEnergy = pressure / (gamma - 1);
            double kineticEnergy = 0.5 * density * (velocityX * velocityX + velocityY * velocityY);
            double totalEnergy = innerEnergy + kineticEnergy;

            // Boundary conditions in PRIMITIVE variables
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#L", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#R", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => pressure);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => pressure);

            // In theory no outflow boundary condition has to be specified as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", XESFVariables.Pressure + "#L", (X, t) => 0.0);
            c.AddBoundaryValue("SupersonicOutlet", XESFVariables.Pressure + "#R", (X, t) => 0.0);
            //c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if(restart == false) {
                // Initial conditions in PRIMITIVE variables
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => density);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => velocityX);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => velocityY);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => pressure);

                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => velocityX);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => velocityY);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => pressure);
            }
            #endregion
            
            c.DgDegree_Start = dgDegreeStart;
            c.ProjectName = "XDGBowShock_TwoLs";
            c.SessionName = string.Format($"XDGBowShock_TwoLs_p{dgDegree}_x{numOfCellsX}_y{numOfCellsY}_agg{c.AgglomerationThreshold}_"
                + $"{c.solRunType.ToString()}_LS2Flux{interfaceFluxLS2.ToString()}");
            if(restart == true) {
                c.SessionName = c.SessionName + "_RESTART";
            }
            
            //## Additional Variables
            c.AddVariable(XESFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESFVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(XESFVariables.Pressure, dgDegree);
            //c.AddVariable(XESFVariables.PressureStep, dgDegree);
            c.AddVariable(XESFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESFVariables.Enthalpy_Error, dgDegree);
            c.ExactEnthalpy = 14.7;
            c.AddVariable(XESFVariables.LocalMachNumber, dgDegree );

            // ### Queries ###
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density));
            c.Queries.Add("L2NormMomentumX", QueryLibrary.L2Norm(CompressibleVariables.Momentum[0]));
            c.Queries.Add("L2NormEnergy", QueryLibrary.L2Norm(CompressibleVariables.Energy));

            return c;
        }
        /// <summary>
        /// COntrol object for Bow shock where initial value and Level Set are seeded as wished. Not in use ATM
        /// </summary>
        /// <param name="optiLSDegree"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegreeEnd"></param>
        /// <param name="dgDegreeStart"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="IVtsNumber"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsTwoDegree"></param>
        /// <param name="lsOneDegree"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="interfaceFluxLS2"></param>
        /// <param name="FluxVersion"></param>
        /// <param name="shockLevelSet_Db"></param>
        /// <param name="shockLevelSet_SessionId"></param>
        /// <param name="pointPath"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="tALNRs"></param>
        /// <param name="TermNs"></param>
        /// <param name="MaxReInits"></param>
        /// <param name="ReInitTols"></param>
        /// <param name="gammaMax"></param>
        /// <param name="gammaMin"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <param name="solverRunType"></param>
        /// <param name="MinPIter"></param>
        /// <param name="terStrat"></param>
        /// <param name="staggeredTS"></param>
        /// <param name="applyReInit"></param>
        /// <param name="getLevelSet"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        public static XESFControl XDGBowShock_TwoLs(int optiLSDegree = 3, int MaxIterations = 200, int dgDegreeEnd = 3, int dgDegreeStart = 0, int numOfCellsX = 10,
       int numOfCellsY = 22, int IVtsNumber = 39, int PlotInterval = -1,
       string dbPath = null, int lsTwoDegree = 3, int lsOneDegree = 4, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.OptimizedHLLC,
       ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var,
       ConvectiveInterfaceFluxes interfaceFluxLS2 = ConvectiveInterfaceFluxes.GodunovInterface,
       Fluxes.FluxVersion FluxVersion = Fluxes.FluxVersion.Optimized, string shockLevelSet_Db = null, string shockLevelSet_SessionId = @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2", string pointPath = null,
        OptiLevelSetType optiLevelSetType = OptiLevelSetType.SplineLevelSet, double[] tALNRs = null, int[] TermNs = null, int[] MaxReInits = null, double[] ReInitTols = null, double gammaMax = 1, double gammaMin = 1e-2,
       GetInitialValue initialValue = GetInitialValue.FromDBSinglePhase, bool iVFromShockRelations = false,
       double agg = 0.2, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch,
       MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit, SolverRunType solverRunType = SolverRunType.PContinuation, int[] MinPIter = null, TerminationStrategy terStrat = TerminationStrategy.Skyline,
       int[] staggeredTS = null, bool applyReInit = true, GetLevelSet getLevelSet = GetLevelSet.FromReconstructionFromPoints) {
            XESFControl c = new XESFControl();

            int dgDegree = dgDegreeEnd;
            c.SolDegree = dgDegreeEnd;
            #region database stuff
            c.DbPath = dbPath;
            c.savetodb = true;
            if(dbPath == null)
                c.savetodb = false;
            c.saveperiod = 1;
            c.ImmediatePlotPeriod = PlotInterval;
            c.SuperSampling = 2;
            c.NoOfTimesteps = MaxIterations;
            bool restart = false;
            if(shockLevelSet_Db == null) {
                shockLevelSet_Db = @"C:\experimental\internal\src\private-mag\XDGShock\Tests\bosss_db_levelSets.zip";
                //shockLevelSet_Db = @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_P0_db";
            }
            c.IVTimestepNumber = IVtsNumber;
            #endregion
            #region Optimization variables
            c.solRunType = solverRunType;
            if(staggeredTS != null)
                c.staggeredTimeSteps = staggeredTS;
            c.MinPIter = MinPIter;
            // Adaptive Regularization
            c.Gamma_Start = 1;
            c.Gamma_Max = gammaMax;
            c.Gamma_Min = gammaMin;
            c.tauGamma = 1.5;
            c.L = 1;
            c.sigma_1 = 1e-2;
            c.sigma_2 = 1e-1;
            //Globalization Parameters
            c.Alpha_Start = 1;
            c.Alpha_Min = 1e-8;

            c.GlobalizationStrategy = globalization;
            c.GetInitialValue = initialValue;
            c.MeritFunctionType = meritFunctionType;
            //shockLevelSet_Db = @"C:\BoSSS-experimental\internal\src\private-mag\XDGShock\Tests\bosss_db_levelSets.zip";

            c.reInit_c1 = -1.5;

            // ### Time config ###
            c.NoOfTimesteps = MaxIterations;

            // Reinitialization
            c.ApplyReiInit = applyReInit;
            c.reInitTols = ReInitTols != null ? ReInitTols : new double[] { -2e-1, -2e-1, 0, 0, 0, 0 };
            //Termination
            c.terStrat = terStrat;
            c.TerminationMinNs = TermNs != null ? TermNs : new int[] { 8, 8, 8, 8, 8, 8 }; ;
            c.tALNRs = tALNRs != null ? tALNRs : new double[] { 1.005, 1.005, 1.001, 1.01, 1.01, 1.01 };
            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            //c.NonlinearQuadratureDegree = 6;
            #endregion
            #region Grid
            double xMin = -2.0;
            double xMax = 0.0;
            double yMin = -4.0;
            double yMax = 4.0;
            // Shift grid (the same way the stored simulation is shifted
            xMin = xMin - 0.025;
            xMax = xMax - 0.025;

            if(restart == true) {
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("3b7ec178-3859-48c6-91c6-13c566ee6246"), new TimestepNumber(1));
                c.GridGuid = new Guid("650eb82f-3642-4150-9cc6-f6925f46ef63");
            } else {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                            if(Math.Abs(X[1]) - 0.1 < 1e-14) { // Right boundary (part of void area)
                                return 2;
                            } else {
                                return 2;
                            }
                            //return 2;
                        } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                            return 1;
                        } else {    // Top and bottom boundary
                            return 3;
                        }
                    });

                    //var gDat = new GridData(grid);
                    //var em1 = gDat.GetBoundaryEdges();
                    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                    return grid;
                };
            }
            #endregion
            #region LevelSets
            /// Species Setting
            c.SpeciesTable[0, 0] = "X"; // 'Forbidden' species
            c.SpeciesTable[0, 1] = "V"; // Void area (inside the blunt body)
            c.SpeciesTable[1, 0] = "L"; // Pre-shock region (on the left of the shock)
            c.SpeciesTable[1, 1] = "R"; // Post-shock region (on the right of the shock)

            c.LsOne_NegSpecies = "V";
            c.LsOne_PosSpecies = "R";
            c.LsTwo_NegSpecies = "L";
            c.LsTwo_PosSpecies = "R";
            c.SpeciesToEvaluate = new string[] { "L", "R" };

            //degree Setting
            c.AddVariable(XESFVariables.LevelSet, lsOneDegree);
            c.AddVariable(XESFVariables.LevelSetTwo, lsTwoDegree);
            c.LevelSetDegree = lsOneDegree;
            c.LevelSetTwoDegree = lsTwoDegree;
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = optiLSDegree;

            //Initial Value for LevelSet
            c.GetLevelSet = getLevelSet;
            //string shockLevelSet_SessionId = @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2";
            //string shockLevelSet_SessionId = @"3a743670-0220-4f90-81b6-b28a027944d5";
            c.ShockLevelSet_Db = shockLevelSet_Db ?? throw new NotSupportedException("Shock level set DB is null.");
            c.ShockLevelSet_Info = new Tuple<Guid, TimestepNumber>(new Guid(shockLevelSet_SessionId), -1);
            c.ShockLevelSet_FieldName = "levelSet_recon_prc_cont";
            c.ShockLevelSet_SeedFromDb = true;

            /// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var;
            c.ConvectiveInterfaceFlux_LsTwo = interfaceFluxLS2;

            //c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;

            c.PointPath = pointPath;
            if(pointPath == null) {
                c.PointPath = @"C:\experimental\internal\src\private-mag\XDGShock\Tests\BowShockPoints.txt";
                //c.PointPath = @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_P0_10x20.txt";
                c.PointPath = @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_P0_10x20.txt";
            }



            if(optiLevelSetType == OptiLevelSetType.GlobalLevelSet) {
                GlobalOptiLevelSet tmp_LS = GlobalOptiLevelSet.CreateONBLevelSet(xMin, xMax, yMin, yMax, optiLSDegree);
                c.OptiLevelSet_ParamNames = new List<string>(tmp_LS.m_ParamNames);
                c.OptiLevelSet_ParamValues = new List<double>(tmp_LS.m_ParamValues);
                c.OptiLevelSet_Param_Functions = new List<Func<double[], double, double>>(tmp_LS.m_phi);
                c.OptiLSIsOrthonormal = true;
                c.OptiLevelSetDegree = optiLSDegree;
            }

            // ### Level set function of the blunt body###
            c.LevelSetOneInitialValue = delegate (double[] X) {
                // Circle 1
                double x0 = 0.0;
                double y0 = 0.5;
                double r0 = 0.5;

                // Circle 2
                double x1 = 0.0;
                double y1 = -0.5;
                double r1 = 0.5;

                // Signed distance formulation
                //if (X[1] >= 0.5) {
                //    return Math.Sqrt((X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0)) - r0;
                //} else if (X[1] <= -0.5) {
                //    return Math.Sqrt((X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1)) - r1;
                //} else {
                //    return -(X[0] + 0.5);
                //}

                // Quadratic formulation
                if(X[1] >= 0.5) {
                    return (X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0) - r0 * r0;
                } else if(X[1] <= -0.5) {
                    return (X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1) - r1 * r1;
                } else {
                    return X[0] * X[0] - 0.5 * 0.5;
                }
            };
            #endregion
            
            #region     Boundary Value stuff

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            //c.AddVariable(XESFVariables.Velocity.xComponent, dgDegree);
            //c.AddVariable(XESFVariables.Velocity.yComponent, dgDegree);
            //c.AddVariable(XESFVariables.Pressure, dgDegree);

            //c.AddVariable(XESFVariables.Enthalpy, dgDegree);
            //c.AddVariable(XESFVariables.LocalMachNumber, dgDegree);

            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;
            double Mach = 4.0;
            double velocityX = Mach * Math.Sqrt(c.EquationOfState.HeatCapacityRatio * pressure / density);
            double velocityY = 0.0;

            double momentumX = density * velocityX;


            double gamma = IdealGas.Air.HeatCapacityRatio;
            double innerEnergy = pressure / (gamma - 1);
            double kineticEnergy = 0.5 * density * (velocityX * velocityX + velocityY * velocityY);
            double totalEnergy = innerEnergy + kineticEnergy;

            // Boundary conditions in PRIMITIVE variables
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#L", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.yComponent + "#R", (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => pressure);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => pressure);

            // In theory no outflow boundary condition has to be specified as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", XESFVariables.Pressure + "#L", (X, t) => 0.0);
            c.AddBoundaryValue("SupersonicOutlet", XESFVariables.Pressure + "#R", (X, t) => 0.0);
            //c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if(restart == false) {
                // Initial conditions in PRIMITIVE variables
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => density);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => velocityX);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#L", X => velocityY);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => pressure);

                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => density);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => velocityX);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.yComponent + "#R", X => velocityY);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => pressure);
            }
            #endregion

            c.DgDegree_Start = dgDegreeStart;
            c.ProjectName = "XDGBowShock_TwoLs";
            c.SessionName = string.Format($"XDGBowShock_TwoLs_p{dgDegree}_x{numOfCellsX}_y{numOfCellsY}_agg{c.AgglomerationThreshold}_"
                + $"{c.solRunType.ToString()}_LS2Flux{interfaceFluxLS2.ToString()}");
            if(restart == true) {
                c.SessionName = c.SessionName + "_RESTART";
            }

            //## Additional Variables
            c.AddVariable(XESFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESFVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(XESFVariables.Pressure, dgDegree);
            //c.AddVariable(XESFVariables.PressureStep, dgDegree);
            c.AddVariable(XESFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESFVariables.Enthalpy_Error, dgDegree);
            c.ExactEnthalpy = 14.7;
            c.AddVariable(XESFVariables.LocalMachNumber, dgDegree);

            // ### Queries ###
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density));
            c.Queries.Add("L2NormMomentumX", QueryLibrary.L2Norm(CompressibleVariables.Momentum[0]));
            c.Queries.Add("L2NormEnergy", QueryLibrary.L2Norm(CompressibleVariables.Energy));

            return c;
        }
    }
}

