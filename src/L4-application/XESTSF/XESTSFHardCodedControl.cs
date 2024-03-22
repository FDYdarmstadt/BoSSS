using ApplicationWithIDT;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using static BoSSS.Solution.CompressibleFlowCommon.CompressibleHelperFunc;
using ilPSP;
using ilPSP.Utils;
using XESF.Fluxes;
using XESF.Variables;
using XESTSF.Variables;

using ApplicationWithIDT.OptiLevelSets;
using System.Reflection.PortableExecutable;

namespace XESTSF {
    public static class XESTSFHardCodedControl {

        /// <summary>
        /// Simulates an acoustic wave in 1D space via space-time method.
        /// </summary>
        /// <param name="MachL">Mach number of the base flow in the left region.</param>
        /// <param name="xMin">Left spacial domain boundary.</param>
        /// <param name="xMax">Right spacial domain boundary.</param>
        /// <param name="shockPosition">Position of the shock wave.</param>
        /// <param name="endTime">Time duration of the simulation.</param>
        /// <param name="numOfCellsX">Number of spatial cells in the x-axis.</param>
        /// <param name="numOfCellsT">Number of temporal cells.</param>
        /// <param name="dgDegree">Degree of the Discontinuous Galerkin method.</param>
        /// <param name="lsDegree">Degree of the Level Set.</param>
        /// <param name="withShock">Flag to include the shock wave in the simulation.</param>
        /// <param name="withLevelSet">Flag to include the Level Set.</param>
        /// <param name="MaxIterations">Maximum number of non-linear solver iterations.</param>
        /// <param name="agg">Aggregation factor.</param>
        /// <param name="PlotInterval">Interval for plotting.</param>
        /// <param name="dbPath">Path for database storage (can be null).</param>
        /// <param name="bulkFlux">Convective bulk flux method.</param>
        /// <param name="interfaceFluxLS1">Convective interface flux for Level Set.</param>
        /// <param name="p_amp_neg">Amplitude of negative acoustic wave.</param>
        /// <param name="p_amp_pos">Amplitude of positive acoustic wave.</param>
        /// <param name="waveform">Type of waveform (e.g., sinusoidal).</param>
        /// <param name="waveLength">Wavelength of the wave.</param>
        /// <param name="wavePosition">Position of the wave in the simulation space.</param>
        /// <returns> XESTSFControl object for controlling the simulation.</returns>
        /// <exception cref="NotSupportedException"></exception>
        public static XESTSFControl AcousticWave1D(double MachL = 4.0, double xMin = 0.0, double xMax = 3.0, double shockPosition = 2.0, double endTime = 1.0,
        int numOfCellsX = 10, int numOfCellsT = 10, int dgDegree = 3, int lsDegree = 3, double scaling = 10000,
        bool withShock = true, bool withLevelSet = true, int MaxIterations = 100, double agg = 0.4, int PlotInterval = 1, string dbPath = null,
        ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.Roe, ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.RoeInterface,
        double p_amp_neg = 0.00001, double p_amp_pos = 0.0,  string waveform = "1sinus", double waveLength = 0.4,  double wavePosition = 0.1) {
            
            var c = new XESTSFControl();
            c.tNref = numOfCellsT;
            #region DBStuff and Plotting
            c.DbPath = dbPath;
            c.savetodb = true;
            if(dbPath == null)
                c.savetodb = false;
            c.NoOfTimesteps = MaxIterations;
            c.ImmediatePlotPeriod = PlotInterval;
            if(dgDegree == 0 && numOfCellsX < 10 && numOfCellsT < 10) {
                c.SuperSampling = 4;
            } else {
                c.SuperSampling = 3;
            }

            #endregion

            #region SQP constants/parameters
            // ### Time config ###
            c.NoOfTimesteps = MaxIterations;
            #endregion

            #region LevelSetStuff
            c.IsTwoLevelSetRun = false;
            c.LevelSetDegree = lsDegree;
            c.OptiLevelSetType = OptiLevelSetType.SplineLevelSet;
            c.GetLevelSet = GetLevelSet.FromFunction;
            c.SpeciesTable[0, 0] = "L";
            c.SpeciesTable[1, 0] = "R";
            c.LsOne_NegSpecies = "L";
            c.LsOne_PosSpecies = "R";
            c.LsOne_SpeciesPairs = new string[,] { { "L", "R" } };
            c.SpeciesToEvaluate = new string[] { "L", "R" };

            if(dgDegree>0) {
                c.solRunType = SolverRunType.PContinuation;
                c.minimalSQPIterations = new int[] { 10, 10, 10, 10, 10,10,10 };
                c.reInitTols = new double[] { 0.0, -1, -1, -1, -1, -1, -1 };
            }
            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.AddVariable(XESFVariables.LevelSet, lsDegree);
            #endregion 

            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.ConvectiveInterfaceFlux_LsTwo = interfaceFluxLS1;

            c.FluxVersion = FluxVersion.NonOptimized;
            #endregion


            #region Grid
            //double xMin = 0.0;
            //double xMax = 3.0;
            double tMin = 0;
            double tMax = endTime;
            
            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] tNodes = GenericBlas.Linspace(tMin, tMax, numOfCellsT + 1);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, tNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.DefineEdgeTags(delegate (double[] X) {
                    if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                        return 2;
                    } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                        return 1;
                    } else if(Math.Abs(X[1] - tMax) < 1e-14) { // top boundary
                        return 2;
                    } else { //bottom boundary
                        return 1;
                    }
                });
                return grid;
            };
            #endregion

            #region Initial Values for Solution + Boundary Values
            //supersonic Base Flow - left values
            double densityL = 1*scaling; double pressureL = 1*scaling;
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double cL = Math.Sqrt(gamma * pressureL / densityL);
            double velocityL = MachL * cL;

            //get rigth values
            (double densityR, double velocityR, double pressureR, double cR, double MachR)
                = ComputeNormalShockWaveRelations(densityL, velocityL, pressureL, MachL, gamma);

            
            Func<double, double> f_waveform = x => 0; //default no pertubation
            //Waveform
            if(waveform == "1sinus") { //one period of sinus
                f_waveform = delegate (double X) {
                    if(wavePosition < X && X < wavePosition + waveLength) {
                        return Math.Sin(2 * Math.PI * X / waveLength);
                    } else {
                        return 0;
                    }
                };
            } else if(waveform == "halfsinus") { //half of period of sinus
                f_waveform = delegate (double X) {
                    if(wavePosition < X && X < wavePosition + waveLength) {
                        return Math.Sin(Math.PI * X / waveLength);
                    } else {
                        return 0;
                    }
                };
            } else if(waveform == "fullsinus") { //full sinus
                f_waveform = delegate (double X) {
                    if(shockPosition < wavePosition) {
                        if(wavePosition<X ) {
                            return Math.Sin(2 * Math.PI * X / waveLength);
                        } else {
                            return 0;
                        }
                    } else {
                        if(wavePosition > X) {
                            return Math.Sin(2 * Math.PI * X / waveLength);
                        } else {
                            return 0;
                        }
                    }
                };
            } else if(waveform == "bump") { //c ifnity bumb (differentiable everywhere)

                double L = waveLength / 2;
                double bumpPos = wavePosition + L;
                Func<double, double> f_base = x => Math.Exp(-1 / (1 - x * x)) * Math.E;

                f_waveform = delegate (double X) {
                    if(bumpPos - L < X && X < bumpPos + L) {
                        return f_base((X - bumpPos) / L);
                    } else {
                        return 0;
                    }
                };
            }

                // Functions for the perturbations
                Func<double[], double> pressure_per = delegate (double[] X) {
                if(X[0] > 0) {
                    if(X[0] < shockPosition) {
                        return p_amp_neg*pressureL * f_waveform(X[0] + (velocityL - cL) * X[1]) + p_amp_pos * pressureL * f_waveform(X[0] - (velocityL + cL) * X[1]);
                    } else {
                        return p_amp_neg * pressureL * f_waveform(X[0] + (velocityR - cR) * X[1]) + p_amp_pos * pressureL * f_waveform(X[0] - (velocityR + cR) * X[1]); ;
                    }
                } else {
                    return 0;
                }
                
            };
            Func<double[], double> density_per = delegate (double[] X) {
                if(X[0] > 0) {
                    if(X[0] < shockPosition) {
                    return pressure_per(X) / (cL * cL);
                } else {
                    return pressure_per(X) / (cR * cR);
                }
                } else {
                    return 0;
                }
            };
            Func<double[], double> velocity_per = delegate (double[] X) {
                if(X[0] > 0) {
                    if(X[0] < shockPosition) {
                    return -p_amp_neg * pressureL / (densityL * cL) * f_waveform(X[0] + (velocityL - cL) * X[1]) + p_amp_pos * pressureL / (densityL * cL) * f_waveform(X[0] - (velocityL + cL) * X[1]);
                } else {
                    return -p_amp_neg * pressureL / (densityR * cR) * f_waveform(X[0] + (velocityR - cR) * X[1]) + p_amp_pos * pressureL / (densityR * cR) * f_waveform(X[0] - (velocityR + cR) * X[1]); ;
                }
                } else {
                    return 0;
                }
            };

            c.LevelSetOneInitialValue = X => 200;
            if(withLevelSet) {
                c.LevelSetOneInitialValue = X => shockPosition;
                c.LevelSetTwoInitialValue = X => shockPosition;
                //c.LevelSetOneInitialValue = X => (velocityL - cL) * X[1];
            }
            //c.LevelSetTwoInitialValue = X => 200;

            //dirichlet boundary vals
            Func<double[], int, double> bnd_vals = delegate (double[] X, int component) {
                if(X[0] < shockPosition) {
                    if(component == 0) {
                        return densityL + density_per(X);
                    } else if(component == 1) {
                        return velocityL + velocity_per(X);
                    } else if(component == 2) {
                        return pressureL + pressure_per(X);
                    } else {
                        throw new NotSupportedException("component not supported");
                    }
                } else {
                    if(component == 0) {
                        return densityR + density_per(X);
                    } else if(component == 1) {
                        return velocityR + velocity_per(X);
                    } else if(component == 2) {
                        return pressureR + pressure_per(X);
                    } else {
                        throw new NotSupportedException("component not supported");
                    }
                }
                
            };

            c.SolDegree = dgDegree;
            c.GetInitialValue = GetInitialValue.FromFunctionPerSpecies;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;
            c.AddVariable(XESTSFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESTSFVariables.Pressure, dgDegree);
            c.AddVariable(XESTSFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESTSFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            if(withShock) {
                c.BaseFlowPerSpeciesAndField = new Dictionary<Tuple<string, string>, double> {
                    { new Tuple<string, string>("p", "L"), pressureL },
                    { new Tuple<string, string>("p", "R"), pressureR  },

                    { new Tuple<string, string>("rho", "L"), densityL  },
                    { new Tuple<string, string>("rho", "R"), densityR  },

                    { new Tuple<string, string>("u0", "L"), velocityL  },
                    { new Tuple<string, string>("u0", "R"), velocityR  },

                    { new Tuple<string, string>("m0", "L"), velocityL*densityL },
                    { new Tuple<string, string>("m0", "R"), velocityR*densityR }
                 };

                c.AddVariable(XESTSFVariables.PertubationPressure, dgDegree);
                c.AddVariable(XESTSFVariables.PertubationDensity, dgDegree);
                c.AddVariable(XESTSFVariables.PertubationVelocityX, dgDegree);
                c.AddVariable(XESTSFVariables.PertubationMomentumX, dgDegree);

            }

            // Boundary conditions in PRIMITIVE variables
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => bnd_vals(new double[] { X[0], t }, 0));
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => bnd_vals(new double[] { X[0], t }, 0));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => bnd_vals(new double[] { X[0], t }, 1));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => bnd_vals(new double[] { X[0], t }, 1));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => bnd_vals(new double[] { X[0], t }, 2));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => bnd_vals(new double[] { X[0], t }, 2));

            //// Initial conditions in PRIMITIVE variables 
            if(withLevelSet) {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", (X) => densityL);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", (X) => velocityL);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", (X) => pressureL);

                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", (X) => densityR);
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", (X) => velocityR);
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", (X) => pressureR);
            } else {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", (X) => bnd_vals(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", (X) => bnd_vals(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", (X) => bnd_vals(X, 2));

                c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", (X) => bnd_vals(X, 0));
                c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", (X) => bnd_vals(X, 1));
                c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", (X) => bnd_vals(X, 2));
            }
            #endregion

            
            #region Project Information
            c.ProjectName = "XESTSF_AcousticWave1D";
            c.SessionName = $"AW_p{dgDegree}_xCells{numOfCellsX}_tCells{numOfCellsT}_tMax{endTime}_{waveform}_Mach{MachL}_sP{shockPosition}_wP{wavePosition}_ampneg{p_amp_neg}_amppos{p_amp_pos}_wL{withLevelSet}_wS{withShock}_wL{waveLength}";

            #endregion

            return c;
        }
        /// <summary>
        /// Simple stationary shock wave at position 0.1 on a 9x10 grid with xMin=-1 and xMax=+1
        /// </summary>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsDegree"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="shocksetup"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="optiLSDegree"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <returns></returns>
        public static XESTSFControl StationaryShockWave(int MaxIterations = 100, int dgDegree = 0, int numOfCellsX = 9,
        int numOfCellsY = 10, int PlotInterval = 1,
        string dbPath = null, int lsDegree = 3, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.Roe,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.RoeInterface,
        GetLevelSet shocksetup = GetLevelSet.FromParams, OptiLevelSetType optiLevelSetType = OptiLevelSetType.SplineLevelSet, int optiLSDegree = 1,
        GetInitialValue initialValue = GetInitialValue.FromFunctionPerSpecies, bool iVFromShockRelations = true, double agg = 0.4, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch, MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit
        ) {
            var c = new XESTSFControl();


            //c.PartiallyFixLevelSetForSpaceTime = false;
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
            c.minimalSQPIterations = new int[] { 40, 20, 20 };
            c.reInitTols = new double[] { -10, -10, -10, -10, -10, -10, -10, -10 };
            #endregion

            #region LevelSetStuff

            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            c.IsTwoLevelSetRun = false;

            c.LevelSetDegree = lsDegree;
            c.AddVariable(XESFVariables.LevelSet, lsDegree);

            c.SpeciesTable[0, 0] = "L";
            c.SpeciesTable[1, 0] = "R";
            c.LsOne_NegSpecies = "L";
            c.LsOne_PosSpecies = "R";
            c.LsOne_SpeciesPairs = new string[,] { { "L", "R" } };
            c.SpeciesToEvaluate = new string[] { "L", "R" };
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = optiLSDegree;

            #endregion 

            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.FluxVersion = FluxVersion.NonOptimized;
            c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;
            #endregion

            #region Grid
            double xMin = -1.0;
            double xMax = 1.0;
            double tMin = 0.0;
            double tMax = 1.0;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] tNodes = GenericBlas.Linspace(tMin, tMax, numOfCellsY + 1);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, tNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");
                //grid.EdgeTagNames.Add(4, "SubsonicOutlet");
                grid.DefineEdgeTags(delegate (double[] X) {
                    if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                        return 2;
                    } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                        return 1;
                    } else if(Math.Abs(X[1] - tMax) < 1e-14) { // top boundary
                        return 2;
                    } else { //bottom boundary
                        return 1;
                    }
                });
                return grid;
            };
            #endregion

            #region Initial Position of LevelSet
            double x_s = 0.1;
            double Mach = 2.0;
            // ### straight evel set function ###
            c.GetLevelSet = GetLevelSet.FromFunction;
            switch(c.OptiLevelSetType) {
                case OptiLevelSetType.SplineLevelSet:
                c.LevelSetOneInitialValue = delegate (double[] X) { return X[1]*0.1+ x_s; };
                c.LevelSetTwoInitialValue = delegate (double[] X) { return X[1] * 0.1 + x_s; };
                break;
                default:
                c.LevelSetOneInitialValue = delegate (double[] X) { return x_s - X[0]; };
                c.LevelSetTwoInitialValue = delegate (double[] X) { return x_s - X[0] ; };
                break;
            }
            
            #endregion

            #region Inital Values for Solution + Boundary Values
            c.SolDegree = dgDegree;
            c.GetInitialValue = initialValue;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;
            c.ExactEnthalpy = 6.3;
            c.AddVariable(XESTSFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESTSFVariables.Pressure, dgDegree);
            c.AddVariable(XESTSFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESTSFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            // #########################
            // Parameters
            // #########################
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double Ms = 1.5;

            // #########################
            // Shock
            // #########################
            double densityLeft = 1;
            double densityRight = (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms) * densityLeft;
            double pressureLeft = 1;
            double pressureRight = (1 + 2 * gamma / (gamma + 1) * (Ms * Ms - 1)) * pressureLeft;
            double velocityXLeft = Ms * Math.Sqrt(gamma * pressureLeft / densityLeft);
            double velocityXRight = (2 + (gamma - 1) * Ms * Ms) / ((gamma + 1) * Ms * Ms) * velocityXLeft;    // (1)
            //double velocityXRight2 = velocityXLeft * densityLeft / densityRight; // equivalent to (1)
            //double MsPostShock = Math.Sqrt((1 + ((gamma - 1) / 2) * Ms * Ms) / (gamma * Ms * Ms - (gamma - 1) / 2));
            //double velocityXRight3 = MsPostShock * Math.Sqrt(gamma * pressureRight / densityRight);     // equivalent to (1)
            double velocityYLeft = 0.0;
            double velocityYRight = 0.0;

            // ### Initial conditions ###
            //c.uAInitPrimitive = new Func<double[], double>[] {
            //        X => densityLeft,
            //        X => velocityXLeft,
            //        X => velocityYLeft,
            //        X => pressureLeft
            //};
            //c.uBInitPrimitive = new Func<double[], double>[] {
            //        X => densityRight,
            //        X => velocityXRight,
            //        X => velocityYRight,
            //        X => pressureRight
            //};

            double momentumXLeft = densityLeft * velocityXLeft;
            double momentumXRight = densityRight * velocityXRight;
            double momentumYLeft = 0.0;
            double momentumYRight = 0.0;

            double innerEnergyLeft = pressureLeft / (gamma - 1);
            double innerEnergyRight = pressureRight / (gamma - 1);
            double kineticEnergyLeft = 0.5 * densityLeft * (velocityXLeft * velocityXLeft + velocityYLeft * velocityYLeft);
            double kineticEnergyRight = 0.5 * densityRight * (velocityXRight * velocityXRight + velocityYRight * velocityYRight);
            double totalEnergyLeft = innerEnergyLeft + kineticEnergyLeft;
            double totalEnergyRight = innerEnergyRight + kineticEnergyRight;

            // ### Initial smoothing ###
            double crossProduct2D(double[] a, double[] b) {
                return a[0] * b[1] - a[1] * b[0];
            }

            // Normal vector of initial shock
            Vector normalVector = new Vector(1, 0);

            // Direction vector of initial shock
            Vector r = new Vector(normalVector.y, -normalVector.x);
            r.NormalizeInPlace();

            //Distance from a point X to the initial shock
            double[] p = new double[] { x_s, 0.0 };

            double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
                double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
                double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

                return distance;
            }


            // Boundary conditions in PRIMITIVE variables
            
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => X[0]<x_s? densityLeft : densityRight);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => X[0] < x_s ? densityLeft : densityRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => X[0] < x_s ? velocityXLeft: velocityXRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => X[0] < x_s ? velocityXLeft : velocityXRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => X[0] < x_s ? pressureLeft:pressureRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => X[0] < x_s ? pressureLeft:pressureRight);


            double exact_sol(double[] X, int component) {
                //if(X[0] > 0.5 + (X[1] / LevelSet1Prime)) { //wedge -X[0] + 0.5 + (X[1] / LevelSet1Prime)
                //    return 0;
                if(X[0] >= x_s) {/*&& (X[0] <= 0.5 + (X[1] / LevelSet1Prime)))*/ //post Shock
                    if(component == 0) {
                        return densityRight;
                    } else if(component == 1) {
                        return velocityXRight;
                    } else if(component == 2) {
                        return velocityYRight;
                    } else if(component == 3) {
                        return pressureRight;
                    }
                } else { //pre Shock
                    if(component == 0) {
                        return densityLeft;
                    } else if(component == 1) {
                        return velocityXLeft;
                    } else if(component == 2) {
                        return 0;
                    } else if(component == 3) {
                        return pressureLeft;
                    }
                }
                return 0;
            }

            //// Initial conditions in PRIMITIVE variables -- Pre Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => exact_sol(X, 0));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L", X => exact_sol(X, 1));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => exact_sol(X, 3));

            //// Initial conditions in PRIMITIVE variables -- Post Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => exact_sol(X, 0));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => exact_sol(X, 1));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => exact_sol(X, 3));
            
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
            c.ProjectName = "XESTSF_Stationary Shock Wave";
            c.SessionName = string.Format("1DStatShockWave_p{0}_xCells{1}_yCells{2}_agg{3}_",
                dgDegree, numOfCellsX, numOfCellsY, agg);

            #endregion

            return c;
        }

        /// <summary>
        /// Stationary shock wave employing two level Sets (only for testing, if it still works with a dummy levelset)
        /// </summary>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsDegree"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="shocksetup"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="optiLSDegree"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <returns></returns>
        public static XESTSFControl StationaryShockWaveTwoLS(int MaxIterations = 100, int dgDegree = 0, int numOfCellsX = 9,
        int numOfCellsY = 10, int PlotInterval = 1,
        string dbPath = null, int lsDegree = 3, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.Roe,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.RoeInterface,
        GetLevelSet shocksetup = GetLevelSet.FromParams, OptiLevelSetType optiLevelSetType = OptiLevelSetType.SplineLevelSet, int optiLSDegree = 1,
        GetInitialValue initialValue = GetInitialValue.FromFunctionPerSpecies, bool iVFromShockRelations = true, double agg = 0.4, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch, MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit
        ) {
            var c = new XESTSFControl();



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
            c.minimalSQPIterations = new int[] { 40, 20, 20 };
            c.reInitTols = new double[] { -10, -10, -10, -10, -10, -10, -10, -10 };
            #endregion

            #region LevelSetStuff

            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            c.IsTwoLevelSetRun = true;

            c.LevelSetDegree = lsDegree;
            c.AddVariable(XESFVariables.LevelSet, 1);
            c.AddVariable(XESFVariables.LevelSetTwo, lsDegree); 

            string RALS = "A";
            string LALS = "B";
            string RARS = "L";
            string LARS = "R";
            c.SpeciesTable[0, 0] = LALS; 
            c.SpeciesTable[1, 0] = RALS; 
            c.SpeciesTable[0, 1] = LARS;
            c.SpeciesTable[1, 1] = RARS;
            c.LsOne_SpeciesPairs = new string[,] { { LALS, RALS }, { LARS,RARS} };
            c.LsTwo_SpeciesPairs = new string[,] { { LALS, LARS }, { RALS,RARS} };

            //c.SpeciesToEvaluate = new string[] { RARS, LALS, RALS, LARS};
            c.SpeciesToEvaluate = new string[] {  LALS,LARS,RARS,RALS };
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = optiLSDegree;

            #endregion 

            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.ConvectiveInterfaceFlux_LsTwo = interfaceFluxLS1;

            c.FluxVersion = FluxVersion.NonOptimized;
            c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;
            #endregion

            #region Grid
            double xMin = -1.0;
            double xMax = 1.0;
            double tMin = 0.0;
            double tMax = 1.0;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] tNodes = GenericBlas.Linspace(tMin, tMax, numOfCellsY + 1);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, tNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");
                //grid.EdgeTagNames.Add(4, "SubsonicOutlet");
                grid.DefineEdgeTags(delegate (double[] X) {
                    if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                        return 2;
                    } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                        return 1;
                    } else if(Math.Abs(X[1] - tMax) < 1e-14) { // top boundary
                        return 2;
                    } else { //bottom boundary
                        return 1;
                    }
                });
                return grid;
            };
            #endregion

            #region Initial Position of LevelSet
            double x_s = 0.0;
            double Mach = 2.0;
            // ### straight level set function ###
            c.GetLevelSet = GetLevelSet.FromFunction;
            switch(c.OptiLevelSetType) {
                case OptiLevelSetType.SplineLevelSet:
                //c.LevelSetOneInitalValue = delegate (double[] X) { return -(X[1] * 0.1+0.1); };
                c.LevelSetOneInitialValue = delegate (double[] X) { return -(x_s - X[0] + 1.5 * X[1]-0.15); };
                c.LevelSetTwoInitialValue = delegate (double[] X) { return x_s; };
                break;
                default:
                c.LevelSetOneInitialValue = delegate (double[] X) { return x_s - X[0]; };
                c.LevelSetTwoInitialValue = delegate (double[] X) { return x_s - X[0]; };
                break;
            }
           

            #endregion

            #region Inital Values for Solution + Boundary Values
            c.SolDegree = dgDegree;
            c.GetInitialValue = initialValue;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;
            c.ExactEnthalpy = 6.3;
            c.AddVariable(XESTSFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESTSFVariables.Pressure, dgDegree);
            c.AddVariable(XESTSFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESTSFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            // #########################
            // Shock
            // #########################
            //Base Flow
            double MachLeft = 1.5; double densityLeft = 1; double pressureLeft = 1;
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double cLeft = Math.Sqrt(gamma * pressureLeft / densityLeft);
            double velocityXLeft = ComputeVelocity(densityLeft, pressureLeft, MachLeft, gamma);
            (double densityRight, double velocityXRight, double pressureRight, double cRight, double MachRight)
                = ComputeNormalShockWaveRelations(densityLeft, velocityXLeft, pressureLeft, MachLeft, gamma);


            double velocityYLeft = 0.0;
            double velocityYRight = 0.0;
            //double velocityXRight2 = velocityXLeft * densityLeft / densityRight; // equivalent to (1)
            //double MsPostShock = Math.Sqrt((1 + ((gamma - 1) / 2) * Ms * Ms) / (gamma * Ms * Ms - (gamma - 1) / 2));
            //double velocityXRight3 = MsPostShock * Math.Sqrt(gamma * pressureRight / densityRight);     // equivalent to (1)

                        //Pertubation wave parameters
            double kRight = 2;
            double omega = (cLeft + velocityXLeft) / kRight;
            double amplitude_p_hat = 0.0;

            // ### Initial conditions ###
            //c.uAInitPrimitive = new Func<double[], double>[] {
            //        X => densityLeft,
            //        X => velocityXLeft,
            //        X => velocityYLeft,
            //        X => pressureLeft
            //};
            //c.uBInitPrimitive = new Func<double[], double>[] {
            //        X => densityRight,
            //        X => velocityXRight,
            //        X => velocityYRight,
            //        X => pressureRight
            //};

            double momentumXLeft = densityLeft * velocityXLeft;
            double momentumXRight = densityRight * velocityXRight;
            double momentumYLeft = 0.0;
            double momentumYRight = 0.0;

            double innerEnergyLeft = pressureLeft / (gamma - 1);
            double innerEnergyRight = pressureRight / (gamma - 1);
            double kineticEnergyLeft = 0.5 * densityLeft * (velocityXLeft * velocityXLeft + velocityYLeft * velocityYLeft);
            double kineticEnergyRight = 0.5 * densityRight * (velocityXRight * velocityXRight + velocityYRight * velocityYRight);
            double totalEnergyLeft = innerEnergyLeft + kineticEnergyLeft;
            double totalEnergyRight = innerEnergyRight + kineticEnergyRight;

            // ### Initial smoothing ###
            double crossProduct2D(double[] a, double[] b) {
                return a[0] * b[1] - a[1] * b[0];
            }

            // Normal vector of initial shock
            Vector normalVector = new Vector(1, 0);

            // Direction vector of initial shock
            Vector r = new Vector(normalVector.y, -normalVector.x);
            r.NormalizeInPlace();

            //Distance from a point X to the initial shock
            double[] p = new double[] { x_s, 0.0 };

            double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
                double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
                double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

                return distance;
            }


            // Boundary conditions in PRIMITIVE variables

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#"+ RALS, (X, t) => X[0] < x_s ? densityLeft : densityRight);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#" + RARS, (X, t) => X[0] < x_s ? densityLeft : densityRight);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#" + LALS, (X, t) => X[0] < x_s ? densityLeft : densityRight);
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#" + LARS, (X, t) => X[0] < x_s ? densityLeft : densityRight);

            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#" + RALS, (X, t) => X[0] < x_s ? velocityXLeft : velocityXRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#" + RARS, (X, t) => X[0] < x_s ? velocityXLeft : velocityXRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#" + LALS, (X, t) => X[0] < x_s ? velocityXLeft : velocityXRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#" + LARS, (X, t) => X[0] < x_s ? velocityXLeft : velocityXRight);
            
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#" + LALS, (X, t) => X[0] < x_s ? pressureLeft : pressureRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#" + LARS, (X, t) => X[0] < x_s ? pressureLeft : pressureRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#" + RALS, (X, t) => X[0] < x_s ? pressureLeft : pressureRight);
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#" + RARS, (X, t) => X[0] < x_s ? pressureLeft : pressureRight);

            //Func<double[], int, double> exact_sol = delegate (double[] X, int component) {
            //    //if(X[0] > 0.5 + (X[1] / LevelSet1Prime)) { //wedge -X[0] + 0.5 + (X[1] / LevelSet1Prime)
            //    //    return 0;
            //    if(X[0] >= x_s) {/*&& (X[0] <= 0.5 + (X[1] / LevelSet1Prime)))*/ //post Shock
            //        if(component == 0) {
            //            return densityRight;
            //        } else if(component == 1) {
            //            return velocityXRight;
            //        } else if(component == 2) {
            //            return velocityYRight;
            //        } else if(component == 3) {
            //            return pressureRight;
            //        }
            //    } else { //pre Shock
            //        if(component == 0) {
            //            return densityLeft;
            //        } else if(component == 1) {
            //            return velocityXLeft;
            //        } else if(component == 2) {
            //            return 0;
            //        } else if(component == 3) {
            //            return pressureLeft;
            //        }
            //    }
            //    return 0;
            //};
            //Func<double[], int, string, double> exact_sol_sd = delegate (double[] X, int component, string subdomain) {

            //    if(subdomain == LALS || subdomain == RALS) { // acoustic is passed & pre shock //LALS
            //        if(component == 0) {
            //            return densityLeft;
            //        } else if(component == 1) {
            //            return velocityXLeft;
            //        } else if(component == 2) {
            //            return pressureLeft;
            //        } else {
            //            throw new NotSupportedException("component: " + component + " not supported");
            //        }
            //    } else if(subdomain == RARS || subdomain == LARS) { // acoustic is not passed & post Shock  //RARS
            //        if(component == 0) {
            //            return densityRight;
            //        } else if(component == 1) {
            //            return velocityXRight;
            //        } else if(component == 2) {
            //            return pressureRight;
            //        } else {
            //            throw new NotSupportedException("component: " + component + " not supported");
            //        }
            //    } else {
            //        throw new NotSupportedException("subdomain: " + subdomain + " not supported");
            //    }
            //};

            Func<double[], int, string, double> exact_sol_sd = delegate (double[] X, int component, string subdomain) {

                if(subdomain == LALS) { // acoustic is passed & pre shock //LALS
                    if(component == 0) {
                        return densityLeft + amplitude_p_hat * Math.Cos(kRight * X[0] - (omega) * (X[1] + ((X[0]) / (cLeft + velocityXLeft)))) / (cLeft * cLeft);
                    } else if(component == 1) {
                        return velocityXLeft - kRight * amplitude_p_hat * Math.Sin(kRight * X[0] - (omega) * (X[1] + ((X[0]) / (cLeft + velocityXLeft)))) / (omega * densityLeft);
                    } else if(component == 2) {
                        return pressureLeft + amplitude_p_hat * Math.Cos(kRight * X[0] - (omega) * (X[1] + ((X[0]) / (cLeft + velocityXLeft))));
                    } else {
                        throw new NotSupportedException("component: " + component + " not supported");
                    }
                } else if(subdomain == RALS) { // acoustic is not passed & pre Shock //RALS
                    if(component == 0) {
                        return densityLeft;
                    } else if(component == 1) {
                        return velocityXLeft;
                    } else if(component == 2) {
                        return pressureLeft;
                    } else {
                        throw new NotSupportedException("component: " + component + " not supported");
                    }
                } else if(subdomain == RARS) { // acoustic is not passed & post Shock  //RARS
                    if(component == 0) {
                        return densityRight;
                    } else if(component == 1) {
                        return velocityXRight;
                    } else if(component == 2) {
                        return pressureRight;
                    } else {
                        throw new NotSupportedException("component: " + component + " not supported");
                    }
                } else if(subdomain == LARS) { // acoustic is passed & post Shock  //LARS
                    if(component == 0) {
                        return densityRight + amplitude_p_hat * Math.Cos((omega / (cRight + velocityXRight)) * X[0] - (omega) * (X[1] + (0.4 / (cLeft + velocityXLeft)) + ((X[0] - 0.4) / (cRight + velocityXRight)))) / (cRight * cRight);
                    } else if(component == 1) {
                        return velocityXRight - (omega / (cRight + velocityXRight)) * amplitude_p_hat * Math.Sin((omega / (cRight + velocityXRight)) * X[0] - (omega) * (X[1] + (0.4 / (cLeft + velocityXLeft)) + ((X[0] - 0.4) / (cRight + velocityXRight)))) / (omega * densityRight);
                    } else if(component == 2) {
                        return pressureRight + amplitude_p_hat * Math.Cos((omega / (cRight + velocityXRight)) * X[0] - (omega) * (X[1] + (0.4 / (cLeft + velocityXLeft)) + ((X[0] - 0.4) / (cRight + velocityXRight))));
                    } else {
                        throw new NotSupportedException("component: " + component + " not supported");
                    }
                } else {
                    throw new NotSupportedException("subdomain: " + subdomain + " not supported");
                }
            };

            //TODO: put the exact solution
            Func<double[], int, double> exact_sol = delegate (double[] X, int component) {

                double tstar = (cLeft + velocityXLeft) * X[1];
                if(X[1] < tstar && X[0] < x_s) { // acoustic is passed & pre shock //LALS
                    return exact_sol_sd(X, component, LALS);

                } else if(X[1] > tstar && X[0] < x_s) { // acoustic is not passed & pre Shock //RALS
                    return exact_sol_sd(X, component, RALS);

                } else if(X[1] > tstar && X[0] > x_s) { // acoustic is not passed & post Shock  //RARS
                    return exact_sol_sd(X, component, RARS);

                } else { // acoustic is passed & post Shock  //LARS
                    return exact_sol_sd(X, component, LARS);
                }
            };

            ////// Initial conditions in PRIMITIVE variables -- Pre Shock
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + RALS, X => exact_sol(X, 0));
            //c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + RALS, X => exact_sol(X, 1));
            //c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#"+ RALS, X => exact_sol(X, 3));

            ////// Initial conditions in PRIMITIVE variables -- Post Shock
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + RARS, X => exact_sol(X, 0));
            //c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + RARS, X => exact_sol(X, 1));
            //c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#" + RARS, X => exact_sol(X, 3));

            ////// Initial conditions in PRIMITIVE variables -- Pre Shock
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + LALS, X => exact_sol(X, 0));
            //c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + LALS, X => exact_sol(X, 1));
            //c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#" +LALS, X => exact_sol(X, 3));

            ////// Initial conditions in PRIMITIVE variables -- Post Shock
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + LARS, X => exact_sol(X, 0));
            //c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + LARS, X => exact_sol(X, 1));
            //c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#" + LARS, X => exact_sol(X, 3));

            //// Initial conditions in PRIMITIVE variables -- Pre Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + RALS, X => exact_sol_sd(X, 0, RALS));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + RALS, X => exact_sol_sd(X, 1, RALS));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#" + RALS, X => exact_sol_sd(X, 2, RALS));

            //// Initial conditions in PRIMITIVE variables -- Post Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + RARS, X => exact_sol_sd(X, 0, RARS));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + RARS, X => exact_sol_sd(X, 1, RARS));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#" + RARS, X => exact_sol_sd(X, 2, RARS));

            //// Initial conditions in PRIMITIVE variables -- Pre Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + LALS, X => exact_sol_sd(X, 0, LALS));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + LALS, X => exact_sol_sd(X, 1, LALS));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#" + LALS, X => exact_sol_sd(X, 2, LALS));

            //// Initial conditions in PRIMITIVE variables -- Post Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#" + LARS, X => exact_sol_sd(X, 0, LARS));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#" + LARS, X => exact_sol_sd(X, 1, LARS));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#" + LARS, X => exact_sol_sd(X, 2, LARS));


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
            c.ProjectName = "XESTSF_Stationary Shock Wave 2LS";
            c.SessionName = string.Format("1DStatShockWave_p{0}_xCells{1}_yCells{2}_agg{3}_",
                dgDegree, numOfCellsX, numOfCellsY, agg);

            #endregion

            return c;
        }

    
        //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
        /// WIP - Dipen - 1D sod shock tube :-
        /// governing equations are euler equations & ideal gas law
        /// initially two regions with different density & pressure (left side high pressure region & right side low pressure region) seperated by a diaphragm
        /// at t=0 diaphragm rupters and create three waves
        /// 1. shock wave which, travells from left to right, starting from the position of the diaphragm
        /// 2. rarefaction wave, which travells from right to left, starting from position of the diaphragm
        /// 3. contact discontinuity, which moves from left to right, with less speed than the shock & chasing it.
        /// two levelsets, one for shock wave and one for contact discontinuity

        public static XESTSFControl XDG_1D_sod_shock_tube_TwoLs_Base(int MaxIterations = 100, int dgDegree = 3, int numOfCellsX = 50,
        int numOfCellsY = 1, int numOfCellsT = 50, int PlotInterval = 1,
        string dbPath = null, int lsDegree = 1, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.Roe,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.RoeInterface,
        ConvectiveInterfaceFluxes interfaceFluxLS2 = ConvectiveInterfaceFluxes.RoeInterface,
        GetLevelSet shocksetup = GetLevelSet.FromParams, OptiLevelSetType optiLevelSetType = OptiLevelSetType.SplineLevelSet, int optiLSDegree = 1,
        GetInitialValue initialValue = GetInitialValue.FromFunctionPerSpecies, double agg = 0.4, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch, MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit
        ) {
            var c = new XESTSFControl();

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

            // Solver Type
            c.solRunType = SolverRunType.PContinuation;
            c.DgDegree_Start = 0;
            c.SolDegree = dgDegree;
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
            c.minimalSQPIterations = new int[] { 20, 20, 20, 20, 20, 20 };
            c.reInitTols = new double[] { 0, -0.6, -1, -2, -2, -2 };
            c.ReiniTMaxIters = new int[] { 0, 15, 15, 15, 15, 15 };

            c.fphiType = FphiType.CurvatureAll;
            #endregion

            #region LevelSetStuff

            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            c.IsTwoLevelSetRun = true;

            c.LevelSetDegree = lsDegree;
            c.AddVariable(XESFVariables.LevelSet, 1); // this is the levelset for the contact discontinuity
            c.AddVariable(XESFVariables.LevelSetTwo, lsDegree); // this is the levelSet used for the shock wave

            ///Left side of the Contact Discontinuity the are rarefaction waves also....???????
            //c.SpeciesTable[0, 0] = "LCLS"; // left side of both the Contact discontinuity && the shock wave
            //c.SpeciesTable[1, 0] = "RCLS"; // right side of the contact discontinuity && left side of the shock wave(post shock region)
            //c.SpeciesTable[1, 1] = "RCRS"; // right side of both the contact discontinuity && the shock wave
            //c.SpeciesTable[0, 1] = "LCRS"; // left side of the contact discontinuity && right side of the shock (empty)
            //c.LsOne_NegSpecies = "LCLS";
            //c.LsOne_PosSpecies = "RCLS";
            //c.LsTwo_NegSpecies = "RCLS";
            //c.LsTwo_PosSpecies = "RCRS";
            //c.SpeciesToEvaluate = new string[] { "LCLS", "RCLS", "RCRS", "LCRS" };
            string RCLS = "A";
            string LCLS = "B";
            string RCRS = "L";
            string LCRS = "R";
            c.SpeciesTable[0, 0] = LCLS;
            c.SpeciesTable[1, 0] = RCLS;
            c.SpeciesTable[0, 1] = LCRS;
            c.SpeciesTable[1, 1] = RCRS;
            c.LsOne_SpeciesPairs = new string[,] { { LCLS, RCLS }, { LCRS, RCRS } };
            c.LsTwo_SpeciesPairs = new string[,] { { LCLS, LCRS }, { RCLS, RCRS } };

            c.LevelSetTwoDegree = lsDegree;
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = optiLSDegree;
            #endregion

            //shockLevelSet_Db = @"C:\BoSSS-experimental\internal\src\private-mag\XESF\Tests\bosss_db_levelSets.zip";


            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.ConvectiveInterfaceFlux_LsTwo = interfaceFluxLS2;
            c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;
            #endregion

            #region Grid
            double xMin = 0;
            double xMax = 1;
            double tMin = 0;
            double tMax = 0.2;

            //todo: check bc
            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] tNodes = GenericBlas.Linspace(tMin, tMax, numOfCellsT + 1);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, tNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SubsonicInlet");// in shock tube left and right side are always at rest untill the waves comes to the endwall and reflect
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");//so left and right are taken as subsonic inlet
                //Also at the top boundary flows past the waves will be supersonic so supersonic outlet
                //And at bottom boundary (at t=0) there is no flow/waves, so subsonic inlet
                grid.DefineEdgeTags(delegate (double[] X) {
                    if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                        return 1;
                    } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                        return 1;
                    } else if(Math.Abs(X[1] - tMax) < 1e-14) { // top boundary
                        return 2;
                    } else { //bottom boundary
                        return 1;
                    }
                });
                return grid;
            };

            #endregion

            #region Initial Position of LevelSets

            ///velocity of the expansion fan-head = -1.18322 m/s
            ///velocity of the expansion fan-tail = -0.07027 m/s
            ///velocity of the contact discontinuity = 0.92745 m/s ( same for the gases both side of the contact wave)
            /// " " " pressure = 0.30313 Pa " " "
            /// left density (region 3) = 0.42632
            /// right density (region 2) = 0.26557
            /// velocity of the shock wave = 1.75216 m/s

            //TODO: set right functions and describe

            //// level set 1

            c.LevelSetOneInitialValue = delegate (double[] X) { return X[0] - 0.5 - 0.92745 * X[1]; };

            //// level set 2
            c.LevelSetTwoInitialValue = X => X[0] - 0.5 - 1.75216 * X[1];
            #endregion

            #region Inital Values for Solution + Boundary Values
            c.SolDegree = dgDegree;
            c.GetInitialValue = initialValue;
            c.EquationOfState = IdealGas.Air;
            //c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            //c.ReynoldsNumber = 1.0;
            //c.PrandtlNumber = 0.71;

            c.AddVariable(XESTSFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESTSFVariables.Pressure, dgDegree);
            c.AddVariable(XESTSFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESTSFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            //TODO: rescale
            // Boundary conditions
            double density_4 = 1; // region 4 is the initial left side of the diaphrgm / high pressure side
            double pressure_4 = 1;
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double speedOfSound = Math.Sqrt(gamma * pressure_4 / density_4);
            double velocity_4 = 0;
            double momentumX_4 = density_4 * velocity_4;
            double innerEnergy_4 = pressure_4 / (gamma - 1);
            double kineticEnergy_4 = 0.5 * density_4 * velocity_4 * velocity_4;
            double totalEnergy_4 = innerEnergy_4 + kineticEnergy_4;

            double density_1 = 0.125; // region 1 is the initial right side of the diaphragm / low pressure side
            double pressure_1 = 0.1;
            double velocity_1 = 0;

            double density_3 = 0.42632; // region 3 is between the tail of expansion fan and the contact discontinuity
            double pressure_3 = 0.30313;
            double velocity_3 = 0.92745;

            double density_2 = 0.26557; // region 2 is between the shock wave and the contact discontinuity
            double pressure_2 = 0.30313;
            double velocity_2 = 0.92745;

            // Boundary conditions in PRIMITIVE variables
            //TODO: adjust for new sub-domain names and correct values
            c.AddBoundaryValue("AdiabaticSlipWall", CompressibleVariables.Density + "#LCLS", (X, t) => exact_sol(new double[] { X[0], t }, 0));
            c.AddBoundaryValue("AdiabaticSlipWall", CompressibleVariables.Density + "#RCLS", (X, t) => exact_sol(new double[] { X[0], t }, 0));
            c.AddBoundaryValue("AdiabaticSlipWall", CompressibleVariables.Density + "#RCRS", (X, t) => exact_sol(new double[] { X[0], t }, 0));
            c.AddBoundaryValue("AdiabaticSlipWall", XESFVariables.Velocity.xComponent + "#LCLS", (X, t) => exact_sol(new double[] { X[0], t }, 1));
            c.AddBoundaryValue("AdiabaticSlipWall", XESFVariables.Velocity.xComponent + "#RCLS", (X, t) => exact_sol(new double[] { X[0], t }, 1));
            c.AddBoundaryValue("AdiabaticSlipWall", XESFVariables.Velocity.xComponent + "#RCRS", (X, t) => exact_sol(new double[] { X[0], t }, 1));
            c.AddBoundaryValue("AdiabaticSlipWall", XESFVariables.Pressure + "#LCLS", (X, t) => exact_sol(new double[] { X[0], t }, 2));
            c.AddBoundaryValue("AdiabaticSlipWall", XESFVariables.Pressure + "#RCLS", (X, t) => exact_sol(new double[] { X[0], t }, 2));
            c.AddBoundaryValue("AdiabaticSlipWall", XESFVariables.Pressure + "#RCRS", (X, t) => exact_sol(new double[] { X[0], t }, 2));


            //   double density_Right = (gamma + 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) / ((gamma - 1.0) * Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) + 2.0);
            //   double pressure_Right = 1.0 + (2.0 * gamma * (Mach * Mach * Math.Sin(shock_angle_radial) * Math.Sin(shock_angle_radial) - 1.0)) / (gamma + 1.0);

            // calculate Velocity via Transformation matrix

            double innerEnergy_1 = pressure_1 / (gamma - 1);
            double kineticEnergy_1 = 0.5 * density_1 * (velocity_1 * velocity_1);
            double totalEnergy_1 = innerEnergy_1 + kineticEnergy_1;


            //TODO: put the exact solution

            double exact_sol(double[] X, int component) {

                if(X[0] < 0.5 - 1.18322 * X[1]) { // region 4
                    exact_sol_sd(X, component, "LCLS");

                } else if(X[0] > 0.5 - 0.07027 * X[1] && X[0] < 0.5 + 0.92745 * X[1]) { // region 3
                    exact_sol_sd(X, component, "LCLS");

                } else if(X[0] > 0.5 + 0.92745 * X[1] && X[0] < 0.5 + 1.75216 * X[1]) { // region 2
                    exact_sol_sd(X, component, "RCLS");

                } else if(X[0] > 0.5 + 1.75216 * X[1]) { // region 1
                    exact_sol_sd(X, component, "RCRS");
                }

                return 0;
            }
            double exact_sol_sd(double[] X, int component, string subdomain) {

                if(subdomain == "LCLS" && X[0] < 0.5 - 1.18322 * X[1]) { // region 4
                    if(component == 0) {
                        return density_4;
                    } else if(component == 1) {
                        return velocity_4;
                    } else if(component == 2) {
                        return pressure_4;
                    }
                } else if(subdomain == "LCLS" && X[0] > 0.5 - 0.07027 * X[1] && X[0] < 0.5 + 0.92745 * X[1]) { // region 3
                    if(component == 0) {
                        return density_3;
                    } else if(component == 1) {
                        return velocity_3;
                    } else if(component == 2) {
                        return pressure_3;
                    }

                } else if(subdomain == "RCLS") { // region 2
                    if(component == 0) {
                        return density_2;
                    } else if(component == 1) {
                        return velocity_2;
                    } else if(component == 2) {
                        return pressure_2;
                    }

                } else if(subdomain == "RCRS") { // region 1
                    if(component == 0) {
                        return density_1;
                    } else if(component == 1) {
                        return velocity_1;
                    } else if(component == 2) {
                        return pressure_1;
                    }
                }

                return 0;
            }
            //TODO: adjust for new sub-domain names and add exact solution

            //// Initial conditions in PRIMITIVE variables -- LCLS
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#LCLS", X => exact_sol_sd(X, 0, "LCLS"));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#LCLS", X => exact_sol_sd(X, 1, "LCLS"));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#LCLS", X => exact_sol_sd(X, 2, "LCLS"));

            //// Initial conditions in PRIMITIVE variables -- RCLS
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#RCLS", X => exact_sol_sd(X, 0, "RCLS"));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#RCLS", X => exact_sol_sd(X, 1, "RCLS"));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#RCLS", X => exact_sol_sd(X, 2, "RCLS"));

            //// Initial conditions in PRIMITIVE variables -- RCRS
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#RCRS", X => exact_sol_sd(X, 0, "RCRS"));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#RCRS", X => exact_sol_sd(X, 1, "RCRS"));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#RCRS", X => exact_sol_sd(X, 2, "RCRS"));

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
            c.ProjectName = "XESTSTSF_SodShockTube_TwoLs";
            c.SessionName = string.Format("XESTSTSF_SodShockTube_TwoLs_p{0}_xCells{1}_yCells{2}_{3}_iterations{4}_angle{5}_lSDegree{6}_Glob{7}_agg{8}_",
                dgDegree, numOfCellsX, numOfCellsY, interfaceFluxLS2.ToString(), c.NoOfTimesteps, optiLSDegree, globalization.ToString(), agg);

            #endregion

            return c;

        }


        /// <summary>
        /// Shu Osher problem.
        /// </summary>
        /// <param name="DgDegree_start"></param>
        /// <param name="endTime"></param>
        /// <param name="MaxIterations"></param>
        /// <param name="dgDegree"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsT"></param>
        /// <param name="PlotInterval"></param>
        /// <param name="withDirichletBoundaryFunc"></param>
        /// <param name="dbPath"></param>
        /// <param name="lsDegree"></param>
        /// <param name="bulkFlux"></param>
        /// <param name="interfaceFluxLS1"></param>
        /// <param name="optiLevelSetType"></param>
        /// <param name="optiLSDegree"></param>
        /// <param name="initialValue"></param>
        /// <param name="iVFromShockRelations"></param>
        /// <param name="agg"></param>
        /// <param name="globalization"></param>
        /// <param name="meritFunctionType"></param>
        /// <returns></returns>
        public static XESTSFControl ShuOsher1D(int DgDegree_start=1,double endTime=1.1, int MaxIterations = 100, int dgDegree = 3, int numOfCellsX = 9,
        int numOfCellsT = 10, int PlotInterval = -1, bool withDirichletBoundaryFunc=false,
        string dbPath = null, int lsDegree = 3, ConvectiveBulkFluxes bulkFlux = ConvectiveBulkFluxes.Roe,
        ConvectiveInterfaceFluxes interfaceFluxLS1 = ConvectiveInterfaceFluxes.RoeInterface,
        OptiLevelSetType optiLevelSetType = OptiLevelSetType.SplineLevelSet, int optiLSDegree = 1,
        GetInitialValue initialValue = GetInitialValue.FromFunctionPerSpecies, bool iVFromShockRelations = true, double agg = 0.4, GlobalizationStrategy globalization = ApplicationWithIDT.GlobalizationStrategy.LineSearch, MeritFunctionType meritFunctionType = MeritFunctionType.ExactMerit
        ) {
            var c = new XESTSFControl();


            //c.PartiallyFixLevelSetForSpaceTime = false;
            #region DBStuff and Plotting
            c.DbPath = dbPath;
            c.savetodb = true;
            if(dbPath == null)
                c.savetodb = false;
            c.NoOfTimesteps = MaxIterations;
            c.ImmediatePlotPeriod = PlotInterval;
            if(dgDegree == 0 && numOfCellsX < 10 && numOfCellsT < 10) {
                c.SuperSampling = 4;
            } else {
                c.SuperSampling = 3;
            }
            c.solRunType = SolverRunType.PContinuation;
            c.DgDegree_Start = DgDegree_start;
            c.GetInitialValue = GetInitialValue.FromFunctionPerSpecies;
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
            c.minimalSQPIterations = new int[] { 40, 20, 20 , 20, 20 };
            c.reInitTols = new double[] { -0.5, -0.5, -0.5, -0.5, -0.5, -10 };
            #endregion

            #region LevelSetStuff

            // ### Agglomeration and quadrature ###
            c.AgglomerationThreshold = agg;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            c.IsTwoLevelSetRun = false;

            c.LevelSetDegree = lsDegree;
            c.AddVariable(XESFVariables.LevelSet, lsDegree);

            c.SpeciesTable[0, 0] = "L";
            c.SpeciesTable[1, 0] = "R";
            c.LsOne_NegSpecies = "L";
            c.LsOne_PosSpecies = "R";
            c.LsOne_SpeciesPairs = new string[,] { { "L", "R" } };
            c.SpeciesToEvaluate = new string[] { "L", "R" };
            c.OptiLevelSetType = optiLevelSetType;
            c.OptiLevelSetDegree = optiLSDegree;

            #endregion 

            #region Fluxes
            //// ### Fluxes ###
            c.ConvectiveBulkFlux = bulkFlux;
            c.ConvectiveInterfaceFlux_LsOne = interfaceFluxLS1;
            c.FluxVersion = FluxVersion.NonOptimized;
            c.MassMatrixShapeandDependence = BoSSS.Solution.XdgTimestepping.MassMatrixShapeandDependence.IsNonIdentity;
            #endregion

            #region Grid
            double xMin = -5.0;
            double xMax = 2.0;
            double tMin = 0.0;
            double tMax = endTime;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] tNodes = GenericBlas.Linspace(tMin, tMax, numOfCellsT + 1);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, tNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");
                //grid.EdgeTagNames.Add(4, "SubsonicOutlet");
                grid.DefineEdgeTags(delegate (double[] X) {
                    if(Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                        return 2;
                    } else if(Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                        return 1;
                    } else if(Math.Abs(X[1] - tMax) < 1e-14) { // top boundary
                        return 2;
                    } else { //bottom boundary
                        return 1;
                    }
                });
                return grid;
            };
            #endregion

            #region Initial Position of LevelSet
            Func<double,double> x_s = t => -4.0 + t*3.5; //approx. shock position
            // ### straight evel set function ###
            c.GetLevelSet = GetLevelSet.FromFunction;
            switch(c.OptiLevelSetType) {
                case OptiLevelSetType.SplineLevelSet:
                c.LevelSetOneInitialValue = delegate (double[] X) { return x_s(X[1]); };
                break;
                default:
                c.LevelSetOneInitialValue = delegate (double[] X) { return x_s(X[1]) - X[0]; };
                break;
            }

            #endregion

            #region Inital Values for Solution + Boundary Values
            c.SolDegree = dgDegree;
            //c.GetInitialValue = initialValue;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;
            c.AddVariable(XESTSFVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(XESTSFVariables.Pressure, dgDegree);
            c.AddVariable(XESTSFVariables.Enthalpy, dgDegree);
            c.AddVariable(XESTSFVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            // #########################
            // Parameters
            // #########################


            Func<double, double> densityLeft = x => 3.857143;
            Func<double, double> velocityXLeft = x => 2.629369;
            Func<double, double> pressureLeft = x => 10.3333;
            Func<double, double> densityRight = x => 1+0.2*Math.Sin(5*x);
            Func<double, double> velocityXRight = x => 0;
            Func<double, double> pressureRight = x => 1;



            if(withDirichletBoundaryFunc) {
                c.SpeciesToEvaluate = new string[] { "L" }; // only evaluate left side of shock
                double gamma = 1.4;
                Func<double, double> MomRight = x => densityRight(x) * velocityXRight(x);
                Func<double, double> EnergyRight = x => pressureRight(x) / (gamma - 1) + 0.5 * densityRight(x) * velocityXRight(x) * velocityXRight(x);
                c.hasDirichletBoundary= true;
                c.DirichletBoundaryFunc = X => new double[] { densityRight(X[0]), MomRight(X[0]), EnergyRight(X[0]) };
            }
            // Boundary conditions in PRIMITIVE variables

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#L", (X, t) => X[0] < x_s(0) ? densityLeft(X[0]) : densityRight(X[0]));
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density + "#R", (X, t) => X[0] < x_s(0) ? densityLeft(X[0]) : densityRight(X[0]));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#L", (X, t) => X[0] < x_s(0) ? velocityXLeft(X[0]) : velocityXRight(X[0]));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Velocity.xComponent + "#R", (X, t) => X[0] < x_s(0) ? velocityXLeft(X[0]) : velocityXRight(X[0]));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#L", (X, t) => X[0] < x_s(0) ? pressureLeft(X[0]) : pressureRight(X[0]));
            c.AddBoundaryValue("SupersonicInlet", XESFVariables.Pressure + "#R", (X, t) => X[0] < x_s(0) ? pressureLeft(X[0]) : pressureRight(X[0]));


            //// Initial conditions in PRIMITIVE variables -- Pre Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#L", X => densityLeft(X[0]));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#L",X=>velocityXLeft(X[0]));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#L", X => pressureLeft(X[0]));

            //// Initial conditions in PRIMITIVE variables -- Post Shock
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density + "#R", X => densityRight(X[0]));
            c.InitialValues_Evaluators.Add(XESFVariables.Velocity.xComponent + "#R", X => velocityXRight(X[0]));
            c.InitialValues_Evaluators.Add(XESFVariables.Pressure + "#R", X => pressureRight(X[0]));

            #endregion


            #region Project Information
            c.ProjectName = "XESTSF_ShuOsher";
            c.SessionName = string.Format($"XESTSF_1DShuOsher_pStart{DgDegree_start}_pEnd{dgDegree}_xCells{numOfCellsX}_tCells{numOfCellsT}_agg{agg}_tMax{endTime}{(withDirichletBoundaryFunc? "_withDBF":"")}");

            #endregion

            return c;
        }
    }
}

