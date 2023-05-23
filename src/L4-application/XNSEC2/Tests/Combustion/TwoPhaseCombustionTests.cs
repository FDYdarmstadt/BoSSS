using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

//using BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases;
//using log4net.Filter;

namespace BoSSS.Application.XNSEC {

    static public partial class FullNSEControlExamples {


        /// <summary>
        /// A cold  droplet in a non-gravity medium, surrounded by air.
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSEC_Control ColdDroplet(int p = 1, int kelem = 16, int AMRlvl = 4) {

            XNSEC_Control C = new XNSEC_Control();

            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;
            bool steadyInterface = false;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            string _DbPath = @"C:\Databases\BoSSS_DB";
            //string _DbPath = null; // @"\\HPCCLUSTER\hpccluster-scratch\smuda\XNSE_studyDB";
            //string _DbPath = @"\\terminal03\Users\smuda\local\terminal03_XNSE_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "TwoPhaseComb";
            C.rhoOne = true;
            C.EnableMassFractions = true;
            C.EnableTemperature = true;
            C.NumberOfChemicalSpecies = 2;
            C.ChemicalReactionActive = false;

            //C.ProjectDescription = "Static droplet";
            //C.SessionName = "SD_meshStudy_Hysing_mesh" + kelem; // "_AMR"+AMRlvl;

            C.ContinueOnIoError = false;
            //C.LogValues = XNSE_Control.LoggingValues.Dropletlike;
            //C.PostprocessingModules.Add(new Dropletlike() { LogPeriod = 10 });

            #endregion


            // DG degrees
            // ==========
            #region degrees
            C.SetDGdegree(p);

         

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            //C.Tags.Add("Hysing");
            C.Tags.Add("La = 5000");
            C.PhysicalParameters.rho_A = 1e4;
            C.PhysicalParameters.rho_B = 1e4;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            double sigma = 1.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.Tags.Add("La = 0.005");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 10;
            //double sigma = 1;
            //C.PhysicalParameters.Sigma = sigma;

            //Air - Water(lenght scale == centimeters, 3D space)
            //C.PhysicalParameters.rho_A = 1e3;      // kg / cm^3
            //C.PhysicalParameters.rho_B = 1.2;    // kg / cm^3
            //C.PhysicalParameters.mu_A = 1e-3;       // kg / cm * sec
            //C.PhysicalParameters.mu_B = 17.1e-6;    // kg / cm * sec
            //double sigma = 72.75e-3;                // kg / sec^2 
            //C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double Lscl = 1.0;
            double xSize = Lscl * 1.0;
            double ySize = Lscl * 1.0;
            double zSize = Lscl * 1.0;

            if(D == 2) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize / 2.0, xSize / 2.0, kelem + 0);
                    double[] Ynodes = GenericBlas.Linspace(-ySize / 2.0, ySize / 2.0, kelem + 0);

                    double[] CutOut1Point1 = new double[2] { 0.001, 0.001 };
                    double[] CutOut1Point2 = new double[2] { -0.001, -0.001 };
                    

                    var CutOut1 = new BoSSS.Platform.Utils.Geom.BoundingBox(2); 
                    CutOut1.AddPoint(CutOut1Point1);
                    CutOut1.AddPoint(CutOut1Point2);

                    var CutOuts = new BoSSS.Platform.Utils.Geom.BoundingBox[] { CutOut1 };
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes,CutOuts: CutOuts);
                    //var grd = Grid2D.UnstructuredTriangleGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "Velocity_Inlet_FarAwayCondition");
                    grd.EdgeTagNames.Add(2, "Velocity_Inlet_Center");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(
                        Math.Abs(X[1] + ySize / 2.0) <= 1.0e-8||
                        Math.Abs(X[1] - ySize / 2.0) <= 1.0e-8 ||
                        Math.Abs(X[0] + xSize / 2.0) <= 1.0e-8 ||                            
                        Math.Abs(X[0] - xSize / 2.0) <= 1.0e-8)
                            et = 1;
                        else
                            et = 2;
                        return et;
                    });
                    var gDat = new GridData(grd);
                    var em1 = gDat.GetBoundaryEdges();
                    em1.SaveToTextFile("alledges2.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);
                    return grd;
                };
            }


            #endregion


            // Initial Values
            // ==============
            #region init

            double r = Lscl * 0.1;

            Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()).Sqrt() - r);         // signed distance
            //Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()) - r.Pow2());         // quadratic

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);

            C.InitialValues_Evaluators.Add("Temperature#A", X => 1.0);
            C.InitialValues_Evaluators.Add("Temperature#B", X => 1.0);

            C.InitialValues_Evaluators.Add("MassFraction0#A", X => 1.0);
            C.InitialValues_Evaluators.Add("MassFraction0#B", X => 0.0);
            C.InitialValues_Evaluators.Add("MassFraction1#A", X => 0.0);
            C.InitialValues_Evaluators.Add("MassFraction1#B", X => 1.0);

            C.InitialValues_Evaluators.Add("Pressure#A", X => 0.0);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);


            #endregion


            // exact solution
            // ==============
            #region exact

            C.Phi = ((X, t) => PhiFunc(X));

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            if(D == 2) {
                C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
                C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            }



            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", (X, t) => 0.0);
            C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("Velocity_Inlet_FarAwayCondition", VariableNames.MassFraction0, (X, t) => 0.00);
            C.AddBoundaryValue("Velocity_Inlet_FarAwayCondition", VariableNames.MassFraction1, (X, t) => 0.23);
            C.AddBoundaryValue("Velocity_Inlet_FarAwayCondition", VariableNames.Temperature, (X, t) => 1.0);


            C.AddBoundaryValue("Velocity_Inlet_Center", VariableNames.MassFraction0, (X, t) => 1.00);
            C.AddBoundaryValue("Velocity_Inlet_Center", VariableNames.MassFraction1, (X, t) => 0.00);
            C.AddBoundaryValue("Velocity_Inlet_Center", VariableNames.Temperature, (X, t) => 1.0);


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;
            C.solveKineticEnergyEquation = false;

            C.CheckJumpConditions = false;
            C.CheckInterfaceProps = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.0;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            if(AMRlvl > 0) {
                C.AdaptiveMeshRefinement = true;
                C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
                C.BaseRefinementLevel = AMRlvl;
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = AMRlvl });
            }
     

            C.BaseRefinementLevel =2;
            C.RefinementLevel = 1;
            C.AMR_startUpSweeps = 2;

            //C.InitSignedDistance = false;
            C.adaptiveReInit = false;

            C.LinearSolver = LinearSolverCode.automatic.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;


            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-7;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            //C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };


            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.Option_LevelSetEvolution = (steadyInterface) ? LevelSetEvolution.None : LevelSetEvolution.Fourier;
            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = 1.0 * sigma;
            //C.PhysicalParameters.lambda_I = 2.0 * sigma;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;



            if(C.Option_LevelSetEvolution == LevelSetEvolution.Fourier) {

                int numSp = 1800;
                double[] FourierP = new double[numSp];
                double[] samplP = new double[numSp];
                for(int sp = 0; sp < numSp; sp++) {
                    FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                    samplP[sp] = r;
                }

                C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
                    center = new double[] { 0.0, 0.0 },
                    FourierEvolve = Fourier_Evolution.MaterialPoints,
                    centerMove = CenterMovement.Reconstructed,
                };
            }


            #endregion


            // Timestepping
            // ============
            #region time

            //switch(p) {
            //    case 1: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //            break;
            //        }
            //    case 2: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //            C.Timestepper_BDFinit = TimeStepperInit.MultiInit;
            //            break;
            //        }
            //    default:
            //        C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //        C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //        break;
            //}

            //if(D == 3) {
            //    C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //    C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //}

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;


            //C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            //C.LSunderrelax = 0.05;
            C.Timestepper_LevelSetHandling = (steadyInterface) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = compMode;
            //C.CompMode = AppControl._CompMode.Transient; 

            double dt = 0.1; //0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1; // 12500; // (int)(125.0 / dt);
            C.saveperiod = 10;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            #endregion


            return C;

        }
       
    
    }
}