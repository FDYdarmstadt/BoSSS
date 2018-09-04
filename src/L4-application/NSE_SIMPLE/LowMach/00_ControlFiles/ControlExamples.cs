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

using ilPSP;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Queries;
using ilPSP.LinSolvers.PARDISO;
using System;
using static System.Math;

namespace NSE_SIMPLE.LowMach {

    public static class ControlExamples {

        public static LowMachSIMPLEControl SteadyCouetteFlowWithTemperatureGradient() {
            LowMachSIMPLEControl c = new LowMachSIMPLEControl();

            c.DbPath = @"..\..\Base\06_ZipDatabases\NUnitTests.zip";
            c.savetodb = false;

            c.GridGuid = new Guid("b3eb0eac-d1a1-440c-9f08-5dae1284607d");
            c.GridPartType = GridPartType.METIS;

            c.ProjectName = "Couette with temperature gradient";
            c.ProjectDescription = "Steady Low Mach SIMPLE";

            // Required fields
            c.FieldOptions.Add(
                VariableNames.VelocityX,
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                VariableNames.VelocityY,
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                VariableNames.Pressure,
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                VariableNames.Temperature,
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                VariableNames.ThermodynamicPressure,
                new FieldOpts() { Degree = 0, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            c.FieldOptions.Add(
                "Density",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                "Eta",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            // Auxiliary fields
            c.FieldOptions.Add(
                "DivB4",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "DivAfter",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorTest_x",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorTest_y",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorAna_x",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorAna_y",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });

            c.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => X[1]);
            c.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
            c.InitialValues_Evaluators.Add(VariableNames.Pressure, X => -0.97809076838538383654 * Log(1.2 * X[1] + 0.4) - 0.06641065188714375468);
            c.InitialValues_Evaluators.Add(VariableNames.Temperature, X => (1.6 - 0.4) * X[1] + 0.4);

            c.Algorithm = SolutionAlgorithms.Steady_SIMPLE;
            c.NoOfTimesteps = 1;
            c.L2NormPressureCorrection = 1.0e-6;
            c.L2NormVelocityResidual = 1.0e-6;
            c.L2NormTemperatureResidual = 1e-6;

            c.PredictorSolverFactory = () => new PARDISOSolver();
            c.CorrectorSolverFactory = () => new PARDISOSolver();
            c.TemperatureSolverFactory = () => new PARDISOSolver();

            c.AnalyticVelocityX = X => -0.33333333333333333332 + 1.2523108062960316903 * (X[1] + 0.11013981986815616002).Pow(3.0 / 5.0);
            c.AnalyticVelocityY = X => 0.0;
            c.AnalyticPressure = X => -1.4882576522041879973 * (1.9716158024185796207 * X[1] + 0.21715340932759252572).Pow(2.0 / 5.0) + 1.5508248735452649660;
            c.AnalyticTemperature = X => (1.9716158024185796207 * X[1] + 0.21715340932759252572).Pow(3.0 / 5.0);

            c.AddBoundaryValue("velocity_inlet_top", VariableNames.VelocityX, X => 1.0);
            c.AddBoundaryValue("velocity_inlet_top", VariableNames.VelocityY, X => 0.0);
            c.AddBoundaryValue("velocity_inlet_top", VariableNames.Temperature, X => 1.6);
            c.AddBoundaryValue("velocity_inlet_left", VariableNames.VelocityX, c.AnalyticVelocityX);
            c.AddBoundaryValue("velocity_inlet_left", VariableNames.VelocityY, c.AnalyticVelocityY);
            c.AddBoundaryValue("velocity_inlet_left", VariableNames.Temperature, c.AnalyticTemperature);
            c.AddBoundaryValue("velocity_inlet_right", VariableNames.VelocityX, c.AnalyticVelocityX);
            c.AddBoundaryValue("velocity_inlet_right", VariableNames.VelocityY, c.AnalyticVelocityY);
            c.AddBoundaryValue("velocity_inlet_right", VariableNames.Temperature, c.AnalyticTemperature);
            c.AddBoundaryValue("wall_bottom", VariableNames.Temperature, X => 0.4);

            c.PhysicsMode = PhysicsMode.LowMach;
            c.ThermodynamicPressureMode = ThermodynamicPressureMode.Constant;
            c.Reynolds = 10.0;
            c.Prandtl = 0.71;
            c.Gamma = 1.4;
            c.EoS = new MaterialLawLowMach(600.0, MaterialParamsMode.PowerLaw);

            c.Froude = 0.92303846;
            c.GravityDirection = new double[] { 0.0, -1.0 };

            c.PressureReferencePoint = new double[] { 0.0, 0.5 };
            c.PressureMeanValue = 0.0;

            c.PredictorApproximation = PredictorApproximations.BlockDiagonal;
            c.PressureStabilizationScaling = 0.0;
            c.PredictorApproximationUpdateCycle = 1;
            c.MaxNoSIMPLEsteps = 1000;
            c.SavePeriodSIMPLE = 1000;
            c.RelaxationFactorPressure = 0.5;
            c.RelexationFactorVelocity = 0.8;
            c.RelexationFactorTemperature = 1.0;
            c.RelaxationModeTemperature = RelaxationTypes.Implicit;
            c.ViscousPenaltyScaling = 1.0;
            c.MaxFlowSolverIterations = 1;
            c.MaxTemperatureSolverIterations = 1;
            c.PrintLinerSolverResults = false;

            int queryQuadOrder = 15;
            c.Queries.Add(
                "SolL2err_u",
                QueryLibrary.L2Error(VariableNames.VelocityX, c.AnalyticVelocityX, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_v",
                QueryLibrary.L2Error(VariableNames.VelocityY, c.AnalyticVelocityY, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_p",
                QueryLibrary.L2Error(VariableNames.Pressure, c.AnalyticPressure, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_T",
                QueryLibrary.L2Error(VariableNames.Temperature, c.AnalyticTemperature, queryQuadOrder));

            return c;
        }
    }
}
