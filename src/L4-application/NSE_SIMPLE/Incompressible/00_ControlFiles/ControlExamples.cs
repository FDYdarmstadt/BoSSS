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

using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Queries;
using ilPSP.LinSolvers.PARDISO;
using System;
using static System.Math;

namespace NSE_SIMPLE.Incompressible {

    public static class ControlExamples {

        public static SIMPLEControl PoiseuilleFlow() {
            SIMPLEControl c = new SIMPLEControl();

            c.DbPath = @"..\..\Base\06_ZipDatabases\NUnitTests.zip";
            c.savetodb = false;

            c.GridGuid = new Guid("0062a338-a8a4-4d52-9b16-bc79379dd4d5");
            c.GridPartType = GridPartType.METIS;

            c.ProjectName = "Steady_SIMPLE 2D Channel";
            c.ProjectDescription = "NUnitTest for Steady_SIMPLE";

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

            // Auxiliary fields
            //c.FieldOptions.Add(
            //    "VelocityX*",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityY*",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityX'",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityY'",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityXRes",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityYRes",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "Pressure'",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "DivB4",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "DivAfter",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });

            c.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => 0.0);
            c.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
            c.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);

            c.Algorithm = SolutionAlgorithms.Steady_SIMPLE;
            c.NoOfTimesteps = 1;
            c.L2NormPressureCorrection = 1.0e-6;
            c.L2NormVelocityResidual = 1.0e-6;

            c.PredictorSolverFactory = () => new PARDISOSolver();
            c.CorrectorSolverFactory = () => new PARDISOSolver();

            c.AddBoundaryCondition("velocity_inlet", VariableNames.VelocityX, X => 1.0 - X[1] * X[1]);
            c.AddBoundaryCondition("velocity_inlet", VariableNames.VelocityY, X => 0.0);
            c.AddBoundaryCondition("pressure_outlet");
            c.AddBoundaryCondition("wall_lower");
            c.AddBoundaryCondition("wall_upper");

            c.PhysicsMode = PhysicsMode.Incompressible;
            c.Reynolds = 100.0;

            c.PredictorApproximation = PredictorApproximations.Identity;
            c.PressureStabilizationScaling = 0.0;
            c.PredictorApproximationUpdateCycle = 500;
            c.MaxNoSIMPLEsteps = 500;
            c.SavePeriodSIMPLE = 500;
            c.RelaxationFactorPressure = 1.0;
            c.RelexationFactorVelocity = 0.2;
            c.ViscousPenaltyScaling = 1.0;
            c.PrintLinerSolverResults = false;

            c.AnalyticVelocityX = X => 1.0 - X[1] * X[1];
            c.AnalyticVelocityY = X => 0.0;
            c.AnalyticPressure = X => -2.0 / 100.0 * X[0] + 0.2;

            int queryQuadOrder = 10;
            c.Queries.Add(
                "SolL2err_u",
                QueryLibrary.L2Error(VariableNames.VelocityX, c.AnalyticVelocityX, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_v",
                QueryLibrary.L2Error(VariableNames.VelocityY, c.AnalyticVelocityY, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_p",
                QueryLibrary.L2Error(VariableNames.Pressure, c.AnalyticPressure, queryQuadOrder));

            return c;
        }

        public static SIMPLEControl UnsteadyTaylorVortex() {
            SIMPLEControl c = new SIMPLEControl();

            c.DbPath = @"..\..\Base\06_ZipDatabases\NUnitTests.zip";
            c.savetodb = false;

            c.GridGuid = new Guid("06e506b2-9cb7-48c5-ba5d-22fc57645aac");
            c.GridPartType = GridPartType.METIS;

            c.ProjectName = "Unsteady_SIMPLE Taylor vortex";
            c.ProjectDescription = "NUnitTest for Unsteady_SIMPLE";

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

            // Auxiliary fields
            //c.FieldOptions.Add(
            //    "VelocityX*",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityY*",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityX'",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityY'",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityXRes",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            //c.FieldOptions.Add(
            //    "VelocityYRes",
            //    new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "Pressure'",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "DivB4",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "DivAfter",
                new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });

            c.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => -Cos(PI * X[0]) * Sin(PI * X[1]));
            c.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => Sin(PI * X[0]) * Cos(PI * X[1]));
            c.InitialValues_Evaluators.Add(VariableNames.Pressure, X => -0.25 * (Cos(2.0 * PI * X[0]) + Cos(2.0 * PI * X[1])));

            c.Algorithm = SolutionAlgorithms.Unsteady_SIMPLE;
            c.TimeOrder = 4;

            c.NoOfTimesteps = 5;
            c.Endtime = 2.0;
            c.L2NormPressureCorrection = 1.0e-5;
            c.L2NormVelocityResidual = 1.0e-5;

            c.PredictorSolverFactory = () => new PARDISOSolver();
            c.CorrectorSolverFactory = () => new PARDISOSolver();

            c.PhysicsMode = PhysicsMode.Incompressible;
            c.Reynolds = 100.0;

            c.PredictorApproximation = PredictorApproximations.Identity;
            c.PressureStabilizationScaling = 0.0;
            c.PredictorApproximationUpdateCycle = 45;
            c.MaxNoSIMPLEsteps = 45;
            c.SavePeriodSIMPLE = 500;
            c.RelaxationFactorPressure = 1.0;
            c.RelexationFactorVelocity = 1.0;
            c.ViscousPenaltyScaling = 1.0;
            c.PrintLinerSolverResults = false;

            c.PressureReferencePoint = new double[] { 0.0, 0.0 };
            c.PressureMeanValue = 0.0;

            c.AnalyticVelocityX = X => -Cos(PI * X[0]) * Sin(PI * X[1]) * Exp(-2.0 * PI * PI * 2.0 / 100.0);
            c.AnalyticVelocityY = X => Sin(PI * X[0]) * Cos(PI * X[1]) * Exp(-2.0 * PI * PI * 2.0 / 100.0);
            c.AnalyticPressure = X => -0.25 * (Cos(2.0 * PI * X[0]) + Cos(2.0 * PI * X[1])) * Exp(-4.0 * PI * PI * 2.0 / 100.0);

            int queryQuadOrder = 10;
            c.Queries.Add(
                "SolL2err_u",
                QueryLibrary.L2Error(VariableNames.VelocityX, c.AnalyticVelocityX, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_v",
                QueryLibrary.L2Error(VariableNames.VelocityY, c.AnalyticVelocityY, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_p",
                QueryLibrary.L2Error(VariableNames.Pressure, c.AnalyticPressure, queryQuadOrder));

            return c;
        }
    }
}
