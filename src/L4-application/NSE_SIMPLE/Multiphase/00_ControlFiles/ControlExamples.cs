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

namespace NSE_SIMPLE.Multiphase {

    public static class ControlExamples {

        public static SIMPLEControl UnsteadyMultiphaseWave() {
            MultiphaseSIMPLEControl c = new MultiphaseSIMPLEControl();

            c.DbPath = @"..\..\Base\06_ZipDatabases\NUnitTests.zip";
            c.savetodb = false;

            c.GridGuid = new Guid("a5135fbb-4243-4677-97f6-562860e4f95c");
            c.GridPartType = GridPartType.METIS;

            c.ProjectName = "MultiphaseWave";
            c.ProjectDescription = "Convected density wave";

            // Required fields
            c.FieldOptions.Add(
                VariableNames.VelocityX,
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                VariableNames.VelocityY,
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                VariableNames.Pressure,
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                "LevelSet",
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                "Density",
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            c.FieldOptions.Add(
                "Eta",
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            // Auxiliary fields
            c.FieldOptions.Add(
                "DivB4",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "DivAfter",
                new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorTest_x",
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorTest_y",
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorAna_x",
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });
            c.FieldOptions.Add(
                "OperatorAna_y",
                new FieldOpts() { Degree = 3, SaveToDB = FieldOpts.SaveToDBOpt.FALSE });

            c.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => 1.0);
            c.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
            c.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
            c.InitialValues_Evaluators.Add("LevelSet", X => 0.5 - 0.5 * Tanh(20.0 * X[0] - 12.0));

            c.Algorithm = SolutionAlgorithms.Unsteady_SIMPLE;
            c.NoOfTimesteps = 2;
            c.Endtime = 0.02;
            c.TimeOrder = 1;

            c.L2NormPressureCorrection = 1.0e-8;
            c.L2NormVelocityResidual = 1.0e-8;
            c.L2NormLevelSetResidual = 1.0E-8;

            c.PredictorSolverFactory = () => new PARDISOSolver();
            c.CorrectorSolverFactory = () => new PARDISOSolver();
            c.LevelSetSolverFactory = () => new PARDISOSolver();

            c.AddBoundaryCondition("velocity_inlet", VariableNames.VelocityX, X => 1.0);
            c.AddBoundaryCondition("velocity_inlet", VariableNames.VelocityY, X => 0.0);
            c.AddBoundaryCondition("velocity_inlet", "LevelSet", X => 1.0);
            c.AddBoundaryCondition("pressure_outlet");

            c.PhysicsMode = PhysicsMode.Multiphase;
            c.Reynolds = 1.0;
            c.EoS = new MaterialLawMultiphase(rho1: 1000.0, rho2: 1.0, mu1: 1.0, mu2: 1.0);

            c.PredictorApproximation = PredictorApproximations.Diagonal;
            c.PressureStabilizationScaling = 0.0;
            c.PredictorApproximationUpdateCycle = 54;
            c.MaxNoSIMPLEsteps = 54;
            c.SavePeriodSIMPLE = 100;
            c.RelaxationFactorPressure = 0.3;
            c.RelexationFactorVelocity = 0.7;
            c.RelaxationFactorLevelSet = 1.0;
            c.LevelSetRelaxationType = RelaxationTypes.Implicit;
            c.ViscousPenaltyScaling = 1.0;
            c.PrintLinerSolverResults = false;

            c.AnalyticVelocityX = X => 1.0;
            c.AnalyticVelocityY = X => 0.0;
            c.AnalyticPressure = X => 0.0;
            c.AnalyticLevelSet = X => 0.5 - 0.5 * Tanh(20.0 * X[0] - 12.4);
            c.AnalyticDensity = X => 1000.0 * (0.5 - 0.5 * Tanh(20.0 * X[0] - 12.4)) + (1.0 - (0.5 - 0.5 * Tanh(20.0 * X[0] - 12.4)));

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
            c.Queries.Add(
                "SolL2err_phi",
                QueryLibrary.L2Error("LevelSet", c.AnalyticLevelSet, queryQuadOrder));
            c.Queries.Add(
                "SolL2err_Rho",
                QueryLibrary.L2Error("Density", c.AnalyticDensity, queryQuadOrder));

            return c;
        }
    }
}
