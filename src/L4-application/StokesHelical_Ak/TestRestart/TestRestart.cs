using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.Utils;
using ilPSP;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using NUnit.Framework.Constraints;
using NUnitLite;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using BoSSS.Application.BoSSSpad;
using BoSSS.Solution.Tecplot;
using ilPSP.Tracing;
using ilPSP.LinSolvers;
using BoSSS.Solution.NSECommon;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Platform;
using MPI.Wrappers;
using System.IO;
using System.Collections;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.Statistic;
using static BoSSS.Application.BoSSSpad.BoSSSshell;
using BoSSS.Solution.Control;
using BoSSS.Solution.GridImport;

namespace StokesHelical_Ak.TestRestart {
    [TestFixture]
    internal class TestRestart {
        // TestRestart for BDF1
        Guid sessionID;
        [Test]
        static public void Restart_Comparison_Regular_Grid_BDF3_with_R0fix() {
            // Initialize the number of timesteps and cell configurations.
            int timeSteps = 6;
            int noOfCells = 4; // Note: noOfCellsXi is aquidistant.
            int deGREE = 2;

            // Create and configure the HelicalControl object for the first simulation run.
            HelicalControl referenceContrl = StokesHelical_Ak.DNS_Hagen_Poiseulle.HagenPoiseulle(noOfCellsR: noOfCells, noOfCellsXi: noOfCells, dtRefining: timeSteps, bdfOrder: 3, degree: deGREE);

            // Set the database path and related settings.

            string restartDB = $"restart_HG_BDF{referenceContrl.TimeSteppingScheme}_degree={ deGREE}_noOfCellsR={ noOfCells}_noOfCellsXi={ noOfCells}_{DateTime.Now.ToString("MMMdd_HHmm")}";
            string restartDBfullPath = Path.Combine(Directory.GetCurrentDirectory(), restartDB);
            DatabaseUtils.CreateDatabase(restartDBfullPath);
            referenceContrl.DbPath = restartDBfullPath;
            referenceContrl.savetodb = referenceContrl.DbPath != null;
            referenceContrl.ProjectName = "DNS_first_try";
            referenceContrl.SessionName = $"degree={deGREE} noOfCellsR={noOfCells} noOfCellsXi={noOfCells}";

            // Initialize and run the solver for the first simulation.
            var solver = new HelicalMain();
            solver.Init(referenceContrl);
            SplittingTimestepper.BackupTimestep = 3;
            solver.RunSolverMode();
            SplittingTimestepper.BackupTimestep = -1;

            // Assert conditions to ensure correct configuration.
            Assert.That(referenceContrl.R0fixOn == true, "R0_fix should be true");
            Assert.That(referenceContrl.PressureReferencePoint == true, "Calculation should proceed without PRP since R0 is on");

            // Issue reminders based on global multipliers and minimum radius conditions.
            if(Globals.activeMult == Globals.Multiplier.one && referenceContrl.rMin < 10e-6) {
                Console.WriteLine("Friendly Reminder: Multiplier One and rMin<10e-6");
            }
            Console.WriteLine($"Remember: R0_fix is {referenceContrl.R0fixOn}");
            Console.WriteLine($"Remember: Calculating with PRP is {referenceContrl.PressureReferencePoint}");

            // Clone the original configuration for a second run with modifications for restart.
            Guid sessionID = solver.CurrentSessionInfo.ID;
            HelicalControl restartControl = referenceContrl.CloneAs();
            restartControl.InitialValues.Clear();
            restartControl.InitialValues_Evaluators.Clear();

            int restart_time_step;

            if(restartControl.TimeSteppingScheme == TimeSteppingScheme.BDF3) {
                restart_time_step = (timeSteps / 2) - 2; // -2 weil ZUSÄTZLICH noch ZWEI weitere weitere Zeitschritte von Nöten sind!
                restartControl.NoOfTimesteps = timeSteps - (restart_time_step + 2); // Wieder Rückgängig wegen BFD3
            } else if(restartControl.TimeSteppingScheme == TimeSteppingScheme.ImplicitEuler) {
                restart_time_step = timeSteps / 2;
                restartControl.NoOfTimesteps = timeSteps - restart_time_step;
            } else {
                throw new ArgumentException("Unsupported BDF scheme: :(( ");
            }


            restartControl.RestartInfo = new Tuple<Guid, TimestepNumber>(sessionID, new TimestepNumber(restart_time_step));
            restartControl.GridFunc = null;

            // Initialize and run the solver for the second simulation.
            var solver2 = new HelicalMain();
            solver2.Init(restartControl);
            solver2.RunSolverMode();
            Guid sessionID2 = solver2.CurrentSessionInfo.ID;

            {
                // Initialize the database connection
                var testDb = DatabaseInfo.CreateOrOpen(referenceContrl.DbPath);
                // Retrieve the session that was restarted using its session ID
                var restartedSession = testDb.Sessions.FirstOrDefault(s => s.ID == sessionID2);
                if(restartedSession == null) {
                    throw new InvalidOperationException("Restarted session not found.");
                }

                // Extract and print the major timestep numbers for the restarted session
                var timestepNumbers = restartedSession.Timesteps.Select(tsi => tsi.TimeStepNumber.MajorNumber).ToArray();
                Console.WriteLine($"tsiNumbers = {string.Join(" ", timestepNumbers)}");

                // Flag to track if any comparisons fail
                bool comparisonFailed = false;

                // Retrieve the reference session from which this session was restarted
                var referenceSession = testDb.Sessions.Single(s => s.ID.Equals(restartedSession.RestartedFrom));
                // Ensure the restarted session and reference session are not the same
                if(restartedSession.ID == referenceSession.ID) {
                    throw new InvalidOperationException("Session IDs are the same.");
                }
                // Iterate through each timestep in the restarted session for comparison
                foreach(var tsiRestart in restartedSession.Timesteps) {
                    // Find the matching timestep in the reference session
                    var tsiReference = referenceSession.Timesteps.Single(t => t.TimeStepNumber.Equals(tsiRestart.TimeStepNumber));

                    // Compare each field within the timestep
                    foreach(var field in tsiRestart.Fields) {
                        var referenceField = tsiReference.Fields.Single(f => f.Identification == field.Identification);
                        // Adjust the reference field coordinates by subtracting the current field coordinates
                        referenceField.Coordinates.Acc(-1.0, field.Coordinates);
                        // OutPut
                        Console.WriteLine($"Loaded data at timestep {tsiRestart.TimeStepNumber.MajorNumber} for field {field.Identification} are not exact: L2-norm = {referenceField.L2Norm()}");
                        // Check if the adjusted coordinates' L2 norm is greater than 0 (indicating a discrepancy)
                        Assert.LessOrEqual(referenceField.L2Norm(), 1E-13, $"ERROR!!!:Loaded data at timestep {tsiRestart.TimeStepNumber.MajorNumber} for field {field.Identification} are not exact: L2-norm = {referenceField.L2Norm()}. Should be less or equal {1E-13}");

                        // Check if the adjusted coordinates' L2 norm is greater than 0 (indicating a discrepancy)
                        if(referenceField.L2Norm() > 1E-13) {
                            Console.WriteLine($"Loaded data at timestep {tsiRestart.TimeStepNumber.MajorNumber} for field {field.Identification} are not exact: L2-norm = {referenceField.L2Norm()}");
                            comparisonFailed = true;
                        }
                    }
                }

                // Assert that there were no comparison failures
                if(comparisonFailed) {
                    Console.WriteLine("Comparison between reference solution and restart solution not equal.");
                    throw new InvalidOperationException("Comparison between reference solution and restart solution not equal.");
                }
                Console.WriteLine("Test successfull!");
            }
        }
    }
}
