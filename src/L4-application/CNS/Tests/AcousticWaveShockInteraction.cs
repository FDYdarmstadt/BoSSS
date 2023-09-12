using CNS.EquationSystem;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace CNS.Tests
{
    public class AcousticWaveShockInteraction
    {
        // The idea was to compare a case where the inital normal shock lies on a cell boundary (deactivating Artificial Viscosity) and a case where it is inside a cell. 
        // The goal was to asses the dampening of pertubations due to AV.
        // Turns out, that the shock move a little bit because of the pertubation, so Oscillations/Dampening occurs for both cases
        // TODO: Compare results with analytical solution
        public static void NonAllignedVSAlligned()
        {
            string dbPath = @"C:\Users\sebastian\Documents\BossDB\StationaryShockWave_Perturbation";
            //BoSSSpad.BoSSSshell.OpenOrCreateDatabase(dbPath);
            var c1_Alligned = ControlExamples_Supersonic.StationaryShockWave_Perturbation(AV:false,shockPosition : 3.15, dbPath:dbPath); 
            c1_Alligned.SessionName = c1_Alligned.SessionName + "Alligned";

            var c2_NonAlligned = ControlExamples_Supersonic.StationaryShockWave_Perturbation(AV: true, shockPosition: 3.14, dbPath: dbPath);
            c2_NonAlligned.SessionName = c2_NonAlligned.SessionName + "NonAlligned";

            var p_Alligned = new CNSProgram();
            p_Alligned.Init(c1_Alligned);
            p_Alligned.RunSolverMode();

            var p_NonAlligned = new CNSProgram();
            p_NonAlligned.Init(c2_NonAlligned);
            p_NonAlligned.RunSolverMode();



        } 
    }
}
