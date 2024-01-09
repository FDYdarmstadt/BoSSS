using BoSSS.Solution;
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
        public static void Test00()
        {
            Application.InitMPI();
            //BoSSSpad.BoSSSshell.OpenOrCreateDatabase(dbPath);
            var c = ControlExamples_Supersonic.AcousticWave(dbPath: null, perStartTime: 0.0, endTime: 32); 
            var p = new CNSProgram();
            p.Init(c);
            p.RunSolverMode();




        } 
    }
}
