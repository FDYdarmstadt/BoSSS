using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Tests {

    [TestFixture]
    public static class ThreePhaseTests {


        [Test]
        public static void ThreePhaseContactLine_TensionBalance([Values(true)] bool UseGravity) {


            var C = ZwoLevelSetSolver.ControlFiles.Droplet.Wiki(p: 4, AMRlvl: 2, SlipLength: 0, useGravity: UseGravity);
            C.NoOfTimesteps = 1;
            
            //if(UseGravity == false) {
            //    // turn off gravity
            //    var toRemove = C.InitialValues.Keys.Where(name => name.StartsWith("Gravity")).ToArray();
            //    foreach(var name in toRemove) 
            //        C.InitialValues.Remove(name);
            //}
           

            using(var S = new ZLS()) {
                S.Init(C);
                S.RunSolverMode();

                double[] L2_vel = new double[2];
                for(int d = 0; d < 2; d++) {
                    L2_vel[d] = S.Velocity[d].L2Norm();
                    Console.WriteLine($"Velocity {d} norm: " + L2_vel[d]);
                }
                // with gravity
                // Velocity 0 norm: 9.070818979018125E-06
                // Velocity 1 norm: 1.5000248940334127E-05

                if(UseGravity) {
                    Assert.Less(L2_vel[0], 5.0e-5, "VelocityX out of bound");
                    Assert.Less(L2_vel[1], 5.0e-5, "VelocityY out of bound");
                }
            }


        }


    }
}
