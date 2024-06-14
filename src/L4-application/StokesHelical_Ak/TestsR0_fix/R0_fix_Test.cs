using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace StokesHelical_Ak.NUnitTestsR0_fix {

    [TestFixture]
    static public class R0_fix_Test {


        /// <summary>
        /// checks if the solution which the solver provides complies with the compatibility conditions,
        /// <see cref="R0fix.CheckSolutionR0Compatibility(SinglePhaseField, SinglePhaseField, SinglePhaseField, SinglePhaseField)"/>.
        /// </summary>
        [Test]
        static public void TestR0_Fix_AtSolution([Values(2, 3, 4, 5, 6)] int pOrder) {
            using(HelicalMain p = new HelicalMain()) {
                var ctrl = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: 32, noOfCellsXi: 32,  rMin: 0.0);
                // overwrite the globals
                Globals.activeMult = Globals.Multiplier.Bsq;
                Globals.pressureStabilEtaMom = false;  // no pressure stabilization for this test
                p.Init(ctrl);
                p.RunSolverMode();
                Assert.That(Globals.activeMult == Globals.Multiplier.Bsq,
                String.Format("Multiplier is not f(r)=B(r)^2"));
                // Note: Fk, 13dec23: I think, `pressureStabilConti`  is essential, therefore it is turned on globally.
                //Assert.That(Globals.pressureStabilConti == false,
                //    String.Format(
                //        "Pressure Stabilization in Conti equation is true. No pressure stabilization for this unit test!")
                //    );
                Assert.That(Globals.pressureStabilEtaMom == false,
                    String.Format(
                        "Pressure Stabilization in Eta momentum equation is true. No pressure stabilization for this unit test!")
                    );
                R0fix myR0fix = new R0fix(p.CurrentSolution.Mapping, p.Control.rMin);
                myR0fix.CheckSolutionR0Compatibility(p.ur, p.uxi, p.ueta, p.Pressure);
                Assert.That(ctrl.R0fixOn == true, "R0_fix have to be true!");
                Assert.That(ctrl.PressureReferencePoint == true, "Pressure Reference Point have to be true!");
                Console.WriteLine("R0_fix is working correctly!");
            }
        }

        /*
        /// <summary>
        /// verifies certain properties which the <see cref="R0fix"/> should comply with:
        /// - Prolongation and subsequent restriction must be an identity operation
        /// </summary>
        [Test]
        static public void TestR0_Fix([Values(2, 3, 4, 5, 6)] int pOrder) {
            using(HelicalMain p = new HelicalMain()) {
                var ctrl = StokesHelical_Ak.DNS_Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 3, noOfCellsXi: 3, dtRefining: 1, _DbPath: null, bdfOrder: 1);
                ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

                // overwrite the globals
                Globals.activeMult = Globals.Multiplier.Bsq;
                Globals.pressureStabilEtaMom = false;  // no pressure stabilization for this test
                p.Init(ctrl);
                p.RunSolverMode();

                p.CurrentSolution.SaveToTextFile("U0.txt");
                
                R0fix myR0fix = new R0fix(p.CurrentSolution.Mapping, p.Control.rMin);
                myR0fix.InternalChecks();
            }
        }
        */
    }
}
