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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using CNS.Diffusion;
using CNS.EquationSystem;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Globalization;
using System.IO;
using System.Linq;

namespace CNS.Tests.DiffusiveFlux {

    /// <summary>
    /// NUnit test, which testes the implementation of the SIPG and optimizedSIPG implementation. It is done
    /// by a comparison of one operator evaluation with a reference implementation stored in txt file, e.g
    /// "SIPGConsistencyRefSolution_pX" for different polynomial degree pX 
    /// </summary>
    public class SIPGConsistency : TestProgram<CNSControl> {


        int noOfCells = 3;
        CoordinateVector output;

        /// <summary>
        /// Performs one operator evaluation
        /// </summary>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            CoordinateMapping output_map = new CoordinateMapping(WorkingSet.ConservativeVariables);
            operatorFactory.GetDiffusiveOperator().ToSpatialOperator(WorkingSet).Evaluate(
                1.0, 0.0, new CoordinateMapping(WorkingSet.ConservativeVariables), null, output_map);
            output = new CoordinateVector(output_map);

            //int p = new int[] { WorkingSet.Density.Basis.Degree, WorkingSet.Momentum[0].Basis.Degree, WorkingSet.Energy.Basis.Degree }.Max();
            //String refSolutionFile = @"..\..\Tests\DiffusiveFlux\SIPGConsistencyRefSolution_p" + p + ".txt";
            //File.WriteAllLines(refSolutionFile, output.Select(d => d.ToString()).ToArray());

            //double[] refSol = File.ReadAllLines(@"..\..\Tests\ViscousFlux\SIPGConsistencyRefSolution_p1.txt").Select(c => Convert.ToDouble(c)).ToArray();

            //int test = 0;

            //for (int i = 0; i < output.Count; i++) {
            //    if (refSol[i] - output[i] > 1e-10) {
            //        test = test + 1;
            //    }
            //}       

            return 1.0;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="physTime"></param>
        /// <param name="timestepNo"></param>
        /// <param name="superSampling"></param>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
        }

        /// <summary>
        /// Performs the consistency test for polynomial degree p=1..3 and 
        /// the SIPG flux
        /// </summary>
        /// <param name="p">polynomial degree</param>
        [Test]
        public static void SIPGconsistencyTest([Range(1, 3)] int p) {
            int numOfCells = 3;

            string controlFile = String.Format(
                @"-c cs:CNS.Tests.DiffusiveFlux.SIPGConsistency.Control({0},{1},{2})",
                numOfCells,
                p,
                1);
            string refSolutionFile = @"..\..\Tests\DiffusiveFlux\SIPGConsistencyRefSolution_p" + p + ".txt";

            SIPGConsistency prog = null;
            Application<CNSControl>._Main(
                new string[] { controlFile },
                false,
                delegate () {
                    prog = new SIPGConsistency();
                    prog.noOfCells = numOfCells;
                    return prog;
                });

            double threshold = 1E-11;
            string[] refSolString = File.ReadAllLines(refSolutionFile);
            double[] refSol = refSolString.Select(c => Convert.ToDouble(c, NumberFormatInfo.InvariantInfo)).ToArray();

            Console.WriteLine("SIPG Consistency Test: [" + numOfCells + "x" + numOfCells + "] and p_dg=" + p);
            for (int i = 0; i < prog.output.Count; i++) {
                Assert.IsTrue(Math.Abs(refSol[i] - prog.output[i]) < threshold, "wrong flux in DOF: {0}, difference is {1}", i, Math.Abs(refSol[i] - prog.output[i]));
            }
        }
        /// <summary>
        /// Performs the consistency test for polynomial degree p=1..3 and 
        /// the optimized version of the SIPG flux
        /// </summary>
        /// <param name="p">polynomial degree</param>
        [Test]
        public static void OptimizedSIPGconsistencyTest([Range(1, 3)] int p) {
            int numOfCells = 3;

            string controlFile = String.Format(
                @"-c cs:CNS.Tests.DiffusiveFlux.SIPGConsistency.Control({0},{1},{2})",
                numOfCells,
                p,
                1);
            string refSolutionFile = @"..\..\Tests\DiffusiveFlux\SIPGConsistencyRefSolution_p" + p + ".txt";

            SIPGConsistency prog = null;
            Application<CNSControl>._Main(
                new string[] { controlFile },
                false,
                delegate () {
                    prog = new SIPGConsistency();
                    prog.noOfCells = numOfCells;
                    return prog;
                });

            double threshold = 1E-11;
            string[] refSolString = File.ReadAllLines(refSolutionFile);
            double[] refSol = refSolString.Select(c => Convert.ToDouble(c, NumberFormatInfo.InvariantInfo)).ToArray();

            Console.WriteLine("SIPG Consistency Test: [" + numOfCells + "x" + numOfCells + "] and p_dg=" + p);
            for (int i = 0; i < prog.output.Count; i++) {
                Assert.IsTrue(Math.Abs(refSol[i] - prog.output[i]) < threshold, "wrong flux in DOF: {0}, difference is {1}", i, Math.Abs(refSol[i] - prog.output[i]));
            }
        }

        /// <summary>
        /// Creates the control file
        /// </summary>
        /// <param name="noOfCells">number of Cells in each direction=3 for the reference solution</param>
        /// <param name="dgDegree">polynomial DG degree</param>
        /// <param name="optimized">SIPG flux: 0=normal, 1=optimized</param>
        /// <returns></returns>
        new public static CNSControl Control(int noOfCells, int dgDegree, int optimized) {
            CNSControl c = new CNSControl();

            c.savetodb = false;
            c.NoOfTimesteps = 1;
            c.CFLFraction = 0.1;

            c.ActiveOperators = Operators.Diffusion;
            switch (optimized) {
                case 0:
                    c.DiffusiveFluxType = DiffusiveFluxTypes.SIPG;
                    break;
                case 1:
                    c.DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;
                    break;
                default:
                    throw new ArgumentException("Wrong diffusive type, only SIPG (0) and OptimizedSIPG (1) exist");
            }
            c.SIPGPenaltyScaling = 1.0;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.EquationOfState = new IdealGas(2.0);
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 1.0;
            c.ViscosityLaw = new SutherlandLaw();
            c.MachNumber = 1 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum[0], dgDegree);
            c.AddVariable(CompressibleVariables.Momentum[1], dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.GridFunc = delegate () {
                GridCommons grid = Grid2D.Cartesian2DGrid(
                    GenericBlas.Linspace(0.0, 2 * Math.PI, noOfCells + 1),
                    GenericBlas.Linspace(0.0, 2 * Math.PI, noOfCells + 1),
                    periodicX: false,
                    periodicY: false);
                grid.EdgeTagNames.Add(1, "supersonicinlet");
                grid.DefineEdgeTags(x => 1);

                return grid;
            };

            Func<double[], double> rho = X => 0.1 * Math.Sin(X[0]) + 0.1 * Math.Cos(X[1]) + 1.0;
            Func<double[], double> u0 = X => Math.Sin(X[0]) + Math.Cos(X[1]) + 1.0;
            Func<double[], double> u1 = X => Math.Cos(X[0]) + Math.Sin(X[1]) + 1.0;
            Func<double[], double> p = X => Math.Cos(X[0]) + Math.Sin(X[0]) + 4.0;

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, rho);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[0], u0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity[1], u1);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, p);

            c.AddBoundaryValue("supersonicinlet", CompressibleVariables.Density, rho);
            c.AddBoundaryValue("supersonicinlet", CNSVariables.Velocity[0], u0);
            c.AddBoundaryValue("supersonicinlet", CNSVariables.Velocity[1], u1);
            c.AddBoundaryValue("supersonicinlet", CNSVariables.Pressure, p);

            return c;
        }
    }
}
