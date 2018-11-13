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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Utils;
using CNS.MaterialProperty;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace CNS {

    /// <summary>
    /// Static class for CNS specific queries
    /// </summary>
    public static class CNSQueries {

        private static EdgeIntegral dragIntegral;
        private static EdgeIntegral liftIntegral;
        private static TextWriter logger;

        private static DGField densityField;
        private static DGField m0Field;
        private static DGField m1Field;
        private static DGField energyField;

        private static SinglePhaseField DrhoDx;
        private static SinglePhaseField DrhoDy;

        private static SinglePhaseField Dm0Dx;
        private static SinglePhaseField Dm0Dy;
        private static SinglePhaseField Dm1Dx;
        private static SinglePhaseField Dm1Dy;
        
        /// <summary>
        /// Calculates the lift and drag force on given edges, e.g edges around an airfoil.
        /// Currently only in 2D!!!
        /// </summary>
        /// <param name="edgeTagName">EdgeTag on which the drag force will be calculated</param>
        /// <returns>total drag force</returns>
        /// <remarks>It assumes that pressure and velocity fields exists</remarks>
        public static Query LiftAndDragForce(String edgeTagName) {
            return delegate (IApplication<AppControl> app, double time) {

                if (dragIntegral == null) {
                    densityField = app.IOFields.SingleOrDefault((f) => f.Identification == "rho");
                    m0Field = app.IOFields.SingleOrDefault((f) => f.Identification == "m0");
                    m1Field = app.IOFields.SingleOrDefault((f) => f.Identification == "m1");
                    energyField = app.IOFields.SingleOrDefault((f) => f.Identification == "rhoE");

                    DrhoDx = new SinglePhaseField(densityField.Basis);
                    DrhoDy = new SinglePhaseField(densityField.Basis);

                    Dm0Dx = new SinglePhaseField(m0Field.Basis);
                    Dm0Dy = new SinglePhaseField(m0Field.Basis);
                    Dm1Dx = new SinglePhaseField(m1Field.Basis);
                    Dm1Dy = new SinglePhaseField(m1Field.Basis);

                    CNSControl c = app.Control as CNSControl;
                    Program prog = app as Program;

                    byte edgeTag = app.Grid.EdgeTagNames.First(item => item.Value.Equals(edgeTagName)).Key;

                    CoordinateMapping mapping = new CoordinateMapping(
                        densityField, m0Field, m1Field, energyField, DrhoDx, DrhoDy, Dm0Dx, Dm0Dy, Dm1Dx, Dm1Dy);

                    dragIntegral = new EdgeIntegral((BoSSS.Foundation.Grid.Classic.GridData) (app.GridData), edgeTag, new ForceFlux(
                        c.ReynoldsNumber, prog.SpeciesMap.GetMaterial(double.NaN), 0), mapping);
                    liftIntegral = new EdgeIntegral((BoSSS.Foundation.Grid.Classic.GridData) (app.GridData), edgeTag, new ForceFlux(
                        c.ReynoldsNumber, prog.SpeciesMap.GetMaterial(double.NaN), 1), mapping);

                    if (logger == null && app.CurrentSessionInfo.ID != Guid.Empty && app.MPIRank == 0) {
                        logger = app.DatabaseDriver.FsDriver.GetNewLog("LiftAndDragForce", app.CurrentSessionInfo.ID);
                        string header = "PhysTime\t LiftForce\t DragForce";
                        logger.WriteLine(header);
                    }
                }

                DrhoDx.Clear();
                DrhoDy.Clear();

                DrhoDx.Derivative(1.0, densityField, 0);
                DrhoDy.Derivative(1.0, densityField, 1);

                Dm0Dx.Clear();
                Dm0Dy.Clear();
                Dm1Dx.Clear();
                Dm1Dy.Clear();

                Dm0Dx.Derivative(1.0, m0Field, 0);
                Dm0Dy.Derivative(1.0, m0Field, 1);
                Dm1Dx.Derivative(1.0, m1Field, 0);
                Dm1Dy.Derivative(1.0, m1Field, 1);

                double lift = liftIntegral.Evaluate();
                double drag = dragIntegral.Evaluate();

                if (logger != null && app.MPIRank == 0) {
                    string line = time + "\t" + lift + "\t" + drag;
                    logger.WriteLine(line);
                    logger.Flush();
                }
                return drag;
            };
        }

        /// <summary>
        /// Force flux required for the construction of an
        /// <see cref="EdgeIntegral"/>.
        /// </summary>
        private class ForceFlux : EdgeIntegral.EdgeFlux {

            private double ReynoldsNumber;
            private Material material;
            private int direction;

            public ForceFlux(double reynolds, Material material, int direction) {
                this.material = material;
                ReynoldsNumber = reynolds;
                this.direction = direction;
            }

            /// <summary>
            /// Evaluates to
            /// \f$ 
            /// p \vec{n} \cdot \vec{e_i} - (\tau \cdot \vec{n}) \cdot \vec{e_i}, where i is the given direction
            /// \f$ 
            /// </summary>
            protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
                StateVector state = new StateVector(
                    material, Uin[0], new Vector(Uin[1], Uin[2], 0.0), Uin[3]);

                double mu = 0.0;
                if (ReynoldsNumber != 0.0) {
                    mu = state.GetViscosity(0) / ReynoldsNumber;
                }

                double[] gradRho = new double[2];
                gradRho[0] = Uin[4]; // "DrhoDx"
                gradRho[1] = Uin[5]; // "DrhoDy"

                double[,] gradM = new double[2, 2];
                gradM[0, 0] = Uin[6]; // "Dm0Dx"
                gradM[0, 1] = Uin[7]; // "Dm0Dy"
                gradM[1, 0] = Uin[8]; // "Dm1Dx"
                gradM[1, 1] = Uin[9]; // "Dm1Dy"

                double[,] gradU = new double[2, 2];
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        // Apply chain rule
                        gradU[i, j] = (gradM[i, j] - state.Momentum[i] / state.Density * gradRho[j]) / state.Density;
                    }
                }
                double divU = gradU[0, 0] + gradU[1, 1];

                switch (direction) {
                    case 0: // x-Direction
                        return state.Pressure * normal[0]
                            - mu * (2.0 * gradU[0, 0] - 2.0 / 3.0 * divU) * normal[0]  //tau_11 * n_1
                            - mu * (gradU[0, 1] + gradU[1, 0]) * normal[1];            //tau_12 * n_2
                    case 1: // y-Direction
                        return state.Pressure * normal[1]
                            - mu * (gradU[0, 1] + gradU[1, 0]) * normal[0]             //tau_12 * n_1
                            - mu * (2.0 * gradU[1, 1] - 2.0 / 3.0 * divU) * normal[1]; //tau_22*n_2
                    default:
                        throw new ArgumentException("Lift and Drag currently only in 2D implemented");
                }

            }

            /// <summary>
            /// Returns { "rho", "m0", "m1", "rhoE", "DrhoDx", "DrhoDy", "Dm0Dx", "Dm0Dy", "Dm1Dx", "Dm1Dy" }.
            /// </summary>
            public override IList<string> ArgumentOrdering {
                get {
                    //                    0      1     2     3       4         5         6        7        8        9
                    return new string[] { "rho", "m0", "m1", "rhoE", "DrhoDx", "DrhoDy", "Dm0Dx", "Dm0Dy", "Dm1Dx", "Dm1Dy" };
                }
            }
        }
    }
}
