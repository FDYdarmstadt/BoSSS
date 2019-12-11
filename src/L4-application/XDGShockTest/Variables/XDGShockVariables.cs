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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using XDGShock.ShockCapturing;
using static BoSSS.Solution.CompressibleFlowCommon.Variable;

namespace XDGShockTest {
    public static class XDGShockVariables {

        /// <summary>
        /// The obligatory level set
        /// </summary>
        public static readonly Variable LevelSet = new Variable("levelSet", VariableTypes.Other);

        /// <summary>
        /// The local sensor value of a shock sensor
        /// </summary>
        public static readonly DerivedVariable<SinglePhaseField> Sensor = new DerivedVariable<SinglePhaseField>(
            "sensor",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XDGShockTestMain program) {
                dgField.Clear();

                // Should be removed later
                //Console.WriteLine("sensor L2-norm: " + sensorField.L2Norm());

                dgField.Acc(1.0, program.Sensor.GetSensorSinglePhaseField());
            });

        /// <summary>
        /// XDG artificial viscosity
        /// </summary>
        public static readonly DerivedVariable<SinglePhaseField> ArtificialViscosity = new DerivedVariable<SinglePhaseField>(
            "artificialViscosity",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XDGShockTestMain program) {
                dgField.Clear();

                dgField.Acc(1.0, program.ArtificialViscosityField);
            });

        /// <summary>
        /// The optional pressure field:
        /// \f[ p = f(\rho, \vec{m}, \rho E) \f], where
        /// \f[ f \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable<XDGField> Pressure = new DerivedVariable<XDGField>(
            "p",
            VariableTypes.Pressure,
            delegate (XDGField dgfield, XDGShockTestMain program) {

                IList<DGField> fields_A = new List<DGField>();
                IList<DGField> fields_B = new List<DGField>();

                foreach (XDGField field in program.ConservativeFields) {
                    fields_A.Add(field.GetSpeciesShadowField("A"));
                    fields_B.Add(field.GetSpeciesShadowField("B"));
                }

                dgfield.Clear();

                CellMask cellMask_A = program.LevelSetTracker.Regions.GetSpeciesMask(program.LevelSetTracker.GetSpeciesId("A"));
                CellMask cellMask_B = program.LevelSetTracker.Regions.GetSpeciesMask(program.LevelSetTracker.GetSpeciesId("B"));

                dgfield.GetSpeciesShadowField("A").ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                    new CellQuadratureScheme(true, cellMask_A),
                    fields_A.ToArray());

                dgfield.GetSpeciesShadowField("B").ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                    new CellQuadratureScheme(true, cellMask_B),
                    fields_B.ToArray());
            }
            );

        /// <summary>
        /// The optional velocity field:
        /// \f[ \vec{u} = \frac{1}{\rho} \vec{m} \f] 
        /// </summary>
        public static readonly Vector<DerivedVariable<XDGField>> Velocity = new Vector<DerivedVariable<XDGField>>(
            d => new DerivedVariable<XDGField>(
                "u" + d,
                VariableTypes.Velocity,
                delegate (XDGField dgField, XDGShockTestMain program) {
                    dgField.Clear();

                    CellMask cellMask_A = program.LevelSetTracker.Regions.GetSpeciesMask(program.LevelSetTracker.GetSpeciesId("A"));
                    CellMask cellMask_B = program.LevelSetTracker.Regions.GetSpeciesMask(program.LevelSetTracker.GetSpeciesId("B"));

                    dgField.GetSpeciesShadowField("A").ProjectQuotient(
                        1.0,
                        program.Momentum[d].GetSpeciesShadowField("A"),
                        program.Density.GetSpeciesShadowField("A"),
                        cellMask_A,
                        accumulateResult: true);

                    dgField.GetSpeciesShadowField("B").ProjectQuotient(
                        1.0,
                        program.Momentum[d].GetSpeciesShadowField("B"),
                        program.Density.GetSpeciesShadowField("B"),
                        cellMask_B,
                        accumulateResult: true);
                }
                ));
    }
}
