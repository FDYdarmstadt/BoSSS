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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Utils;
using CNS.MaterialProperty;
using ilPSP;
using System;
using System.Diagnostics;
using System.Linq;

namespace CNS.IBM {

    /// <summary>
    /// Queries able to cope with the peculiarities of IBM applications
    /// </summary>
    public static class IBMQueries {

        /// <summary>
        /// L2 error with respect to given reference solution. The quadrature
        /// is determined from the settings in <see cref="IBMControl"/>
        /// </summary>
        /// <param name="fieldName"></param>
        /// <param name="referenceSolution"></param>
        /// <returns></returns>
        public static Query L2Error(string fieldName, Func<double[], double, double> referenceSolution) {
            return delegate (IApplication<AppControl> app, double time) {
                IProgram<CNSControl> program = app as IProgram<CNSControl>;
                if (program == null) {
                    throw new Exception();
                }

                ImmersedSpeciesMap speciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                IBMControl control = program.Control as IBMControl;
                if (speciesMap == null || control == null) {
                    throw new Exception(
                        "Query is only valid for immersed boundary runs");
                }

                SpeciesId species = speciesMap.Tracker.GetSpeciesId(control.FluidSpeciesName);
                int order = control.LevelSetQuadratureOrder;
                CellQuadratureScheme scheme = speciesMap.QuadSchemeHelper.GetVolumeQuadScheme(
                    species, true, speciesMap.SubGrid.VolumeMask);

                DGField dgField = app.IOFields.Single(f => f.Identification == fieldName);

                return dgField.L2Error(referenceSolution.Vectorize(time), order, scheme);
            };
        }

        public static Query Integral(string fieldName) {
            return delegate (IApplication<AppControl> app, double time) {
                IProgram<CNSControl> program = app as IProgram<CNSControl>;
                if (program == null) {
                    throw new Exception();
                }

                ImmersedSpeciesMap speciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                IBMControl control = program.Control as IBMControl;
                if (speciesMap == null || control == null) {
                    throw new Exception(
                        "Query is only valid for immersed boundary runs");
                }

                SpeciesId species = speciesMap.Tracker.GetSpeciesId(control.FluidSpeciesName);
                int order = control.LevelSetQuadratureOrder;
                CellQuadratureScheme scheme = speciesMap.QuadSchemeHelper.GetVolumeQuadScheme(
                    species, true, speciesMap.SubGrid.VolumeMask);

                DGField dgField = app.IOFields.Single(f => f.Identification == fieldName);

                return DGField.IntegralOverEx(scheme, new Func((X,U,j) => (U[0])), 2 , dgField);
            };
        }


        static double refIntegral = 0.0;
        static bool firstCall = false;

        public static Query IntegralError(string fieldName) {
            var subQuery = Integral(fieldName);
            
            return delegate (IApplication<AppControl> app, double time) {
                IProgram<CNSControl> program = app as IProgram<CNSControl>;
                if (program == null) {
                    throw new Exception();
                }

                if (!firstCall) {
                    refIntegral = subQuery(app, time);
                    firstCall = true;
                } 
                double subResult = subQuery(app, time);

                
                return subResult - refIntegral;
            };
        }

        /// <summary>
        /// L2 error of some quantity derived from the state vector (e.g.,
        /// entropy) with respect to given reference solution. The quadrature
        /// is determined from the settings in <see cref="IBMControl"/>
        /// </summary>
        /// <param name="quantityOfInterest"></param>
        /// <param name="referenceSolution"></param>
        /// <returns></returns>
        public static Query L2Error(Func<StateVector, double> quantityOfInterest, Func<double[], double, double> referenceSolution) {
            return delegate (IApplication<AppControl> app, double time) {
                IProgram<CNSControl> program = app as IProgram<CNSControl>;
                if (program == null) {
                    throw new Exception();
                }

                ImmersedSpeciesMap speciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                IBMControl control = program.Control as IBMControl;
                if (speciesMap == null || control == null) {
                    throw new Exception(
                        "Query is only valid for immersed boundary runs");
                }

                SpeciesId species = speciesMap.Tracker.GetSpeciesId(control.FluidSpeciesName);
                int order = control.LevelSetQuadratureOrder;
                CellQuadratureScheme scheme = speciesMap.QuadSchemeHelper.GetVolumeQuadScheme(
                    species, true, speciesMap.SubGrid.VolumeMask);
                var composititeRule = scheme.Compile(program.GridData, order);
                IChunkRulePair<QuadRule>[] chunkRulePairs = composititeRule.ToArray();

                DGField density = program.WorkingSet.Density;
                VectorField<DGField> momentum = program.WorkingSet.Momentum;
                DGField energy = program.WorkingSet.Energy;

                // Construct dummy field since L2Error is currently only supported
                // for Field's; However, _avoid_ a projection.
                DGField dummy = new SinglePhaseField(new Basis(program.GridData, 0));
                Material material = speciesMap.GetMaterial(double.NaN);
                int index = 0;
                double value = dummy.LxError(
                    (ScalarFunctionEx)delegate (int j0, int Len, NodeSet nodes, MultidimensionalArray result) {
                        MultidimensionalArray input = program.GridData.GlobalNodes.GetValue_Cell(nodes, j0, Len);

                        Chunk chunk = chunkRulePairs[index].Chunk;
                        QuadRule rule = chunkRulePairs[index].Rule;

                        if (chunk.i0 != j0 || chunk.Len != Len) {
                            throw new Exception();
                        }

                        if (rule.NoOfNodes != nodes.GetLength(0)) {
                            throw new Exception();
                        }

                        MultidimensionalArray rho = MultidimensionalArray.Create(chunk.Len, rule.NoOfNodes);
                        density.Evaluate(chunk.i0, chunk.Len, nodes, rho);

                        MultidimensionalArray[] m = new MultidimensionalArray[CNSEnvironment.NumberOfDimensions];
                        for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                            m[d] = MultidimensionalArray.Create(chunk.Len, rule.NoOfNodes);
                            momentum[d].Evaluate(chunk.i0, chunk.Len, nodes, m[d]);
                        }

                        MultidimensionalArray rhoE = MultidimensionalArray.Create(chunk.Len, rule.NoOfNodes);
                        energy.Evaluate(chunk.i0, chunk.Len, nodes, rhoE);

                        double[] X = new double[CNSEnvironment.NumberOfDimensions];
                        Vector3D mVec = new Vector3D();
                        for (int i = 0; i < chunk.Len; i++) {
                            for (int j = 0; j < rule.NoOfNodes; j++) {
                                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                                    X[d] = input[i, j, d];
                                    mVec[d] = m[d][i, j];
                                }

                                StateVector state = new StateVector(material, rho[i, j], mVec, rhoE[i, j]);
                                double qoi = quantityOfInterest(state);

                                Debug.Assert(
                                    !double.IsNaN(qoi),
                                    "Encountered node with unphysical state"
                                        + " (not able to determine quantity of interest)");

                                result[i, j] = qoi - referenceSolution(X, time);
                            }
                        }

                        index++;
                    },
                    (X, a, b) => (a - b) * (a - b),
                    composititeRule);

                // No value is NaN, but the results. How can this be?
                // => All values around 0, but values in void region are a little
                // farther away from the exact solution
                // => weights in the void zone sum up to something slightly negative
                Debug.Assert(
                    value >= 0,
                    "Encountered unphysical norm even though individual values where valid."
                    + " This indicates a problem with cut-cell quadrature.");

                return Math.Sqrt(value);
            };
        }


        /// <summary>
        /// Calculates a force along the zero iso surface of the Level-Set
        /// </summary>
        /// <param name="direction">Direction of the force projection, e.g. 0=x-axis, 1=y-axis</param>
        /// <returns></returns>
        static Query LiftOrDragForce(int direction) {
            return delegate (IApplication<AppControl> app, double time) {
                IProgram<CNSControl> program = app as IProgram<CNSControl>;
                if (program == null) {
                    throw new Exception();
                }

                ImmersedSpeciesMap speciesMap = program.SpeciesMap as ImmersedSpeciesMap;
                IBMControl control = program.Control as IBMControl;
                if (speciesMap == null || control == null) {
                    throw new Exception(
                        "Query is only valid for immersed boundary runs");
                }
                DGField density = program.WorkingSet.Density;
                VectorField<DGField> momentum = program.WorkingSet.Momentum;
                DGField energy = program.WorkingSet.Energy;

                return XDGUtils.GetIntegralOverZeroLevelSet(
                    speciesMap.Tracker,
                    GetSurfaceForce(
                        density,
                        momentum,
                        energy,
                        speciesMap,
                        direction,
                        speciesMap.Tracker._Regions.GetCutCellMask()),
                    control.MomentFittingVariant,
                    control.LevelSetQuadratureOrder,
                    speciesMap.QuadSchemeHelper);
            };
        }

        /// <summary>
        /// Gives the lift force of an immersed object, e.g airfoil or cylinder
        /// </summary>
        /// <returns></returns>
        public static Query IBMLiftForce() {
            return LiftOrDragForce(1);
        }

        /// <summary>
        /// Gives the lift force of an immersed object, e.g airfoil or cylinder
        /// </summary>
        /// <returns></returns>
        public static Query IBMDragForce() {
            return LiftOrDragForce(0);
        }

        /// <summary>
        /// Defines the force which is integrated over an immersed boundary, 
        /// called by <see cref="IBMQueries.LiftOrDragForce"/> 
        /// </summary>
        /// <param name="density"></param>
        /// <param name="momentum"></param>
        /// <param name="energy"></param>
        /// <param name="speciesMap"></param>
        /// <param name="direction">Direction of the force projection, e.g. 0=x-axis, 1=y-axis</param>
        /// <param name="cutCellMask">Cells intersected by the interface</param>
        /// <returns></returns>
        static ScalarFunctionEx GetSurfaceForce(
            DGField density, VectorField<DGField> momentum, DGField energy, ImmersedSpeciesMap speciesMap, int direction, CellMask cutCellMask) {

            return delegate (int j0, int Len, NodeSet nodes, MultidimensionalArray result) {
                int noOfNodes = nodes.GetLength(0);
                int D = nodes.SpatialDimension;

                double Reynolds = speciesMap.Control.ReynoldsNumber;
                double Mach = speciesMap.Control.MachNumber;
                double gamma = speciesMap.Control.EquationOfState.HeatCapacityRatio;

                double MachScaling = gamma * Mach * Mach;


                MultidimensionalArray rho = MultidimensionalArray.Create(Len, noOfNodes);
                density.Evaluate(j0, Len, nodes, rho);

                MultidimensionalArray[] m = new MultidimensionalArray[CNSEnvironment.NumberOfDimensions];
                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                    m[d] = MultidimensionalArray.Create(Len, noOfNodes);
                    momentum[d].Evaluate(j0, Len, nodes, m[d]);
                }

                MultidimensionalArray rhoE = MultidimensionalArray.Create(Len, noOfNodes);
                energy.Evaluate(j0, Len, nodes, rhoE);

                MultidimensionalArray gradRho = MultidimensionalArray.Create(Len, noOfNodes, D);
                density.EvaluateGradient(j0, Len, nodes, gradRho);

                MultidimensionalArray gradM = MultidimensionalArray.Create(Len, noOfNodes, D, D);
                for (int d = 0; d < D; d++) {
                    momentum[d].EvaluateGradient(
                        j0,
                        Len,
                        nodes,
                        gradM.ExtractSubArrayShallow(-1, -1, d, -1),
                        0,
                        0.0);
                }

                MultidimensionalArray normals = speciesMap.Tracker.GetLevelSetNormals(0, nodes, j0, Len);

                Vector3D mVec = new Vector3D();
                for (int i = 0; i < Len; i++) {
                    for (int j = 0; j < noOfNodes; j++) {
                        for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                            mVec[d] = m[d][i, j];
                        }

                        Material material = speciesMap.GetMaterial(double.NaN);
                        StateVector state = new StateVector(material, rho[i, j], mVec, rhoE[i, j]);

                        double mu = 0.0;
                        if (Reynolds != 0.0) {
                            mu = state.GetViscosity(j0 + j) / Reynolds;
                        }

                        double[,] gradU = new double[D, D];
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                // Apply chain rule
                                gradU[d1, d2] = (gradM[i, j, d1, d2] - state.Momentum[d1] / state.Density * gradRho[i, j, d2]) / state.Density;
                            }
                        }
                        double divU = gradU[0, 0] + gradU[1, 1];

                        switch (direction) {
                            // Attention: Changed sign, because normal vector is pointing inwards, not outwards!
                            case 0: // x-Direction
                                result[i, j, 0] = -state.Pressure/MachScaling * normals[i, j, 0]
                                    + mu * (2.0 * gradU[0, 0] - 2.0 / 3.0 * divU) * normals[i, j, 0]  //tau_11 * n_1
                                    + mu * (gradU[0, 1] + gradU[1, 0]) * normals[i, j, 1];            //tau_12 * n_2
                                break;

                            case 1: // y-Direction
                                result[i, j, 0] = -state.Pressure/MachScaling * normals[i, j, 1]
                                    + mu * (gradU[0, 1] + gradU[1, 0]) * normals[i, j, 0]             //tau_12 * n_1
                                    + mu * (2.0 * gradU[1, 1] - 2.0 / 3.0 * divU) * normals[i, j, 1]; //tau_22*n_2
                                break;

                            default:
                                throw new ArgumentException("Lift and Drag currently only in 2D implemented");
                        }
                    }
                }
            };
        }
    }
}
