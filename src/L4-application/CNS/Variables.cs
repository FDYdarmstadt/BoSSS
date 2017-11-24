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
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Timestepping;
using CNS.IBM;
using CNS.ShockCapturing;
using ilPSP;
using System;
using System.IO;
using System.Linq;
using static CNS.Variable;

namespace CNS {

    /// <summary>
    /// Common variables that are often used in control files. Note that
    /// additional custom variables (e.g., for debugging purposes) can
    /// analogously be specified in the control file itself.
    /// </summary>
    public static class Variables {

        /// <summary>
        /// <see cref="VariableTypes.Density"/>
        /// </summary>
        public static readonly Variable Density = new Variable("rho", VariableTypes.Density);

        /// <summary>
        /// <see cref="VariableTypes.Momentum"/>
        /// </summary>
        public static readonly Vector<Variable> Momentum = new Vector<Variable>(
            d => new Variable("m" + d, VariableTypes.Momentum));

        /// <summary>
        /// <see cref="VariableTypes.Energy"/>
        /// </summary>
        public static readonly Variable Energy = new Variable("rhoE", VariableTypes.Energy);

        /// <summary>
        /// The optional velocity field:
        /// \f[ \vec{u} = \frac{1}{\rho} \vec{m} \f] 
        /// </summary>
        public static readonly Vector<DerivedVariable> Velocity = new Vector<DerivedVariable>(
            d => new DerivedVariable(
                "u" + d,
                VariableTypes.Velocity,
                delegate (DGField u, CellMask cellMask, IProgram<CNSControl> program) {
                    u.Clear();
                    u.ProjectQuotient(
                        1.0,
                        program.WorkingSet.Momentum[d],
                        program.WorkingSet.Density,
                        cellMask,
                        accumulateResult: true);
                }));

        /// <summary>
        /// The optional inner energy field:
        /// \f[ e = E - \frac{1}{2} u_i u_i \f] 
        /// </summary>
        public static readonly DerivedVariable SpecificInnerEnergy = new DerivedVariable(
            "e",
            VariableTypes.Other,
            delegate (DGField e, CellMask cellMask, IProgram<CNSControl> program) {
                e.Clear();
                e.ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.SpeciesMap.GetMaterial(double.NaN)).SpecificInnerEnergy,
                    new CellQuadratureScheme(true, cellMask),
                    program.WorkingSet.ConservativeVariables);
            });

        /// <summary>
        /// The optional temperature field:
        /// \f[ (\kappa - 1) \kappa Ma_\infty e\f] 
        /// </summary>
        public static readonly DerivedVariable Temperature = new DerivedVariable(
            "T",
            VariableTypes.Other,
            delegate (DGField T, CellMask cellMask, IProgram<CNSControl> program) {
                T.Clear();
                T.ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.SpeciesMap.GetMaterial(double.NaN)).Temperature,
                    new CellQuadratureScheme(true, cellMask),
                    program.WorkingSet.ConservativeVariables);
            });

        /// <summary>
        /// The optional pressure field:
        /// \f[ p = f(\rho, \vec{m}, \rho E) \f], where
        /// \f[ f \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable Pressure = new DerivedVariable(
            "p",
            VariableTypes.Pressure,
            delegate (DGField p, CellMask cellMask, IProgram<CNSControl> program) {
                p.Clear();
                p.ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.SpeciesMap.GetMaterial(double.NaN)).Pressure,
                    new CellQuadratureScheme(true, cellMask),
                    program.WorkingSet.ConservativeVariables);
            });

        /// <summary>
        /// The optional local Mach number field
        /// \f[ Ma = \frac{\Vert u \Vert}{a} \f], where the formula for the
        /// local speed of sound \f[ a \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable LocalMachNumber = new DerivedVariable(
            "Ma",
            VariableTypes.Other,
            delegate (DGField Ma, CellMask cellMask, IProgram<CNSControl> program) {
                Ma.Clear();
                Ma.ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.SpeciesMap.GetMaterial(double.NaN)).LocalMachNumber,
                    new CellQuadratureScheme(true, cellMask),
                    program.WorkingSet.ConservativeVariables);
            });

        /// <summary>
        /// The optional kinetic energy field:
        /// \f[ K = \rho u_i u_i \f] 
        /// </summary>
        public static readonly DerivedVariable KineticEnergy = new DerivedVariable(
            "K",
            VariableTypes.Other,
            delegate (DGField K, CellMask cellMask, IProgram<CNSControl> program) {
                K.Clear();
                K.ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.SpeciesMap.GetMaterial(double.NaN)).KineticEnergy,
                    new CellQuadratureScheme(true, cellMask),
                    program.WorkingSet.ConservativeVariables);
            });

        /// <summary>
        /// The optional entropy field (see
        /// <see cref="MaterialProperty.IEquationOfState.GetEntropy"/>).
        /// </summary>
        public static readonly DerivedVariable Entropy = new DerivedVariable(
            "S",
            VariableTypes.Other,
            delegate (DGField S, CellMask cellMask, IProgram<CNSControl> program) {
                S.Clear();
                S.ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.SpeciesMap.GetMaterial(double.NaN)).Entropy,
                    new CellQuadratureScheme(true, cellMask),
                    program.WorkingSet.ConservativeVariables);
            });

        /// <summary>
        /// The local CFL restriction
        /// </summary>
        public static readonly DerivedVariable CFL = new DerivedVariable(
            "cfl",
            VariableTypes.Other,
            delegate (DGField cfl, CellMask cellMask, IProgram<CNSControl> P) {
                if (P.FullOperator == null) {
                    return;
                }

                // Query each cell individually so we get local results
                for (int i = 0; i < P.Grid.NoOfUpdateCells; i++) {
                    double localCFL = P.FullOperator.CFLConstraints.Min(c => c.GetLocalStepSize(i, 1));
                    cfl.SetMeanValue(i, localCFL);
                }
            });

        /// <summary>
        /// The local convective CFL restriction
        /// </summary>
        public static readonly DerivedVariable CFLConvective = new DerivedVariable(
            "cflConvective",
            VariableTypes.Other,
            delegate (DGField cfl, CellMask cellMask, IProgram<CNSControl> P) {
                if (P.FullOperator == null) {
                    return;
                }

                // Query each cell individually so we get local results
                TimeStepConstraint cflConstraint = P.FullOperator.CFLConstraints.OfType<Convection.ConvectiveCFLConstraint>().Single();
                for (int i = 0; i < P.Grid.NoOfUpdateCells; i++) {
                    double localCFL = cflConstraint.GetLocalStepSize(i, 1);
                    cfl.SetMeanValue(i, localCFL);
                }
            });

        /// <summary>
        /// The local diffusive CFL restriction based on artificial viscosity
        /// </summary>
        public static readonly DerivedVariable CFLArtificialViscosity = new DerivedVariable(
            "cflArtificialViscosity",
            VariableTypes.Other,
            delegate (DGField cfl, CellMask cellMask, IProgram<CNSControl> P) {
                if (P.FullOperator == null) {
                    return;
                }

                // Query each cell individually so we get local results
                TimeStepConstraint cflConstraint = P.FullOperator.CFLConstraints.OfType<ArtificialViscosityCFLConstraint>().Single();
                for (int i = 0; i < P.Grid.NoOfUpdateCells; i++) {
                    double localCFL = cflConstraint.GetLocalStepSize(i, 1);
                    cfl.SetMeanValue(i, localCFL);
                }
            });

        /// <summary>
        /// The local residual
        /// </summary>
        public static readonly DerivedVariable Residual = new DerivedVariable(
            "residual",
            VariableTypes.Other,
            delegate (DGField residual, CellMask cellMask, IProgram<CNSControl> program) {
                residual.Clear();

                throw new NotImplementedException();

                //SpatialOperator op = program.FullOperator.ToSpatialOperator();

                // Problem: Not clear how to evaluate (needs a lot of
                // parameters which are only known to the time-stepper...)
                // Relation to residual logger?
            });

        /// <summary>
        /// The optional vorticity field:
        /// \f[ \vec{\omega} = \nabla \times \vec{u} \f]
        /// </summary>
        public static readonly Vector<DerivedVariable> Vorticity = new Vector<DerivedVariable>(
            d => new DerivedVariable(
                "omega" + d,
                VariableTypes.Other,
                delegate (DGField omega, CellMask cellMask, IProgram<CNSControl> program) {
                    throw new NotImplementedException();

                    // Problem: Order of updates is not defined, so relying on
                    // the fact that velocity has already been upated is EVIL

                    //if (Velocity == null) {
                    //    throw new Exception(
                    //        "Currently, computing vorticity requires calculating the velocity first");
                    //}

                    //switch (CNSEnvironment.NumberOfDimensions) {
                    //    case 1:
                    //        throw new Exception(
                    //            "The concept of vorticity does not make sense for"
                    //            + " one-dimensional flows");

                    //    case 2:
                    //        Debug.Assert(
                    //            Vorticity.Count == 1,
                    //            "In 2D, vorticity is scalar quantity");

                    //        Vorticity[0].Clear();
                    //        Vorticity[0].Curl2D(1.0, Velocity);
                    //        break;

                    //    case 3:
                    //        Vorticity.Clear();
                    //        Vorticity.Curl3D(1.0, Velocity);
                    //        break;

                    //    default:
                    //        throw new Exception();
                    //}
                }));

        /// <summary>
        /// The MPI rank of the current process
        /// </summary>
        public static readonly DerivedVariable Rank = new DerivedVariable(
            "mpiRank",
            VariableTypes.Other,
            delegate (DGField rank, CellMask cellMask, IProgram<CNSControl> program) {
                // Ignore cell mask, we want it everywhere!
                foreach (Chunk chunk in CellMask.GetFullMask(program.GridData)) {
                    foreach (int cell in chunk.Elements) {
                        rank.SetMeanValue(cell, program.GridData.MpiRank);
                    }
                }
            });

        /// <summary>
        /// The local non-dimensional kinemtic viscosity
        /// </summary>
        public static readonly DerivedVariable Viscosity = new DerivedVariable(
            "nu",
            VariableTypes.Other,
            delegate (DGField nu, CellMask cellMask, IProgram<CNSControl> program) {
                foreach (Chunk chunk in cellMask) {
                    foreach (int cell in chunk.Elements) {
                        nu.SetMeanValue(cell, program.Control.ViscosityLaw.GetViscosity(0.0, cell));
                    }
                }
            });

        /// <summary>
        /// The local non-dimensional artifical viscosity
        /// </summary>
        /// 
        //######################################################################
        //IMPORTANT: UPDATE ONLY POSSIBLE AFTER SENSOR FIELD WAS UPDATED/CREATED
        //depends on the order of the variables in the variable list
        //######################################################################
        private static SpecFemBasis avSpecFEMBasis;
        public static readonly DerivedVariable ArtificialViscosity = new DerivedVariable(
            "artificialViscosity",
            VariableTypes.Other,
            delegate (DGField artificialViscosity, CellMask cellMask, IProgram<CNSControl> program) {
                ConventionalDGField avField = artificialViscosity as ConventionalDGField;
                int D = cellMask.GridData.SpatialDimension;

                // Determine piecewise constant viscosity
                avField.Clear();
                foreach (Chunk chunk in cellMask) {
                    // Compute local mean state in cell
                    MultidimensionalArray meanValues = MultidimensionalArray.Create(chunk.Len, D + 2);
                    program.WorkingSet.Density.EvaluateMean(
                        chunk.i0,
                        chunk.Len,
                        meanValues.ExtractSubArrayShallow(-1, 0));
                    for (int d = 0; d < D; d++) {
                        program.WorkingSet.Momentum[d].EvaluateMean(
                            chunk.i0,
                            chunk.Len,
                            meanValues.ExtractSubArrayShallow(-1, d + 1));
                    }
                    program.WorkingSet.Energy.EvaluateMean(
                        chunk.i0,
                        chunk.Len,
                        meanValues.ExtractSubArrayShallow(-1, D + 1));

                    for (int i = 0; i < chunk.Len; i++) {
                        int cell = chunk.i0 + i;
                        StateVector state = new StateVector(
                            meanValues.ExtractSubArrayShallow(i, -1).To1DArray(),
                            program.SpeciesMap.GetMaterial(double.NaN));

                        double localViscosity = program.Control.ArtificialViscosityLaw.GetViscosity(
                           cell, program.GridData.Cells.h_min[cell], state);
                        avField.SetMeanValue(cell, localViscosity);
                    }
                }

                // Project visocsity onto continuous, multilinear space
                if (D < 3) {
                    // Standard version
                    if (avSpecFEMBasis == null || !avField.Basis.Equals(avSpecFEMBasis.ContainingDGBasis)) {
                        avSpecFEMBasis = new SpecFemBasis(program.GridData, 2);
                    }
                    SpecFemField specFemField = new SpecFemField(avSpecFEMBasis);
                    specFemField.ProjectDGFieldMaximum(1.0, avField);
                    avField.Clear();
                    specFemField.AccToDGField(1.0, avField);
                } else {
                    //MultidimensionalArray verticeCoordinates = MultidimensionalArray.Create(
                    //    NoOfCells, verticesPerCell, dimension);
                    //context.TransformLocal2Global(
                    //    localVerticeCoordinates,
                    //    cnk.i0,
                    //    cnk.Len,
                    //    verticeCoordinates,
                    //    cnt);
                    //PlotDriver.ZoneDriver.InitializeVertice2(
                    //    ;

                    if (program.GridData.MpiSize > 1) {
                        throw new NotImplementedException();
                    }

                    // Version that should finally also work in 3D
                    RefElement refElement = program.Grid.RefElements[0];
                    int N = program.GridData.Cells.NoOfLocalUpdatedCells;
                    int V = refElement.NoOfVertices;

                    // Sample maximum at vertices
                    MultidimensionalArray vertexValues = avField.Evaluate(refElement.Vertices);
                    for (int cell = 0; cell < N; cell++) {
                        for (int node = 0; node < vertexValues.GetLength(1); node++) {
                            double x = vertexValues[cell, node, 0];
                            double y = vertexValues[cell, node, 1];
                            double z = vertexValues[cell, node, 2];
                            double value = vertexValues[cell, node, D];

                            for (int coCell = 0; coCell < cell; coCell++) {
                                for (int coNode = 0; coNode < vertexValues.GetLength(1); coNode++) {
                                    double coX = vertexValues[coCell, coNode, 0];
                                    double coY = vertexValues[coCell, coNode, 1];
                                    double coZ = vertexValues[coCell, coNode, 2];

                                    if (Math.Abs(x - coX) < 1e-15
                                        && Math.Abs(y - coY) < 1e-15
                                        && Math.Abs(z - coZ) < 1e-15) {
                                        // Same point in space!

                                        double coValue = vertexValues[coCell, coNode, D];
                                        double max = Math.Max(value, coValue);
                                        vertexValues[cell, node, D] = max;
                                        vertexValues[coCell, coNode, D] = max;
                                        break;
                                    }
                                }
                            }
                        }

                        // Nodal projection using tri-linear coefficients only
                        avField.Clear();
                        ScalarFunctionEx func = delegate (int j0, int Len, NodeSet nodes, MultidimensionalArray result) {
                            for (int i = j0; i < Len; i++) {
                                for (int node = 0; node < vertexValues.GetLength(1); node++) {
                                    result[i, node] = vertexValues[i, node, D];
                                }
                            }
                        };
                        avField.ProjectNodalMultilinear(1.0, func, refElement.Vertices);
                    }
                }
            });

        /// <summary>
        /// The local sensor value of a shock sensor
        /// </summary>
        public static readonly DerivedVariable Sensor = new DerivedVariable(
            "sensor",
            VariableTypes.Other,
            delegate (DGField s, CellMask cellMask, IProgram<CNSControl> program) {
                IShockSensor sensor = program.Control.ShockSensor;
                foreach (Chunk chunk in cellMask) {
                    foreach (int cell in chunk.Elements) {
                        s.SetMeanValue(cell, sensor.GetSensorValue(cell));
                    }
                }
            });

        /// <summary>
        /// The clusters when using local time stepping
        /// </summary>
        public static readonly DerivedVariable LTSClusters = new DerivedVariable(
            "clusterLTS",
            VariableTypes.Other,
            delegate (DGField ClusterVisualizationField, CellMask cellMask, IProgram<CNSControl> program) {
                AdamsBashforthLTS LTSTimeStepper = (AdamsBashforthLTS)program.TimeStepper;

                if (LTSTimeStepper != null) {
                    for (int i = 0; i < LTSTimeStepper.CurrentClustering.NumberOfClusters; i++) {
                        SubGrid currentCluster = LTSTimeStepper.CurrentClustering.Clusters[i];
                        for (int j = 0; j < currentCluster.LocalNoOfCells; j++) {
                            foreach (Chunk chunk in currentCluster.VolumeMask) {
                                foreach (int cell in chunk.Elements) {
                                    ClusterVisualizationField.SetMeanValue(cell, i);
                                }
                            }
                        }
                    }
                }
            });
    }
}
