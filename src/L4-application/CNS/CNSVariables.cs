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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.Timestepping;
using CNS.Convection;
using CNS.ShockCapturing;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Diagnostics;
using System.Linq;
using static BoSSS.Foundation.Grid.Classic.GridData;
using static BoSSS.Solution.CompressibleFlowCommon.Variable;

namespace CNS {

    /// <summary>
    /// Common variables that are often used in control files. Note that
    /// additional custom variables (e.g., for debugging purposes) can
    /// analogously be specified in the control file itself.
    /// </summary>
    public static class CNSVariables {

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
        /// <see cref="BoSSS.Solution.CompressibleFlowCommon.MaterialProperty.IEquationOfState.GetEntropy"/>).
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
            delegate (DGField cfl, CellMask cellMask, IProgram<CNSControl> program) {
                // CFL sizes cannnot be plotted for the initial step, as the operator does not exist yet
                if (program.FullOperator == null) {
                    return;
                }

                // Query each cell individually so we get local results
                int J = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                for (int i = 0; i < J; i++) {
                    // Use "harmonic sum" of individual step sizes - see ExplicitEuler
                    double localCFL = 1.0 / program.FullOperator.CFLConstraints.Sum(c => 1.0 / (program.Control.CFLFraction * c.GetLocalStepSize(i, 1)));
                    cfl.SetMeanValue(i, localCFL);
                }
            });

        /// <summary>
        /// The local convective CFL restriction
        /// </summary>
        public static readonly DerivedVariable CFLConvective = new DerivedVariable(
            "cflConvective",
            VariableTypes.Other,
            delegate (DGField cfl, CellMask cellMask, IProgram<CNSControl> program) {
                // CFL sizes cannnot be plotted for the initial step, as the operator does not exist yet
                if (program.FullOperator == null) {
                    return;
                }

                TimeStepConstraint cflConstraint = program.FullOperator.CFLConstraints.OfType<ConvectiveCFLConstraint>().Single();

                // Query each cell individually so we get local results
                int J = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                for (int i = 0; i < J; i++) {
                    double localCFL = program.Control.CFLFraction * cflConstraint.GetLocalStepSize(i, 1);
                    cfl.SetMeanValue(i, localCFL);
                }
            });

        /// <summary>
        /// The local diffusive CFL restriction based on artificial viscosity
        /// </summary>
        public static readonly DerivedVariable CFLArtificialViscosity = new DerivedVariable(
            "cflArtificialViscosity",
            VariableTypes.Other,
            delegate (DGField cfl, CellMask cellMask, IProgram<CNSControl> program) {
                // CFL sizes cannnot be plotted for the initial step, as the operator does not exist yet
                if (program.FullOperator == null) {
                    return;
                }

                TimeStepConstraint cflConstraint = program.FullOperator.CFLConstraints.OfType<ArtificialViscosityCFLConstraint>().Single();

                // Query each cell individually so we get local results
                int J = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                for (int i = 0; i < J; i++) {
                    double localCFL = program.Control.CFLFraction * cflConstraint.GetLocalStepSize(i, 1);
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

                    //switch (CompressibleEnvironment.NumberOfDimensions) {
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
                foreach (int cell in cellMask.ItemEnum) {
                    nu.SetMeanValue(cell, program.Control.ViscosityLaw.GetViscosity(0.0, cell));
                }
            });

        /// <summary>
        /// The local non-dimensional specific enthalpy which is given by
        /// \f[h = \frac{\rho E + p}{\rho} \f], where \f[\rho E\f] is the total energy,
        /// \f[p\f] is the pressure, and \f[\rho\f] is the density.
        /// </summary>
        public static readonly DerivedVariable Enthalpy = new DerivedVariable(
            "h",
            VariableTypes.Other,
            delegate (DGField h, CellMask cellMask, IProgram<CNSControl> program) {
                h.Clear();
                h.ProjectFunction(
                    1.0,
                    (X, U, j) => new StateVector(U, program.SpeciesMap.GetMaterial(double.NaN)).Enthalpy,
                    new CellQuadratureScheme(true, cellMask),
                    program.WorkingSet.ConservativeVariables);
            });

        /// <summary>
        /// The local non-dimensional artificial viscosity
        /// </summary>
        /// <remarks>
        /// IMPORTANT: UPDATE ONLY POSSIBLE AFTER SENSOR FIELD WAS UPDATED/CREATED
        /// </remarks>
        private static SpecFemBasis avSpecFEMBasis;
        //private static Basis avContinuousDGBasis;
        public static readonly DerivedVariable ArtificialViscosity = new DerivedVariable(
            "artificialViscosity",
            VariableTypes.Other,
            delegate (DGField artificialViscosity, CellMask cellMask, IProgram<CNSControl> program) {
                using (var tr = new FuncTrace()) {
                    ConventionalDGField avField = artificialViscosity as ConventionalDGField;
                    int D = cellMask.GridData.SpatialDimension;
                    CellData cells = ((BoSSS.Foundation.Grid.Classic.GridData)program.GridData).Cells;

                    // h = shortest edge (not valid for true cut cells --> volume/surface)
                    var h_min = cells.h_min;

                    // h = volume / surface --> trying this for curved elements 
                    //var h_min = cells.h_min;
                    //h_min.Clear();
                    //foreach (int j in cellMask.ItemEnum) {
                    //    if (!cells.IsCellAffineLinear(j)) {
                    //        h_min[j] = cells.GetCellVolume(j) / cells.CellSurfaceArea[j];
                    //    }
                    //}

                    //var h_min = cells.CellLengthScale.CloneAs();
                    //h_min.ApplyAll(delegate (double x) {
                    //    return x * 2;
                    //});

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
                               cell, h_min[cell], state);

                            Debug.Assert(localViscosity >= 0.0);

                            avField.SetMeanValue(cell, localViscosity);
                        }
                    }


                    // Project visocsity onto continuous, multilinear space
                    if (D < 3) {
                        // Standard version
                        if (avSpecFEMBasis == null || !avField.Basis.Equals(avSpecFEMBasis.ContainingDGBasis)) {
                            avSpecFEMBasis = new SpecFemBasis((BoSSS.Foundation.Grid.Classic.GridData)program.GridData, 2);
                        }
                        SpecFemField specFemField = new SpecFemField(avSpecFEMBasis);
                        specFemField.ProjectDGFieldMaximum(1.0, avField);
                        avField.Clear();
                        specFemField.AccToDGField(1.0, avField);
                        avField.Clear(CellMask.GetFullMask(program.GridData).Except(program.SpeciesMap.SubGrid.VolumeMask));

                        // Continuous DG version
                        //Basis continuousDGBasis = new Basis(program.GridData, 2);
                        //ContinuousDGField continuousDGField = new ContinuousDGField(continuousDGBasis);
                        //continuousDGField.ProjectDGField(1.0, avField);
                        //avField.Clear();
                        //continuousDGField.AccToDGField(1.0, avField);
                        //avField.Clear(CellMask.GetFullMask(program.GridData).Except(program.SpeciesMap.SubGrid.VolumeMask));
                    } else {
                        if (program.GridData.MpiSize > 1) {
                            throw new NotImplementedException();
                        }


                        // Version that should finally also work in 3D
                        RefElement refElement = ((BoSSS.Foundation.Grid.Classic.GridCommons)(program.Grid)).RefElements[0];
                        int N = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
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
                }
            });

        /// <summary>
        /// The local sensor value of a shock sensor
        /// </summary>
        public static readonly DerivedVariable ShockSensor = new DerivedVariable(
            "sensor",
            VariableTypes.Other,
            delegate (DGField s, CellMask cellMask, IProgram<CNSControl> program) {
                s.Clear();  // Just to be sure that we get new values only

                ICNSShockSensor sensor = program.Control.CNSShockSensor;
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
                // Don't fail, just ignore
                if (LTSTimeStepper == null) {
                    return;
                }

                for (int i = 0; i < LTSTimeStepper.CurrentClustering.NumberOfClusters; i++) {
                    SubGrid currentCluster = LTSTimeStepper.CurrentClustering.Clusters[i];
                    //for (int j = 0; j < currentCluster.LocalNoOfCells; j++) {
                    foreach (int cell in currentCluster.VolumeMask.ItemEnum) {
                        //foreach (Chunk chunk in currentCluster.VolumeMask) {
                        //    foreach (int cell in chunk.Elements) {
                        ClusterVisualizationField.SetMeanValue(cell, i);
                        //    }
                    }
                    //}
                }
            });

        /// <summary>
        /// The so-called Schlieren variables is based on the magnitude of the density gradient
        /// The implementation is based on the definition for the shock-vortex interaction test case
        /// of the HiOCFD5 workshop
        /// </summary>
        public static readonly DerivedVariable Schlieren = new DerivedVariable(
            "schlieren",
            VariableTypes.Other,
            delegate (DGField schlierenField, CellMask cellMask, IProgram<CNSControl> program) {
                schlierenField.Clear();
                int D = program.GridData.SpatialDimension;

                // Old version by Björn
                /*
                // Calculate the magnitude of the density gradient
                SinglePhaseField derivative = new SinglePhaseField(schlierenField.Basis, "derivative");

                for (int d = 0; d < D; d++) {
                    derivative.Derivative(1.0, program.WorkingSet.Density, d);
                    foreach (int cell in cellMask.ItemEnum) {
                        double updateValue = schlierenField.GetMeanValue(cell) + Math.Pow(derivative.GetMeanValue(cell), 2);
                        if (d == (D - 1)) {
                            schlierenField.SetMeanValue(cell, Math.Sqrt(updateValue));
                        } else {
                            schlierenField.SetMeanValue(cell, updateValue);
                        }

                    }
                }
                */

                // New version by Florian
                DGField Density = program.WorkingSet.Density;
                VectorField<SinglePhaseField> DensityGradient = new VectorField<SinglePhaseField>(D, Density.Basis, (b, S) => new SinglePhaseField(b, S));
                DensityGradient.Gradient(1.0, Density, cellMask);

                schlierenField.ProjectFunction(1.0,
                    delegate (double[] X, double[] U, int jCell) {
                        double R;
                        R = 1.0;
                        if (D == 2) {
                            R += Math.Sqrt(U[0] * U[0] + U[1] * U[1]);
                        } else if (D == 3) {
                            R += Math.Sqrt(U[0] * U[0] + U[1] * U[1] + U[2] * U[2]);
                        } else {
                            throw new NotSupportedException();
                        }
                        R = Math.Log(R) / Math.Log(10);
                        return R;
                    },
                    new CellQuadratureScheme(true, cellMask),
                    DensityGradient.ToArray());

            });

        public static readonly DerivedVariable EmptyField = new DerivedVariable(
            "emptyField",
            VariableTypes.Other,
            delegate (DGField field, CellMask cellMask, IProgram<CNSControl> program) {
            });

        public static readonly DerivedVariable RiemannDensity = new DerivedVariable(
            "riemannDensity",
            VariableTypes.Other,
            delegate (DGField field, CellMask cellMask, IProgram<CNSControl> program) {
                if (program.TimeStepper != null) {
                    // ### Initial conditions ###
                    double densityLeft = 1.0;
                    double densityRight = 0.125;
                    double pressureLeft = 1.0;
                    double pressureRight = 0.1;
                    double velocityLeft = 0.0;
                    double velocityRight = 0.0;
                    double discontinuityPosition = 0.5;

                    var control = program.Control;
                    Material material = new Material(control.EquationOfState, control.ViscosityLaw, control.MachNumber, control.ReynoldsNumber, control.PrandtlNumber, control.FroudeNumber, control.ViscosityRatio);
                    StateVector stateLeft = StateVector.FromPrimitiveQuantities(
                        material, densityLeft, new Vector(velocityLeft, 0.0, 0.0), pressureLeft);
                    StateVector stateRight = StateVector.FromPrimitiveQuantities(
                        material, densityRight, new Vector(velocityRight, 0.0, 0.0), pressureRight);

                    var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0, 0.0, 0.0));
                    double pStar, uStar;
                    riemannSolver.GetStarRegionValues(out pStar, out uStar);
                    field.Clear();
                    field.ProjectFunction(1.0,
                        (X, U, j) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, program.TimeStepper.Time).Density,
                        new CellQuadratureScheme(true, cellMask),
                        program.WorkingSet.ConservativeVariables
                        );
                }
            });
    }
}
