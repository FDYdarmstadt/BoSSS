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
using ilPSP;
using System;
using static BoSSS.Solution.CompressibleFlowCommon.Variable;

namespace XESF.Variables {
    public static class XESFVariables {

        /// <summary>
        /// The obligatory level set
        /// </summary>
        /// <summary>
        /// The obligatory level set
        /// </summary>
        public static readonly Variable LevelSet = new Variable("levelSet", VariableTypes.Other);

        public static readonly Variable LevelSetTwo = new Variable("levelSetTwo", VariableTypes.Other);


        public static readonly DerivedVariable<XDGField> CutCells = new DerivedVariable<XDGField>(
            "cutCells",
            VariableTypes.Other,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();

                CellMask cutCellMask = program.LsTrk.Regions.GetCutCellMask();

                foreach(Chunk chunk in cutCellMask) {
                    foreach(int cell in chunk.Elements) {
                        foreach(SpeciesId id in program.LsTrk.SpeciesIdS) {
                            dgField.GetSpeciesShadowField(id).SetMeanValue(cell, 1);
                        }
                    }
                }
            });

        public static readonly DerivedVariable<XDGField> CutCellsWithOutSourceCells = new DerivedVariable<XDGField>(
            "cutCellsWithoutSourceCells",
            VariableTypes.Other,
            delegate (XDGField dgField, XESFMain program) {
                // Do not do this in the first time step, because the MultiphaseAgglomerator does not exist yet
                if(program.XSpatialOperator == null) {
                    return;
                }

                dgField.Clear();

                CellMask cutCellMask = program.LsTrk.Regions.GetCutCellMask();

                foreach(SpeciesId id in program.SpeciesToEvaluate_Ids) {
                    foreach(Chunk chunk in cutCellMask.Except(program.MultiphaseAgglomerator.GetAgglomerator(id).AggInfo.SourceCells)) {
                        foreach(int cell in chunk.Elements) {
                            dgField.GetSpeciesShadowField(id).SetMeanValue(cell, 1);
                        }
                    }
                }
            });







        public static readonly DerivedVariable<double> residual = new DerivedVariable<double>("Residual",
            VariableTypes.Other,
            delegate (double residual, XESFMain program) {

                //CellMask cellMask = program.LsTrk.Regions.GetCutCellMask();
                LevelSetTracker LsTrk = program.LsTrk;
                XDGBasis[] basis = new XDGBasis[program.ConservativeFields.Length];
                XDGBasis tmp_basis;
                //
                for(int i = 0; i < program.ConservativeFields.Length; i++) {
                    tmp_basis = new XDGBasis(LsTrk, program.ConservativeFields[i].Basis.Degree);
                    basis[i] = tmp_basis;
                }
                var mapping2 = new CoordinateMapping(program.ConservativeFields);
                var eval2 = program.XSpatialOperator.GetEvaluatorEx(LsTrk, program.ConservativeFields, null, mapping2);
                double[] residual_vector2 = new double[mapping2.Ntotal];
                eval2.Evaluate(1, 1, residual_vector2);


                residual = residual_vector2.L2Norm() / 2;

                //Console.WriteLine("original residual = " + residual);


            });
        /// <summary>
        /// The local sensor value of a shock sensor
        /// </summary>
        //public static readonly DerivedVariable<SinglePhaseField> Sensor = new DerivedVariable<SinglePhaseField>(
        //    "sensor",
        //    VariableTypes.Other,
        //    delegate (SinglePhaseField dgField, XESFMain program) {
        //        dgField.Clear();

        //        // Should be removed later
        //        //Console.WriteLine("sensor L2-norm: " + sensorField.L2Norm());

        //        dgField.Acc(1.0, program.Sensor.GetSensorSinglePhaseField());
        //    });

        ///// <summary>
        ///// XDG artificial viscosity
        ///// </summary>
        //public static readonly DerivedVariable<SinglePhaseField> ArtificialViscosity = new DerivedVariable<SinglePhaseField>(
        //    "artificialViscosity",
        //    VariableTypes.Other,
        //    delegate (SinglePhaseField dgField, XESFMain program) {
        //        dgField.Clear();

        //        dgField.Acc(1.0, program.ArtificialViscosityField);
        //    });

        /// <summary>
        /// The optional pressure field:
        /// \f[ p = f(\rho, \vec{m}, \rho E) \f], where
        /// \f[ f \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable<XDGField> Pressure = new DerivedVariable<XDGField>(
            "p",
            VariableTypes.Pressure,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);
                //foreach(string species in program.Control.SpeciesToEvaluate) {
                //    CellMask cellMask = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(species));

        //    IList<DGField> fields = new List<DGField>();
        //    foreach(XDGField field in program.ConservativeFields) {
        //        fields.Add(field.GetSpeciesShadowField(species));
        //    }

        //    dgField.GetSpeciesShadowField(species).ProjectFunction(
        //        1.0,
        //        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
        //        new CellQuadratureScheme(true, cellMask),
        //        fields.ToArray());
        //    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
        //}
            }
            );
        /// <summary>
        /// The optional pressure- step -field :
        /// \f[ p = f(\rho, \vec{m}, \rho E) \f], where
        /// \f[ f \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable<XDGField> PressureStep = new DerivedVariable<XDGField>(
            "p_step",
            VariableTypes.Pressure,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                    program.SpeciesToEvaluate_Ids,
                    program.originalStep);
                //foreach(string species in program.Control.SpeciesToEvaluate) {
                //    CellMask cellMask = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(species));

                //    IList<DGField> fields = new List<DGField>();
                //    foreach(XDGField field in program.originalStep) {
                //        fields.Add(field.GetSpeciesShadowField(species));
                //    }

                //    dgField.GetSpeciesShadowField(species).ProjectFunction(
                //        1.0,
                //        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                //        new CellQuadratureScheme(true, cellMask),
                //        fields.ToArray());
                //    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                //}
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
                delegate (XDGField dgField, XESFMain program) {
                    dgField.Clear();
                                    dgField.ProjectFunctionXDG(1.0,
                        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Velocity[d],
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);
                    //foreach(string species in program.Control.SpeciesToEvaluate) {
                    //    CellMask cellMask = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(species));
                    //    dgField.GetSpeciesShadowField(species).ProjectQuotient(
                    //        1.0,
                    //        program.Momentum[d].GetSpeciesShadowField(species),
                    //        program.Density.GetSpeciesShadowField(species),
                    //        cellMask,
                    //        accumulateResult: true);
                    //}
                }
                ));

        /// <summary>
        /// The local non-dimensional specific enthalpy which is given by
        /// \f[h = \frac{\rho E + p}{\rho} \f], where \f[\rho E\f] is the total energy,
        /// \f[p\f] is the pressure, and \f[\rho\f] is the density.
        /// </summary>
        public static readonly DerivedVariable<XDGField> Enthalpy = new DerivedVariable<XDGField>(
            "h",
            VariableTypes.Other,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Enthalpy, 
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);
                //foreach(string species in program.Control.SpeciesToEvaluate) {
                //    IList<DGField> fields = new List<DGField>();
                //    foreach(XDGField field in program.ConservativeFields) {
                //        fields.Add(field.GetSpeciesShadowField(species));
                //    }

                //    CellMask cellMask = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(species));

                //    dgField.GetSpeciesShadowField(species).ProjectFunction(
                //        1.0,
                //        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Enthalpy,
                //        new CellQuadratureScheme(true, cellMask),
                //        fields.ToArray());
                //}
            });
        public static readonly DerivedVariable<XDGField> Enthalpy_Error = new DerivedVariable<XDGField>(
           "h_err",
           VariableTypes.Other,
           delegate (XDGField dgField, XESFMain program) {
               dgField.Clear();
               dgField.ProjectFunctionXDG(1.0,
                       (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Enthalpy - program.Control.ExactEnthalpy,
                       program.SpeciesToEvaluate_Ids,
                       program.ConservativeFields);
               //foreach(string species in program.Control.SpeciesToEvaluate) {
               //    IList<DGField> fields = new List<DGField>();
               //    foreach(XDGField field in program.ConservativeFields) {
               //        fields.Add(field.GetSpeciesShadowField(species));
               //    }

               //    CellMask cellMask = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(species));

               //    dgField.GetSpeciesShadowField(species).ProjectFunction(
               //        1.0,
               //        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Enthalpy,
               //        new CellQuadratureScheme(true, cellMask),
               //        fields.ToArray());
               //}
           });
        /// <summary>
        /// The optional local Mach number field
        /// \f[ Ma = \frac{\Vert u \Vert}{a} \f], where the formula for the
        /// local speed of sound \f[ a \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable<XDGField> LocalMachNumber = new DerivedVariable<XDGField>(
            "Ma",
            VariableTypes.Other,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(
                        1.0,
                        delegate (Vector X, double[] U, int j) {
                            var state = new StateVector(U, program.Control.GetMaterial());
                            double r = state.LocalMachNumber;
                            //Debug.Assert(r.IsNaN() == false, "NAN in local Mach number");
                            //Debug.Assert(r.IsInfinity() == false, "INF in local Mach number");
                            return r;
                        },
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);
                //foreach(string species in program.Control.SpeciesToEvaluate) {
                //    IList<DGField> fields = new List<DGField>();
                //    foreach(XDGField field in program.ConservativeFields) {
                //        fields.Add(field.GetSpeciesShadowField(species));
                //    }

                //    CellMask cellMask = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(species));

                //    dgField.GetSpeciesShadowField(species).ProjectFunction(
                //        1.0,
                //        delegate (Vector X, double[] U, int j) {
                //            var state = new StateVector(U, program.Control.GetMaterial());
                //            double r = state.LocalMachNumber;
                //            Debug.Assert(r.IsNaN() == false, "NAN in local Mach number");
                //            Debug.Assert(r.IsInfinity() == false, "INF in local Mach number");
                //            return r;
                //        },
                //        new CellQuadratureScheme(true, cellMask),
                //        fields.ToArray());
                //}
            });

        /// <summary>
        /// The optional entropy field (see
        /// <see cref="BoSSS.Solution.CompressibleFlowCommon.MaterialProperty.IEquationOfState.GetEntropy"/>).
        /// </summary>
        public static readonly DerivedVariable<XDGField> Entropy = new DerivedVariable<XDGField>(
            "S",
            VariableTypes.Other,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Entropy,
                    program.SpeciesToEvaluate_Ids,
                    program.ConservativeFields);
            });


        /// <summary>
        /// The local convective CFL restriction
        /// </summary>
        //public static readonly DerivedVariable<XDGField> CFLConvective = new DerivedVariable<XDGField>(
        //    "cflConvective",
        //    VariableTypes.Other,
        //    delegate (XDGField dgField, XESFMain program) {
        //        // CFL sizes cannot be plotted for the initial step, as the operator does not exist yet
        //        if(program.XSpatialOperator == null) {
        //            return;
        //        }

        //        foreach(ConvectiveXDGCFLConstraint constraint in program.TimeStepConstraints) {
        //            // Query each cell individually so we get local results
        //            foreach(Chunk chunk in program.LsTrk.Regions.GetSpeciesSubGrid(constraint.FluidSpeciesId).VolumeMask) {
        //                foreach(int cell in chunk.Elements) {
        //                    double localCFL = program.Control.CFLFraction * constraint.GetLocalStepSize(cell, 1);
        //                    dgField.GetSpeciesShadowField(constraint.FluidSpeciesId).SetMeanValue(cell, localCFL);
        //                }
        //            }
        //        }
        //    });

        

        public static readonly DerivedVariable<SinglePhaseField> JumpNormAtShock = new DerivedVariable<SinglePhaseField>(
            "jump_rho",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XESFMain program) {
                dgField.Clear();

                var LsTrk = program.LsTrk;
                int D = LsTrk.GridDat.SpatialDimension;
                int quadDregree = program.Control.NonlinearQuadratureDegree;

                var field_L = program.Density.GetSpeciesShadowField("L");
                var field_R = program.Density.GetSpeciesShadowField("R");

                void ErrFunc(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of nodes
                    int _D = D; // local var may be a bit faster
                    MultidimensionalArray field_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray field_R_Res = MultidimensionalArray.Create(Len, K);

                    field_L.Evaluate(j0, Len, NS, field_L_Res);
                    field_R.Evaluate(j0, Len, NS, field_R_Res);

                    for(int j = 0; j < Len; j++) {
                        for(int k = 0; k < K; k++) {
                            double q = field_R_Res[j, k] - field_L_Res[j, k];
                            result[j, k] = q * q;
                        }
                    }
                }

                double Jump_NORM = 0.0;

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId("L") }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: 1, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: 1));

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, quadDregree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            Jump_NORM += ResultsOfIntegration[i, 0];
                            dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0].Sqrt());
                        }
                    }
                ).Execute();

                Console.WriteLine("JumpNormAtShock: " + Math.Sqrt(Jump_NORM));
            });


        public static readonly DerivedVariable<SinglePhaseField> rh_conti_indicator_JS = new DerivedVariable<SinglePhaseField>(
            "rh_conti_indicator_JS",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XESFMain program) {
                dgField.Clear();

                // Quadrature
                int quadDregree = program.Control.NonlinearQuadratureDegree;

                // Level set
                LevelSetTracker LsTrk = program.LsTrk;
                int D = LsTrk.GridDat.SpatialDimension;
                string left = program.Control.LsOne_NegSpecies;
                string right = program.Control.LsOne_PosSpecies;
                int levelSetIndex = 0;

                // Density (rho_1, rho_2)
                XDGField.SpeciesShadowField momentum_L = program.Momentum[0].GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField momentum_R = program.Momentum[0].GetSpeciesShadowField(right);


                CellMask cellMask_L = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(left));
                CellMask cellMask_R = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(right));


                // Evaluate RH-condition at the interface
                void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of nodes
                    int _D = D; // local var may be a bit faster

                    MultidimensionalArray momentum_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray momentum_R_Res = MultidimensionalArray.Create(Len, K);

                    momentum_L.Evaluate(j0, Len, NS, momentum_L_Res);
                    momentum_R.Evaluate(j0, Len, NS, momentum_R_Res);


                    // RH (first eq.), stationary case
                    for(int j = 0; j < Len; j++) {
                        for(int k = 0; k < K; k++) {
                            // right side - left side
                            double q = momentum_L_Res[j, k] - momentum_R_Res[j, k];
                            result[j, k] = q;
                        }
                    }
                }

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, quadDregree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                        }
                    }
                ).Execute();
            });
        public static readonly DerivedVariable<SinglePhaseField> rh_conti_indicator_LsOne = new DerivedVariable<SinglePhaseField>(
            "rh_conti_indicator_LsOne",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XESFMain program) {
                dgField.Clear();

                // Quadrature
                int quadDregree = program.Control.NonlinearQuadratureDegree;

                // Level set
                LevelSetTracker LsTrk = program.LsTrk;
                int D = LsTrk.GridDat.SpatialDimension;
                string left = program.Control.LsOne_NegSpecies;
                string right = program.Control.LsOne_PosSpecies;
                int levelSetIndex = 0;

                // Density (rho_1, rho_2)
                XDGField.SpeciesShadowField density_L = program.Density.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField density_R = program.Density.GetSpeciesShadowField(right);

                // Velocity (u_1, u_2)
                XDGField velocity = new XDGField(program.Density.Basis);
                XDGField.SpeciesShadowField velocity_L = velocity.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField velocity_R = velocity.GetSpeciesShadowField(right);

                CellMask cellMask_L = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(left));
                velocity_L.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(left),
                    density_L,
                    cellMask_L,
                    accumulateResult: true);

                CellMask cellMask_R = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(right));
                velocity_R.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(right),
                    density_R,
                    cellMask_R,
                    accumulateResult: true);

                // Evaluate RH-condition at the interface
                void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of nodes
                    int _D = D; // local var may be a bit faster

                    MultidimensionalArray density_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray density_R_Res = MultidimensionalArray.Create(Len, K);

                    MultidimensionalArray velocity_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray velocity_R_Res = MultidimensionalArray.Create(Len, K);

                    density_L.Evaluate(j0, Len, NS, density_L_Res);
                    density_R.Evaluate(j0, Len, NS, density_R_Res);

                    velocity_L.Evaluate(j0, Len, NS, velocity_L_Res);
                    velocity_R.Evaluate(j0, Len, NS, velocity_R_Res);

                    // RH (first eq.), stationary case
                    for(int j = 0; j < Len; j++) {
                        for(int k = 0; k < K; k++) {
                            // right side - left side
                            double q = density_L_Res[j, k] * velocity_L_Res[j, k] - density_R_Res[j, k] * velocity_R_Res[j, k];
                            result[j, k] = q;
                        }
                    }
                }

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, quadDregree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                        }
                    }
                ).Execute();
            });
        public static readonly DerivedVariable<SinglePhaseField> rh_conti_indicator_LsTwo = new DerivedVariable<SinglePhaseField>(
            "rh_conti_indicator_LsTwo",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XESFMain program) {
                dgField.Clear();

                // Quadrature
                int quadDregree = program.Control.NonlinearQuadratureDegree;

                // Level set
                LevelSetTracker LsTrk = program.LsTrk;
                int D = LsTrk.GridDat.SpatialDimension;
                string left = program.Control.LsTwo_NegSpecies;
                string right = program.Control.LsTwo_PosSpecies;
                int levelSetIndex = 1;

                // Density (rho_1, rho_2)
                XDGField.SpeciesShadowField density_L = program.Density.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField density_R = program.Density.GetSpeciesShadowField(right);

                // Velocity (u_1, u_2)
                XDGField velocity = new XDGField(program.Density.Basis);
                XDGField.SpeciesShadowField velocity_L = velocity.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField velocity_R = velocity.GetSpeciesShadowField(right);

                CellMask cellMask_L = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(left));
                velocity_L.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(left),
                    density_L,
                    cellMask_L,
                    accumulateResult: true);

                CellMask cellMask_R = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(right));
                velocity_R.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(right),
                    density_R,
                    cellMask_R,
                    accumulateResult: true);

                // Evaluate RH-condition at the interface
                void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of nodes
                    int _D = D; // local var may be a bit faster

                    MultidimensionalArray density_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray density_R_Res = MultidimensionalArray.Create(Len, K);

                    MultidimensionalArray velocity_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray velocity_R_Res = MultidimensionalArray.Create(Len, K);

                    density_L.Evaluate(j0, Len, NS, density_L_Res);
                    density_R.Evaluate(j0, Len, NS, density_R_Res);

                    velocity_L.Evaluate(j0, Len, NS, velocity_L_Res);
                    velocity_R.Evaluate(j0, Len, NS, velocity_R_Res);

                    // RH (first eq.), stationary case
                    for(int j = 0; j < Len; j++) {
                        for(int k = 0; k < K; k++) {
                            // right side - left side
                            double q = density_L_Res[j, k] * velocity_L_Res[j, k] - density_R_Res[j, k] * velocity_R_Res[j, k];
                            result[j, k] = q;
                        }
                    }
                }

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, quadDregree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                        }
                    }
                ).Execute();
            });

        public static readonly DerivedVariable<SinglePhaseField> rh_mom_indicator_LsOne = new DerivedVariable<SinglePhaseField>(
            "rh_mom_indicator_LsOne",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XESFMain program) {
                dgField.Clear();

                // Quadrature
                int quadDregree = program.Control.NonlinearQuadratureDegree;

                // Level set
                LevelSetTracker LsTrk = program.LsTrk;
                int D = LsTrk.GridDat.SpatialDimension;
                string left = program.Control.LsOne_NegSpecies;
                string right = program.Control.LsOne_PosSpecies;
                int levelSetIndex = 0;

                // Density (rho_1, rho_2)
                XDGField.SpeciesShadowField density_L = program.Density.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField density_R = program.Density.GetSpeciesShadowField(right);

                // Velocity (u_1, u_2)
                XDGField velocity = new XDGField(program.Density.Basis);
                XDGField.SpeciesShadowField velocity_L = velocity.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField velocity_R = velocity.GetSpeciesShadowField(right);

                CellMask cellMask_L = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(left));
                velocity_L.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(left),
                    density_L,
                    cellMask_L,
                    accumulateResult: true);

                CellMask cellMask_R = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(right));
                velocity_R.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(right),
                    density_R,
                    cellMask_R,
                    accumulateResult: true);

                // Evaluate RH-condition at the interface
                void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of nodes
                    int _D = D; // local var may be a bit faster

                    MultidimensionalArray density_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray density_R_Res = MultidimensionalArray.Create(Len, K);

                    MultidimensionalArray velocity_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray velocity_R_Res = MultidimensionalArray.Create(Len, K);

                    density_L.Evaluate(j0, Len, NS, density_L_Res);
                    density_R.Evaluate(j0, Len, NS, density_R_Res);

                    velocity_L.Evaluate(j0, Len, NS, velocity_L_Res);
                    velocity_R.Evaluate(j0, Len, NS, velocity_R_Res);

                    // RH (first eq.), stationary case
                    for(int j = 0; j < Len; j++) {
                        for(int k = 0; k < K; k++) {
                            // right side - left side
                            double q = density_L_Res[j, k] * (density_L_Res[j, k] * velocity_L_Res[j, k] * velocity_L_Res[j, k]) - density_R_Res[j, k] * (density_R_Res[j, k] * velocity_R_Res[j, k] * velocity_R_Res[j, k]);
                            result[j, k] = q;
                        }
                    }
                }

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, quadDregree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                        }
                    }
                ).Execute();
            });

        public static readonly DerivedVariable<SinglePhaseField> rh_mom_indicator_LsTwo = new DerivedVariable<SinglePhaseField>(
            "rh_mom_indicator_LsTwo",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XESFMain program) {
                dgField.Clear();

                // Quadrature
                int quadDregree = program.Control.NonlinearQuadratureDegree;

                // Level set
                LevelSetTracker LsTrk = program.LsTrk;
                int D = LsTrk.GridDat.SpatialDimension;
                string left = program.Control.LsTwo_NegSpecies;
                string right = program.Control.LsTwo_PosSpecies;
                int levelSetIndex = 1;

                // Density (rho_1, rho_2)
                XDGField.SpeciesShadowField density_L = program.Density.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField density_R = program.Density.GetSpeciesShadowField(right);

                // Velocity (u_1, u_2)
                XDGField velocity = new XDGField(program.Density.Basis);
                XDGField.SpeciesShadowField velocity_L = velocity.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField velocity_R = velocity.GetSpeciesShadowField(right);

                CellMask cellMask_L = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(left));
                velocity_L.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(left),
                    density_L,
                    cellMask_L,
                    accumulateResult: true);

                CellMask cellMask_R = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(right));
                velocity_R.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(right),
                    density_R,
                    cellMask_R,
                    accumulateResult: true);

                // Evaluate RH-condition at the interface
                void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of nodes
                    int _D = D; // local var may be a bit faster

                    MultidimensionalArray density_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray density_R_Res = MultidimensionalArray.Create(Len, K);

                    MultidimensionalArray velocity_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray velocity_R_Res = MultidimensionalArray.Create(Len, K);

                    density_L.Evaluate(j0, Len, NS, density_L_Res);
                    density_R.Evaluate(j0, Len, NS, density_R_Res);

                    velocity_L.Evaluate(j0, Len, NS, velocity_L_Res);
                    velocity_R.Evaluate(j0, Len, NS, velocity_R_Res);

                    // RH (first eq.), stationary case
                    for(int j = 0; j < Len; j++) {
                        for(int k = 0; k < K; k++) {
                            // right side - left side
                            double q = density_L_Res[j, k] * (density_L_Res[j, k] * velocity_L_Res[j, k] * velocity_L_Res[j, k]) - density_R_Res[j, k] * (density_R_Res[j, k] * velocity_R_Res[j, k] * velocity_R_Res[j, k]);
                            result[j, k] = q;
                        }
                    }
                }

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, quadDregree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                        }
                    }
                ).Execute();
            });

        public static readonly DerivedVariable<SinglePhaseField> rh_energy_indicator_LsOne = new DerivedVariable<SinglePhaseField>(
            "rh_energy_indicator_LsOne",
            VariableTypes.Other,
            delegate (SinglePhaseField dgField, XESFMain program) {
                dgField.Clear();

                // Quadrature
                int quadDregree = program.Control.NonlinearQuadratureDegree;

                // Level set
                LevelSetTracker LsTrk = program.LsTrk;
                int D = LsTrk.GridDat.SpatialDimension;
                string left = program.Control.LsOne_NegSpecies;
                string right = program.Control.LsOne_PosSpecies;
                int levelSetIndex = 0;

                // Density (rho_1, rho_2)
                XDGField.SpeciesShadowField density_L = program.Density.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField density_R = program.Density.GetSpeciesShadowField(right);

                // Velocity (u_1, u_2)
                XDGField velocity = new XDGField(program.Density.Basis);
                XDGField.SpeciesShadowField velocity_L = velocity.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField velocity_R = velocity.GetSpeciesShadowField(right);

                CellMask cellMask_L = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(left));
                velocity_L.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(left),
                    density_L,
                    cellMask_L,
                    accumulateResult: true);

                CellMask cellMask_R = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(right));
                velocity_R.ProjectQuotient(
                    1.0,
                    program.Momentum[0].GetSpeciesShadowField(right),
                    density_R,
                    cellMask_R,
                    accumulateResult: true);

                // Energy (rho_E_1, rho_E_2 )
                XDGField.SpeciesShadowField energy_L = program.Energy.GetSpeciesShadowField(left);
                XDGField.SpeciesShadowField energy_R = program.Energy.GetSpeciesShadowField(right);

                // Evaluate RH-condition at the interface
                void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No of nodes
                    int _D = D; // local var may be a bit faster

                    MultidimensionalArray density_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray density_R_Res = MultidimensionalArray.Create(Len, K);

                    MultidimensionalArray velocity_L_Res = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray velocity_R_Res = MultidimensionalArray.Create(Len, K);

                    density_L.Evaluate(j0, Len, NS, density_L_Res);
                    density_R.Evaluate(j0, Len, NS, density_R_Res);

                    velocity_L.Evaluate(j0, Len, NS, velocity_L_Res);
                    velocity_R.Evaluate(j0, Len, NS, velocity_R_Res);

                    // RH (first eq.), stationary case
                    for(int j = 0; j < Len; j++) {
                        for(int k = 0; k < K; k++) {
                            // right side - left side
                            double q = density_L_Res[j, k] * (density_L_Res[j, k] * velocity_L_Res[j, k] * velocity_L_Res[j, k]) - density_R_Res[j, k] * (density_R_Res[j, k] * velocity_R_Res[j, k] * velocity_R_Res[j, k]);
                            result[j, k] = q;
                        }
                    }
                }

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, quadDregree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                        }
                    }
                ).Execute();
            });

        // not implemented
        //public static readonly DerivedVariable<SinglePhaseField> rh_energy_indicator_LsTwo = new DerivedVariable<SinglePhaseField>(
        //    "rh_mom_indicator_LsTwo",
        //    VariableTypes.Other,
        //    delegate (SinglePhaseField dgField, XESFMain program) {
        //        dgField.Clear();

        //        // Quadrature
        //        int quadDregree = program.Control.NonlinearQuadratureDegree;

        //        // Level set
        //        LevelSetTracker LsTrk = program.LsTrk;
        //        int D = LsTrk.GridDat.SpatialDimension;
        //        string left = program.Control.LsTwo_NegSpecies;
        //        string right = program.Control.LsTwo_PosSpecies;
        //        int levelSetIndex = 0;

        //        // Density (rho_1, rho_2)
        //        XDGField.SpeciesShadowField density_L = program.Density.GetSpeciesShadowField(left);
        //        XDGField.SpeciesShadowField density_R = program.Density.GetSpeciesShadowField(right);

        //        // Velocity (u_1, u_2)
        //        XDGField velocity = new XDGField(program.Density.Basis);
        //        XDGField.SpeciesShadowField velocity_L = velocity.GetSpeciesShadowField(left);
        //        XDGField.SpeciesShadowField velocity_R = velocity.GetSpeciesShadowField(right);

        //        CellMask cellMask_L = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(left));
        //        velocity_L.ProjectQuotient(
        //            1.0,
        //            program.Momentum[0].GetSpeciesShadowField(left),
        //            density_L,
        //            cellMask_L,
        //            accumulateResult: true);

        //        CellMask cellMask_R = program.LsTrk.Regions.GetSpeciesMask(program.LsTrk.GetSpeciesId(right));
        //        velocity_R.ProjectQuotient(
        //            1.0,
        //            program.Momentum[0].GetSpeciesShadowField(right),
        //            density_R,
        //            cellMask_R,
        //            accumulateResult: true);


        //        // Energy (rho_E_1, rho_E_2 )
        //        XDGField.SpeciesShadowField energy_L = program.Energy.GetSpeciesShadowField(left);
        //        XDGField.SpeciesShadowField energy_R = program.Energy.GetSpeciesShadowField(right);

        //        // Evaluate RH-condition at the interface
        //        void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
        //            int K = result.GetLength(1); // No of nodes
        //            int _D = D; // local var may be a bit faster

        //            MultidimensionalArray energy_L_Res = MultidimensionalArray.Create(Len, K);
        //            MultidimensionalArray energy_R_Res = MultidimensionalArray.Create(Len, K);

        //            MultidimensionalArray velocity_L_Res = MultidimensionalArray.Create(Len, K);
        //            MultidimensionalArray velocity_R_Res = MultidimensionalArray.Create(Len, K);

        //            energy_L.Evaluate(j0, Len, NS, energy_L_Res);
        //            energy_R.Evaluate(j0, Len, NS, energy_R_Res);

        //            velocity_L.Evaluate(j0, Len, NS, velocity_L_Res);
        //            velocity_R.Evaluate(j0, Len, NS, velocity_R_Res);

        //            // RH (first eq.), stationary case
        //            for(int j = 0; j < Len; j++) {
        //                for(int k = 0; k < K; k++) {
        //                    // right side - left side
        //                    //double q = density_L_Res[j, k] * (density_L_Res[j, k] * velocity_L_Res[j, k] * velocity_L_Res[j, k]) - density_R_Res[j, k] * (density_R_Res[j, k] * velocity_R_Res[j, k] * velocity_R_Res[j, k]);
        //                    //result[j, k] = q;
        //                }
        //            }
        //        }

        //        var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
        //        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

        //        CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
        //            cqs.Compile(LsTrk.GridDat, quadDregree),
        //            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
        //                RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
        //            },
        //            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
        //                for(int i = 0; i < Length; i++) {
        //                    dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
        //                }
        //            }
        //        ).Execute();
        //    });

        public static readonly DerivedVariable<XDGField> rho_p0 = new DerivedVariable<XDGField>(
            "rho_p0",
            VariableTypes.Other,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();

                XDGField density = program.Density;
                XDGField tmp = new XDGField(density.Basis);
                CellMask cellMask = CellMask.GetFullMask(program.GridData);

                // Copy only the coordinates that belong to P0
                foreach(string s in program.Control.SpeciesToEvaluate) {
                    foreach(int cell in cellMask.ItemEnum) {
                        foreach(int coordinate in density.GetSpeciesShadowField(s).Basis.GetPolynomialIndicesForDegree(cell, degree: 0)) {
                            tmp.GetSpeciesShadowField(s).Coordinates[cell, coordinate] = density.GetSpeciesShadowField(s).Coordinates[cell, coordinate];
                        }
                    }
                }

                // Save mean value for output
                foreach(string s in program.Control.SpeciesToEvaluate) {
                    foreach(int cell in cellMask.ItemEnum) {
                        double mean = tmp.GetSpeciesShadowField(s).GetMeanValue(cell);

                        dgField.GetSpeciesShadowField(s).SetMeanValue(cell, mean);
                    }
                }
            });

        public static readonly DerivedVariable<XDGField> rho_p0_indicator = new DerivedVariable<XDGField>(
            "rho_p0_indicator",
            VariableTypes.Other,
            delegate (XDGField dgField, XESFMain program) {
                dgField.Clear();

                XDGField density = program.Density;
                XDGField tmp = new XDGField(density.Basis);

                CellMask cellMask = program.LsTrk.Regions.GetCutCellMask();

                // Copy only the coordinates that belong to P0
                foreach(string s in program.Control.SpeciesToEvaluate) {
                    foreach(int cell in cellMask.ItemEnum) {
                        foreach(int coordinate in density.GetSpeciesShadowField(s).Basis.GetPolynomialIndicesForDegree(cell, degree: 0)) {
                            tmp.GetSpeciesShadowField(s).Coordinates[cell, coordinate] = density.GetSpeciesShadowField(s).Coordinates[cell, coordinate];
                        }
                    }
                    //density.GetSpeciesShadowField(s).Coordinates.SaveToTextFile("Density_" + s + ".txt");
                    //dgField.GetSpeciesShadowField(s).Coordinates.SaveToTextFile("dgField_" + s + ".txt");
                }

                // Density on the right side of the shock
                double Ms = 1.5;
                double gamma = 1.4;
                double densityLeft = 1;
                double densityRight = (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms) * densityLeft;
                double meanShock = (densityRight + densityLeft) / 2;

                // Copy only the coordinates that belong to P0
                foreach(string s in program.Control.SpeciesToEvaluate) {
                    foreach(int cell in cellMask.ItemEnum) {
                        double mean_a = tmp.GetSpeciesShadowField(program.Control.SpeciesToEvaluate[0]).GetMeanValue(cell);
                        double mean_b = tmp.GetSpeciesShadowField(program.Control.SpeciesToEvaluate[1]).GetMeanValue(cell);

                        double val = meanShock - (mean_b + mean_a) / 2;
                        dgField.GetSpeciesShadowField(s).SetMeanValue(cell, val);
                    }
                }
            });
        public static readonly DerivedVariable<double> enrichedResidual = new DerivedVariable<double>("enrichedResidual",
            VariableTypes.Other,
            delegate (double enriched_residual, XESFMain program) {
                //CellMask cellMask = program.LsTrk.Regions.GetCutCellMask();
                LevelSetTracker LsTrk = program.LsTrk;
                XDGField[] enriched_fields = new XDGField[program.ConservativeFields.Length];
                XDGBasis[] enriched_basis = new XDGBasis[program.ConservativeFields.Length];
                XDGBasis tmp_basis;
                XDGField tmp_field;

                for(int i = 0; i < program.ConservativeFields.Length; i++) {
                    tmp_basis = new XDGBasis(LsTrk, program.ConservativeFields[i].Basis.Degree + 1);
                    tmp_field = new XDGField(tmp_basis);
                    tmp_field.AccLaidBack(1, program.ConservativeFields[i]);
                    enriched_fields[i] = tmp_field;
                    enriched_basis[i] = tmp_basis;
                }

                var mapping = new CoordinateMapping(enriched_fields);
                //UnsetteledCoordinateMapping mapping = new UnsetteledCoordinateMapping(enriched_basis);
                var eval = program.XSpatialOperator.GetEvaluatorEx(enriched_fields, null, mapping);
                double[] residual_vector = new double[mapping.Ntotal];
                eval.Evaluate(1, 0, residual_vector);


                enriched_residual = residual_vector.L2Norm() / 2;
                //Console.WriteLine("Enriched Residual = " + enriched_residual);



                //// test the original residual
                ////
                //for(int i = 0; i < program.ConservativeFields.Length; i++) {
                //    tmp_basis = new XDGBasis(LsTrk, program.ConservativeFields[i].Basis.Degree);
                //    enriched_basis[i] = tmp_basis;
                //}
                //var mapping2 = new CoordinateMapping(program.ConservativeFields);
                //var eval2 = program.XSpatialOperator.GetEvaluatorEx(program.ConservativeFields, null, mapping2);
                //double[] residual_vector2 = new double[mapping2.Ntotal];
                //eval2.Evaluate(1, 1, residual_vector2);


                //double residual2 = residual_vector2.L2Norm() / 2;

                //Console.WriteLine("original Residual = " + residual2);


            });
        

        public static readonly DerivedVariable<SinglePhaseField> Rankine_Hugnoit_Indicator_LsOne = new DerivedVariable<SinglePhaseField>(
        "RH_LS_One",
        VariableTypes.Other,
        delegate (SinglePhaseField dgField, XESFMain program) {
            var LsTrk = program.LsTrk;
            var left = program.Control.LsOne_PosSpecies;
            var right = program.Control.LsOne_NegSpecies;
            int D = LsTrk.GridDat.SpatialDimension;
            int quadDregree = program.Control.NonlinearQuadratureDegree;
            int levelSetIndex = 0;

            var density_L = program.Density.GetSpeciesShadowField(left);
            var density_R = program.Density.GetSpeciesShadowField(right);
            var momentum_x_L = program.Momentum[0].GetSpeciesShadowField(left);
            var momentum_x_R = program.Momentum[0].GetSpeciesShadowField(right);
            var momentum_y_L = program.Momentum[1].GetSpeciesShadowField(left);
            var momentum_y_R = program.Momentum[1].GetSpeciesShadowField(right);
            var energy_L = program.Energy.GetSpeciesShadowField(left);
            var energy_R = program.Energy.GetSpeciesShadowField(right);

            XDGField pressure = new(program.ConservativeFields[0].Basis);
            XESFVariables.Pressure.UpdateFunction(pressure, program);
            var pressure_L = pressure.GetSpeciesShadowField(left);
            var pressure_R = pressure.GetSpeciesShadowField(right);

            XDGField velocity_x = new(program.ConservativeFields[0].Basis);
            XESFVariables.Velocity[0].UpdateFunction(velocity_x, program);
            var velocity_x_L = velocity_x.GetSpeciesShadowField(left);
            var velocity_x_R = velocity_x.GetSpeciesShadowField(right);

            XDGField velocity_y = new(program.ConservativeFields[0].Basis);
            XESFVariables.Velocity[1].UpdateFunction(velocity_y, program);
            var velocity_y_L = velocity_y.GetSpeciesShadowField(left);
            var velocity_y_R = velocity_y.GetSpeciesShadowField(right);

            void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No of nodes

                MultidimensionalArray density_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray density_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pressure_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pressure_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_x_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_x_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_y_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_y_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray energy_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray energy_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_x_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_x_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_y_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_y_R_Res = MultidimensionalArray.Create(Len, K);

                density_L.Evaluate(j0, Len, NS, density_L_Res);
                density_R.Evaluate(j0, Len, NS, density_R_Res);
                pressure_L.Evaluate(j0, Len, NS, pressure_L_Res);
                pressure_R.Evaluate(j0, Len, NS, pressure_R_Res);
                momentum_x_L.Evaluate(j0, Len, NS, momentum_x_L_Res);
                momentum_x_R.Evaluate(j0, Len, NS, momentum_x_R_Res);
                momentum_y_L.Evaluate(j0, Len, NS, momentum_y_L_Res);
                momentum_y_R.Evaluate(j0, Len, NS, momentum_y_R_Res);
                energy_R.Evaluate(j0, Len, NS, energy_R_Res);
                energy_L.Evaluate(j0, Len, NS, energy_L_Res);
                velocity_x_L.Evaluate(j0, Len, NS, velocity_x_L_Res);
                velocity_x_R.Evaluate(j0, Len, NS, velocity_x_R_Res);
                velocity_y_L.Evaluate(j0, Len, NS, velocity_y_L_Res);
                velocity_y_R.Evaluate(j0, Len, NS, velocity_y_R_Res);

                for(int j = 0; j < Len; j++) {
                    for(int k = 0; k < K; k++) {
                        double a = momentum_x_R_Res[j, k] + momentum_y_R_Res[j, k] - momentum_x_L_Res[j, k] - momentum_y_L_Res[j, k];
                        double b = momentum_x_R_Res[j, k] * (velocity_x_R_Res[j, k] + velocity_y_R_Res[j, k]) + pressure_R_Res[j, k] - momentum_x_L_Res[j, k] * (velocity_x_L_Res[j, k] + velocity_y_L_Res[j, k]) - pressure_L_Res[j, k];
                        double c = momentum_y_R_Res[j, k] * (velocity_x_R_Res[j, k] + velocity_y_R_Res[j, k]) + pressure_R_Res[j, k] - momentum_y_L_Res[j, k] * (velocity_x_L_Res[j, k] + velocity_y_L_Res[j, k]) - pressure_L_Res[j, k];
                        double d = (energy_R_Res[j, k] + pressure_R_Res[j, k]) * (velocity_x_R_Res[j, k] + velocity_y_R_Res[j, k]) - (energy_L_Res[j, k] + pressure_L_Res[j, k]) * (velocity_x_L_Res[j, k] + velocity_y_L_Res[j, k]);
                        result[j, k] = Math.Sqrt(a * a + b * b + c * c + d * d);
                    }
                }
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, quadDregree),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++) {
                        dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                    }
                }
            ).Execute();

            //Console.WriteLine("JumpNormAtShock: " + Math.Sqrt(Jump_NORM));


        });
        public static readonly DerivedVariable<SinglePhaseField> Rankine_Hugnoit_Indicator_LsTwo = new DerivedVariable<SinglePhaseField>(
        "RH_LS_Two",
        VariableTypes.Other,
        delegate (SinglePhaseField dgField, XESFMain program) {
            var LsTrk = program.LsTrk;
            var left = program.Control.LsTwo_PosSpecies;
            var right = program.Control.LsTwo_NegSpecies;
            int D = LsTrk.GridDat.SpatialDimension;
            int quadDregree = program.Control.NonlinearQuadratureDegree;
            int levelSetIndex = 1;

            var density_L = program.Density.GetSpeciesShadowField(left);
            var density_R = program.Density.GetSpeciesShadowField(right);
            var momentum_x_L = program.Momentum[0].GetSpeciesShadowField(left);
            var momentum_x_R = program.Momentum[0].GetSpeciesShadowField(right);
            var momentum_y_L = program.Momentum[1].GetSpeciesShadowField(left);
            var momentum_y_R = program.Momentum[1].GetSpeciesShadowField(right);
            var energy_L = program.Energy.GetSpeciesShadowField(left);
            var energy_R = program.Energy.GetSpeciesShadowField(right);

            XDGField pressure = new(program.ConservativeFields[0].Basis);
            XESFVariables.Pressure.UpdateFunction(pressure, program);
            var pressure_L = pressure.GetSpeciesShadowField(left);
            var pressure_R = pressure.GetSpeciesShadowField(right);

            XDGField velocity_x = new(program.ConservativeFields[0].Basis);
            XESFVariables.Velocity[0].UpdateFunction(velocity_x, program);
            var velocity_x_L = velocity_x.GetSpeciesShadowField(left);
            var velocity_x_R = velocity_x.GetSpeciesShadowField(right);

            XDGField velocity_y = new(program.ConservativeFields[0].Basis);
            XESFVariables.Velocity[1].UpdateFunction(velocity_y, program);
            var velocity_y_L = velocity_y.GetSpeciesShadowField(left);
            var velocity_y_R = velocity_y.GetSpeciesShadowField(right);

            void RH_Condition(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No of nodes

                MultidimensionalArray density_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray density_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pressure_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pressure_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_x_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_x_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_y_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray momentum_y_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray energy_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray energy_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_x_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_x_R_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_y_L_Res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray velocity_y_R_Res = MultidimensionalArray.Create(Len, K);

                density_L.Evaluate(j0, Len, NS, density_L_Res);
                density_R.Evaluate(j0, Len, NS, density_R_Res);
                pressure_L.Evaluate(j0, Len, NS, pressure_L_Res);
                pressure_R.Evaluate(j0, Len, NS, pressure_R_Res);
                momentum_x_L.Evaluate(j0, Len, NS, momentum_x_L_Res);
                momentum_x_R.Evaluate(j0, Len, NS, momentum_x_R_Res);
                momentum_y_L.Evaluate(j0, Len, NS, momentum_y_L_Res);
                momentum_y_R.Evaluate(j0, Len, NS, momentum_y_R_Res);
                energy_R.Evaluate(j0, Len, NS, energy_R_Res);
                energy_L.Evaluate(j0, Len, NS, energy_L_Res);
                velocity_x_L.Evaluate(j0, Len, NS, velocity_x_L_Res);
                velocity_x_R.Evaluate(j0, Len, NS, velocity_x_R_Res);
                velocity_y_L.Evaluate(j0, Len, NS, velocity_y_L_Res);
                velocity_y_R.Evaluate(j0, Len, NS, velocity_y_R_Res);

                for(int j = 0; j < Len; j++) {
                    for(int k = 0; k < K; k++) {
                        double a = momentum_x_R_Res[j, k] + momentum_y_R_Res[j, k] - momentum_x_L_Res[j, k] - momentum_y_L_Res[j, k];
                        double b = momentum_x_R_Res[j, k] * (velocity_x_R_Res[j, k] + velocity_y_R_Res[j, k]) + pressure_R_Res[j, k] - momentum_x_L_Res[j, k] * (velocity_x_L_Res[j, k] + velocity_y_L_Res[j, k]) - pressure_L_Res[j, k];
                        double c = momentum_y_R_Res[j, k] * (velocity_x_R_Res[j, k] + velocity_y_R_Res[j, k]) + pressure_R_Res[j, k] - momentum_y_L_Res[j, k] * (velocity_x_L_Res[j, k] + velocity_y_L_Res[j, k]) - pressure_L_Res[j, k];
                        double d = (energy_R_Res[j, k] + pressure_R_Res[j, k]) * (velocity_x_R_Res[j, k] + velocity_y_R_Res[j, k]) - (energy_L_Res[j, k] + pressure_L_Res[j, k]) * (velocity_x_L_Res[j, k] + velocity_y_L_Res[j, k]);
                        result[j, k] = Math.Sqrt(a * a + b * b + c * c + d * d);
                    }
                }
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId(left), LsTrk.GetSpeciesId(right) }, quadDregree, HistoryIndex: 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet: levelSetIndex, IntegrationDom: LsTrk.Regions.GetCutCellMask4LevSet(LevSetIdx: levelSetIndex));

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, quadDregree),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    RH_Condition(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++) {
                        dgField.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0] * 10);
                    }
                }
            ).Execute();

            //Console.WriteLine("JumpNormAtShock: " + Math.Sqrt(Jump_NORM));


        });


    }
}