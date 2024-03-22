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
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.CompressibleFlowCommon;
using ilPSP;
using MathNet.Numerics;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using XDGShock.TimeStepping;
using static BoSSS.Solution.CompressibleFlowCommon.Variable;

namespace XESTSF.Variables {
    public static class XESTSFVariables {

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
            delegate (XDGField dgField, XESTSFMain program) {
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
            delegate (XDGField dgField, XESTSFMain program) {
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
            delegate (double residual, XESTSFMain program) {

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
        /// The optional pressure field:
        /// \f[ p = f(\rho, \vec{m}, \rho E) \f], where
        /// \f[ f \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable<XDGField> Pressure = new DerivedVariable<XDGField>(
            "p",
            VariableTypes.Pressure,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);

            }
            );

        /// <summary>
        /// The optional Pertubation pressure field
        /// </summary>
        public static readonly DerivedVariable<XDGField> PertubationPressure = new DerivedVariable<XDGField>(
            "p_per",
            VariableTypes.Pressure,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                // get pressure
                dgField.ProjectFunctionXDG(1.0,
                        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);
                //compute the normalized pertubation
                foreach(string spc in program.Control.SpeciesToEvaluate) {
                    var sf = dgField.GetSpeciesShadowField(spc);
                    var tp = new Tuple<string, string>("p", spc);
                    program.Control.BaseFlowPerSpeciesAndField.TryGetValue(tp, out double baseFlow);
                    sf.AccConstant(-baseFlow);
                    
                }
            }
            );
        /// <summary>
        /// The optional Pertubation density field
        /// </summary>
        public static readonly DerivedVariable<XDGField> PertubationDensity = new DerivedVariable<XDGField>(
            "rho_per",
            VariableTypes.Density,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                // get pressure
                dgField.AccLaidBack(1.0, program.ConservativeFields[0]);
                //compute the normalized pertubation
                foreach(string spc in program.Control.SpeciesToEvaluate) {
                    var sf = dgField.GetSpeciesShadowField(spc);
                    var tp = new Tuple<string, string>("rho", spc);
                    program.Control.BaseFlowPerSpeciesAndField.TryGetValue(tp, out double baseFlow);
                    sf.AccConstant(-baseFlow);
                }
            }
            );
        /// <summary>
        /// The optional Pertubation x-velocity field
        /// </summary>
        public static readonly DerivedVariable<XDGField> PertubationVelocityX = new DerivedVariable<XDGField>(
            "u0_per",
            VariableTypes.Velocity,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                // get pressure
                dgField.ProjectFunctionXDG(1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Velocity[0],
                    program.SpeciesToEvaluate_Ids,
                    program.ConservativeFields);
                //compute the normalized pertubation
                foreach(string spc in program.Control.SpeciesToEvaluate) {
                    var sf = dgField.GetSpeciesShadowField(spc);
                    var tp = new Tuple<string, string>("u0", spc);
                    program.Control.BaseFlowPerSpeciesAndField.TryGetValue(tp, out double baseFlow);
                    sf.AccConstant(-baseFlow);
                }
            }
            );
        /// <summary>
        /// The optional Pertubation y-velocity field
        /// </summary>
        public static readonly DerivedVariable<XDGField> PertubationVelocityY = new DerivedVariable<XDGField>(
            "u1_per",
            VariableTypes.Velocity,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                // get pressure
                dgField.ProjectFunctionXDG(1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Velocity[1],
                    program.SpeciesToEvaluate_Ids,
                    program.ConservativeFields);
                //compute the normalized pertubation
                foreach(string spc in program.Control.SpeciesToEvaluate) {
                    var sf = dgField.GetSpeciesShadowField(spc);
                    var tp = new Tuple<string, string>("u1", spc);
                    program.Control.BaseFlowPerSpeciesAndField.TryGetValue(tp, out double baseFlow);
                    sf.AccConstant(-baseFlow);
                }
            }
            );
        /// <summary>
        /// The optional Pertubation x momentum field
        /// </summary>
        public static readonly DerivedVariable<XDGField> PertubationMomentumX = new DerivedVariable<XDGField>(
            "m0_per",
            VariableTypes.Momentum,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                // get pressure
                dgField.AccLaidBack(1.0, program.ConservativeFields[1]);
                //compute the normalized pertubation
                foreach(string spc in program.Control.SpeciesToEvaluate) {
                    var sf = dgField.GetSpeciesShadowField(spc);
                    var tp = new Tuple<string, string>("m0", spc);
                    program.Control.BaseFlowPerSpeciesAndField.TryGetValue(tp, out double baseFlow);
                    sf.AccConstant(-baseFlow);
                }
            }
            );
        /// <summary>
        /// The optional Pertubation x momentum field
        /// </summary>
        public static readonly DerivedVariable<XDGField> PertubationMomentumY = new DerivedVariable<XDGField>(
            "m1_per",
            VariableTypes.Momentum,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                // get pressure
                dgField.AccLaidBack(1.0, program.ConservativeFields[2]);
                //compute the normalized pertubation
                foreach(string spc in program.Control.SpeciesToEvaluate) {
                    var sf = dgField.GetSpeciesShadowField(spc);
                    var tp = new Tuple<string, string>("m1", spc);
                    program.Control.BaseFlowPerSpeciesAndField.TryGetValue(tp, out double baseFlow);
                    sf.AccConstant(-baseFlow);
                }
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
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Pressure,
                    program.SpeciesToEvaluate_Ids,
                    program.originalStep);
               
            }
         );

        /// <summary>
        /// The optional Pertubation density field
        /// </summary>
        public static readonly DerivedVariable<XDGField> PertubationEnergy = new DerivedVariable<XDGField>(
            "rhoE_per",
            VariableTypes.Energy,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                // get pressure
                dgField.AccLaidBack(1.0, program.ConservativeFields[program.ConservativeFields.Length-1]);
                //compute the normalized pertubation
                foreach(string spc in program.Control.SpeciesToEvaluate) {
                    var sf = dgField.GetSpeciesShadowField(spc);
                    var tp = new Tuple<string, string>("rhoE", spc);
                    program.Control.BaseFlowPerSpeciesAndField.TryGetValue(tp, out double baseFlow);
                    sf.AccConstant(-baseFlow);
                }
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
                delegate (XDGField dgField, XESTSFMain program) {
                    dgField.Clear();
                    dgField.ProjectFunctionXDG(1.0,
        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Velocity[d],
        program.SpeciesToEvaluate_Ids,
        program.ConservativeFields);
                  
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
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                        (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Enthalpy,
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);
               
            });
        /// <summary>
        /// The optional local Mach number field
        /// \f[ Ma = \frac{\Vert u \Vert}{a} \f], where the formula for the
        /// local speed of sound \f[ a \f] depends on the equation of state
        /// </summary>
        public static readonly DerivedVariable<XDGField> LocalMachNumber = new DerivedVariable<XDGField>(
            "Ma",
            VariableTypes.Other,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(
                        1.0,
                        delegate (Vector X, double[] U, int j) {
                            var state = new StateVector(U, program.Control.GetMaterial());
                            double r = state.LocalMachNumber;
                            //Debug.Assert(r.IsNaN() == false, "NAN in local mach number");
                            //Debug.Assert(r.IsInfinity() == false, "INF in local mach number");
                            return r;
                        },
                        program.SpeciesToEvaluate_Ids,
                        program.ConservativeFields);
            });

        /// <summary>
        /// The optional entropy field (see
        /// <see cref="BoSSS.Solution.CompressibleFlowCommon.MaterialProperty.IEquationOfState.GetEntropy"/>).
        /// </summary>
        public static readonly DerivedVariable<XDGField> Entropy = new DerivedVariable<XDGField>(
            "S",
            VariableTypes.Other,
            delegate (XDGField dgField, XESTSFMain program) {
                dgField.Clear();
                dgField.ProjectFunctionXDG(1.0,
                    (X, U, j) => new StateVector(U, program.Control.GetMaterial()).Entropy,
                    program.SpeciesToEvaluate_Ids,
                    program.ConservativeFields);
            });

    }
}