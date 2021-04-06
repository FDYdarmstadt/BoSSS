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

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Numerics;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.Tracing;
using ilPSP.LinSolvers;

using BoSSS.Platform;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.XDG;

using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.Grid.Aggregation;
using NUnit.Framework;
using MPI.Wrappers;
using System.Collections;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XNSE_Solver.Legacy {

    /// <summary>
    /// Solver for Incompressible Multiphase flows; 
    /// </summary>
    public partial class XNSE_SolverMain : BoSSS.Solution.Application<XNSE_Control> {

        //=============================================
        // partial file for heat equation related code
        //=============================================


        //=====================================
        // Field declaration and instantiation
        //=====================================
        #region fields

#pragma warning disable 649

        /// <summary>
        /// Temperature
        /// </summary>
        XDGField Temperature;

        /// <summary>
        /// Heat flux; either derived via temperature gradient
        /// or directly solved via LDG
        /// </summary>
        VectorField<XDGField> HeatFlux;

        /// <summary>
        /// Residual of the heat equation
        /// </summary>
        XDGField ResidualHeat;

        /// <summary>
        /// Residual of the auxiliary heat flux equations for LDG
        /// </summary>
        VectorField<XDGField> ResidualAuxHeatFlux;

        /// <summary>
        /// 
        /// </summary>
        SinglePhaseField MassFluxField;

        /// <summary>
        /// 
        /// </summary>
        SinglePhaseField MassFluxExtension;

        /// <summary>
        /// prescribed disjoining pressure field for evaporation near wall 
        /// </summary>
        SinglePhaseField DisjoiningPressure;

#pragma warning restore 649

        /// <summary>
        /// creates heat equation related fields
        /// </summary>
        public void CreateHeatFields() {

            int D = this.GridData.SpatialDimension;

            this.Temperature = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Temperature].Degree), VariableNames.Temperature);
            base.RegisterField(this.Temperature);
            this.ResidualHeat = new XDGField(this.Temperature.Basis, "ResidualHeat");
            base.RegisterField(this.ResidualHeat);

            this.HeatFlux = new VectorField<XDGField>(D.ForLoop(d => new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Temperature].Degree), VariableNames.HeatFluxVectorComponent(d))));
            base.RegisterField(this.HeatFlux);
            if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                this.ResidualAuxHeatFlux = new VectorField<XDGField>(D.ForLoop(d => new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.HeatFluxVectorComponent(d)].Degree), VariableNames.ResidualAuxHeatFluxVectorComponent(d))));
                base.RegisterField(this.ResidualAuxHeatFlux);
            }

            this.MassFluxField = new SinglePhaseField(new Basis(this.GridData, this.Control.FieldOptions[VariableNames.VelocityX].Degree), "MassFluxField");
            base.RegisterField(MassFluxField);
            this.MassFluxExtension = new SinglePhaseField(new Basis(this.GridData, this.Control.FieldOptions[VariableNames.VelocityX].Degree), VariableNames.MassFluxExtension);
            base.RegisterField(MassFluxExtension);

            this.DisjoiningPressure = new SinglePhaseField(new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree), "DisjoiningPressure");
            if (this.Control.DisjoiningPressureFunc != null) {
                DisjoiningPressure.ProjectField(this.Control.DisjoiningPressureFunc);
            }
            base.RegisterField(this.DisjoiningPressure);

        }


        #endregion



        ThermalMultiphaseBoundaryCondMap m_thermBcMap;

        /// <summary>
        /// Boundary conditions.
        /// </summary>
        ThermalMultiphaseBoundaryCondMap thermBcMap {
            get {
                if (m_thermBcMap == null && this.Control.solveCoupledHeatEquation) {
                    m_thermBcMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
                }
                return m_thermBcMap;
            }
        }
   


        public void ComputeHeatFlux() {
            using (FuncTrace ft = new FuncTrace()) {

                int D = this.LsTrk.GridDat.SpatialDimension;
                for (int d = 0; d < D; d++) {

                    this.HeatFlux[d].Clear();

                    foreach (var Spc in LsTrk.SpeciesNames) { // loop over species...
                                                              // shadow fields
                        DGField Temp_Spc = (this.Temperature.GetSpeciesShadowField(Spc));

                        double kSpc = 0.0;
                        switch (Spc) {
                            case "A": kSpc = -this.Control.ThermalParameters.k_A; break;
                            case "B": kSpc = -this.Control.ThermalParameters.k_B; break;
                            default: throw new NotSupportedException("Unknown species name '" + Spc + "'");
                        }

                        this.HeatFlux[d].GetSpeciesShadowField(Spc).DerivativeByFlux(kSpc, Temp_Spc, d, this.LsTrk.Regions.GetSpeciesSubGrid(Spc));
                    }
                }

                this.HeatFlux.ForEach(F => F.CheckForNanOrInf(true, true, true));

            }
        }



        public void ComputeMassFluxField() {

            double kA = this.Control.ThermalParameters.k_A;
            double kB = this.Control.ThermalParameters.k_B;

            int order = MassFluxExtension.Basis.Degree * MassFluxExtension.Basis.Degree + 2;

            MassFluxField.Clear();

            MassFluxField.ProjectField(1.0,
                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    int D = this.LsTrk.GridDat.SpatialDimension;

                    MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                    this.Temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                    this.Temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

                    MultidimensionalArray HeatFluxA_Res = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray HeatFluxB_Res = MultidimensionalArray.Create(Len, K, D);
                    if (XOpConfig.getConductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                        for (int dd = 0; dd < D; dd++) {
                            this.HeatFlux[dd].GetSpeciesShadowField("A").Evaluate(j0, Len, NS, HeatFluxA_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                            this.HeatFlux[dd].GetSpeciesShadowField("B").Evaluate(j0, Len, NS, HeatFluxB_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                        }
                    }

                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                    for (int j = 0; j < Len; j++) {

                        MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                        this.GridData.TransformLocal2Global(NS, globCoord, j);

                        for (int k = 0; k < K; k++) {

                            double qEvap = 0.0;
                            //macro region
                            for (int dd = 0; dd < D; dd++) {
                                if (XOpConfig.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                                    qEvap += ((-kB) * GradTempB_Res[j, k, dd] - (-kA) * GradTempA_Res[j, k, dd]) * Normals[j, k, dd];
                                } else {
                                    qEvap += (HeatFluxB_Res[j, k, dd] - HeatFluxA_Res[j, k, dd]) * Normals[j, k, dd];
                                }
                            }

                            //Console.WriteLine("qEvap delUpdateLevelSet = {0}", qEvap);
                            double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };
                            double mEvap = (this.XOpConfig.prescribedMassflux != null) ? this.XOpConfig.prescribedMassflux(globX, hack_Phystime) 
                                                                                            : qEvap / this.Control.ThermalParameters.hVap; // mass flux

                            //Console.WriteLine("mEvap - delUpdateLevelSet = {0}", mEvap);
                            result[j, k] = mEvap;  
                        }
                    }
                }, (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order));

        }


        public void ConstructMassFluxExtensionField() {

            SubGrid CCgrid = LsTrk.Regions.GetCutCellSubGrid();
            CellMask CC = LsTrk.Regions.GetCutCellMask();
            CellMask NEAR = LsTrk.Regions.GetNearFieldMask(1);
            int J = this.LsTrk.GridDat.Cells.NoOfLocalUpdatedCells;
            double[][] MassFluxMin = new double[1][];
            double[][] MassFluxMax = new double[1][];
            MassFluxMin[0] = new double[J];
            MassFluxMax[0] = new double[J];

            NarrowMarchingBand.ConstructExtVel_PDE(this.LsTrk, CCgrid, new SinglePhaseField[] { MassFluxExtension }, new SinglePhaseField[] { MassFluxField },
                this.DGLevSet.Current, this.DGLevSetGradient, MassFluxMin, MassFluxMax, this.m_HMForder);

            var marcher = new FastMarchReinit(this.DGLevSet.Current.Basis);
            marcher.ConstructExtension(this.DGLevSet.Current, NEAR.Except(CC), CC, new SinglePhaseField[] { MassFluxExtension },
                MassFluxMin, MassFluxMax, this.DGLevSetGradient, this.hack_TimestepIndex);

            MassFluxExtension.CheckForNanOrInf(true, true, true);
            
        }

    }

}
