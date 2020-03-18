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

namespace BoSSS.Application.XNSE_Solver {

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



    }

}
