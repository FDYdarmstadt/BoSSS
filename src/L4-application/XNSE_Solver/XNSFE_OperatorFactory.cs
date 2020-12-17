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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;

using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.EnergyCommon;
using BoSSS.Solution.XheatCommon;
using System.Collections;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// class for defining the equation components of the XNSFE_Operator (extended Navier-Stokes-Fourier equations) 
    /// and assembly of the corresponding matrix 
    /// </summary>
    public class XNSFE_OperatorFactory : XOperatorFactoryBase {


        string[] CodNameSelected = new string[0];
        string[] DomNameSelected = new string[0];

        string[] Params;

        int HMFDegree;
        XNSFE_OperatorConfiguration config;

        /// <summary>
        /// ctor for the operator factory, where the equation compnents are set
        /// </summary>
        /// <param name="_config"></param>
        /// <param name="_LsTrk"></param>
        /// <param name="_HMFdegree"></param>
        /// <param name="BcMap"></param>
        /// <param name="thermBcMap"></param>
        /// <param name="degU"></param>
        public XNSFE_OperatorFactory(XNSFE_OperatorConfiguration _config, LevelSetTracker _LsTrk, int _HMFdegree, 
            IncompressibleMultiphaseBoundaryCondMap BcMap, ThermalMultiphaseBoundaryCondMap thermBcMap, int degU, IDictionary<SpeciesId, IEnumerable<double>> MassScale) {

            this.config = _config;
            this.LsTrk = _LsTrk;
            this.D = _LsTrk.GridDat.SpatialDimension;

            this.HMFDegree = _HMFdegree;

            this.physParams = config.getPhysParams;
            this.thermParams = config.getThermParams;
            this.dntParams = config.getDntParams;


            // test input
            // ==========
            {
                if ((config.getPhysParams.mu_A <= 0) && (config.getPhysParams.mu_B <= 0)) {
                    config.isViscous = false;
                } else {
                    if ((config.getPhysParams.mu_A <= 0) || (config.getPhysParams.mu_B <= 0))
                        throw new ArgumentException("Viscosity does not fall within expected range");
                }

                if ((config.getPhysParams.rho_A <= 0) || (config.getPhysParams.rho_B <= 0))
                    throw new ArgumentException("Density does not fall within expected range");

                if (_LsTrk.SpeciesNames.Count != 2)
                    throw new ArgumentException("Species name count does not fall within expected range");

                if (!(_LsTrk.SpeciesNames.Contains("A") && _LsTrk.SpeciesNames.Contains("B")))
                    throw new ArgumentException("Missing species name");
            }


            // full operator:
            // ==============
            CodName = ArrayTools.Cat(EquationNames.MomentumEquations(D), EquationNames.ContinuityEquation);
            Params = ArrayTools.Cat(
                VariableNames.Velocity0Vector(D),
                VariableNames.Velocity0MeanVector(D),
                VariableNames.NormalVector(D),
                VariableNames.Curvature,
                VariableNames.SurfaceForceVector(D)
                );
            DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure);

            //if (config.solveEnergy) {
            //    CodName = ArrayTools.Cat(CodName, EquationNames.KineticEnergyEquation);
            //    Params = ArrayTools.Cat(Params,
            //        VariableNames.VelocityX_GradientVector(),
            //        VariableNames.VelocityY_GradientVector(),
            //        new string[] { "VelocityXGradX_GradientX", "VelocityXGradX_GradientY" },
            //        new string[] { "VelocityXGradY_GradientX", "VelocityXGradY_GradientY" },
            //        new string[] { "VelocityYGradX_GradientX", "VelocityYGradX_GradientY" },
            //        new string[] { "VelocityYGradY_GradientX", "VelocityYGradY_GradientY" },
            //        VariableNames.Pressure0,
            //        VariableNames.PressureGradient(D),
            //        VariableNames.GravityVector(D)
            //        );
            //    DomName = ArrayTools.Cat(DomName, VariableNames.KineticEnergy);
            //}

            if (config.solveHeat) {
                Params = ArrayTools.Cat(Params,
                    VariableNames.Temperature0,
                    VariableNames.HeatFlux0Vector(D),
                    VariableNames.MassFluxExtension
                    );
                if (config.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                    CodName = ArrayTools.Cat(CodName, EquationNames.HeatEquation);
                    DomName = ArrayTools.Cat(DomName, VariableNames.Temperature);
                } else {
                    CodName = ArrayTools.Cat(CodName, EquationNames.HeatEquation, EquationNames.AuxHeatFlux(D));
                    DomName = ArrayTools.Cat(DomName, VariableNames.Temperature, VariableNames.HeatFluxVector(D));
                }
            }

            storedParams = new DGField[Params.Length];

            // selected part:
            if (config.getCodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, D));
            if (config.getCodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D, 1));

            if (config.getDomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, D));
            if (config.getDomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D, 1));

            int nBlocks = 2;
            //if (config.solveEnergy) {
            //    nBlocks = 3;
            //    if (config.getCodBlocks[2])
            //        CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D + 1, 1));
            //    if (config.getDomBlocks[2])
            //        DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D + 1, 1));
            //}

            if (config.solveHeat) {
                if (config.getCodBlocks[nBlocks])
                    CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D + 1, 1));
                if (config.getDomBlocks[nBlocks])
                    DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D + 1, 1));

                if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                    if (config.getCodBlocks[nBlocks + 1])
                        CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D + 2, D));
                    if (config.getDomBlocks[nBlocks + 1])
                        DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D + 2, D));
                }
            }


            // create Operator
            // ===============

            m_XOp = new XSpatialOperatorMk2(DomNameSelected, Params, CodNameSelected, (A, B, C) => _HMFdegree, this.LsTrk.SpeciesNames);

            // add Navier-Stokes components
            // ============================

            // species bulk components
            for (int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
                // Navier Stokes equations
                Solution.XNSECommon.XOperatorComponentsFactory.AddSpeciesNSE(m_XOp, config, D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, LsTrk, out U0meanrequired);

                // continuity equation
                if (config.isContinuity)
                    Solution.XNSECommon.XOperatorComponentsFactory.AddSpeciesContinuityEq(m_XOp, config, D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap);
            }

            // interface components
            Solution.XNSECommon.XOperatorComponentsFactory.AddInterfaceNSE(m_XOp, config, D, BcMap, LsTrk);     // surface stress tensor
            Solution.XNSECommon.XOperatorComponentsFactory.AddSurfaceTensionForce(m_XOp, config, D, BcMap, LsTrk, degU, out NormalsRequired, out CurvatureRequired);     // surface tension force

            if (config.isContinuity)
                Solution.XNSECommon.XOperatorComponentsFactory.AddInterfaceContinuityEq(m_XOp, config, D, LsTrk);       // continuity equation


            // add kinetic energy equation components
            // ======================================
            //if (config.solveEnergy) {

            //    // species bulk components
            //    for (int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
            //        Solution.EnergyCommon.XOperatorComponentsFactory.AddSpeciesKineticEnergyEquation(m_XOp, config, D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, LsTrk);
            //    }

            //    // interface components
            //    Solution.EnergyCommon.XOperatorComponentsFactory.AddInterfaceKineticEnergyEquation(m_XOp, config, D, BcMap, LsTrk, degU);
            //    CurvatureRequired = true;
            //}


            // add Heat equation components
            // ============================
            if (config.solveHeat) {

                // species bulk components
                for (int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
                    Solution.XheatCommon.XOperatorComponentsFactory.AddSpeciesHeatEq(m_XOp, config,
                         D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], thermBcMap, LsTrk);
                }

                // interface components
                Solution.XheatCommon.XOperatorComponentsFactory.AddInterfaceHeatEq(m_XOp, config, D, thermBcMap, LsTrk);
            }


            // add Evaporation interface components
            // ====================================

            if (config.isEvaporation) {

                XOperatorComponentsFactory.AddInterfaceNSE_withEvaporation(m_XOp, config, D, LsTrk);
                if (config.isContinuity)
                    XOperatorComponentsFactory.AddInterfaceContinuityEq_withEvaporation(m_XOp, config, D, LsTrk);

            }

            // create temporal operator
            // ========================
            var TempOp = new ConstantXTemporalOperator(m_XOp, 0.0);
            m_XOp.TemporalOperator = TempOp;
            foreach(var kv in MassScale) {
                TempOp.DiagonalScale[LsTrk.GetSpeciesName(kv.Key)].SetV(kv.Value.ToArray());
            }

            // Finalize
            // ========

            m_XOp.Commit();
        }

        PhysicalParameters physParams;
        ThermalParameters thermParams;
        DoNotTouchParameters dntParams;

        bool NormalsRequired;

        bool CurvatureRequired;

        bool U0meanrequired;

        /// <summary>
        /// 
        /// </summary>
        public void AssembleMatrix<T>(BlockMsrMatrix OpMatrix, double[] OpAffine,
            UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping,
            IEnumerable<T> CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time, double dt,
            int CutCellQuadOrder, VectorField<SinglePhaseField> SurfaceForce,
            VectorField<SinglePhaseField> LevelSetGradient, SinglePhaseField ExternalyProvidedCurvature,
            bool[] updateSolutionParams = null, DGField[] ExtParams = null) where T : DGField {
            //IEnumerable<T> CoupledCurrentState = null, IEnumerable<T> CoupledParams = null) where T : DGField {

            // checks:
            if (ColMapping.BasisS.Count != this.m_XOp.DomainVar.Count)
                throw new ArgumentException();
            if (RowMapping.BasisS.Count != this.m_XOp.CodomainVar.Count)
                throw new ArgumentException();

            int D = this.LsTrk.GridDat.SpatialDimension;
            if (CurrentState != null && !config.solveEnergy && !config.solveHeat && CurrentState.Count() != (D + 1))
                throw new ArgumentException();

            if (OpMatrix == null && CurrentState == null)
                throw new ArgumentException();

            DGField[] U0;
            if (CurrentState != null)
                U0 = CurrentState.Take(D).ToArray();
            else
                U0 = null;



            // parameter assembly
            // ==================
            #region param assembly

            LevelSet Phi = (LevelSet)(this.LsTrk.LevelSets[0]);
            SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();


            // linearization velocity:
            DGField[] U0_U0mean;
            if (this.U0meanrequired || this.config.isEvaporation) {
                XDGBasis U0meanBasis = new XDGBasis(this.LsTrk, 0);
                VectorField<XDGField> U0mean = new VectorField<XDGField>(D, U0meanBasis, "U0mean_", XDGField.Factory);
                U0mean.Clear();
                if (this.physParams.IncludeConvection)
                    ComputeAverageU(U0, U0mean, CutCellQuadOrder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, CutCellQuadOrder, 1).XQuadSchemeHelper);

                U0_U0mean = ArrayTools.Cat<DGField>(U0, U0mean);
            } else {
                U0_U0mean = new DGField[2 * D];
            }

            // linearization velocity:
            //if (this.U0meanrequired) {
            //    VectorField<XDGField> U0mean = new VectorField<XDGField>(U0_U0mean.Skip(D).Take(D).Select(f => ((XDGField)f)).ToArray());

            //    U0mean.Clear();
            //    if (this.physParams.IncludeConvection)
            //        ComputeAverageU(U0, U0mean, CutCellQuadOrder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, CutCellQuadOrder, 1).XQuadSchemeHelper);
            //}


            // normals:
            SinglePhaseField[] Normals; // Normal vectors: length not normalized - will be normalized at each quad node within the flux functions.
            if (this.NormalsRequired) {
                if (LevelSetGradient == null) {
                    LevelSetGradient = new VectorField<SinglePhaseField>(D, Phi.Basis, SinglePhaseField.Factory);
                    LevelSetGradient.Gradient(1.0, Phi);
                }
                Normals = LevelSetGradient.ToArray();
            } else {
                Normals = new SinglePhaseField[D];
            }

            // curvature:
            SinglePhaseField Curvature;
            if (this.CurvatureRequired) {
                Curvature = ExternalyProvidedCurvature;
            } else {
                Curvature = null;
            }


            // velocity gradient vectors
            var VelMap = new CoordinateMapping(U0);
            DGField[] VelParam = VelMap.Fields.ToArray();

            VectorField<DGField> GradVelX = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((VelParam[0] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (GradVelX[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            GradVelX.ForEach(F => F.CheckForNanOrInf(true, true, true));

            VectorField<DGField> GradVelY = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((VelParam[1] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (GradVelY[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            GradVelY.ForEach(F => F.CheckForNanOrInf(true, true, true));


            VectorField<DGField> GradVelXGradX = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradX_Gradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((GradVelX[0] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (GradVelXGradX[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            GradVelXGradX.ForEach(F => F.CheckForNanOrInf(true, true, true));

            VectorField<DGField> GradVelXGradY = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradY_Gradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((GradVelX[1] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (GradVelXGradY[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            GradVelXGradY.ForEach(F => F.CheckForNanOrInf(true, true, true));

            VectorField<DGField> GradVelYGradX = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradX_Gradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((GradVelY[0] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (GradVelYGradX[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            GradVelYGradX.ForEach(F => F.CheckForNanOrInf(true, true, true));

            VectorField<DGField> GradVelYGradY = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradY_Gradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((GradVelY[1] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (GradVelYGradY[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            GradVelYGradY.ForEach(F => F.CheckForNanOrInf(true, true, true));



            // pressure and gradient
            var PressMap = new CoordinateMapping(CurrentState.ToArray()[D]);
            DGField[] PressParam = PressMap.Fields.ToArray();

            VectorField<DGField> PressGrad = new VectorField<DGField>(D, PressParam[0].Basis, "PressureGrad", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((PressParam[0] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (PressGrad[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            PressGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

            // gravity
            //var GravMap = new CoordinateMapping(ExtParams);
            //DGField[] GravParam = GravMap.Fields.ToArray();


            // heat flux for evaporation
            DGField[] HeatFluxParam = new DGField[D];
            if (config.solveHeat) {
                if (config.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP && updateSolutionParams[D+1]) {
                    HeatFluxParam = new VectorField<XDGField>(D, CurrentState.ToArray()[D + 1].Basis, "HeatFlux0_", XDGField.Factory).ToArray();
                    Dictionary<string, double> kSpc = new Dictionary<string, double>();
                    kSpc.Add("A", -thermParams.k_A);
                    kSpc.Add("B", -thermParams.k_B);
                    XNSEUtils.ComputeGradientForParam(CurrentState.ToArray()[D + 1], HeatFluxParam, this.LsTrk, kSpc, this.LsTrk.Regions.GetCutCellSubGrid());
                } else if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP && updateSolutionParams[D+2]) {
                    var HeatFluxMap = new CoordinateMapping(CurrentState.ToArray().GetSubVector(D + 2, D));
                    HeatFluxParam = HeatFluxMap.Fields.ToArray();
                } else {
                    HeatFluxParam = storedParams.GetSubVector(2 * D + 4, D);
                }
            }
            //if (ExtParams != null) {
            //    HeatFluxParam = ExtParams;
            //}

            #endregion


            // concatenate everything
            var Params = ArrayTools.Cat<DGField>(
                U0_U0mean,
                Normals,
                Curvature,
                ((SurfaceForce != null) ? SurfaceForce.ToArray() : new SinglePhaseField[D]));

            //if (config.solveEnergy) {
            //    Params = ArrayTools.Cat<DGField>(Params.ToArray<DGField>(),
            //        GradVelX,
            //        GradVelY,
            //        GradVelXGradX,
            //        GradVelXGradY,
            //        GradVelYGradX,
            //        GradVelYGradY,
            //        PressParam,
            //        PressGrad,
            //        GravMap);
            //}

            if (config.solveHeat) {
                Params = ArrayTools.Cat<DGField>(Params.ToArray<DGField>(),
                    CurrentState.ToArray<DGField>().GetSubVector(D + 1, 1),
                    HeatFluxParam,
                    ExtParams);
            }


            // store old params
            for (int p = 0; p < Params.Length; p++) {
                if (Params[p] != null)
                    storedParams[p] = Params[p].CloneAs();
            }



            // advanced settings for the navier slip boundary condition
            // ========================================================


            CellMask SlipArea;
            switch (this.dntParams.GNBC_Localization) {
                case NavierSlip_Localization.Bulk: {
                        SlipArea = this.LsTrk.GridDat.BoundaryCells.VolumeMask;
                        break;
                    }
                case NavierSlip_Localization.ContactLine: {
                        SlipArea = null;
                        break;
                    }
                case NavierSlip_Localization.Nearband: {
                        SlipArea = this.LsTrk.GridDat.BoundaryCells.VolumeMask.Intersect(this.LsTrk.Regions.GetNearFieldMask(this.LsTrk.NearRegionWidth));
                        break;
                    }
                case NavierSlip_Localization.Prescribed: {
                        throw new NotImplementedException();
                    }
                default:
                    throw new ArgumentException();
            }


            MultidimensionalArray SlipLengths;
            SlipLengths = this.LsTrk.GridDat.Cells.h_min.CloneAs();
            SlipLengths.Clear();
            //SlipLengths.AccConstant(-1.0);
            if (SlipArea != null) {
                foreach (Chunk cnk in SlipArea) {
                    for (int i = cnk.i0; i < cnk.JE; i++) {
                        switch (this.dntParams.GNBC_SlipLength) {
                            case NavierSlip_SlipLength.hmin_DG: {
                                    int degU = ColMapping.BasisS.ToArray()[0].Degree;
                                    SlipLengths[i] = this.LsTrk.GridDat.Cells.h_min[i] / (degU + 1);
                                    break;
                                }
                            case NavierSlip_SlipLength.hmin_Grid: {
                                    SlipLengths[i] = SlipLengths[i] = this.LsTrk.GridDat.Cells.h_min[i];
                                    break;
                                }
                            case NavierSlip_SlipLength.Prescribed_SlipLength: {
                                    SlipLengths[i] = this.physParams.sliplength;
                                    break;
                                }
                            case NavierSlip_SlipLength.Prescribed_Beta: {
                                    SlipLengths[i] = -1.0;
                                    break;
                                }
                        }
                    }
                }

            }


            // interface coefficients 
            // ======================

            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengths = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), CutCellQuadOrder).CutCellMetrics.InterfaceArea;

            MultidimensionalArray sigmaMaxValue;
            sigmaMaxValue = this.LsTrk.GridDat.Cells.h_min.CloneAs();
            sigmaMaxValue.Clear();

            //double LevSet_Deg = ((LevelSet)this.LsTrk.LevelSets[0]).Basis.Degree + 1;

            foreach (Chunk cnk in this.LsTrk.Regions.GetCutCellMask()) {
                for (int i = cnk.i0; i < cnk.JE; i++) {

                    double ILen = InterfaceLengths.ElementAt(0).Value[i];
                    //ILen /= LevSet_Deg;
                    double sigmaILen_Max = (this.physParams.rho_A + this.physParams.rho_B)
                           * Math.Pow(ILen, 3) / (2 * Math.PI * dt.Pow2());

                    if (dntParams.SetSurfaceTensionMaxValue && (physParams.Sigma > sigmaILen_Max)) {
                        sigmaMaxValue[i] = sigmaILen_Max * 0.5;
                        //Console.WriteLine("set new sigma value: {0}; {1}", sigmaILen_Max, sigmaILen_Max/physParams.Sigma);
                    } else {
                        sigmaMaxValue[i] = this.physParams.Sigma * 0.5;
                    }

                }
            }

            // dissipative interface model / local stabilization

            MultidimensionalArray lambdaI, muI;
            //lambdaI = SlipLengths.CloneAs();
            //lambdaI.Clear();
            muI = SlipLengths.CloneAs();
            muI.Clear();
            //SinglePhaseField muI_DGField = new SinglePhaseField(new Basis(this.LsTrk.GridDat, 0));

            foreach (Chunk cnk in this.LsTrk.Regions.GetCutCellMask()) {
                for (int i = cnk.i0; i < cnk.JE; i++) {

                    double ILen = InterfaceLengths.ElementAt(0).Value[i];
                    double h_sigma = this.LsTrk.GridDat.Cells.h_min[i];
                    int LevSet_Deg = ((LevelSet)this.LsTrk.LevelSets[0]).Basis.Degree;
                    h_sigma /= LevSet_Deg;

                    if (ILen < h_sigma) {
                        double dt_sigmaILen = (this.physParams.rho_A + this.physParams.rho_B)
                           * Math.Pow(ILen, 3) / (2 * Math.PI * physParams.Sigma);
                        double dt_sigma = (this.physParams.rho_A + this.physParams.rho_B)
                           * Math.Pow(h_sigma, 3) / (2 * Math.PI * physParams.Sigma);

                        muI[i] = 0.5 * (physParams.mu_A + physParams.mu_B) * (1 - (dt_sigmaILen / dt_sigma));
                        //muI_DGField.SetMeanValue(i, muI[i]);
                    }
                }
            }


            // assemble the matrix & affine vector
            // ===================================

            BitArray EvapMicroRegion = this.LsTrk.GridDat.GetBoundaryCells().GetBitMask();
            EvapMicroRegion.SetAll(false);

            // compute matrix
            if (OpMatrix != null) {

                XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = this.m_XOp.GetMatrixBuilder(LsTrk, ColMapping, Params, RowMapping);

                foreach (var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                    this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["SlipLengths"] = SlipLengths;
                    this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["EvapMicroRegion"] = EvapMicroRegion;
                    if (config.prescribedMassflux != null) {
                        double[] dummyX = new double[] { 0.0, 0.0 };
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["prescribedMassflux"] = config.prescribedMassflux(dummyX, time);
                    }
                }

                if (this.m_XOp.SurfaceElementOperator.TotalNoOfComponents > 0) {
                    foreach (var kv in InterfaceLengths) {
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["InterfaceLengths"] = kv.Value;
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["sigmaMaxValue"] = sigmaMaxValue;
                        //this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["lambda_interface"] = lambdaI;
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["mu_interface"] = muI;
                    }
                }

                mtxBuilder.time = time;

                mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);

            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin eval = this.m_XOp.GetEvaluatorEx(this.LsTrk,
                    CurrentState.ToArray(), Params, RowMapping);

                foreach (var kv in AgglomeratedCellLengthScales) {
                    eval.CellLengthScales[kv.Key] = kv.Value;
                    this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["SlipLengths"] = SlipLengths;
                    this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["EvapMicroRegion"] = EvapMicroRegion;
                    if (config.prescribedMassflux != null) {
                        double[] dummyX = new double[] { 0.0, 0.0 };
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["prescribedMassflux"] = config.prescribedMassflux(dummyX, time);
                    }
                }

                if (this.m_XOp.SurfaceElementOperator.TotalNoOfComponents > 0) {
                    foreach (var kv in InterfaceLengths) {
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["InterfaceLengths"] = kv.Value;
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["sigmaMaxValue"] = sigmaMaxValue;
                        //this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["lambda_interface"] = lambdaI;
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["mu_interface"] = muI;
                    }
                }

                eval.time = time;

                eval.Evaluate(1.0, 1.0, OpAffine);

            }

        }


        DGField[] storedParams;



        private void ComputeAverageU<T>(IEnumerable<T> U0, VectorField<XDGField> U0mean, int order, XQuadSchemeHelper qh) where T : DGField {
            using (FuncTrace ft = new FuncTrace()) {

                var CC = this.LsTrk.Regions.GetCutCellMask();
                int D = this.LsTrk.GridDat.SpatialDimension;
                double minvol = Math.Pow(this.LsTrk.GridDat.Cells.h_minGlobal, D);


                //var qh = new XQuadSchemeHelper(agg);
                foreach (var Spc in this.LsTrk.SpeciesIdS) { // loop over species...
                    //var Spc = this.LsTrk.GetSpeciesId("B"); {
                    // shadow fields
                    DGField[] U0_Spc = U0.Select(U0_d => (U0_d is XDGField) ? ((DGField)((U0_d as XDGField).GetSpeciesShadowField(Spc))) : ((DGField)U0_d)).ToArray();
                    var U0mean_Spc = U0mean.Select(U0mean_d => U0mean_d.GetSpeciesShadowField(Spc)).ToArray();


                    // normal cells:
                    for (int d = 0; d < D; d++) {
                        U0mean_Spc[d].AccLaidBack(1.0, U0_Spc[d], this.LsTrk.Regions.GetSpeciesMask(Spc));
                    }

                    // cut cells
                    var scheme = qh.GetVolumeQuadScheme(Spc, IntegrationDomain: this.LsTrk.Regions.GetCutCellMask());
                    var rule = scheme.Compile(this.LsTrk.GridDat, order);
                    CellQuadrature.GetQuadrature(new int[] { D + 1 }, // vector components: ( avg_vel[0], ... , avg_vel[D-1], cell_volume )
                        this.LsTrk.GridDat,
                        rule,
                        delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            EvalResult.Clear();
                            for (int d = 0; d < D; d++)
                                U0_Spc[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                            var Vol = EvalResult.ExtractSubArrayShallow(-1, -1, D);
                            Vol.SetAll(1.0);
                        },
                        delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                            for (int i = 0; i < Length; i++) {
                                int jCell = i + i0;

                                double Volume = ResultsOfIntegration[i, D];
                                if (Math.Abs(Volume) < minvol * 1.0e-12) {
                                    // keep current value
                                    // since the volume of species 'Spc' in cell 'jCell' is 0.0, the value in this cell should have no effect
                                } else {
                                    for (int d = 0; d < D; d++) {
                                        double IntVal = ResultsOfIntegration[i, d];
                                        U0mean_Spc[d].SetMeanValue(jCell, IntVal / Volume);
                                    }
                                }

                            }
                        }).Execute();

                }

#if DEBUG
                {
                    var Uncut = LsTrk.Regions.GetCutCellMask().Complement();


                    VectorField<SinglePhaseField> U0mean_check = new VectorField<SinglePhaseField>(D, new Basis(LsTrk.GridDat, 0), SinglePhaseField.Factory);
                    for (int d = 0; d < D; d++) {
                        U0mean_check[d].ProjectField(1.0, U0.ElementAt(d).Evaluate,
                            new CellQuadratureScheme(false, Uncut).AddFixedOrderRules(LsTrk.GridDat, U0.ElementAt(d).Basis.Degree + 1));
                    }

                    foreach (var _Spc in this.LsTrk.SpeciesIdS) { // loop over species...
                        for (int d = 0; d < D; d++) {
                            U0mean_check[d].AccLaidBack(-1.0, U0mean[d].GetSpeciesShadowField(_Spc), Uncut.Intersect(LsTrk.Regions.GetSpeciesMask(_Spc)));
                        }
                    }

                    double checkNorm = U0mean_check.L2Norm();
                    Debug.Assert(checkNorm < 1.0e-6);
                }
#endif


                U0mean.ForEach(F => F.CheckForNanOrInf(true, true, true));

            }
        }

    }
}

