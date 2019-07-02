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

        /// <summary>
        /// ctor for the operator factory, where the equation compnents are set
        /// </summary>
        /// <param name="config"></param>
        /// <param name="_LsTrk"></param>
        /// <param name="_HMFdegree"></param>
        /// <param name="BcMap"></param>
        /// <param name="degU"></param>
        public XNSFE_OperatorFactory(XNSFE_OperatorConfiguration config, LevelSetTracker _LsTrk, int _HMFdegree, IncompressibleMultiphaseBoundaryCondMap BcMap, int degU) {

            this.LsTrk = _LsTrk;
            this.D = _LsTrk.GridDat.SpatialDimension;

            this.HMFDegree = _HMFdegree;

            this.physParams = config.getPhysParams;
            this.dntParams = config.getDntParams;


            // test input
            // ==========
            {
                if (config.getDomBlocks.GetLength(0) != 2 || config.getCodBlocks.GetLength(0) != 2)
                    throw new ArgumentException();

                if ((config.getPhysParams.mu_A <= 0) && (config.getPhysParams.mu_B <= 0)) {
                    config.isViscous = false;
                } else {
                    if ((config.getPhysParams.mu_A <= 0) || (config.getPhysParams.mu_B <= 0))
                        throw new ArgumentException();
                }

                if ((config.getPhysParams.rho_A <= 0) || (config.getPhysParams.rho_B <= 0))
                    throw new ArgumentException();

                if (_LsTrk.SpeciesNames.Count != 2)
                    throw new ArgumentException();
                if (!(_LsTrk.SpeciesNames.Contains("A") && _LsTrk.SpeciesNames.Contains("B")))
                    throw new ArgumentException();
            }

            // full operator:
            // ==============
            CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D)).Cat("div");
            Params = ArrayTools.Cat(
                VariableNames.Velocity0Vector(D),
                VariableNames.Velocity0MeanVector(D),
                "Curvature",
                (new string[] { "surfForceX", "surfForceY", "surfForceZ" }).GetSubVector(0, D),
                (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D),
                (new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D)),
                VariableNames.Temperature,
                "DisjoiningPressure");
            DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure);

            // selected part:
            if (config.getCodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, D));
            if (config.getCodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D, 1));

            if (config.getDomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, D));
            if (config.getDomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D, 1));


            // create Operator
            // ===============
            m_XOp = new XSpatialOperatorMk2(DomNameSelected, Params, CodNameSelected, (A, B, C) => _HMFdegree, this.LsTrk.SpeciesIdS.ToArray());

            // add Navier-Stokes components
            // ============================

            // species bulk components
            for (int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
                // Navier Stokes equations
                Solution.XNSECommon.XOperatorComponentsFactory.AddSpeciesNSE(m_XOp, CodName.GetSubVector(0, D), LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, config, LsTrk, out U0meanrequired);

                // continuity equation
                if (config.isContinuity)
                    Solution.XNSECommon.XOperatorComponentsFactory.AddSpeciesContinuityEq(m_XOp, CodName[D], D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, config, LsTrk);
            }

            // interface components
            Solution.XNSECommon.XOperatorComponentsFactory.AddInterfaceNSE(m_XOp, CodName.GetSubVector(0, D), BcMap, config, LsTrk);    // surface stress tensor
            Solution.XNSECommon.XOperatorComponentsFactory.AddSurfaceTensionForce(m_XOp, CodName.GetSubVector(0, D), BcMap, config, LsTrk, degU, out NormalsRequired, out CurvatureRequired);     // surface tension force

            if (config.isContinuity)
                Solution.XNSECommon.XOperatorComponentsFactory.AddInterfaceContinuityEq(m_XOp, CodName[D], D, BcMap, config, LsTrk);       // continuity equation


            // add Evaporation interface components
            // ====================================

            if (config.isEvaporation) {

                XOperatorComponentsFactory.AddInterfaceNSE_withEvaporation(m_XOp, CodName.GetSubVector(0, D), config, LsTrk);

                if (config.isContinuity)
                    XOperatorComponentsFactory.AddInterfaceContinuityEq_withEvaporation(m_XOp, CodName[D], D, config, LsTrk);

            }


            m_XOp.Commit();
        }

        PhysicalParameters physParams;
        DoNotTouchParameters dntParams;

        bool NormalsRequired;

        bool CurvatureRequired;

        bool U0meanrequired;

        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="OpMatrix"></param>
        /// <param name="OpAffine"></param>
        /// <param name="RowMapping"></param>
        /// <param name="ColMapping"></param>
        /// <param name="CurrentState"></param>
        /// <param name="AgglomeratedCellLengthScales"></param>
        /// <param name="time"></param>
        /// <param name="CutCellQuadOrder"></param>
        /// <param name="SurfaceForce"></param>
        /// <param name="LevelSetGradient"></param>
        /// <param name="ExternalyProvidedCurvature"></param>
        public void AssembleMatrix<T>(BlockMsrMatrix OpMatrix, double[] OpAffine,
            UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping,
            IEnumerable<T> CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time,
            int CutCellQuadOrder, VectorField<SinglePhaseField> SurfaceForce,
            VectorField<SinglePhaseField> LevelSetGradient, SinglePhaseField ExternalyProvidedCurvature,
            IEnumerable<T> CoupledCurrentState = null, IEnumerable<T> CoupledParams = null) where T : DGField {

            // checks:
            if (ColMapping.BasisS.Count != this.m_XOp.DomainVar.Count)
                throw new ArgumentException();
            if (RowMapping.BasisS.Count != this.m_XOp.CodomainVar.Count)
                throw new ArgumentException();

            int D = this.LsTrk.GridDat.SpatialDimension;
            if (CurrentState != null && CurrentState.Count() != (D + 1))
                throw new ArgumentException();

            if (OpMatrix == null && CurrentState == null)
                throw new ArgumentException();

            DGField[] U0;
            if (CurrentState != null)
                U0 = CurrentState.Take(D).ToArray();
            else
                U0 = null;



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


            // parameter assembly
            // ==================

            LevelSet Phi = (LevelSet)(this.LsTrk.LevelSets[0]);
            SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();

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

            // linearization velocity:
            DGField[] U0_U0mean;
            if (this.U0meanrequired) {
                XDGBasis U0meanBasis = new XDGBasis(this.LsTrk, 0);
                VectorField<XDGField> U0mean = new VectorField<XDGField>(D, U0meanBasis, "U0mean_", XDGField.Factory);

                U0_U0mean = ArrayTools.Cat<DGField>(U0, U0mean);
            } else {
                U0_U0mean = new DGField[2 * D];
            }

            // Temperature gradient for evaporation
            VectorField<DGField> GradTemp = new VectorField<DGField>(D, new XDGBasis(LsTrk, 0), XDGField.Factory);
            if (CoupledCurrentState != null) {
                DGField Temp = CoupledCurrentState.ToArray()[0];
                GradTemp = new VectorField<DGField>(D, Temp.Basis, "GradTemp", XDGField.Factory);
                XNSEUtils.ComputeGradientForParam(Temp, GradTemp, this.LsTrk);
            }

            // concatenate everything
            var Params = ArrayTools.Cat<DGField>(
                U0_U0mean,
                Curvature,
                ((SurfaceForce != null) ? SurfaceForce.ToArray() : new SinglePhaseField[D]),
                Normals,
                ((CoupledCurrentState != null) ? GradTemp.ToArray() : new SinglePhaseField[D]),
                ((CoupledCurrentState != null) ? CoupledCurrentState.ToArray<DGField>() : new SinglePhaseField[1]),
                ((CoupledCurrentState != null) ? CoupledParams.ToArray<DGField>() : new SinglePhaseField[1]));  //((evaporation) ? GradTemp.ToArray() : new SinglePhaseField[D]));

            // linearization velocity:
            if (this.U0meanrequired) {
                VectorField<XDGField> U0mean = new VectorField<XDGField>(U0_U0mean.Skip(D).Take(D).Select(f => ((XDGField)f)).ToArray());

                U0mean.Clear();
                if (this.physParams.IncludeConvection)
                    ComputeAverageU(U0, U0mean, CutCellQuadOrder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, CutCellQuadOrder, 1).XQuadSchemeHelper);
            }



            // assemble the matrix & affine vector
            // ===================================

            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengths = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), CutCellQuadOrder).CutCellMetrics.InterfaceArea;

            BitArray EvapMicroRegion = this.LsTrk.GridDat.GetBoundaryCells().GetBitMask();
            EvapMicroRegion.SetAll(false);

            // compute matrix
            if (OpMatrix != null) {

                XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = this.m_XOp.GetMatrixBuilder(LsTrk, ColMapping, Params, RowMapping, SpcToCompute);

                foreach (var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("SlipLengths", SlipLengths);
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("EvapMicroRegion", EvapMicroRegion);
                }

                if (this.m_XOp.SurfaceElementOperator.TotalNoOfComponents > 0) {
                    foreach (var kv in InterfaceLengths) {
                        mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                    }
                }

                mtxBuilder.time = time;

                mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);

            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin eval = this.m_XOp.GetEvaluatorEx(this.LsTrk,
                    CurrentState.ToArray(), Params, RowMapping,
                    SpcToCompute);

                foreach (var kv in AgglomeratedCellLengthScales) {
                    eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                    eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("SlipLengths", SlipLengths);
                    eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("EvapMicroRegion", EvapMicroRegion);
                }

                if (this.m_XOp.SurfaceElementOperator.TotalNoOfComponents > 0) {
                    foreach (var kv in InterfaceLengths) {
                        eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                    }
                }

                eval.time = time;

                eval.Evaluate(1.0, 1.0, OpAffine);

            }

        }


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

