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
using BoSSS.Solution.RheologyCommon;
using System.Collections;
using BoSSS.Solution.Statistic;

namespace BoSSS.Application.XRheology_Solver {

    /// <summary>
    /// class for defining the equation components of the XRheology_Operator 
    /// and assembly of the corresponding matrix 
    /// </summary>
    public class XRheology_OperatorFactory : XOperatorFactoryBase {


        string[] CodNameSelected = new string[0];
        string[] DomNameSelected = new string[0];

        string[] Params;

        int HMFDegree;
        //int compRow;
        //int compCol;
        int stressDegree;

        PhysicalParameters physParams;
        DoNotTouchParameters dntParams;

        // Parameters: Velocity Gradient
        public VectorField<XDGField> VelocityXGradient;
        public VectorField<XDGField> VelocityYGradient;

        bool useJacobianForOperatorMatrix;

        bool useArtificialDiffusion;

        bool NormalsRequired;

        bool CurvatureRequired;

        bool U0meanrequired;

        /// <summary>
        /// ctor for the operator factory, where the equation compnents are set
        /// </summary>
        /// <param name="config"></param>
        /// <param name="_LsTrk"></param>
        /// <param name="_HMFdegree"></param>
        /// <param name="BcMap"></param>
        /// <param name="degU"></param>
        public XRheology_OperatorFactory(XRheology_OperatorConfiguration config, LevelSetTracker _LsTrk, int _HMFdegree, IncompressibleMultiphaseBoundaryCondMap BcMap, int _stressDegree, int degU) {


            this.LsTrk = _LsTrk;
            this.D = _LsTrk.GridDat.SpatialDimension;

            this.HMFDegree = _HMFdegree;

            this.physParams = config.getPhysParams;
            this.dntParams = config.getDntParams;
            this.useJacobianForOperatorMatrix = config.isUseJacobian;
            this.useArtificialDiffusion = config.isUseArtificialDiffusion;


            // test input
            // ==========
            {
                if (config.getDomBlocks.GetLength(0) != 3 || config.getCodBlocks.GetLength(0) != 3)
                    throw new ArgumentException();

                //if ((config.getPhysParams.Weissenberg_a <= 0) && (config.getPhysParams.Weissenberg_a <= 0)) {
                //    config.isOldroydB = false;
                //} else {
                //    if ((config.getPhysParams.mu_A <= 0) || (config.getPhysParams.mu_B <= 0))
                //        throw new ArgumentException();
                //}

                //if ((config.getPhysParams.reynolds_A <= 0) && (config.getPhysParams.reynolds_B <= 0)) {
                //    config.isViscous = false;
                //} else {  
                //    if ((config.getPhysParams.reynolds_A <= 0) || (config.getPhysParams.reynolds_B <= 0))
                //    throw new ArgumentException();

                //if ((config.getPhysParams.rho_A <= 0) || (config.getPhysParams.rho_B <= 0))
                //    throw new ArgumentException();

                if (_LsTrk.SpeciesNames.Count != 2)
                    throw new ArgumentException();
                if (!(_LsTrk.SpeciesNames.Contains("A") && _LsTrk.SpeciesNames.Contains("B")))
                    throw new ArgumentException();
            }

            // full operator:
            // ==============
            CodName = ArrayTools.Cat(EquationNames.MomentumEquations(this.D), EquationNames.ContinuityEquation, EquationNames.Constitutive(this.D));
            Params = ArrayTools.Cat(
                VariableNames.Velocity0Vector(this.D),
                VariableNames.Velocity0MeanVector(this.D),
                VariableNames.VelocityX_GradientVector(),
                VariableNames.VelocityY_GradientVector(),
                VariableNames.StressXXP,
                VariableNames.StressXYP,
                VariableNames.StressYYP,
                // "artificialViscosity",
                VariableNames.NormalVector(this.D),
                VariableNames.Curvature,
                VariableNames.SurfaceForceVector(this.D)
                );
            DomName = ArrayTools.Cat(VariableNames.VelocityVector(this.D), VariableNames.Pressure, VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressYY);


            // selected part:
            if (config.getCodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, this.D));
            if (config.getCodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(this.D, 1));
            if (config.getCodBlocks[2])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(this.D + 1, 3));

            if (config.getDomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, this.D));
            if (config.getDomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(this.D, 1));
            if (config.getDomBlocks[2])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(this.D + 1, 3));


            // create Operator
            // ===============
            m_XOp = new XDifferentialOperatorMk2(DomNameSelected, Params, CodNameSelected, (A, B, C) => _HMFdegree, this.LsTrk.SpeciesNames);

            // add components
            // ============================

            // species bulk components
            for (int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
                // Navier Stokes equations
                XOperatorComponentsFactory.AddSpeciesNSE(m_XOp, config, this.D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, LsTrk, out U0meanrequired);

                // continuity equation
                if (config.isContinuity)
                    XOperatorComponentsFactory.AddSpeciesContinuityEq(m_XOp, config, this.D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap);

                // constitutive equation
                XConstitutiveOperatorComponentsFactory.AddSpeciesConstitutive(m_XOp, config, this.D, stressDegree, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, LsTrk, out U0meanrequired);

            }

            //// interface components
            XOperatorComponentsFactory.AddInterfaceNSE(m_XOp, config, this.D, BcMap, LsTrk);    // surface stress tensor
            XOperatorComponentsFactory.AddSurfaceTensionForce(m_XOp, config, this.D, BcMap, LsTrk, degU, out NormalsRequired, out CurvatureRequired);     // surface tension force
            XConstitutiveOperatorComponentsFactory.AddInterfaceConstitutive(m_XOp, config, this.D, BcMap, LsTrk, out U0meanrequired); //constitutive eq at interfeac

            if (config.isContinuity)
                XOperatorComponentsFactory.AddInterfaceContinuityEq(m_XOp, config, this.D, LsTrk);       // continuity equation

            m_XOp.Commit();
        }


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
            VectorField<SinglePhaseField> LevelSetGradient, SinglePhaseField ExternalyProvidedCurvature, double[] currentWeissenberg,
            IEnumerable<T> CoupledCurrentState = null, IEnumerable<T> CoupledParams = null) where T : DGField {

            // checks:
            if (ColMapping.BasisS.Count != this.m_XOp.DomainVar.Count)
                throw new ArgumentException();
            if (RowMapping.BasisS.Count != this.m_XOp.CodomainVar.Count)
                throw new ArgumentException();

            int D = this.LsTrk.GridDat.SpatialDimension;
            if (CurrentState != null && CurrentState.Count() != (D + 4))
                throw new ArgumentException();

            if (OpMatrix == null && CurrentState == null)
                throw new ArgumentException();

            DGField[] U0;
            if (CurrentState != null)
                U0 = CurrentState.Take(D).ToArray();
            else
                U0 = null;

            DGField[] Stress0;
            Stress0 = CurrentState.Skip(D + 1).Take(3).ToArray();

            if (U0.Count() != D)
                throw new ArgumentException("Spatial dimesion and number of velocity parameter components does not match!");

            if (Stress0.Count() != D + 1)
                throw new ArgumentException("Spatial dimesion and number of stress parameter components does not match!");



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


            // Velocity and stresses for linearization

            // linearization velocity:
            DGField[] U0_U0mean;
            if (this.U0meanrequired) {
                XDGBasis U0meanBasis = new XDGBasis(this.LsTrk, 0);
                VectorField<XDGField> U0mean = new VectorField<XDGField>(D, U0meanBasis, "U0mean_", XDGField.Factory);

                U0_U0mean = ArrayTools.Cat<DGField>(U0, U0mean);
            } else {
                U0_U0mean = new DGField[2 * D];
            }


            if (VelocityXGradient == null) {
                VelocityXGradient = new VectorField<XDGField>(D, CurrentState.ElementAt(0).Basis, "VelocityX_Gradient", XDGField.Factory);
            }
            if (VelocityYGradient == null) {
                VelocityYGradient = new VectorField<XDGField>(D, CurrentState.ElementAt(1).Basis, "VelocityY_Gradient", XDGField.Factory);
            }


            // concatenate everything
            var Params = ArrayTools.Cat<DGField>(
                U0_U0mean,
                VelocityXGradient, 
                VelocityYGradient, 
                Stress0, 
                //artificialViscosity,
                Normals,
                Curvature,
                ((SurfaceForce != null) ? SurfaceForce.ToArray() : new SinglePhaseField[D]));


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

                if (!useJacobianForOperatorMatrix) {
                    XDifferentialOperatorMk2.XEvaluatorLinear mtxBuilder = this.m_XOp.GetMatrixBuilder(LsTrk, ColMapping, Params, RowMapping);
                    this.ParameterUpdate(CurrentState, Params, CutCellQuadOrder, AgglomeratedCellLengthScales);

                    foreach (var kv in AgglomeratedCellLengthScales) {
                        mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                        //mtxBuilder.CellLengthScales[kv.Key].EdgeLengthScales = kv.Value; // this.LsTrk.GridDat.Edges.h_max_Edge;
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["SlipLengths"] = SlipLengths;
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["EvapMicroRegion"] = EvapMicroRegion;
                    }

                    if (this.m_XOp.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0) {
                        foreach (var kv in InterfaceLengths) {
                            this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["InterfaceLengths"] = kv.Value;
                        }
                    }

                    foreach (var kv in AgglomeratedCellLengthScales) {
                        int id;
                        if (kv.Key == LsTrk.SpeciesIdS[0])
                            id = 0;
                        else 
                            id = 1;
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["Weissenbergnumber"] = currentWeissenberg[id];
                    }

                    mtxBuilder.time = time;
                    mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);


                } else {

                    throw new NotImplementedException("The FDbuilder for the Jacobian is missing for the XSpatial Opearator!");
                    // Finite Difference Linearization
                    //XSpatialOperatorMk2.XEvaluatorLinear FDBuilder = this.m_XOp.GetFDJacobianBuilder(domMap, Params, codMap, this.ParameterUpdate);

                    //foreach (var kv in AgglomeratedCellLengthScales) {
                    //    FDbuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                    //    FDbuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("SlipLengths", SlipLengths);
                    //    FDbuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("EvapMicroRegion", EvapMicroRegion);
                    //}

                    //if (this.m_XOp.SurfaceElementOperator.TotalNoOfComponents > 0) {
                    //    foreach (var kv in InterfaceLengths) {
                    //        FDbuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                    //    }
                    //}
                    //FDbuilder.time = time;

                    //FDbuilder.ComputeMatrix(OpMatrix, OpAffine);

                    //// FDJacobian has (Mx +b) as RHS, for unsteady calc. we must subtract Mx for real affine Vector!
                    //OpMatrix.SpMV(-1.0, new CoordinateVector(CurrentState), 1.0, OpAffine);

                    //foreach (var kv in AgglomeratedCellLengthScales) {
                    //    FDbuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("Weissenbergnumber", currentWeissenberg);
                    //}
                }

            } else {
                XDifferentialOperatorMk2.XEvaluatorNonlin eval = this.m_XOp.GetEvaluatorEx(this.LsTrk,
                    CurrentState.ToArray(), Params, RowMapping);

                foreach (var kv in AgglomeratedCellLengthScales) {
                    eval.CellLengthScales[kv.Key] = kv.Value;
                    //eval.SpeciesOperatorCoefficients[kv.Key].EdgeLengthScales = kv.Value; //this.LsTrk.GridDat.Edges.h_max_Edge;
                    this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["SlipLengths"] = SlipLengths;
                }

                if (this.m_XOp.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0) {
                    foreach (var kv in InterfaceLengths) {
                        this.m_XOp.UserDefinedValues[this.LsTrk.GetSpeciesName(kv.Key)]["InterfaceLengths"] = kv.Value;
                    }
                }

                foreach (var kv in AgglomeratedCellLengthScales) {
                    int id;
                    if (kv.Key == LsTrk.SpeciesIdS[0])
                        id = 0;
                    else
                        id = 1;
                    this.m_XOp.UserDefinedValues[LsTrk.GetSpeciesName(kv.Key)]["Weissenbergnumber"] = currentWeissenberg[id];
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

        void ParameterUpdate(IEnumerable<DGField> CurrentState, IEnumerable<DGField> ParameterVar, int CutCellQuadOrder, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales) {

            LevelSet Phi = (LevelSet)(this.LsTrk.LevelSets[0]);
            SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();

            DGField[] U0;
            if (CurrentState != null)
                U0 = CurrentState.Take(D).ToArray();
            else
                U0 = null;

            DGField[] Stress0;
            Stress0 = CurrentState.Skip(D + 1).Take(3).ToArray();

            if (U0.Count() != D)
                throw new ArgumentException("Spatial dimesion and number of velocity parameter components does not match!");

            if (Stress0.Count() != D + 1)
                throw new ArgumentException("Spatial dimesion and number of stress parameter components does not match!");

            // linearization velocity:
            DGField[] U0_U0mean;
            if (this.U0meanrequired) {
                XDGBasis U0meanBasis = new XDGBasis(this.LsTrk, 0);
                VectorField<XDGField> U0mean = new VectorField<XDGField>(D, U0meanBasis, "U0mean_", XDGField.Factory);

                U0_U0mean = ArrayTools.Cat<DGField>(U0, U0mean);

                U0mean.Clear();
                ComputeAverageU(U0, U0mean, CutCellQuadOrder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, CutCellQuadOrder, 1).XQuadSchemeHelper);
            } else {
                U0_U0mean = new DGField[2 * D];
            }


            //if (this.Control.SetParamsAnalyticalSol == false) {
            //SinglePhaseField[] __VelocityXGradient = ParameterVar.Skip(2 * D).Take(D).Select(f => f as SinglePhaseField).ToArray();
            //SinglePhaseField[] __VelocityYGradient = ParameterVar.Skip(3 * D).Take(D).Select(f => f as SinglePhaseField).ToArray();
            //Debug.Assert(ArrayTools.AreEqual(__VelocityXGradient, VelocityXGradient.ToArray(), (fa, fb) => object.ReferenceEquals(fa, fb)));
            //Debug.Assert(ArrayTools.AreEqual(__VelocityYGradient, VelocityYGradient.ToArray(), (fa, fb) => object.ReferenceEquals(fa, fb)));

            //if (VelocityXGradient == null) {
            //    VelocityXGradient = new VectorField<XDGField>(D, CurrentState.ElementAt(0).Basis, "VelocityX_Gradient", XDGField.Factory);
            //}
            VelocityXGradient.Clear();
            VelocityXGradient.GradientByFlux(1.0, U0[0]);
            //if (VelocityYGradient == null) {
            //    VelocityYGradient = new VectorField<XDGField>(D, CurrentState.ElementAt(1).Basis, "VelocityY_Gradient", XDGField.Factory);
            //}
            VelocityYGradient.Clear();
            VelocityYGradient.GradientByFlux(1.0, U0[1]);
            //}

            if (this.useArtificialDiffusion == true) {

                throw new NotImplementedException("artificial diffusion not jet fully implemented...");
                //SinglePhaseField __ArtificialViscosity = ParameterVar.Skip(5 * D + 1).Take(1).Select(f => f as SinglePhaseField).ToArray()[0];
                //if (!object.ReferenceEquals(this.artificalViscosity, __ArtificialViscosity))
                //    throw new ApplicationException();

                //ArtificialViscosity.ProjectArtificalViscosityToDGField(__ArtificialViscosity, perssonsensor, this.Control.SensorLimit, artificialMaxViscosity);
            }
        }

    }
}

