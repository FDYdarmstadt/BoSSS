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
using System.Linq;
using System.Text;
using BoSSS.Foundation.XDG;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon.Operator;
using BoSSS.Foundation;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Tracing;
using MPI.Wrappers;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using ilPSP;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using NUnit.Framework;

namespace BoSSS.Solution.XNSECommon {

    public class OperatorFactory {

        LevelSetTracker LsTrk;

        DoNotTouchParameters dntParams;

        PhysicalParameters physParams;

        int D;

        bool muIs0 = false;

        string[] DomName;
        string[] CodName;
        string[] Params;

        string[] CodNameSelected = new string[0];
        string[] DomNameSelected = new string[0];

        double muA;
        double muB;
        double rhoA;
        double rhoB;
        double sigma;
        bool MatInt;
        bool UseExtendedVelocity;

        public bool movingmesh;

        //XQuadFactoryHelper.MomentFittingVariants momentFittingVariant;
        int HMFDegree;

        public OperatorFactory(
            OperatorConfiguration config,
            LevelSetTracker _LsTrk,
            int _HMFdegree,
            int degU,
            IncompressibleMultiphaseBoundaryCondMap BcMap,
            bool _movingmesh) {

            // variable names
            // ==============
            D = _LsTrk.GridDat.SpatialDimension;
            this.LsTrk = _LsTrk;
            //this.momentFittingVariant = momentFittingVariant;
            this.HMFDegree = _HMFdegree;
            this.dntParams = config.dntParams.CloneAs();
            this.physParams = config.physParams.CloneAs();
            this.UseExtendedVelocity = config.UseXDG4Velocity;
            this.movingmesh = _movingmesh;
            // test input
            // ==========
            {
                if (config.DomBlocks.GetLength(0) != 2 || config.CodBlocks.GetLength(0) != 2)
                    throw new ArgumentException();
                if (config.physParams.mu_A == 0.0 && config.physParams.mu_B == 0.0) {
                    muIs0 = true;
                } else {
                    if (config.physParams.mu_A <= 0)
                        throw new ArgumentException();
                    if (config.physParams.mu_B <= 0)
                        throw new ArgumentException();
                }

                if (config.physParams.rho_A <= 0)
                    throw new ArgumentException();
                if (config.physParams.rho_B <= 0)
                    throw new ArgumentException();

                if (_LsTrk.SpeciesNames.Count != 2)
                    throw new ArgumentException();
                if (!(_LsTrk.SpeciesNames.Contains("A") && _LsTrk.SpeciesNames.Contains("B")))
                    throw new ArgumentException();
            }


            // full operator:
            CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D)).Cat("div");
            Params = ArrayTools.Cat(
                VariableNames.Velocity0Vector(D),
                VariableNames.Velocity0MeanVector(D),
                "Curvature",
                (new string[] { "surfForceX", "surfForceY", "surfForceZ" }).GetSubVector(0, D),
                (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D)
                );
            DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure);

            // selected part:
            if (config.CodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, D));
            if (config.CodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D, 1));

            if (config.DomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, D));
            if (config.DomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D, 1));

            muA = config.physParams.mu_A;
            muB = config.physParams.mu_B;
            rhoA = config.physParams.rho_A;
            rhoB = config.physParams.rho_B;
            sigma = config.physParams.Sigma;

            MatInt = config.physParams.Material;

            double MassFlux = -0.1;

            //if (!MatInt)
            //    throw new NotSupportedException("Non-Material interface is NOT tested!");

            // create Operator
            // ===============
            m_OP = new XSpatialOperator(DomNameSelected, Params, CodNameSelected, (A,B,C) => _HMFdegree);

            // build the operator
            // ==================
            {

                // Momentum equation
                // =================

                if (config.physParams.IncludeConvection && config.Transport) {

                    for (int d = 0; d < D; d++) {
                        var comps = m_OP.EquationComponents[CodName[d]];

                        // convective part:
                        // variante 1: 

                        double LFFA = config.dntParams.LFFA;
                        double LFFB = config.dntParams.LFFB;


                        var conv = new Operator.Convection.ConvectionInBulk_LLF(D, BcMap, d, rhoA, rhoB, LFFA, LFFB, LsTrk);
                        comps.Add(conv); // Bulk component
                        comps.Add(new Operator.Convection.ConvectionAtLevelSet_LLF(d, D, LsTrk, rhoA, rhoB, LFFA, LFFB, config.physParams.Material, BcMap, movingmesh));       // LevelSet component


                        // variante 3:
                        //var convA = new LocalConvection(D, d, rhoA, rhoB, this.config.varMode, LsTrk);
                        //XOP.OnIntegratingBulk += convA.SetParameter;
                        //comps.Add(convA);
                        //////var convB = new LocalConvection2(D, d, rhoA, rhoB, varMode, LsTrk); // macht Bum-Bum!
                        ////XOP.OnIntegratingBulk += convB.SetParameter;
                        ////comps.Add(convB);
                    }

                    this.U0meanrequired = true;
                }

                // pressure gradient
                // =================

                if (config.PressureGradient) {

                    for (int d = 0; d < D; d++) {
                        var comps = m_OP.EquationComponents[CodName[d]];


                        var pres = new Operator.Pressure.PressureInBulk(d, BcMap);
                        comps.Add(pres);

                        //if (!MatInt)
                        //    throw new NotSupportedException("New Style pressure coupling does not support non-material interface.");
                        var presLs = new Operator.Pressure.PressureFormAtLevelSet(d, D, LsTrk);
                        comps.Add(presLs);
                    }
                }

                // viscous operator
                // ================

                if (config.Viscous && !muIs0) {


                    for (int d = 0; d < D; d++) {
                        var comps = m_OP.EquationComponents[CodName[d]];
                        // viscous part:
                        //double _D = D;
                        //double penalty_mul = dntParams.PenaltySafety;
                        //double _p = degU;
                        //double penalty_base = (_p + 1) * (_p + _D) / _D;
                        //double penalty = penalty_base * penalty_mul;
                        double penalty = dntParams.PenaltySafety;
                        switch (dntParams.ViscosityMode) {
                            case ViscosityMode.Standard: {
                                    // Bulk operator:
                                    var Visc = new Operator.Viscosity.ViscosityInBulk_GradUTerm(
                                        dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                        BcMap, d, D, muA, muB); // , _betaA: this.physParams.betaS_A, _betaB: this.physParams.betaS_B);

                                    comps.Add(Visc);

                                    if (dntParams.UseGhostPenalties) {
                                        var ViscPenalty = new Operator.Viscosity.ViscosityInBulk_GradUTerm(penalty * 1.0, 0.0, BcMap, d, D, muA, muB);
                                        m_OP.GhostEdgesOperator.EquationComponents[CodName[d]].Add(ViscPenalty);
                                    }
                                    // Level-Set operator:
                                    comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(LsTrk, muA, muB, penalty * 1.0, d, true));

                                    break;
                                }
                            case ViscosityMode.TransposeTermMissing: {
                                    // Bulk operator:
                                    var Visc = new Operator.Viscosity.ViscosityInBulk_GradUTerm(
                                        dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                        BcMap, d, D, muA, muB);
                                    comps.Add(Visc);
                                    
                                    if (dntParams.UseGhostPenalties) {
                                        var ViscPenalty = new Operator.Viscosity.ViscosityInBulk_GradUTerm(penalty * 1.0, 0.0, BcMap, d, D, muA, muB);
                                        m_OP.GhostEdgesOperator.EquationComponents[CodName[d]].Add(ViscPenalty);
                                    }
                                    // Level-Set operator:
                                    comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(LsTrk, muA, muB, penalty * 1.0, d, false));

                                    break;
                                }
                            case ViscosityMode.ExplicitTransformation: {
                                    // Bulk operator
                                    var Visc = new Operator.Viscosity.ViscosityInBulk_GradUTerm(
                                        dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                        BcMap, d, D, muA, muB);
                                    comps.Add(Visc);
                                    if (dntParams.UseGhostPenalties) {
                                        var ViscPenalty = new Operator.Viscosity.ViscosityInBulk_GradUTerm(penalty * 1.0, 0.0, BcMap, d, D, muA, muB);
                                        m_OP.GhostEdgesOperator.EquationComponents[CodName[d]].Add(ViscPenalty);
                                    }

                                    //Level-Set operator:
                                    //comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Explicit(d, D, LsTrk, penalty, muA, muB));
                                    throw new NotSupportedException("Beim refact rausgeflogen, braucht eh kein Mensch. fk, 08jan16.");

                                    //break;
                                }

                            case ViscosityMode.FullySymmetric: {
                                    // Bulk operator
                                    var Visc1 = new Operator.Viscosity.ViscosityInBulk_GradUTerm(
                                        dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                        BcMap, d, D, muA, muB, _betaA: this.physParams.betaS_A, _betaB: this.physParams.betaS_B);
                                    var Visc2 = new Operator.Viscosity.ViscosityInBulk_GradUtranspTerm(
                                        dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                        BcMap, d, D, muA, muB, _betaA: this.physParams.betaS_A, _betaB: this.physParams.betaS_B);
                                    //var Visc3 = new Operator.Viscosity.ViscosityInBulk_divTerm(dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0, BcMap, d, D, muA, muB);


                                    comps.Add(Visc1);
                                    comps.Add(Visc2);
                                    //comps.Add(Visc3);

                                    if (dntParams.UseGhostPenalties) {
                                        var Visc1Penalty = new Operator.Viscosity.ViscosityInBulk_GradUTerm(
                                            penalty, 0.0,
                                            BcMap, d, D, muA, muB);
                                        var Visc2Penalty = new Operator.Viscosity.ViscosityInBulk_GradUtranspTerm(
                                            penalty, 0.0,
                                            BcMap, d, D, muA, muB);
                                        //var Visc3Penalty = new Operator.Viscosity.ViscosityInBulk_divTerm(
                                        //    penalty, 0.0,
                                        //    BcMap, d, D, muA, muB);
                                        //m_OP.OnIntegratingBulk += Visc3Penalty.SetParameter;
                                        m_OP.GhostEdgesOperator.EquationComponents[CodName[d]].Add(Visc1Penalty);
                                        m_OP.GhostEdgesOperator.EquationComponents[CodName[d]].Add(Visc2Penalty);
                                        //m_OP.AndresHint.EquationComponents[CodName[d]].Add(Visc3Penalty);
                                    }

                                    // Level-Set operator
                                    comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_FullySymmetric(LsTrk, muA, muB, penalty, d));

                                    //comps.Add(new Operator.Viscosity.GeneralizedViscosityAtLevelSet_FullySymmetric(LsTrk, muA, muB, penalty, d, rhoA, rhoB, MassFlux));

                                    break;
                                }

                            default:
                                throw new NotImplementedException();
                        }
                    }

                }

                // Continuum equation
                // ==================

                if (config.continuity) {

                    for (int d = 0; d < D; d++) {
                        var src = new Operator.Continuity.DivergenceInBulk_Volume(d, D, rhoA, rhoB, config.dntParams.ContiSign, config.dntParams.RescaleConti);
                        var flx = new Operator.Continuity.DivergenceInBulk_Edge(d, BcMap, rhoA, rhoB, config.dntParams.ContiSign, config.dntParams.RescaleConti);
                        m_OP.EquationComponents["div"].Add(flx);
                        m_OP.EquationComponents["div"].Add(src);
                    }

                    var divPen = new Operator.Continuity.DivergenceAtLevelSet(D, LsTrk, rhoA, rhoB, MatInt, config.dntParams.ContiSign, config.dntParams.RescaleConti);
                    m_OP.EquationComponents["div"].Add(divPen);

                    //var divPenGen = new Operator.Continuity.GeneralizedDivergenceAtLevelSet(D, LsTrk, rhoA, rhoB, MatInt, config.dntParams.ContiSign, config.dntParams.RescaleConti, MassFlux);
                    //m_OP.EquationComponents["div"].Add(divPenGen);

                    //// pressure stabilization
                    //if (this.config.PressureStab) {
                    //    Console.WriteLine("Pressure Stabilization active.");
                    //var pStabi = new PressureStabilization(0.0001, 0.0001);
                    //    m_OP.OnIntegratingBulk += pStabi.SetParameter;
                    //    m_OP.EquationComponents["div"].Add(pStabi);
                    ////var pStabiLS = new PressureStabilizationAtLevelSet(D, LsTrk, rhoA, muA, rhoB, muB, sigma, this.config.varMode, MatInt);
                    //    //XOP.EquationComponents["div"].Add(pStabiLS);
                    //} else {
                    //    Console.WriteLine("Pressure Stabilization INACTIVE.");
                    //}

                }

                // surface tension
                // ===============

                if(config.PressureGradient && config.physParams.Sigma != 0.0) {

                    // isotropic part of the surface stress tensor
                    if(config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux
                     || config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local
                     || config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {

                        for(int d = 0; d < D; d++) {

                            if(config.dntParams.SST_isotropicMode != SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                                IEquationComponent G = new SurfaceTension_LaplaceBeltrami_Surface(d, config.physParams.Sigma * 0.5);
                                IEquationComponent H = new SurfaceTension_LaplaceBeltrami_BndLine(d, config.physParams.Sigma * 0.5, config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);
                                m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(G);
                                m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(H);
                            } else {
                                //G = new SurfaceTension_LaplaceBeltrami2_Surface(d, config.physParams.Sigma * 0.5);
                                //H = new SurfaceTension_LaplaceBeltrami2_BndLine(d, config.physParams.Sigma * 0.5, config.physParams.Theta_e, config.physParams.betaL);
                                IEquationComponent isoSurfT = new IsotropicSurfaceTension_LaplaceBeltrami(d, D, config.physParams.Sigma * 0.5, BcMap.EdgeTag2Type, config.physParams.theta_e, config.physParams.betaL);
                                m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(isoSurfT);
                            }

                        }

                        this.NormalsRequired = true;


                    } else if(config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Projected
                            || config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint
                            || config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean
                            || config.dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Fourier) {

                        for(int d = 0; d < D; d++) {
                            m_OP.EquationComponents[CodName[d]].Add(new CurvatureBasedSurfaceTension(d, D, LsTrk, config.physParams.Sigma));
                        }

                        this.CurvatureRequired = true;

                        /*
                        Console.WriteLine("REM: hack in Operator factory");
                        for(int d = 0; d < D; d++) {
                            //var G = new SurfaceTension_LaplaceBeltrami_Surface(d, config.physParams.Sigma * 0.5);
                            var H = new SurfaceTension_LaplaceBeltrami_BndLine(d, config.physParams.Sigma * 0.5, config.dntParams.surfTensionMode == SurfaceTensionMode.LaplaceBeltrami_Flux);

                            //m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(G);
                            m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(H);
                        }

                        this.NormalsRequired = true;
                        */

                    } else {
                        throw new NotImplementedException("Not implemented.");
                    }


                    // dynamic part
                    if(config.dntParams.SurfStressTensor != SurfaceSressTensor.Isotropic) {

                        double muI = config.physParams.mu_I;
                        double lamI = config.physParams.lambda_I;

                        double penalty_base = (degU + 1) * (degU + D) / D;
                        double penalty = penalty_base * dntParams.PenaltySafety;

                        // surface shear viscosity 
                        if(config.dntParams.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation ||
                            config.dntParams.SurfStressTensor == SurfaceSressTensor.SemiImplicit ||
                            config.dntParams.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                            for(int d = 0; d < D; d++) {
                                var surfDeformRate = new BoussinesqScriven_SurfaceDeformationRate_GradU(d, muI * 0.5, penalty);
                                m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(surfDeformRate);
                                //m_OP.OnIntegratingSurfaceElement += surfDeformRate.SetParameter;

                                if(config.dntParams.SurfStressTensor != SurfaceSressTensor.SemiImplicit) {
                                    var surfDeformRateT = new BoussinesqScriven_SurfaceDeformationRate_GradUTranspose(d, muI * 0.5, penalty);
                                    m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(surfDeformRateT);
                                    //m_OP.OnIntegratingSurfaceElement += surfDeformRateT.SetParameter;
                                }
                            }

                        }
                        // surface dilatational viscosity
                        if(config.dntParams.SurfStressTensor == SurfaceSressTensor.SurfaceVelocityDivergence ||
                            config.dntParams.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                            for(int d = 0; d < D; d++) {
                                var surfVelocDiv = new BoussinesqScriven_SurfaceVelocityDivergence(d, muI * 0.5, lamI * 0.5, penalty, BcMap.EdgeTag2Type);
                                m_OP.SurfaceElementOperator.EquationComponents[CodName[d]].Add(surfVelocDiv);
                                //m_OP.OnIntegratingSurfaceElement += surfVelocDiv.SetParameter;
                            }

                        }
                    }


                    // stabilization
                    if(config.dntParams.UseLevelSetStabilization) {

                        for(int d = 0; d < D; d++) {
                            m_OP.EquationComponents[CodName[d]].Add(new LevelSetStabilization(d, D, LsTrk));
                        }
                    }

                }


                // surface force term
                // ==================


                if (config.PressureGradient && config.physParams.useArtificialSurfaceForce) {
                    for (int d = 0; d < D; d++) {
                        m_OP.EquationComponents[CodName[d]].Add(new SurfaceTension_ArfForceSrc(d, D, LsTrk));
                    }
                }


                // evaporation (mass flux)
                // =======================

                //for(int d = 0; d < D; d++) {
                //    m_OP.EquationComponents[CodName[d]].Add(new Operator.DynamicInterfaceConditions.PrescribedMassFlux(d, D, LsTrk, rhoA, rhoB, MassFlux));
                //}


            }

            // Finalize
            // ========

            m_OP.Commit();
        }
        

       
        
        
        XSpatialOperator m_OP;

        /// <summary>
        /// The incompressible XNSE-Operator 
        /// </summary>
        public XSpatialOperator Op {
            get {
                return m_OP;
            }
        }
        
        bool NormalsRequired = false;

        bool CurvatureRequired = false;

        bool U0meanrequired = false;

        

        /// <summary>
        /// 
        /// </summary>
        public void AssembleMatrix_Timestepper<T>(
            int CutCellQuadOrder,
            BlockMsrMatrix OpMatrix, double[] OpAffine, 
            Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales,
            IEnumerable<T> CurrentState,
            VectorField<SinglePhaseField> SurfaceForce,
            VectorField<SinglePhaseField> LevelSetGradient, SinglePhaseField ExternalyProvidedCurvature,
            UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping,
            double time) where T : DGField {

            if (ColMapping.BasisS.Count != this.Op.DomainVar.Count)
                throw new ArgumentException();
            if (RowMapping.BasisS.Count != this.Op.CodomainVar.Count)
                throw new ArgumentException();

            // check:
            var Tracker = this.LsTrk;
            int D = Tracker.GridDat.SpatialDimension;
            if (CurrentState != null && CurrentState.Count() != (D + 1))
                throw new ArgumentException();
            if (OpMatrix == null && CurrentState == null)
                throw new ArgumentException();
            DGField[] U0;
            if (CurrentState != null)
                U0 = CurrentState.Take(D).ToArray();
            else
                U0 = null;

            
                       
            LevelSet Phi = (LevelSet)(Tracker.LevelSets[0]);

            SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();

            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengths = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), CutCellQuadOrder).CutCellMetrics.InterfaceArea;


            // advanced settings for the navier slip boundary condition
            // ========================================================

            CellMask SlipArea;
            switch(this.dntParams.GNBC_Localization) {
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

            if(SlipArea != null) {
                foreach(Chunk cnk in SlipArea) {
                    for(int i = cnk.i0; i < cnk.JE; i++) {
                        switch(this.dntParams.GNBC_SlipLength) {
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
                XDGBasis U0meanBasis = new XDGBasis(Tracker, 0);
                VectorField<XDGField> U0mean = new VectorField<XDGField>(D, U0meanBasis, "U0mean_", XDGField.Factory);

                U0_U0mean = ArrayTools.Cat<DGField>(U0, U0mean);
            } else {
                U0_U0mean = new DGField[2 * D];
            }

            // concatenate everything
            var Params = ArrayTools.Cat<DGField>(
                U0_U0mean,
                Curvature,
                ((SurfaceForce != null) ? SurfaceForce.ToArray() : new SinglePhaseField[D]),
                Normals);

            // linearization velocity:
            if (this.U0meanrequired) {
                VectorField<XDGField> U0mean = new VectorField<XDGField>(U0_U0mean.Skip(D).Take(D).Select(f => ((XDGField)f)).ToArray());

                U0mean.Clear();
                if (this.physParams.IncludeConvection)
                    ComputeAverageU(U0, U0mean, CutCellQuadOrder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, CutCellQuadOrder, 1).XQuadSchemeHelper );
            }

            // assemble the matrix & affine vector
            // ===================================

            // compute matrix
            if (OpMatrix != null) {
                //Op.ComputeMatrixEx(Tracker,
                //    ColMapping, Params, RowMapping,
                //    OpMatrix, OpAffine, false, time, true,
                //    AgglomeratedCellLengthScales,
                //    InterfaceLengths, SlipLengths,
                //    SpcToCompute);

                XSpatialOperator.XEvaluatorLinear mtxBuilder = Op.GetMatrixBuilder(LsTrk, ColMapping, Params, RowMapping, SpcToCompute);

                foreach(var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                    mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("SlipLengths", SlipLengths);
                }

                if(Op.SurfaceElementOperator.TotalNoOfComponents > 0) {
                    foreach(var kv in InterfaceLengths)
                        mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                }

                mtxBuilder.time = time;

                mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);

            } else {
                XSpatialOperator.XEvaluatorNonlin eval = Op.GetEvaluatorEx(Tracker,
                    CurrentState.ToArray(), Params, RowMapping,
                    SpcToCompute);

                foreach(var kv in AgglomeratedCellLengthScales) {
                    eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
                    eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("SlipLengths", SlipLengths);
                }

                if(Op.SurfaceElementOperator.TotalNoOfComponents > 0) {
                    foreach(var kv in InterfaceLengths)
                        eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                }

                eval.time = time;

                eval.Evaluate(1.0, 1.0, OpAffine);

#if DEBUG
                // remark: remove this piece in a few months from now on (09may18) if no problems occur
                {

                    BlockMsrMatrix checkOpMatrix = new BlockMsrMatrix(RowMapping, ColMapping);
                    double[] checkAffine = new double[OpAffine.Length];

                    Op.ComputeMatrixEx(Tracker,
                    ColMapping, Params, RowMapping,
                    OpMatrix, OpAffine, false, time, true,
                    AgglomeratedCellLengthScales,
                    InterfaceLengths, SlipLengths,
                    SpcToCompute);


                    double[] checkResult = checkAffine.CloneAs();
                    var currentVec = new CoordinateVector(CurrentState.ToArray());
                    checkOpMatrix.SpMV(1.0, new CoordinateVector(CurrentState.ToArray()), 1.0, checkResult);

                    double L2_dist = GenericBlas.L2DistPow2(checkResult, OpAffine).MPISum().Sqrt();
                    double RefNorm = (new double[] { checkResult.L2NormPow2(), OpAffine.L2NormPow2(), currentVec.L2NormPow2() }).MPISum().Max().Sqrt();

                    Assert.LessOrEqual(L2_dist, RefNorm * 1.0e-6);
                    Debug.Assert(L2_dist < RefNorm * 1.0e-6);
                }
#endif
            }


            // check
            // =====

            /*
            {
                DGField[] testDomainFieldS = ColMapping.BasisS.Select(bb => new XDGField(bb as XDGBasis)).ToArray();
                CoordinateVector test = new CoordinateVector(testDomainFieldS);

                DGField[] errFieldS = ColMapping.BasisS.Select(bb => new XDGField(bb as XDGBasis)).ToArray();
                CoordinateVector Err = new CoordinateVector(errFieldS);

                var eval = Op.GetEvaluatorEx(LsTrk,
                    testDomainFieldS, Params, RowMapping);

                foreach (var s in this.LsTrk.SpeciesIdS)
                    eval.SpeciesOperatorCoefficients[s].CellLengthScales = AgglomeratedCellLengthScales[s];
                
                eval.time = time;
                int L = test.Count;
                Random r = new Random();
                for(int i = 0; i < L; i++) {
                    test[i] = r.NextDouble();
                }
                


                double[] R1 = new double[L];
                double[] R2 = new double[L];
                eval.Evaluate(1.0, 1.0, R1);

                R2.AccV(1.0, OpAffine);
                OpMatrix.SpMV(1.0, test, 1.0, R2);

                Err.AccV(+1.0, R1);
                Err.AccV(-1.0, R2);

                double ErrDist = GenericBlas.L2DistPow2(R1, R2).MPISum().Sqrt();

                double Ref = test.L2NormPow2().MPISum().Sqrt();

                Debug.Assert(ErrDist <= Ref*1.0e-5, "Mismatch between explicit evaluation of XDG operator and matrix.");
            }
            */

        }


        private void ComputeAverageU<T>(IEnumerable<T> U0, VectorField<XDGField> U0mean, int order, XQuadSchemeHelper qh) where T: DGField {
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
                    var rule = scheme.Compile(this.LsTrk.GridDat,  order);
                    CellQuadrature.GetQuadrature(new int[] { D + 1 }, // vector components: ( avg_vel[0], ... , avg_vel[D-1], cell_volume )
                        this.LsTrk.GridDat,
                        rule,
                        delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            EvalResult.Clear();
                            for(int d = 0; d < D; d++)
                                U0_Spc[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                            var Vol = EvalResult.ExtractSubArrayShallow(-1, -1, D);
                            Vol.SetAll(1.0);
                        },
                        delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
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
                    for(int d = 0; d < D; d++) {
                        U0mean_check[d].ProjectField(1.0, U0.ElementAt(d).Evaluate,
                            new CellQuadratureScheme(false, Uncut).AddFixedOrderRules(LsTrk.GridDat, U0.ElementAt(d).Basis.Degree + 1));
                    }

                    foreach(var _Spc in this.LsTrk.SpeciesIdS) { // loop over species...
                        for(int d = 0; d < D; d++) {
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
