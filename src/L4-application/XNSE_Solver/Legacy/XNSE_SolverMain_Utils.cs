﻿/* =======================================================================
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

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.Tracing;
using ilPSP.LinSolvers;

using BoSSS.Platform;

using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;

using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XdgTimestepping;
using NUnit.Framework;
using MPI.Wrappers;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XNSE_Solver.Legacy {

    /// <summary>
    /// Solver for Incompressible Multiphase flows; 
    /// </summary>
    public partial class XNSE_SolverMain : BoSSS.Solution.Application<XNSE_Control> {

        //=======================================
        // partial file for various utility code
        //=======================================


        //=====================================
        // Field declaration and instantiation
        //=====================================
        #region fields

#pragma warning disable 649

        /// <summary>
        /// Divergence of velocity ->
        /// Conservation of mass for incompressibility
        /// </summary>
        XDGField divVelocity;


        SinglePhaseField MassBalanceAtInterface;

        VectorField<SinglePhaseField> MomentumBalanceAtInterface;

        VectorField<SinglePhaseField> SurfaceTensionForce;

        SinglePhaseField EnergyBalanceAtInterface;

        SinglePhaseField InterfaceDivergence;


#pragma warning restore 649

        /// <summary>
        /// creates utility fields
        /// </summary>
        public void CreateUtilityFields() {

            int D = this.GridData.SpatialDimension;

            IOListOption register = (this.Control.RegisterUtilitiesToIOFields) ? IOListOption.Always : IOListOption.ControlFileDetermined;

            XDGBasis b = new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree);
            this.divVelocity = new XDGField(b, "DivergenceVelocity");
            base.RegisterField(this.divVelocity, register);


            if (this.Control.CheckJumpConditions) {

                Basis basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.VelocityX].Degree + (this.Control.FieldOptions[VariableNames.VelocityX].Degree - 1));
                this.MassBalanceAtInterface = new SinglePhaseField(basis, "MassBalanceAtInterface");
                base.RegisterField(this.MassBalanceAtInterface, register);

                basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree + (this.Control.FieldOptions[VariableNames.VelocityX].Degree - 1));
                this.MomentumBalanceAtInterface = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(basis, d + "-MomentumBalanceAtInterface")));
                base.RegisterField(this.MomentumBalanceAtInterface, register);

                basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree + (this.Control.FieldOptions[VariableNames.VelocityX].Degree - 1));
                this.SurfaceTensionForce = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(basis, d + "-SurfaceTensionForce")));
                base.RegisterField(this.SurfaceTensionForce, register);

                //basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree + this.Control.FieldOptions[VariableNames.VelocityX].Degree + (this.Control.FieldOptions[VariableNames.VelocityX].Degree - 1));
                //this.EnergyBalanceAtInterface = new SinglePhaseField(basis, "EnergyBalanceAtInterface");
                //base.RegisterField(this.EnergyBalanceAtInterface, register);

                basis = new Basis(this.GridData, this.Control.FieldOptions[VariableNames.Pressure].Degree * 2);
                this.InterfaceDivergence = new SinglePhaseField(basis, "InterfaceDivergence");
                base.RegisterField(this.InterfaceDivergence, register);
            }
        }


        #endregion


/*
        //=====================
        // pre-/postprocessing
        //=====================

        /// <summary>
        /// 
        /// </summary>
        /// <param name="TimestepInt"></param>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <param name="TimestepNo"></param>
        private void Preprocessing(int TimestepInt, double phystime, double dt, TimestepNumber TimestepNo) {


            //if (this.Control.CheckInterfaceProps) {
            //    double CL_length = this.GetContactLineLength();
            //    Console.WriteLine("contact line length = {0}", CL_length);

            //    double[] props = this.ComputeSphericalPorperties();
            //    Console.WriteLine("volume = {0}", props[0]);
            //    Console.WriteLine("surface = {0}", props[1]);
            //}


            //if (this.Control.solveKineticEnergyEquation) {
            //    double[] rhoS = new double[] { this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B };
            //    EnergyUtils.ProjectKineticEnergy(this.KineticEnergy, this.LsTrk, this.XDGvelocity.Velocity.ToArray(), rhoS, this.m_HMForder);
            //}

        }

    */
        /// <summary>
        /// 
        /// </summary>
        /// <param name="TimestepInt"></param>
        /// <param name="phystime"></param>
        /// <param name="dt"></param>
        /// <param name="TimestepNo"></param>
        private void Postprocessing(int TimestepInt, double phystime, double dt, TimestepNumber TimestepNo) {

            // =================
            // plot shadowfields
            // =================


            //DGField[] HeatFluxParam = new DGField[this.GridData.SpatialDimension];
            //HeatFluxParam = new VectorField<XDGField>(this.GridData.SpatialDimension, this.Temperature.Basis, "HeatFlux0_", XDGField.Factory).ToArray();
            //Dictionary<string, double> kSpc = new Dictionary<string, double>();
            //kSpc.Add("A", -this.Control.ThermalParameters.k_A);
            //kSpc.Add("B", -this.Control.ThermalParameters.k_B);
            //XNSEUtils.ComputeGradientForParam(this.Temperature, HeatFluxParam, this.LsTrk, kSpc);

            //DGField TempA = this.Temperature.GetSpeciesShadowField("A");
            //DGField TempB = this.Temperature.GetSpeciesShadowField("B");
            //DGField HeatFluxYA = (HeatFluxParam[1] as XDGField).GetSpeciesShadowField("A");
            //DGField HeatFluxYB = (HeatFluxParam[1] as XDGField).GetSpeciesShadowField("B");

            //Tecplot.PlotFields(new DGField[] { TempA, TempB, HeatFluxYA, HeatFluxYB, this.DGLevSet.Current, this.LevSet}, "HeatFluxParam" + hack_TimestepIndex, hack_Phystime, 4);


            // ======================================
            // Check jump conditions at the interface 
            // ======================================

            #region check jump conditions

            if (this.Control.CheckJumpConditions) {

                this.MassBalanceAtInterface.Clear();
                for (int d = 0; d < this.Grid.SpatialDimension; d++) {
                    this.MomentumBalanceAtInterface[d].Clear();
                    this.SurfaceTensionForce[d].Clear();
                }
                //this.EnergyBalanceAtInterface.Clear();
                this.InterfaceDivergence.Clear();


                CurvatureAlgorithms.CurvatureDriver(
                    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                    CurvatureAlgorithms.FilterConfiguration.NoFilter,
                    this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                    this.m_HMForder, this.DGLevSet.Current);

                ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(this.XDGvelocity.Velocity.ToArray());


                // mass balance
                XNSEUtils.ProjectMassBalanceAtInterface(this.MassBalanceAtInterface, 1.0, this.XDGvelocity.Velocity);

                //double velJump_Norm = XNSEUtils.VelocityJumpNorm(this.XDGvelocity.Velocity, false, this.m_HMForder);
                //Console.WriteLine("velocity jump norm = {0}", velJump_Norm);


                // momentum balance
                XNSEUtils.ProjectMomentumBalanceAtInterface(this.MomentumBalanceAtInterface, 1.0, this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                    this.Control.PhysicalParameters);
                //XNSEUtils.ProjectMomentumBalanceNorm(this.MomentumBalanceAtInterface, this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                //   this.Control.PhysicalParameters, true);

                //double[] momBal_Norm = XNSEUtils.MomentumBalanceNormAtInterface(this.Pressure, this.XDGvelocity.Velocity, this.Curvature,
                //    this.Control.PhysicalParameters, this.Control.AdvancedDiscretizationOptions.SurfStressTensor, this.m_HMForder);
                //Console.WriteLine("x-momentum balance norm = {0}", momBal_Norm[0]);
                //Console.WriteLine("y-momentum balance norm = {0}", momBal_Norm[1]);


                // surface tension force
                XNSEUtils.ProjectSurfaceTensionForce(this.SurfaceTensionForce, 1.0, meanVelocity, this.Curvature, this.LsTrk, 
                    this.Control.PhysicalParameters, this.Control.AdvancedDiscretizationOptions.SurfStressTensor);


                // energy balance
                //EnergyUtils.ProjectEnergyBalanceAtInterface(this.EnergyBalanceAtInterface, 1.0, this.Pressure, this.XDGvelocity.Velocity, meanVelocity, this.Curvature,
                //    this.Control.PhysicalParameters);

                //double energyBal_Norm = XNSEUtils.EnergyBalanceNormAtInterface(this.Pressure, this.XDGvelocity.Velocity, meanVelocity, this.Curvature,
                //    this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, this.Control.PhysicalParameters.Sigma, this.m_HMForder);
                //Console.WriteLine("energy balance norm = {0}", energyBal_Norm);



                // interface divergence
                EnergyUtils.ProjectInterfaceDivergence(this.InterfaceDivergence, 1.0, meanVelocity, this.LsTrk, this.Control.PhysicalParameters);



            }

            #endregion


            // ====================================================
            // Compute integral energy of the system
            // ====================================================

            #region integral energy computation 

            if (this.Control.ComputeEnergyProperties) {

                // compute current energies
                double[] rhoS = new double[] { this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B };
                double currentKinEnergy = EnergyUtils.GetKineticEnergy(this.LsTrk, this.XDGvelocity.Velocity.ToArray(), rhoS, this.m_HMForder);
                double currentSurfEnergy = EnergyUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma, this.m_HMForder);
                EnergyUtils.ProjectKineticEnergy(this.DerivedKineticEnergy, this.LsTrk, this.XDGvelocity.Velocity.ToArray(), rhoS, this.m_HMForder);

                // compute changerates (kinetic, surface)
                double CR_KinEnergy = 0.0;
                double CR_SurfEnergy = 0.0;
                if (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
                    double prevKinEnergy = EnergyUtils.GetKineticEnergy(this.LsTrk, this.prevVel, rhoS, this.m_HMForder, 0);
                    CR_KinEnergy = (currentKinEnergy - prevKinEnergy) / dt;

                    double prevSurfEnergy = EnergyUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma, this.m_HMForder, 0);
                    CR_SurfEnergy = (currentSurfEnergy - prevSurfEnergy) / dt;

                    //Console.WriteLine("current kinetic energy = {0}; actual changerate = {1}", currentKinEnergy, CR_KinEnergy);
                    //Console.WriteLine("current surface energy = {0}; actual changerate = {1}", currentSurfEnergy, CR_SurfEnergy);
                }

                // changerate of kinetic energy from discretization
                double[] muS = new double[] { this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B };
                DGField[] extForce = new DGField[this.LsTrk.GridDat.SpatialDimension];
                for (int i = 0; i < this.LsTrk.GridDat.SpatialDimension; i++) {
                    extForce[i] = this.XDGvelocity.Gravity[i].CloneAs();
                    for (int iSpc = 0; iSpc < LsTrk.SpeciesIdS.Count; iSpc++) {
                        SpeciesId spcId = LsTrk.SpeciesIdS[iSpc];
                        switch (LsTrk.GetSpeciesName(spcId)) {
                            case "A": {
                                    (extForce[i] as XDGField).GetSpeciesShadowField("A").Scale(this.Control.PhysicalParameters.rho_A);
                                    break;
                                }
                            case "B": {
                                    (extForce[i] as XDGField).GetSpeciesShadowField("B").Scale(this.Control.PhysicalParameters.rho_B);
                                    break;
                                }
                        }
                    }
                }
                double kineticDissipationBulk = EnergyUtils.GetKineticDissipation(this.LsTrk, this.XDGvelocity.Velocity.ToArray(), muS, this.m_HMForder);
                //EnergyUtils.ProjectKineticDissipation(this.KineticDissipation, this.LsTrk, this.XDGvelocity.Velocity.ToArray(), muS, this.m_HMForder);
                //EnergyUtils.ProjectPowerOfStresses(this.PowerOfStresses, this.LsTrk, this.Pressure, this.XDGvelocity.Velocity.ToArray(), muS, this.m_HMForder);

                // changerate of surface energy form discretization
                ConventionalDGField[] meanVelocity = XNSEUtils.GetMeanVelocity(this.XDGvelocity.Velocity, this.LsTrk,
                    this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B);
                double SurfDivergence = EnergyUtils.GetSurfaceChangerate(this.LsTrk, meanVelocity, this.m_HMForder);

                //EnergyUtils.ProjectEnergyBalanceNorm(this.EnergyJumpCondition, 1.0, this.Pressure, this.XDGvelocity.Velocity, meanVelocity, this.Curvature, muS[0], muS[1], this.Control.PhysicalParameters.Sigma, this.m_HMForder);


                // logging
                this.EnergyLogger.TimeStep = TimestepInt;
                this.EnergyLogger.CustomValue(phystime + dt, "PhysicalTime");
                this.EnergyLogger.CustomValue(currentKinEnergy, "KineticEnergy");
                this.EnergyLogger.CustomValue(currentSurfEnergy, "SurfaceEnergy");
                this.EnergyLogger.CustomValue(CR_KinEnergy, "ChangerateKineticEnergy");
                this.EnergyLogger.CustomValue(CR_SurfEnergy, "ChangerateSurfaceEnergy");
                this.EnergyLogger.CustomValue(SurfDivergence, "SurfaceDivergence");
                this.EnergyLogger.CustomValue(kineticDissipationBulk, "KineticDissipationBulk");


                // surface viscosity parts
                if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor != SurfaceSressTensor.Isotropic) {

                    double shearViscEnergyCR = 0.0;
                    double dilViscEnergyCR = 0.0;

                    // surface shear viscosity energy
                    if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation
                        || this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        shearViscEnergyCR = EnergyUtils.GetInterfaceShearViscosityEnergyCR(this.LsTrk, meanVelocity, this.Control.PhysicalParameters.mu_I, this.m_HMForder);
                    }

                    // surface dilatational viscosity energy
                    if (this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation
                        || this.Control.AdvancedDiscretizationOptions.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        dilViscEnergyCR = EnergyUtils.GetInterfaceDilatationalViscosityEnergyCR(this.LsTrk, meanVelocity, this.Control.PhysicalParameters.lambda_I, this.m_HMForder);
                    }


                    this.EnergyLogger.CustomValue(shearViscEnergyCR, "ShearViscosityDR");
                    this.EnergyLogger.CustomValue(dilViscEnergyCR, "DilatationalViscosityDR");

                    Console.WriteLine("current kinetic energy dissipation from discretization = {0}", kineticDissipationBulk + shearViscEnergyCR + dilViscEnergyCR);

                } else {

                    Console.WriteLine("current kinetic energy dissipation from discretization = {0}", kineticDissipationBulk);
                }


                // logging
                // =======

                this.EnergyLogger.NextTimestep(true);


                //double[] RhoS = new double[] { this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B };
                //double newKinEnergy = XNSEUtils.GetKineticEnergy(this.LsTrk, this.CurrentVel, RhoS, this.m_HMForder);
                //double oldKinEnergy;
                //if (base.Control.CompMode == AppControl._CompMode.Transient) {
                //    DGField[] prevVel;
                //    if (this.XDGvelocity != null)
                //        prevVel = this.XDGvelocity.Velocity.ToArray();//<DGField>();
                //    else
                //        prevVel = this.DGvelocity.Velocity.ToArray();
                //    oldKinEnergy = XNSEUtils.GetKineticEnergy(this.LsTrk, prevVel, RhoS, this.m_HMForder);
                //} else if (base.Control.CompMode == AppControl._CompMode.Steady) {
                //    oldKinEnergy = newKinEnergy;
                //} else {
                //    throw new NotSupportedException();
                //}
                //double surfEnergy = XNSEUtils.GetSurfaceEnergy(this.LsTrk, this.Control.PhysicalParameters.Sigma.Abs(), this.m_HMForder);


                //// Logging and Console Output
                //// ===========================

                //this.EnergyLogger.TimeStep = TimestepInt;
                //this.EnergyLogger.CustomValue(phystime + dt, "PhysicalTime");
                //this.EnergyLogger.CustomValue(oldKinEnergy, "OldKineticEnergy");
                //this.EnergyLogger.CustomValue(newKinEnergy, "NewKineticEnergy");
                //this.EnergyLogger.CustomValue(surfEnergy, "SurfaceEnergy");

                //this.EnergyLogger.NextTimestep(true);

            }



            #endregion


            // ====================================
            // divergence of velocity and curvature
            // ====================================

            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            this.divVelocity.Clear();
            this.divVelocity.Divergence(1.0, this.XDGvelocity.Velocity);
            //ilPSP.Environment.StdoutOnlyOnRank0 = true;


            //CurvatureAlgorithms.CurvatureDriver(
            //    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
            //    CurvatureAlgorithms.FilterConfiguration.NoFilter,   // fromC0
            //    this.Curvature, out this.LevSetGradient, this.LsTrk,
            //    this.m_HMForder, this.DGLevSet.Current);


            //CurvatureAlgorithms.CurvatureDriver(
            //    SurfaceStressTensor_IsotropicMode.Curvature_Projected,
            //    new CurvatureAlgorithms.FilterConfiguration() {
            //        gradOpt = CurvatureAlgorithms.GradientOption.LevSet,
            //        hessOpt = CurvatureAlgorithms.HessianOption.LevSetGrad,
            //        useFiltLevSetGrad = false,
            //        useFiltLevSetHess = false,
            //        FilterCurvatureCycles = 0,
            //        LevelSetSource = CurvatureAlgorithms.LevelSetSource.fromDG,
            //        PatchRecoveryDomWidth = 0,
            //        NoOfPatchRecoverySweeps = 0,
            //        CurvatureLimiting = false
            //    },
            //    this.DGCurvature, out this.DGLevSetGradient, this.LsTrk,
            //    this.m_HMForder, this.DGLevSet.Current);


            // ====================================
            // L2 error against exact solution
            // ====================================

            //this.ComputeL2Error(phystime + dt);

            // =========== 
            // check area
            // ===========

            //SpeciesId spcId = LsTrk.SpeciesIdS[0];
            //double area = XNSEUtils.GetSpeciesArea(LsTrk, spcId, MomentFittingVariant);
            //Console.WriteLine("Area of species 'A' = {0}", area);


            //double[] props = this.ComputeSphericalPorperties();
            //Console.WriteLine("volume = {0}", props[0]);
            //Console.WriteLine("surface = {0}", props[1]);

            //double CL_length = this.GetContactLineLength();
            //Console.WriteLine("contact line length = {0}", CL_length);

            //double CapHeight = GetCapillaryHeight();
            //Console.WriteLine("Capillary height = {0}", CapHeight);

            //ContinuityEnforcer = new ContinuityProjection(
            //        ContBasis: this.LevSet.Basis,
            //        DGBasis: this.DGLevSet.Current.Basis,
            //        gridData: GridData,
            //        Option: Control.LSContiProjectionMethod
            //        );


            //var returnValues = this.ComputeContactLineQuantities();
            //List<double> contactangles = returnValues.Item3;
            //Console.WriteLine("contact angles: {0},{1}", contactangles.ToArray()[0], contactangles.ToArray()[1]);

            // ====================================
            // IO related to Fourier level set
            // ====================================

            if (Log_FourierLS != null) {
                // save restart infos for FLS
                Guid vecSamplP_id = this.DatabaseDriver.SaveVector<double>(Fourier_LevSet.getRestartInfo());
                if (base.MPIRank == 0) {
                    Log_FourierLS.WriteLine(vecSamplP_id);
                    Log_FourierLS.Flush();
                }
                // Log_files for FLS
                //if (this.Control.FourierLevSetControl.WriteFLSdata) {
                //    Fourier_LevSet.saveToLogFiles(TimestepNo.MajorNumber, phystime + dt);
                //}
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }


            // ====================================================================== 
            // IO for further external postprocessing/ Query handling for Testprogram
            // ======================================================================

            //if (this.Control.TestMode == true) {
            //    //LogQueryValue(phystime + dt);
            //} else {
            //    if (Log != null && this.Control.LogValues != XNSE_Control.LoggingValues.None && base.MPIRank == 0 && (TimestepNo.MajorNumber % this.Control.LogPeriod == 0))
            //        try {
            //            WriteLogLine(TimestepNo, phystime + dt);
            //        } catch (Exception e) {
            //            Console.WriteLine("An error occurred during WriteLogLine: '{0}'", e);
            //        }
            //    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            //}

            //Console.WriteLine("Pause");

            //=======================
            //var jmpNorm = XNSEUtils.VelocityJumpNorm(this.XDGvelocity.Velocity, true, MomentFittingVariant, -1);
            //Console.WriteLine("Velocity Jump Norm: " + jmpNorm);
            //var jmpStressNorm = XNSEUtils.MomentumJumpNorm(this.XDGvelocity.Velocity, this.Pressure, this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, MomentFittingVariant, -1);
            //Console.WriteLine("Stress Jump Norm [0]: " + jmpStressNorm[0]);
            //Console.WriteLine("Stress Jump Norm [1]: " + jmpStressNorm[1]);
            //PrintVelocityAtLevSet(TimestepNo.MajorNumber);

        }


        //==============================
        // various property computation
        //==============================
        #region property computation

        /*

        public double GetContactLineLength() {

            double CL_length = 0.0;

            if (this.LsTrk.GridDat.SpatialDimension == 3) {

                var metrics = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder);

                XQuadSchemeHelper SchemeHelper = metrics.XQuadSchemeHelper;
                EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"));

                var QuadDom = SurfaceElement_Edge.Domain;
                var boundaryCutEdge = QuadDom.Intersect(this.GridData.GetBoundaryEdgeMask());

                var innerDom = QuadDom.Except(this.GridData.GetBoundaryEdgeMask());

                System.Collections.BitArray lowerBits = new System.Collections.BitArray(((GridData)this.GridData).Edges.Count);
                foreach (Chunk cnk in boundaryCutEdge) {
                    for (int iE = cnk.i0; iE < cnk.JE; iE++) {
                        if (((GridData)this.GridData).Edges.EdgeTags[iE] == 1) {
                            lowerBits[iE] = true;
                        }
                    }
                }
                EdgeMask lowerDom = new EdgeMask(this.GridData, lowerBits);

                EdgeMask dom = lowerDom;

                var factory = metrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, LsTrk.GridDat.Grid.RefElements[0]);
                SurfaceElement_Edge = new EdgeQuadratureScheme(factory, dom);

                EdgeQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    SurfaceElement_Edge.Compile(LsTrk.GridDat, this.m_HMForder),
                    delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {
                        EvalResult.SetAll(1.0);
                    },
                    delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < length; i++)
                            CL_length += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
            }

            return CL_length;

        }

        */
        


        public double[] ComputeBenchmarkQuantities_LineInterface() {

            // interface length
            double length = 0.0;
            length = XNSEUtils.GetInterfaceLength(LsTrk);

            // species area
            double area = 0.0;
            area = XNSEUtils.GetSpeciesArea(LsTrk, LsTrk.SpeciesIdS[0]);

            // interface mean angle

            return new double[] { length, area };
        }


        


        #endregion


        //=========
        // logging
        //=========
        #region logging

        
        /// <summary>
        /// saves the vector Guid for the sample points 
        /// </summary>
        TextWriter Log_FourierLS;

        /*
        /// <summary>
        /// testcase specific LogFile
        /// </summary>
        TextWriter Log;

        /// <summary>
        /// saves interface points
        /// </summary>
        TextWriter LogInterfaceP;
        */

        /*
        /// <summary>
        /// initializes the format of the Log File
        /// </summary>
        /// <param name="sessionID"></param>
        public void InitLogFile(Guid sessionID) {

            string header;


            switch (this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        // File for physical data
                        TextWriter setUpData = base.DatabaseDriver.FsDriver.GetNewLog("SetUpData", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", "lambda", "H0", "rho1", "rho2", "mu1", "mu2", "sigma", "g");
                        setUpData.WriteLine(header);
                        setUpData.Flush();
                        string data = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", this.Control.AdditionalParameters[1], this.Control.AdditionalParameters[2], this.Control.PhysicalParameters.rho_A, this.Control.PhysicalParameters.rho_B,
                            this.Control.PhysicalParameters.mu_A, this.Control.PhysicalParameters.mu_B, this.Control.PhysicalParameters.Sigma, this.Control.AdditionalParameters[3]);
                        setUpData.WriteLine(data);
                        setUpData.Flush();


                        
                    }
                
            }

            
        }
        */


        
        /*
        /// <summary>
        /// writes one line to the Log File
        /// </summary>
        public void WriteLogLine(TimestepNumber TimestepNo, double phystime) {

            

            switch (this.Control.LogValues) {


                case XNSE_Control.LoggingValues.ChannelFlow: {

                        

                        break;
                    }
                default:
                    throw new ArgumentException("No specified LogFormat");
            }

        }

        double EvapVelocMean;

        /*
        /// <summary>
        /// encapsulated handling of query values
        /// </summary>
        public void LogQueryValue(double phystime) {

           

            switch (this.Control.LogValues) {
                
                


                case XNSE_Control.LoggingValues.CirclelikeLS: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        base.QueryResultTable.LogValue("area", BmQ_RB[0]);
                        base.QueryResultTable.LogValue("xM", BmQ_RB[1]);
                        base.QueryResultTable.LogValue("yM", BmQ_RB[2]);
                        base.QueryResultTable.LogValue("circ", BmQ_RB[3]);
                        base.QueryResultTable.LogValue("vM_x", BmQ_RB[4]);
                        base.QueryResultTable.LogValue("vM_y", BmQ_RB[5]);

                        break;
                    }
                default:
                    return;
            }

        }
        //*/

        #endregion


        //==================================
        // spatial operator matrix analysis
        //==================================
        #region analysis


        UnsetteledCoordinateMapping SaddlePointProblemMapping {
            get {
                return this.CurrentSolution.Mapping;
            }
        }

        int RepairZeroRows(MsrMatrix Mtx) {
            int NoOfZeroRows = 0;
            for (long iRow = Mtx.RowPartitioning.i0; iRow < Mtx.RowPartitioning.iE; iRow++) {
                if (Mtx.GetNoOfNonZerosPerRow(iRow) == 0) {
                    Mtx[iRow, iRow] = +1.0;
                    NoOfZeroRows++;
                }
            }
            return NoOfZeroRows;
        }

        /// <summary>
        /// Computes condition number, etc. of the current system matrix.
        /// </summary>
        /// <param name="CheckAssertions"></param>
        /// <param name="AnalysisLevel">
        /// - equal 0: check that pressure gradient and velocity divergence are transpose
        /// - equal 1: in addition, positive definiteness test.
        /// - equal 2: in addition, check condition number and eigenvalues using MATLAB
        /// </param>
        public void SpatialOperatorMatrixAnalysis(bool CheckAssertions, int AnalysisLevel) {
            using (new FuncTrace()) {
                int D = this.Grid.SpatialDimension;

                if (AnalysisLevel < 0 || AnalysisLevel > 2)
                    throw new ArgumentException();

                // ========================================================
                // compute agglomeration & operator (saddle-point) matrix
                // ========================================================

                //var CurrentAgglomeration = new MultiphaseCellAgglomerator(
                //        this.DelUpdateCutCellMetrics(),
                //        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                //        AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);
                var CurrentAgglomeration = this.LsTrk.GetAgglomerator(
                    this.LsTrk.SpeciesIdS.ToArray(), m_HMForder,
                    this.Control.AgglomerationThreshold,
                    AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);

                BlockMsrMatrix SaddlePointMatrix = new BlockMsrMatrix(this.SaddlePointProblemMapping);
                double[] AffineDummy = new double[this.SaddlePointProblemMapping.LocalLength];



                DelComputeOperatorMatrix(SaddlePointMatrix, AffineDummy, this.SaddlePointProblemMapping,
                    this.CurrentSolution.Mapping.Fields.ToArray(), CurrentAgglomeration.CellLengthScales, 0.0, 1);

                // =============================
                // AnalysisLevel 0
                // =============================
                {
                    var SaddlePointMatrixT = SaddlePointMatrix.Transpose();

                    CoordinateVector TestVec = new CoordinateVector(this.CurrentSolution.Mapping.Fields.Select(f => f.CloneAs()).ToArray());

                    double testsumPos = 0.0;
                    double testsumNeg = 0.0;
                    for (int rnd_seed = 0; rnd_seed < 20; rnd_seed++) {

                        // fill the pressure components of the test vector
                        TestVec.Clear();
                        Random rnd = new Random(rnd_seed);
                        XDGField Pressack = TestVec.Mapping.Fields[D] as XDGField;
                        int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                        for (int j = 0; j < J; j++) {
                            int N = Pressack.Basis.GetLength(j);

                            for (int n = 0; n < N; n++)
                                Pressack.Coordinates[j, n] = rnd.NextDouble();
                        }

                        //Pressack.Clear(H.Complement());
                        //Pressack.GetSpeciesShadowField("A").Clear();

                        // Gradient times P:
                        double[] R1 = new double[TestVec.Count];
                        SaddlePointMatrix.SpMV(1.0, TestVec, 0.0, R1);       // R1 = Grad * P
                                                                             //Console.WriteLine("L2 of 'Grad * P': " + R1.L2Norm());

                        // transpose of Divergence times P: 
                        double[] R2 = new double[TestVec.Count];
                        SaddlePointMatrixT.SpMV(1.0, TestVec, 0.0, R2);      // R2 = divT * P
                                                                             //Console.WriteLine("L2 of 'divT * P': " + R2.L2Norm());


                        TestVec.Clear();
                        TestVec.Acc(1.0, R1);
                        TestVec.Acc(1.0, R2);


                        // analyze!
                        testsumNeg += GenericBlas.L2Dist(R1, R2);

                        R2.ScaleV(-1.0);
                        testsumPos += GenericBlas.L2Dist(R1, R2);

                    }

                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (+): " + testsumPos);
                    Console.WriteLine("Pressure/Divergence Symmetry error in all tests (-): " + testsumNeg);

                    if (CheckAssertions)
                        Assert.LessOrEqual(Math.Abs(testsumNeg), testsumPos * 1.0e-13);
                }


                // =============================
                // AnalysisLevel 1 and 2
                // =============================

                if (AnalysisLevel > 0) {
                    AggregationGridBasis[][] MgBasis = AggregationGridBasis.CreateSequence(this.MultigridSequence, this.CurrentSolution.Mapping.BasisS);
                    //todo: AsyncCallback update
                    MgBasis.UpdateXdgAggregationBasis(CurrentAgglomeration);
                    MultigridOperator mgOp = new MultigridOperator(MgBasis, this.SaddlePointProblemMapping,
                        SaddlePointMatrix, this.MassFact.GetMassMatrix(this.SaddlePointProblemMapping, false),
                        this.MultigridOperatorConfig, null);

                    // extract
                    ////////////

                    MsrMatrix FullMatrix = mgOp.OperatorMatrix.ToMsrMatrix();

                    MsrMatrix DiffMatrix;
                    {
                        int[] VelVarIdx = D.ForLoop(d => d);

                        long[] USubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        long[] USubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int L = USubMatrixIdx_Row.Length;

                        DiffMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(long[]), USubMatrixIdx_Col, default(long[]));
                    }

                    int Zeros_FullMatrix = RepairZeroRows(FullMatrix);
                    int Zeros_DiffMatrix = RepairZeroRows(DiffMatrix);

                    Console.WriteLine("Indefinite Basis elements in Diffusion matrix:\t" + Zeros_DiffMatrix);
                    Console.WriteLine("Indefinite Basis elements in Saddle-point matrix:\t" + Zeros_FullMatrix);

                    base.QueryHandler.ValueQuery("NoOfIndef_DiffMtx", Zeros_DiffMatrix, false);
                    base.QueryHandler.ValueQuery("NoOfIndef_FullMtx", Zeros_FullMatrix, false);

                    // operator analysis
                    //////////////////////

                    bool posDef;
                    if (AnalysisLevel > 1) {
                        // +++++++++++++++++++++++++++++++
                        // check condition number, etc
                        // +++++++++++++++++++++++++++++++

                        MultidimensionalArray ret = MultidimensionalArray.Create(1, 5);
                        Console.WriteLine("Calling MATLAB/Octave...");
                        using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                            bmc.PutSparseMatrix(FullMatrix, "FullMatrix");
                            bmc.PutSparseMatrix(DiffMatrix, "DiffMatrix");
                            bmc.Cmd("DiffMatrix = 0.5*(DiffMatrix + DiffMatrix');");
                            bmc.Cmd("condNoDiffMatrix = condest(DiffMatrix);");
                            bmc.Cmd("condNoFullMatrix = condest(FullMatrix);");
                            bmc.Cmd("eigiMaxi = eigs(DiffMatrix,1,'lm')");
                            bmc.Cmd("eigiMini = eigs(DiffMatrix,1,'sm')");
                            bmc.Cmd("lasterr");
                            bmc.Cmd("[V,r]=chol(DiffMatrix);");
                            bmc.Cmd("ret = [condNoDiffMatrix, condNoFullMatrix, eigiMaxi, eigiMini, r]");
                            bmc.GetMatrix(ret, "ret");

                            bmc.Execute(false);
                        }

                        double condNoDiffMatrix = ret[0, 0];
                        double condNoFullMatrix = ret[0, 1];
                        double eigiMaxi = ret[0, 2];
                        double eigiMini = ret[0, 3];
                        posDef = ret[0, 4] == 0;

                        Console.WriteLine("Eigenvalue range of diffusion matrix: {0} to {1}", eigiMini, eigiMaxi);

                        Console.WriteLine("Condition number diffusion operator: {0:0.####E-00}", condNoDiffMatrix);
                        Console.WriteLine("Condition number full operator: {0:0.####E-00}", condNoFullMatrix);

                    } else {
                        // +++++++++++++++++++++++++++++++++++++++
                        // test only for positive definiteness
                        // +++++++++++++++++++++++++++++++++++++++

                        var DiffMatrixFull = DiffMatrix.ToFullMatrixOnProc0();

                        posDef = true;
                        try {
                            DiffMatrixFull.Cholesky();
                        } catch (ArithmeticException) {
                            posDef = false;
                        }
                    }

                    double DiffSymm = DiffMatrix.SymmetryDeviation();
                    Console.WriteLine("Symmetry deviation of diffusion matrix: " + DiffSymm);

                    if (posDef)
                        Console.WriteLine("Good news: Diffusion operator matrix seems to be positive definite.");
                    else
                        Console.WriteLine("WARNING: Diffusion operator matrix is not positive definite.");

                    if (CheckAssertions) {
                        if (Control.AdvancedDiscretizationOptions.ViscosityMode == ViscosityMode.FullySymmetric && Control.PhysicalParameters.IncludeConvection == false) {
                            double compVal = DiffMatrix.InfNorm() * 1e-13;
                            Assert.LessOrEqual(DiffSymm, compVal, "Diffusion matrix seems to be non-symmetric.");
                            Assert.IsTrue(posDef, "Positive definiteness test failed.");
                        }
                    }
                }
            }
        }


        /// <summary>
        /// makes direct use of <see cref="XdgTimesteppingBase.OperatorAnalysis"/>; aids the condition number scaling analysis
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis() {
            return this.m_BDF_Timestepper.OperatorAnalysis(plotStencilCondNumV: true);
        }


        #endregion

    }

}
