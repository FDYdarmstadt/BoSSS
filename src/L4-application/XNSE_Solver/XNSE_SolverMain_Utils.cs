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
        /// creates utilitiy fields
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


            CurvatureAlgorithms.CurvatureDriver(
                SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                CurvatureAlgorithms.FilterConfiguration.NoFilter,   // fromC0
                this.Curvature, out this.LevSetGradient, this.LsTrk,
                this.m_HMForder, this.DGLevSet.Current);


            CurvatureAlgorithms.CurvatureDriver(
                SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                new CurvatureAlgorithms.FilterConfiguration() {
                    gradOpt = CurvatureAlgorithms.GradientOption.LevSet,
                    hessOpt = CurvatureAlgorithms.HessianOption.LevSetGrad,
                    useFiltLevSetGrad = false,
                    useFiltLevSetHess = false,
                    FilterCurvatureCycles = 0,
                    LevelSetSource = CurvatureAlgorithms.LevelSetSource.fromDG,
                    PatchRecoveryDomWidth = 0,
                    NoOfPatchRecoverySweeps = 0,
                    CurvatureLimiting = false
                },
                this.DGCurvature, out this.DGLevSetGradient, this.LsTrk,
                this.m_HMForder, this.DGLevSet.Current);



            // ====================================
            // L2 error against exact solution
            // ====================================

            this.ComputeL2Error(phystime + dt);

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

            //if (Log_FourierLS != null) {
            //    // save restart infos for FLS
            //    Guid vecSamplP_id = this.DatabaseDriver.SaveVector<double>(Fourier_LevSet.getRestartInfo());
            //    if (base.MPIRank == 0) {
            //        Log_FourierLS.WriteLine(vecSamplP_id);
            //        Log_FourierLS.Flush();
            //    }
            //    // Log_files for FLS
            //    //if (this.Control.FourierLevSetControl.WriteFLSdata) {
            //    //    Fourier_LevSet.saveToLogFiles(TimestepNo.MajorNumber, phystime + dt);
            //    //}
            //    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            //}


            // ====================================================================== 
            // IO for further external postprocessing/ Query handling for Testprogram
            // ======================================================================

            if (this.Control.TestMode == true) {
                LogQueryValue(phystime + dt);
            } else {
                if (Log != null && this.Control.LogValues != XNSE_Control.LoggingValues.None && base.MPIRank == 0 && (TimestepNo.MajorNumber % this.Control.LogPeriod == 0))
                    try {
                        WriteLogLine(TimestepNo, phystime + dt);
                    } catch (Exception e) {
                        Console.WriteLine("An error occured during WriteLogLine: '{0}'", e);
                    }
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }

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


        public double[] ComputeSphericalPorperties() {

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), this.m_HMForder, 1).XQuadSchemeHelper;

            // area/volume
            double volume = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        volume += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            // surface
            double surface = 0.0;
            //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            var surfElemVol = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                surfElemVol.Compile(LsTrk.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return new double[] { volume, surface };

        }


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


        public double[] ComputeBenchmarkQuantities_RisingBubble() {

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // area of bubble
            double area = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        area += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            // center of mass/geometric center (for incommpressible fluid)
            int D = this.Grid.SpatialDimension;
            MultidimensionalArray center = MultidimensionalArray.Create(1, D);
            CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    NodeSet nodes_global = QR.Nodes.CloneAs();
                    for (int i = i0; i < i0 + Length; i++) {
                        LsTrk.GridDat.TransformLocal2Global(QR.Nodes, nodes_global, i);
                        EvalResult.AccSubArray(1.0, nodes_global, new int[] { i - i0, -1, -1 });
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            center[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();

            center.Scale(1.0 / area);

            // rise velocity
            MultidimensionalArray VelocityAtCenter = MultidimensionalArray.Create(1, D);

            // integral computation
            CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    for (int d = 0; d < D; d++) {
                        this.CurrentVel[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            VelocityAtCenter[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();
            VelocityAtCenter.Scale(1.0 / area);

            double v_rise = VelocityAtCenter[0, 1];

            //Console.WriteLine("rise velocity = " + v_rise);


            // circularity
            double diamtr_c = Math.Sqrt(4 * area / Math.PI);
            double perimtr_b = 0.0;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        perimtr_b += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            double circ = Math.PI * diamtr_c / perimtr_b;

            return new double[] { area, center[0, 0], center[0, 1], circ, VelocityAtCenter[0, 0], VelocityAtCenter[0, 1] };
        }


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


        public Tuple<List<double[]>, List<double[]>, List<double>> ComputeContactLineQuantities() {

            List<double[]> contactPoints = new List<double[]>();
            List<double[]> contactVelocities = new List<double[]>();
            List<double> contactAngles = new List<double>();

            ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(this.CurrentVel);

            var Phi = (LevelSet)LsTrk.LevelSets[0];
            var LevelSetGradient = new VectorField<SinglePhaseField>(Grid.SpatialDimension, Phi.Basis, SinglePhaseField.Factory);
            LevelSetGradient.Gradient(1.0, (SinglePhaseField)LsTrk.LevelSets[0]);
            SinglePhaseField[] Normals = LevelSetGradient.ToArray();

            XQuadSchemeHelper SchemeHelper = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadSchemeHelper;
            EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"));

            var QuadDom = SurfaceElement_Edge.Domain;
            var boundaryEdge = ((GridData)this.GridData).GetBoundaryEdgeMask().GetBitMask();
            var boundaryCutEdge = QuadDom.Intersect(new EdgeMask((GridData)this.GridData, boundaryEdge, MaskType.Geometrical));

            var factory = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, LsTrk.GridDat.Grid.RefElements[0]);
            SurfaceElement_Edge = new EdgeQuadratureScheme(factory, boundaryCutEdge);

            EdgeQuadrature.GetQuadrature(new int[] { 5 }, LsTrk.GridDat,
                SurfaceElement_Edge.Compile(LsTrk.GridDat, 0),
                delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {

                                // contact point
                                NodeSet Enode_l = QR.Nodes;
                    int trf = LsTrk.GridDat.Edges.Edge2CellTrafoIndex[i0, 0];
                    NodeSet Vnode_l = Enode_l.GetVolumeNodeSet(LsTrk.GridDat, trf);
                    NodeSet Vnode_g = Vnode_l.CloneAs();
                    int cell = LsTrk.GridDat.Edges.CellIndices[i0, 0];
                    LsTrk.GridDat.TransformLocal2Global(Vnode_l, Vnode_g, cell);
                                //Console.WriteLine("contact point: ({0},{1})", Vnode_g[0, 0], Vnode_g[0, 1]);

                                int D = Grid.SpatialDimension;
                    for (int d = 0; d < D; d++) {
                        EvalResult[0, 0, d] = Vnode_g[0, d];
                    }

                                // contact line velocity
                                MultidimensionalArray U_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    MultidimensionalArray U_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    for (int d = 0; d < D; d++) {
                        (meanVelocity[d] as SinglePhaseField).EvaluateEdge(i0, length, QR.Nodes, U_IN.ExtractSubArrayShallow(-1, -1, d), U_OUT.ExtractSubArrayShallow(-1, -1, d));
                    }

                    for (int d = 0; d < D; d++) {
                        EvalResult[0, 0, 2 + d] = U_IN[0, 0, d];
                    }

                                // contact angle
                                MultidimensionalArray normal_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    MultidimensionalArray normal_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    for (int d = 0; d < D; d++) {
                        Normals[d].EvaluateEdge(i0, length, QR.Nodes, normal_IN.ExtractSubArrayShallow(-1, -1, d), normal_OUT.ExtractSubArrayShallow(-1, -1, d));
                    }

                    double theta_surf = Math.Atan2(normal_IN[0, 0, 1], normal_IN[0, 0, 0]);
                    double theta_edge = Math.Atan2(LsTrk.GridDat.Edges.NormalsForAffine[i0, 1], LsTrk.GridDat.Edges.NormalsForAffine[i0, 0]);
                    double theta = (theta_surf - theta_edge) * (180 / Math.PI);

                    EvalResult[0, 0, 2 * D] = (theta > 180) ? theta - 180 : theta;
                                //Console.WriteLine("contact angle = {0}", EvalResult[0, 0, 2]);

                            },
                delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                    int D = Grid.SpatialDimension;
                    for (int i = 0; i < length; i++) {
                        if (ResultsOfIntegration[i, 2 * D] != 0.0) {
                            contactAngles.Add(Math.Abs(ResultsOfIntegration[i, 2 * D]));
                            double[] cp = new double[D];
                            double[] cpV = new double[D];
                            for (int d = 0; d < D; d++) {
                                cp[d] = ResultsOfIntegration[i, d];
                                cpV[d] = ResultsOfIntegration[i, 2 + d];
                            }
                            contactPoints.Add(cp);
                            contactVelocities.Add(cpV);
                        }
                    }
                }
            ).Execute();

            return new Tuple<List<double[]>, List<double[]>, List<double>>(contactPoints, contactVelocities, contactAngles);

        }


        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionVelocity"/> and <see cref="XNSE_Control.ExactSolutionPressure"/>).
        /// </summary>
        internal double[] ComputeL2Error(double time) {
            int D = this.GridData.SpatialDimension;
            double[] Ret = new double[D + 1];

            if (this.Control.ExactSolutionVelocity == null && this.Control.ExactSolutionPressure == null)
                // nothing to do
                return Ret;

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;


            // Velocity error
            // ==============
            if (this.Control.ExactSolutionVelocity != null) {
                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[D];

                foreach (var spc in this.LsTrk.SpeciesNames) {
                    L2Error_Species.Add(spc, new double[D]);

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                    for (int d = 0; d < D; d++) {
                        ConventionalDGField Vel_d = ((XDGField)this.CurrentVel[d]).GetSpeciesShadowField(spc);

                        L2Error_Species[spc][d] = Vel_d.L2Error(this.Control.ExactSolutionVelocity[spc][d].Vectorize(time), order, scheme);
                        L2Error[d] += L2Error_Species[spc][d].Pow2();

                        base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);
                    }
                }
                L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

                for (int d = 0; d < D; d++) {
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                    Ret[d] = L2Error[d];
                }

            }


            // pressure error
            // ==============
            if (this.Control.ExactSolutionPressure != null) {

                // pass 1: mean value of pressure difference
                double DiffInt = 0;
                foreach (var spc in this.LsTrk.SpeciesNames) {

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.GridData, order);

                    DiffInt += this.Pressure.GetSpeciesShadowField(spc).LxError(this.Control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b), rule);
                    //Volume +=  this.Pressure.GetSpeciesShadowField(spc).LxError(null, (a, b) => (1.0), rule);
                }
                double Volume2 = (new SubGrid(CellMask.GetFullMask(this.GridData))).Volume;
                double PressureDiffMean = DiffInt / Volume2;


                double L2Error = 0;
                Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                foreach (var spc in this.LsTrk.SpeciesNames) {

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.GridData, order);

                    double IdV = this.Pressure.GetSpeciesShadowField(spc).LxError(this.Control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b - PressureDiffMean).Pow2(), rule);
                    L2Error += IdV;
                    L2Error_Species.Add(spc, IdV.Sqrt());

                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure + "#" + spc, L2Error_Species[spc], true);
                }


                L2Error = L2Error.Sqrt();
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
                Ret[D] = L2Error;

            } //*/

            return Ret;
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


        /// <summary>
        /// testcase specific LogFile
        /// </summary>
        TextWriter Log;

        /// <summary>
        /// saves interface points
        /// </summary>
        TextWriter LogInterfaceP;


        /// <summary>
        /// initializes the format of the Log File
        /// </summary>
        /// <param name="sessionID"></param>
        public void InitLogFile(Guid sessionID) {

            string header;

            if (this.Control.WriteInterfaceP) {
                LogInterfaceP = base.DatabaseDriver.FsDriver.GetNewLog("InterfaceP", sessionID);
                header = String.Format("{0}\t{1}\t{2}", "#timestep", "#time", "interfacePoints");
                LogInterfaceP.WriteLine(header);
                LogInterfaceP.Flush();
            }

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


                        // Log file for the interface height
                        Log = base.DatabaseDriver.FsDriver.GetNewLog("Amplitude", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#timestep", "time", "magnitude", "real", "imaginary");

                        break;
                    }
                case XNSE_Control.LoggingValues.Dropletlike: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("SemiAxis", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", "#timestep", "time", "semi axis x", "semi axis y", "area", "perimeter");

                        break;
                    }
                case XNSE_Control.LoggingValues.RisingBubble: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("BenchmarkQuantities_RisingBubble", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "time", "area", "center of mass - x", "center of mass - y", "circularity", "rise velocity");

                        break;
                    }
                case XNSE_Control.LoggingValues.MovingContactLine: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("ContactAngle", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-angle");

                        break;
                    }
                case XNSE_Control.LoggingValues.DropletOnWall: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("DropletOnWall", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", "#timestep", "time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", 
                            "contact-angle", "droplet-height");

                        break;
                    }
                case XNSE_Control.LoggingValues.CapillaryHeight: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("CapillaryHeight", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}", "#timestep", "time", "capillary-height", "at-PositionX");

                        break;
                    }
                case XNSE_Control.LoggingValues.EvaporationL:
                case XNSE_Control.LoggingValues.EvaporationC: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("Evaporation", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "'#timestep", "time", "interfacePosition", "meanInterfaceVelocity", "meanMassFlux", "Temperatureprofile");

                        break;
                    }
                case XNSE_Control.LoggingValues.ChannelFlow: {

                        Log = base.DatabaseDriver.FsDriver.GetNewLog("ChannelFlow", sessionID);
                        header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#timestep", "time", "U_max", "derivedKinEnergy_max", "KinEnergy_max");

                        break;
                    }
                case XNSE_Control.LoggingValues.LinelikeLS:
                case XNSE_Control.LoggingValues.CirclelikeLS: {
                        throw new ArgumentException("specified LogFormat only valid for queries");
                    }
                default:
                    throw new ArgumentException("No specified LogFormat");
            }

            Log.WriteLine(header);
            Log.Flush();
        }


        List<double[]> contactPointsRef;

        
        /// <summary>
        /// writes one line to the Log File
        /// </summary>
        public void WriteLogLine(TimestepNumber TimestepNo, double phystime) {

            if (this.Control.WriteInterfaceP) {
                double[] interfaceP;
                if (Fourier_LevSet != null) {
                    interfaceP = Fourier_LevSet.current_interfaceP.To1DArray();
                } else {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                    interfaceP = interP.ResizeShallow(interP.Length).To1DArray();
                }
                string logline = String.Format("{0}\t{1}", TimestepNo, phystime);
                //for (int ip = 0; ip < interfaceP.Length; ip++) {
                //    logline = logline + "\t" + interfaceP[ip].ToString();
                //}
                logline = logline + "\t" + String.Join("\t", interfaceP.Select(ip => ip.ToString()).ToArray());
                LogInterfaceP.WriteLine(logline);
                LogInterfaceP.Flush();
            }

            switch (this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        Complex DFT_k;
                        int numP;
                        if (Fourier_LevSet != null) {
                            //amplitude = 2.0 * (Fourier_LevSet.DFT_coeff[1].Magnitude / Fourier_LevSet.current_samplP.Length);
                            DFT_k = Fourier_LevSet.DFT_coeff[(int)this.Control.AdditionalParameters[0]];
                            numP = Fourier_LevSet.current_samplP.Length;
                        } else {
                            MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                            Complex[] DFT_coeff = DiscreteFourierTransformation.TransformForward_nonequidistant(interP, this.Control.AdditionalParameters[1]);
                            DFT_k = DFT_coeff[(int)this.Control.AdditionalParameters[0]];
                            numP = interP.Lengths[0];
                            //amplitude = -2.0 * DFT_coeff[1].Imaginary / (double)interP.Lengths[0];
                            //amplitude = DiscreteFourierTransformation.SingleSidedPowerSpectrum(DFT_coeff)[1];
                        }
                        string logline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, 2.0 * DFT_k.Magnitude / numP, 2.0 * DFT_k.Real / numP, -2.0 * DFT_k.Imaginary / numP);
                        Log.WriteLine(logline);
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.Dropletlike: {

                        MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                        int numP = interP.Lengths[0];
                        double[] Xcoord = new double[numP];
                        double[] Ycoord = new double[numP];
                        for (int i = 0; i < numP; i++) {
                            Xcoord[i] = interP[i, 0];
                            Ycoord[i] = interP[i, 1];
                        }
                        double semiAxisX = Xcoord.Max() - Xcoord.Min();
                        double semiAxisY = Ycoord.Max() - Ycoord.Min();

                        double[] sphereProps = this.ComputeSphericalPorperties();

                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", TimestepNo, phystime, semiAxisX, semiAxisY, sphereProps[0], sphereProps[1]);
                        Log.WriteLine(line);
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.RisingBubble: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime, BmQ_RB[0], BmQ_RB[1], BmQ_RB[2], BmQ_RB[3], BmQ_RB[5]);
                        Log.WriteLine(line);
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.MovingContactLine: {

                        // contact angles at contact points
                        //=================================

                        List<double[]> contactPoints = new List<double[]>();
                        List<double[]> contactVelocities = new List<double[]>();
                        List<double> contactAngles = new List<double>();

                        var returnValues = this.ComputeContactLineQuantities();
                        contactPoints = returnValues.Item1;
                        contactVelocities = returnValues.Item2;
                        contactAngles = returnValues.Item3;


                        List<double[]> contactPointsSorted = new List<double[]>();
                        List<double[]> contactVelocitiesSorted = new List<double[]>();
                        List<double> contactAnglesSorted = new List<double>();

                        if (contactPointsRef == null) {
                            contactPointsRef = contactPoints;
                            contactPointsSorted = contactPoints;
                            contactVelocitiesSorted = contactVelocities; ;
                            contactAnglesSorted = contactAngles;
                        } else {
                            // sort
                            double eps = 1e-7;

                            do {
                                contactPointsSorted.Clear();
                                contactVelocitiesSorted.Clear();
                                contactAnglesSorted.Clear();

                                eps *= 10;
                                //Console.WriteLine("sorting contact line points");
                                foreach (var ptR in contactPointsRef) {
                                    //Console.WriteLine("ref point: ({0}, {1})", ptR[0], ptR[1]);
                                    for (int i = 0; i < contactAngles.Count(); i++) {
                                        double[] pt = contactPoints.ElementAt(i);
                                        //Console.WriteLine("sorting point: ({0}, {1})", pt[0], pt[1]);
                                        double xDiff = Math.Abs(pt[0] - ptR[0]);
                                        //Console.WriteLine("x diff = {0}", xDiff);
                                        double yDiff = Math.Abs(pt[1] - ptR[1]);
                                        //Console.WriteLine("y diff = {0}", yDiff);
                                        if (xDiff < eps && yDiff < eps) {
                                            //Console.WriteLine("sorted");
                                            contactPointsSorted.Add(pt.CloneAs());
                                            contactVelocitiesSorted.Add(contactVelocities.ElementAt(i).CloneAs());
                                            contactAnglesSorted.Add(contactAngles.ElementAt(i));
                                            break;
                                        }
                                    }
                                }

                            } while (contactPointsSorted.Count != contactPointsRef.Count && eps < 1.0);

                            contactPointsRef = contactPointsSorted;
                        }

                        MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                        double[] height = interP.ExtractSubArrayShallow(new int[] { -1, 1 }).To1DArray();
                        //Console.WriteLine("Droplet height: {0}", height.Max());

                        for (int p = 0; p < contactAnglesSorted.Count; p++) {
                            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", TimestepNo, phystime,
                                contactPointsSorted.ElementAt(p)[0], contactPointsSorted.ElementAt(p)[1],
                                contactVelocitiesSorted.ElementAt(p)[0], contactVelocitiesSorted.ElementAt(p)[1], contactAnglesSorted.ElementAt(p), height.Max());
                            Log.WriteLine(line);
                        }
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.DropletOnWall: {

                        // contact angles at contact points
                        //=================================

                        List<double[]> contactPoints = new List<double[]>();
                        List<double[]> contactVelocities = new List<double[]>();
                        List<double> contactAngles = new List<double>();

                        var returnValues = this.ComputeContactLineQuantities();
                        contactPoints = returnValues.Item1;
                        contactVelocities = returnValues.Item2;
                        contactAngles = returnValues.Item3;


                        List<double[]> contactPointsSorted = new List<double[]>();
                        List<double[]> contactVelocitiesSorted = new List<double[]>();
                        List<double> contactAnglesSorted = new List<double>();

                        if (contactPointsRef == null) {
                            contactPointsRef = contactPoints;
                            contactPointsSorted = contactPoints;
                            contactVelocitiesSorted = contactVelocities; ;
                            contactAnglesSorted = contactAngles;
                        } else {
                            // sort
                            double eps = 1e-7;

                            do {
                                contactPointsSorted.Clear();
                                contactVelocitiesSorted.Clear();
                                contactAnglesSorted.Clear();

                                eps *= 10;
                                //Console.WriteLine("sorting contact line points");
                                foreach (var ptR in contactPointsRef) {
                                    //Console.WriteLine("ref point: ({0}, {1})", ptR[0], ptR[1]);
                                    for (int i = 0; i < contactAngles.Count(); i++) {
                                        double[] pt = contactPoints.ElementAt(i);
                                        //Console.WriteLine("sorting point: ({0}, {1})", pt[0], pt[1]);
                                        double xDiff = Math.Abs(pt[0] - ptR[0]);
                                        //Console.WriteLine("x diff = {0}", xDiff);
                                        double yDiff = Math.Abs(pt[1] - ptR[1]);
                                        //Console.WriteLine("y diff = {0}", yDiff);
                                        if (xDiff < eps && yDiff < eps) {
                                            //Console.WriteLine("sorted");
                                            contactPointsSorted.Add(pt.CloneAs());
                                            contactVelocitiesSorted.Add(contactVelocities.ElementAt(i).CloneAs());
                                            contactAnglesSorted.Add(contactAngles.ElementAt(i));
                                            break;
                                        }
                                    }
                                }

                            } while (contactPointsSorted.Count != contactPointsRef.Count && eps < 1.0);

                            contactPointsRef = contactPointsSorted;
                        }


                        for (int p = 0; p < contactAnglesSorted.Count; p++) {
                            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime,
                                contactPointsSorted.ElementAt(p)[0], contactPointsSorted.ElementAt(p)[1],
                                contactVelocitiesSorted.ElementAt(p)[0], contactVelocitiesSorted.ElementAt(p)[1], contactAnglesSorted.ElementAt(p));
                            Log.WriteLine(line);
                        }
                        Log.Flush();

                        return;
                    }
                case XNSE_Control.LoggingValues.CapillaryHeight: {

                        MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);

                        double h_min = double.MaxValue, x_pos = 0.0;
                        for (int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                            if (InterfacePoints[i, 1] < h_min) {
                                h_min = InterfacePoints[i, 1];
                                x_pos = InterfacePoints[i, 0];
                            }
                        }

                        string line = String.Format("{0}\t{1}\t{2}\t{3}", TimestepNo, phystime, h_min, x_pos);
                        Log.WriteLine(line);
                        Log.Flush();

                        break;
                    }
                case XNSE_Control.LoggingValues.EvaporationL:
                case XNSE_Control.LoggingValues.EvaporationC: {

                        MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                        double nNodes = InterfacePoints.GetLength(0);

                        double posI = 0.0;
                        if (this.Control.LogValues == XNSE_Control.LoggingValues.EvaporationL)
                            posI = InterfacePoints.ExtractSubArrayShallow(-1, 1).To1DArray().Sum() / nNodes;

                        double hVap = this.Control.ThermalParameters.hVap;
                        double MassFlux = EvapVelocMean * hVap;

                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, posI, EvapVelocMean, MassFlux);

                        // temperature profile
                        int N = 100;
                        double[] tempP = new double[N + 1];

                        double[] probe = new double[N + 1];
                        //if (this.Control.LogValues == XNSE_Control.LoggingValues.EvaporationL) {
                        double L = this.Control.AdditionalParameters[0];
                        double x_probe = this.Control.AdditionalParameters[1];
                        for (int i = 0; i <= N; i++) {
                            probe = new double[] { x_probe, i * (L / (double)N) };
                            try {
                                tempP[i] = this.Temperature.ProbeAt(probe);
                            } catch {
                                tempP[i] = 0.0;
                            }
                        }
                        //}

                        line = line + "\t" + String.Join("\t", tempP.Select(ip => ip.ToString()).ToArray());

                        Log.WriteLine(line);
                        Log.Flush();

                        break;
                    }
                case XNSE_Control.LoggingValues.ChannelFlow: {

                        double C_halfw = 0.5;

                        double Umax = this.XDGvelocity.Velocity[0].ProbeAt(new double[]{ 1.01, C_halfw });
                        double derivedKinEmax = 0.5 * Umax.Pow2();
                        double KinEmax = this.KineticEnergy.ProbeAt(new double[] { 1.01, C_halfw });

                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, Umax, derivedKinEmax, KinEmax);

                        Log.WriteLine(line);
                        Log.Flush();

                        break;
                    }
                default:
                    throw new ArgumentException("No specified LogFormat");
            }

        }

        double EvapVelocMean;

        /// <summary>
        /// encapsulated handling of query values
        /// </summary>
        public void LogQueryValue(double phystime) {

            base.QueryResultTable.LogValue("time", phystime);

            if (this.Control.WriteInterfaceP) {
                double[] interfaceP;
                if (Fourier_LevSet != null) {
                    interfaceP = Fourier_LevSet.current_interfaceP.To1DArray();
                } else {
                    MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                    interfaceP = interP.ResizeShallow(interP.Length).To1DArray();
                }

                base.QueryResultTable.LogValue("interfaceP", interfaceP);

            }

            switch (this.Control.LogValues) {
                case XNSE_Control.LoggingValues.Wavelike: {

                        double amplitude;
                        if (Fourier_LevSet != null) {
                            amplitude = 2.0 * (Fourier_LevSet.DFT_coeff[1].Magnitude / Fourier_LevSet.current_samplP.Length);
                        } else {
                            MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                            Complex[] DFT_coeff = DiscreteFourierTransformation.TransformForward_nonequidistant(interP, this.Control.AdditionalParameters[1]);
                            amplitude = -2.0 * DFT_coeff[1].Imaginary / (double)interP.Lengths[0];
                            //amplitude = DiscreteFourierTransformation.SingleSidedPowerSpectrum(DFT_coeff)[1];
                        }

                        base.QueryResultTable.LogValue("amplitude", amplitude);
                        return;
                    }
                case XNSE_Control.LoggingValues.RisingBubble: {

                        double[] BmQ_RB = this.ComputeBenchmarkQuantities_RisingBubble();

                        base.QueryResultTable.LogValue("area", BmQ_RB[0]);
                        base.QueryResultTable.LogValue("yCM", BmQ_RB[2]);
                        base.QueryResultTable.LogValue("circ", BmQ_RB[3]);
                        base.QueryResultTable.LogValue("riseV", BmQ_RB[5]);

                        return;
                    }
                case XNSE_Control.LoggingValues.LinelikeLS: {
                        break;
                    }
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
            for (int iRow = Mtx.RowPartitioning.i0; iRow < Mtx.RowPartitioning.iE; iRow++) {
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
                    this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                    AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);

                BlockMsrMatrix SaddlePointMatrix = new BlockMsrMatrix(this.SaddlePointProblemMapping);
                double[] AffineDummy = new double[this.SaddlePointProblemMapping.LocalLength];



                DelComputeOperatorMatrix(SaddlePointMatrix, AffineDummy, this.SaddlePointProblemMapping,
                    this.CurrentSolution.Mapping.Fields.ToArray(), CurrentAgglomeration.CellLengthScales, 0.0);

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
                        this.MultigridOperatorConfig);

                    // extract
                    ////////////

                    MsrMatrix FullMatrix = mgOp.OperatorMatrix.ToMsrMatrix();

                    MsrMatrix DiffMatrix;
                    {
                        int[] VelVarIdx = D.ForLoop(d => d);

                        int[] USubMatrixIdx_Row = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int[] USubMatrixIdx_Col = mgOp.Mapping.GetSubvectorIndices(VelVarIdx);
                        int L = USubMatrixIdx_Row.Length;

                        DiffMatrix = new MsrMatrix(L, L, 1, 1);
                        FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(int[]), USubMatrixIdx_Col, default(int[]));
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
