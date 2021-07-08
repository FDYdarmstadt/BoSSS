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

        //======================================
        // partial file for energy related code
        //======================================


        //=====================================
        // Field declaration and instantiation
        //=====================================
        #region fields

#pragma warning disable 649

        XDGField[] prevVel;

        /// <summary>
        /// kinetic energy derived via \f$ \rho \frac{vec{u} \cdot \vec{u}}{ 2 } \f$
        /// </summary>
        XDGField DerivedKineticEnergy;

        //XDGField DerivedKineticEnergyChangerate;

        XDGField GeneratedKineticEnergy;

        //XDGField KineticEnergyNSE;

        //XDGField KineticEnergyNSEchangerate;

        /// <summary>
        /// kinetic energy computed via <see cref="KineticEnergyBalanceOperator"/>
        /// </summary>
        internal XDGField KineticEnergy {
            get;
            private set;
        }

        /// <summary>
        /// Residual of the kinetic energy balance
        /// </summary>
        XDGField ResidualKineticEnergy;


        //XDGField prevKineticEnergyChangerate;

        //XDGField KineticEnergyChangerate;

        //XDGField prevKineticEnergy;
        //XDGField pprevKineticEnergy;


        /// <summary>
        /// source term for the kinetic energy
        /// </summary>
        XDGField KineticDissipation;

        //XDGField PowerOfStresses;


#pragma warning restore 649

        /// <summary>
        /// creates energy related fields
        /// </summary>
        public void CreateEnergyFields() {

            int D = this.GridData.SpatialDimension;

            IOListOption register = (this.Control.RegisterUtilitiesToIOFields) ? IOListOption.Always : IOListOption.ControlFileDetermined;

            if (this.Control.solveKineticEnergyEquation) {

                this.KineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), VariableNames.KineticEnergy);
                base.RegisterField(this.KineticEnergy);
                this.ResidualKineticEnergy = new XDGField(this.KineticEnergy.Basis, "ResidualKineticEnergy");
                base.RegisterField(this.ResidualKineticEnergy);

                //this.prevKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions["KineticEnergy"].Degree)));

                this.GeneratedKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "GeneratedKineticEnergy");
                base.RegisterField(this.GeneratedKineticEnergy, register);

                //this.KineticEnergyChangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticEnergyChangerate");
                //base.RegisterField(this.KineticEnergyChangerate, register);

                //this.prevKineticEnergyChangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)));

            }

            if (this.Control.ComputeEnergyProperties) {

                if (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
                    prevVel = new XDGField[D];
                    for (int d = 0; d < D; d++) {
                        prevVel[d] = new XDGField(this.XDGvelocity.Velocity[d].Basis);
                    }
                }

                //this.prevKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "previousKineticEnergy");
                //base.RegisterField(this.prevKineticEnergy);
                //this.pprevKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "ppreviousKineticEnergy");
                //base.RegisterField(this.pprevKineticEnergy);

                //this.KineticEnergyNSE = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticEnergyNSE");
                //base.RegisterField(this.KineticEnergyNSE);

                //this.KineticEnergyNSEchangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticEnergyNSEchangerate");
                //base.RegisterField(this.KineticEnergyNSEchangerate);

                this.DerivedKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "DerivedKineticEnergy");
                base.RegisterField(this.DerivedKineticEnergy, register);

                //this.DerivedKineticEnergyChangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "DerivedKineticEnergyChangerate");
                //base.RegisterField(this.DerivedKineticEnergyChangerate, register);

                this.KineticDissipation = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticDissipation");
                base.RegisterField(this.KineticDissipation, register);

                //this.PowerOfStresses = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "PowerOfStresses");
                //base.RegisterField(this.PowerOfStresses, register);

            }

        }


        #endregion



        // =========================================================
        // related stuff for property tracking (e.g. kinetic energy)
        // =========================================================
        #region tracking

        ResidualLogger m_EnergyLogger;

        /// <summary>
        /// Logger for kinetic and surface energy.
        /// </summary>
        ResidualLogger EnergyLogger {
            get {
                if (!this.Control.ComputeEnergyProperties)
                    return null;

                if (m_EnergyLogger == null) {
                    m_EnergyLogger = new ResidualLogger(base.MPIRank, base.DatabaseDriver, base.CurrentSessionInfo.ID);
                    m_EnergyLogger.WriteResidualsToConsole = false;
                    m_EnergyLogger.WriteResidualsToTextFile = true;
                    m_EnergyLogger.TextFileFileName = "Energy";
                }

                return m_EnergyLogger;
            }
        }


        /// <summary>
        /// spatial Operator for the kinetic energy balance
        /// </summary>
        XSpatialOperatorMk2 KineticEnergyOperator;

        int m_HMForder_kinE;

        public void generateKinEnergyOperator(XNSFE_OperatorConfiguration config) {

            int degK = this.KineticEnergy.Basis.Degree;

            int D = this.GridData.SpatialDimension;

            string[] CodName = new string[] { EquationNames.KineticEnergyEquation };
            string[] Params = ArrayTools.Cat(
                VariableNames.VelocityVector(D),
                //(new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0,D),
                //VariableNames.NormalVector(D),
                //VariableNames.Curvature,
                //VariableNames.SurfaceForceVector(D),
                VariableNames.VelocityX_GradientVector(),
                VariableNames.VelocityY_GradientVector(),
                //new string[] { "VelocityXGradX_GradientX", "VelocityXGradX_GradientY" },
                //new string[] { "VelocityXGradY_GradientX", "VelocityXGradY_GradientY" },
                //new string[] { "VelocityYGradX_GradientX", "VelocityYGradX_GradientY" },
                //new string[] { "VelocityYGradY_GradientX", "VelocityYGradY_GradientY" },
                VariableNames.Pressure,
                VariableNames.PressureGradient(D),
                VariableNames.GravityVector(D)
                );

            string[] DomName = new string[] { VariableNames.KineticEnergy };


            // create operator
            // ===============

            this.m_HMForder_kinE = degK * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);

            KineticEnergyOperator = new XSpatialOperatorMk2(DomName, Params, CodName, (A, B, C) => this.m_HMForder_kinE,  this.LsTrk.SpeciesNames);
            KineticEnergyOperator.AgglomerationThreshold = this.Control.AgglomerationThreshold;

            // build the operator
            // ==================
            {

                // species bulk components
                for (int spc = 0; spc < LsTrk.TotalNoOfSpecies; spc++) {
                    Solution.EnergyCommon.XOperatorComponentsFactory.AddSpeciesKineticEnergyEquation(KineticEnergyOperator, config, D, LsTrk.SpeciesNames[spc], LsTrk.SpeciesIdS[spc], BcMap, LsTrk);
                }

                // interface components
                //Solution.EnergyCommon.XOperatorComponentsFactory.AddInterfaceKineticEnergyEquation(KineticEnergyOperator, config, D, BcMap, LsTrk, degK);
                //CurvatureRequired = true;


                // finalize
                // ========

                KineticEnergyOperator.Commit();

            }

        }


        /// <summary>
        /// timestepper for solving the kinetic energy equation postprocessing
        /// </summary>
        XdgBDFTimestepping KineticEnergyTimestepper;

        public void CreateTimestepperForKineticEnergySolve() {

            KineticEnergyTimestepper = new XdgBDFTimestepping(
                    this.KineticEnergy.Mapping.Fields,
                    this.ResidualKineticEnergy.Mapping.Fields,
                    KineticEnergyOperator.InvokeParameterFactory(this.KineticEnergy.Mapping.Fields),
                    LsTrk,
                    true,
                    DelComputeEnergyOperatorMatrix, KineticEnergyOperator, DelUpdateLevelSet_EnergyOperator,
                    (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? bdfOrder + 0 : 1,
                    (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? LevelSetHandling.Coupled_Once : this.Control.Timestepper_LevelSetHandling,
                    this.XOpConfig.mmsd,
                    SpatialOperatorType.LinearTimeDependent, //(this.Control.PhysicalParameters.IncludeConvection) ? SpatialOperatorType.Nonlinear : SpatialOperatorType.LinearTimeDependent,
                    this.MultigridOperatorConfig_kinE, base.MultigridSequence,
                    this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder_kinE,
                    this.Control.AgglomerationThreshold,
                    true,
                    this.Control.NonLinearSolver,
                    this.Control.LinearSolver
                    );
            KineticEnergyTimestepper.m_ResLogger = m_EnergyResLogger;
            KineticEnergyTimestepper.m_ResidualNames = this.ResidualKineticEnergy.Mapping.Fields.Select(f => f.Identification).ToArray();
            KineticEnergyTimestepper.Timestepper_Init = (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) ? this.Control.Timestepper_BDFinit : TimeStepperInit.SingleInit;
            KineticEnergyTimestepper.incrementTimesteps = this.Control.incrementTimesteps;
            KineticEnergyTimestepper.PushLevelSet = delegate() { }; // dummy push does nothing
            KineticEnergyTimestepper.IterUnderrelax = this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? this.Control.LSunderrelax : 1.0;

            KineticEnergyTimestepper.Config_LevelSetConvergenceCriterion = this.Control.LevelSet_ConvergenceCriterion;

            // solver 
            this.Control.NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.NonLinearSolver.MinSolverIterations; 
            //m_BDF_Timestepper.config_NonLinearSolver.MinSolverIterations = (this.Control.Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) ? 1 : this.Control.Solver_MinIterations;


        }


        //public void EvaluateKineticEnergy(double physTime, double[] fluxEval) {


        //    // parameter assembly
        //    // ==================
        //    #region param assembly

        //    int D = this.Grid.SpatialDimension;
        //    //LevelSet Phi = (LevelSet)(this.LsTrk.LevelSets[0]);


        //    //// linearization velocity:
        //    //DGField[] U0_U0mean;
        //    //if (this.U0meanrequired) {
        //    //    XDGBasis U0meanBasis = new XDGBasis(this.LsTrk, 0);
        //    //    VectorField<XDGField> U0mean = new VectorField<XDGField>(D, U0meanBasis, "U0mean_", XDGField.Factory);
        //    //    U0mean.Clear();
        //    //    if (this.physParams.IncludeConvection)
        //    //        ComputeAverageU(U0, U0mean, CutCellQuadOrder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, CutCellQuadOrder, 1).XQuadSchemeHelper);

        //    //    U0_U0mean = ArrayTools.Cat<DGField>(U0, U0mean);
        //    //} else {
        //    //    U0_U0mean = new DGField[2 * D];
        //    //}


        //    //// normals:
        //    //SinglePhaseField[] Normals; // Normal vectors: length not normalized - will be normalized at each quad node within the flux functions.
        //    //if (this.NormalsRequired) {
        //    //    if (LevelSetGradient == null) {
        //    //        LevelSetGradient = new VectorField<SinglePhaseField>(D, Phi.Basis, SinglePhaseField.Factory);
        //    //        LevelSetGradient.Gradient(1.0, Phi);
        //    //    }
        //    //    Normals = LevelSetGradient.ToArray();
        //    //} else {
        //    //    Normals = new SinglePhaseField[D];
        //    //}

        //    //// curvature:
        //    //SinglePhaseField Curvature;
        //    //if (this.CurvatureRequired) {
        //    //    Curvature = ExternalyProvidedCurvature;
        //    //} else {
        //    //    Curvature = null;
        //    //}


        //    // velocity gradient vectors
        //    var VelMap = new CoordinateMapping(this.CurrentSolution.Fields.GetSubVector(0, D));
        //    DGField[] VelParam = VelMap.Fields.ToArray();

        //    VectorField<DGField> GradVelX = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradient", XDGField.Factory);
        //    for (int d = 0; d < D; d++) {
        //        foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //            DGField f_Spc = ((VelParam[0] as XDGField).GetSpeciesShadowField(Spc));
        //            SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
        //            (GradVelX[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
        //        }
        //    }
        //    GradVelX.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    VectorField<DGField> GradVelY = new VectorField<DGField>(D, VelParam[1].Basis, "VelocityYGradient", XDGField.Factory);
        //    for (int d = 0; d < D; d++) {
        //        foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //            DGField f_Spc = ((VelParam[1] as XDGField).GetSpeciesShadowField(Spc));
        //            SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
        //            (GradVelY[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
        //        }
        //    }
        //    GradVelY.ForEach(F => F.CheckForNanOrInf(true, true, true));


        //    //VectorField<DGField> GradXVelXGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradX_Gradient", XDGField.Factory);
        //    //for (int d = 0; d < D; d++) {
        //    //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //    //        DGField f_Spc = ((GradVelX[0] as XDGField).GetSpeciesShadowField(Spc));
        //    //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
        //    //        (GradXVelXGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
        //    //    }
        //    //}
        //    //GradXVelXGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    //VectorField<DGField> GradYVelXGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradY_Gradient", XDGField.Factory);
        //    //for (int d = 0; d < D; d++) {
        //    //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //    //        DGField f_Spc = ((GradVelX[1] as XDGField).GetSpeciesShadowField(Spc));
        //    //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
        //    //        (GradYVelXGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
        //    //    }
        //    //}
        //    //GradYVelXGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    //VectorField<DGField> GradXVelYGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradX_Gradient", XDGField.Factory);
        //    //for (int d = 0; d < D; d++) {
        //    //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //    //        DGField f_Spc = ((GradVelY[0] as XDGField).GetSpeciesShadowField(Spc));
        //    //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
        //    //        (GradXVelYGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
        //    //    }
        //    //}
        //    //GradXVelYGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    //VectorField<DGField> GradYVelYGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradY_Gradient", XDGField.Factory);
        //    //for (int d = 0; d < D; d++) {
        //    //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //    //        DGField f_Spc = ((GradVelY[1] as XDGField).GetSpeciesShadowField(Spc));
        //    //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
        //    //        (GradYVelYGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
        //    //    }
        //    //}
        //    //GradYVelYGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));


        //    //Tecplot.PlotFields(ArrayTools.Cat<DGField>(GradVelX, GradVelY, GradXVelXGrad, GradXVelYGrad, GradYVelXGrad, GradYVelYGrad), 
        //    //    "GradientParam", physTime, 3);


        //    // pressure and gradient
        //    var PressMap = new CoordinateMapping(this.CurrentSolution.Fields.ToArray()[D]);
        //    DGField[] PressParam = PressMap.Fields.ToArray();

        //    VectorField<DGField> PressGrad = new VectorField<DGField>(D, PressParam[0].Basis, "PressureGrad", XDGField.Factory);
        //    for (int d = 0; d < D; d++) {
        //        foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //            DGField f_Spc = ((PressParam[0] as XDGField).GetSpeciesShadowField(Spc));
        //            SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
        //            (PressGrad[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
        //        }
        //    }
        //    PressGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    // gravity
        //    var GravMap = new CoordinateMapping(this.XDGvelocity.Gravity.ToArray());
        //    DGField[] GravParam = GravMap.Fields.ToArray();


        //    #endregion


        //    // concatenate everything
        //    var Params = ArrayTools.Cat<DGField>(
        //        this.CurrentSolution.Fields.GetSubVector(0, D),
        //        GradVelX,
        //        GradVelY,
        //        //GradXVelXGrad,
        //        //GradYVelXGrad,
        //        //GradXVelYGrad,
        //        //GradYVelYGrad,
        //        PressParam,
        //        PressGrad,
        //        GravMap);



        //    XSpatialOperatorMk2.XEvaluatorNonlin eval = this.KineticEnergyOperator.GetEvaluatorEx(this.LsTrk,
        //        this.KineticEnergy.ToEnumerable().ToArray(), Params, this.KineticEnergy.Mapping,
        //        this.LsTrk.SpeciesIdS.ToArray());

        //    var agg = this.LsTrk.GetAgglomerator(this.LsTrk.SpeciesIdS.ToArray(),
        //        this.KineticEnergy.Basis.Degree * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2),
        //        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold);

        //    foreach (var kv in agg.CellLengthScales) {
        //        eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
        //    }

        //    //if (this.KineticEnergyOperator.SurfaceElementOperator.TotalNoOfComponents > 0) {
        //    //    foreach (var kv in InterfaceLengths) {
        //    //        eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
        //    //        //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("lambda_interface", lambdaI);
        //    //        //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("mu_interface", muI);
        //    //    }
        //    //}

        //    eval.time = physTime;

        //    eval.Evaluate(1.0, 1.0, fluxEval);

        //}


        /// <summary>
        /// 
        /// </summary>
        IDictionary<SpeciesId, IEnumerable<double>> KinEnergyMassScale {
            get {
                double rho_A = this.Control.PhysicalParameters.rho_A,
                    rho_B = this.Control.PhysicalParameters.rho_B;

                int D = this.GridData.SpatialDimension;


                double[] scale_A = new double[1];
                double[] scale_B = new double[1];

                scale_A.SetAll(rho_A); // mass matrix in momentum equation (kinetic energy equation)
                scale_B.SetAll(rho_B); // mass matrix in momentum equation (kinetic energy equation)

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), scale_A);
                R.Add(this.LsTrk.GetSpeciesId("B"), scale_B);


                return R;
            }
        }


        void DelComputeEnergyOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime, int iLsTrkHist) {
            if(iLsTrkHist != 1)
                throw new NotSupportedException();
            
            // parameter assembly
            // ==================
            #region param assembly

            int D = this.Grid.SpatialDimension;
            //LevelSet Phi = (LevelSet)(this.LsTrk.LevelSets[0]);


            //// linearization velocity:
            //DGField[] U0_U0mean;
            //if (this.U0meanrequired) {
            //    XDGBasis U0meanBasis = new XDGBasis(this.LsTrk, 0);
            //    VectorField<XDGField> U0mean = new VectorField<XDGField>(D, U0meanBasis, "U0mean_", XDGField.Factory);
            //    U0mean.Clear();
            //    if (this.physParams.IncludeConvection)
            //        ComputeAverageU(U0, U0mean, CutCellQuadOrder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, CutCellQuadOrder, 1).XQuadSchemeHelper);

            //    U0_U0mean = ArrayTools.Cat<DGField>(U0, U0mean);
            //} else {
            //    U0_U0mean = new DGField[2 * D];
            //}


            //// normals:
            //SinglePhaseField[] Normals; // Normal vectors: length not normalized - will be normalized at each quad node within the flux functions.
            //if (this.NormalsRequired) {
            //    if (LevelSetGradient == null) {
            //        LevelSetGradient = new VectorField<SinglePhaseField>(D, Phi.Basis, SinglePhaseField.Factory);
            //        LevelSetGradient.Gradient(1.0, Phi);
            //    }
            //    Normals = LevelSetGradient.ToArray();
            //} else {
            //    Normals = new SinglePhaseField[D];
            //}

            //// curvature:
            //SinglePhaseField Curvature;
            //if (this.CurvatureRequired) {
            //    Curvature = ExternalyProvidedCurvature;
            //} else {
            //    Curvature = null;
            //}


            // velocity gradient vectors
            var VelMap = new CoordinateMapping(this.CurrentSolution.Fields.GetSubVector(0, D));
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

            VectorField<DGField> GradVelY = new VectorField<DGField>(D, VelParam[1].Basis, "VelocityYGradient", XDGField.Factory);
            for (int d = 0; d < D; d++) {
                foreach (var Spc in this.LsTrk.SpeciesIdS) {
                    DGField f_Spc = ((VelParam[1] as XDGField).GetSpeciesShadowField(Spc));
                    SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
                    (GradVelY[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                }
            }
            GradVelY.ForEach(F => F.CheckForNanOrInf(true, true, true));


            //VectorField<DGField> GradXVelXGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradX_Gradient", XDGField.Factory);
            //for (int d = 0; d < D; d++) {
            //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
            //        DGField f_Spc = ((GradVelX[0] as XDGField).GetSpeciesShadowField(Spc));
            //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
            //        (GradXVelXGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
            //    }
            //}
            //GradXVelXGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

            //VectorField<DGField> GradYVelXGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradY_Gradient", XDGField.Factory);
            //for (int d = 0; d < D; d++) {
            //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
            //        DGField f_Spc = ((GradVelX[1] as XDGField).GetSpeciesShadowField(Spc));
            //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
            //        (GradYVelXGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
            //    }
            //}
            //GradYVelXGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

            //VectorField<DGField> GradXVelYGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradX_Gradient", XDGField.Factory);
            //for (int d = 0; d < D; d++) {
            //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
            //        DGField f_Spc = ((GradVelY[0] as XDGField).GetSpeciesShadowField(Spc));
            //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
            //        (GradXVelYGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
            //    }
            //}
            //GradXVelYGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

            //VectorField<DGField> GradYVelYGrad = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradY_Gradient", XDGField.Factory);
            //for (int d = 0; d < D; d++) {
            //    foreach (var Spc in this.LsTrk.SpeciesIdS) {
            //        DGField f_Spc = ((GradVelY[1] as XDGField).GetSpeciesShadowField(Spc));
            //        SubGrid sf = this.LsTrk.Regions.GetSpeciesSubGrid(Spc);
            //        (GradYVelYGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d); //, optionalSubGrid: sf);
            //    }
            //}
            //GradYVelYGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));


            //Tecplot.PlotFields(ArrayTools.Cat<DGField>(GradVelX, GradVelY, GradXVelXGrad, GradXVelYGrad, GradYVelXGrad, GradYVelYGrad), 
            //    "GradientParam", physTime, 3);


            // pressure and gradient
            var PressMap = new CoordinateMapping(this.CurrentSolution.Fields.ToArray()[D]);
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
            var GravMap = new CoordinateMapping(this.XDGvelocity.Gravity.ToArray());
            DGField[] GravParam = GravMap.Fields.ToArray();


            #endregion


            // concatenate everything
            var Params = ArrayTools.Cat<DGField>(
                this.CurrentSolution.Fields.GetSubVector(0, D),
                GradVelX,
                GradVelY,
                //GradXVelXGrad,
                //GradYVelXGrad,
                //GradXVelYGrad,
                //GradYVelYGrad,
                PressParam,
                PressGrad,
                GravMap);



            // assemble the matrix & affine vector
            // ===================================



            // compute matrix
            if (OpMtx != null) {

                var mtxBuilder = KineticEnergyOperator.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping);

                mtxBuilder.time = phystime;

                foreach (var kv in AgglomeratedCellLengthScales) {
                    mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                }

                //if (this.KineticEnergyOperator.SurfaceElementOperator.TotalNoOfComponents > 0) {
                //    foreach (var kv in InterfaceLengths) {
                //        mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                //        //mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("lambda_interface", lambdaI);
                //        //mtxBuilder.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("mu_interface", muI);
                //    }
                //}

                mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

            } else {
                XSpatialOperatorMk2.XEvaluatorNonlin eval = this.KineticEnergyOperator.GetEvaluatorEx(this.LsTrk,
                    this.KineticEnergy.ToEnumerable().ToArray(), Params, this.KineticEnergy.Mapping);

                foreach (var kv in AgglomeratedCellLengthScales) {
                    eval.CellLengthScales[kv.Key] = kv.Value;
                }

                //if (this.KineticEnergyOperator.SurfaceElementOperator.TotalNoOfComponents > 0) {
                //    foreach (var kv in InterfaceLengths) {
                //        eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
                //        //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("lambda_interface", lambdaI);
                //        //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("mu_interface", muI);
                //    }
                //}

                eval.time = phystime;

                eval.Evaluate(1.0, 1.0, OpAffine);

            }

        }

        /// <summary>
        /// dummy delegate for coupled operators
        /// </summary>
        /// <param name="CurrentState"></param>
        /// <param name="Phystime"></param>
        /// <param name="dt"></param>
        /// <param name="underrelax"></param>
        /// <param name="incremental"></param>
        /// <returns></returns>
        double DelUpdateLevelSet_EnergyOperator(DGField[] CurrentState, double Phystime, double dt, double underrelax, bool incremental) {
            // do nothing
            return 0.0;
        }


        /// <summary>
        /// The residual logger for this application.
        /// </summary>
        public ResidualLogger EnergyResLogger {
            get {
                return m_EnergyResLogger;
            }
        }

        ResidualLogger m_EnergyResLogger;


        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig_kinE {
            get {
                int pKinE = this.KineticEnergy.Basis.Degree;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {

                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[1];

                    // configuration for kinetic energy
                    configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Math.Max(1, pKinE - iLevel) },
                        mode = this.Control.KineticEnergyeBlockPrecondMode,
                        VarIndex = new int[] { 0 }
                    };

                }


                return configs;
            }
        }


        #endregion



    }

}
