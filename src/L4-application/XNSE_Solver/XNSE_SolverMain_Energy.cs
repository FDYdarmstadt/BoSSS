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

        XDGField DerivedKineticEnergyChangerate;

        XDGField GeneratedKineticEnergy;

        XDGField KineticEnergyNSE;

        XDGField KineticEnergyNSEchangerate;

        /// <summary>
        /// kinetic energy computed via <see cref="KineticEnergyBalanceOperator"/>
        /// </summary>
        XDGField KineticEnergy;

        XDGField prevKineticEnergyChangerate;

        XDGField KineticEnergyChangerate;

        XDGField prevKineticEnergy;

        /// <summary>
        /// Residual of the kinetic energy balance
        /// </summary>
        XDGField ResidualKineticEnergy;

        /// <summary>
        /// source term for the kinetic energy
        /// </summary>
        XDGField KineticDissipation;

        XDGField PowerOfStresses;


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

                this.KineticEnergyChangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticEnergyChangerate");
                base.RegisterField(this.KineticEnergyChangerate, register);

                this.prevKineticEnergyChangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)));

            }

            if (this.Control.ComputeEnergyProperties) {

                if (this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
                    prevVel = new XDGField[D];
                    for (int d = 0; d < D; d++) {
                        prevVel[d] = new XDGField(this.XDGvelocity.Velocity[d].Basis);
                    }
                }

                this.prevKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "previousKineticEnergy");
                base.RegisterField(this.prevKineticEnergy);

                this.KineticEnergyNSE = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticEnergyNSE");
                base.RegisterField(this.KineticEnergyNSE);

                this.KineticEnergyNSEchangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticEnergyNSEchangerate");
                base.RegisterField(this.KineticEnergyNSEchangerate);

                this.DerivedKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "DerivedKineticEnergy");
                base.RegisterField(this.DerivedKineticEnergy, register);

                this.DerivedKineticEnergyChangerate = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "DerivedKineticEnergyChangerate");
                base.RegisterField(this.DerivedKineticEnergyChangerate, register);

                this.KineticDissipation = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticDissipation");
                base.RegisterField(this.KineticDissipation, register);

                this.PowerOfStresses = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "PowerOfStresses");
                base.RegisterField(this.PowerOfStresses, register);


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
            KineticEnergyOperator = new XSpatialOperatorMk2(DomName, Params, CodName, (A, B, C) => degK * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2), this.LsTrk.SpeciesIdS.ToArray());


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



        public void EvaluateKineticEnergy(double physTime, double[] fluxEval) {


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



            XSpatialOperatorMk2.XEvaluatorNonlin eval = this.KineticEnergyOperator.GetEvaluatorEx(this.LsTrk,
                this.KineticEnergy.ToEnumerable().ToArray(), Params, this.KineticEnergy.Mapping,
                this.LsTrk.SpeciesIdS.ToArray());

            var agg = this.LsTrk.GetAgglomerator(this.LsTrk.SpeciesIdS.ToArray(),
                this.KineticEnergy.Basis.Degree * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2),
                this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold);

            foreach (var kv in agg.CellLengthScales) {
                eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
            }

            //if (this.KineticEnergyOperator.SurfaceElementOperator.TotalNoOfComponents > 0) {
            //    foreach (var kv in InterfaceLengths) {
            //        eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("InterfaceLengths", kv.Value);
            //        //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("lambda_interface", lambdaI);
            //        //eval.SpeciesOperatorCoefficients[kv.Key].UserDefinedValues.Add("mu_interface", muI);
            //    }
            //}

            eval.time = physTime;

            eval.Evaluate(1.0, 1.0, fluxEval);

        }






        //void DelComputeEnergyOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {

        //    int D = this.GridData.SpatialDimension;

        //    SpeciesId[] SpcToCompute = AgglomeratedCellLengthScales.Keys.ToArray();

        //    // parameter assembly
        //    // ==================    

        //    // velocity
        //    var VelMap = new CoordinateMapping(this.XDGvelocity.Velocity.ToArray());
        //    DGField[] VelParam = VelMap.Fields.ToArray();

        //    // velocity mean
        //    VectorField<XDGField> VelMeanParam = new VectorField<XDGField>(D, new XDGBasis(LsTrk, 0), "VelMean_", XDGField.Factory);
        //    XheatUtils.ComputeAverageU(VelParam, VelMeanParam, m_HMForder, LsTrk.GetXDGSpaceMetrics(SpcToCompute, m_HMForder, 1).XQuadSchemeHelper, this.LsTrk);

        //    // velocity gradient vectors
        //    VectorField<DGField> GradVelX = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityXGradient", XDGField.Factory);
        //    for (int d = 0; d < D; d++) {
        //        foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //            DGField f_Spc = ((VelParam[0] as XDGField).GetSpeciesShadowField(Spc));
        //            (GradVelX[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
        //        }
        //    }
        //    GradVelX.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    VectorField<DGField> GradVelY = new VectorField<DGField>(D, VelParam[0].Basis, "VelocityYGradient", XDGField.Factory);
        //    for (int d = 0; d < D; d++) {
        //        foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //            DGField f_Spc = ((VelParam[1] as XDGField).GetSpeciesShadowField(Spc));
        //            (GradVelY[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
        //        }
        //    }
        //    GradVelY.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    // pressure and gradient
        //    var PressMap = new CoordinateMapping(this.Pressure);
        //    DGField[] PressParam = PressMap.Fields.ToArray();

        //    VectorField<DGField> PressGrad = new VectorField<DGField>(D, PressParam[0].Basis, "PressureGrad", XDGField.Factory);
        //    for (int d = 0; d < D; d++) {
        //        foreach (var Spc in this.LsTrk.SpeciesIdS) {
        //            DGField f_Spc = ((PressParam[0] as XDGField).GetSpeciesShadowField(Spc));
        //            (PressGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
        //        }
        //    }
        //    PressGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

        //    // gravity
        //    var GravMap = new CoordinateMapping(this.XDGvelocity.Gravity.ToArray());
        //    DGField[] GravParam = GravMap.Fields.ToArray();

        //    // normals:
        //    SinglePhaseField[] Normals; // Normal vectors: length not normalized - will be normalized at each quad node within the flux functions.
        //    var LevelSetGradient = new VectorField<SinglePhaseField>(D, LevSet.Basis, SinglePhaseField.Factory);
        //    LevelSetGradient.Gradient(1.0, LevSet);
        //    Normals = LevelSetGradient.ToArray();

        //    // Curvature
        //    CurvatureAlgorithms.CurvatureDriver(
        //        SurfaceStressTensor_IsotropicMode.Curvature_Projected,
        //        CurvatureAlgorithms.FilterConfiguration.Default,
        //        this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
        //        this.m_HMForder, this.DGLevSet.Current);


        //    // concatenate everything
        //    var Params = ArrayTools.Cat<DGField>(
        //        VelParam,
        //        VelMeanParam,
        //        GradVelX,
        //        GradVelY,
        //        PressParam,
        //        PressGrad,
        //        GravMap,
        //        Normals,
        //        this.Curvature);



        //    // assemble the matrix & affine vector
        //    // ===================================


        //    // compute matrix
        //    if (OpMtx != null) {

        //        var mtxBuilder = KineticEnergyBalanceOperator.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping, SpcToCompute);

        //        mtxBuilder.time = phystime;

        //        foreach (var kv in AgglomeratedCellLengthScales) {
        //            mtxBuilder.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
        //        }

        //        mtxBuilder.ComputeMatrix(OpMtx, OpAffine);

        //    } else {
        //        XSpatialOperatorMk2.XEvaluatorNonlin eval = KineticEnergyBalanceOperator.GetEvaluatorEx(LsTrk,
        //            CurrentState.ToArray(), Params, Mapping,
        //            SpcToCompute);

        //        foreach (var kv in AgglomeratedCellLengthScales) {
        //            eval.SpeciesOperatorCoefficients[kv.Key].CellLengthScales = kv.Value;
        //        }

        //        eval.time = phystime;

        //        eval.Evaluate(1.0, 1.0, OpAffine);

        //    }


        //OpAffine.ScaleV(-1.0);

        //// mass matrix factory
        //MassFact = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), m_HMForder, 1).MassMatrixFactory;// new MassMatrixFactory(maxB, CurrentAgg);
        //var WholeMassMatrix = MassFact.GetMassMatrix(Mapping, MassScale); // mass matrix scaled with density rho

        //// add power of gravity forces
        //var WholeGravity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity.Gravity.ToArray<DGField>()));
        //var WholeVelocity = new CoordinateVector(ArrayTools.Cat<DGField>(this.XDGvelocity.Velocity.ToArray<DGField>()));
        //WholeMassMatrix.SpMV(1.0, WholeVelocity, 1.0, OpAffine);

        //// transform from RHS to Affine
        //OpAffine.ScaleV(-1.0);

        //}

        /// <summary>
        /// dummy delegate for coupled operators
        /// </summary>
        /// <param name="CurrentState"></param>
        /// <param name="Phystime"></param>
        /// <param name="dt"></param>
        /// <param name="underrelax"></param>
        /// <param name="incremental"></param>
        /// <returns></returns>
        //double DelUpdateLevelSet_EnergyOperator(DGField[] CurrentState, double Phystime, double dt, double underrelax, bool incremental) {
        //    // do nothing
        //    return 0.0;
        //}


        /// <summary>
        /// The residual logger for this application.
        /// </summary>
        //public ResidualLogger EnergyResLogger {
        //    get {
        //        return m_EnergyResLogger;
        //    }
        //}

        //ResidualLogger m_EnergyResLogger;


        #endregion



    }

}
