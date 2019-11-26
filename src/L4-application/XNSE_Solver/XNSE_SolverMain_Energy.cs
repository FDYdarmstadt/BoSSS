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

        XDGField GeneratedKineticEnergy;

        /// <summary>
        /// kinetic energy computed via <see cref="KineticEnergyBalanceOperator"/>
        /// </summary>
        XDGField KineticEnergy;

        //XDGField prevKineticEnergy;

        /// <summary>
        /// Residual of the kinetic energy balance
        /// </summary>
        XDGField ResidualKineticEnergy;

        /// <summary>
        /// source term for the kinetic energy
        /// </summary>
        XDGField KineticDissipation;

#pragma warning restore 649

        /// <summary>
        /// creates energy related fields
        /// </summary>
        public void CreateEnergyFields() {

            int D = this.GridData.SpatialDimension;


            if (this.Control.solveKineticEnergyEquation) {

                this.KineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), VariableNames.KineticEnergy);
                base.RegisterField(this.KineticEnergy);
                this.ResidualKineticEnergy = new XDGField(this.KineticEnergy.Basis, "ResidualKineticEnergy");
                base.RegisterField(this.ResidualKineticEnergy);

                //this.prevKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions["KineticEnergy"].Degree)));

                this.DerivedKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "DerivedKineticEnergy");
                base.RegisterField(this.DerivedKineticEnergy);

                this.GeneratedKineticEnergy = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticEnergyProduction");
                base.RegisterField(this.GeneratedKineticEnergy);

            }

            if (this.Control.ComputeEnergyProperties) {

                if (this.Control.CompMode == AppControl._CompMode.Transient) {
                    prevVel = new XDGField[D];
                    for (int d = 0; d < D; d++) {
                        prevVel[d] = new XDGField(this.XDGvelocity.Velocity[d].Basis);
                    }
                }


                this.KineticDissipation = new XDGField(new XDGBasis(this.LsTrk, (this.Control.FieldOptions[VariableNames.KineticEnergy].Degree)), "KineticDissipation");
                base.RegisterField(this.KineticDissipation);

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
        //XSpatialOperatorMk2 KineticEnergyBalanceOperator;

        //IDictionary<SpeciesId, IEnumerable<double>> MassScaleForEnergy {
        //    get {
        //        double rho_A = this.Control.PhysicalParameters.rho_A,
        //            rho_B = this.Control.PhysicalParameters.rho_B;

        //        double[] _rho_A = new double[1];
        //        _rho_A[0] = rho_A;
        //        double[] _rho_B = new double[1];
        //        _rho_B[0] = rho_B;

        //        Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
        //        R.Add(this.LsTrk.GetSpeciesId("A"), _rho_A);
        //        R.Add(this.LsTrk.GetSpeciesId("B"), _rho_B);

        //        return R;
        //    }
        //}

        //MultigridOperator.ChangeOfBasisConfig[][] MultigridEnergyOperatorConfig {
        //    get {
        //        int pEnergy = this.KineticEnergy.Basis.Degree;

        //        // set the MultigridOperator configuration for each level:
        //        // it is not necessary to have exactly as many configurations as actual multigrid levels:
        //        // the last configuration enty will be used for all higher level
        //        MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[1][];
        //        for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
        //            configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[1];

        //            // configuration for Temperature
        //            configs[iLevel][0] = new MultigridOperator.ChangeOfBasisConfig() {
        //                Degree = Math.Max(0, pEnergy - iLevel),
        //                mode = MultigridOperator.Mode.Eye,
        //                VarIndex = new int[] { 0 }
        //            };
        //        }

        //        return configs;
        //    }
        //}


        //EnergyMultiphaseBoundaryCondMap m_energyBcMap;

        ///// <summary>
        ///// Boundary conditions.
        ///// </summary>
        //EnergyMultiphaseBoundaryCondMap energyBcMap {
        //    get {
        //        if (m_energyBcMap == null) {
        //            m_energyBcMap = new EnergyMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
        //        }
        //        return m_energyBcMap;
        //    }
        //}

        //CoordinateVector m_CurrentEnergySolution;

        ///// <summary>
        ///// Current temperature;
        ///// </summary>
        //internal CoordinateVector CurrentEnergySolution {
        //    get {
        //        if (m_CurrentEnergySolution == null) {
        //            m_CurrentEnergySolution = new CoordinateVector(this.KineticEnergy);
        //        }
        //        return m_CurrentEnergySolution;
        //    }
        //}

        //CoordinateVector m_CurrentEnergyResidual;

        ///// <summary>
        ///// Current residual for coupled heat equation.
        ///// </summary>
        //internal CoordinateVector CurrentEnergyResidual {
        //    get {
        //        if (m_CurrentEnergyResidual == null) {
        //            m_CurrentEnergyResidual = new CoordinateVector(this.ResidualKineticEnergy);
        //        }
        //        return m_CurrentEnergyResidual;
        //    }
        //}


        /// <summary>
        /// Implicit timestepping using Backward-Differentiation-Formulas (BDF),
        /// specialized for XDG applications.
        /// </summary>
        //XdgBDFTimestepping m_BDF_energyTimestepper;


        //public void generateKinEnergyOperator() {

        //    int degK = this.KineticEnergy.Basis.Degree;

        //    int D = this.GridData.SpatialDimension;

        //    string[] CodName = new string[] { "kinBalance" };
        //    string[] Params = ArrayTools.Cat(
        //         VariableNames.VelocityVector(D),
        //         (new string[] { "VelocityX_Mean", "VelocityY_Mean", "VelocityZ_Mean" }).GetSubVector(0, D),
        //         VariableNames.VelocityX_GradientVector(),
        //         VariableNames.VelocityY_GradientVector(),
        //         VariableNames.Pressure,
        //         (new string[] { "PressureGradX", "PressureGradY", "PressureGradZ" }).GetSubVector(0, D),
        //         (new string[] { "GravityX", "GravityY", "GravityZ" }).GetSubVector(0, D),
        //         (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D),
        //         "Curvature");
        //    string[] DomName = new string[] { "KineticEnergy" };

        //    double rhoA = this.Control.PhysicalParameters.rho_A;
        //    double rhoB = this.Control.PhysicalParameters.rho_B;
        //    double muA = this.Control.PhysicalParameters.mu_A;
        //    double muB = this.Control.PhysicalParameters.mu_B;
        //    double sigma = this.Control.PhysicalParameters.Sigma;

        //    double LFFA = this.Control.AdvancedDiscretizationOptions.LFFA;
        //    double LFFB = this.Control.AdvancedDiscretizationOptions.LFFB;


        //    var dntParams = this.Control.AdvancedDiscretizationOptions;

        //    // create operator
        //    // ===============
        //    KineticEnergyBalanceOperator = new XSpatialOperatorMk2(DomName, Params, CodName, (A, B, C) => degK * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2), null);


        //    // build the operator
        //    // ==================
        //    {

        //        // convective part
        //        // ================
        //        {
        //            if (this.Control.PhysicalParameters.IncludeConvection) {

        //                var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];

        //                // kinetic energy
        //                var convK = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyConvectionInBulk(D, energyBcMap, rhoA, rhoB, LFFA, LFFB, LsTrk);
        //                comps.Add(convK); // Bulk component


        //                bool movingmesh;
        //                switch (this.Control.Timestepper_LevelSetHandling) {
        //                    case LevelSetHandling.Coupled_Once:
        //                        movingmesh = true;
        //                        break;
        //                    case LevelSetHandling.LieSplitting:
        //                    case LevelSetHandling.StrangSplitting:
        //                    case LevelSetHandling.None:
        //                        movingmesh = false;
        //                        break;
        //                    case LevelSetHandling.Coupled_Iterative:
        //                    default:
        //                        throw new NotImplementedException();
        //                }

        //                //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyConvectionAtLevelSet(D, LsTrk, rhoA, rhoB, LFFA, LFFB, this.Control.PhysicalParameters.Material, energyBcMap, movingmesh));       // LevelSet component
        //            }
        //        }

        //        // Laplace of kinetic energy
        //        // =========================
        //        {
        //            var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];

        //            double penalty = dntParams.PenaltySafety;

        //            var Visc = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyLaplace(
        //                dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
        //                energyBcMap, D, muA, muB);

        //            comps.Add(Visc);

        //            if (dntParams.UseGhostPenalties) {
        //                var ViscPenalty = new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergyLaplace(penalty * 1.0, 0.0, energyBcMap, D, muA, muB);
        //                KineticEnergyBalanceOperator.GhostEdgesOperator.EquationComponents[CodName[0]].Add(ViscPenalty);
        //            }

        //            // Level-Set operator:
        //            //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.KineticEnergylaplceAtLevelSet(LsTrk, muA, muB, penalty * 1.0));
        //        }

        //        // Divergence of stress tensor
        //        // ===========================
        //        {
        //            var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
        //            comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.StressDivergence(D, energyBcMap, muA, muB));

        //            // Level-Set operator:
        //            //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.StressDivergenceAtLevelSet(LsTrk, muA, muB));
        //        }

        //        // surface energy (surface tension)
        //        // ================================
        //        {
        //            //var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
        //            //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.SurfaceEnergy(D, LsTrk, sigma));
        //        }

        //        // pressure term
        //        // =============
        //        {
        //            var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
        //            comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.DivergencePressureEnergy(D, energyBcMap));
        //            //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.PressureConvectionInBulk(D, energyBcMap, LFFA, LFFB, LsTrk));
        //            //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.PressureGradientConvection(D));

        //            // Level-Set operator:
        //            //comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.DivergencePressureEnergyAtLevelSet(LsTrk));
        //        }

        //        // dissipation
        //        // ===========
        //        {
        //            var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
        //            comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.Dissipation(D, muA, muB));
        //        }

        //        // gravity (volume forces)
        //        // =======================
        //        {
        //            var comps = KineticEnergyBalanceOperator.EquationComponents[CodName[0]];
        //            comps.Add(new BoSSS.Solution.XNSECommon.Operator.Energy.PowerofGravity(D, rhoA, rhoB));
        //        }


        //        // finalize
        //        // ========

        //        KineticEnergyBalanceOperator.Commit();

        //    }

        //}


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
