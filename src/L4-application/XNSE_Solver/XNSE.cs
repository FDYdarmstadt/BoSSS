﻿using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver {

    public class XNSE : SolverWithLevelSetUpdater<XNSE_Control> {
        protected IncompressibleMultiphaseBoundaryCondMap boundaryMap;

        public int QuadOrder() {
            //QuadOrder
            int degU;
            if (Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts field)) {
                degU = field.Degree;
            } else if (Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts field1)) {
                degU = field1.Degree;
            } else {
                throw new Exception("Velocity not found!");
            }
            int quadOrder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if (this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                quadOrder *= 2;
                quadOrder += 1;
            }
            return quadOrder;
        }

        protected int VelocityDegree() {
            int pVel;
            if (this.Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v)) {
                pVel = v.Degree;
            } else if (this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1)) {
                pVel = v1.Degree;
            } else {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pVel;
        }

        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;

        protected override LevelSetUpdater InstantiateLevelSetUpdater() {
            int levelSetDegree = Control.FieldOptions["Phi"].Degree;

            LevelSetUpdater lsUpdater;
            ILevelSetParameter levelSetVelocity = new LevelSetVelocity("Phi", GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
            switch (Control.Option_LevelSetEvolution) {
                case LevelSetEvolution.Fourier:
                if (Control.EnforceLevelSetConservation) {
                    throw new NotSupportedException("mass conservation correction currently not supported");
                }
                FourrierLevelSet fourrierLevelSet = new FourrierLevelSet(Control.FourierLevSetControl, new Basis(GridData, levelSetDegree), "Phi");
                fourrierLevelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, fourrierLevelSet);
                lsUpdater.AddLevelSetParameter("Phi", levelSetVelocity);
                break;

                case LevelSetEvolution.FastMarching:
                LevelSet levelSet = new LevelSet(new Basis(GridData, levelSetDegree), "Phi");
                levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet);
                var fastMarcher = new FastMarchingEvolver("Phi", QuadOrder(), levelSet.GridDat.SpatialDimension);
                lsUpdater.AddEvolver("Phi", fastMarcher);
                lsUpdater.AddLevelSetParameter("Phi", levelSetVelocity);
                break;

                case LevelSetEvolution.SplineLS:
                int nodeCount = 30;
                Console.WriteLine("Achtung, Spline node count ist hart gesetzt. Was soll hier hin?");
                SplineLevelSet SplineLevelSet = new SplineLevelSet(Control.Phi0Initial, new Basis(GridData, levelSetDegree), "Phi", nodeCount);
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, SplineLevelSet);
                var SplineEvolver = new SplineLevelSetEvolver("Phi", (GridData)SplineLevelSet.GridDat);
                lsUpdater.AddEvolver("Phi", SplineEvolver);
                lsUpdater.AddLevelSetParameter("Phi", levelSetVelocity);
                break;

                case LevelSetEvolution.None:
                LevelSet levelSet1 = new LevelSet(new Basis(GridData, levelSetDegree), "Phi");
                levelSet1.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet1);
                break;

                default:
                throw new NotImplementedException();
            }
            return lsUpdater;
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            int pVel = VelocityDegree();
            int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
            int D = this.GridData.SpatialDimension;

            if (this.Control.UseSchurBlockPrec) {
                // using a Schur complement for velocity & pressure
                var confMomConti = new MultigridOperator.ChangeOfBasisConfig();
                for (int d = 0; d < D; d++) {
                    d.AddToArray(ref confMomConti.VarIndex);
                    //Math.Max(1, pVel - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                    pVel.AddToArray(ref confMomConti.DegreeS);
                }
                D.AddToArray(ref confMomConti.VarIndex);
                //Math.Max(0, pPrs - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                pPrs.AddToArray(ref confMomConti.DegreeS);

                confMomConti.mode = MultigridOperator.Mode.SchurComplement;

                configsLevel.Add(confMomConti);
            } else {
                // configurations for velocity
                for (int d = 0; d < D; d++) {
                    var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pVel },
                        //DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                        mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.VelocityVector(D)[d]) }
                    };
                    configsLevel.Add(configVel_d);
                }
                // configuration for pressure
                var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pPrs },
                    //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                    mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Pressure) }
                };
                configsLevel.Add(configPres);
            }
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {
            OperatorFactory opFactory = new OperatorFactory();
            boundaryMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            XOP.FreeMeanValue[BoSSS.Solution.NSECommon.VariableNames.Pressure] = !boundaryMap.DirichletPressureBoundary;
            XOP.LinearizationHint = LinearizationHint.AdHoc;
            XOP.Commit();

            return XOP;
        }

        protected virtual void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            int quadOrder = QuadOrder();
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            for (int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NavierStokes("A", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A));
                opFactory.AddEquation(new NavierStokes("B", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B));
                opFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, LsTrk, config, config.isMovingMesh));
                opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, LsTrk, config));
            }
            opFactory.AddCoefficient(new SlipLengths(config, VelocityDegree()));
            opFactory.AddParameter(new Velocity0(D));
            Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
            opFactory.AddParameter(v0Mean);

            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);

            if (config.isContinuity) {
                opFactory.AddEquation(new Continuity(config, D, "A", LsTrk.GetSpeciesId("A"), boundaryMap));
                opFactory.AddEquation(new Continuity(config, D, "B", LsTrk.GetSpeciesId("B"), boundaryMap));
                opFactory.AddEquation(new InterfaceContinuity(config, D, LsTrk, config.isMatInt));
            }

            lsUpdater.AddLevelSetParameter("Phi", v0Mean);
            lsUpdater.AddLevelSetParameter("Phi", normalsParameter);
            switch (Control.AdvancedDiscretizationOptions.SST_isotropicMode) {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine:
                MaxSigma maxSigmaParameter = new MaxSigma(Control.PhysicalParameters, Control.AdvancedDiscretizationOptions, QuadOrder(), Control.dtFixed);
                opFactory.AddParameter(maxSigmaParameter);
                lsUpdater.AddLevelSetParameter("Phi", maxSigmaParameter);
                BeltramiGradient lsBGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                lsUpdater.AddLevelSetParameter("Phi", lsBGradient);
                break;

                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                BeltramiGradient lsGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                lsUpdater.AddLevelSetParameter("Phi", lsGradient);
                break;

                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                BeltramiGradientAndCurvature lsGradientAndCurvature =
                    FromControl.BeltramiGradientAndCurvature(Control, "Phi", quadOrder, D);
                opFactory.AddParameter(lsGradientAndCurvature);
                lsUpdater.AddLevelSetParameter("Phi", lsGradientAndCurvature);
                break;

                case SurfaceStressTensor_IsotropicMode.Curvature_Fourier:
                FourrierLevelSet ls = (FourrierLevelSet)lsUpdater.LevelSets["Phi"].DGLevelSet;
                var fourrier = new FourierEvolver(
                    ls,
                    Control.FourierLevSetControl,
                    Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Curvature].Degree);
                lsUpdater.AddLevelSetParameter("Phi", fourrier);
                lsUpdater.AddEvolver("Phi", fourrier);
                opFactory.AddParameter(fourrier);
                break;

                default:
                break;
            }
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {
            Tecplot.PlotFields(this.m_RegisteredFields, "XNSE_Solver" + timestepNo, physTime, superSampling);
            if (Timestepping?.Parameters != null) {
                Tecplot.PlotFields(Timestepping.Parameters, "XNSE_Solver_Time_Params" + timestepNo, physTime, superSampling);
            }
            if (LsUpdater?.Parameters != null) {
                Tecplot.PlotFields(Timestepping.Parameters, "XNSE_Solver_LevelSet_Params" + timestepNo, physTime, superSampling);
            }
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetFixedTimestep();
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            Console.WriteLine($"done with timestep {TimestepNo}");
            return dt;
        }
    }
}