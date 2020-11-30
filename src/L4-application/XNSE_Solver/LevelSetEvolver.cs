﻿using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver
{
    interface ILevelSetParameter
    {
        IList<string> ParameterNames { get; }

        (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields);

        void LevelSetParameterUpdate(
            DualLevelSet levelSet,
            LevelSetTracker lsTrkr,
            double time,
            IReadOnlyDictionary<string, DGField> ParameterVarFields);
    }

    interface ILevelSetEvolver
    {
        IList<string> ParameterNames { get; }

        IList<string> VariableNames { get; }

        void MovePhaseInterface(
            DualLevelSet levelSet,
            LevelSetTracker lsTrkr,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields);
    }

    class FastMarcher : ILevelSetEvolver
    {
        SinglePhaseField[] extensionVelocity;

        int m_HMForder;

        XNSE_Control Control;

        string[] parameters;

        string[] variables;

        public FastMarcher(XNSE_Control control, int hMForder, int D)
        {
            Control = control;
            this.m_HMForder = hMForder;
            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            variables = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
        }

        public IList<string> VariableNames => variables;

        public IList<string> ParameterNames => parameters;

        public void MovePhaseInterface(
            DualLevelSet phaseInterface,
            LevelSetTracker lsTrkr,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            //Mean Velocity
            XDGField[] EvoVelocity;
            try {
                EvoVelocity = new XDGField[]
                {
                    (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX],
                    (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityY],
                };
            } catch (KeyNotFoundException e) {
                Console.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                EvoVelocity = new XDGField[]
                {                    
                    (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0X],
                    (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0Y],
                };
            }

            ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(EvoVelocity, lsTrkr, Control);

            SinglePhaseField[] filtLevSetGradientArrray = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient0],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient1],
            };

            VectorField<SinglePhaseField> filtLevSetGradient = new VectorField<SinglePhaseField>(filtLevSetGradientArrray);

            //Extension Velocity: Was macht man damit?
            if (extensionVelocity == null)
            {
                int D = lsTrkr.GridDat.SpatialDimension;
                extensionVelocity = new SinglePhaseField[D];
                Basis basis;
                try {
                    basis = new Basis(lsTrkr.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
                } catch (KeyNotFoundException e) {
                    Console.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                    basis = new Basis(lsTrkr.GridDat, ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0X].Basis.Degree);
                }
                for (int d = 0; d < D; ++d)
                {
                    extensionVelocity[d] = new SinglePhaseField(basis, "ExtensionVelocity" + d);
                }
            }
            
            //Move LevelSet
            SinglePhaseField lsBuffer = phaseInterface.DGLevelSet.CloneAs();

            NarrowMarchingBand.Evolve_Mk2(
                dt, lsTrkr, lsBuffer, phaseInterface.DGLevelSet, filtLevSetGradient,
                meanVelocity, extensionVelocity,
                m_HMForder);
        }

        public static ConventionalDGField[] GetMeanVelocityFromXDGField(XDGField[] EvoVelocity, LevelSetTracker lsTrkr, XNSE_Control control)
        {
            int D = EvoVelocity.Length;
            ConventionalDGField[] meanVelocity;

            meanVelocity = new ConventionalDGField[D];

            double rho_A = control.PhysicalParameters.rho_A, rho_B = control.PhysicalParameters.rho_B;
            double mu_A = control.PhysicalParameters.mu_A, mu_B = control.PhysicalParameters.mu_B;
            CellMask CC = lsTrkr.Regions.GetCutCellMask4LevSet(0);
            CellMask Neg = lsTrkr.Regions.GetLevelSetWing(0, -1).VolumeMask;
            CellMask Pos = lsTrkr.Regions.GetLevelSetWing(0, +1).VolumeMask;
            CellMask posNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Neg);
            CellMask negNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Pos);

            for (int d = 0; d < D; d++)
            {
                Basis b = EvoVelocity[d].Basis.NonX_Basis;
                meanVelocity[d] = new SinglePhaseField(b);


                foreach (string spc in lsTrkr.SpeciesNames)
                {
                    double rhoSpc;
                    double muSpc;
                    switch (spc)
                    {
                        case "A": rhoSpc = rho_A; muSpc = mu_A; break;
                        case "B": rhoSpc = rho_B; muSpc = mu_B; break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }

                    double scale = 1.0;
                    switch (control.InterVelocAverage)
                    {
                        case XNSE_Control.InterfaceVelocityAveraging.mean:
                            {
                                scale = 0.5;
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.density:
                            {
                                scale = rhoSpc / (rho_A + rho_B);
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.viscosity:
                            {
                                scale = muSpc / (mu_A + mu_B);
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.phaseA:
                            {
                                scale = (spc == "A") ? 1.0 : 0.0;
                                break;
                            }
                        case XNSE_Control.InterfaceVelocityAveraging.phaseB:
                            {
                                scale = (spc == "B") ? 1.0 : 0.0;
                                break;
                            }
                    }

                    meanVelocity[d].Acc(scale, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), CC);
                    switch (spc)
                    {
                        //case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), Neg.Except(CC)); break;
                        case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), negNear); break;
                        case "B": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), posNear); break;
                        default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                    }
                }
            }
            return meanVelocity;
        }
    }

    class FourierEvolver : Parameter, ILevelSetParameter, ILevelSetEvolver
    {
        FourierLevSetBase Fourier_LevSet;

        /// <summary>
        /// specialized timestepper (Runge-Kutta-based) for the evoultion of the Fourier-LS
        /// </summary>
        FourierLevSetTimestepper Fourier_Timestepper;

        XNSE_Control Control;

        int m_HMForder;

        int curvatureDegree;

        public IList<string> VariableNames => new[] { BoSSS.Solution.NSECommon.VariableNames.VelocityX, BoSSS.Solution.NSECommon.VariableNames.VelocityY };

        public override IList<string> ParameterNames => new[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };

        public override DelParameterFactory Factory => ParameterFactory;

        public FourierEvolver(XNSE_Control Control, int hMForder, int curvatureDegree)
        {
            m_HMForder = hMForder;
            if (Control.FourierLevSetControl == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of FourierLevSetControl!");

            Fourier_LevSet = FourierLevelSetFactory.Build(Control.FourierLevSetControl);
            //create specialized fourier timestepper
            Fourier_Timestepper = FourierLevelSetFactory.Build_Timestepper(Control.FourierLevSetControl, Fourier_LevSet.GetFLSproperty(),
                                                        Fourier_LevSet.ComputeChangerate, Fourier_LevSet.EvolveFourierLS);
            if (Control.EnforceLevelSetConservation)
            {
                throw new NotSupportedException("mass conservation correction currently not supported");
            }
            this.Control = Control;
            this.curvatureDegree = curvatureDegree;
        }

        public void LevelSetParameterUpdate(
            DualLevelSet levelSet,
            LevelSetTracker lsTrkr,
            double time,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            Fourier_Timestepper.updateFourierLevSet();
            VectorField<SinglePhaseField> filtLevSetGradient;
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature];
            Fourier_LevSet.ProjectToDGCurvature(Curvature, out filtLevSetGradient, lsTrkr.Regions.GetCutCellMask4LevSet(levelSet.LevelSetIndex));
        }

        public void MovePhaseInterface(
            DualLevelSet levelSet,
            LevelSetTracker lsTrkr,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            //Mean Velocity
            XDGField[] EvoVelocity = new XDGField[]
            {
                (XDGField) DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX],
                (XDGField) DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityY],
            };
            ConventionalDGField[] meanVelocity = FastMarcher.GetMeanVelocityFromXDGField(EvoVelocity, lsTrkr, Control);

            Fourier_Timestepper.moveLevelSet(dt, meanVelocity, lsTrkr.Regions.GetNearFieldMask(1));
            if (incremental)
                Fourier_Timestepper.updateFourierLevSet();
            Fourier_LevSet.ProjectToDGLevelSet(levelSet.DGLevelSet, lsTrkr);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            //Curvature
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[1];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis curvatureBasis = new Basis(gridData, curvatureDegree);
            string curvatureName = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            fields[0] = (curvatureName, new SinglePhaseField(curvatureBasis, curvatureName));
            return fields;
        }
    }
}
