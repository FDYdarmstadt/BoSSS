using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver
{
    interface ILevelSetEvolver
    {
        void UpdatePhaseInterface(
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

        public FastMarcher(XNSE_Control control, int hMForder)
        {
            this.Control = control;
            this.m_HMForder = hMForder;
        }

        public void UpdatePhaseInterface(
            DualLevelSet phaseInterface,
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
            ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(EvoVelocity, lsTrkr, Control);
            
            //Extension Velocity
            if(extensionVelocity == null)
            {
                int D = lsTrkr.GridDat.SpatialDimension;
                extensionVelocity = new SinglePhaseField[D];
                Basis basis = new Basis(lsTrkr.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
                for(int d = 0; d < D; ++d)
                {
                    extensionVelocity[d] = new SinglePhaseField(basis, "ExtensionVelocity" + d);
                }
            }

            //dgLevSetGradient and update curvature
            VectorField<SinglePhaseField> filtLevSetGradient;
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature];

            switch (Control.AdvancedDiscretizationOptions.SST_isotropicMode)
            {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine:
                    {
                        CurvatureAlgorithms.LaplaceBeltramiDriver(
                            Control.AdvancedDiscretizationOptions.SST_isotropicMode,
                            Control.AdvancedDiscretizationOptions.FilterConfiguration,
                            out filtLevSetGradient, 
                            lsTrkr,
                            phaseInterface.DGLevelSet);
                        break;
                    }
                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                    CurvatureAlgorithms.CurvatureDriver(
                        Control.AdvancedDiscretizationOptions.SST_isotropicMode,
                        Control.AdvancedDiscretizationOptions.FilterConfiguration,
                        Curvature, 
                        out filtLevSetGradient, 
                        lsTrkr,
                        this.m_HMForder,
                        phaseInterface.DGLevelSet);
                    //CurvatureAlgorithms.MakeItConservative(LsTrk, this.Curvature, this.Control.PhysicalParameters.Sigma, this.SurfaceForce, filtLevSetGradient, MomentFittingVariant, this.m_HMForder);
                    break;
                default: throw new NotImplementedException("Unknown SurfaceTensionMode");
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
                    switch (control.InterAverage)
                    {
                        case XNSE_Control.InterfaceAveraging.mean:
                            {
                                scale = 0.5;
                                break;
                            }
                        case XNSE_Control.InterfaceAveraging.density:
                            {
                                scale = rhoSpc / (rho_A + rho_B);
                                break;
                            }
                        case XNSE_Control.InterfaceAveraging.viscosity:
                            {
                                scale = muSpc / (mu_A + mu_B);
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

    class FourierEvolver : ILevelSetEvolver
    {

        FourierLevSetBase Fourier_LevSet;

        /// <summary>
        /// specialized timestepper (Runge-Kutta-based) for the evoultion of the Fourier-LS
        /// </summary>
        FourierLevSetTimestepper Fourier_Timestepper;

        XNSE_Control Control;

        int m_HMForder;

        public FourierEvolver(XNSE_Control Control, int hMForder)
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
        }

        public void UpdatePhaseInterface(
            DualLevelSet levelSet, 
            LevelSetTracker lsTrkr, 
            double time, 
            double dt, 
            bool incremental, 
            IReadOnlyDictionary<string, DGField> DomainVarFields, 
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {

            VectorField<SinglePhaseField> filtLevSetGradient;
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature];

            Fourier_LevSet.ProjectToDGCurvature(Curvature, out filtLevSetGradient, lsTrkr.Regions.GetCutCellMask4LevSet(levelSet.LevelSetIndex));

            Fourier_Timestepper.updateFourierLevSet();
            Fourier_LevSet.ProjectToDGLevelSet(levelSet.DGLevelSet, lsTrkr);

            //Mean Velocity
            XDGField[] EvoVelocity = new XDGField[]
            {
                (XDGField) DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX],
                (XDGField) DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityY],
            };
            ConventionalDGField[] meanVelocity = FastMarcher.GetMeanVelocityFromXDGField(EvoVelocity, lsTrkr, Control);

            Fourier_Timestepper.moveLevelSet(dt, meanVelocity);
            //if (incremental)
            //    Fourier_Timestepper.updateFourierLevSet();
            Fourier_LevSet.ProjectToDGLevelSet(levelSet.DGLevelSet, lsTrkr);
        }
    }
}
