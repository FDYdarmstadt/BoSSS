using BoSSS.Application.XNSE_Solver;
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
            double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields);
    }

    interface ILevelSetEvolver
    {
        IList<string> ParameterNames { get; }

        IList<string> VariableNames { get; }

        void MovePhaseInterface(
            DualLevelSet levelSet,
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

        IList<string> parameters;

        string[] variables;

        string levelSetName;

        public FastMarcher(string levelSetName, int hMForder, int D)
        {
            this.m_HMForder = hMForder;
            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            parameters = parameters.Cat(BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)));
        }

        public IList<string> VariableNames => null;

        public IList<string> ParameterNames => parameters;

        public void MovePhaseInterface(
            DualLevelSet phaseInterface,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            SinglePhaseField[] meanVelocity = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityX)],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName,BoSSS.Solution.NSECommon.VariableNames.VelocityY)],
            };

            SinglePhaseField[] filtLevSetGradientArrray = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient0],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient1],
            };

            VectorField<SinglePhaseField> filtLevSetGradient = new VectorField<SinglePhaseField>(filtLevSetGradientArrray);

            //Extension Velocity: Was macht man damit?
            if (extensionVelocity == null)
            {
                int D = phaseInterface.Tracker.GridDat.SpatialDimension;
                extensionVelocity = new SinglePhaseField[D];
                Basis basis = new Basis(phaseInterface.Tracker.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
                for (int d = 0; d < D; ++d)
                {
                    extensionVelocity[d] = new SinglePhaseField(basis, "ExtensionVelocity" + d);
                }
            }
            //Move LevelSet
            SinglePhaseField lsBuffer = phaseInterface.DGLevelSet.CloneAs();

            NarrowMarchingBand.Evolve_Mk2(
                dt, phaseInterface.Tracker, lsBuffer, phaseInterface.DGLevelSet, filtLevSetGradient,
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

        int curvatureDegree;

        string levelSetName;

        IList<string> parameters;

        public IList<string> VariableNames => null;

        public override IList<string> ParameterNames => parameters;

        IList<string> ILevelSetParameter.ParameterNames
        {
            get
            {
                return new[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };
            }
        }

        IList<string> ILevelSetEvolver.ParameterNames
        {
            get
            {
                return BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(2));
            }
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public FourierEvolver(string levelSetName, FourierLevSetControl control, int curvatureDegree)
        {
            parameters = new[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };
            this.levelSetName = levelSetName;

            if (control == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of FourierLevSetControl!");

            Fourier_LevSet = FourierLevelSetFactory.Build(control);
            //create specialized fourier timestepper
            Fourier_Timestepper = FourierLevelSetFactory.Build_Timestepper(control, Fourier_LevSet.GetFLSproperty(),
                Fourier_LevSet.ComputeChangerate, Fourier_LevSet.EvolveFourierLS);
            
            this.curvatureDegree = curvatureDegree;
        }

        public void LevelSetParameterUpdate(
            DualLevelSet levelSet,
            double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            Fourier_Timestepper.updateFourierLevSet();
            VectorField<SinglePhaseField> filtLevSetGradient;
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature];
            Fourier_LevSet.ProjectToDGCurvature(Curvature, out filtLevSetGradient, levelSet.Tracker.Regions.GetCutCellMask4LevSet(levelSet.LevelSetIndex));
        }

        public void MovePhaseInterface(
            DualLevelSet levelSet,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            SinglePhaseField[] meanVelocity = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityX)],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName,BoSSS.Solution.NSECommon.VariableNames.VelocityY)],
            };

            Fourier_Timestepper.moveLevelSet(dt, meanVelocity, levelSet.Tracker.Regions.GetNearFieldMask(1));
            if (incremental)
                Fourier_Timestepper.updateFourierLevSet();
            Fourier_LevSet.ProjectToDGLevelSet(levelSet.DGLevelSet, levelSet.Tracker);
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
