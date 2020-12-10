using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    public class FastMarchingEvolver : ILevelSetEvolver
    {
        SinglePhaseField[] extensionVelocity;

        int m_HMForder;

        IList<string> parameters;

        string[] variables;

        string levelSetName;

        public FastMarchingEvolver(string levelSetName, int hMForder, int D)
        {
            this.m_HMForder = hMForder;
            this.levelSetName = levelSetName;
            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            parameters = parameters.Cat(BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(this.levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)));
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

            VectorField<SinglePhaseField> filtLevSetGradient = new VectorField<SinglePhaseField>(
                phaseInterface.DGLevelSet.GridDat.SpatialDimension.ForLoop(d => new SinglePhaseField(phaseInterface.DGLevelSet.Basis)));

            if (extensionVelocity == null)
            {
                int D = phaseInterface.Tracker.GridDat.SpatialDimension;
                extensionVelocity = new SinglePhaseField[D];
                Basis basis;
                try {
                    basis = new Basis(phaseInterface.Tracker.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
                } catch {
                    Console.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                    basis = new Basis(phaseInterface.Tracker.GridDat, ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0X].Basis.Degree);
                }

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
    }

    public class FourierEvolver : Parameter, ILevelSetParameter, ILevelSetEvolver
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
