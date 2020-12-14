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
    public interface ILevelSetEvolver {
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
}
