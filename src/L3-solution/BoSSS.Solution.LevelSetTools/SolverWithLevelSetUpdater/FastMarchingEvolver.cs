using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {



    /// <summary>
    /// Driver class for the <see cref="NarrowMarchingBand.Evolve_Mk2"/>
    /// </summary>
    /// <remarks>
    /// In contrast to other methods, this should be reliable for general interface topologies
    /// and it is also known to work for contact line setups.
    /// The fast marching is a level-set evolution technique which 
    /// Development history and literature:
    /// - parallelized and extended to 3D by M. Smuda, Q4 2020.
    /// - used extensively in the thesis om M. Smuda, 2020
    /// - also described in the FDY Annual Report 2015, F. Kummer
    /// - in 2015/2016 the cell-wise solvers were re-formulated 
    ///   based upon the Reinitialization and Extension ideas from T. Utz
    /// - initial implementation from F Kummer, 2013 & 2014
    /// </remarks>
    public class FastMarchingEvolver : ILevelSetEvolver {
        SinglePhaseField[] extensionVelocity;

        int m_HMForder;

        IList<string> parameters;

        //string[] variables;

        string levelSetName;

        public FastMarchingEvolver(string levelSetName, int hMForder, int D, int ReInitPeriod) {
            this.m_HMForder = hMForder;
            this.levelSetName = levelSetName;

            this.ReInit_Period = ReInitPeriod;
            ReInit_Control = new EllipticReInitAlgoControl();

            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D);
            parameters = parameters.Cat(BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(this.levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)));
        }

        public IList<string> VariableNames => null;

        public IList<string> ParameterNames => parameters;

        // reinitialization
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => Reinitialize;

        /// <summary>
        /// Provides access to the internally constructed extension velocity.
        /// <see cref="ILevelSetEvolver.InternalFields"/>
        /// </summary>
        public IDictionary<string, DGField> InternalFields {
            get {
                var Ret = new Dictionary<string, DGField>();

                if (extensionVelocity != null) {
                    foreach (var f in extensionVelocity)
                        Ret.Add(f.Identification, f);
                }

                return Ret;
            }
        }


        public void MovePhaseInterface(
            DualLevelSet phaseInterface,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            int D = phaseInterface.Tracker.GridDat.SpatialDimension;


            SinglePhaseField[] meanVelocity = D.ForLoop(
                d => (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d))]
                );

            VectorField<SinglePhaseField> filtLevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(
                d => new SinglePhaseField(phaseInterface.DGLevelSet.Basis)
                ));

            if (extensionVelocity == null) {
                extensionVelocity = new SinglePhaseField[D];
                Basis basis;
                try {
                    basis = new Basis(phaseInterface.Tracker.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
                } catch {
                    Console.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                    basis = new Basis(phaseInterface.Tracker.GridDat, ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0X].Basis.Degree);
                }

                for (int d = 0; d < D; ++d) {
                    extensionVelocity[d] = new SinglePhaseField(basis, "ExtensionVelocity[" + d + "]");
                }
            }
            //Move LevelSet
            SinglePhaseField lsBuffer = phaseInterface.DGLevelSet.CloneAs();

            int TimestepNo = (int)(time / dt);
            //TimestepNumber tsn = new TimestepNumber(new int[] { TimestepNo, 0 });
            //DGField[] plotFields = ArrayTools.Cat<DGField>(meanVelocity, extensionVelocity);
            //Tecplot.Tecplot.PlotFields(plotFields, "NarrowMarchingBand" + tsn, time, 2);

            NarrowMarchingBand.Evolve_Mk2(
                dt, phaseInterface.Tracker, lsBuffer, phaseInterface.DGLevelSet, filtLevSetGradient,
                meanVelocity, extensionVelocity,
                m_HMForder,
                TimestepNo, false);

            //tsn = new TimestepNumber(new int[] { TimestepNo, 1 });
            //Tecplot.Tecplot.PlotFields(plotFields, "NarrowMarchingBand" + tsn, 0.0, 2);
            //Tecplot.Tecplot.PlotFields(plotFields, this.GetType().ToString().Split('.').Last() + "-" + TimestepNo, (double)TimestepNo, 2);            
        }

        private EllipticReInitAlgoControl ReInit_Control;
        private int ReInit_TimestepIndex = 0;
        private int ReInit_Period = 0;
        public void Reinitialize(
            DualLevelSet phaseInterface,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            // after level-set evolution and for initializing non-signed-distance level set fields
            if (ReInit_Period > 0 && ReInit_TimestepIndex % ReInit_Period == 0) {
                
                Console.WriteLine("Performing ReInit");
                ReInit_Control.Potential = ReInitPotential.BastingSingleWell;
                EllipticReInit.EllipticReInit ReInitPDE = new EllipticReInit.EllipticReInit(phaseInterface.Tracker, ReInit_Control, phaseInterface.DGLevelSet);
                ReInitPDE.ReInitialize(Restriction: phaseInterface.Tracker.Regions.GetCutCellSubGrid());

                FastMarchReinit FastMarchReinitSolver = new FastMarchReinit(phaseInterface.DGLevelSet.Basis);
                CellMask Accepted = phaseInterface.Tracker.Regions.GetCutCellMask();
                CellMask ActiveField = phaseInterface.Tracker.Regions.GetNearFieldMask(1);
                CellMask NegativeField = phaseInterface.Tracker.Regions.GetSpeciesMask("A");
                FastMarchReinitSolver.FirstOrderReinit(phaseInterface.DGLevelSet, Accepted, NegativeField, ActiveField);
            }

            ReInit_TimestepIndex++;
        }
    }
}
