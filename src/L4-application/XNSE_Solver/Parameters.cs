using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using ilPSP.Tracing;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.XNSE_Solver {
    public static class FromControl {
        public static BeltramiGradient BeltramiGradient(XNSE_Control control, string levelSetName, int D) {
            string levelSet = levelSetName;
            int levelSetDegree;
            if (control.FieldOptions.TryGetValue(levelSet, out FieldOpts lsOpts)) {
                var levelSetSource = control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource;
                levelSetDegree = (levelSetSource == CurvatureAlgorithms.LevelSetSource.fromDG) ? lsOpts.Degree : lsOpts.Degree + 1;
            } else {
                throw new Exception("Level set options not found in FieldOptions");
            }

            DoNotTouchParameters AdvancedDiscretizationOptions = control.AdvancedDiscretizationOptions;
            return new BeltramiGradient(D, AdvancedDiscretizationOptions, levelSetDegree);
        }

        public static BeltramiGradientAndCurvature BeltramiGradientAndCurvature(XNSE_Control control, string levelSetName, int m_HMForder, int D) {
            string curvature = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            int curvatureDegree;
            if (control.FieldOptions.TryGetValue(curvature, out FieldOpts opts)) {
                curvatureDegree = opts.Degree;
            } else {
                throw new Exception("Curvature options not found in FieldOptions");
            }
            
            string levelSet = levelSetName;
            int levelSetDegree;
            if (control.FieldOptions.TryGetValue(levelSet, out FieldOpts lsOpts)) {
                var levelSetSource = control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource;
                levelSetDegree = (levelSetSource == CurvatureAlgorithms.LevelSetSource.fromDG) ? lsOpts.Degree : lsOpts.Degree + 1;
            } else {
                levelSetDegree = 1;
            }
            DoNotTouchParameters AdvancedDiscretizationOptions = control.AdvancedDiscretizationOptions;
            return new BeltramiGradientAndCurvature(curvatureDegree, levelSetDegree, m_HMForder, AdvancedDiscretizationOptions, D);
        }
    }

    /// <summary>
    /// Computation of the fluid interface velocity for material interfaces:
    /// Due to the discontinuous approximation at the interface, 
    /// and the weak enforcement of the velocity jump condition `$ [[\vec{u}]] = 0 `$
    /// the velocities of both phases do not match exactly in the discrete setting.
    /// (They are equal in the continuous setting, however.)
    /// 
    /// Therefore, the phase velocities are averaged according to <see cref="XNSE_Control.InterfaceVelocityAveraging"/>
    /// to obtain a single interface velocity.
    /// </summary>
    public class LevelSetVelocity : ILevelSetParameter {
        protected int D;

        protected IList<string> parameters;

        protected int degree;

        public IList<string> ParameterNames => parameters;

        protected XNSE_Control.InterfaceVelocityAveraging averagingMode;

        protected PhysicalParameters physicalParameters;

        /// <summary>
        /// ctor.
        /// </summary>
        public LevelSetVelocity(string levelSetName, int D, int degree, XNSE_Control.InterfaceVelocityAveraging averagingMode, PhysicalParameters physicalParameters) {
            this.averagingMode = averagingMode;
            this.physicalParameters = physicalParameters;
            this.D = D;
            this.degree = degree;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
        }

        /// <summary>
        /// averaging velocity at interface
        /// </summary>
        public virtual void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(new FuncTrace()) {
                int D = levelSet.Tracker.GridDat.SpatialDimension;

                //Mean Velocity
                XDGField[] EvoVelocity; // = new XDGField[]
                try {
                    EvoVelocity = D.ForLoop(
                        d => (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d)]
                        );
                } catch {
                    Console.Error.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                    EvoVelocity = D.ForLoop(
                        d => (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0_d(d)]
                        );
                }

                DGField[] meanVelocity = new ConventionalDGField[D];

                double rho_A = physicalParameters.rho_A, rho_B = physicalParameters.rho_B;
                double mu_A = physicalParameters.mu_A, mu_B = physicalParameters.mu_B;
                LevelSetTracker lsTrkr = levelSet.Tracker;
                CellMask CC = lsTrkr.Regions.GetCutCellMask4LevSet(0);
                CellMask Neg = lsTrkr.Regions.GetLevelSetWing(0, -1).VolumeMask;
                CellMask Pos = lsTrkr.Regions.GetLevelSetWing(0, +1).VolumeMask;
                CellMask posNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Neg);
                CellMask negNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Pos);

                for(int d = 0; d < D; d++) {
                    //Basis b = EvoVelocity[d].Basis.NonX_Basis;
                    meanVelocity[d] = ParameterVarFields[ParameterNames[d]];
                    meanVelocity[d].Clear();


                    foreach(string spc in lsTrkr.SpeciesNames) {
                        double rhoSpc;
                        double muSpc;
                        switch(spc) {
                            case "A": rhoSpc = rho_A; muSpc = mu_A; break;
                            case "B": rhoSpc = rho_B; muSpc = mu_B; break;
                            case "C": continue; // solid phase; ignore
                            default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                        }

                        double scale = 1.0;
                        switch(averagingMode) {
                            case XNSE_Control.InterfaceVelocityAveraging.mean: {
                                scale = 0.5;
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.density: {
                                scale = rhoSpc / (rho_A + rho_B);
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.viscosity: {
                                scale = muSpc / (mu_A + mu_B);
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.phaseA: {
                                scale = (spc == "A") ? 1.0 : 0.0;
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.phaseB: {
                                scale = (spc == "B") ? 1.0 : 0.0;
                                break;
                            }
                        }

                        meanVelocity[d].Acc(scale, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), CC);
                        switch(spc) {
                            //case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), Neg.Except(CC)); break;
                            case "A": {
                                if(averagingMode != XNSE_Control.InterfaceVelocityAveraging.phaseB)
                                    meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), negNear);
                                break;
                            }
                            case "B": {
                                if(averagingMode != XNSE_Control.InterfaceVelocityAveraging.phaseA)
                                    meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), posNear);
                                break;
                            }
                            default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                        }
                    }
                }
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocties = new (string, DGField)[D];
            for(int d = 0; d < D; ++d) {
                Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(basis, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }


    /// <summary>
    /// Level set velocity, i.e. parameters with name <see cref="BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(string, IList{string})"/>
    /// </summary>
    public class ExplicitLevelSetVelocity :  ParameterS, ILevelSetParameter {

        string[] m_ParameterNames;

        public override IList<string> ParameterNames {
            get {
                return m_ParameterNames;
            }
        }

        // delegate (string ParameterName, DGField ParamField)[] DelParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        public override DelParameterFactory Factory => ParameterFactory;

        public ExplicitLevelSetVelocity(string levelSetName, ScalarFunctionTimeDep[] components) : base() {
            this.D = components.Length;
            m_ParameterNames = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
            m_components = components.CloneAs();
            
        }

        ScalarFunctionTimeDep[] m_components;
        int D;

        public override DelPartialParameterUpdate Update {
            get {
                return InternalParameterUpdate;
            }
        }

        void InternalParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(new FuncTrace()) {

                DGField[] meanVelocity = new ConventionalDGField[D];
                for(int d = 0; d < D; d++) {
                    //Basis b = EvoVelocity[d].Basis.NonX_Basis;
                    meanVelocity[d] = ParameterVarFields[ParameterNames[d]];
                    meanVelocity[d].Clear();
                    if(m_components[d] != null)
                        meanVelocity[d].ProjectField(m_components[d].SetTime(t));
                }
            }
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            InternalParameterUpdate(time, DomainVarFields, ParameterVarFields);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocties = new (string, DGField)[D];
            var bv = DomainVarFields[VariableNames.VelocityX].Basis;
            var b = new Basis(bv.GridDat, bv.Degree);

            for(int d = 0; d < D; ++d) {
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(b, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }
}

