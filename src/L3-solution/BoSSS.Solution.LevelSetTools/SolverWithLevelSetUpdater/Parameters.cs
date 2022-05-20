using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// Computation of normals from the level-set;
    /// - computed normals are typically **not of unit length**, i.e. the vectors must be normalized before use!
    /// - computed from broken derivatives, i.e. un-filtered
    /// </summary>
    public class Normals : ParameterS, ILevelSetParameter {

        int D;
        int degree;
        IList<string> parameterNames;

        public Normals(string LsName, int D, int degree) {
            this.D = D;
            this.degree = degree;
            if (LsName == VariableNames.LevelSetCG) {
                VariableNames.NormalVector(D);
            } else {
                parameterNames = VariableNames.AsLevelSetVariable(LsName, VariableNames.NormalVector(D));
            }
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => parameterNames;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            LevelSet Phi = levelSet.CGLevelSet;
            DGField[] Normals = new SinglePhaseField[D];
            for (int i = 0; i < D; ++i) {
                Normals[i] = ParameterVarFields[parameterNames[i]];
            }
            VectorField<DGField> normalVector = new VectorField<DGField>(Normals);
            normalVector.Clear();
            normalVector.Gradient(1.0, Phi);
        }

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis basis = new Basis(gridData, degree);
            VectorField<SinglePhaseField> Normals = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(basis, ParameterNames[d])));

            (string, DGField)[] normals = new (string, DGField)[D];
            for (int d = 0; d < D; ++d) {
                normals[d] = (parameterNames[d], Normals[d]);
            }
            return normals;
        }
    }

    public class GradientAndCurvature : ParameterS, ILevelSetParameter {

        int m_HMForder;

        string[] lsParameters;

        int gradientDegree;

        int curvatureDegree;

        public GradientAndCurvature(string LsName, int curvatureDegree, int gradientDegree, int m_HMForder, int D) : this(LsName, curvatureDegree, gradientDegree, m_HMForder, D, CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode.Curvature_Projected, CurvatureAlgorithmsForLevelSet.FilterConfiguration.NoFilter) {            
        }

        CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode m_mode;
        CurvatureAlgorithmsForLevelSet.FilterConfiguration m_filter;

        public GradientAndCurvature(string LsName, int curvatureDegree, int gradientDegree, int m_HMForder, int D, CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode mode, CurvatureAlgorithmsForLevelSet.FilterConfiguration filter) {
            if (LsName == VariableNames.LevelSetCG) {
                lsParameters = VariableNames.LevelSetGradient(D).Cat(VariableNames.Curvature);
            } else {
                lsParameters = VariableNames.AsLevelSetVariable(LsName, VariableNames.LevelSetGradient(D)).Cat(VariableNames.AsLevelSetVariable(LsName, VariableNames.Curvature));
            }
            this.m_HMForder = m_HMForder;
            this.gradientDegree = gradientDegree;
            this.curvatureDegree = curvatureDegree;

            m_mode = mode;
            m_filter = filter;
        }

        IList<string> ILevelSetParameter.ParameterNames => lsParameters;

        public override IList<string> ParameterNames => new string[] { VariableNames.Curvature };

        public override DelParameterFactory Factory => CurvatureFactory;

        (string ParameterName, DGField ParamField)[] CurvatureFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            //Curvature
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[1];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis curvatureBasis = new Basis(gridData, curvatureDegree);
            string curvatureName = VariableNames.Curvature;
            fields[0] = (curvatureName, new SinglePhaseField(curvatureBasis, curvatureName));
            return fields;
        }

        (string ParameterName, DGField ParamField)[] ILevelSetParameter.ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            int paramCount = lsParameters.Length;
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[lsParameters.Length];
            IGridData gridData = DomainVarFields.First().Value.GridDat;

            Basis basis = new Basis(gridData, gradientDegree);
            for (int i = 0; i < paramCount - 1; ++i) {
                fields[i] = (lsParameters[i], new SinglePhaseField(basis, lsParameters[i]));
            }
            Basis curvatureBasis = new Basis(gridData, curvatureDegree);
            fields[2] = (lsParameters[2], new SinglePhaseField(curvatureBasis, lsParameters[2]));
            return fields;
        }

        public void LevelSetParameterUpdate(
           DualLevelSet phaseInterface,
           double time,
           IReadOnlyDictionary<string, DGField> DomainVarFields,
           IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[lsParameters[2]];
            VectorField<SinglePhaseField> filtLevSetGradient;
            CurvatureAlgorithmsForLevelSet.CurvatureDriver(
                m_mode,
                m_filter,
                Curvature,
                out filtLevSetGradient,
                phaseInterface.Tracker,
                m_HMForder,
                phaseInterface.DGLevelSet);
            for (int i = 0; i < lsParameters.Length - 1; ++i) {
                ParameterVarFields[lsParameters[i]].Clear();
                if (filtLevSetGradient != null)
                    ParameterVarFields[lsParameters[i]].AccLaidBack(1.0, filtLevSetGradient[i]);
            }
        }
    }

}
