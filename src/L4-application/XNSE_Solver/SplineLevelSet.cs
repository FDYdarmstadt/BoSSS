using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using MathNet.Numerics.Interpolation;

namespace XNSE_Solver
{
    class SplineLevelSet : LevelSet
    {
        CubicSpline spline;

        public SplineLevelSet(Basis basis, string name) : base(basis, name)
        {

        }

        void Embedd(Func<Vector, double> intial, int numberOfNodes)
        {
            double[] x = new double[numberOfNodes];
            double[] y = new double[numberOfNodes];

            spline = CubicSpline.InterpolateBoundariesSorted();
        }

        void Embedd(CubicSpline spline)
        {

        }

        CubicSpline GetSpline()
        {
            return spline;
        }
    }

    class SplineLevelSetEvolver : ILevelSetEvolver
    {
        int curvatureDegree;

        int m_HMForder;

        SplineLevelSetTimeStepper timeStepper;

        string[] parameters;

        string[] variables;

        string levelSetName;

        public SplineLevelSetEvolver(string levelSetName, int D)
        {
            this.levelSetName = levelSetName;
            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            variables = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
        }

        public IList<string> ParameterNames => parameters;

        public IList<string> VariableNames => variables;

        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            SplineLevelSet splineLs = levelSet.DGLevelSet as SplineLevelSet;

            SinglePhaseField[] meanVelocity = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityX)],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName,BoSSS.Solution.NSECommon.VariableNames.VelocityY)],
            };
            timeStepper.MoveSpline(splineLs, meanVelocity);
        }
    }

    class SplineLevelSetTimeStepper
    {
        public void MoveSpline(SplineLevelSet levelSet, IList<SinglePhaseField> velocity)
        {

        }
    }
}
