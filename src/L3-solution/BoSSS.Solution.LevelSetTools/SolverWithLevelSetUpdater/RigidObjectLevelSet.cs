using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    public class RigidObjectLevelSet : LevelSet {

        /// <summary>
        /// ctor
        /// </summary>
        public RigidObjectLevelSet(Func<double[], double, double>[] ObjectsLevelSet, double GridLengthParameter, CellMask Region, Basis basis, string name) : base(basis, name) {
            this.ObjectsLevelSet = ObjectsLevelSet;
            this.GridLengthParameter = GridLengthParameter;
            SetLevelSet(LevelSetFunction, this, Region);
        }

        private readonly Func<double[], double, double>[] ObjectsLevelSet;
        private readonly double GridLengthParameter;

        double LevelSetFunction(double[] X, double t) {
            double levelSetFunction = int.MinValue;
            for(int i = 0; i < ObjectsLevelSet.Length; i++) {
                if (levelSetFunction < ObjectsLevelSet[i](X, GridLengthParameter))
                    levelSetFunction = ObjectsLevelSet[i](X, GridLengthParameter);
            }
            return levelSetFunction;
        }

        static void SetLevelSet(Func<double[], double, double> ObjectsLevelSet, SinglePhaseField levelSet, CellMask Region = null) {
            if (Region == null) {
                Region = CellMask.GetFullMask(levelSet.GridDat);
            }
            ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(ObjectsLevelSet, 0.0);
            levelSet.Clear(Region);
            levelSet.ProjectField(1.0, Function, new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, Region));
        }
    }

    /// <summary>
    /// A <see cref="ILevelSetEvolver"/>-Driver for multiple rigid objects, i.e. the level set evolution depends on center of mass and orientation angle of the objects.
    /// </summary>
    /// <remarks>
    /// implemented by B. Deußen, Jan. 2021.
    /// </remarks>
    public class RigidObjectLevelSetEvolver : ILevelSetEvolver {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="LevSetFunction"></param>
        public RigidObjectLevelSetEvolver(Func<double[], double, double>[] ObjectsLevelSet, double GridLengthParameter) {
            this.ObjectsLevelSet = ObjectsLevelSet;
            this.GridLengthParameter = GridLengthParameter;
        }

        private readonly Func<double[], double, double>[] ObjectsLevelSet;
        private readonly double GridLengthParameter;

        double LevelSetFunction(double[] X, double t) {
            double levelSetFunction = int.MinValue;
            for (int i = 0; i < ObjectsLevelSet.Length; i++) {
                if (levelSetFunction < ObjectsLevelSet[i](X, GridLengthParameter))
                    levelSetFunction = ObjectsLevelSet[i](X, GridLengthParameter);
            }
            return levelSetFunction;
        }

        public IList<string> ParameterNames => null;

        public IList<string> VariableNames => null;

        // nothing to do
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => null;

        public IDictionary<string, DGField> InternalFields => null;

        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var ls = levelSet.DGLevelSet;
            ls.Clear();
            ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(LevelSetFunction, time);
            ls.ProjectField(1.0, Function, new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, CellMask.GetFullMask(ls.GridDat)));
        }
    }
}
