using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;
using ilPSP;
using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;


namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    
    
    /// <summary>
    /// Conversion of an explicit interface representation into a Level-Set 
    /// - in a 2D setting;
    /// - the internal spline represents the y-position in dependence of the x-coordinate.
    /// </summary>
    public class SplineLevelSet : LevelSet {
        public double[] x { get; private set; }

        public double[] y { get; private set; }

        /// <summary>
        /// the y-position of the interface in dependence of the x-coordinate
        /// </summary>
        public CubicSpline Spline { get; private set; }

        /// <summary>
        /// ctor
        /// </summary>
        public SplineLevelSet(Basis basis, string name, int numberOfNodes) : base(basis, name) {
            this.numberOfNodes = numberOfNodes;
            x = new double[numberOfNodes];
            y = new double[numberOfNodes];
            Nodes = MultidimensionalArray.Create(numberOfNodes, 2);
            if(basis.GridDat is GridData grid) {
                BoundingBox bbox = grid.LocalBoundingBox;
                double left = bbox.Min[0];
                double right = bbox.Max[0];
                double increment = (right - left) / numberOfNodes;

                left += increment / 2;
                for(int i = 0; i < numberOfNodes; ++i) {
                    Nodes[i, 0] = left + i * increment;
                }
            } else {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        public SplineLevelSet(Func<double, double> initial, Basis basis, string name, int numberOfNodes) : this(basis, name, numberOfNodes) {
            Interpolate(initial);
        }

        public int numberOfNodes { get; private set; }

        public MultidimensionalArray Nodes { get; private set; }

        /// <summary>
        /// Sets a node
        /// </summary>
        /// <param name="i">number of node</param>
        /// <param name="j">x or y node</param>
        /// <param name="value"></param>
        public void SetNode(int i, int j, double value )
        {
            Nodes[i, j] = value;
        }

        public void Interpolate(Func<double, double> initial, CellMask region = null) {
            for(int i = 0; i < numberOfNodes; ++i) {
                Nodes[i, 1] = initial(Nodes[i, 0]);
            }
            Interpolate(region);
        }

        public void Interpolate(MultidimensionalArray yValues, CellMask region = null) {
            if(yValues.Dimension != 1 || yValues.Length != numberOfNodes) {
                throw new NotSupportedException();
            }
            for(int i = 0; i < numberOfNodes; ++i) {
                Nodes[i, 1] = yValues[i];
            }
            Interpolate(region);
        }

        public void Interpolate(CellMask region = null) {
            for(int i = 0; i < numberOfNodes; ++i) {
                x[i] = Nodes[i, 0];
                y[i] = Nodes[i, 1];
            }
            Spline = CubicSpline.InterpolateNaturalSorted(x, y);
            EmbeddInLevelSet(Spline, this, region);
        }

        static void EmbeddInLevelSet(CubicSpline spline, SinglePhaseField levelSet, CellMask region = null) {
            if(region == null) {
                region = CellMask.GetFullMask(levelSet.GridDat);
            }
            levelSet.Clear(region);
            levelSet.ProjectField(
                1.0,
                LevelSet(spline),
                new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, region));
        }

        static ScalarFunction LevelSet(CubicSpline spline) {
            return (MultidimensionalArray nodes, MultidimensionalArray results) => {
                for(int i = 0; i < nodes.Lengths[0]; ++i) {
                    //phi(x,y) = position(x) - y
                    results[i] = spline.Interpolate(nodes[i, 0]);
                    results[i] -= nodes[i, 1];
                    results[i] *= -1; // negative y is <0
                };
            };
        }
    }

    public class SplineLevelSetEvolver : ILevelSetEvolver {
        SplineLevelSetTimeStepper timeStepper;

        IList<string> parameters;

        string levelSetName;

        public SplineLevelSetEvolver(string levelSetName, GridData gridData) {
            this.levelSetName = levelSetName;
            int D = gridData.SpatialDimension;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            timeStepper = new SplineLevelSetTimeStepper(gridData);
        }

        public IList<string> ParameterNames => parameters;

        public IList<string> VariableNames => new string[] { };

        // nothing to do
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => null;

        /// <summary>
        /// <see cref="ILevelSetEvolver.InternalFields"/>; here, empty;
        /// </summary>
        public IDictionary<string, DGField> InternalFields { 
            get { return null; }
        }


        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            SplineLevelSet splineLs = levelSet.DGLevelSet as SplineLevelSet;

            SinglePhaseField[] meanVelocity = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityX)],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName,BoSSS.Solution.NSECommon.VariableNames.VelocityY)],
            };

            CellMask near = levelSet.Tracker.Regions.GetNearMask4LevSet(levelSet.LevelSetIndex, 1);

            timeStepper.MoveLevelSet(dt, splineLs, meanVelocity, near);

            CellMask posFar = levelSet.Tracker.Regions.GetLevelSetWing(levelSet.LevelSetIndex, +1).VolumeMask.Except(near);
            CellMask negFar = levelSet.Tracker.Regions.GetLevelSetWing(levelSet.LevelSetIndex, -1).VolumeMask.Except(near);
            splineLs.Clear(posFar);
            splineLs.AccConstant(1, posFar);
            splineLs.Clear(negFar);
            splineLs.AccConstant(-1, negFar);
        }
    }

    class SplineLevelSetTimeStepper {
        MultidimensionalArray splineVelocity;

        FieldEvaluation evaluator;

        public SplineLevelSetTimeStepper(GridData gridData) {
            evaluator = new FieldEvaluation(gridData);
        }

        public void MoveLevelSet(
            double dt, SplineLevelSet levelSet,
            IList<SinglePhaseField> velocity,
            CellMask near) {

            splineVelocity = MultidimensionalArray.Create(levelSet.Nodes.GetLength(0), velocity.Count);
            evaluator.Evaluate(1.0, velocity, levelSet.Nodes, 0.0, splineVelocity);

            CubicSpline spline = levelSet.Spline;
            for(int i = 0; i < levelSet.Nodes.GetLength(0); ++i) {
                double x = levelSet.Nodes[i, 0];
                double f = -splineVelocity[i, 0] * spline.Differentiate(x) + splineVelocity[i, 1];
                levelSet.Nodes[i, 1] += dt * f;
            }
            levelSet.Interpolate(near);
        }
    }
}
