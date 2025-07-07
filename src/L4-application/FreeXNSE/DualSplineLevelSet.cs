using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Statistic;
using ilPSP;
using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;


namespace FreeXNSE {
    
    
    /// <summary>
    /// Conversion of an explicit interface representation into a Level-Set 
    /// - in a 2D setting;
    /// - the internal spline represents the x-position in dependence of the y-coordinate.
    /// </summary>
    public class DualSplineLevelSet : LevelSet {
        public double[] xL { get; private set; }
        public double[] xR { get; private set; }

        public double[] yL { get; private set; }
        public double[] yR { get; private set; }


        /// <summary>
        /// the y-position of the interface in dependence of the x-coordinate
        /// </summary>
        public CubicSpline SplineL { get; private set; }

        /// <summary>
        /// the y-position of the interface in dependence of the x-coordinate
        /// </summary>
        public CubicSpline SplineR { get; private set; }

        FreeXNSE_Control ctrl;

        /// <summary>
        /// ctor
        /// </summary>
        public DualSplineLevelSet(Basis basis, string name, int numberOfNodes) : base(basis, name) {
            this.numberOfNodes = numberOfNodes;
            xL = new double[numberOfNodes];
            xR = new double[numberOfNodes];
            yL = new double[numberOfNodes];
            yR = new double[numberOfNodes];
            NodesL = MultidimensionalArray.Create(numberOfNodes, 2);
            NodesR = MultidimensionalArray.Create(numberOfNodes, 2);
            if(basis.GridDat is GridData grid) {
                BoundingBox bbox = grid.LocalBoundingBox;
                double bottom = bbox.Min[1];
                double top = bbox.Max[1];
                double increment = (top - bottom) / (numberOfNodes - 1);

                //bottom += increment / 2;
                for(int i = 0; i < numberOfNodes - 1; ++i) {
                    NodesL[i, 1] = bottom + i * increment;
                    NodesR[i, 1] = bottom + i * increment;
                }
                NodesL[numberOfNodes - 1, 1] = top;
                NodesR[numberOfNodes - 1, 1] = top;
            } else {
                throw new NotImplementedException();
            }
        }


        CubicSpline[] m_DualSplinePhi0Initial;

        /// <summary>
        /// ctor, assuming special form of initial, (y0, y1) => f0(y0) bzw.(+) f1(y1)
        /// </summary>
        public DualSplineLevelSet(FreeXNSE_Control ctrl, Basis basis, string name, int numberOfNodes) : this(basis, name, numberOfNodes) {
            this.ctrl = ctrl;
            m_DualSplinePhi0Initial = new CubicSpline[2];
            m_DualSplinePhi0Initial[0] = CubicSpline.InterpolateNaturalSorted(ctrl.DualSplinePhi0Initial[0].Item1, ctrl.DualSplinePhi0Initial[0].Item2);
            m_DualSplinePhi0Initial[1] = CubicSpline.InterpolateNaturalSorted(ctrl.DualSplinePhi0Initial[1].Item1, ctrl.DualSplinePhi0Initial[1].Item2);
            Interpolate(m_DualSplinePhi0Initial[0], true);
            Interpolate(m_DualSplinePhi0Initial[1], false);

        }

        public int numberOfNodes { get; private set; }

        public MultidimensionalArray NodesL { get; private set; }
        public MultidimensionalArray NodesR { get; private set; }


        /// <summary>
        /// Sets a node
        /// </summary>
        /// <param name="i">number of node</param>
        /// <param name="j">x or y node</param>
        /// <param name="value"></param>
        public void SetNode(int i, int j, double value, bool left)
        {
            if (left)
                NodesL[i, j] = value;
            else
                NodesR[i, j] = value;
        }

        public void Interpolate(CubicSpline initial, bool left, CellMask region = null) {
            if(left) {
                for(int i = 0; i < numberOfNodes; ++i) {
                    NodesL[i, 0] = initial.Interpolate(NodesL[i, 1]);
                }
                Interpolate(left, region);
            } else {
                for(int i = 0; i < numberOfNodes; ++i) {
                    NodesR[i, 0] = initial.Interpolate(NodesR[i, 1]);
                }
                Interpolate(left, region);
            }
        }

        public void Interpolate(MultidimensionalArray xValues, bool left, CellMask region = null) {
            if(left) {
                if(xValues.Dimension != 1 || xValues.Length != numberOfNodes) {
                    throw new NotSupportedException();
                }
                for(int i = 0; i < numberOfNodes; ++i) {
                    NodesL[i, 0] = xValues[i];
                }
                Interpolate(left, region);
            } else {
                if(xValues.Dimension != 1 || xValues.Length != numberOfNodes) {
                    throw new NotSupportedException();
                }
                for(int i = 0; i < numberOfNodes; ++i) {
                    NodesR[i, 0] = xValues[i];
                }
                Interpolate(left, region);
            }
        }

        public void Interpolate(bool left, CellMask region = null) {
            if(left) {
                for(int i = 0; i < numberOfNodes; ++i) {
                    xL[i] = NodesL[i, 0];
                    yL[i] = NodesL[i, 1];
                }
                SplineL = CubicSpline.InterpolateNaturalSorted(yL, xL);
                m_DualSplinePhi0Initial[0] = SplineL;
                ctrl.DualSplinePhi0Initial[0]= Tuple.Create(yL,xL);
            } else {
                for(int i = 0; i < numberOfNodes; ++i) {
                    xR[i] = NodesR[i, 0];
                    yR[i] = NodesR[i, 1];
                }
                SplineR = CubicSpline.InterpolateNaturalSorted(yR, xR);
                m_DualSplinePhi0Initial[1] = SplineR;
                ctrl.DualSplinePhi0Initial[1] = Tuple.Create(yR, xR);
            }
            if(SplineL != null && SplineR != null) EmbeddInLevelSet(SplineL, SplineR, this, region);
        }

        static void EmbeddInLevelSet(CubicSpline splineL, CubicSpline splineR, SinglePhaseField levelSet, CellMask region = null) {
            if(region == null) {
                region = CellMask.GetFullMask(levelSet.GridDat);
            }
            levelSet.Clear(region);
            levelSet.ProjectField(
                1.0,
                LevelSet(splineL, splineR),
                new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, region));
        }

        static ScalarFunction LevelSet(CubicSpline splineL, CubicSpline splineR) {
            return (MultidimensionalArray nodes, MultidimensionalArray results) => {
                for(int i = 0; i < nodes.Lengths[0]; ++i) {
                    //phi(x,y) = position(x) - y
                    double pL = splineL.Interpolate(nodes[i, 1]);
                    pL -= nodes[i, 0];
                    pL *= 1;
                    double pR = splineR.Interpolate(nodes[i, 1]);
                    pR -= nodes[i, 0];
                    pR *= -1;

                    results[i] = Math.Max(pL, pR);
                };
            };
        }
    }

    public class DualSplineLevelSetEvolver : ILevelSetEvolver {
        DualSplineLevelSetTimeStepper timeStepper;

        IList<string> parameters;

        string levelSetName;

        public DualSplineLevelSetEvolver(string levelSetName, GridData gridData) {
            this.levelSetName = levelSetName;
            int D = gridData.SpatialDimension;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            timeStepper = new DualSplineLevelSetTimeStepper(gridData);
        }

        public IList<string> ParameterNames => parameters;

        public IList<string> VariableNames => new string[] { };

        // nothing to do
        public Func<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>, bool> AfterMovePhaseInterface => null;

        /// <summary>
        /// <see cref="ILevelSetEvolver.InternalFields"/>; here, empty;
        /// </summary>
        public IDictionary<string, DGField> InternalFields { 
            get { return null; }
        }


        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            DualSplineLevelSet splineLs = levelSet.DGLevelSet as DualSplineLevelSet;

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

    class DualSplineLevelSetTimeStepper {
        MultidimensionalArray splineVelocityL;
        MultidimensionalArray splineVelocityR;

        FieldEvaluation evaluator;

        public DualSplineLevelSetTimeStepper(GridData gridData) {
            evaluator = new FieldEvaluation(gridData);
        }

        public void MoveLevelSet(
            double dt, DualSplineLevelSet levelSet,
            IList<SinglePhaseField> velocity,
            CellMask near) {
            {
                splineVelocityL = MultidimensionalArray.Create(levelSet.NodesL.GetLength(0), velocity.Count);
                evaluator.Evaluate(1.0, velocity, levelSet.NodesL, 0.0, splineVelocityL);

                CubicSpline spline = levelSet.SplineL;
                for(int i = 0; i < levelSet.NodesL.GetLength(0); ++i) {
                    double y = levelSet.NodesL[i, 1];
                    double f = -splineVelocityL[i, 1] * spline.Differentiate(y) + splineVelocityL[i, 0];
                    levelSet.NodesL[i, 0] += dt * f;
                }
                levelSet.Interpolate(true, near);
            }
            {
                splineVelocityR = MultidimensionalArray.Create(levelSet.NodesR.GetLength(0), velocity.Count);
                evaluator.Evaluate(1.0, velocity, levelSet.NodesR, 0.0, splineVelocityR);

                CubicSpline spline = levelSet.SplineR;
                for(int i = 0; i < levelSet.NodesR.GetLength(0); ++i) {
                    double y = levelSet.NodesR[i, 1];
                    double f = -splineVelocityR[i, 1] * spline.Differentiate(y) + splineVelocityR[i, 0];
                    levelSet.NodesR[i, 0] += dt * f;
                }
                levelSet.Interpolate(false, near);
            }
        }
    }
}
