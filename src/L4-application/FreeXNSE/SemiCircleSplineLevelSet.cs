using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Statistic;
using ilPSP;
using ilPSP.Utils;
using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;


namespace FreeXNSE {


    /// <summary>
    /// Conversion of an explicit interface representation into a Level-Set 
    /// - in a 2D setting;
    /// - the internal spline represents the x- and y-positions over an angular parametric coordinate.
    /// </summary>
    public class SemiCircleSplineLevelSet : LevelSet {
        public double[] x { get; private set; }
        public double[] y { get; private set; }
        public double[] a { get; private set; }

        public double[] center { get; private set; } = new double[2];

        /// <summary>
        /// the x-position of the interface in dependence of the a-coordinate
        /// </summary>
        public CubicSpline SplineX { get; private set; }

        /// <summary>
        /// the y-position of the interface in dependence of the a-coordinate
        /// </summary>
        public CubicSpline SplineY { get; private set; }

        FreeXNSE_Control ctrl;

        /// <summary>
        /// ctor
        /// </summary>
        public SemiCircleSplineLevelSet(Basis basis, string name, int numberOfNodes) : base(basis, name) {
            this.numberOfNodes = numberOfNodes;
            a = new double[numberOfNodes];
            x = new double[numberOfNodes];
            y = new double[numberOfNodes];
            Nodes = MultidimensionalArray.Create(numberOfNodes, 2); // x, y

            if(basis.GridDat is GridData grid) {
                double bottom = 0.0;
                double top = Math.PI;
                double increment = (top - bottom) / (numberOfNodes - 3);
                double eps = increment / 1;
                a[0] = bottom;
                a[1] = bottom + eps;
                for(int i = 2; i < numberOfNodes - 1; ++i) {
                    a[i] = bottom + i * increment;
                }
                a[numberOfNodes - 2] = top - eps;
                a[numberOfNodes - 1] = top;
            } else {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// ctor, assuming special form of initial, (y0, y1) => f0(y0) bzw.(+) f1(y1)
        /// </summary>
        public SemiCircleSplineLevelSet(FreeXNSE_Control ctrl, Basis basis, string name, int numberOfNodes) : this(basis, name, numberOfNodes) {
            this.ctrl = ctrl;
            Interpolate(ctrl.SemiCircleSplinePhi0Initial);
        }

        public int numberOfNodes { get; private set; }

        public MultidimensionalArray Nodes { get; private set; }

        /// <summary>
        /// Sets a node
        /// </summary>
        /// <param name="i">number of node</param>
        /// <param name="j">a or x or y node</param>
        /// <param name="value"></param>
        public void SetNode(int i, int j, double value, bool left) {
            Nodes[i, j] = value;
        }

        public void Interpolate(Func<double, double[]> initial, CellMask region = null) {
            for(int i = 0; i < numberOfNodes; ++i) {
                double[] point = initial(a[i]);
                Nodes[i, 0] = point[0];
                Nodes[i, 1] = point[1];
            }
            Interpolate(region);
        }

        public void Interpolate(MultidimensionalArray xyValues, CellMask region = null) {
            if(xyValues.Dimension != 2 || xyValues.Length != 2 * numberOfNodes) {
                throw new NotSupportedException();
            }
            for(int i = 0; i < numberOfNodes; ++i) {
                Nodes[i, 0] = xyValues[i, 0];
                Nodes[i, 1] = xyValues[i, 1];
            }
            Interpolate(region);
        }

        public void Interpolate(CellMask region = null) {
            for(int i = 0; i < numberOfNodes; ++i) {
                x[i] = Nodes[i, 0];
                y[i] = Nodes[i, 1];
            }
            SplineX = CubicSpline.InterpolateNaturalSorted(a, x);
            SplineY = CubicSpline.InterpolateNaturalSorted(a, y);

            double r0 = new double[] { x[0] - center[0], y[0] - center[1] }.L2Norm();
            double rE = new double[] { x[numberOfNodes - 1] - center[0], y[numberOfNodes - 1] - center[1] }.L2Norm();

            double imbalance = 1.5;
            if(r0 / rE > imbalance || rE / r0 > imbalance) {
                Console.Write("attempting to move centerpoint ...");
                try {
                    UpdateCenter();
                    SplineX = CubicSpline.InterpolateNaturalSorted(a, x);
                    SplineY = CubicSpline.InterpolateNaturalSorted(a, y);
                    Console.WriteLine("success!");
                    Console.WriteLine("new center : ({0}|{1})", center[0], center[1]);
                } catch(MathNet.Numerics.NonConvergenceException e) {
                    Console.WriteLine("unable to move centerpoint, continuing");
                }
            }
            ctrl.SemiCircleSplinePhi0Initial = a => new double[] { SplineX.Interpolate(a), SplineY.Interpolate(a) };
            //Console.WriteLine("({0}|{1})", center[0], center[1]);
            //Console.WriteLine("({0}|{1})",x[0], y[0]);
            //Console.WriteLine("({0}|{1})", x[numberOfNodes - 1], y[numberOfNodes - 1]);
            EmbeddInLevelSet(SplineX, SplineY, this, region);
        }

        void UpdateCenter() {
            double[] newcenter = new double[2];
            newcenter[0] = 0.5 * (x[0] + x[numberOfNodes - 1]);
            newcenter[1] = 0.5 * (y[0] + y[numberOfNodes - 1]);

            // calculate new Nodes
            for(int i = 0; i < numberOfNodes; ++i) {
                double oldangle = MathNet.Numerics.RootFinding.RobustNewtonRaphson.FindRoot(b => Math.Atan2(SplineY.Interpolate(b) - newcenter[1] ,SplineX.Interpolate(b) - newcenter[0]) - a[i], b => (SplineX.Differentiate(b) * (SplineY.Interpolate(b) - newcenter[0]) - SplineY.Differentiate(b) * (SplineY.Interpolate(b) - newcenter[0])) / (Math.Pow(SplineX.Interpolate(b) - newcenter[0], 2) + Math.Pow(SplineY.Interpolate(b) - newcenter[1], 2)), 0, Math.PI);
                //Console.WriteLine("New:{0}, old:{1}", a[i], oldangle);
                Nodes[i, 0] = SplineX.Interpolate(oldangle); 
                Nodes[i, 1] = SplineY.Interpolate(oldangle);
            }

            // reset x and y
            for(int i = 0; i < numberOfNodes; ++i) {
                x[i] = Nodes[i, 0];
                y[i] = Nodes[i, 1];
            }

            center = newcenter;
        }

        void EmbeddInLevelSet(CubicSpline splineL, CubicSpline splineR, SinglePhaseField levelSet, CellMask region = null) {
            if(region == null) {
                region = CellMask.GetFullMask(levelSet.GridDat);
            }
            levelSet.Clear(region);
            levelSet.ProjectField(
                1.0,
                LevelSet(splineL, splineR),
                new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, region));
        }

        ScalarFunction LevelSet(CubicSpline splineX, CubicSpline splineY) {
            return (MultidimensionalArray nodes, MultidimensionalArray results) => {
                for(int i = 0; i < nodes.Lengths[0]; ++i) {
                    double[] point = new double[2];
                    point[0] = nodes[i, 0] - center[0];
                    point[1] = nodes[i, 1] - center[1];

                    double r = point.L2Norm();
                    double angle = Math.Atan(point[1] / point[0]);
                    if(point[0] < 0.0) angle = angle + Math.PI;

                    double[] point0 = new double[] { splineX.Interpolate(angle) - center[0], splineY.Interpolate(angle) - center[1] };
                    double r0 = point0.L2Norm();
                    results[i] = r - r0;
                };
            };
        }

        static int No = 0;
        public void WriteToText() {
            if(No == 0) {
                using(var file = new StreamWriter("SemicircleLevelSet-" + No++ + ".txt")) {
                    for(int i = 0; i < numberOfNodes; ++i) {
                        file.WriteLine($"{i}\t{Nodes[i, 0]}\t{Nodes[i, 1]}");
                    }
                }
            }
            using(var file = new StreamWriter("SemicircleLevelSet-" + No++ + ".txt")) {
                for(int i = 0; i < numberOfNodes; ++i) {
                    file.WriteLine($"{i}\t{Nodes[i, 0]}\t{Nodes[i, 1]}");
                }
            }            
        }
    }

    public class SemiCircleSplineLevelSetEvolver : ILevelSetEvolver {
        SemiCircleSplineLevelSetTimeStepper timeStepper;

        IList<string> parameters;

        string levelSetName;

        public SemiCircleSplineLevelSetEvolver(string levelSetName, GridData gridData) {
            this.levelSetName = levelSetName;
            int D = gridData.SpatialDimension;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            timeStepper = new SemiCircleSplineLevelSetTimeStepper(gridData);
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
            SemiCircleSplineLevelSet splineLs = levelSet.DGLevelSet as SemiCircleSplineLevelSet;

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

            splineLs.WriteToText();
        }
    }

    class SemiCircleSplineLevelSetTimeStepper {
        MultidimensionalArray splineVelocity;

        FieldEvaluation evaluator;

        public SemiCircleSplineLevelSetTimeStepper(GridData gridData) {
            evaluator = new FieldEvaluation(gridData);
        }

        public void MoveLevelSet(
            double dt, SemiCircleSplineLevelSet levelSet,
            IList<SinglePhaseField> velocity,
            CellMask near) {
            
            splineVelocity = MultidimensionalArray.Create(levelSet.Nodes.GetLength(0), velocity.Count);
            evaluator.Evaluate(1.0, velocity, levelSet.Nodes, 0.0, splineVelocity);

            CubicSpline splineX = levelSet.SplineX;
            CubicSpline splineY = levelSet.SplineY;

            for(int i = 0; i < levelSet.Nodes.GetLength(0); ++i) {
                double x = levelSet.Nodes[i, 0] - levelSet.center[0];
                double y = levelSet.Nodes[i, 1] - levelSet.center[0];
                double r = new[] { x, y }.L2Norm();
                double angle = levelSet.a[i];

                double ur = splineVelocity[i,0] * Math.Cos(angle) + splineVelocity[i,1] * Math.Sin(angle);
                double ua = -splineVelocity[i, 0] * Math.Sin(angle) + splineVelocity[i, 1] * Math.Cos(angle);

                double dr = dt * (ur - ua / r * (Math.Cos(angle) * splineX.Differentiate(angle) + Math.Sin(angle) * splineY.Differentiate(angle)));

                levelSet.Nodes[i, 0] += dr * Math.Cos(angle);
                levelSet.Nodes[i, 1] += dr * Math.Sin(angle);               

            }
            levelSet.Interpolate(near);
                        
        }
    }
}
