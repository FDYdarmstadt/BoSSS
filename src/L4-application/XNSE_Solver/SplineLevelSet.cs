using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;
using ilPSP;
using MathNet.Numerics.Interpolation;

namespace BoSSS.Application.XNSE_Solver
{
    class SplineLevelSet : LevelSet
    {
        double[] x;

        double[] y;

        public CubicSpline Spline { get; private set; }

        public SplineLevelSet(Basis basis, string name, int numberOfNodes) : base(basis, name)
        {
            this.numberOfNodes = numberOfNodes;
            x = new double[numberOfNodes];
            y = new double[numberOfNodes];
            Nodes = MultidimensionalArray.Create(2, numberOfNodes);
            if (basis.GridDat is GridData grid) 
            {
                BoundingBox bbox = grid.LocalBoundingBox;
                double left = bbox.Min[0];
                double right = bbox.Max[0];
                double increment = (right - left) / numberOfNodes;

                left += increment / 2;
                for(int i = 0; i < numberOfNodes; ++i)
                {
                    Nodes[0, i] = left + i * increment;
                }
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        public SplineLevelSet(Func<double, double> initial, Basis basis, string name, int numberOfNodes) : this(basis, name, numberOfNodes)
        {
            Interpolate(initial);
        }

        public int numberOfNodes { get; private set; }

        public MultidimensionalArray Nodes { get; private set; }

        public void Interpolate(Func<double, double> initial, CellMask region = null)
        {
            for(int i = 0; i < numberOfNodes; ++i)
            {
                Nodes[1, i] = initial(Nodes[0, i]); 
            }
            Interpolate(region);
        }

        public void Interpolate(MultidimensionalArray yValues, CellMask region = null)
        {
            if( yValues.Dimension != 1 || yValues.Length != numberOfNodes)
            {
                throw new NotSupportedException();
            }
            for(int i = 0; i < numberOfNodes; ++i)
            {
                Nodes[1,i] = yValues[i];
            }
            Interpolate(region);
        }

        public void Interpolate(CellMask region = null)
        {
            for (int i = 0; i < numberOfNodes; ++i)
            {
                x[i] = Nodes[0, i];
                y[i] = Nodes[1, i];
            }
            Spline = CubicSpline.InterpolateNaturalSorted(x, y);
            EmbeddInLevelSet(Spline, this, region);
        }

        static void EmbeddInLevelSet(CubicSpline spline, SinglePhaseField levelSet, CellMask region = null)
        {
            if(region == null)
            {
                region = CellMask.GetFullMask(levelSet.GridDat);
            }
            levelSet.Clear(region);
            levelSet.ProjectField(
                1.0,
                LevelSet(spline),
                new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, region));
        }

        static ScalarFunction LevelSet(CubicSpline spline)
        {
            return (MultidimensionalArray nodes, MultidimensionalArray results) =>
            {
                for (int i = 0; i < nodes.Lengths[0]; ++i)
                {
                    //phi(x,y) = position(x) - y
                    results[i] = spline.Interpolate(nodes[i, 0]);
                    results[i] -= nodes[i, 1];
                };
            };
        }
    }

    class SplineLevelSetEvolver : ILevelSetEvolver
    {
        SplineLevelSetTimeStepper timeStepper;

        IList<string> parameters;

        string levelSetName;

        public SplineLevelSetEvolver(string levelSetName, GridData gridData)
        {
            this.levelSetName = levelSetName;
            int D = gridData.SpatialDimension;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            timeStepper = new SplineLevelSetTimeStepper(gridData);
        }

        public IList<string> ParameterNames => parameters;

        public IList<string> VariableNames => new string[] { };

        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
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

    class SplineLevelSetTimeStepper
    {
        MultidimensionalArray splineVelocity;

        FieldEvaluation evaluator; 

        public SplineLevelSetTimeStepper(GridData gridData)
        {
            evaluator = new FieldEvaluation(gridData);
        }

        public void MoveLevelSet(
            double dt, SplineLevelSet levelSet, 
            IList<SinglePhaseField> velocity,
            CellMask near)
        {
            evaluator.Evaluate(1.0, velocity, levelSet.Nodes, 1.0, splineVelocity);

            CubicSpline spline = levelSet.Spline;
            for (int i = 0; i < levelSet.Nodes.GetLength(1); ++i)
            {
                double x = levelSet.Nodes[0, i];
                double f = - splineVelocity[0, i] * spline.Differentiate(x) + splineVelocity[1, i];
                levelSet.Nodes[1, i] += dt * f;
            }
            levelSet.Interpolate(near);
        }
    }
}
