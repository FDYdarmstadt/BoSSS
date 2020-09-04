using System;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.IO;

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.Timestepping;
using System.Diagnostics;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTool.SemiLagrangianLevelSet;

namespace BoSSS.Application.SemiLagrangianLevelSetTestSuite
{
    public class MergingCircles : TestCase
    {
        private int CellNbr;

        private double CircleRadius;
        private double CircleIniDist;

        private Func<double [], double, double> u = (X, t) => -(Math.Abs(X[0])/X[0]);
        private Func<double [], double, double> v = (X, t) => 0;

        public MergingCircles (LagrangianMode Mode, double CircleRadius = 0.75, double CircleIniDist = 0.1, int CellNbr = 40)
            :base(0.01,30,Mode)
        {
            this.CircleRadius = CircleRadius;
            this.CircleIniDist = CircleIniDist;
            this.CellNbr = CellNbr;
            PlotInterval = 1;
        }

        public override GridCommons CreateGrid ()
        {
            double xMin = -2.5;
            double xMax = 2.5;
            int xCellNum = CellNbr*4;

            double yMin = -1.0;
            double yMax = 1.0;
            int yCellNum = CellNbr*2;

            double [] x_nodes = GenericBlas.Linspace (xMin, xMax, xCellNum - 1);
            double [] y_nodes = GenericBlas.Linspace (yMin, yMax, yCellNum - 1);
            GridCom = Grid2D.Cartesian2DGrid (x_nodes, y_nodes, type: CellType.Square_Linear);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet ()
        {
            ScalarFunctionEx MergingCircles = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                Func<double [], double> Circle = (X) => Math.Sqrt (X [0].Pow2 () + X [1].Pow2 ()) - CircleRadius;

                MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell (Ns, cell0, Len);

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    for (int i = 0; i < Ns.NoOfNodes; i++)
                    {
                        double [] X = new double [2];

                        if (GlobalNodes[c,i,0] > 0)
                        {
                            X [0] = GlobalNodes [c, i, 0] - 0.5*CircleIniDist-CircleRadius;
                            X [1] = GlobalNodes [c, i, 1];
                            result [c, i] = Circle (X);
                        }
                        else
                        {
                            X [0] = GlobalNodes [c, i, 0] + 0.5 * CircleIniDist + CircleRadius;
                            X [1] = GlobalNodes [c, i, 1];
                            result [c, i] = Circle (X);
                        }
                    }
                }
            };
            LevelSet.ProjectField (MergingCircles);
        }

        protected override void ProjectVelocity (double phystime, double dt)
        {
            Velocity_Current [0].ProjectField (X => u (X, phystime));
            Velocity_Current [1].ProjectField (X => v (X, phystime));
            Velocity_Next [0].ProjectField (X => u (X, phystime + dt));
            Velocity_Next [1].ProjectField (X => v (X, phystime + dt));
        }
    }

    public class MergingRectangles : TestCase
    {
        private int CellNbr;

        private double RecCenterDistY = 0.4;
        private double RecLen = 0.7;
        private double RecWidth = 0.5;

        public MergingRectangles (LagrangianMode Mode, int CellNbr = 100, double dt = 0.01, int NbrofTimesteps = 10)
            : base (dt, NbrofTimesteps, Mode)
        {
            this.CellNbr = CellNbr;
        }
        public override GridCommons CreateGrid ()
        {
            double xMin = -1.0;
            double xMax = 1.0;
            int xCellNum = CellNbr;

            double yMin = -1.0;
            double yMax = 1.0;
            int yCellNum = CellNbr;

            double [] x_nodes = GenericBlas.Linspace (xMin, xMax, xCellNum - 1);
            double [] y_nodes = GenericBlas.Linspace (yMin, yMax, yCellNum - 1);
            GridCom = Grid2D.Cartesian2DGrid (x_nodes, y_nodes, type: CellType.Square_Linear);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet ()
        {
            ScalarFunctionEx TwoRectangles = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                double LeftXCorner,RightXCorner;
                if (RecLen <= RecWidth)
                {
                    LeftXCorner = 0;
                    RightXCorner = 0;
                }
                else
                {
                    LeftXCorner = -(RecLen / 2) + 0.5 * RecWidth;
                    RightXCorner = (RecLen / 2) - 0.5 * RecWidth;
                }

                double UpperRecUpLeftX, UpperRecUpLeftY, UpperRecUpRightX, UpperRecUpRightY; 
                double UpperRecDownLeftX, UpperRecDownLeftY, UpperRecDownRightX,UpperRecDownRightY;

                double LowerRecUpLeftX, LowerRecUpLeftY, LowerRecUpRightX, LowerRecUpRightY;
                double LowerRecDownLeftX, LowerRecDownLeftY, LowerRecDownRightX, LowerRecDownRightY;

                double Radius = (RecLen < RecWidth) ? 0.25 * RecLen : 0.25 * RecWidth;

                UpperRecUpLeftX = UpperRecDownLeftX = LowerRecUpLeftX = LowerRecDownLeftX = -0.5 * RecLen + Radius;
                UpperRecUpRightX = UpperRecDownRightX = LowerRecUpRightX = LowerRecDownRightX = 0.5 * RecLen - Radius;

                UpperRecUpLeftY = UpperRecUpRightY = RecCenterDistY + 0.5 * RecWidth - Radius;
                UpperRecDownLeftY = UpperRecDownRightY = RecCenterDistY - 0.5 * RecWidth + Radius;

                LowerRecUpLeftY = LowerRecUpRightY = -RecCenterDistY + 0.5 * RecWidth - Radius;
                LowerRecDownLeftY = LowerRecDownRightY = -RecCenterDistY - 0.5 * RecWidth + Radius;

                MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell (Ns, cell0, Len);

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    for (int i = 0; i < Ns.NoOfNodes; i++)
                    {                        
                        if(GlobalNodes[c,i,0] < UpperRecUpLeftX && GlobalNodes[c,i,0] > UpperRecUpLeftX-2*Radius)
                        {
                            if(GlobalNodes[c,i,1] > UpperRecUpLeftY && GlobalNodes[c,i,1] < UpperRecUpLeftY+2*Radius)
                            {
                                result [c, i] = -Radius + Math.Sqrt ((GlobalNodes [c, i, 0] - UpperRecUpLeftX).Pow2 () + (GlobalNodes [c, i, 1] - UpperRecUpLeftY).Pow2 ());
                                continue;
                            }
                            if (GlobalNodes [c, i, 1] < UpperRecDownLeftY && GlobalNodes [c, i, 1] > UpperRecDownLeftY - 2 * Radius)
                            {
                                result [c, i] = -Radius + Math.Sqrt ((GlobalNodes [c, i, 0] - UpperRecDownLeftX).Pow2 () + (GlobalNodes [c, i, 1] - UpperRecDownLeftY).Pow2 ());
                                continue;
                            }
                            if (GlobalNodes[c, i, 1] > LowerRecUpLeftY && GlobalNodes[c, i, 1] < LowerRecUpLeftY + 2 * Radius)
                            {
                                result[c, i] = -Radius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerRecUpLeftX).Pow2() + (GlobalNodes[c, i, 1] - LowerRecUpLeftY).Pow2());
                                continue;
                            }
                            if (GlobalNodes[c, i, 1] < LowerRecDownLeftY && GlobalNodes[c, i, 1] > LowerRecDownLeftY - 2 * Radius)
                            {
                                result[c, i] = -Radius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerRecDownLeftX).Pow2() + (GlobalNodes[c, i, 1] - LowerRecDownLeftY).Pow2());
                                continue;
                            }
                        }
                        if (GlobalNodes [c, i, 0] > UpperRecUpRightX && GlobalNodes [c, i, 0] < UpperRecUpRightX + 2 * Radius)
                        {
                            if (GlobalNodes[c, i, 1] > UpperRecUpRightY && GlobalNodes[c, i, 1] < UpperRecUpRightY + 2 * Radius)
                            {
                                result[c, i] = -Radius + Math.Sqrt((GlobalNodes[c, i, 0] - UpperRecUpRightX).Pow2() + (GlobalNodes[c, i, 1] - UpperRecUpRightY).Pow2());
                                continue;
                            }
                            if (GlobalNodes[c, i, 1] < UpperRecDownRightY && GlobalNodes[c, i, 1] > UpperRecDownRightY - 2 * Radius)
                            {
                                result[c, i] = -Radius + Math.Sqrt((GlobalNodes[c, i, 0] - UpperRecDownRightX).Pow2() + (GlobalNodes[c, i, 1] - UpperRecDownRightY).Pow2());
                                continue;
                            }
                            if (GlobalNodes[c, i, 1] > LowerRecUpRightY && GlobalNodes[c, i, 1] < LowerRecUpRightY + 2 * Radius)
                            {
                                result[c, i] = -Radius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerRecUpRightX).Pow2() + (GlobalNodes[c, i, 1] - LowerRecUpRightY).Pow2());
                                continue;
                            }
                            if (GlobalNodes[c, i, 1] < LowerRecDownRightY && GlobalNodes[c, i, 1] > LowerRecDownRightY - 2 * Radius)
                            {
                                result[c, i] = -Radius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerRecDownRightX).Pow2() + (GlobalNodes[c, i, 1] - LowerRecDownRightY).Pow2());
                                continue;
                            }
                        }                        

                        if (GlobalNodes [c, i, 1] > 0)
                        {
                            if (GlobalNodes[c,i,0]<LeftXCorner && Math.Abs(GlobalNodes[c,i,1]-RecCenterDistY)<Math.Abs(GlobalNodes[c,i,0]-LeftXCorner))
                            {
                                result[c, i] = Math.Abs (GlobalNodes [c, i, 0] - LeftXCorner) - 0.5 * RecWidth;
                            }
                            else 
                            if(GlobalNodes[c,i,0]>RightXCorner && Math.Abs(GlobalNodes[c,i,1]-RecCenterDistY)<Math.Abs(GlobalNodes[c,i,0]-RightXCorner))
                            {
                                result[c, i] = Math.Abs (GlobalNodes [c, i, 0] - RightXCorner) - 0.5 * RecWidth;
                            }
                            else
                            {
                                result[c, i] = Math.Abs (GlobalNodes [c, i, 1] - RecCenterDistY) - 0.5 * RecWidth;
                            }
                        }                        
                        else
                        {
                            if (GlobalNodes [c, i, 0] < LeftXCorner && Math.Abs (GlobalNodes [c, i, 1] - -RecCenterDistY) < Math.Abs (GlobalNodes [c, i, 0] - LeftXCorner))
                            {
                                result [c, i] = Math.Abs (GlobalNodes [c, i, 0] - LeftXCorner) - 0.5 * RecWidth;
                            }
                            else
                            if (GlobalNodes [c, i, 0] > RightXCorner && Math.Abs (GlobalNodes [c, i, 1] - -RecCenterDistY) < Math.Abs (GlobalNodes [c, i, 0] - RightXCorner))
                            {
                                result [c, i] = Math.Abs (GlobalNodes [c, i, 0] - RightXCorner) - 0.5 * RecWidth;
                            }
                            else
                            {
                                result [c, i] = Math.Abs (GlobalNodes [c, i, 1] - -RecCenterDistY) - 0.5 * RecWidth;
                            }
                        }                        
                    }
                }
            };
            LevelSet.ProjectField (TwoRectangles);
        }

        protected override void ProjectVelocity(double phystime, double dt)
        {
            ScalarFunctionEx CurvatureVelocity_X = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                LevelSet.EvaluateTotalCurvature(cell0, Len, Ns, result);

                MultidimensionalArray[] Gradient = new MultidimensionalArray[Grid.SpatialDimension];
                for (int dim = 0; dim < Grid.SpatialDimension; dim++)
                {
                    Gradient[dim] = result.CloneAs();
                    LevelSetGradient[dim].Evaluate(cell0, Len, Ns, Gradient[dim]);
                }

                MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                double GradientNorm;
                double Gradient_Y;

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    for (int i = 0; i < Ns.NoOfNodes; i++)
                    {
                        GradientNorm = 0;
                        for (int dim = 0; dim < Grid.SpatialDimension; dim++)
                        {
                            GradientNorm += Gradient[dim][c, i].Pow2();
                        }
                        GradientNorm = Math.Sqrt(GradientNorm);
                        Gradient_Y = Gradient[0][c, i];
                        result[c, i] = result[c, i] * (Gradient[0][c, i] / GradientNorm);
                    }

                }
            };
            ScalarFunctionEx CurvatureVelocity_Y = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                LevelSet.EvaluateTotalCurvature(cell0, Len, Ns, result);

                MultidimensionalArray[] Gradient = new MultidimensionalArray[Grid.SpatialDimension];
                for (int dim = 0; dim < Grid.SpatialDimension; dim++)
                {
                    Gradient[dim] = result.CloneAs();
                    LevelSetGradient[dim].Evaluate(cell0, Len, Ns, Gradient[dim]);
                }

                MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                double GradientNorm;
                double Gradient_Y;

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    for (int i = 0; i < Ns.NoOfNodes; i++)
                    {
                        GradientNorm = 0;
                        for (int dim = 0; dim < Grid.SpatialDimension; dim++)
                        {
                            GradientNorm += Gradient[dim][c, i].Pow2();
                        }
                        GradientNorm = Math.Sqrt(GradientNorm);
                        Gradient_Y = Gradient[1][c, i];
                        result[c, i] = result[c, i] * (Gradient[1][c, i] / GradientNorm);                        
                    }

                }
            };

            Velocity_Current[0].ProjectField(CurvatureVelocity_X);
            Velocity_Current[1].ProjectField(CurvatureVelocity_Y);
            Velocity_Next = Velocity_Current;
        }
    }

    public class CircleToRotatingStar : TestCase
    {
        private int CellNbr;
        private double CircleCenterX = 0;
        private double CircleCenterY = 0;
        private double CircleRadius = 0.3;
        private Func<double [], double, double> F;
        private Func<double [], double, double> u_expansion;
        private Func<double [], double, double> v_expansion;
        private double AngularSpeed = Math.PI;
        private Func<double [], double, double> u_rotation;
        private Func<double [], double, double> v_rotation;

        public CircleToRotatingStar(LagrangianMode Mode,int CellNbr = 100, double dt = 0.01, int NbrofTimesteps = 295)
            :base(dt,NbrofTimesteps,Mode)
        {
            this.CellNbr = CellNbr;
            F = (X, t) => ((Math.Sqrt (X [0].Pow2 () + X [1].Pow2 ()) - CircleRadius).Pow2 () + 1) * (2 + Math.Sin (4 * Math.Atan (X [1] / X [0])));
            u_expansion = (X, t) => F (X, t) * Math.Cos (Math.Atan (X [1] / X [0]));
            v_expansion = (X, t) => F (X, t) * Math.Sin (Math.Atan (X [1] / X [0]));
            u_rotation = (X, t) => -Math.Sin (Math.Atan (X [1] / X [0])) * Math.Sqrt ((X [0] - CircleCenterX).Pow2 () + (X [1] - CircleCenterY).Pow2 ()) * AngularSpeed;
            v_rotation = (X, t) =>  Math.Cos (Math.Atan (X [1] / X [0])) * Math.Sqrt ((X [0] - CircleCenterX).Pow2 () + (X [1] - CircleCenterY).Pow2 ()) * AngularSpeed;
        }
        public override GridCommons CreateGrid ()
        {
            double xMin = 0.0;
            double xMax = 1.0;
            int xCellNum = CellNbr;

            double yMin = 0.0;
            double yMax = 1.0;
            int yCellNum = CellNbr;

            double [] x_nodes = GenericBlas.Linspace (xMin, xMax, xCellNum - 1);
            double [] y_nodes = GenericBlas.Linspace (yMin, yMax, yCellNum - 1);
            GridCom = Grid2D.Cartesian2DGrid (x_nodes, y_nodes, type: CellType.Square_Linear);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet ()
        {
            Func<double [], double> Phi_0 = (X) => Math.Sqrt ((X [0] - CircleCenterX).Pow2 () + (X [1] - CircleCenterY).Pow2 ()) - CircleRadius;
            LevelSet.ProjectField (X => Phi_0 (X));
        }

        protected override void ProjectVelocity (double phystime, double dt)
        {
            if (phystime < 0.95)
            {
                Velocity_Current [0].ProjectField (X => u_expansion (X, phystime));
                Velocity_Current [1].ProjectField (X => v_expansion (X, phystime));
                Velocity_Next [0].ProjectField (X => u_expansion (X, phystime + dt));
                Velocity_Next [1].ProjectField (X => v_expansion (X, phystime + dt));
            }
            else if(phystime >= 0.95)
            {
                Velocity_Current [0].ProjectField (X => u_rotation (X, phystime));
                Velocity_Current [1].ProjectField (X => v_rotation (X, phystime));
                Velocity_Next [0].ProjectField (X => u_rotation (X, phystime + dt));
                Velocity_Next [1].ProjectField (X => v_rotation (X, phystime + dt));
            }
        }
    }

    public class TheKartoffel : TestCase
    {
        private int CellNbr;

        private Func<double[], double, double> F;
        private Func<double[], double, double> u;
        private Func<double[], double, double> v;


        public TheKartoffel(LagrangianMode Mode, int CellNbr = 20, double dt = 0.01, int NbrofTimesteps = 295)
            : base(dt, NbrofTimesteps, Mode)
        {
            this.CellNbr = CellNbr;
            u = (X, t) => (- X[1]);
            v = (X, t) => (X[0]);
            PlotInterval = 1;
        }
        public override GridCommons CreateGrid()
        {
            double xMin = -1.0;
            double xMax = 1.0;
            int xCellNum = CellNbr;

            double yMin = -1.0;
            double yMax = 1.0;
            int yCellNum = CellNbr;

            double[] x_nodes = GenericBlas.Linspace(xMin, xMax, xCellNum - 1);
            double[] y_nodes = GenericBlas.Linspace(yMin, yMax, yCellNum - 1);
            GridCom = Grid2D.Cartesian2DGrid(x_nodes, y_nodes, type: CellType.Square_Linear);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet()
        {
            Func<double[], double> Phi_0 = X => (X[0].Pow2() / (0.8 * 0.8) * 1.25 + X[1].Pow2() / (0.8 * 0.8) * 0.8) - 1.0 + 0.2 * Math.Sin(10 * X[0] * X[1]);
            LevelSet.ProjectField(X => Phi_0(X));
        }

        protected override void ProjectVelocity(double phystime, double dt)
        {
            Velocity_Current[0].ProjectField(X => u(X, phystime));
            Velocity_Current[1].ProjectField(X => v(X, phystime));
            Velocity_Next[0].ProjectField(X => u(X, phystime + dt));
            Velocity_Next[1].ProjectField(X => v(X, phystime + dt));
        }
    }

    public class ZalesaksDisk : TestCase
    {
        private readonly int CellNbr;
        private Func<double [], double, double> u;
        private Func<double [], double, double> v;

        private double DiskCenterX = 50;
        private double DiskCenterY = 75;
        private double DiskRadius = 15;
        private double SlotWidth = 5;
        private double SlotLenfromCenter = 5;

        public ZalesaksDisk(LagrangianMode Mode, int CellNbr = 100, double dt = 1, int NbrofTimesteps = 1256)
            : base(dt, NbrofTimesteps, Mode)
        {
            this.CellNbr = CellNbr;
            u = (X, t) => (Math.PI / 314) * (DiskCenterX - X [1]);
            v = (X, t) => (Math.PI / 314) * (X [0] - DiskCenterX);
            PlotInterval = 10;
        }

        public override GridCommons CreateGrid()
        {
            double xMin = 0.0;
            double xMax = 100.0;
            int xCellNum = CellNbr;

            double yMin = 0.0;
            double yMax = 100.0;
            int yCellNum = CellNbr;

            double[] x_nodes = GenericBlas.Linspace(xMin, xMax, xCellNum - 1);
            double[] y_nodes = GenericBlas.Linspace(yMin, yMax, yCellNum - 1);
            GridCom = Grid2D.Cartesian2DGrid(x_nodes, y_nodes, type: CellType.Square_Linear, periodicX: true, periodicY: true);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet()
        {
            ScalarFunctionEx ZalesaksDisk = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                Func<double [], double> Circle = (X) => Math.Sqrt ((X [0] - DiskCenterX).Pow2 () + (X [1] - DiskCenterY).Pow2 ()) - DiskRadius;
                double SlotTopYCorner = DiskCenterY + SlotLenfromCenter - 0.5 * SlotWidth;
                double InnerCornerLeftX, InnerCornerLeftY, InnerCornerRightX, InnerCornerRightY;
                double LowerCornerLeftX, LowerCornerLeftY, LowerCornerRightX, LowerCornerRightY;

                double InnerCornerRadius, LowerCornerRadius;
                InnerCornerRadius = 0.25 * SlotWidth;
                LowerCornerRadius = 0.25 * SlotWidth;

                InnerCornerLeftX = DiskCenterX - 0.5*SlotWidth + InnerCornerRadius;
                InnerCornerRightX = DiskCenterX + 0.5*SlotWidth - InnerCornerRadius;
                InnerCornerLeftY = InnerCornerRightY = DiskCenterY + SlotLenfromCenter - InnerCornerRadius;
                
                LowerCornerLeftX = DiskCenterX - 0.5*SlotWidth - LowerCornerRadius;
                LowerCornerRightX = DiskCenterX + 0.5*SlotWidth + LowerCornerRadius;
                LowerCornerLeftY = LowerCornerRightY = -Math.Sqrt((DiskRadius - LowerCornerRadius).Pow2() - (LowerCornerRightX-DiskCenterX).Pow2())+DiskCenterY;

                MultidimensionalArray VectorToLowerCornerLeftPoint = MultidimensionalArray.Create(1, 2);
                VectorToLowerCornerLeftPoint[0, 0] = LowerCornerLeftX - DiskCenterX;
                VectorToLowerCornerLeftPoint[0, 1] = LowerCornerLeftY - DiskCenterY;
                MultidimensionalArray VectorToLowerCornerRightPoint = MultidimensionalArray.Create(1, 2);
                VectorToLowerCornerRightPoint[0, 0] = LowerCornerRightX - DiskCenterX;
                VectorToLowerCornerRightPoint[0, 1] = LowerCornerRightY - DiskCenterY;
                MultidimensionalArray VectorMinusY = MultidimensionalArray.Create(1, 2);
                VectorMinusY[0, 0] = 0;
                VectorMinusY[0, 1] = -1;

                double LowerCornerLeftAngle = SemiLagrangianLevelSetMethod<MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.VectorAngle(VectorMinusY, VectorToLowerCornerLeftPoint);
                double LowerCornerRightAngle = SemiLagrangianLevelSetMethod<MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.VectorAngle(VectorMinusY, VectorToLowerCornerRightPoint);

                MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell (Ns, cell0, Len);

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    for (int i = 0; i < Ns.NoOfNodes; i++)
                    {
                        double[] X = new double[2];
                        X[0] = GlobalNodes[c, i, 0];
                        X[1] = GlobalNodes[c, i, 1];

                        if (GlobalNodes[c, i, 1] > InnerCornerLeftY && GlobalNodes[c, i, 1] < InnerCornerLeftY + 2*InnerCornerRadius)
                        {
                            if (GlobalNodes[c, i, 0] < InnerCornerLeftX && GlobalNodes[c, i, 0] > InnerCornerLeftX - 2*InnerCornerRadius)
                            {
                                result[c, i] = InnerCornerRadius - Math.Sqrt((GlobalNodes[c, i, 0] - InnerCornerLeftX).Pow2() + (GlobalNodes[c, i, 1] - InnerCornerLeftY).Pow2());
                                continue;
                            }
                            if (GlobalNodes[c, i, 0] > InnerCornerRightX && GlobalNodes[c, i, 0] < InnerCornerRightX + 2*InnerCornerRadius)
                            {
                                result[c, i] = InnerCornerRadius - Math.Sqrt((GlobalNodes[c, i, 0] - InnerCornerRightX).Pow2() + (GlobalNodes[c, i, 1] - InnerCornerRightY).Pow2());
                                continue;
                            }
                        }
                        if (GlobalNodes[c, i, 1] < LowerCornerLeftY && GlobalNodes[c, i, 1] > LowerCornerLeftY - 0.5*(SlotWidth/LowerCornerRadius)*LowerCornerRadius)
                        {
                            if (GlobalNodes[c, i, 0] > LowerCornerLeftX && GlobalNodes[c, i, 0] < LowerCornerLeftX + 0.5*(SlotWidth / LowerCornerRadius) * LowerCornerRadius)
                            {
                                result[c, i] = -LowerCornerRadius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerCornerLeftX).Pow2() + (GlobalNodes[c, i, 1] - LowerCornerLeftY).Pow2());
                                continue;
                            }                            
                            if (GlobalNodes[c, i, 0] < LowerCornerRightX && GlobalNodes[c, i, 0] > LowerCornerRightX - 0.5*(SlotWidth / LowerCornerRadius) * LowerCornerRadius)
                            {
                                result[c, i] = -LowerCornerRadius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerCornerRightX).Pow2() + (GlobalNodes[c, i, 1] - LowerCornerRightY).Pow2());
                                continue;
                            }
                            MultidimensionalArray VectorToGlobalNode = MultidimensionalArray.Create(1, 2);
                            VectorToGlobalNode[0, 0] = GlobalNodes[c,i,0] - DiskCenterX;
                            VectorToGlobalNode[0, 1] = GlobalNodes[c,i,1] - DiskCenterY;
                            double GlobalNodeAngle = SemiLagrangianLevelSetMethod<MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.VectorAngle(VectorMinusY, VectorToGlobalNode);

                            if (GlobalNodes[c, i, 0] < LowerCornerLeftX && GlobalNodeAngle <= LowerCornerLeftAngle)
                            {
                                result[c, i] = -LowerCornerRadius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerCornerLeftX).Pow2() + (GlobalNodes[c, i, 1] - LowerCornerLeftY).Pow2());
                                continue;
                            }
                            if (GlobalNodes[c, i, 0] > LowerCornerRightX && GlobalNodeAngle <= LowerCornerRightAngle)
                            {
                                result[c, i] = -LowerCornerRadius + Math.Sqrt((GlobalNodes[c, i, 0] - LowerCornerRightX).Pow2() + (GlobalNodes[c, i, 1] - LowerCornerRightY).Pow2());
                                continue;
                            }

                        }
                        if (Math.Sqrt((GlobalNodes[c, i, 0] - DiskCenterX).Pow2() + (GlobalNodes[c, i, 1] - DiskCenterY).Pow2()) > DiskRadius)
                        {
                            double RecRes = 0.5 * SlotWidth - Math.Abs(GlobalNodes[c, i, 0] - DiskCenterX);
                            double CircRes = Circle(X);
                            if (Math.Abs(GlobalNodes[c, i, 0] - DiskCenterX) <= SlotWidth / 2 && GlobalNodes[c, i, 1] < DiskCenterY)
                                result[c, i] = (Math.Abs(RecRes) < Math.Abs(CircRes)) ? CircRes : RecRes;
                            else
                                result[c, i] = CircRes;
                        }
                        else
                        {
                            if (Math.Abs(GlobalNodes[c, i, 0] - DiskCenterX) < (SlotWidth / 2) && GlobalNodes[c, i, 1] < (DiskCenterY + SlotLenfromCenter))
                            {
                                if (GlobalNodes[c, i, 1] > SlotTopYCorner && Math.Abs(GlobalNodes[c, i, 0] - DiskCenterX) <= Math.Abs(GlobalNodes[c, i, 1] - SlotTopYCorner))
                                {
                                    result[c, i] = (SlotWidth / 2) - Math.Abs(GlobalNodes[c, i, 1] - SlotTopYCorner);
                                }
                                else
                                {
                                    result[c, i] = 0.5 * SlotWidth - Math.Abs(GlobalNodes[c, i, 0] - DiskCenterX);
                                }
                            }
                            else
                            {
                                double RecRes;
                                double Circ = Circle(X);
                                if (GlobalNodes[c, i, 1] > SlotTopYCorner && Math.Abs(GlobalNodes[c, i, 0] - DiskCenterX) <= Math.Abs(GlobalNodes[c, i, 1] - SlotTopYCorner))
                                    RecRes = (SlotWidth / 2) - Math.Abs(GlobalNodes[c, i, 1] - SlotTopYCorner);
                                else
                                {
                                    RecRes = 0.5 * SlotWidth - Math.Abs(GlobalNodes[c, i, 0] - DiskCenterX);
                                }
                                result[c, i] = (Math.Abs(RecRes) < Math.Abs(Circ)) ? RecRes : Circ;
                            }
                        }
                    }
                }
            };
            LevelSet.ProjectField(ZalesaksDisk);
        }

        protected override void ProjectVelocity(double phystime, double dt)
        {
            Velocity_Current[0].ProjectField(X => u(X, phystime));
            Velocity_Current[1].ProjectField(X => v(X, phystime));
            Velocity_Next[0].ProjectField(X => u(X, phystime + dt));
            Velocity_Next[1].ProjectField(X => v(X, phystime + dt));
        }
    }

    public class ZalesaksSphere : TestCase
    {
        private readonly int CellNbr;
        private Func<double[], double, double> u = (X, t) => (Math.PI / 314) * (50 - X[1]);
        private Func<double[], double, double> v = (X, t) => (Math.PI / 314) * (X[0] - 50);
        private Func<double[], double, double> w = (X, t) => 0;

        private double DiskCenterX;
        private double DiskCenterY;
        private double DiskCenterZ;
        private double DiskRadius;
        private double SlotWidth;
        private double SlotLenfromCenter;

        public ZalesaksSphere(LagrangianMode Mode, int CellNbr = 100, double dt = 1, int NbrofTimesteps = 1000)
            : base(dt, NbrofTimesteps, Mode)
        {
            this.CellNbr = CellNbr;
        }

        public override GridCommons CreateGrid()
        {
            double xMin = 0.0;
            double xMax = 100.0;
            int xCellNum = CellNbr;

            double yMin = 0.0;
            double yMax = 100.0;
            int yCellNum = CellNbr;

            double zMin = 0.0;
            double zMax = 100.0;
            int zCellNum = CellNbr;

            double[] x_nodes = GenericBlas.Linspace(xMin, xMax, xCellNum - 1);
            double[] y_nodes = GenericBlas.Linspace(yMin, yMax, yCellNum - 1);
            double[] z_nodes = GenericBlas.Linspace(zMin, zMax, zCellNum - 1);
            GridCom = Grid3D.Cartesian3DGrid(x_nodes, y_nodes, z_nodes, periodicX: false, periodicY: false,
                periodicZ: false, _CellType: CellType.Cube_Linear);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet()
        {
            ScalarFunctionEx ZalesaksSphere = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                Func<double [], double> Circle = (X) => Math.Sqrt ((X [0] - DiskCenterX).Pow2 () + (X [1] - DiskCenterY).Pow2 () + (X [2] - DiskCenterZ).Pow2 ()) - DiskRadius;
                double SlotTopYCorner = DiskCenterY + SlotLenfromCenter - 0.5 * SlotWidth;

                MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell (Ns, cell0, Len);

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    for (int i = 0; i < Ns.NoOfNodes; i++)
                    {
                        double [] X = new double [3];
                        X [0] = GlobalNodes [c, i, 0];
                        X [1] = GlobalNodes [c, i, 1];
                        X [2] = GlobalNodes [c, i, 2];
                        if (Math.Sqrt ((GlobalNodes [c, i, 0] - DiskCenterX).Pow2 () + (GlobalNodes [c, i, 1] - DiskCenterY).Pow2 ()) > DiskRadius)
                        {
                            if (Math.Abs (GlobalNodes [c, i, 0] - DiskCenterX) <= SlotWidth / 2)
                            {
                                double RecRes = 0.5 * SlotWidth - Math.Abs (GlobalNodes [c, i, 0] - DiskCenterX);
                                double CircRes = Circle (X);
                                result [c, i] = (Math.Abs (RecRes) < Math.Abs (CircRes)) ? CircRes : RecRes;
                            }
                            else
                            {
                                double RecRes;
                                if (GlobalNodes [c, i, 1] > SlotTopYCorner && Math.Abs (GlobalNodes [c, i, 0] - DiskCenterX) <= Math.Abs (GlobalNodes [c, i, 1] - SlotTopYCorner))
                                    RecRes = (SlotWidth / 2) - Math.Abs (GlobalNodes [c, i, 1] - SlotTopYCorner);
                                else
                                    RecRes = (SlotWidth / 2) - Math.Abs (GlobalNodes [c, i, 1] - DiskCenterY);

                                double CircRes = Circle (X);

                                if (Math.Abs (GlobalNodes [c, i, 0] - DiskCenterX) < SlotWidth / 2 && GlobalNodes [c, i, 1] < (DiskCenterY + SlotLenfromCenter))
                                    result [c, i] = RecRes;
                                else
                                    result [c, i] = (Math.Abs (RecRes) < Math.Abs (CircRes)) ? RecRes : CircRes;
                            }
                        }
                    }
                }
            };
        }

        protected override void ProjectVelocity(double phystime, double dt)
        {
            Velocity_Current[0].ProjectField(X => u(X, phystime));
            Velocity_Current[1].ProjectField(X => v(X, phystime));
            Velocity_Current[2].ProjectField(X => w(X, phystime));
            Velocity_Next[0].ProjectField(X => u(X, phystime + dt));
            Velocity_Next[1].ProjectField(X => v(X, phystime + dt));
            Velocity_Next[2].ProjectField(X => w(X, phystime));
        }
    }

    public class CircleInSingleVortex : TestCase
    {
        private readonly int CellNbr;
        /*
        * Streamfunction
        * SIGMA = (1/Math.PI)*(Math.Sin(Math.PI*X[0]).Pow2())*(Math.Sin(Math.PI*X[1]).Pow2())
        */
        private Func<double[], double, double> u = (X, t) => 2 * Math.Cos(Math.PI * X[1]) * Math.Sin(Math.PI * X[1]) * (Math.Sin(Math.PI * X[0])).Pow2();
        private Func<double[], double, double> v = (X, t) =>-2 * Math.Cos(Math.PI * X[0]) * Math.Sin(Math.PI * X[0]) * (Math.Sin(Math.PI * X[1])).Pow2();

        public CircleInSingleVortex(LagrangianMode Mode, int CellNbr = 400, double dt = 0.0025, int NbrofTimesteps=1000)
            :base(dt,NbrofTimesteps,Mode)
        {
            this.CellNbr = CellNbr;
            PlotInterval = 25;
        }

        public override GridCommons CreateGrid()
        {
            double xMin = 0.0;
            double xMax = 1.0;
            int xCellNum = CellNbr;

            double yMin = 0.0;
            double yMax = 1.0;
            int yCellNum = CellNbr;

            double[] x_nodes = GenericBlas.Linspace(xMin, xMax, xCellNum-1);
            double[] y_nodes = GenericBlas.Linspace(yMin, yMax, yCellNum-1);
            GridCom = Grid2D.Cartesian2DGrid(x_nodes, y_nodes, type: CellType.Square_Linear);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet()
        {
            double CircleCenterX = 0.5;
            double CircleCenterY = 0.75;
            double CircleRadius = 0.15;
            Func<double[], double> Phi_0 = (X) => Math.Sqrt((X[0] - CircleCenterX).Pow2() + (X[1] - CircleCenterY).Pow2()) - CircleRadius;
            LevelSet.ProjectField(X => Phi_0(X));
        }

        protected override void ProjectVelocity(double phystime, double dt)
        {
            if (phystime < 1.25)
            {
                Velocity_Current[0].ProjectField(X => u(X, phystime));
                Velocity_Current[1].ProjectField(X => v(X, phystime));
                Velocity_Next[0].ProjectField(X => u(X, phystime + dt));
                Velocity_Next[1].ProjectField(X => v(X, phystime + dt));
            }
            else if(phystime >= 1.25)
            {
                Velocity_Current[0].ProjectField(X => -u(X, phystime));
                Velocity_Current[1].ProjectField(X => -v(X, phystime));
                Velocity_Next[0].ProjectField(X => -u(X, phystime + dt));
                Velocity_Next[1].ProjectField(X => -v(X, phystime + dt));
            }
        }
    }

    public class MovingSphere : TestCase
    {
        private Func<double[], double, double> u = (X, t) => 1;
        private Func<double[], double, double> v = (X, t) => 0;
        private Func<double[], double, double> w = (X, t) => 0;

        public MovingSphere(double dt,int NbrofTimesteps,LagrangianMode Mode)
            :base(dt,NbrofTimesteps,Mode)
        { }

        public override GridCommons CreateGrid()
        {
            double xMin = -2.0;
            double xMax = 8.0;
            int xCellNum = 64;

            double yMin = -1.0;
            double yMax = 1.0;
            int yCellNum = 16;

            double zMin = -1.0;
            double zMax = 1.0;
            int zCellNum = 16;

            double[] x_nodes = GenericBlas.Linspace(xMin, xMax, xCellNum);
            double[] y_nodes = GenericBlas.Linspace(yMin, yMax, yCellNum);
            double[] z_nodes = GenericBlas.Linspace(zMin, zMax, zCellNum);
            GridCom = Grid3D.Cartesian3DGrid(x_nodes, y_nodes, z_nodes, periodicX: false, periodicY: false,
                periodicZ: false,_CellType:CellType.Cube_Linear);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet()
        {
            double CircleCenterX = 0;
            double CircleCenterY = 0;
            double CircleCenterZ = 0;
            double CircleRadius = 0.7;
            Func<double[], double> Phi_0 = (X) => Math.Sqrt((X[0] - CircleCenterX).Pow2() + (X[1] - CircleCenterY).Pow2() + (X[2] - CircleCenterZ).Pow2()) - CircleRadius;
            LevelSet.ProjectField(X => Phi_0(X));
        }

        protected override void ProjectVelocity(double phystime, double dt)
        {
            Velocity_Current[0].ProjectField(X => u(X, phystime));
            Velocity_Current[1].ProjectField(X => v(X, phystime));
            Velocity_Current[2].ProjectField(X => w(X, phystime));
            Velocity_Next[0].ProjectField(X => u(X, phystime + dt));
            Velocity_Next[1].ProjectField(X => v(X, phystime + dt));
            Velocity_Next[2].ProjectField(X => w(X, phystime));
        }
    }

    public class MovingCircle : TestCase
    {
        public MovingCircle(LagrangianMode Mode, double dt = 0.1, int NbrofTimesteps = 1)
            : base(dt, NbrofTimesteps, Mode)
        {
            PlotInterval = 1;
        }

        public override GridCommons CreateGrid()
        {
            double xMin = -2.0;
            double xMax = 8.0;
            int xCellNum = 50;

            double yMin = -1.5;
            double yMax = 1.5;
            int yCellNum = 15;

            double[] x_nodes = GenericBlas.Linspace(xMin, xMax, xCellNum);
            double[] y_nodes = GenericBlas.Linspace(yMin, yMax, yCellNum);
            GridCom = Grid2D.Cartesian2DGrid(x_nodes, y_nodes, type: CellType.Square_Linear, periodicX: false, periodicY: false);
            return GridCom;
        }

        protected override void ProjectInitialLevelSet()
        {
            double CircleCenterX = 0;
            double CircleCenterY = 0;
            double CircleRadius = 1;
            Func<double[], double> Phi_0 = (X) => Math.Sqrt((X[0] - CircleCenterX).Pow2() + (X[1] - CircleCenterY).Pow2()) - CircleRadius;
            LevelSet.ProjectField(X => Phi_0(X));
        }

        protected override void ProjectVelocity(double phystime, double dt)
        {
            Func<double[], double, double> u = (X, t) => 0;
            Func<double[], double, double> v = (X, t) => 0;
            Velocity_Current[0].ProjectField(X => u(X, phystime));
            Velocity_Current[1].ProjectField(X => v(X, phystime));
            Velocity_Next[0].ProjectField(X => u(X, phystime + dt));
            Velocity_Next[1].ProjectField(X => v(X, phystime + dt));
        }
    }    

    public abstract class TestCase
    {
        protected LevelSet LevelSet;
        protected VectorField<SinglePhaseField> LevelSetGradient;
        protected LevelSetTracker LevelSetTrck;
        protected VectorField<SinglePhaseField> Velocity_Next;
        protected VectorField<SinglePhaseField> Velocity_Current;
        protected IGridData Grid;
        protected GridCommons GridCom;
        protected DGField[] PlotFields;
        protected Tecplot plt;
        protected int dgDegree = 3;
        protected int PlotInterval = 10;

        protected Stopwatch sw;
        public readonly double dt;
        public readonly int NbrofTimesteps;
        private readonly LagrangianCorrectors Corrector;

        public TestCase(double dt, int NbrofTimesteps, LagrangianMode Mode)
        {
            this.dt = dt;
            this.NbrofTimesteps = NbrofTimesteps;
            sw = new Stopwatch();
            Corrector = new LagrangianCorrectors(Mode);
        }

        public abstract GridCommons CreateGrid();

        public virtual void CreateFields(IGridData Grid)
        {
            this.Grid = Grid;
            LevelSet = new LevelSet(new Basis(this.Grid, dgDegree), "LevelSet");
            LevelSetTrck = new LevelSetTracker((GridData)this.Grid, XQuadFactoryHelper.MomentFittingVariants.Saye, 4, new string[] { "A", "B" }, LevelSet);
            LevelSetGradient = new VectorField<SinglePhaseField>(this.Grid.SpatialDimension, new Basis(this.Grid, dgDegree), "LevelSetGradient", SinglePhaseField.Factory);
            Velocity_Current = new VectorField<SinglePhaseField>(this.Grid.SpatialDimension, new Basis(this.Grid, dgDegree), "VelocityCurrentStep", SinglePhaseField.Factory);
            Velocity_Next = new VectorField<SinglePhaseField>(this.Grid.SpatialDimension, new Basis(this.Grid, dgDegree), "VelocityNextStep", SinglePhaseField.Factory);
            plt = new Tecplot(this.Grid, 3);
            PlotFields = new DGField[2 * this.Grid.SpatialDimension + 1];
            PlotFields[0] = LevelSet;
            for (int dim = 0; dim < this.Grid.SpatialDimension; dim++)
            {
                PlotFields[dim + 1] = Velocity_Current[dim];
                PlotFields[dim + 1 + this.Grid.SpatialDimension] = LevelSetGradient[dim];
            }
            Corrector.Constructor(LevelSet,LevelSetGradient,LevelSetTrck,Velocity_Next,Velocity_Current,Grid);
        }

        protected abstract void ProjectInitialLevelSet();

        protected abstract void ProjectVelocity(double phystime, double dt);

        public virtual void InitialWork()
        {
            ProjectInitialLevelSet();
            Console.WriteLine("Projected");
            LevelSetGradient.Gradient(1.0, LevelSet);
            LevelSetTrck.UpdateTracker();
            ProjectVelocity(0.0, 0.0);
            Corrector.Initialize();
        }

        public virtual void TimestepWork(double phystime, double dt, int Timestep)
        {
            ProjectVelocity(phystime, dt);
            LevelSetTrck.UpdateTracker();
            LevelSetGradient.Gradient(1.0, LevelSet);
            Corrector.Timestep(dt, 1, Timestep);
        }

        public virtual void Plot(TimestepNumber timestepNo, double physTime)
        {
            if (timestepNo.MajorNumber % PlotInterval == 0)
            {
                sw.Start();
                plt.PlotFields("LevelSet_" + timestepNo, physTime, PlotFields);
                sw.Stop();
                Console.WriteLine("Print .plt Elapsed={0}", sw.Elapsed);
                sw.Reset();

                sw.Start();
                Corrector.PlotPoints(timestepNo);
                sw.Stop();
                Console.WriteLine("Print .csv Elapsed={0}", sw.Elapsed);
                sw.Reset();
            }
        }
    }

    public enum LagrangianMode { Particle, Marker, FrontTracking };

    public static class Curvature
    {
        public static VectorField<DGField> CurvatureVelocity(VectorField<SinglePhaseField> LevelSetGradient)
        {
            IGridData Grid = LevelSetGradient.GridDat;
            Basis Basis = LevelSetGradient[0].Basis;
            SinglePhaseField Curvature = new SinglePhaseField(Basis, "Curvature");
            VectorField<SinglePhaseField> LevelSetLaplacian = new VectorField<SinglePhaseField>(Grid.SpatialDimension, Basis, "Intermediate", SinglePhaseField.Factory);

            for (int dim = 0; dim < Grid.SpatialDimension; dim++)
            {
                SinglePhaseField LevelSetGradientDIMClone = LevelSetGradient[dim].CloneAs();
                ScalarFunctionEx GradientNorm = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
                {
                    LevelSetGradientDIMClone.Evaluate(cell0, Len, Ns, result);
                    MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                    int c = 0;
                    for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                    {
                        for (int i = 0; i < Ns.NoOfNodes; i++)
                        {
                            result[c, i] = result[c, i] / Math.Abs(result[c, i]);
                        }
                    }
                };
                LevelSetGradient[dim].ProjectField(GradientNorm);
                VectorField<SinglePhaseField> LevelSetLaplacianInter = new VectorField<SinglePhaseField>(Grid.SpatialDimension, Basis, "Intermediate2", SinglePhaseField.Factory);
                LevelSetLaplacianInter.Gradient(1.0, LevelSetGradient[dim]);
                LevelSetLaplacian[dim] = LevelSetLaplacianInter[dim];
            }

            ScalarFunctionEx CalcCurvature = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                MultidimensionalArray[] result_dim = new MultidimensionalArray[Grid.SpatialDimension];
                for (int dim = 0; dim < Grid.SpatialDimension; dim++)
                {
                    result_dim[dim] = result.CloneAs();
                    LevelSetLaplacian[dim].Evaluate(cell0, Len, Ns, result_dim[dim]);
                }
                MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    for (int i = 0; i < Ns.NoOfNodes; i++)
                    {
                        result[c, i] = 1;
                        for (int dim = 0; dim < Grid.SpatialDimension; dim++)
                        {
                            result[c, i] *= result_dim[dim][c, i];
                        }
                    }
                }
            };

            Curvature.ProjectField(CalcCurvature);
            VectorField<DGField> CurvatureVelocity = new VectorField<DGField>(Grid.SpatialDimension, Basis, "CurvatureVelocity",SinglePhaseField.Factory);

            for (int dim=0;dim<Grid.SpatialDimension;dim++)
            {
                ScalarFunctionEx CurvatureVelo_dim = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
                {
                    MultidimensionalArray result_Curv = result.CloneAs();
                    MultidimensionalArray result_normal = result.CloneAs();
                    Curvature.Evaluate(cell0, Len, Ns, result_Curv);
                    LevelSetGradient[dim].Evaluate(cell0, Len, Ns, result_normal);

                    MultidimensionalArray GlobalNodes = Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                    int c = 0;
                    for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                    {
                        for (int i = 0; i < Ns.NoOfNodes; i++)
                        {
                            result[c, i] = result_Curv[c,i]*result_normal[c,i];
                        }
                    }
                };
                CurvatureVelocity[dim].ProjectField(CalcCurvature);
            }

            return CurvatureVelocity;
        }
    }

    public class LagrangianCorrectors
    {
        private FrontTrackingLevelSet2D FrontTracking = null;
        private MarkerLevelSet Marker = null;
        private ParticleLevelSet Particle = null;

        private readonly LagrangianMode Mode;

        public LagrangianCorrectors(LagrangianMode Mode)
        {
            this.Mode = Mode;
        }
        public void Constructor(SinglePhaseField LevelSet, VectorField<SinglePhaseField> LevelSetGradient,
            LevelSetTracker LevelSetTrck, VectorField<SinglePhaseField> Velocity_Next, VectorField<SinglePhaseField> Velocity_Current,
            IGridData GridData)
        {
            switch (Mode)
            {
                case LagrangianMode.Particle:
                    Particle = new ParticleLevelSet(Velocity_Next, Velocity_Current, LevelSet, LevelSetTrck, LevelSetGradient,
                        40, 45, 35, 0.2, 0.001, 0.2, 0.001, 10, 2, GridData, true, 101,
                        ParticleLevelSet.MinimalDistanceSearchMode.FullSearch,
                        ParticleLevelSet.TopologyMergingMode.Off,ParticleLevelSet.NormalVectorDampingMode.Off);
                    break;
                case LagrangianMode.Marker:
                    Marker = new MarkerLevelSet(Velocity_Next, Velocity_Current, LevelSet, LevelSetTrck, LevelSetGradient,
                        //6,8,5,10, 1, GridData, 50, MarkerLevelSet.MinimalDistanceSearchMode.FullSearch,
                        15, 20, 10, 10, 1, GridData, 50, MarkerLevelSet.MinimalDistanceSearchMode.FullSearch,
                        MarkerLevelSet.MarkerCorrectionMode.MinDistanceWithNormal, MarkerLevelSet.TopologyMergingMode.On,
                        MarkerLevelSet.NormalVectorDampingMode.MaxAngle45_avg);
                    break;
                case LagrangianMode.FrontTracking:
                    FrontTracking = new FrontTrackingLevelSet2D(Velocity_Next, Velocity_Current, LevelSet, LevelSetTrck,
                        LevelSetGradient, 10, 6, GridData, true, 10,
                        FrontTrackingLevelSet2D.MinimalDistanceSearchMode.FullSearch, FrontTrackingLevelSet2D.TopologyMergingMode.Off,
                        FrontTrackingLevelSet2D.NormalVectorDampingMode.Off, 0.1);
                    break;
            }
        }

        public void Initialize()
        {
            switch (Mode)
            {
                case LagrangianMode.Particle:
                    Particle.Initialize();
                    break;
                case LagrangianMode.Marker:
                    Marker.Initialize();
                    break;
                case LagrangianMode.FrontTracking:
                    FrontTracking.Initialize();
                    break;
            }
        }

        public void Timestep(double dt, int substeps, int Timestep)
        {
            switch (Mode)
            {
                case LagrangianMode.Particle:
                    Particle.PerformTimestep(dt, substeps, Timestep);
                    break;
                case LagrangianMode.Marker:
                    Marker.PerformTimestep(dt, substeps, Timestep);                
                    break;
                case LagrangianMode.FrontTracking:
                    FrontTracking.PerformTimestep(dt, substeps, Timestep);
                    break;
            }
        }

        public void PlotPoints(TimestepNumber timestepNo)
        {
            switch (Mode)
            {
                case LagrangianMode.Particle:
                    Particle.PrintCSV("Particle_" + timestepNo);
                    break;
                case LagrangianMode.Marker:
                    Marker.PrintCSV("Marker_" + timestepNo);
                    break;
                case LagrangianMode.FrontTracking:
                    FrontTracking.PrintCSV("FrontTracking_" + timestepNo);
                    break;
            }
        }
    }
}