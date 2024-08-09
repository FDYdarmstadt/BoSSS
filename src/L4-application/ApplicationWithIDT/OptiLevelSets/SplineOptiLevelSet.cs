using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using ilPSP;
using System;
using System.Linq;
using BoSSS.Solution.Statistic;
using System.Collections.Generic;
using ilPSP.LinSolvers;
using MathNet.Numerics.Interpolation;
using ilPSP.Utils;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using System.IO;

namespace ApplicationWithIDT.OptiLevelSets {
    /// <summary>
    /// Spline Level Set for 2D
    /// Level Set defined by an explicit spline representation of the interface, i.e. phi(x,y)=x-S(y), where S:R->R is th spline
    /// Conversion of an explicit interface representation into a Level-Set 
    /// - in a 2D setting;
    /// - the internal spline represents the x-position in dependence of the y-coordinate.
    /// </summary>
    public class SplineOptiLevelSet : LevelSet, IOptiLevelSet {
        public double[] x { get; private set; }

        public double[] y { get; private set; }

        public Mode m_Mode;
        /// <summary>
        /// the y-position of the interface in dependence of the x-coordinate
        /// </summary>
        public IInterpolation Spline { get; private set; }

        /// <summary>
        /// made up of Node information and 1st derivative information
        /// </summary>
        public MultidimensionalArray m_AllParams;
        public LevelSetTracker m_thisTracker;
        public LevelSet m_LevelSet;
        public int m_LevelSetNumber;
        public double xMax;
        public double xMin;
        public LevelSetTracker m_Tracker;
        /// <summary>
        /// ctor
        /// </summary>
        public SplineOptiLevelSet(Basis basis, LevelSet LsTBO, LevelSetTracker Tracker, int LevelSetNumber) : base(basis, "splineLevelSet") {

            //we extract (hopefully) all possible y values for the given grid
            if(GridDat is GridData grid) {
                var Distinct_y_values = grid.Vertices.Coordinates.ExtractSubArrayShallow(-1, 1).CloneAs().To1DArray().ToList().Distinct().ToList();
                Distinct_y_values.Sort();
                // in practice the given values will feature double points, so next they are sorted out
                for(int i = 0; i < Distinct_y_values.Count - 1; i++) {
                    //remove a node if it is the same as the next
                    if(Distinct_y_values[i + 1] - Distinct_y_values[i] < 1e-14) {
                        Distinct_y_values.RemoveAt(i);
                        i = i - 1;
                    }
                }
                numberOfNodes = Distinct_y_values.Count;
                var Distinct_x_values = grid.Vertices.Coordinates.ExtractSubArrayShallow(-1, 0).CloneAs().To1DArray();
                xMax = Distinct_x_values.Max();
                xMin = Distinct_x_values.Min();
                //dependent on the degree we have more free Parameters
                if(Basis.Degree == 0) {
                    throw new NotSupportedException("degree 0 not allowed");
                } else if(Basis.Degree == 1) {
                    m_AllParams = MultidimensionalArray.Create(Distinct_y_values.Count);
                    m_Mode = Mode.Deg1;
                } else if(Basis.Degree == 2) {
                    m_AllParams = MultidimensionalArray.Create(Distinct_y_values.Count * 2); //+ (int)(((double)Distinct_y_values.Count) / 2.0), 1);
                    m_Mode = Mode.Deg2;
                } else {
                    m_AllParams = MultidimensionalArray.Create(numberOfNodes * 2);
                    m_Mode = Mode.Deg3;
                }

                y = Distinct_y_values.ToArray();
            } else {
                throw new NotImplementedException();
            }
            m_LevelSet = LsTBO;
            m_Tracker = Tracker;
            m_LevelSetNumber = LevelSetNumber;
        }

        int numberOfNodes;


        public void Interpolate(CellMask region = null) {
            GetSpline();
            EmbeddInLevelSet(Spline, this, region);
        }
        /// <summary>
        /// changes the Member Spline with respect to the member Points x, y
        /// </summary>
        public void GetSpline() {
            x = m_AllParams.ExtractSubArrayShallow(new int[] { 0 }, new int[] { numberOfNodes - 1 }).To1DArray();
            switch(m_Mode) {
                case Mode.Deg1:
                Spline = LinearSpline.InterpolateSorted(y, x);
                break;
                case Mode.Deg2:
                // here one should normally calculate c3 from c2 and x
                // We could also make it C2- Continuous which would give us one free parameter alongside x
                double[] c2 = m_AllParams.ExtractSubArrayShallow(new int[] { numberOfNodes }, new int[] { numberOfNodes * 2 - 1 }).To1DArray();

                double[] c3 = new double[numberOfNodes];
                Spline = CubicSpline.InterpolateHermiteSorted(y, x, c2);
                break;
                case Mode.Deg3:
                c2 = m_AllParams.ExtractSubArrayShallow(new int[] { numberOfNodes }, new int[] { numberOfNodes * 2 - 1 }).To1DArray();

                Spline = CubicSpline.InterpolateHermiteSorted(y, x, c2);
                break;
            }
        }
        /// <summary>
        /// changes the member Spline with respect to some given points
        /// </summary>
        /// <param name="xPoints"></param>
        /// <param name="yPoints"></param>
        /// <param name="cPoints"></param>
        public void GetSpline(double[] xPoints, double[] yPoints, double[] cPoints) {
            switch(m_Mode) {
                case Mode.Deg1:
                Spline = LinearSpline.Interpolate(yPoints, xPoints);
                break;
                case Mode.Deg2:
                Spline = CubicSpline.InterpolateHermite(yPoints, xPoints, cPoints);
                break;
                case Mode.Deg3:
                Spline = CubicSpline.InterpolateHermite(yPoints, xPoints, cPoints);
                break;
            }
        }

        /// <summary>
        /// Get the spline from a big chunk of points, using that spline we also determine x
        /// points(:,0) are the x - Values
        /// points(:,1) are the y - Points
        /// points(:,2) are the first derivatives c <- optional, here the lengths must be the same as the SPline has yPoints
        /// </summary>
        /// <param name="points"></param>
        public void GetSplineOverDetermined(MultidimensionalArray points) {
            //GetSplineOverDeterminedV2(points);
            GetSplineOverDeterminedOld(points);

            GetSpline();
        }

        public void GetSplineOverDeterminedV2(MultidimensionalArray points) {
            #region new version
            #region Mapping CellToPoints
            /// the purpose if this part is to obtain mapping between the y-Cells and the Points (e.g Cell 0 (= [ y[0],y[1] ]) -> (0,1,2,3) (= [yPoins[0],yPoins[1],yPoins[2],yPoins[3]]) )
            Dictionary<int, List<int>> CellToPoints = new Dictionary<int, List<int>>();
            for(int j = 0; j < y.Length - 1; j++) {
                CellToPoints.Add(j, new List<int>());
            }
            var xPoints = points.ExtractSubArrayShallow(-1, 0).To1DArray();
            var yPoints = points.ExtractSubArrayShallow(-1, 1).To1DArray();
            for(int i = 0; i < yPoints.Length; i++) {
                for(int j = 0; j < y.Length - 1; j++) {
                    if(yPoints[i] > y[j] && yPoints[i] < y[j + 1]) {
                        CellToPoints.TryGetValue(j, out List<int> list);
                        list.Add(i);
                    }
                }
            }
            #endregion

            #region Find the indices of the two Minimum and Maximum y-Points (so far unused)
            //find two maximum y-Points and two Minimum y-Points (so {yMax = Max(y-Points),Max(y-Points \ yMax ) }, {yMin = Min(y-Points),Min(y-Points \ yMin ) }
            //int[] twoBiggestY;
            //int[] twoSmallestY;

            //if(yPoints[0] > yPoints[1]) {
            //    twoBiggestY = new int[] { 0, 1 };
            //    twoSmallestY = new int[] { 1, 0 };
            //} else {
            //    twoBiggestY = new int[] { 1, 0 };
            //    twoSmallestY = new int[] { 0, 1 };
            //}
            //double yMax = yPoints[twoBiggestY[0]];
            //double yMin = yPoints[twoSmallestY[0]];
            //for(int i = 1; i < yPoints.Length; i++) {
            //    yMax = yPoints[twoBiggestY[0]];
            //    yMin = yPoints[twoSmallestY[0]];
            //    if(yPoints[i] > yMax) {
            //        twoBiggestY[1] = twoBiggestY[0];
            //        twoBiggestY[0] = i;
            //    }
            //    if(yMin > yPoints[i]) {
            //        twoSmallestY[1] = twoSmallestY[0];
            //        twoSmallestY[0] = i;
            //    }
            //}
            //yMax = yPoints[twoBiggestY[0]];
            //yMin = yPoints[twoSmallestY[0]];



            #endregion

            #region PointsInCell (unused)
            //here we want to obtain a mapping between the Spline y points and a boolean indicating whether the point is in the range of the given y-Points, 
            // yMin < y[i] < yMax -> true
            //      else          -> false
            // by this we know in which cells the Spline is outside of the domain
            //Dictionary<int, bool> PointsInCell = new Dictionary<int, bool>();
            //for(int i = 0; i < y.Length; i++) {
            //    if(y[i] < yMin) {
            //        PointsInCell.Add(i, false);
            //    } else if(y[i] > yMax) {
            //        PointsInCell.Add(i, false);
            //    } else {
            //        PointsInCell.Add(i, true);
            //    }
            //}
            #endregion



            if(points.Lengths[1] == 2) {
                #region determine cPoints (1st Derivative)
                //lastly we want to determine the first derivatives
                var cPoints = new double[yPoints.Length];
                {
                    //forward FD on first point
                    var eps = (yPoints[1] - yPoints[0]);
                    var diff = xPoints[1] - xPoints[0];
                    cPoints[0] = diff / eps;
                    //backward FD on the rest
                    for(int i = 1; i < cPoints.Length; i++) {
                        eps = (yPoints[i] - yPoints[i - 1]);
                        diff = xPoints[i] - xPoints[i - 1];
                        cPoints[i] = diff / eps;
                    }
                }
                #endregion
                #region determine cPoints (2nd Derivative) (so far unused)
                //var cPoints = new double[yPoints.Length];
                //{
                //    //right/forward FD on first point
                //    var eps = ((yPoints[1] - yPoints[0]) + (yPoints[2] - yPoints[1])) / 2;
                //    var diff = xPoints[2] - 2 * xPoints[1] + xPoints[0];
                //    cPoints[0] = diff / (eps * eps);
                //    //don't include first and last point, use central FDs in the middle-points
                //    for(int i = 1; i < cPoints.Length - 1; i++) {
                //        eps = ((yPoints[i] - yPoints[i - 1]) + (yPoints[i + 1] - yPoints[i])) / 2;
                //        diff = xPoints[i + 1] - 2 * xPoints[i] + xPoints[i - 1];
                //        cPoints[i] = diff / (eps * eps);
                //    }
                //    //do backward FDS on the last one
                //    var iEnd = yPoints.Length - 1;
                //    eps = ((yPoints[iEnd - 1] - yPoints[iEnd - 2]) + (yPoints[iEnd] - yPoints[iEnd - 1])) / 2;
                //    diff = xPoints[iEnd] - 2 * xPoints[iEnd - 1] + xPoints[iEnd - 2];
                //    cPoints[iEnd] = diff / (eps * eps);
                //}
                #endregion 
                GetSpline(xPoints, yPoints, cPoints);
            } else if(points.Lengths[1] == 3) {
                var cPoints = points.ExtractSubArrayShallow(-1, 2).To1DArray();
                GetSpline(xPoints, yPoints, cPoints);
            } else {
                throw new ArgumentException($"array has {points.Lengths[1]} columns but must have 2 or 3");
            }
            //after the Spline is constructed we can get the values for x,y,c
            this.x = new double[y.Length];


            bool CheckifPointHasNeighbours(int iY) {
                bool exbiggerpoint = false;
                bool exsmallerpoint = false;
                //Checks if the point lies outside, so if it has a smaller and bigger neighbor
                for(int iYPoint = 0; iYPoint < yPoints.Length; iYPoint++) {
                    if(yPoints[iYPoint] - y[iY] < 1e-06) {
                        exsmallerpoint = true;
                    } else if(y[iY] - yPoints[iYPoint] < 1e-06) {
                        exbiggerpoint = true;
                    }
                }

                //if(exbiggerpoint && exsmallerpoint) {

                //    for(int iYPoint = 0; iYPoint < yPoints.Length; iYPoint++) {
                //        if(Math.Abs(yPoints[iYPoint] - y[iY]) < 1e-01) {
                //            return true;
                //        }
                //     }
                //}
                if(exbiggerpoint && exsmallerpoint) {
                    return true;
                } else {
                    return false;
                }
            }



            //loop over all interpolatory points y, check if point has a close neighbor in yPoints, if so use the spline, if not use the linear Extrapolation Func
            List<int> CellsLeft = new List<int>();
            List<int> CellsWithEnoughPoints = new List<int>();
            for(int iY = 0; iY < y.Length; iY++) {
                if(CheckifPointHasNeighbours(iY)) {
                    this.x[iY] = this.Spline.Interpolate(y[iY]);
                    m_AllParams[iY] = this.Spline.Interpolate(y[iY]);
                    m_AllParams[iY + y.Length] = this.Spline.Differentiate(y[iY]);
                    CellsWithEnoughPoints.Add(iY);
                } else {
                    CellsLeft.Add(iY);
                }
                //else {
                //    this.x[iY] = ExtrapolationFunc(y[iY]);
                //    m_AllParams[iY] = ExtrapolationFunc(y[iY]);
                //    m_AllParams[iY + y.Length] = ExtrapolationFuncGetA(y[iY]);
                //}
            }
            //for the remaining y's we do an linear interpolation using the nearest neighbor cell with points 
            int[] twoBiggestY = new int[] { CellsWithEnoughPoints[CellsWithEnoughPoints.Count - 2], CellsWithEnoughPoints.Last() };
            int[] twoSmallestY = new int[] { CellsWithEnoughPoints[0], CellsWithEnoughPoints[1] };
            double yMin = y[CellsWithEnoughPoints[0]];
            double yMax = y[CellsWithEnoughPoints.Last()];

            double ExtrapolationFuncGetA(double t) {
                if(t > (yMin + yMax) / 2) {
                    var x0 = x[twoBiggestY[0]];
                    var x1 = x[twoBiggestY[1]];
                    var y0 = y[twoBiggestY[0]];
                    var y1 = y[twoBiggestY[1]];
                    double a = (x0 - x1) / (y0 - y1);
                    return a;
                } else {
                    var x0 = x[twoSmallestY[0]];
                    var x1 = x[twoSmallestY[1]];
                    var y0 = y[twoSmallestY[0]];
                    var y1 = y[twoSmallestY[1]];
                    double a = (x0 - x1) / (y0 - y1);
                    return a;
                }
            };
            double ExtrapolationFunc(double t) {
                if(t > (yMin + yMax) / 2) {
                    var x0 = x[twoBiggestY[0]];
                    var x1 = x[twoBiggestY[1]];
                    var y0 = y[twoBiggestY[0]];
                    var y1 = y[twoBiggestY[1]];
                    double a = (x0 - x1) / (y0 - y1);
                    double b = x0 - a * y0;
                    return a * t + b;
                } else {
                    var x0 = x[twoSmallestY[0]];
                    var x1 = x[twoSmallestY[1]];
                    var y0 = y[twoSmallestY[0]];
                    var y1 = y[twoSmallestY[1]];
                    double a = (x0 - x1) / (y0 - y1);
                    double b = x0 - a * y0;
                    return a * t + b;
                }

            }


            foreach(int iY in CellsLeft) {

                this.x[iY] = ExtrapolationFunc(y[iY]);
                m_AllParams[iY] = ExtrapolationFunc(y[iY]);
                m_AllParams[iY + y.Length] = ExtrapolationFuncGetA(y[iY]);
            }



            //foreach(KeyValuePair<int, bool> pair in PointsInCell) {
            //    if(pair.Value) { // if there are y-Points in the Cell [y[i],y[i+1]) 
            //        CellToPoints.TryGetValue(pair.Key, out List<int> list);
            //        bool is_regular_cell = false; 
            //        for(int i = 0; i < list.Count; i++) { //we call a cell regular if there are "enough" close points to the left boundary [y[i],y[i+1]) 
            //            is_regular_cell = is_regular_cell || Math.Abs(yPoints[list[i]] - y[pair.Key]) < 1e-01;
            //        }
            //        if(is_regular_cell) { // if its regular we interpolate the value x[i]
            //            this.x[pair.Key] = this.Spline.Interpolate(y[pair.Key]);
            //        } else { //if not e set the value to something under minX
            //            this.x[pair.Key] = ExtrapolationFunc(y[pair.Key]); 
            //        }

            //        m_AllParams[pair.Key] = x[pair.Key];
            //    } else {
            //        this.x[pair.Key] = 3.0;
            //        m_AllParams[pair.Key] = x[pair.Key];
            //    }
            //}


            //for(int i = 0; i <  this.y.Length; i++) {

            //    this.x[i] = this.Spline.Interpolate(y[i]);
            //    m_AllParams[i] = x[i];

            //}
            //if(m_AllParams.Length > this.y.Length) {
            //    for(int i = 0; i < this.y.Length; i++) {
            //        double diff = Spline.Differentiate(y[i]);
            //        if(Math.Abs(diff) < 2) {
            //            m_AllParams[i + y.Length] = diff;
            //        } else {
            //            m_AllParams[i + y.Length] = Math.Sign(diff)*2;
            //        }

            //    }
            //}
            #endregion


            GetSpline();
        }

        public void GetSplineOverDeterminedOld(MultidimensionalArray points) {
            Dictionary<int, bool> PointsInCell = new Dictionary<int, bool>();
            Dictionary<int, List<int>> CellToPoints = new Dictionary<int, List<int>>();
            for(int j = 0; j < y.Length - 1; j++) {
                CellToPoints.Add(j, new List<int>());
            }
            var xPoints = points.ExtractSubArrayShallow(-1, 0).To1DArray();
            var yPoints = points.ExtractSubArrayShallow(-1, 1).To1DArray();
            for(int i = 0; i < yPoints.Length; i++) {
                for(int j = 0; j < y.Length - 1; j++) {
                    if(yPoints[i] > y[j] && yPoints[i] < y[j + 1]) {
                        CellToPoints.TryGetValue(j, out List<int> list);
                        list.Add(i);
                    }
                }
            }
            if(points.Lengths[1] == 2) {

                //find maximum points in Minimum Points
                double[] twoBiggestY;
                double[] twoSmallestY;
                if(yPoints[0] > yPoints[1]) {
                    twoBiggestY = new double[] { yPoints[0], yPoints[1] };
                    twoSmallestY = new double[] { yPoints[1], yPoints[0] };
                } else {
                    twoBiggestY = new double[] { yPoints[1], yPoints[0] };
                    twoSmallestY = new double[] { yPoints[0], yPoints[1] };
                }
                double yMax = twoBiggestY[0];
                double yMin = twoSmallestY[0];
                for(int i = 1; i < yPoints.Length; i++) {
                    yMax = twoBiggestY[0];
                    yMin = twoSmallestY[0];
                    if(yPoints[i] > yMax) {
                        twoBiggestY[1] = twoBiggestY[0];
                        twoBiggestY[0] = yPoints[i];
                    }
                    if(yMin > yPoints[i]) {
                        twoSmallestY[1] = twoSmallestY[0];
                        twoSmallestY[0] = yPoints[i];
                    }
                }
                yMax = twoBiggestY[0];
                yMin = twoSmallestY[0];

                for(int i = 0; i < y.Length; i++) {
                    if(y[i] < yMin) {
                        PointsInCell.Add(i, false);
                    } else if(y[i] > yMax) {
                        PointsInCell.Add(i, false);
                    } else {
                        PointsInCell.Add(i, true);
                    }
                }

                //List<double> interpolationY = new List<double>();
                //List<double> interpolationX = new List<double>();
                //for(int iY = 0; iY < y.Length-1; iY++) {
                //    interpolationY.Add(y[iY]);
                //    interpolationY.Add((y[iY]+ y[iY+1])/2);

                //    double closestPointiX = points[0, 0];
                //    double closestPointiX2 = points[0, 0];
                //    double deltaiY = Math.Abs(y[iY] - points[0, 1]);
                //    double deltaiY2 = Math.Abs((y[iY] + y[iY + 1])/ 2 - points[0, 1]);

                //    for(int i = 1; i < points.Lengths[0]; i++) {

                //        double deltaiY_test = Math.Abs(y[iY] - points[i, 1]);
                //        double deltaiY2_test = Math.Abs((y[iY] + y[iY + 1]) / 2 - points[i, 1]);
                //        if(deltaiY_test< deltaiY) {
                //            closestPointiX = points[i, 0];
                //        }
                //        if(deltaiY2_test < deltaiY2) {
                //            closestPointiX2 = points[i, 0];
                //        }

                //    }
                //    interpolationX.Add(closestPointiX);
                //    interpolationX.Add(closestPointiX2);
                //}
                //var xPoints = interpolationX.ToArray();
                //var yPoints = interpolationY.ToArray();
                #region determine cPoints (1st Derivative)
                var cPoints = new double[yPoints.Length];
                {
                    //forward FD on first point
                    var eps = (yPoints[1] - yPoints[0]);
                    var diff = xPoints[1] - xPoints[0];
                    cPoints[0] = diff / eps;
                    //backward FD on the rest
                    for(int i = 1; i < cPoints.Length; i++) {
                        eps = (yPoints[i] - yPoints[i - 1]);
                        diff = xPoints[i] - xPoints[i - 1];
                        cPoints[i] = diff / eps;
                    }
                }
                #endregion
                #region determine cPoints (2nd Derivative)
                //var cPoints = new double[yPoints.Length];
                //{
                //    //right/forward FD on first point
                //    var eps = ((yPoints[1] - yPoints[0]) + (yPoints[2] - yPoints[1])) / 2;
                //    var diff = xPoints[2] - 2 * xPoints[1] + xPoints[0];
                //    cPoints[0] = diff / (eps * eps);
                //    //don't include first and last point, use central FDs in the middle-points
                //    for(int i = 1; i < cPoints.Length - 1; i++) {
                //        eps = ((yPoints[i] - yPoints[i - 1]) + (yPoints[i + 1] - yPoints[i])) / 2;
                //        diff = xPoints[i + 1] - 2 * xPoints[i] + xPoints[i - 1];
                //        cPoints[i] = diff / (eps * eps);
                //    }
                //    //do backward FDS on the last one
                //    var iEnd = yPoints.Length - 1;
                //    eps = ((yPoints[iEnd - 1] - yPoints[iEnd - 2]) + (yPoints[iEnd] - yPoints[iEnd - 1])) / 2;
                //    diff = xPoints[iEnd] - 2 * xPoints[iEnd - 1] + xPoints[iEnd - 2];
                //    cPoints[iEnd] = diff / (eps * eps);
                //}
                #endregion
                GetSpline(xPoints, yPoints, cPoints);
            } else if(points.Lengths[1] == 3) {
                var cPoints = points.ExtractSubArrayShallow(-1, 2).To1DArray();
                GetSpline(xPoints, yPoints, cPoints);

            } else {
                throw new ArgumentException("array must have 2-3 columns");
            }
            //after the Spline is constructed we can get the values for x,y,c
            this.x = new double[y.Length];
            if(PointsInCell.Count == y.Length) {
                foreach(KeyValuePair<int, bool> pair in PointsInCell) {
                    if(pair.Value) {
                        CellToPoints.TryGetValue(pair.Key, out List<int> list);
                        bool is_regular_cell = false;
                        for(int i = 0; i < list.Count; i++) {
                            is_regular_cell = is_regular_cell || Math.Abs(yPoints[list[i]] - y[pair.Key]) < 1e-01;
                        }
                        if(is_regular_cell) {
                            this.x[pair.Key] = this.Spline.Interpolate(y[pair.Key]);
                        } else {
                            this.x[pair.Key] = 0.0;
                        }

                        m_AllParams[pair.Key] = x[pair.Key];
                    } else {
                        this.x[pair.Key] = 3.0;
                        m_AllParams[pair.Key] = x[pair.Key];
                    }
                }

            }
            //for(int i = 0; i <  this.y.Length; i++) {

            //    this.x[i] = this.Spline.Interpolate(y[i]);
            //    m_AllParams[i] = x[i];

            //}
            if(m_AllParams.Length > this.y.Length) {
                for(int i = 0; i < this.y.Length; i++) {
                    double diff = Spline.Differentiate(y[i]);
                    if(Math.Abs(diff) < 2) {
                        m_AllParams[i + y.Length] = diff;
                    } else {
                        m_AllParams[i + y.Length] = Math.Sign(diff) * 2;
                    }

                }
            }
            GetSpline();
        }
        public SplineOptiLevelSet CloneAs() {
            return (SplineOptiLevelSet)Clone();
        }
        public override object Clone() {
            var ret = new SplineOptiLevelSet(Basis, m_LevelSet, m_Tracker, m_LevelSetNumber);
            y.CopyTo(ret.y, 0);
            for(int i = 0; i < m_AllParams.Length; i++) {
                ret.m_AllParams[i] = m_AllParams[i];
            }
            ret.Interpolate();
            return ret;
        }

        public static void EmbeddInLevelSet(IInterpolation spline, SinglePhaseField levelSet, CellMask region = null) {
            if(region == null) {
                region = CellMask.GetFullMask(levelSet.GridDat);
            }
            levelSet.Clear(region);
            levelSet.ProjectField(
                1.0,
                LevelSet(spline),
                new BoSSS.Foundation.Quadrature.CellQuadratureScheme(true, region));
        }

        static ScalarFunction LevelSet(IInterpolation spline) {
            return (nodes, results) => {
                for(int i = 0; i < nodes.Lengths[0]; ++i) {
                    //phi(x,y) = position(y) - x
                    results[i] = spline.Interpolate(nodes[i, 1]);
                    results[i] -= nodes[i, 0];
                    results[i] *= -1; // negative x is <0
                };
            };
        }
        #region IOptiLevelSet members
        public IGridData GetGrid() {
            return Basis.GridDat;
        }

        public int GetLength() {
            return m_AllParams.Lengths[0];
        }

        public virtual void SetParam(int index, double val) {
            m_AllParams[index] = val;
        }
        /// <summary>
        /// gets the param
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public double GetParam(int index) {
            return m_AllParams[index];
        }

        /// <summary>
        /// need still to think about it
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public bool IsInNearBand(int index) {
            throw new NotImplementedException();
        }

        public double[] GetParamsAsArray() {
            return m_AllParams.To1DArray();
        }

        public void CopyParamsFrom(IOptiLevelSet source) {
            for(int i = 0; i < m_AllParams.Length; i++) {
                m_AllParams[i] = source.GetParam(i);
            };
        }

        public void AccToParam(int index, double acc) {
            m_AllParams[index] += acc;
        }
        /// <summary>
        /// returns name of param
        /// </summary>
        /// <param name="index">index of param</param>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public string GetParamName(int index) {
            if(index < numberOfNodes) {
                return "y[" + index + "]" + m_AllParams[index];
            } else {
                return "c2[" + index + "]" + m_AllParams[index];
            }

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="degree"></param>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public SinglePhaseField ToSinglePhaseField(int degree) {
            var ret_LS = new SinglePhaseField(new Basis(GridDat, degree), "splineLevelSet");

            ret_LS.AccLaidBack(1.0, this);
            return ret_LS;
        }
        /// <summary>
        /// no Trans Mat needed
        /// </summary>
        /// <param name="targetLS"></param>
        public void AssembleTransMat(LevelSet targetLS) {
            ;
        }

        public void ProjectOntoLevelSet(LevelSet targetLS) {
            var p = 4;
            //m_AllParams.To1DArray().SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\mAP{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"bImAP{p}*").Length}.txt");
            //this.CoordinateVector.SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\tOLS{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"bItOLS{p}*").Length}.txt");
            //targetLS.Clear();
            //targetLS.ProjectFromForeignGrid(1.0, this);
            Interpolate();
            //m_AllParams.To1DArray().SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\mAP{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"mAP{p}*").Length}.txt");
            //this.CoordinateVector.SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\tOLS{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"tOLS{p}*").Length}.txt");

            //targetLS.CoordinateVector.SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\tLS{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"tLS{p}*").Length}.txt");


            if(Basis == targetLS.Basis) {
                targetLS.CoordinateVector.CopyFrom(CoordinateVector, 0);
            } else {
                targetLS.Clear();
                targetLS.ProjectFromForeignGrid(1.0, this);
            }
            //targetLS.CoordinateVector.SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\tLS{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"tLS{p}*").Length}.txt");
            
            //Interpolate();
            //targetLS.Clear();
            //targetLS.ProjectFromForeignGrid(1.0, this);
            //targetLS.CoordinateVector.SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\tLS{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"tLS{p}*").Length}.txt");
            //var clone = this.Clone();
            //targetLS.Clear();
            //Interpolate();
            //targetLS.ProjectFromForeignGrid(1.0, this);
            //targetLS.CoordinateVector.SaveToTextFile($@"C:\Users\sebastian\Documents\Forschung\tLS{p}_{Directory.GetFiles(@"C:\Users\sebastian\Documents\Forschung", $"tLS{p}*").Length}.txt");
            
        }

    public void Print() {

            for(int i = 0; i < m_AllParams.Length; i++) {
                Console.Write(GetParamName(i) + ", ");
            }

        }
        /// <summary>
        /// here will need to do some root finding stuff
        /// </summary>
        /// <param name="sourceLS"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void ProjectFromForeignLevelSet(SinglePhaseField sourceLS) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// here will need to do some root finding stuff
        /// </summary>
        /// <param name="sourceLS"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void ProjectFromLevelSet(ConventionalDGField sourceLS) {
            ProjectFromForeignLevelSet((SinglePhaseField)sourceLS);
        }

        public bool TestOrthonormality() {
            return true;
        }
        /// <summary>
        /// 
        /// if there is a isolated cut-cell pair (that is if 
        /// 
        ///         x_{i_i},x_{i+1} \notin [x_{min},x_{max}] 
        ///         but  x_{i} \in [x_{min},x_{max}]
        ///         
        /// we then want to get rid of it by setting x_{i} \notin [x_{min},x_{max}] and the derivatives to zero
        /// 
        /// </summary>
        public void Reinitialize(double L, double kappa_S) {
            //List<int> CPs = new List<int>();
            //{
            //    for(int iY = 1; iY < y.Length - 2; iY++) {
            //        if ((m_AllParams[iY] - m_AllParams[iY-1])* (m_AllParams[iY] - m_AllParams[iY +1]) > 0) {
            //            CPs.Add(iY);
            //        }
            //    }
            //}
            //List<int> Ps = new List<int>();
            //{
            //    for(int iY = 1; iY < CPs.Count - 1; iY++) {
            //        if(Math.Abs((m_AllParams[iY] - m_AllParams[iY - 1])) > L*kappa_S && Math.Abs((m_AllParams[iY] - m_AllParams[iY +1])) > L * kappa_S) {
            //            Ps.Add(iY);
            //        }

            //    }
            //}

            //Interpolate();

            var metrics = this.m_Tracker.GetXDGSpaceMetrics(this.m_Tracker.SpeciesIdS, this.GetDegree());
            var CutCellVolumes = metrics.CutCellMetrics.CutCellVolumes;
            var cutCells = m_Tracker.Regions.GetCutCellMask().ItemEnum;

            var agglomerator = this.m_Tracker.GetAgglomerator(this.m_Tracker.SpeciesIdS.ToArray(), 3, 0.1); // maybe parameterize 
            var CutCellVolumesAgg = agglomerator.CutCellVolumes;

            var grid = this.m_Tracker.GridDat;
            double maxCellVol = grid.Cells.GetCellVolume(0); //in a cartesian mesh all cells have the same volume
            var eps = Math.Abs(y[0] - y[1]) / 10;

            for(int iY = 0; iY < y.Length; iY++) { // loop over all y-Cells
                //first we check whether the cell y_i,y_i+1 is isolated,
                //that means that the x_i and x_i+1 are not inside the domain 
                //Then we check whether the respective neighbors are isolated to
                // a special treaty for the first and the last cell has to be done
                // if a neighbor is isolated we can "freely" manipulate its parameters
                bool hasRightIsolatedNeighbor = true;
                bool hasLeftIsolatedNeighbor = true;
                bool isIsolatedCell = false;
                if(iY == 0) {
                    hasRightIsolatedNeighbor = (x[iY + 2] < xMin || x[iY + 2] > xMax) && (x[iY + 1] < xMin || x[iY + 1] > xMax);
                    isIsolatedCell = (x[iY] < xMin || x[iY] > xMax) && (x[iY + 1] < xMin || x[iY + 1] > xMax);
                }else if(iY==y.Length-1){
                    hasLeftIsolatedNeighbor = (x[iY - 2] < xMin || x[iY - 2] > xMax) && (x[iY - 1] < xMin || x[iY - 1] > xMax);
                    isIsolatedCell = (x[iY] < xMin || x[iY] > xMax) && (x[iY - 1] < xMin || x[iY - 1] > xMax);
                } else {
                    hasLeftIsolatedNeighbor = (x[iY - 2] < xMin || x[iY - 2] > xMax) && (x[iY - 1] < xMin || x[iY - 1] > xMax);
                    isIsolatedCell = (x[iY] < xMin || x[iY] > xMax) && (x[iY - 1] < xMin || x[iY - 1] > xMax);
                    hasRightIsolatedNeighbor = (x[iY + 2] < xMin || x[iY + 2] > xMax) && (x[iY + 1] < xMin || x[iY + 1] > xMax);
                }
                if(isIsolatedCell) {
                    grid.LocatePoint(new double[] { xMin + eps, y[iY] + eps }, out long GlobIdRight, out long GlobIndexRight, out bool isInsideRight, out bool onProcRight);
                    grid.LocatePoint(new double[] { xMax - eps, y[iY] + eps }, out long GlobIdLeft, out long GlobIndexLeft, out bool isInsideLeft, out bool onProcLeft);

                    bool isLeftBoundaryCutCell = cutCells.ToList().Contains((int)GlobIndexLeft);
                    bool isRightBoundaryCutCell = cutCells.ToList().Contains((int)GlobIndexRight);

                    double a = (x[iY + 1] - x[iY]) / (y[iY + 1] - y[iY]);
                    if(isLeftBoundaryCutCell || isRightBoundaryCutCell) { //only do something if its a cut cell
                        if(!hasLeftIsolatedNeighbor && !hasRightIsolatedNeighbor) {
                        //here it is not clear what to do, it means that an isolated cell has two non-isolated neighbors
                    } else if(hasLeftIsolatedNeighbor && hasRightIsolatedNeighbor) {
                        
                            m_AllParams[iY + y.Length] = a;
                            m_AllParams[iY + 1 + y.Length] = a;
                        
                    } else {
                        if(hasRightIsolatedNeighbor) {
                            if(isLeftBoundaryCutCell || isRightBoundaryCutCell) {
                                
                            }
                        } else if(hasLeftIsolatedNeighbor) {

                        } else {
                            //should not be reached
                        }
                    }
                }
            }
            }
            Interpolate();
        }

        public void ProjectFromFunction(Func<double[], double> initialShockPostion) {
            for(int i = 0; i < y.Length; ++i) {
                m_AllParams[i] = initialShockPostion(new double[] { 0, y[i] });
            }
            var eps = 1e-08;
            //This will be the first derivative for Quadratic and Cubic Splines
            for(int i = y.Length; i < m_AllParams.Length; ++i) {
                m_AllParams[i] = (initialShockPostion(new double[] { 0, y[i- y.Length] +eps }) - initialShockPostion(new double[] { 0, y[i - y.Length] }) )/ eps;
                //m_AllParams[i] = 0;
            }
            Interpolate();
        }

        public int GetDegree() {
            return Basis.Degree > 3 ? Basis.Degree : 3;
        }

        public LevelSet ToLevelSet(int v) {
            var ret_LS = new LevelSet(new Basis(GridDat, v), "splineLevelSet");
            return this;
        }

        public void AssembleTracker() {
            var dummy_LS = new LevelSet(Basis, "dummy_LS");
            m_thisTracker = new LevelSetTracker((GridData)Basis.GridDat, m_Tracker.CutCellQuadratureType, GetDegree(), m_Tracker.GetSpeciesSeparatedByLevSet(m_LevelSetNumber).ToArray(), dummy_LS);
        }
        public double Norm(double[] levelSetStepCoordinates) {
            //Old version of calculation of norm (works better for ECCOMAS Test Case Scalar Advection)
            double norm = 0;
            for(int i = 0; i < this.GetLength(); i++) {
                norm += levelSetStepCoordinates[i]* levelSetStepCoordinates[i];
            }
            return norm.Sqrt();
        }
        /// <summary>
        /// Experimental Regularization 
        /// </summary>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public MsrMatrix GetRegMatrix() {
            var length_l = this.GetLength();
            var reg = new MsrMatrix(length_l, length_l, 1, 1);
            //identity
            IMutuableMatrixEx_Extensions.AccEyeSp(reg, 1.0);
            return reg;

            //something else
            //for(int iY = 0; iY < y.Length; iY++) {
            //    reg[iY, iY] = 1; // this is inversely proportional to the size of the minimal cut Cell
            //    reg[iY + y.Length, iY + y.Length] = 2;  // this is a 1st derivative so lets square it.
            //}
            //return reg;

            //Experimental
            //return ExperimentalReg();
        }

        public MsrMatrix ExperimentalReg() {
            var length_l = this.GetLength();
            var reg = new MsrMatrix(length_l, length_l, 1, 1);
            //identity
            //IMutuableMatrixEx_Extensions.AccEyeSp(reg, 1.0);

            //Experimental Reg
            var metrics = this.m_Tracker.GetXDGSpaceMetrics(this.m_Tracker.SpeciesIdS, this.GetDegree());
            var CutCellVolumes = metrics.CutCellMetrics.CutCellVolumes;

            //Function that determines the minimal CutCell Volume for 2 cells across all species
            double MinCutCellVolume(long cell1, long cell2) {
                double ret = double.MaxValue;
                foreach(SpeciesId spc in this.m_Tracker.SpeciesIdS) {
                    var volumesForSpc1 = CutCellVolumes[spc][(int)cell1];
                    var volumesForSpc2 = CutCellVolumes[spc][(int)cell2];
                    if(volumesForSpc1 != 0 && volumesForSpc2 != 0) {
                        ret = volumesForSpc1 < ret ? volumesForSpc1 : ret;
                        ret = volumesForSpc2 < ret ? volumesForSpc2 : ret;
                    }
                }
                if(ret == double.MaxValue) {
                    throw new Exception("something went wrong");
                }
                return ret;
            }

            var agglomerator = this.m_Tracker.GetAgglomerator(this.m_Tracker.SpeciesIdS.ToArray(), 3, 0.1); // maybe parameterize 
            var CutCellVolumesAgg = agglomerator.CutCellVolumes;

            var grid = this.m_Tracker.GridDat;
            double maxCellVol = grid.Cells.GetCellVolume(0); //in a cartesian mesh all cells have the same volume

            var eps = Math.Abs(y[0] - y[1]) / 10;
            var diag = new double[y.Length];
            for(int iY = 0; iY < y.Length; iY++) {

                //this should give indices for 2 neighboring cells
                grid.LocatePoint(new double[] { this.GetParam(iY), y[iY] + eps }, out long GlobIdUp, out long GlobIndexUp, out bool isInsideUp, out bool onProcUp);
                grid.LocatePoint(new double[] { this.GetParam(iY), y[iY] - eps }, out long GlobIdDown, out long GlobIndexDown, out bool isInsideDown, out bool onProcDown);
                double minCCV;
                if(isInsideUp) {
                    if(isInsideDown) {
                        minCCV = MinCutCellVolume(GlobIndexUp, GlobIndexDown);
                    } else {
                        minCCV = MinCutCellVolume(GlobIndexUp, GlobIndexUp); //lower boundary of grid
                    }
                } else {
                    if(isInsideDown) {
                        minCCV = MinCutCellVolume(GlobIndexDown, GlobIndexDown); //upper boundary of grid
                    } else {
                        minCCV = maxCellVol; //point outside grid <----- we want a zero here ? or better a one? or should the result be infinity?
                    }
                }
                diag[iY] = Math.Sqrt(maxCellVol / minCCV);
                reg[iY, iY] = Math.Sqrt(maxCellVol / minCCV); // this is inversely proportional to the size of the minimal cut Cell
                reg[iY + y.Length, iY + y.Length] = maxCellVol / minCCV;  // this is a 1st derivative so lets square it.
            }
            return reg;
        }

        #endregion
    }
    public enum Mode {
        Deg1,
        Deg2,
        Deg3

    }
}
