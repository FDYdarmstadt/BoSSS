/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;

namespace BoSSS.Solution.LevelSetTools.TestCases {
    public class ZalesaksDisk {

        double Radius;
        double YCutout;
        double[] XCutout;

        double yEdgePointFanLeft ;
        double yEdgePointFanRight;

        Func<double[], double> ExpansionFanLeft ;
        Func<double[], double> ExpansionFanRight ;

        Func<double[], double> PointLeft ;
        Func<double[], double> PointRight ;

        Func<double[], double> LineLeft ;
        Func<double[], double> LineRight ;
        Func<double[], double> LineBottom ;
        Func<double[], double> Circle ;

        public ZalesaksDisk(double[] XCutout, double YCutout, double Radius) {
            this.XCutout = XCutout;
            this.YCutout = YCutout;
            this.Radius = Radius;

            yEdgePointFanLeft = Math.Sqrt(Radius.Pow2() - XCutout[0].Pow2());
            yEdgePointFanRight = Math.Sqrt(Radius.Pow2() - XCutout[1].Pow2());

            ExpansionFanLeft = (x => Math.Sqrt((x[0] - XCutout[0]).Pow2() + (x[1] - yEdgePointFanLeft).Pow2()));
            ExpansionFanRight = (x => Math.Sqrt((x[0] - XCutout[1]).Pow2() + (x[1] - yEdgePointFanRight).Pow2()));

            PointLeft = (x => -Math.Sqrt((x[0] - XCutout[0]).Pow2() + (x[1] - YCutout).Pow2()));
            PointRight = (x => -Math.Sqrt((x[0] - XCutout[1]).Pow2() + (x[1] - YCutout).Pow2()));

            LineLeft = x => x[0] - XCutout[0];
            LineRight = x => -(x[0] - XCutout[1]);
            LineBottom = x => (x[1] - YCutout);
            Circle = x => (Math.Sqrt(x[0].Pow2() + x[1].Pow2()) - Radius);
        }


        public double NonSignedDistance(double[] X) {
            //NonSigned-Distance, No Expansion Fans:
            double CircleDistance = (Math.Sqrt(X[0].Pow2() + X[1].Pow2()) - Radius);
            double CutoutDistance = Math.Min(Math.Min((X[0] - XCutout[0]), (-X[0] + XCutout[1])), X[1] - YCutout);
            return Math.Max(CircleDistance, CutoutDistance);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="X">Spatial Coordinate</param>
        /// <param name="initialFunction">Function to extend</param>
        /// <returns>
        /// index[0]: Signed-distance level-set value
        /// index[1]: Extension, based on signed distance field
        /// </returns>
        public double[] SignedDistance(double[] X, Func<double[], double> initialFunction) {

            double levelSet = double.NaN;
            double AnalyticalSolution = double.NaN;
          
            Func<double[], double> CircleSolution = x => initialFunction(new double[] { Radius * Math.Cos(Math.Atan2(X[1], X[0])), Radius * Math.Sin(Math.Atan2(X[1], X[0])) });
            Func<double, double, double> soln = (x, y) => initialFunction(new double[] { x, y });

            //Full-on:
            //Outside Circle
            if (Circle(X) >= 0 || X[1] >= yEdgePointFanLeft || X[1] >= yEdgePointFanRight) {
                // Expansion Fans Top Side
                if (X[1] > YCutout && X[0] > XCutout[0] && X[0] < XCutout[1]) {
                    if (ExpansionFanLeft(X) < ExpansionFanRight(X)) {
                        levelSet = ExpansionFanLeft(X);
                        AnalyticalSolution = soln(XCutout[0], yEdgePointFanLeft);
                    }
                    else {
                        levelSet = ExpansionFanRight(X);
                        AnalyticalSolution = soln(XCutout[1], yEdgePointFanRight);
                    }
                }
                // 
                else {
                    levelSet = Circle(X);
                    AnalyticalSolution = CircleSolution(X);
                }
            }
            //InsideCircle, Positive side of Cutout
            else if (X[0] > XCutout[0] && X[0] < XCutout[1] && X[1] > YCutout) {
                double centerline = (XCutout[0] + XCutout[1]) / 2;
                // levelSet = Math.Min(LineBottom(X), Math.Min(LineRight(X), LineLeft(X)));
                if (X[0] < centerline) {
                    if (LineBottom(X) < LineLeft(X)) {
                        levelSet = LineBottom(X);
                        AnalyticalSolution = soln(X[0], YCutout);
                    }
                    else {
                        levelSet = LineLeft(X);
                        AnalyticalSolution = soln(XCutout[0], X[1]);
                    }
                }
                else {
                    if (LineBottom(X) < LineRight(X)) {
                        levelSet = LineBottom(X);
                        AnalyticalSolution = soln(X[0], YCutout);
                    }
                    else {
                        levelSet = LineRight(X);
                        AnalyticalSolution = soln(XCutout[1], X[1]);
                    }
                }


            }
            //Inside
            // Bottom Left
            else if (X[0] <= XCutout[0] && X[1] <= YCutout) {
                if (PointLeft(X) > Circle(X)) {
                    levelSet = PointLeft(X);
                    AnalyticalSolution = soln(XCutout[0], YCutout);
                }
                else {
                    levelSet = Circle(X);
                    AnalyticalSolution = CircleSolution(X);
                }
            }
            // Bottom Right
            else if (X[0] >= XCutout[1] && X[1] <= YCutout) {
                if (PointRight(X) > Circle(X)) {
                    levelSet = PointRight(X);
                    AnalyticalSolution = soln(XCutout[1], YCutout);
                }
                else {
                    levelSet = Circle(X);
                    AnalyticalSolution = CircleSolution(X);
                }
            }
            //Bottom Center
            else if (X[1] <= YCutout) {
                if (LineBottom(X) > Circle(X)) {
                    levelSet = LineBottom(X);
                    AnalyticalSolution = soln(X[0], YCutout);
                }
                else {
                    levelSet = Circle(X);
                    AnalyticalSolution = CircleSolution(X);
                }
            }
            //Top Left
            else if (X[0] < XCutout[0]) {
                if (LineLeft(X) > Circle(X)) {
                    levelSet = LineLeft(X);
                    AnalyticalSolution = soln(XCutout[0], X[1]);
                }
                else {
                    levelSet = Circle(X);
                    AnalyticalSolution = CircleSolution(X);
                }
            }
            //Top Right
            else {
                if (LineRight(X) > Circle(X)) {
                    levelSet = LineRight(X);
                    AnalyticalSolution = soln(XCutout[1], X[1]);
                }
                else {
                    levelSet = Circle(X);
                    AnalyticalSolution = CircleSolution(X);
                }
            }
            return new double[] { levelSet, AnalyticalSolution };

        }


        public double SignedDistanceLevelSet(double[] X) {
            return SignedDistance(X, new Func<double[],double>(x=>0))[0];
        }
        public double SignedDistanceExtension(double[] X, Func<double[], double> initialFunction) {
            return SignedDistance(X, initialFunction)[1];
        }


        public double GetArea() {
            double Circle = Math.PI * Radius * Radius;
            double TriangleLeft = -yEdgePointFanLeft * XCutout[0] * 0.5;
            double TriangleRight = yEdgePointFanRight * XCutout[1] * 0.5;
            double LowerRectangle = YCutout * (XCutout[1] - XCutout[0]);
            double arc = (Math.Atan2(yEdgePointFanRight, XCutout[1]) - Math.Atan2(yEdgePointFanLeft, XCutout[0])) * Radius*Radius*0.5;

            return Circle - TriangleLeft - TriangleRight - arc + LowerRectangle;
        }
    }
}
