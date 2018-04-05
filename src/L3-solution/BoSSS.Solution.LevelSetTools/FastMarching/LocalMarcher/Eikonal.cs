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
using static System.Math;

namespace BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher {
    static class Eikonal {

        public static double approximate(double u_l, double u_r, double u_b, double u_t, double h_l, double h_r, double h_b, double h_t, double R) {

            //Square R
            if (R != 1) {
                R = Pow(R, 2);
            }

            //create List of points of interest
            point U_l = new point( u_l, h_l, "x");
            point U_r = new point(u_r, h_r, "x");
            point U_t = new point(u_t, h_t, "y");
            point U_b = new point(u_b, h_b, "y");
            point U_lr = U_l.solveIntersectionOfRespLine( U_r);
            point U_tb = U_t.solveIntersectionOfRespLine( U_b);

            point[] points = new point[] { U_l, U_r, U_t, U_b, U_lr, U_tb };
            Array.Sort(points);

#if DEBUG
            //Throw exception when all Points are undefined; 
            if (points[0].pos == double.MaxValue) {
                throw new ArgumentException("At least one u_i must be defined");
            }
            if (h_l <= 0 || h_r <= 0 || h_b <= 0 || h_t <= 0) {
                throw new ArgumentException("H must be > 0");
            }
#endif

            //Create function that selects the active conditions in each max() 
            function Eikonal = new function(R);
      
            //Check points
            for(int i = 0; i < points.Length-1; ++i) {
                Eikonal.update(points[i]);
                if (Eikonal.value(points[i + 1]) > R) {
                    return Eikonal.solve();
                }
            }
            Eikonal.update(points[points.Length-1]);
            return Eikonal.solve();
        }

        class point : IComparable<point> {
            public double pos; 
            public double u;
            public double h;
            public string dimension;
            public Boolean intersection;
            
            public point( double u_i, double h_i, string Equation) {
                pos = u_i;
                u = u_i;
                h = h_i;
                dimension = Equation;
                intersection = false;
            }

            public point(double Pos, double u_i, double h_i, string Equation) {
                pos = Pos;
                u = u_i;
                h = h_i;
                dimension = Equation;
                intersection = true;
            }

            public int CompareTo(point Point) {
                if (Point.pos == this.pos) {
                    return 0;
                }
                if (Point.pos < this.pos) {
                    return 1;
                }
                return -1;
            }

            public point solveIntersectionOfRespLine(point V) {
                if (dimension != V.dimension) {
                    throw new ArgumentException("Points belong to different Dimensions");
                }
                //If one line is infinity, set pos to infinity
                if (V.u == double.MaxValue || u == double.MaxValue) {
                    return new point(double.MaxValue, Double.NaN, Double.NaN, dimension);
                }

                //If the lines are parallel, set pos to infinity
                if (Abs(h - V.h) < 1e-13) {
                    return new point(double.MaxValue, Double.NaN, Double.NaN, dimension);
                }
                double u_intersect = (V.h * u - h * V.u) / (V.h - h);

                //If the intersection is below 0 in R direction, set pos to infinity
                if ((u_intersect - V.u) / V.h < 0 || double.IsInfinity(u_intersect)) {
                    return new point(double.MaxValue, Double.NaN, Double.NaN, dimension);
                }

                //Attach the data of the active function for u > u_intersect
                double respU;
                double respH;
                if (h < V.h) {
                    respU = u;
                    respH = h;
                } else {
                    respU = V.u;
                    respH = V.h;
                }
                point intersection = new point(u_intersect, respU, respH, dimension);

                return intersection;
            }

        }

        class function {
            double R;
            point X;
            point Y; 

            public function(double r) {
                R = r;                     
            }

            public void update(point u) {
                if (u.pos != double.MaxValue) {
                    if (u.dimension == "x" && X == null) {
                        X = u;
                    }
                    if (u.dimension == "y" && Y == null) {
                        Y = u;
                    }
                    if (u.intersection == true) {
                        if (u.dimension == "x") {
                            X = u;
                        }
                        if (u.dimension == "y") {
                            Y = u;
                        }
                    }
                }
            }

            public double value(point U) {
                if (U.pos == double.MaxValue) {
                    return R  + 1;
                }
                if (X == null && Y == null) {
                    throw new ArgumentException("Insufficient data");
                }
                if (X == null) {
                    return valueSingle(U.pos, Y);
                }
                if (Y == null) {
                    return valueSingle(U.pos, X);
                }
                return valueDouble(U.pos, X, Y);
            }

            double valueSingle(double u, point U) {
                return Pow(((u - U.u) / U.h), 2);
            }

            double valueDouble(double u, point X, point Y) {
                return (Pow((u - X.u) / X.h, 2) + Pow( (u - Y.u) / Y.h, 2)); 
            }

            public double solve() {
                if (X == null && Y == null) {
                    throw new ArgumentException("Insufficient data");
                }
                if (X == null) {
                    return solveSingle(Y); 
                }
                if (Y == null) {
                    return solveSingle(X);
                }
                return solveDouble(X,Y);
            }

            double solveSingle(point U) {
                double u;
                if (R == 1) {
                    u = U.h + U.u;
                } else {
                    u = U.u + Sqrt(R * Pow(U.h, 2));
                }
                return u;
            }

            double solveDouble(point U, point V) {
                double a = (U.u / Pow(U.h, 2) + V.u / Pow(V.h, 2));
                double b = (1 / Pow(U.h, 2) + 1 / Pow(V.h, 2));
                double u = a / b + Sqrt(Pow(a / b, 2) - (Pow(U.u / U.h, 2) + Pow(V.u / V.h, 2) - R) / b);
                return u;
            }

        }

    }
}
