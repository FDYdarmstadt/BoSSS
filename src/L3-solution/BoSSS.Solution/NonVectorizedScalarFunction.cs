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
using BoSSS.Platform;
using BoSSS.Foundation;
using ilPSP;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// 1D function f(x), non-vectorized
    /// </summary>
    public delegate double _1D(double x0);

    /// <summary>
    /// 2D function f(x,y) or (1+1)D function f(t,x), non-vectorized
    /// </summary>
    public delegate double _2D(double x0, double x1);

    /// <summary>
    /// 2D function f(x,y,z) or (1+2)D function f(t,x,y), non-vectorized
    /// </summary>
    public delegate double _3D(double x0, double x1, double x3);

    /// <summary>
    /// (1+3)D - function f(t,x,y,z), non-vectorized
    /// </summary>
    public delegate double _4D(double t, double x1, double x2, double x3);

    /// <summary>
    /// extension functions (vectorization) of <see cref="_1D"/>, <see cref="_2D"/>, <see cref="_3D"/> and <see cref="_4D"/>;
    /// </summary>
    public static class NonVectorizedScalarFunction {

        /// <summary>
        /// Projection of a 1D-function <paramref name="f"/> onto a DG-Field
        /// </summary>
        public static void ProjectField(this DGField u, _1D f) {
            if (u.Basis.GridDat.SpatialDimension != 1)
                throw new ArgumentException("mismatch in spatial dimension");
            u.ProjectField(f.Vectorize());
        }

        /// <summary>
        /// Projection of a 2D-function <paramref name="f"/> onto a DG-Field
        /// </summary>
        public static void ProjectField(this DGField u, _2D f) {
            if (u.Basis.GridDat.SpatialDimension != 2)
                throw new ArgumentException("mismatch in spatial dimension");
            u.ProjectField(f.Vectorize());
        }

        /// <summary>
        /// Projection of a 3D-function <paramref name="f"/> onto a DG-Field
        /// </summary>
        public static void ProjectField(this DGField u, _3D f) {
            if (u.Basis.GridDat.SpatialDimension != 3)
                throw new ArgumentException("mismatch in spatial dimension");
            u.ProjectField(f.Vectorize());
        }

        /// <summary>
        /// Initializes a DG field with random values.
        /// </summary>
        public static void InitRandom(this DGField u, Random r = null) {
            if (r == null)
                r = new Random();

            int J = u.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < J; j++) {
                int N = u.Basis.GetLength(j);
                for (int n = 0; n < N; n++) {
                    u.Coordinates[j, n] = r.NextDouble();
                }
            }
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double> Convert_Xt2X(this Func<double[],double,double> f, double time) {
            return (X => f(X, time));
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double,double> Convert_X2Xt(this Func<double[],double> f) {
            return (X,t) => f(X);
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double> Convert_x2X(this _1D f) {
            return (double[] X) => f(X[0]);
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[], double> Convert_x2X(this Func<double, double> f) {
            return (double[] X) => f(X[0]);
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double> Convert_xy2X(this _2D f) {
            return (double[] X) => f(X[0], X[1]);
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double> Convert_xyz2X(this _3D f) {
            return (double[] X) => f(X[0], X[1], X[2]);
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double,double> Convert_tx2Xt(this _2D f) {
            return (double[] X, double t) => f(t, X[0]);
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double,double> Convert_txy2Xt(this _3D f) {
            return (double[] X, double t) => f(t, X[0], X[1]);
        }

        /// <summary>
        /// Scalar function conversion.
        /// </summary>
        public static Func<double[],double,double> Convert_txyz2Xt(this _4D f) {
            return (double[] X, double t) => f(t, X[0], X[1], X[2]);
        }


        /// <summary>
        /// Projection of a function <paramref name="f"/> onto a DG-Field
        /// </summary>
        public static void ProjectField(this DGField u, Func<double[], double> f) {
            u.ProjectField(f.Vectorize());
        }

        /// <summary>
        /// Projection of a function <paramref name="f"/> onto a DG-Field
        /// </summary>
        public static void ProjectField(this DGField u, double alpha, Func<double[], double> f, CellQuadratureScheme scheme = null) {
            u.ProjectField(alpha, f.Vectorize(), scheme);
        }

        /// <summary>
        /// L2Error w.r.t. a 1D-function <paramref name="f"/> of a DG-Field
        /// </summary>
        public static double L2Error(this DGField u, _1D f) {
            if (u.Basis.GridDat.SpatialDimension != 1)
                throw new ArgumentException("mismatch in spatial dimension");
            return u.L2Error(f.Vectorize());
        }

        /// <summary>
        /// L2Error w.r.t. a 2D-function <paramref name="f"/> of a DG-Field
        /// </summary>
        public static double L2Error(this DGField u, _2D f) {
            if (u.Basis.GridDat.SpatialDimension != 2)
                throw new ArgumentException("mismatch in spatial dimension");
            return u.L2Error(f.Vectorize());
        }

        /// <summary>
        /// L2Error w.r.t. a 3D-function <paramref name="f"/> of a DG-Field
        /// </summary>
        public static double L2Error(this DGField u, _3D f) {
            if (u.Basis.GridDat.SpatialDimension != 3)
                throw new ArgumentException("mismatch in spatial dimension");
            return u.L2Error(f.Vectorize());
        }


        /// <summary>
        /// L2Error w.r.t. a function <paramref name="f"/> of a DG-Field
        /// </summary>
        public static double L2Error(this DGField u, Func<double[], double> f) {
            return u.L2Error(f.Vectorize());
        }


        /// <summary>
        /// Vectorized 1D function (<see cref="ScalarFunction"/>) from a scalar implementation
        /// </summary>
        /// <param name="f">calling sequence: f(x)</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this _1D f) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                if (D != 1)
                    throw new ArgumentException("wrong spatial dimension.");

                for (int i = 0; i < inp.GetLength(0); i++) {
                    double x = inp[i, 0];

                    res[i] = f(x);
                }
            });
        }

        /// <summary>
        /// Vectorized 1D function (<see cref="ScalarFunction"/>) from a scalar implementation, with fixed time
        /// </summary>
        /// <param name="f">calling sequence: f(<paramref name="time"/>,x)</param>
        /// <param name="time">fixed time</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this _2D f, double time) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                if (D != 1)
                    throw new ArgumentException("wrong spatial dimension.");

                for (int i = 0; i < inp.GetLength(0); i++) {
                    double x = inp[i, 0];

                    res[i] = f(time, x);
                }
            });
        }

        /// <summary>
        /// Vectorized 2D function (<see cref="ScalarFunction"/>) from a scalar implementation
        /// </summary>
        /// <param name="f">calling sequence: f(x,y)</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this _2D f) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                if (D != 2)
                    throw new ArgumentException("wrong spatial dimension.");

                for (int i = 0; i < inp.GetLength(0); i++) {
                    double x = inp[i, 0];
                    double y = inp[i, 1];

                    res[i] = f(x, y);
                }
            });
        }

        /// <summary>
        /// Vectorized 2D function (<see cref="ScalarFunction"/>) from a scalar implementation
        /// </summary>
        /// <param name="f">calling sequence: f(x,y)</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this Func<double, double, double> f) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                if (D != 2)
                    throw new ArgumentException("wrong spatial dimension.");

                for (int i = 0; i < inp.GetLength(0); i++) {
                    double x = inp[i, 0];
                    double y = inp[i, 1];

                    res[i] = f(x, y);
                }
            });
        }

        /// <summary>
        /// Vectorized 2D function (<see cref="ScalarFunction"/>) from a scalar implementation, with fixed time
        /// </summary>
        /// <param name="f">calling sequence: f(<paramref name="time"/>,x,y)</param>
        /// <param name="time">fixed time</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this _3D f, double time) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                if (D != 2)
                    throw new ArgumentException("wrong spatial dimension.");

                for (int i = 0; i < inp.GetLength(0); i++) {
                    double x = inp[i, 0];
                    double y = inp[i, 1];

                    res[i] = f(time, x, y);
                }
            });
        }

        /// <summary>
        /// Vectorized 3D function (<see cref="ScalarFunction"/>) from a scalar implementation
        /// </summary>
        /// <param name="f">calling sequence: f(x,y,z)</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this _3D f) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                if (D != 3)
                    throw new ArgumentException("wrong spatial dimension.");

                for (int i = 0; i < inp.GetLength(0); i++) {
                    double x = inp[i, 0];
                    double y = inp[i, 1];
                    double z = inp[i, 2];

                    res[i] = f(x, y, z);
                }
            });
        }

        /// <summary>
        /// Vectorized function (<see cref="ScalarFunction"/>) from a scalar implementation
        /// </summary>
        /// <param name="f">calling sequence: f(x,y,z)</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this Func<double[], double> f) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                double[] X = new double[D];

                for (int i = 0; i < inp.GetLength(0); i++) {
                    for (int d = 0; d < D; d++)
                        X[d] = inp[i, d];

                    res[i] = f(X);
#if DEBUG
                    if (res[i].IsNaN())
                        throw new ArithmeticException("Vectorizing returns invalid values");
#endif
                }
            });
        }

        /// <summary>
        /// Vectorized function (<see cref="ScalarFunction"/>) from a scala
        /// implementation
        /// </summary>
        /// <param name="f">calling sequence: f(x,y,z,t)</param>
        /// <param name="time"></param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this Func<double[], double, double> f, double time) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                double[] X = new double[D];

                for (int i = 0; i < inp.GetLength(0); i++) {
                    for (int d = 0; d < D; d++)
                        X[d] = inp[i, d];

                    res[i] = f(X, time);
                }
            });
        }

        /// <summary>
        /// Vectorized 3D function (<see cref="ScalarFunction"/>) from a scalar implementation, with fixed time
        /// </summary>
        /// <param name="f">calling sequence: f(<paramref name="time"/>,x,y,z)</param>
        /// <param name="time">fixed time</param>
        /// <returns></returns>
        public static ScalarFunction Vectorize(this _4D f, double time) {
            return (delegate(MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                if (D != 3)
                    throw new ArgumentException("wrong spatial dimension.");

                for (int i = 0; i < inp.GetLength(0); i++) {
                    double x = inp[i, 0];
                    double y = inp[i, 1];
                    double z = inp[i, 2];

                    res[i] = f(time, x, y, z);
                }
            });
        }
    }
}
