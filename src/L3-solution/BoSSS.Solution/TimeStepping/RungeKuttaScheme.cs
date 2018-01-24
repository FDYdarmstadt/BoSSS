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
using System.Linq;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Structure for storing the Butcher Tableau of a Runge-Kutta Method;
    /// The Runge-Kutta Method for the ODE 
    /// \f[
    ///    y' = f(t,y), y(0) = y0 
    /// \f] 
    /// is defined as:
    /// \f{align*}{
    /// k[0] & = f(0,y0), \\
    /// k[r] & = f(c[r],y0 + k[0]*a[r,0] + ... + k[r-1]*a[r,r-1]); \\
    /// y1   & = y0 + (b[0]*k[0] + ... + b[s-1]*k[s-1]);
    /// \f}
    /// </summary>
    public class RungeKuttaScheme : ICloneable {

        /// <summary>
        /// See class summary.
        /// </summary>
        public double[,] a;

        /// <summary>
        /// See class summary.
        /// </summary>
        public double[] b;

        /// <summary>
        /// See class summary.
        /// </summary>
        public double[] c;

        /// <summary>
        /// performs some basic consistency checks on this Runge-Kutta rule.
        /// </summary>
        public void Verify() {
            if (a == null)
                throw new ArgumentNullException();
            if (b == null)
                throw new ArgumentNullException();
            if (c == null)
                throw new ArgumentNullException();


            int S = c.Length;
            if (a.GetLength(0) != a.GetLength(1))
                throw new ArgumentException();
            if (a.GetLength(0) != S)
                throw new ArgumentException();
            if (b.Length != S)
                throw new ArgumentException();

            if (Math.Abs(b.Sum() - 1.0) > 1.0e-12)
                throw new ArgumentException();

            for (int s = 0; s < S; s++) {

                if (Math.Abs(this.a.GetRow(s).Sum() - c[s]) > 1.0e-12)
                    throw new ArgumentException();

            }
        }

        /// <summary>
        /// Guess what?
        /// </summary>
        public object Clone() {
            var C = new RungeKuttaScheme();
            C.a = this.a.CloneAs();
            C.b = this.b.CloneAs();
            C.c = this.c.CloneAs();
            return C;
        }

        /// <summary>
        /// True, if this scheme is explicit.
        /// </summary>
        public bool IsExplicit {
            get {
                this.Verify();
                int S = this.Stages;
                for (int s = 0; s < S; s++) {
                    for (int k = s; k < S; k++)
                        if (a[s, k] != 0.0)
                            return false;
                }
                return true;
            }
        }

        /// <summary>
        /// True, if this scheme is a diagonally implicit Runge Kutta (DIRK) formula.
        /// </summary>
        public bool IsDiagonallyImplicit {
            get {
                this.Verify();
                int S = this.Stages;
                for (int s = 0; s < S; s++) {
                    for (int k = s + 1; k < S; k++)
                        if (a[s, k] != 0.0)
                            return false;
                }
                return true;
            }
        }

        /// <summary>
        /// number of stages of this Runge-Kutta Scheme
        /// </summary>
        public int Stages {
            get {
                return c.Length;
            }
        }

        /// <summary>
        /// 3rd order TVD scheme;<br/>
        /// </summary>
        /// <remarks>
        /// From Mathematics of Computation, Vol. 67, No 221, 1998, Pages 73-85:<br/>
        /// Total Variation diminishing Runge-Kutta schemes,<br/>
        /// by Sigal Gottlieb and Chi-Wang Shu;
        /// </remarks>
        public static RungeKuttaScheme TVD3 {
            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[3, 3];
                ret.b = new double[3];
                ret.c = new double[3];

                ret.a[1, 0] = 1;
                ret.a[2, 0] = 0.25;
                ret.a[2, 1] = 0.25;

                ret.c[0] = 0;
                ret.c[1] = 1;
                ret.c[2] = 0.5;

                ret.b[0] = 1.0 / 6.0;
                ret.b[1] = 1.0 / 6.0;
                ret.b[2] = 2.0 / 3.0;

                return ret;
            }
        }

        /// <summary>
        /// Middlepoint - Rule, Order 2, 2 Stages;
        /// </summary>
        public static RungeKuttaScheme Middlepoint {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[2, 2];
                ret.b = new double[2];
                ret.c = new double[2];

                ret.a[1, 0] = 0.5;

                ret.c[0] = 0;
                ret.c[1] = 0.5;

                ret.b[0] = 0;
                ret.b[1] = 1;

                return ret;
            }
        }

        /// <summary>
        /// Explicit Euler - Rule, Order 1, 2 Stages, for testing purpose;
        /// Equivalent of doing explicit Euler two times with the half time-step;
        /// </summary>
        public static RungeKuttaScheme ExplicitEuler2 {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[2, 2];
                ret.b = new double[2];
                ret.c = new double[2];

                ret.a[1, 0] = 0.5;

                ret.c[0] = 0;
                ret.c[1] = 0.5;

                ret.b[0] = 0.5;
                ret.b[1] = 0.5;


                return ret;
            }
        }

        /// <summary>
        /// Explicit Euler - Rule, Order 1, 1 Stage;
        /// </summary>
        public static RungeKuttaScheme ExplicitEuler {
            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[1, 1];
                ret.b = new double[1];
                ret.c = new double[1];


                ret.c[0] = 0;

                ret.b[0] = 1;


                return ret;
            }
        }

        /// <summary>
        /// classical Method of Runge and Kutta, anno 1901, order 4, 4 stages;
        /// </summary>
        public static RungeKuttaScheme RungeKutta1901 {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[4, 4];
                ret.b = new double[4];
                ret.c = new double[4];

                ret.a[1, 0] = 0.5;
                ret.a[2, 1] = 0.5;
                ret.a[3, 2] = 1;

                ret.c[0] = 0;
                ret.c[1] = 0.5;
                ret.c[2] = 0.5;
                ret.c[3] = 1;

                ret.b[0] = (1.0 / 6.0);
                ret.b[1] = (2.0 / 6.0);
                ret.b[2] = (2.0 / 6.0);
                ret.b[3] = (1.0 / 6.0);

                return ret;
            }
        }

        /// <summary>
        /// 3/8 - Rule, order 4, 4 stages
        /// </summary>
        public static RungeKuttaScheme ThreeOverEight {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[4, 4];
                ret.b = new double[4];
                ret.c = new double[4];

                ret.a[1, 0] = 1.0 / 3.0;
                ret.a[2, 0] = -1.0 / 3.0;
                ret.a[2, 1] = 1.0;
                ret.a[3, 0] = 1;
                ret.a[3, 1] = -1.0;
                ret.a[3, 2] = 1.0;

                ret.c[0] = 0;
                ret.c[1] = 1.0 / 3.0;
                ret.c[2] = 2.0 / 3.0;
                ret.c[3] = 1.0;

                ret.b[0] = 1.0 / 8.0;
                ret.b[1] = 3.0 / 8.0;
                ret.b[2] = 3.0 / 8.0;
                ret.b[3] = 1.0 / 8.0;

                return ret;
            }
        }

        /// <summary>
        /// Heun - Rule, order 3, 3 stages
        /// </summary>
        public static RungeKuttaScheme Heun {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[3, 3];
                ret.b = new double[3];
                ret.c = new double[3];

                ret.a[1, 0] = 1.0 / 3.0;
                ret.a[2, 1] = 2.0 / 3.0;

                ret.c[0] = 0;
                ret.c[1] = 1.0 / 3.0;
                ret.c[2] = 2.0 / 3.0;

                ret.b[0] = 1.0 / 4.0;
                ret.b[1] = 0;
                ret.b[2] = 3.0 / 4.0;

                return ret;
            }
        }

        /// <summary>
        /// Runge Kutta - Rule, order 3, 3 stages
        /// </summary>
        public static RungeKuttaScheme RungeKutta3 {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[3, 3];
                ret.b = new double[3];
                ret.c = new double[3];

                ret.a[1, 0] = 1.0 / 2.0;
                ret.a[2, 0] = -1.0;
                ret.a[2, 1] = 2.0;

                ret.c[0] = 0;
                ret.c[1] = 1.0 / 2.0;
                ret.c[2] = 1.0;

                ret.b[0] = 1.0 / 6.0;
                ret.b[1] = 2.0 / 3.0;
                ret.b[2] = 1.0 / 6.0;

                return ret;
            }
        }


        /// <summary>
        /// Heun - Rule, order 2, 2 stages
        /// </summary>
        public static RungeKuttaScheme Heun2 {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[2, 2];
                ret.b = new double[2];
                ret.c = new double[2];

                ret.a[1, 0] = 1.0;

                ret.c[0] = 0;
                ret.c[1] = 1.0;

                ret.b[0] = 0.5;
                ret.b[1] = 0.5;

                return ret;
            }
        }


        /// <summary>
        /// Kraaijevanger, 1991: Contractivity of Runge-Kutta methods
        /// Strong stability-preserving 5-stage scheme of 4th order
        /// </summary>
        public static RungeKuttaScheme SSP54 {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[5, 5];
                ret.b = new double[5];
                ret.c = new double[5];

                ret.a[1, 0] = 0.39175222657188905833;
                ret.a[2, 0] = 0.21766909626116921036;
                ret.a[2, 1] = 0.36841059305037202075;
                ret.a[3, 0] = 0.08269208665781075441;
                ret.a[3, 1] = 0.13995850219189573938;
                ret.a[3, 2] = 0.25189177427169263984;
                ret.a[4, 0] = 0.06796628363711496324;
                ret.a[4, 1] = 0.11503469850463199467;
                ret.a[4, 2] = 0.20703489859738471851;
                ret.a[4, 3] = 0.54497475022851992204;

                ret.c[0] = 0.0;
                ret.c[1] = 0.39175222657188905833;
                ret.c[2] = 0.58607968931154123111;
                ret.c[3] = 0.47454236312139913362;
                ret.c[4] = 0.93501063096765159845;

                ret.b[0] = 0.14681187608478644956;
                ret.b[1] = 0.24848290944497614757;
                ret.b[2] = 0.10425883033198029567;
                ret.b[3] = 0.27443890090134945681;
                ret.b[4] = 0.22600748323690765039;

                return ret;
            }
        }

        /// <summary>
        /// Toulorge & Desmet 2012: Optimal Runge–Kutta schemes for
        /// discontinuous Galerkin space discretizations applied to wave
        /// propagation problems
        /// Scheme optimized towards size of the stability region and
        /// accuracy
        /// </summary>
        public static RungeKuttaScheme RKC84 {

            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();
                int numberOfStages = 8;

                ret.a = new double[numberOfStages, numberOfStages];
                ret.b = new double[numberOfStages];
                ret.c = new double[numberOfStages];

                //Low - storage RK data; Table A.21
                double[] A = new double[] {
                        0.0,
                        -0.7212962482279240,
                        -0.01077336571612980,
                        -0.5162584698930970,
                        -1.730100286632201,
                        -5.200129304403076,
                        0.7837058945416420,
                        -0.5445836094332190
                    };

                double[] B = new double[] {
                        0.2165936736758085,
                        0.1773950826411583,
                        0.01802538611623290,
                        0.08473476372541490,
                        0.8129106974622483,
                        1.903416030422760,
                        0.1314841743399048,
                        0.2082583170674149
                    };

                ret.c[0] = 0.0;
                ret.c[1] = 0.2165936736758085;
                ret.c[2] = 0.2660343487538170;
                ret.c[3] = 0.2840056122522720;
                ret.c[4] = 0.3251266843788570;
                ret.c[5] = 0.4555149599187530;
                ret.c[6] = 0.7713219317101170;
                ret.c[7] = 0.9199028964538660;

                // Recurrence relation (16)
                //for (int i = 0; i < numberOfStages - 2; i++) {
                for (int i = 0; i < numberOfStages - 1; i++) {
                    ret.a[i + 1, i] = B[i];
                }

                ret.b[numberOfStages - 1] = B[numberOfStages - 1];

                for (int i = numberOfStages - 1; i >= 1; i--) {
                    ret.b[i - 1] = A[i] * ret.b[i] + B[i - 1];
                }

                // This is how I read the recurrence relation, but this
                // gives inconsistent sums of a[:,i], so I used the code
                // below which seems to do the trick
                //for (int i = 1; i < numberOfStages - 1; i++) {
                //    ret.a[i + 1, i - 1] = A[i] * B[i] + ret.c[i];
                //}

                for (int i = 1; i < numberOfStages - 1; i++) {
                    // Some per row must be c and there are only two
                    // entries per row (low storage scheme!)
                    ret.a[i + 1, i - 1] = ret.c[i + 1] - ret.a[i + 1, i];
                }

                return ret;
            }
        }


        /// <summary>
        /// The Crank-Nicolson (aka. Implicit Trapezoidal)
        /// scheme written as diagonally implicit Runge-Kutta method.
        /// </summary>
        public static RungeKuttaScheme CrankNicolson {
            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[2, 2];
                ret.b = new double[2];
                ret.c = new double[2];

                ret.a[1, 0] = 0.5;
                ret.a[1, 1] = 0.5;

                ret.c[0] = 0;
                ret.c[1] = 1.0;

                ret.b[0] = 0.5;
                ret.b[1] = 0.5;

                return ret;
            }
        }


        /// <summary>
        /// The implicit Euler 
        /// scheme written as diagonally implicit Runge-Kutta method.
        /// </summary>
        public static RungeKuttaScheme ImplicitEuler {
            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[1, 1];
                ret.b = new double[1];
                ret.c = new double[1];

                ret.a[0, 0] = 1.0;

                ret.c[0] = 1;

                ret.b[0] = 1;

                return ret;
            }
        }

        /// <summary>
        /// Diagonally implicit scheme, 3rd order.
        /// </summary>
        /// <remarks>
        /// See:
        /// Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations,
        /// Uri M. Ascher and Steven J. Ruuth and Raymond J. Spiteri,
        /// Applied Numerical Mathematics, volume 25, number 2, pages 151 - 167, year 1997, 
        /// </remarks>
        public static RungeKuttaScheme IMEX3 {
            get {
                RungeKuttaScheme ret = new RungeKuttaScheme();

                ret.a = new double[3, 3];
                ret.b = new double[3];
                ret.c = new double[3];

                ret.c[0] = 0.435866521508458999416019451194;
                ret.a[0, 0] = 0.435866521508458999416019451194;

                ret.c[1] = 0.717933260754229499708009725597;
                ret.a[1, 0] = 0.282066739245770500291990274403;
                ret.a[1, 1] = 0.435866521508458999416019451194;

                ret.c[2] = 1.0;
                ret.a[2, 0] = 1.20849664917601007033647768407;
                ret.a[2, 1] = -0.64436317068446906975249713526;
                ret.a[2, 2] = 0.435866521508458999416019451194;

                ret.b[0] = 1.20849664917601007033647768407;
                ret.b[1] = -0.64436317068446906975249713526;
                ret.b[2] = 0.435866521508458999416019451194;

                return ret;
            }
        }

    }
}
