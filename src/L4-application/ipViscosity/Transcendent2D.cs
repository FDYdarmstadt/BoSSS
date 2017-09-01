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

namespace BoSSS.Application.ipViscosity {
    

    class Transcendent2D : TestSolution {

        public Transcendent2D() {
            base.D = 2;
        }

        public override double U(int d, double[] X) {
            double x = X[0], y = X[1];
            switch (d) {
                case 0:
                return Math.Sin(x) + Math.Cos( y);;

                case 1:
                return Math.Exp(-x*x)*y;

                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        public override double dU(int d, double[] X, int i) {
            double x = X[0], y = X[1];
            switch (d) {
                case 0:
                switch (i) {
                    case 0: return Math.Cos(x);
                    case 1: return -Math.Sin(y);
                    default: throw new ArgumentOutOfRangeException();
                }

                case 1:
                switch (i) {
                    case 0: return -2*x*Math.Exp(-x*x)*y;
                    case 1: return Math.Exp(-x*x);
                    default: throw new ArgumentOutOfRangeException();
                }

                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        public override double ddU(int d, double[] X, int i, int j) {
            double x = X[0], y = X[1];
            switch (d) {
                case 0:
                switch (i) {
                    case 0:
                    switch (j) {
                        case 0: return -Math.Sin(x);
                        case 1: return 0;
                        default: throw new ArgumentOutOfRangeException();
                    }

                    case 1:
                    switch (j) {
                        case 0: return 0;
                        case 1: return -Math.Cos(y);
                        default: throw new ArgumentOutOfRangeException();
                    }
                    default: throw new ArgumentOutOfRangeException();
                }

                case 1:
                switch (i) {
                    case 0:
                    switch (j) {
                        case 0: return -2*Math.Exp(-x*x)*y+4*x*x*Math.Exp(-x*x)*y;
                        case 1: return -2*x*Math.Exp(-x*x);
                        default: throw new ArgumentOutOfRangeException();
                    }

                    case 1:
                    switch (j) {
                        case 0: return -2*x*Math.Exp(-x*x);
                        case 1: return 0;
                        default: throw new ArgumentOutOfRangeException();
                    }

                    default: throw new ArgumentOutOfRangeException();
                }

                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        public override double mu(double[] X) {
            double x = X[0], y = X[1];

            return 1.5 + Math.Cos(x)*Math.Sin(y);
        }

        public override double dmu(double[] X, int i) {
            double x = X[0], y = X[1];

            switch (i) {
                case 0:
                return -Math.Sin(x)*Math.Sin(y);

                case 1:
                return Math.Cos(x)*Math.Cos(y);

                default:
                throw new ArgumentOutOfRangeException();
            }
        }
    }

}
