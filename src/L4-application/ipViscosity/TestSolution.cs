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
using ilPSP.Utils;

namespace BoSSS.Application.ipViscosity {
    public abstract class TestSolution {

        protected int D;

        /// <summary>
        /// \f$ 
        /// u_d(\vec{X})
        /// \f$ 
        /// </summary>
        public abstract double U(int d, double[] X);

        /// \f$ 
        /// \partial_{i} u_d(\vec{X})
        /// \f$ 
        public abstract double dU(int d, double[] X, int i);

        /// \f$ 
        /// \partial_{i} \partial_{j} u_d(\vec{X})
        /// \f$ 
        public abstract double ddU(int d, double[] X, int i, int j);

        /// <summary>
        /// \f$ 
        /// \mu(\vec{X})
        /// \f$ 
        /// </summary>
        public abstract double mu(double[] X);

        /// <summary>
        /// \f$ 
        /// \partial_{i} \mu(\vec{X})
        /// \f$ 
        /// </summary>
        public abstract double dmu(double[] X, int i);

        /// <summary>
        /// \f[ 
        ///   - \operatorname{div} \left( \mu \nabla u_d \right)
        /// \f]
        /// </summary>
        public double Term1(int d, double[] X) {
            if (d < 0 || d >= D)
                throw new ArgumentOutOfRangeException();
            if (X.Length != D)
                throw new ArgumentException();

            double _mu = this.mu(X);
            double[] _dmu = new double[D];
            double[] _dU = new double[D];
            double _ddU = 0;

            for (int j = 0; j < D; j++) {
                _dmu[j] = this.dmu(X, j);
                _dU[j] = this.dU(d, X, j);
                _ddU += this.ddU(d, X, j, j);
            }

            double ret = 0;
            ret += GenericBlas.InnerProd(_dmu, _dU);
            ret += _mu*_ddU;

            return -ret;            
        }

        /// <summary>
        /// \f[ 
        ///   - \operatorname{div} \left( \mu (\partial_d \vec{u}) \right)
        /// \f]
        /// </summary>
        public double Term2(int d, double[] X) {
            if (d < 0 || d >= D)
                throw new ArgumentOutOfRangeException();
            if (X.Length != D)
                throw new ArgumentException();

            double _mu = this.mu(X);
            double[] _dmu = new double[D];
            double[] _dU = new double[D];
            double _ddU = 0;


            for (int j = 0; j < D; j++) {
                _dmu[j] = this.dmu(X, j);
                _dU[j] = this.dU(j, X, d);                
                _ddU += this.ddU(j, X, d, j);
            }

            double ret = 0;
            ret += GenericBlas.InnerProd(_dmu, _dU);
            ret += _mu*_ddU;

            return -ret;
        }

        /// <summary>
        /// \f[ 
        ///   \frac{2}{3} \operatorname{div} \left( \mu \myMatrix{I} \operatorname{div} ( \vec{u} )  \right)
        /// \f]
        /// </summary>
        public double Term3(int d, double[] X) {
            if (d < 0 || d >= D)
                throw new ArgumentOutOfRangeException();
            if (X.Length != D)
                throw new ArgumentException();

            double _mu = this.mu(X);
            double _dmu = this.dmu(X, d); 
            double _dU = 0;
            double _ddU = 0;


            for (int j = 0; j < D; j++) {
                _dU += this.dU(j, X, j);
                _ddU += this.ddU(j, X, d, j);
            }

            double ret = 0;
            ret += _dmu*_dU;
            ret += _mu*_ddU;

            return (2.0/3.0)*ret;
        }
    }
}
