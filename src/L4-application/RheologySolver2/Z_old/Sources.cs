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
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Hard-coded implementation of source for testing manufactured solution
    /// </summary>

    class MomentumXSource : LinearSource {

        public MomentumXSource() {

        }

        public override IList<string> ArgumentOrdering
        {
            get { return new string[0]; }
        }

        protected override double Source(double[] x, double[] parameters, double[] U) {

            return -4 * Math.Exp(x[0]) * Math.Sin(x[1]);
        }
    }

    class MomentumYSource : LinearSource {


        public MomentumYSource() {

        }

        protected override double Source(double[] x, double[] parameters, double[] U) {
            //return -Math.Sin(x[1]) - 2.0 * eta * Math.Sin(x[1]) * Math.Cos(x[1]);
            //without ext. pressure
            return -2 * Math.Exp(x[0]) * (-Math.Cos(x[1]) + x[1] * Math.Sin(x[1]));
        }

        public override IList<string> ArgumentOrdering
        {
            get { return VariableNames.VelocityVector(2); }
        }
    }
}
