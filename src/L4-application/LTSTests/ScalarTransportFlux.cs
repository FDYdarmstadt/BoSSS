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
using BoSSS.Solution.Utils;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using System.Collections.Generic;

namespace LTSTests {
    /// <summary>
    /// Flux implementation for the LTS NUnit-test
    /// </summary>
    class ScalarTransportFlux : NonlinearFlux {

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" };
            }
        }

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            Vector2D n;
            n.x = normal[0];
            n.y = normal[1];
            if (n * FlowField(x) >= 0) {
                return (FlowField(x) * Uin[0]) * n;
            } else {
                return (FlowField(x) * Inflow(time)) * n;
            }
        }

        protected override void Flux(double time, double[] x, double[] U, double[] output) {
            Vector2D o;
            o = FlowField(x) * U[0];
            output[0] = o.x;
            output[1] = o.y;
        }

        protected override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            Vector2D n;
            n.x = normal[0];
            n.y = normal[1];
            if (FlowField(x) * n > 0)
                return (FlowField(x) * Uin[0]) * n;
            else
                return (FlowField(x) * Uout[0]) * n;
        }

        Vector2D FlowField(double[] x) {
            Vector2D u;
            u.x = 1;
            u.y = 0;
            return u;
        }

        double Inflow(double time) {
            return 0.0;
        }
    }
}
