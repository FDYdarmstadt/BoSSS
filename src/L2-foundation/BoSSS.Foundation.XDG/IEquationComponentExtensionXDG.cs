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

namespace BoSSS.Foundation.XDG {
    /// <summary>
    /// extension methods for the <see cref="IEquationComponent"/>-interface
    /// </summary>
    public static class IEquationComponentExtension {

        /// <summary>
        /// creates the spatial operator that consists only of component <paramref name="c"/>
        /// </summary>
        public static XSpatialOperator XOperator(this IEquationComponent c, int DegreeOfNonlinearity = 1) {
            return XOperator(c, QuadOrderFunc.NonLinear(DegreeOfNonlinearity));
        }


        /// <summary>
        /// creates the spatial operator that consists only of component <paramref name="c"/>
        /// </summary>
        public static XSpatialOperator XOperator(this IEquationComponent c, Func<int[], int[], int[], int> quadOrderFunc) {

            string[] Codomain = new string[] { "v1" };
            string[] Domain = c.ArgumentOrdering.ToArray();
            string[] Param = (c.ParameterOrdering != null) ? c.ParameterOrdering.ToArray() : new string[0];

            XSpatialOperator ret = new XSpatialOperator(Domain, Param, Codomain, quadOrderFunc);
            ret.EquationComponents[Codomain[0]].Add(c);
            ret.Commit();

            return ret;
        }

    }
}
