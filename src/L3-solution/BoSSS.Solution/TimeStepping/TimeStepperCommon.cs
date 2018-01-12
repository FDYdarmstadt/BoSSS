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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;


namespace BoSSS.Solution.TimeStepping {
    /// <summary>
    /// A helper class 
    /// </summary>
    public static class TimeStepperCommon {

        /// <summary>
        /// Check whether Operator, Fields and Parameters fit toghether
        /// </summary>
        /// <param name="spatialOp"></param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters"></param>
        public static void VerifyInput(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters) {
            if (!spatialOp.ContainsNonlinear && !(spatialOp.ContainsLinear()))
                throw new ArgumentException("spatial differential operator seems to contain no components.", "spatialOp");
            if (spatialOp.DomainVar.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("spatial differential operator must have the same number of domain and codomain variables.", "spatialOp");
            if (Fieldsmap.Fields.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("the number of fields in the coordinate mapping must be equal to the number of domain/codomain variables of the spatial differential operator", "fields");
            if (Parameters == null) {
                if (spatialOp.ParameterVar.Count != 0)
                    throw new ArgumentException("the number of fields in the parameter mapping must be equal to the number of parameter variables of the spatial differential operator", "Parameters");
            }
            else {
                if (Parameters.Fields.Count != spatialOp.ParameterVar.Count)
                    throw new ArgumentException("the number of fields in the parameter mapping must be equal to the number of parameter variables of the spatial differential operator", "Parameters");
            }
        }

    }
}
