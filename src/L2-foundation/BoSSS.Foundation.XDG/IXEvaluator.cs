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

using ilPSP.LinSolvers;
using System.Collections.Generic;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Evaluation and matrix assembly of XSpatial operators
    /// </summary>
    public interface IXEvaluator : IEvaluator {

        /// <summary>
        /// The pointer to a owner object, which totally contradicts the original idea of object-orientation. Hehe.
        /// </summary>
        XSpatialOperatorMk2 Owner { get; }

       
        /// <summary>
        /// Operator coefficients for each species
        /// </summary>
        Dictionary<SpeciesId, CoefficientSet> SpeciesOperatorCoefficients { get; }
    }


    /// <summary>
    /// Explicit evaluation of  Xspatial operators, i.e. matrix assembly.
    /// </summary>
    public interface IXEvaluatorNonLin : IXEvaluator, IEvaluatorNonLin {

        

    }
   
    /// <summary>
    /// Evaluation of linear/linearized Xspatial operators, i.e. matrix assembly.
    /// </summary>
    public interface IXEvaluatorLinear : IXEvaluator, IEvaluatorLinear {

        
    }
    
}
