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

namespace CNS.LoadBalancing {

    /// <summary>
    /// Classifies cells by assigning them a performance class (a number
    /// greater than or equal to zero)
    /// </summary>
    public interface ICellClassifier {

        /// <summary>
        /// Sorts each cell updated by the local process into a performance
        /// class
        /// </summary>
        /// <param name="program"></param>
        /// <returns>
        /// The total number of performance classes (across all processes) and,
        /// for each cell, its assigned performance class
        /// </returns>
        (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program);
    }
}
