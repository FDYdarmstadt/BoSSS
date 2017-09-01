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
using ilPSP;

namespace BoSSS.Foundation {
    
    /// <summary>
    /// Evaluates the scalar function
    /// at points specified in node set <paramref name="nodes"/>
    /// in all cells from <paramref name="j0"/> to
    /// <paramref name="j0"/>+<paramref name="Len"/>-1 and writes them
    /// to <paramref name="result"/>.
    /// </summary>
    /// <param name="j0">local index of the first cell to evaluate</param>
    /// <param name="Len">Number of cells to evaluate</param>
    /// <param name="result">
    /// The output: 
    /// On exit, the value of the DG field at the given nodes are saved there.<br/>
    /// The array is 2-dimensional:
    /// <list type="bullet">
    ///   <item>1st index: cell index <i>j</i> - <paramref name="j0"/>;</item>
    ///   <item>2nd index: node index <i>k</i></item>
    /// </list>
    /// </param>
    /// <param name="nodes">
    /// nodes in reference coordinates, at which the scalar function should be evaluated.
    /// </param>
    /// <remarks>
    /// This method is vectorized: Here, it means that the Points at which the DG field should be evaluated,
    /// are given for one cell in reference coordinates, but
    /// the evaluation is performed for <paramref name="Len"/> cells at once.
    /// </remarks>
    public delegate void ScalarFunctionEx(int j0, int Len, NodeSet nodes, MultidimensionalArray result);
}
