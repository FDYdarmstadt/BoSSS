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

using BoSSS.Platform;
using ilPSP;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// This interface defines which methods any level set must
    /// provide to be used by BoSSS. The level-set itself may be calculated
    /// by BoSSS or come from some external object which implements this
    /// interface.
    /// </summary>
    public interface ILevelSet {

        /// <summary>
        /// Evaluates the level set 
        /// at nodes <paramref name="NS"/>
        /// in all cells from <paramref name="j0"/> to
        /// <paramref name="j0"/>+<paramref name="Len"/>-1.
        /// </summary>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="result">
        /// Output: 
        /// On exit, the value of the level set.
        /// The array is 2-dimensional:
        /// <list type="bullet">
        ///   <item>1st index: cell index <i>j</i> - <paramref name="j0"/>;</item>
        ///   <item>2nd index: node index <i>k</i>, corresponds with 1st index of the node set;</item>
        /// </list>
        /// </param>
        /// <param name="NS">
        /// nodes at wihich to evaluate.
        /// </param>
        /// <remarks>
        /// This method is vectorized: Here, it means that the Points at which the level set should be evaluated,
        /// are given for one cell in reference coordinates, but
        /// the evaluation is performed for <paramref name="Len"/> cells at once.
        /// </remarks>
        void Evaluate(int j0, int Len, NodeSet NS, MultidimensionalArray result);

        /// <summary>
        /// Evaluates the gradient of the level set 
        /// at nodes  <paramref name="NS"/>
        /// in all cells from <paramref name="j0"/> to
        /// <paramref name="j0"/>+<paramref name="Len"/>-1.
        /// </summary>
        /// <remarks>
        /// The returned vectors don't need to be of unit length;
        /// </remarks>
        void EvaluateGradient(int j0, int Len, NodeSet NS, MultidimensionalArray result);

        /// <summary>
        /// Evaluates all 2nd derivatives (by cell-local analytic derivation of the basis polynomials) of this field;
        /// this method my move to the <see cref="DGField"/>-class in future;
        /// </summary>
        /// <param name="j0"></param>
        /// <param name="Len"></param>
        /// <param name="NS"></param>
        /// <param name="result">
        /// <list type="bullet">
        ///   <item>1st index: cell index <em>j</em></item>
        ///   <item>2nd index: node index <em>m</em> into nodeset #<paramref name="NodeSetIndex"/></item>
        ///   <item>3rd index: spatial direction of 1st derivation, <em>k</em></item>
        ///   <item>4th index: spatial direction of 2nd derivation, <em>l</em></item>
        /// </list>
        /// So, the entry [j,m,k,l] = \f$ \frac{\partial}{\partial x_k} \frac{\partial}{\partial x_l} \varphi (\vec{\xi}_m)\f$ 
        /// where \f$ \vec{xi}_m\f$  is the <em>m</em>-th vector in the nodeset #<paramref name="NodeSetIndex"/>,
        /// in the <em>j</em>-th cell.
        /// </param>
        /// <remarks>
        /// Because of 2 derivatives taken, this field needs to be at least of DG degree 2 to get a non-zero result
        /// from this method.
        /// </remarks>
        void EvaluateHessian(int j0, int Len, NodeSet NS, MultidimensionalArray result);

        /// <summary>
        /// Evaluates the total curvature of this field, i.e.
        /// \f$  \mathrm{div} \left( \frac{ \varphi }{ \nabla_h \varphi } \right) \f$.
        /// </summary>
        /// <param name="j0"></param>
        /// <param name="Len"></param>
        /// <param name="NS"></param>
        /// <param name="result">
        /// <list type="bullet">
        ///   <item>1st index: cell index <em>j</em></item>
        ///   <item>2nd index: node index <em>m</em> into nodeset #<paramref name="NodeSetIndex"/></item>
        ///   <item>3rd index: spatial direction of 1st derivation, <em>k</em></item>
        ///   <item>4th index: spatial direction of 2nd derivation, <em>l</em></item>
        /// </list>
        /// So, the entry [j,m,k,l] = \f$ \frac{\partial}{\partial x_k} \frac{\partial}{\partial x_l} \varphi (\vec{\xi}_m)\f$ 
        /// where \f$ \vec{xi}_m\f$  is the <em>m</em>-th vector in the nodeset #<paramref name="NodeSetIndex"/>,
        /// in the <em>j</em>-th cell.
        /// </param>
        /// <remarks>
        /// Because of 2 derivatives taken, this field needs to be at least of DG degree 2 to get a non-zero result
        /// from this method.
        /// </remarks>
        void EvaluateTotalCurvature(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result);
    }
}
