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
using System.Collections;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Utils;
using System.Diagnostics;
using System.Linq;
using ilPSP;
using BoSSS.Foundation.Comm;
using System.Collections.Generic;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.XDG {

    partial class LevelSetTracker {



        public Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper> m_QuadFactoryHelpers
            = new Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper>();

        /*
        /// <summary>
        /// Central 'factory' for creating Level Set - related quadrature.
        /// </summary>
        /// <remarks>
        /// The centralized approach should avoid multiple creation of the same quadrature rule.
        /// </remarks>
        public XQuadFactoryHelper GetXQuadFactoryHelper(XQuadFactoryHelper.MomentFittingVariants variant) {
            if (!m_QuadFactoryHelpers.ContainsKey(variant)) {
                m_QuadFactoryHelpers[variant] = new XQuadFactoryHelper(this, variant);
            }
            return m_QuadFactoryHelpers[variant];
        }
        
      
        */


        public XDGSpaceMetrics GetXDGSpaceMetrics(XQuadFactoryHelper.MomentFittingVariants variant, int quadorder, int stackindex) {
            if(!m_QuadFactoryHelpers.ContainsKey(variant)) {
                m_QuadFactoryHelpers[variant] = new XQuadFactoryHelper(this, variant);
            }


        }



        /// <summary>
        /// Base class for all types of histories,
        /// i.e. scalar histories (constant in space)
        /// and scalar- and vector-fields histories (not constant in space).
        /// </summary>
        /// <typeparam name="T"></typeparam>
        sealed public class HistoryStack<T> {

            /// <summary>
            /// ctor
            /// </summary>
            internal HistoryStack(T curr) {
                if(!(curr.GetType().IsValueType) && (curr == null))
                    throw new ArgumentNullException();

                m_Current = curr;
            }


            /// <summary>
            /// <see cref="PushCount"/>
            /// </summary>
            protected int m_PushCount = 0;

            /// <summary>
            /// counts every call to <see cref="Push"/>
            /// </summary>
            public int PushCount {
                get {
                    return m_PushCount;
                }
            }

            /// <summary>
            /// either a double for time-dependent scalars, which are constant in space (e.g. ambient pressure in Low-Mach solver)
            /// or a <see cref="BoSSS.Foundation.DGField"/> or a <see cref="BoSSS.Foundation.VectorField{T}"/>.
            /// </summary>
            protected T m_Current;

            /// <summary>
            /// see <see cref="HistoryLength"/>
            /// </summary>
            protected int m_HistoryLength = 0;

            /// <summary>
            /// the number of Timesteps that is stored <em>in addition</em> to the current one
            /// </summary>
            public int HistoryLength {
                get { return m_HistoryLength; }
            }

            /// <summary>
            /// returns the minimum of the two numbers:
            /// <list type="bullet">
            ///   <item>number of <see cref="Push"/> - operations carried out on this object</item>
            ///   <item><see cref="m_HistoryLength"/></item>
            /// </list>
            /// If this value is, e.g. 2, this means that timesteps 1, 0, -1 are available;
            /// </summary>
            /// <returns></returns>
            public int GetPopulatedLength() {
                return History.Count;
            }

            /// <summary>
            /// container for previous values: entries are either of type <see cref="DGField"/> or <see cref="VectorField{t}"/>
            /// </summary>
            protected ArrayList History = new ArrayList();

            /// <summary>
            /// pushes a new variable set onto the top of the stack
            /// </summary>
            /// <remarks>
            /// An example of the Push-operation, for <see cref="HistoryLength"/>==3:
            /// <code>
            /// 
            ///    Timestep Index ->  1  0 -1 -2 
            ///   -----------------------------------
            ///    Before Push(..): | a  b  c  d 
            ///                     | 
            ///    After Push(..):  | a  a  b  c
            /// </code>
            /// Because of this behavior, any difference norm (sometimes incorrectly called 'residual')
            /// should be computed prior to the <see cref="Push"/>-operation.
            /// </remarks>
            abstract public void Push();



            /// <summary>
            /// Gets access to previous (and actual) timesteps.
            /// </summary>
            /// <param name="i">
            /// Either 1 (equal to <see cref="Current"/>) or a non-positive value:
            /// the first item in the stack is 0, the second item is -1, etc..
            /// </param>
            public T this[int i] {
                get {
                    if(i > 1)
                        throw new IndexOutOfRangeException();

                    if(i == 1)
                        return (T)m_Current;
                    else
                        return (T)History[-i];
                }
            }

            /// <summary>
            /// The top item of the stack.
            /// </summary>
            public T Current {
                get {
                    return this[1];
                }
            }

        }

    }
}
