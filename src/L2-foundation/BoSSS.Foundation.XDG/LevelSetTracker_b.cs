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
using MPI.Wrappers;

namespace BoSSS.Foundation.XDG {

    partial class LevelSetTracker {



        HistoryStack<Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper>> m_QuadFactoryHelpersHistory = null;
        
        /// <summary>
        /// Central 'factory' for creating Level Set - related quadrature.
        /// </summary>
        /// <remarks>
        /// The centralized approach should avoid multiple creation of the same quadrature rule.
        /// </remarks>
        XQuadFactoryHelper GetXQuadFactoryHelper(XQuadFactoryHelper.MomentFittingVariants variant, int HistoryIndex = 1) {
            var dict = m_QuadFactoryHelpersHistory[HistoryIndex];
            
            if (!dict.ContainsKey(variant)) {
                dict[variant] = new XQuadFactoryHelper(
                    this.DataHistories.Select(hist => hist[HistoryIndex]).ToArray(),
                    variant);
            }

            return dict[variant];
        }
        
      
        HistoryStack<Dictionary<Tuple<SpeciesId[], XQuadFactoryHelper.MomentFittingVariants, int>, XDGSpaceMetrics>> m_XDGSpaceMetricsHistory = null;
        
        Dictionary<Tuple<SpeciesId[], XQuadFactoryHelper.MomentFittingVariants, int>, XDGSpaceMetrics> NewXDGSpaceMetricsCache() {
            return new Dictionary<Tuple<SpeciesId[], XQuadFactoryHelper.MomentFittingVariants, int>, XDGSpaceMetrics>(
                new FuncEqualityComparer<Tuple<SpeciesId[], XQuadFactoryHelper.MomentFittingVariants, int>>(
                    delegate(Tuple<SpeciesId[], XQuadFactoryHelper.MomentFittingVariants, int> A, Tuple<SpeciesId[], XQuadFactoryHelper.MomentFittingVariants, int> B) {
                        //if(.)
                        if(!ArrayTools.ListEquals(A.Item1, B.Item1, (a, b) => a.cntnt == b.cntnt))
                            return false;
                        if(A.Item2 != B.Item2)
                            return false;
                        if(A.Item3 != B.Item3)
                            return false;

                        return true;
                    }));

        }

        /// <summary>
        /// The order of all cut-cell quadrature rules which are present in the cache.
        /// </summary>
        public ISet<int> GetCachedOrders() {
            HashSet<int> hs = new HashSet<int>();
            foreach(var t in m_XDGSpaceMetricsHistory.Current.Keys) {
                hs.Add(t.Item3);
            }

            return hs;
        }

        /// <summary>
        /// Cut Cell and Cut Edge metrics before agglomeration
        /// </summary>
        public XDGSpaceMetrics GetXDGSpaceMetrics(SpeciesId Spc, int CutCellsQuadOrder, int HistoryIndex = 1) {
            return GetXDGSpaceMetrics(new SpeciesId[] { Spc }, CutCellsQuadOrder, HistoryIndex);
        }


        /// <summary>
        /// Cut Cell and Cut Edge metrics before agglomeration
        /// </summary>
        public XDGSpaceMetrics GetXDGSpaceMetrics(IEnumerable<SpeciesId> Spc, int CutCellsQuadOrder, int HistoryIndex = 1) {
            //if(!m_QuadFactoryHelpers.ContainsKey(variant)) {
            //    m_QuadFactoryHelpers[variant] = new XQuadFactoryHelper(this, variant);
            //}
            var _Spc = Spc.ToArray();
#if TEST
            MPICollectiveWatchDog.WatchAtRelease();
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
#endif
            var CutCellsQuadType = this.CutCellQuadratureType;

            //throw new NotImplementedException("todo");
            var dict = m_XDGSpaceMetricsHistory[HistoryIndex];
            var key = Tuple.Create(_Spc, this.CutCellQuadratureType, CutCellsQuadOrder);
            if (!dict.ContainsKey(key)) {
                dict.Add(key,
                    new XDGSpaceMetrics(this,
                        GetXQuadFactoryHelper(CutCellsQuadType, HistoryIndex),
                        CutCellsQuadOrder,
                        _Spc,
                        HistoryIndex
                        )
                );
            }

            return dict[key];
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Spc">all species, for which agglomeration should be performed</param>
        /// <param name="CutCellsQuadOrder">
        /// cut-cell quadrature order for the quadrature rule that is used to determine cell volumes;
        /// this should typically be the same order which is used to evaluate the XDG operator matrix.
        /// </param>
        /// <param name="__AgglomerationTreshold">
        /// Volume fraction, which triggers cell agglomeration;
        /// see <see cref="MultiphaseCellAgglomerator.AgglomerationThreshold"/></param>
        /// <param name="oldTs__AgglomerationTreshold">
        /// Agglomeration thresholds for   _previous_ timesteps.
        /// The number of entries in this array determines how many previous timesteps are considered
        /// (<see cref="LevelSetTracker.HistoryLength"/>).
        /// </param>
        /// <param name="AgglomerateNewborn">
        /// 'Newborn' cut-cells 
        /// are agglomerated to cells which already exist in the previous  timestep.
        /// </param>
        /// <param name="AgglomerateDecased">
        /// 'Deceased' cut-cells are agglomerated to cells which exist in the next timestep.
        /// </param>
        /// <param name="ExceptionOnFailedAgglomeration">
        /// If true, an exception is thrown for 
        /// any cell which should be agglomerated, if no neighbour is found.
        /// </param>
        /// <param name="NewbornAndDecasedThreshold">
        /// Volume fraction threshold at which a cut-cell counts as newborn, resp. deceased, see <paramref name="AgglomerateNewborn"/>, <paramref name="AgglomerateDecased"/>;
        /// </param>       
        /// <param name="Tag">
        /// A string value to pass tags for debugs (e.g. LevelSetAgg)
        /// </param>     
        /// <returns></returns>
        public MultiphaseCellAgglomerator GetAgglomerator(
            SpeciesId[] Spc, int CutCellsQuadOrder,
            double __AgglomerationTreshold, 
            bool AgglomerateNewborn = false, bool AgglomerateDecased = false, bool ExceptionOnFailedAgglomeration = true, 
            double[] oldTs__AgglomerationTreshold = null,
            double NewbornAndDecasedThreshold = 1.0e-6,
            string Tag = null
            ) {
            MPICollectiveWatchDog.Watch(token: 169);
            return new MultiphaseCellAgglomerator(this, Spc, CutCellsQuadOrder,
                __AgglomerationTreshold, AgglomerateNewborn, AgglomerateDecased, ExceptionOnFailedAgglomeration, oldTs__AgglomerationTreshold, NewbornAndDecasedThreshold, Tag);
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
            int m_PushCount = 0;

            /// <summary>
            /// Counts every call to <see cref="Push"/>
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
            T m_Current;

            /// <summary>
            /// see <see cref="HistoryLength"/>
            /// </summary>
            int m_HistoryLength = 0;

            /// <summary>
            /// the number of Timesteps that should stored <em>in addition</em> to the current one, at maximum, before elements are dropped.
            /// </summary>
            public int HistoryLength {
                get {
                    return m_HistoryLength;
                }
                set {
                    m_HistoryLength = value;
                    if(History.Count > m_HistoryLength) {
                        History.RemoveRange(m_HistoryLength, History.Count - m_HistoryLength);
                        Debug.Assert(History.Count == m_HistoryLength);
                    }
                }
            }

            /// <summary>
            /// 1 (including) to -<see cref="GetPopulatedLength()"/> (including) 
            /// </summary>
            public int[] AvailableIndices {
                get {
                    var R = new int[GetPopulatedLength() + 1];
                    for(int i = 0; i < R.Length; i++) {
                        R[i] = 1 - i;
                    }
                    return R;
                }
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
            /// indexes of all available times, i.e. all valid indexes to access <see cref="this[int]"/>
            /// </summary>
            public int[] PopulatedIndices {
                get {
                    int[] ret = new int[GetPopulatedLength() + 1];
                    for(int i = 0; i < ret.Length; i++)
                        ret[i] = 1 - i;
                    return ret;
                }
            }



            /// <summary>
            /// Sontainer for previous states.
            /// </summary>
            List<T> History = new List<T>();

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
            /// </remarks>
            internal void Push(Func<T,T> Replicator1, Func<T,T,T> Replicator0) {
                if(History.Count < HistoryLength) {
                    History.Add(default(T));
                }
                if(History.Count > HistoryLength) {
                    History.RemoveRange(HistoryLength, History.Count - HistoryLength);
                    Debug.Assert(History.Count == HistoryLength);
                }

                for(int i = History.Count - 1; i >= 1; i--) {
                    History[i] = History[i - 1];
                }

                if(History.Count > 0) {
                    History[0] = Replicator0(m_Current, History[0]);
                }
                m_Current = Replicator1(m_Current);
                m_PushCount++;
            }

            internal T Pop(Func<T, T, T> Replacor1) {
                T ret = m_Current;
                m_Current = Replacor1(m_Current, History[0]);
                for(int i = 0; i <= History.Count - 2; i++) {
                    History[i] = History[i + 1];
                }
                History.RemoveAt(History.Count - 1);
                m_PushCount--;
                return ret;
            }



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
                    if(i <= -m_HistoryLength)
                        throw new IndexOutOfRangeException();
                    if(i <= -GetPopulatedLength())
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
