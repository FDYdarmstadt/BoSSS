﻿/* =======================================================================
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

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// thrown by <see cref="LevelSetTracker.UpdateTracker"/>;
    /// </summary>
    public class LevelSetTopologyException : Exception {

        static string ComposeMessage(IEnumerable<(int iLevSet, int j, int Neigh, int dist_j, int dist_neigh)> problems) {
            using (StringWriter stw = new StringWriter()) {
                stw.Write("failed on Level Set Topology: contact of purely positive/negative domain across an edge without a cut cell in between;");
                foreach (var tttt in problems) {
                    stw.Write($"cell {tttt.j}, dist = {tttt.dist_j}, distance of neighbor cell {tttt.Neigh} is {tttt.dist_neigh}; ");
                }

                return stw.ToString();
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        internal LevelSetTopologyException(IEnumerable<(int iLevSet, int j, int Neigh, int dist_j, int dist_neigh)> problems)
            : base(ComposeMessage(problems)) {
            m_Problems = problems.ToArray();
        }

        (int iLevSet, int j, int Neigh, int dist_j, int dist_neigh)[] m_Problems;
    }
}
