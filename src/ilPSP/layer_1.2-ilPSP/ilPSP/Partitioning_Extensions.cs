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

namespace ilPSP {

    /// <summary>
    /// Extension methods for <see cref="IPartitioning"/>.
    /// </summary>
    static public class Partitioning_Extensions {


        /// <summary>
        /// Throws an exception, if <see cref="IsInLocalRange(int)"/>(<paramref name="i"/>) evaluates to false;
        /// </summary>
        /// <param name="i"></param>
        static public void TestIfInLocalRange(this IPartitioning p, int i) {
            if (!p.IsInLocalRange(i)) {
                int i0 = p.i0;
                int iE = p.iE;
                throw new ArgumentOutOfRangeException("i", "index is not within local range of partition: expecting " + i0 + " <= i < " + iE + ", but got i = " + i + ";");
            }
        }

        /// <summary>
        /// True, if a block index is local.
        /// </summary>
        public static bool IsLocalBlock(this IBlockPartitioning p, int iBlock) {
            Debug.Assert(iBlock >= 0, "Block index out of range");
            Debug.Assert(iBlock < p.TotalNoOfBlocks, "Block index out of range");
            return (p.FirstBlock <= iBlock) && (iBlock < p.FirstBlock + p.LocalNoOfBlocks);
        }

        /// <summary>
        /// Subtracts <see cref="IPartitioning.i0"/> from <paramref name="iGlob"/>.
        /// </summary>
        static public int TransformIndexToLocal(this IPartitioning p, int iGlob) {
            return iGlob - p.i0;
        }
        
        /// <summary>
        /// first indices for each process; size is equal to number of processors (<see cref="IPartitioning.MpiSize"/> plus 1;
        /// first entry is always 0, last entry is equal to <see cref="m_TotalLength"/>;
        /// </summary>
        static public int[] GetI0s(this IPartitioning p) {
            int sz = p.MpiSize;
            int[] R = new int[sz + 1];
            for (int i = 0; i < sz; i++) {
                R[i] = p.GetI0Offest(i);
            }
            R[sz] = p.TotalLength;
            return R;
        }

        /// <summary>
        /// Equality of partitioning.
        /// </summary>
        static public bool EqualsPartition(this IPartitioning t, IPartitioning o) {
            if (t == null && o == null)
                return true;
            if (o == null)
                return false;
            if (t == null)
                return false;
            if (object.ReferenceEquals(t, o))
                return true;

            if (o.MPI_Comm != t.MPI_Comm)
                return false;

            if (o.MpiSize != t.MpiSize)
                return false;
            if (o.MpiRank != t.MpiRank)
                return false;
            if (o.TotalLength != t.TotalLength)
                return false;
            if (o.i0 != t.i0)
                return false;
            if (o.LocalLength != t.LocalLength)
                return false;
            //if (o.IsMutuable != t.IsMutuable)
            //    return false;

            return true;
        }
    }
}
