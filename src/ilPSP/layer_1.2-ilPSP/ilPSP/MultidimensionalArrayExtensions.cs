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

using ilPSP.Utils;
using System;

namespace ilPSP {

    /// <summary>
    /// a few extensions for <see cref="MultidimensionalArray"/>
    /// </summary>
    public static class MultidimensionalArrayExtensions {

        /// <summary>
        /// Invokes a transform function on each element of a sequence and
        /// returns the maximum value
        /// </summary>
        static public double Max(this MultidimensionalArray mda, Func<double, double> selector) {
            if (mda.Length <= 0)
                throw new InvalidOperationException();

            double ret = double.MinValue;

            mda.ApplyAll(delegate(double entry) {
                ret = Math.Max(ret, selector(entry));
            });

            return ret;
        }

        /// <summary>
        /// Invokes a transform function on each element of a sequence and
        /// returns the minimum value
        /// </summary>
        static public double Min(this MultidimensionalArray mda, Func<double, double> selector) {
            if (mda.Length <= 0)
                throw new InvalidOperationException();

            double ret = double.MaxValue;

            mda.ApplyAll(delegate(double entry) {
                ret = Math.Min(ret, selector(entry));
            });

            return ret;
        }

        /// <summary>
        /// sum over all entries
        /// </summary>
        static public double Sum(this MultidimensionalArray mda) {
            double ret = 0;

            mda.ApplyAll(delegate(double entry) {
                ret += entry;
            });

            return ret;
        }


        /// <summary>
        /// sum over the absolute value of all entries (i.e. the vector \f$ l_1\f$ --norm);
        /// </summary>
        static public double AbsSum(this MultidimensionalArray mda) {
            double ret = 0;

            mda.ApplyAll(delegate(double entry) {
                ret += Math.Abs(entry);
            });

            return ret;
        }

        /// <summary>
        /// L2-norm over all entries
        /// </summary>
        static public double L2Norm(this MultidimensionalArray mda) {
            double ret = 0;

            mda.ApplyAll(delegate(double entry) {
                ret += entry * entry;
            });

            return Math.Sqrt(ret);
        }

        /// <summary>
        /// L2-norm over all entries
        /// </summary>
        static public double L2Dist(this MultidimensionalArray mda, MultidimensionalArray mdb) {
            if (!ArrayTools.Equals(mda.Lengths, mdb.Lengths, (La, Lb) => La == Lb))
                throw new ArgumentException("Arrays must have the same length.");

            double ret = 0;

            mda.ApplyAll(delegate (int[] idx, double entry_a) {
                double entry_b = mdb[idx];
                double dist = entry_a - entry_b;
                ret += dist * dist;
            });

            return Math.Sqrt(ret);
        }


        /// <summary>
        /// minimum over all entries
        /// </summary>
        static public double Min(this MultidimensionalArray mda) {
            double ret = double.MaxValue;

            mda.ApplyAll(delegate(double entry) {
                ret = Math.Min(ret, entry);
            });

            return ret;
        }

        /// <summary>
        /// maximum over all entries
        /// </summary>
        static public double Max(this MultidimensionalArray mda) {
            double ret = -double.MaxValue;

            mda.ApplyAll(delegate(double entry) {
                ret = Math.Max(ret, entry);
            });

            return ret;
        }
    }
}
