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
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using ilPSP.Utils;
using Newtonsoft.Json;
using Newtonsoft.Json.Converters;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Represents a time-step number. The classical time-step index/number is
    /// represented by <see cref="MajorNumber"/>. In order to handle inner
    /// iterations (and possibly inner iterations of inner iterations, and so
    /// on), an array of numbers can be used. For example, time-step 12 with
    /// solver iteration 2323 and sub-iteration 369 would be time-step number
    /// [12, 2323, 369].
    /// </summary>
    [Serializable]
    [DataContract]
    [JsonObject]
    public struct TimestepNumber : IEnumerable<int>, IEquatable<TimestepNumber>, IComparable<TimestepNumber>, IComparable<int> {

        /// <summary>
        /// The time-step numbers
        /// </summary>
        [DataMember]
        private int[] numbers;

        /// <summary>
        /// Constructs an immutable time-step number
        /// </summary>
        /// <param name="numbers">
        /// An array containing the represented time-step number on each nesting
        /// level, i.e. <paramref name="numbers"/>[0] is the major time-step
        /// number, <paramref name="numbers"/>[1] is one level below, and so on
        /// </param>
        public TimestepNumber(params int[] numbers) {
            this.numbers = numbers;
        }

        /// <summary>
        /// Constructs an immutable time-step number from a string consisting
        /// of dot-separated integers.
        /// </summary>
        /// <param name="numberString">
        /// A string consisting solely of integers separated by dots.
        /// </param>
        public TimestepNumber(string numberString) {
            string[] chunks = (numberString ?? "0").Split('.');
            this.numbers = new int[chunks.Length];
            for (int i = 0; i < chunks.Length; i++) {
                int number;
                if (!int.TryParse(chunks[i], out number)) {
                    throw new ArgumentException(
                        "Given string is not a valid TimestepNumber (integers separated by '.')");
                }

                this.numbers[i] = number;
            }
        }

        /// <summary>
        /// Converts the given <paramref name="numberString"/> to its
        /// <see cref="TimestepNumber"/> equivalent.
        /// </summary>
        /// <param name="numberString">
        /// The string, which should be converted.
        /// </param>
        /// <param name="Timestep">
        /// The return value, if the conversion was successful.
        /// </param>
        /// <returns>
        /// Indicates, whether the conversion was successful.
        /// </returns>
        public static bool TryParse(string numberString, out TimestepNumber Timestep) {
            bool IsTimestepNumber = true;

            try {
                Timestep = new TimestepNumber(numberString);
            } catch (ArgumentException) {
                Timestep = null;
                IsTimestepNumber = false;
            }

            return IsTimestepNumber;
        }

        /// <summary>
        /// The major/classical time-step number
        /// </summary>
        public int MajorNumber {
            get {
                return numbers[0];
            }
        }

        /// <summary>
        /// Allows to access the i-th component of the time-step number. The
        /// significance of the number is sorted in decreasing order, i.e. the
        /// significance of this[0] is larger than of this[1] and so on.
        /// </summary>
        /// <param name="index">
        /// The index into this number.
        /// </param>
        /// <returns>
        /// The i-th component of the time-step number
        /// </returns>
        public int this[int index] {
            get {
                if (numbers.Length < index) {
                    return 0;
                } else {
                    return numbers[index];
                }
            }
        }

        /// <summary>
        /// The number of levels of this time-step number.
        /// </summary>
        public int Length {
            get {
                return numbers.Length;
            }
        }

        /// <summary>
        /// Returns the next number if level <paramref name="level"/> is
        /// advanced
        /// </summary>
        /// <param name="level">
        /// The time-step level to be advanced (e.g., 0 for time-step, 1 for
        /// an iteration within a time-step, and 2 for an inner iteration
        /// within an outer iteration within a time-step)
        /// </param>
        /// <returns>
        /// A new time-step number where the <paramref name="level"/>-th number
        /// in <see cref="P:Item(System.Int32)"/> is incremented by one.
        /// </returns>
        public TimestepNumber NextNumber(uint level) {
            int[] newNumbers = new int[Math.Max(numbers.Length, level + 1)];

            // Only copy part up $level
            Array.Copy(numbers, newNumbers, Math.Min(numbers.Length, level + 1));
            newNumbers[level]++;
            return new TimestepNumber(newNumbers);
        }

        /// <summary>
        /// Returns the number of the next time-step.
        /// </summary>
        /// <returns></returns>
        public TimestepNumber NextTimestep() {
            return NextNumber(0);
        }

        /// <summary>
        /// Returns the number of the next iteration within the given time-step
        /// </summary>
        /// <returns></returns>
        public TimestepNumber NextIteration() {
            return NextNumber(1);
        }

        /// <summary>
        /// Returns the number of the next sub-iteration within the current
        /// (outer) iteration of the current time-step
        /// </summary>
        /// <returns></returns>
        public TimestepNumber NextSubIteration() {
            return NextNumber(2);
        }

        /// <summary>
        /// Appends the numbers in the order of their significance.
        /// </summary>
        /// <returns>
        /// Given a number [1, 2, 3], returns '1.2.3'
        /// </returns>
        public override string ToString() {
            return numbers.Skip(1).Aggregate(MajorNumber.ToString(),
                (s, t) => s + "." + t);
        }

        #region IEnumerable<int> Members

        /// <summary>
        /// See <see cref="IEnumerable{Int32}.GetEnumerator"/>
        /// </summary>
        /// <returns></returns>
        public IEnumerator<int> GetEnumerator() {
            return numbers.AsEnumerable().GetEnumerator();
        }

        #endregion

        #region IEnumerable Members

        /// <summary>
        /// See <see cref="T:int[].GetEnumerator"/>
        /// </summary>
        /// <returns></returns>
        IEnumerator IEnumerable.GetEnumerator() {
            return numbers.GetEnumerator();
        }

        #endregion

        #region IEquatable<TimeStepIndex> Members

        /// <summary>
        /// See <see cref="Equals(TimestepNumber)"/>
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            return obj is TimestepNumber ? Equals((TimestepNumber)obj) : false;
        }

        /// <summary>
        /// See <see cref="T:int[].GetHashCode"/>
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode() {
            return numbers.GetHashCode();
        }

        /// <summary>
        /// See <see cref="T:int[].Equals"/>
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(TimestepNumber other) {
            return ArrayTools.Equals(numbers, other.numbers);
        }

        #endregion

        #region IComparable<TimeStepIndex> Members

        /// <summary>
        /// Compares two time-steps based on each individual number (cf.
        /// <see cref="TimestepNumber"/>) in order of their
        /// significance.
        /// </summary>
        /// <param name="other">
        /// The time-step to be compared to
        /// </param>
        /// <returns>
        /// The result of the integer comparison of the sub-number with the
        /// highest significance for which the result is not zero.
        /// </returns>
        public int CompareTo(TimestepNumber other) {
            int length = Math.Min(this.numbers.Length, other.numbers.Length);

            // Decision based on existing entries
            for (int i = 0; i < length; i++) {
                int result = this.numbers[i].CompareTo(other.numbers[i]);
                if (result != 0) {
                    return result;
                }
            }

            // Still undecided? -> More indices = larger
            return this.numbers.Length.CompareTo(other.numbers.Length);
        }

        #endregion

        #region IComparable<int> Members

        /// <summary>
        /// Compares this time-step number to <paramref name="other"/> by
        /// converting <paramref name="other"/> to a
        /// <see cref="TimestepNumber"/> and using
        /// <see cref="CompareTo(TimestepNumber)"/>
        /// </summary>
        /// <param name="other">
        /// The time-step number to be compared
        /// </param>
        /// <returns>
        /// See <see cref="CompareTo(TimestepNumber)"/>
        /// </returns>
        public int CompareTo(int other) {
            return this.CompareTo(new TimestepNumber(other));
        }

        #endregion

        /// <summary>
        /// Implicit conversion for time-steps with just a single index.
        /// </summary>
        /// <param name="number">
        /// The number of the time-step to be represented
        /// </param>
        /// <returns>
        /// The number <paramref name="number"/> as
        /// <see cref="TimestepNumber"/>
        /// </returns>
        public static implicit operator TimestepNumber(int number) {
            return new TimestepNumber(number);
        }

        /// <summary>
        /// Implicit conversion time-steps given as strings (integers separated
        /// by dots).
        /// </summary>
        /// <param name="numberAsString">
        /// A string following the pattern "1", "1.2", "1.2.3" and the like.
        /// </param>
        /// <returns>
        /// The parsed time-step number.
        /// </returns>
        /// <remarks>
        /// In principle, this shouldn't be an implicit conversion, since it
        /// can fail. However, it greatly simplifies the usage of the command
        /// line which should justify this code smell.
        /// </remarks>
        public static implicit operator TimestepNumber(string numberAsString) {
            return new TimestepNumber(numberAsString);
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for equality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator ==(TimestepNumber x, int y) {
            return x.CompareTo(y) == 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for inequality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator !=(TimestepNumber x, int y) {
            return x.CompareTo(y) != 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for equality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator ==(int x, TimestepNumber y) {
            return y == x;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for inequality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator !=(int x, TimestepNumber y) {
            return y == x;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator <(TimestepNumber x, int y) {
            return x.CompareTo(y) < 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator >(TimestepNumber x, int y) {
            return x.CompareTo(y) > 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator <(int x, TimestepNumber y) {
            return y > x;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator >(int x, TimestepNumber y) {
            return y < x;
        }






        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for equality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator ==(TimestepNumber x, string y) {
            return x.CompareTo(y) == 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for inequality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator !=(TimestepNumber x, string y) {
            return x.CompareTo(y) != 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for equality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator ==(string x, TimestepNumber y) {
            return y == x;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to check <paramref name="x"/>
        /// and <paramref name="y"/> for inequality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator !=(string x, TimestepNumber y) {
            return y == x;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator <(TimestepNumber x, string y) {
            return x.CompareTo(y) < 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator >(TimestepNumber x, string y) {
            return x.CompareTo(y) > 0;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator <(string x, TimestepNumber y) {
            return y > x;
        }

        /// <summary>
        /// Uses <see cref="CompareTo(int)"/> to compare <paramref name="x"/>
        /// and <paramref name="y"/>.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator >(string x, TimestepNumber y) {
            return y < x;
        }
    }
}
