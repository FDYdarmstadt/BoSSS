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
using System.Linq;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// The first part of a full <see cref="Guid"/>. Useful for the selection
    /// of objects without having to specify the full-length
    /// <see cref="Guid"/>.
    /// </summary>
    public struct PartialGuid : IComparable, IComparable<PartialGuid>, IEquatable<PartialGuid> {

        /// <summary>
        /// The length of a guid <b>by definition</b>
        /// </summary>
        private const int GUID_LENGTH = 32;

        /// <summary>
        /// Use a full guid as internal storage so that validation routines are
        /// automatically reused
        /// </summary>
        private Guid guid;

        /// <summary>
        /// The number of digits of <see cref="guid"/> that are significant for
        /// this guid (i.e., the length of this partial guid)
        /// </summary>
        private int significantLength;

        /// <summary>
        /// Initializes a new instance of the <see cref="PartialGuid"/> struct.
        /// </summary>
        /// <param name="partialGuid">
        /// A portion of a guid, i.e. a substring of a string in a format that
        /// defines a valid <see cref="Guid"/>.
        /// </param>
        public PartialGuid(string partialGuid) {
            partialGuid = partialGuid.Replace("-", "");
            if (partialGuid.Length > GUID_LENGTH) {
                throw new ArgumentException(String.Format(
                    "A guid must not be longer than {0} characters",
                    GUID_LENGTH));
            }

            significantLength = partialGuid.Length;
            guid = Guid.Parse(partialGuid.PadRight(GUID_LENGTH, '0'));
        }

        #region IComparable Members

        /// <summary>
        /// See <see cref="CompareTo(PartialGuid)"/>
        /// </summary>
        /// <param name="obj">
        /// See <see cref="CompareTo(PartialGuid)"/>
        /// </param>
        /// <returns>
        /// See <see cref="CompareTo(PartialGuid)"/>
        /// </returns>
        public int CompareTo(object obj) {
            if (obj is PartialGuid) {
                return CompareTo((PartialGuid)obj);
            } else if (obj is Guid) {
                return CompareTo((Guid)obj);
            } else {
                throw new ArgumentException(
                    "Cannot compare to non-Guids");
            }
        }

        #endregion

        #region IComparable<PartialGuid> Members

        /// <summary>
        /// Compares this guid to the given <paramref name="other"/> by
        /// comparing the respective string representations.
        /// </summary>
        /// <param name="other">
        /// The partial guid to be compared to.
        /// </param>
        /// <returns>
        /// See <see cref="IComparable{T}.CompareTo"/>.
        /// </returns>
        public int CompareTo(PartialGuid other) {
            return this.ToString().CompareTo(other.ToString());
        }

        #endregion

        #region IEquatable<PartialGuid> Members

        /// <summary>
        /// Checks of this guid is equal to the given
        /// <paramref name="guid"/> by comparing only the respective string
        /// representation <b>up to the length of the shorter guid</b>. As a
        /// result, this is a diffuse equality comparison.
        /// </summary>
        /// <param name="guid">
        /// The guid to be compared to.
        /// </param>
        /// <returns>
        /// True, if both strings are <b>approximately</b> equal.
        /// </returns>
        public bool Equals(PartialGuid guid) {
            if (significantLength <= guid.significantLength) {
                return this.ToString().Equals(
                    guid.ToString().Substring(0, this.significantLength));
            } else {
                return guid.ToString().Equals(
                    ToString().Substring(0, guid.significantLength));
            }
        }

        #endregion

        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this
        /// instance.
        /// </summary>
        /// <returns>
        /// A <see cref="System.String" /> that represents this instance.
        /// </returns>
        public override string ToString() {
            return guid.ToString("N", null).Substring(0, significantLength);
        }

        /// <summary>
        /// See <see cref="Equals(PartialGuid)"/>
        /// </summary>
        /// <param name="obj">
        /// See <see cref="Equals(PartialGuid)"/>
        /// </param>
        /// <returns>
        /// See <see cref="Equals(PartialGuid)"/>
        /// </returns>
        public override bool Equals(object obj) {
            if (obj is PartialGuid) {
                return Equals((PartialGuid)obj);
            } else if (obj is Guid) {
                return Equals((Guid)obj);
            } else {
                return false;
            }
        }

        /// <summary>
        /// Returns a hash code for this instance.
        /// </summary>
        /// <returns>
        /// A hash code for this instance, suitable for use in hashing
        /// algorithms and data structures like a hash table. 
        /// </returns>
        public override int GetHashCode() {
            return guid.ToByteArray().Take(significantLength).GetHashCode();
        }

        /// <summary>
        /// See <see cref="Equals(PartialGuid)"/>
        /// </summary>
        /// <param name="a">
        /// See <see cref="Equals(PartialGuid)"/>
        /// </param>
        /// <param name="b">
        /// See <see cref="Equals(PartialGuid)"/>
        /// </param>
        /// <returns>
        /// See <see cref="Equals(PartialGuid)"/>
        /// </returns>
        public static bool operator ==(PartialGuid a, PartialGuid b) {
            return a.Equals(b);
        }

        /// <summary>
        /// Negation of <see cref="Equals(PartialGuid)"/>
        /// </summary>
        /// <param name="a">
        /// See <see cref="Equals(PartialGuid)"/>
        /// </param>
        /// <param name="b">
        /// See <see cref="Equals(PartialGuid)"/>
        /// </param>
        /// <returns>
        /// Negation of <see cref="Equals(PartialGuid)"/>
        /// </returns>
        public static bool operator !=(PartialGuid a, PartialGuid b) {
            return !a.Equals(b);
        }

        /// <summary>
        /// Implicit conversion from a string to a partial guid
        /// </summary>
        /// <param name="partialGuid">
        /// A string representing a portion of a <see cref="Guid"/>.
        /// </param>
        /// <returns>
        /// A partial guid.
        /// </returns>
        /// <remarks>
        /// In some sense, this operator violates the contract for an implicit
        /// conversion since it may fail if <paramref name="partialGuid"/> is
        /// not in a suitable format. However, the additional benefit when
        /// using a console interface outweigh this disadvantage.
        /// </remarks>
        public static implicit operator PartialGuid(string partialGuid) {
            return new PartialGuid(partialGuid);
        }

        /// <summary>
        /// Implicit conversion from a <see cref="Guid"/> to a partial guid.
        /// </summary>
        /// <param name="guid">
        /// A string representing a portion of a <see cref="Guid"/>.
        /// </param>
        /// <returns>
        /// A partial guid.
        /// </returns>
        public static implicit operator PartialGuid(Guid guid) {
            return new PartialGuid(guid.ToString("N", null));
        }
    }
}
