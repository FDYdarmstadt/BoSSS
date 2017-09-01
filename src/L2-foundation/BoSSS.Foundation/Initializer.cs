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

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Represents a type that is initializable during IO operations.
    /// </summary>
    /// <typeparam name="T">
    /// The type of the initialized object.
    /// </typeparam>
    public interface IInitializer<out T> {

        /// <summary>
        /// Uses the given initialization context <paramref name="c"/> in order
        /// create an object of the given type T.
        /// </summary>
        /// <param name="c">
        /// Context containing already initialized objects and information
        /// about the grid
        /// </param>
        /// <returns>
        /// A properly initialized object instance.
        /// </returns>
        T Initialize(IInitializationContext c);
    }
    
    /// <summary>
    /// Base class of all initializers that provides a default behavior for
    /// comparison operations.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    [Serializable]
    public abstract class Initializer<T> : IInitializer<T>, IEquatable<Initializer<T>> {

        /// <summary>
        /// See <see cref="IInitializer{T}.Initialize"/>
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public abstract T Initialize(IInitializationContext c);

        /// <summary>
        /// Uses <see cref="Equals(Initializer{T})"/> to check for equality
        /// </summary>
        /// <param name="obj">
        /// The object to be compared to
        /// </param>
        /// <returns>
        /// True, if <paramref name="obj"/> refers to the same instance as this
        /// object.
        /// </returns>
        public override bool Equals(object obj) {
            if (obj is Initializer<T> == false) {
                return false;
            }

            return this.Equals((Initializer<T>)obj);
        }

        /// <summary>
        /// See <see cref="object.GetHashCode"/>
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        #region IEquatable<Initializer<T>> Members

        /// <summary>
        /// Checks for reference equality. May be overridden in sub-classes to
        /// provide a different notion of equality
        /// </summary>
        /// <param name="other">
        /// The object to be compared to
        /// </param>
        /// <returns>
        /// True, if <paramref name="other"/> refers to the same instance as
        /// this object.
        /// </returns>
        public virtual bool Equals(Initializer<T> other) {
            return object.ReferenceEquals(other, this);
        }

        #endregion
    }
}
