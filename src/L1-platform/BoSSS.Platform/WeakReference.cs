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
using System.Runtime.InteropServices;

namespace BoSSS.Platform {

    /// <summary>
    /// Implementation of <see cref="WeakReference"/> using generics.
    /// </summary>
    /// <typeparam name="T">
    /// The type of the object that should be referenced.
    /// </typeparam>
    /// <remarks>
    /// Adapted from
    /// http://blogs.microsoft.co.il/blogs/pavely/archive/2011/06/04/a-proper-weakreference-class.aspx
    /// </remarks>
    public sealed class WeakReference<T> : IDisposable where T : class {

        /// <summary>
        /// Weak reference to the object.
        /// </summary>
        private GCHandle m_Handle;

        /// <summary>
        /// Constructs the weak reference.
        /// </summary>
        /// <param name="referencedObject">
        /// The object to be referenced weakly.
        /// </param>
        public WeakReference(T referencedObject) {
            if (referencedObject == null) {
                throw new ArgumentException("Referenced object must not be null");
            }

            m_Handle = GCHandle.Alloc(referencedObject, GCHandleType.Weak);
        }

        /// <summary>
        /// The referenced object if is still alive. Otherwise, null is
        /// returned.
        /// </summary>
        public T Target {
            get {
                if (m_Handle.IsAllocated) {
                    return m_Handle.Target as T;
                } else {
                    return null;
                }
            }
        }

        /// <summary>
        /// True, if the referenced object has not yet been collected by the
        /// garbage collector. False, otherwise.
        /// </summary>
        public bool IsAlive {
            get {
                return m_Handle.Target != null;
            }
        }

        /// <summary>
        /// Checks whether this object is equal to the given object.
        /// </summary>
        /// <param name="otherReference">
        /// The reference in question.
        /// </param>
        /// <returns>
        /// True, if the given object references the same object as this
        /// object. False, otherwise.
        /// </returns>
        public bool Equals(WeakReference<T> otherReference) {
            if (otherReference == null) {
                return false;
            } else {
                return m_Handle.Equals(otherReference.m_Handle);
            }
        }

        /// <summary>
        /// Checks whether this object is equal to the given object.
        /// </summary>
        /// <param name="obj">
        /// The object in question.
        /// </param>
        /// <returns>
        /// True, if the given object is a weak reference referencing the same
        /// object as this object. False, otherwise.
        /// </returns>
        public override bool Equals(object obj) {
            WeakReference<T> otherReference = obj as WeakReference<T>;
            if (otherReference == null) {
                return false;
            } else {
                return Equals(otherReference);
            }
        }

        /// <summary>
        /// <see cref="Object.GetHashCode"/>
        /// </summary>
        /// <returns>
        /// The hash code of the referenced object.
        /// </returns>
        public override int GetHashCode() {
            return m_Handle.GetHashCode();
        }

        /// <summary>
        /// Frees the reference to the referenced object.
        /// </summary>
        ~WeakReference() {
            m_Handle.Free();
        }

        /// <summary>
        /// <see cref="Equals(WeakReference{T})"/>
        /// </summary>
        public static bool operator ==(WeakReference<T> left, WeakReference<T> right) {
            return left.Equals(right);
        }

        /// <summary>
        /// not <see cref="Equals(WeakReference{T})"/>
        /// </summary>
        /// <returns>
        /// The negation of the result of <see cref="Equals(WeakReference{T})"/>.
        /// </returns>
        public static bool operator !=(WeakReference<T> left, WeakReference<T> right) {
            return !left.Equals(right);
        }

        #region IDisposable Members

        /// <summary>
        /// Frees the reference to the referenced object.
        /// </summary>
        public void Dispose() {
            m_Handle.Free();
            GC.SuppressFinalize(this);
        }

        #endregion
    }
}
