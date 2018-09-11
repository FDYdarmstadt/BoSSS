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

using BoSSS.Foundation.Grid;
using System.Collections.Generic;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// A context for the initialization of objects for IO purposes. Its main
    /// purpose is the caching of objects such that they're not recreated or
    /// deserialized hundreds of times.
    /// </summary>
    public interface IInitializationContext {

        /// <summary>
        /// Information about the grid.
        /// </summary>
        IGridData GridData {
            get;
        }

        /// <summary>
        /// Adds object <paramref name="value"/> to the cache
        /// </summary>
        /// <typeparam name="T">
        /// The type of the object to be cached. Must match the generic type
        /// argument of the initializer
        /// </typeparam>
        /// <param name="initializer">
        /// The initializer that has been used to create
        /// <paramref name="value"/>.
        /// </param>
        /// <param name="value">
        /// The object to be cached
        /// </param>
        void Add<T>(IInitializer<object> initializer, T value);

        /// <summary>
        /// Tries to retrieve a value for the given
        /// <paramref name="initializer"/>.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the object to be cached. Must match the generic type
        /// argument of the initializer
        /// </typeparam>
        /// <param name="initializer">
        /// The initializer that has been used to create
        /// <paramref name="value"/>.
        /// </param>
        /// <param name="value">
        /// On exit: If the object has already been created by the given
        /// <paramref name="initializer"/>, the corresponding cached object.
        /// Otherwise, null is returned.
        /// </param>
        /// <returns>
        /// True, if the object could be retrieved from the cache; false
        /// otherwise.
        /// </returns>
        bool TryGetValue<T>(IInitializer<object> initializer, out T value);
    }

    /// <summary>
    /// A context for the initialization of objects for IO purposes. Its main
    /// purpose is the caching of objects such that they're not recreated or
    /// deserialized hundreds of times. This particular implementation has the
    /// scope of a grid object, i.e. it can be considered as a global cache for
    /// all objects living on the same grid.
    /// </summary>
    public sealed class GridInitializationContext : IInitializationContext {

        /// <summary>
        /// A directory of objects that have already been deserialized within
        /// this context.
        /// </summary>
        /// <remarks>
        /// Now uses the default equality comparer instead of reference
        /// comparisons since it is almost impossible to ensure reference
        /// equality for the initializers after deserialization (i.e., this
        /// whole construct does not have a point anymore since reference
        /// equality never holds).
        /// </remarks>
        private readonly Dictionary<IInitializer<object>, object> objectDirectory =
            new Dictionary<IInitializer<object>, object>();

        /// <summary>
        /// Constructs a new context.
        /// </summary>
        public GridInitializationContext(IGridData gridData) {
            this.GridData = gridData;
        }

        #region IInitializationContext Members

        /// <summary>
        /// Information about the grid.
        /// </summary>
        public IGridData GridData {
            get;
            private set;
        }

        /// <summary>
        /// See <see cref="IInitializationContext"/>
        /// </summary>
        public void Add<T>(IInitializer<object> initializer, T value) {
            objectDirectory.Add(initializer, value);
        }

        /// <summary>
        /// See <see cref="IInitializationContext"/>
        /// </summary>
        public bool TryGetValue<T>(IInitializer<object> initializer, out T value) {
            object val;
            bool success = objectDirectory.TryGetValue(initializer, out val);
            value = (T)val;
            return success;
        }

        #endregion
    }

    /// <summary>
    /// A context for the initialization of objects for IO purposes. Its main
    /// purpose is the caching of objects such that they're not recreated or
    /// deserialized hundreds of times. This particular implementation has the
    /// scope of a grid object, i.e. it can be considered as a global cache for
    /// all objects living on the same grid.
    /// </summary>
    public sealed class TimestepInitializationContext : IInitializationContext {

        /// <summary>
        /// Grid-global context for things like basis objects which don't
        /// depend on the actual time step.
        /// </summary>
        private GridInitializationContext parentContext;

        /// <summary>
        /// A directory of objects that have already been deserialized within
        /// this context.
        /// </summary>
        /// <remarks>
        /// Now uses the default equality comparer instead of reference
        /// comparisons since it is almost impossible to ensure reference
        /// equality for the initializers after deserialization (i.e., this
        /// whole construct does not have a point anymore since reference
        /// equality never holds).
        /// </remarks>
        private readonly Dictionary<IInitializer<object>, object> objectDirectory =
            new Dictionary<IInitializer<object>, object>();

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="parentContext">
        /// The context of the grid where this time-step context lives.
        /// </param>
        public TimestepInitializationContext(GridInitializationContext parentContext) {
            this.parentContext = parentContext;
        }

        #region IInitializationContext Members

        /// <summary>
        /// Information about the grid.
        /// </summary>
        public IGridData GridData {
            get {
                return parentContext.GridData;
            }
        }

        /// <summary>
        /// Adds the given <paramref name="value"/> to the list of cached
        /// objects. Depending <typeparamref name="T"/>, the object will either
        /// go to the local scope or to the grid-global scope (cf.
        /// <see cref="TimestepInitializationContext.TimestepInitializationContext"/>)
        /// </summary>
        /// <typeparam name="T">
        /// The of the <paramref name="value"/>
        /// </typeparam>
        /// <param name="initializer">
        /// The initializer that has been used to created
        /// <paramref name="value"/>
        /// </param>
        /// <param name="value">
        /// The value to be cached.
        /// </param>
        public void Add<T>(IInitializer<object> initializer, T value) {
            if (typeof(T) == typeof(Basis)) {
                parentContext.Add(initializer, value);
            } else {
                // Deactivated for now to avoid memory overflows. 
                objectDirectory.Add(initializer, value);
            }
        }

        /// <summary>
        /// Tries to retrieve a value for the given
        /// <paramref name="initializer"/>. If the value is not in the local
        /// scope of this object, it will be tried to retrieve it from the
        /// grid-global scope (cf.
        /// <see cref="TimestepInitializationContext.TimestepInitializationContext"/>)
        /// </summary>
        /// <typeparam name="T">
        /// The type of the object to be cached. Must match the generic type
        /// argument of the initializer
        /// </typeparam>
        /// <param name="initializer">
        /// The initializer that has been used to create
        /// <paramref name="value"/>.
        /// </param>
        /// <param name="value">
        /// On exit: If the object has already been created by the given
        /// <paramref name="initializer"/>, the corresponding cached object.
        /// Otherwise, null is returned.
        /// </param>
        /// <returns>
        /// True, if the object could be retrieved from the cache; false
        /// otherwise.
        /// </returns>
        public bool TryGetValue<T>(IInitializer<object> initializer, out T value) {
            object val;
            bool success = objectDirectory.TryGetValue(initializer, out val);
            if (!success) {
                success = parentContext.TryGetValue(initializer, out val);
            }
            
            value = (T)val;
            return success;
        }

        #endregion
    }
}
