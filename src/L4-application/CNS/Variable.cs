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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using CNS.Exception;
using System;
using System.Collections;
using System.Collections.Generic;

namespace CNS {

    /// <summary>
    /// Some variables that are recognized in context of the flow state of the
    /// compressible Navier-Stokes equations.
    /// </summary>
    [Flags]
    public enum VariableTypes {

        /// <summary>
        /// The empty state
        /// </summary>
        None = 0x0,

        /// <summary>
        /// The volume density
        /// </summary>
        Density = 0x1,

        /// <summary>
        /// The momentum field
        /// </summary>
        Momentum = 0x2,

        /// <summary>
        /// The velocity field
        /// </summary>
        Velocity = 0x4,

        /// <summary>
        /// The total energy per volume
        /// </summary>
        Energy = 0x8,
        
        /// <summary>
        /// The standard pressure
        /// </summary>
        Pressure = 0x10,

        /// <summary>
        /// A custom field
        /// </summary>
        Other = 0x20,

        /// <summary>
        /// Combination of density, momentum and energy
        /// </summary>
        ConservativeVariables = Density | Momentum | Energy,

        /// <summary>
        /// Combination of density, velocity and pressure
        /// </summary>
        PrimitiveVariables = Density | Velocity | Pressure,
    }

    /// <summary>
    /// Class representing primal variables known to CNS.
    /// </summary>
    [Serializable]
    public class Variable {

        /// <summary>
        /// The unique name of the variable that is used for importing and
        /// exporting data
        /// </summary>
        public readonly string Name;

        /// <summary>
        /// The type of the variable (important for vector-valued variables)
        /// </summary>
        public readonly VariableTypes Type;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name"></param>
        /// <param name="type"></param>
        public Variable(string name, VariableTypes type) {
            this.Name = name;
            this.Type = type;
        }

        /// <summary>
        /// Returns <see cref="Name"/>
        /// </summary>
        /// <param name="variable"></param>
        /// <returns></returns>
        public static implicit operator string(Variable variable) {
            return variable.Name;
        }

        /// <summary>
        /// Returns <see cref="Type"/>
        /// </summary>
        /// <param name="variable"></param>
        /// <returns></returns>
        public static implicit operator VariableTypes(Variable variable) {
            return variable.Type;
        }

        /// <summary>
        /// Returns <see cref="Name"/>
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return Name;
        }

        /// <summary>
        /// Base class for all vector-valued variables, such as the momentum
        /// </summary>
        [Serializable]
        public class Vector<T> : IEnumerable<T> where T : Variable {

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="variableFunc"></param>
            public Vector(Func<int, T> variableFunc) {
                this.xComponent = variableFunc(0);
                this.yComponent = variableFunc(1);
                this.zComponent = variableFunc(2);
            }

            /// <summary>
            /// x-component
            /// </summary>
            public T xComponent {
                get;
                private set;
            }

            /// <summary>
            /// y-component (if D>1)
            /// </summary>
            public T yComponent {
                get;
                private set;
            }

            /// <summary>
            /// z-component (if D=3)
            /// </summary>
            public T zComponent {
                get;
                private set;
            }

            /// <summary>
            /// Access to <see cref="xComponent"/>, <see cref="yComponent"/>
            /// and <see cref="zComponent"/> via indices.
            /// </summary>
            /// <param name="dimension">
            /// The spatial component of the vector
            /// </param>
            /// <returns>
            /// The correct component
            /// </returns>
            public T this[int dimension] {
                get {
                    switch (dimension) {
                        case 0:
                            return xComponent;

                        case 1:
                            return yComponent;

                        case 2:
                            return zComponent;

                        default:
                            throw new InternalErrorException();
                    }
                }
            }

            #region IEnumerable<Variable> Members

            /// <summary>
            /// Enumerator over <see cref="xComponent"/>,
            /// <see cref="yComponent"/> and <see cref="zComponent"/>
            /// </summary>
            /// <returns></returns>
            public IEnumerator<T> GetEnumerator() {
                yield return xComponent;
                yield return yComponent;
                yield return zComponent;
            }

            #endregion

            #region IEnumerable Members

            /// <summary>
            /// Enumerator over <see cref="xComponent"/>,
            /// <see cref="yComponent"/> and <see cref="zComponent"/>
            /// </summary>
            /// <returns></returns>
            IEnumerator IEnumerable.GetEnumerator() {
                return GetEnumerator();
            }

            #endregion
        }
    }

    /// <summary>
    /// Represents auxiliary variables that are updated via some user-defined
    /// update function
    /// </summary>
    public class DerivedVariable : Variable {

        /// <summary>
        /// The update function to be invoked after each (sub-)time-step
        /// </summary>
        public Action<DGField, CellMask, IProgram<CNSControl>> UpdateFunction;

        /// <summary>
        /// See <see cref="Variable"/> and <see cref="UpdateFunction"/>
        /// </summary>
        /// <param name="name"></param>
        /// <param name="type"></param>
        /// <param name="updateFunction"></param>
        public DerivedVariable(string name, VariableTypes type, Action<DGField, CellMask, IProgram<CNSControl>> updateFunction)
            : base(name, type) {
            this.UpdateFunction = updateFunction;
        }
    }
}
