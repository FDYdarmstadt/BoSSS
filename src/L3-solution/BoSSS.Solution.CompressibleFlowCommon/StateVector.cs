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

using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Solution.CompressibleFlowCommon {

    /// <summary>
    /// Struct representing the flow state in a single point which is defined by
    /// the <see cref="Density"/>, the <see cref="Momentum"/> and the
    /// <see cref="Energy"/>. The state represented by an instance of this
    /// class is immutable.
    /// </summary>
    public struct StateVector {

        /// <summary>
        /// Common basis for all constructors
        /// </summary>
        /// <param name="material">
        /// The equation of state associated with the represented state.
        /// </param>
        private StateVector(Material material)
            : this() {
            this.Material = material;
        }

        /// <summary>
        /// Constructs a state with the given data.
        /// </summary>
        /// <param name="material">
        /// The material/fluid associated with the represented state.
        /// </param>
        /// <param name="density"><see cref="Density"/></param>
        /// <param name="momentum"><see cref="Momentum"/></param>
        /// <param name="energy"><see cref="Energy"/></param>
        public StateVector(Material material, double density, Vector3D momentum, double energy)
            : this(material) {
            this.Density = density;
            this.Momentum = momentum;
            this.Energy = energy;
        }

        /// <summary>
        /// Takes an array representation <paramref name="stateVectorAsArray"/>
        /// and saves the values in <see cref="Density"/> and
        /// <see cref="Momentum"/> and <see cref="Energy"/>.
        /// </summary>
        /// <param name="stateVectorAsArray">
        /// The state to be represented
        /// </param>
        /// <param name="material">
        /// The material/fluid associated with the represented state.
        /// </param>
        public StateVector(IList<double> stateVectorAsArray, Material material)
            : this(material) {
            if (stateVectorAsArray.Count < CNSEnvironment.NumberOfDimensions + 2) {
                throw new ArgumentException(
                    "The given state vector has an invalid length. In n dimensions,"
                        + " the length should at least be n + 2",
                    "stateVectorAsArray");
            }

            this.Density = stateVectorAsArray[0];
            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                this.Momentum[d] = stateVectorAsArray[d + 1];
            }
            this.Energy = stateVectorAsArray[CNSEnvironment.NumberOfDimensions + 1];
        }

        /// <summary>
        /// Constructs a new state vector using data given in primitive
        /// variables
        /// </summary>
        /// <param name="material">
        /// The material/fluid associated with the represented state.
        /// </param>
        /// <param name="density"><see cref="Density"/></param>
        /// <param name="velocity"><see cref="Velocity"/></param>
        /// <param name="pressure"><see cref="Pressure"/></param>
        /// <returns></returns>
        public static StateVector FromPrimitiveQuantities(Material material, double density, Vector3D velocity, double pressure) {
            double MachScaling = material.EquationOfState.HeatCapacityRatio * material.Control.MachNumber * material.Control.MachNumber;
            StateVector state = new StateVector(
                material,
                density,
                density * velocity,
                material.EquationOfState.GetInnerEnergy(density, pressure) + 0.5 * MachScaling * density * velocity.AbsSquare());
            return state;
        }

        /// <summary>
        /// Simplifies the construction of a state vector in fluxes where the
        /// conserved quantities are typically given as an array of
        /// two-dimensional arrays.
        /// </summary>
        /// <param name="material">
        /// The material/fluid associated with the represented state.
        /// </param>
        /// <param name="stateAsArray">
        /// An array of two-dimensional instances of
        /// <see cref="MultidimensionalArray"/> containing the n-th conserved
        /// variable value at U[n][<paramref name="i"/>, <paramref name="j"/>]
        /// </param>
        /// <param name="i">
        /// First index into the entries of <paramref name="stateAsArray"/>
        /// </param>
        /// <param name="j">
        /// Second index into the entries of <paramref name="stateAsArray"/>
        /// </param>
        public StateVector(Material material, MultidimensionalArray[] stateAsArray, int i, int j)
            : this(material) {
            int D = CNSEnvironment.NumberOfDimensions;
            if (stateAsArray.Length < D + 2) {
                throw new ArgumentException(
                    "The given state vector has an invalid length. In n dimensions,"
                        + " the length should at least be n + 2",
                    "stateVectorAsArray");
            }

            this.Density = stateAsArray[0][i, j];
            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                this.Momentum[d] = stateAsArray[d + 1][i, j];
            }
            this.Energy = stateAsArray[D + 1][i, j];
        }

        /// <summary>
        /// The mass density
        /// </summary>
        public readonly double Density;

        /// <summary>
        /// The momentum vector $\rho \vec{u}$. If
        /// <see cref="CNSEnvironment.NumberOfDimensions"/> is less than
        /// three, surplus components are set to zero.
        /// </summary>
        public readonly Vector3D Momentum;

        /// <summary>
        /// The total energy per volume $\rho E$
        /// </summary>
        public readonly double Energy;

        /// <summary>
        /// The material/fluid associated with the represented state.
        /// </summary>
        public Material Material {
            get;
            private set;
        }

        /// <summary>
        /// Calculates the velocity from <see cref="Momentum"/> and
        /// <see cref="Density"/>.
        /// </summary>
        /// <returns>$\vec{m} / \rho$</returns>
        public Vector3D Velocity {
            get {
                return Momentum / Density;
            }
        }

        /// <summary>
        /// Calculates the kinetic energy $K$ via
        /// $K = \frac{1}{2} \frac{\vec{m} \cdot \vec{m}}{\rho}$
        /// </summary>
        public double KineticEnergy {
            get {
                double kineticEnergy = 0.5 * Momentum.AbsSquare() / Density;
                // needs scaling according to Nondimensionalization
                double dimensionalScaling = Material.EquationOfState.HeatCapacityRatio * Material.Control.MachNumber * Material.Control.MachNumber;
                return dimensionalScaling * kineticEnergy;
            }
        }

        /// <summary>
        /// Calculates the specific kinetic energy \f$ k \f$ via
        ///  \f$ k = \frac{1}{2} \vec{u} \cdot \vec{u} \f$
        /// </summary>
        public double SpecificKineticEnergy {
            get {
                return KineticEnergy / Density;
            }
        }

        /// <summary>
        /// Calculates the inner energy  \f \rho e  \f via  \f \rho e = (\rho E - K)  \f$
        /// </summary>
        public double InnerEnergy {
            get {
                double innerEnergy = Energy - KineticEnergy;
                return innerEnergy;
            }
        }

        /// <summary>
        /// Calculates the specific inner energy  \f$e \f$ via  \f$e = (\rho E - K) / \rho \f$
        /// </summary>
        public double SpecificInnerEnergy {
            get {
                return InnerEnergy / Density;
            }
        }

        /// <summary>
        /// Calculate the pressure  \f$p \f$ via
        /// <see cref="Material"/>.<see cref="IEquationOfState.GetPressure"/>.
        /// </summary>
        public double Pressure {
            get {
                double pressure = Material.EquationOfState.GetPressure(this);
                return pressure;
            }
        }

        /// <summary>
        /// Calculate the pressure  \f$ T \f$ via
        /// <see cref="Material"/>.<see cref="IEquationOfState.GetTemperature"/>.
        /// </summary>
        public double Temperature {
            get {
                double temperature = Material.EquationOfState.GetTemperature(this);
                return temperature;
            }
        }

        /// <summary>
        /// Calculates the dynamics viscosity  \f$ \nu \f$ using the Sutherland model
        /// which can be written as  \f$ \nu = T^\frac{3}{2} \frac{1 + S}{T + S}  \f$
        /// </summary>
        /// <remarks>
        /// The Sutherland temperature  \f$ S \f$ is currently assuming room
        /// temperature in the far field
        /// </remarks>
        public double GetViscosity(int cellIndex) {
            double viscosity = Material.ViscosityLaw.GetViscosity(Temperature, cellIndex);
            return viscosity;
        }

        /// <summary>
        /// Calculates the local speed of sound  \f$ a  \f$ via
        ///  \f$\mathrm{a} = \sqrt{\kappa \frac{p}{\rho}} \f$
        /// </summary>
        public double SpeedOfSound {
            get {
                //Dimensional scaling is inside GetSpeedOfSound
                double sos = Material.EquationOfState.GetSpeedOfSound(this);
                return sos;
            }
        }

        /// <summary>
        /// Calculates the Mach number Ma via
        /// \f$ \mathrm{Ma} = \frac{\sqrt{\vec{u} \cdot \vec{u}}}{a}\f$ 
        /// where a is the speed of sound (<see cref="SpeedOfSound"/>).
        /// </summary>
        public double LocalMachNumber {
            get {
                return Velocity.Abs() / SpeedOfSound;
            }
        }

        /// <summary>
        /// Calculates the specific enthalpy h via
        /// \f$ h = \frac{\rho E + p}{\rho}\f$ 
        /// </summary>
        public double Enthalpy {
            get {
                return (Energy + Pressure) / Density;
            }
        }

        /// <summary>
        /// Calculates the local entropy using <see cref="Material"/>
        /// (see <see cref="IEquationOfState.GetEntropy"/>).
        /// </summary>
        public double Entropy {
            get {
                double entropy = Material.EquationOfState.GetEntropy(this);
                return entropy;
            }
        }

        /// <summary>
        /// Returns the variables in the fixed order:
        /// rho, m1, [m2, [m3, ]] rhoE
        /// </summary>
        /// <returns></returns>
        public double[] ToArray() {
            int D = CNSEnvironment.NumberOfDimensions;

            Debug.Assert(CNSEnvironment.PrimalArgumentToIndexMap[Variables.Density] == 0);
            for (int d = 0; d < D; d++) {
                Debug.Assert(CNSEnvironment.PrimalArgumentToIndexMap[Variables.Momentum[d]] == d + 1);
            }
            Debug.Assert(CNSEnvironment.PrimalArgumentToIndexMap[Variables.Energy] == D + 1);

            double[] ret = new double[D + 2];
            ret[0] = Density;
            for (int d = 0; d < D; d++) {
                ret[d + 1] = Momentum[d];
            }
            ret[D + 1] = Energy;
            return ret;
        }

        /// <summary>
        /// Returns true if the represented state defines a physically
        /// reasonable configuration.
        /// </summary>
        public bool IsValid {
            get {
                return Density > 0.0
                    && InnerEnergy > 0.0;
            }
        }
    }
}
