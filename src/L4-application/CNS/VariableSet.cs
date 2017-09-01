using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Platform.LinAlg;

namespace CNS {

    public class VariableSet {

        /// <summary>
        /// { "rho", "p1" [,"p2" [, "p3"]], "E" }
        /// </summary>
        public readonly string[] ArgumentOrdering;

        /// <summary>
        /// A mapping between the variable names and their position in
        /// <see cref="ArgumentOrdering"/>
        /// </summary>
        public readonly Dictionary<string, int> ArgumentToIndexMap;

        /// <summary>
        /// Configuration options
        /// </summary>
        private CNSControl.CNSConfig config;

        /// <summary>
        /// Initializes <see cref="ArgumentOrdering"/>,
        /// <see cref="ArgumentToIndexMap"/> and <see cref="config"/>
        /// </summary>
        /// <param name="config">Configuration options</param>
        public VariableSet(CNSControl.CNSConfig config) {
            this.config = config;
            ArgumentToIndexMap = new Dictionary<string, int>(
                config.NumberOfDimensions + 2);

            ArgumentToIndexMap.Add("rho", 0);
            for (int i = 0; i < config.NumberOfDimensions; i++) {
                ArgumentToIndexMap.Add("p" + i, i + 1);
            }
            ArgumentToIndexMap.Add("E", config.NumberOfDimensions + 1);

            ArgumentOrdering = new string[config.NumberOfDimensions + 2];
            ArgumentToIndexMap.Keys.CopyTo(ArgumentOrdering, 0);
        }

        /// <summary>
        /// Extracts the density value from <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector">
        /// A list of doubles representing the values of the arguments in the
        /// ordering defined by <see cref="ArgumentOrdering"/>
        /// </param>
        /// <returns>The dimensionless density</returns>
        public double GetDensity(double[] argumentVector) {
            CheckVectorLength(argumentVector);
            return argumentVector[ArgumentToIndexMap["rho"]];
        }

        public void SetDensity(double[] argumentVector, double newDensity) {
            CheckVectorLength(argumentVector);
            argumentVector[ArgumentToIndexMap["rho"]] = newDensity;
        }

        /// <summary>
        /// Extracts the values of the momentum components from
        /// <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>The dimensionless momentum vector</returns>
        public double[] GetMomentum(double[] argumentVector) {
            CheckVectorLength(argumentVector);
            
            double[] result = new double[config.NumberOfDimensions];
            for (int i = 0; i < config.NumberOfDimensions; i++) {
                result[i] = argumentVector[ArgumentToIndexMap["p" + i]];
            }

            return result;
        }

        public void SetMomentum(double[] argumentVector, Vector3D newMomentum) {
            CheckVectorLength(argumentVector);

            for (int i = 0; i < config.NumberOfDimensions; i++) {
                argumentVector[ArgumentToIndexMap["p" + i]] = newMomentum[i];
            }
        }

        /// <summary>
        /// Extracts the values of the momentum components from
        /// <paramref name="argumentVector"/> and returns it as a 3D vector
        /// (independent of the spatial dimension of the problem!). Components
        /// that are not needed if the spatial dimensions is smaller than three
        /// will be set to zero.
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>The dimensionless momentum vector</returns>
        public Vector3D GetMomentumVector(double[] argumentVector) {
            CheckVectorLength(argumentVector);
            
            Vector3D result = new Vector3D(0.0, 0.0, 0.0);

            int numberOfDimensions = ArgumentToIndexMap.Count - 2;
            for (int i = 0; i < numberOfDimensions; i++) {
                result[i] = argumentVector[ArgumentToIndexMap["p" + i]];
            }

            return result;
        }

        /// <summary>
        /// Extracts the energy part from <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>The dimensionless total energy per volume</returns>
        public double GetEnergy(double[] argumentVector) {
            CheckVectorLength(argumentVector);
            return argumentVector[ArgumentToIndexMap["E"]];
        }

        public void SetEnergy(double[] argumentVector, double newEnergy) {
            CheckVectorLength(argumentVector);
            argumentVector[ArgumentToIndexMap["E"]] = newEnergy;
        }

        /// <summary>
        /// Calculates the velocity from <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>$\vec{p} / \rho$</returns>
        public Vector3D GetVelocityVector(double[] argumentVector) {
            Vector3D result = GetMomentumVector(argumentVector);
            result.Scale(1 / GetDensity(argumentVector));
            return result;
        }

        public double GetKineticEnergy(double[] argumentVector) {
            Vector3D momentumVector = GetMomentumVector(argumentVector);
            return 0.5 * (momentumVector * momentumVector)
                / GetDensity(argumentVector);
        }

        /// <summary>
        /// Calculates the local inner energy from
        /// <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>$e = \rho E - \fract{1}{2} \fract{p_i p_i}{\rho}$</returns>
        public double GetInnerEnergy(double[] argumentVector) {
            Vector3D momentumVector = GetMomentumVector(argumentVector);
            return GetEnergy(argumentVector) - GetKineticEnergy(argumentVector);
        }

        /// <summary>
        /// Calculates the local pressure from
        /// <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>$(\kappa - 1) * \rho * e$</returns>
        public double GetPressure(double[] argumentVector) {
            return (config.HeatCapacityRatio - 1)
                * GetDensity(argumentVector)
                * GetInnerEnergy(argumentVector);
        }

        /// <summary>
        /// Calculates the local temperature from
        /// <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>$(\kappa - 1) \kappa \mathrm{Ma}^2 e$</returns>
        public double GetTemperature(double[] argumentVector) {
            return (config.HeatCapacityRatio - 1)
                * config.HeatCapacityRatio
                * config.MachNumber
                * config.MachNumber
                * GetInnerEnergy(argumentVector);
        }

        /// <summary>
        /// Calculates the dynamics viscosity from
        /// <paramref name="argumentVector"/> using Sutherland's law.
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>$T^\frac{3}{2} \frac{1 + S}{T + S}$</returns>
        /// <remarks>
        /// The Sutherland temperature S is currently assuming room temperature
        /// in the far field
        /// </remarks>
        public double GetDynamicViscosity(double[] argumentVector) {
            // For now, assume room temperature in the far field
            double S = 110.4 / 293.15;
            double T = GetTemperature(argumentVector);
            return Math.Pow(T, 1.5) * (1 + S) / (1 + T);
        }

        /// <summary>
        /// Calculates the local Mach number from
        /// <paramref name="argumentVector"/>
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        /// <returns>$\kappa \sqrt{(\kappa - 1) e}$</returns>
        public double GetMachNumber(double[] argumentVector) {
            return config.HeatCapacityRatio * Math.Sqrt(
                (config.HeatCapacityRatio - 1) * GetEnergy(argumentVector));
        }

        /// <summary>
        /// Checks whether the given list of doubles is compatible with the
        /// argument ordering defined by this class
        /// </summary>
        /// <param name="argumentVector"><see cref="GetDensity"/></param>
        private void CheckVectorLength(double[] argumentVector) {
            if (argumentVector.Length != ArgumentToIndexMap.Count) {
                throw new ArgumentException("Given vector has an incorrect length", "vector");
            }
        }
    }
}
