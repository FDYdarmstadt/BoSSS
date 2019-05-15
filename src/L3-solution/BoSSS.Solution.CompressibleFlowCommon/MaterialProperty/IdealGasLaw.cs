using System;

namespace CNS.EquationOfState {

    /// <summary>
    /// Ideal gas law for gases like air at moderate condtitions.
    /// </summary>
    public class IdealGasLaw : IEquationOfState {

        /// <summary>
        /// The heat capacity ratio.
        /// </summary>
        private double kappa;

        /// <summary>
        /// <see cref="ICNSOptions.MachNumber"/>
        /// </summary>
        private double referenceMachNumber;

        /// <summary>
        /// Constructs a new ideal gas law
        /// </summary>
        /// <param name="kappa">
        /// The heat capacity ratio
        /// </param>
        /// <param name="referenceMachNumber">
        /// <see cref="ICNSOptions.MachNumber"/>
        /// </param>
        public IdealGasLaw(double kappa, double referenceMachNumber) {
            this.kappa = kappa;
            this.referenceMachNumber = referenceMachNumber;
        }

        /// <summary>
        /// Constructs an ideal gas law for standard air
        /// </summary>
        public static IdealGasLaw Air {
            get {
                return new IdealGasLaw(1.4, double.NaN);
            }
        }

        /// <summary>
        /// Constructs an ideal gas law for Helium
        /// </summary>
        public static IdealGasLaw Helium {
            get {
                return new IdealGasLaw(1.6, double.NaN);
            }
        }

        #region IEquationOfState Members

        /// <summary>
        /// Calculate the pressure $p$ via
        /// $p = (\kappa - 1) \rho e$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetPressure"/>
        /// </param>
        /// <returns>
        /// <latex>\rho e (\kappa - 1)</latex>
        /// where
        /// <latex mode="inline">\rho</latex> = <see name="StateVector.Density"/>,
        /// <latex mode="inline">e</latex> = <see name="StateVector.InnerEnergy"/>
        /// and <latex mode="inline">\kappa</latex> is the heat capacity ratio
        /// supplied to <see cref="IdealGasLaw.IdealGasLaw"/>.
        /// </returns>
        public double GetPressure(StateVector state) {
            return state.Density * (kappa - 1.0) * state.InnerEnergy;
        }

        /// <summary>
        /// Calculates the temperature $T$ via
        /// $T = (\kappa - 1) \kappa \mathrm{Ma}^2 e$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetTemperature"/>
        /// </param>
        /// <returns>
        /// <latex>(\kappa - 1) \kappa \mathrm{Ma}^2 e</latex>
        /// where
        /// <latex mode="inline">e</latex> = <see name="StateVector.InnerEnergy"/>
        /// and <latex mode="inline">\kappa</latex> and
        /// <latex mode="inline">\mathrm{Ma}</latex> are the heat capacity
        /// ratio and the reference Mach number supplied to
        /// <see cref="IdealGasLaw.IdealGasLaw"/>, respectively.
        /// </returns>
        public double GetTemperature(StateVector state) {
            return (kappa - 1.0) * kappa * referenceMachNumber * referenceMachNumber * state.InnerEnergy;
        }

        /// <summary>
        /// Calculates the local speed of sound $a$ via
        /// $\mathrm{a} = \sqrt{\kappa \frac{p}{\rho}}$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetSpeedOfSound"/>
        /// </param>
        /// <returns>
        /// <latex>\sqrt{\kappa \frac{p}{\rho}}</latex>
        /// where
        /// <latex mode="inline">\rho</latex> = <see name="StateVector.Density"/>
        /// <latex mode="inline">p</latex> = <see name="StateVector.Pressure"/>
        /// and <latex mode="inline">\kappa</latex> is the heat capacity ratio
        /// supplied to <see cref="IdealGasLaw.IdealGasLaw"/>.
        /// </returns>
        public double GetSpeedOfSound(StateVector state) {
            return Math.Sqrt(kappa * state.Pressure / state.Density);
        }

        /// <summary>
        /// Calculates the local entropy $S$ via
        /// $S = \frac{p}{\rho^\kappa}$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetEntropy"/>
        /// </param>
        /// <returns>
        /// <latex>S = \frac{p}{\rho^\kappa}</latex>
        /// where
        /// <latex mode="inline">\rho</latex> = <see name="StateVector.Density"/>
        /// <latex mode="inline">p</latex> = <see name="StateVector.Pressure"/>
        /// and <latex mode="inline">\kappa</latex> is the heat capacity ratio
        /// supplied to <see cref="IdealGasLaw.IdealGasLaw"/>.
        /// </returns>
        public double GetEntropy(StateVector state) {
            return state.Pressure / Math.Pow(state.Density, kappa);
        }

        #endregion
    }
}
