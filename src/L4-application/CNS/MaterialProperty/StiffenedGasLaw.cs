using System;

namespace CNS.EquationOfState {

    /// <summary>
    /// Model equation for stiff fluids (e.g. water) where an ideal gas at a
    /// very pressure is assumed. For details, e.g. see ChangLiou2007. The
    /// generic form of the equation is
    /// <latex>
    /// p = (\kappa - 1) \rho e - \pi
    /// </latex>
    /// with determining constants <latex mode="inline">\kappa</latex> and
    /// <latex mode="inline">\pi</latex>.
    /// </summary>
    public class StiffenedGasLaw : IEquationOfState {

        /// <summary>
        /// <see cref="StiffenedGasLaw.StiffenedGasLaw"/>
        /// </summary>
        private readonly double kappa;

        /// <summary>
        /// <see cref="StiffenedGasLaw.StiffenedGasLaw"/>
        /// </summary>
        private readonly double pi;

        /// <summary>
        /// <see cref="StiffenedGasLaw.StiffenedGasLaw"/>
        /// </summary>
        private readonly double C_p;

        /// <summary>
        /// Constructs a stiff fluid.
        /// </summary>
        /// <param name="kappa">
        /// The heat capacity ratio.
        /// </param>
        /// <param name="pi">
        /// The reference pressure. A value of zero corresponds to an ideal
        /// gas.
        /// </param>
        /// <param name="C_p">
        /// The heat capacity at constant pressure. Only affects temperature.
        /// </param>
        public StiffenedGasLaw(double kappa, double pi, double C_p) {
            this.kappa = kappa;
            this.pi = pi;
            this.C_p = C_p;
        }

        /// <summary>
        /// Instantiates the specific equation of state water.
        /// </summary>
        public static StiffenedGasLaw Water {
            get {
                return new StiffenedGasLaw(1.932, 1.1645e9, 8095.08);
            }
        }

        #region IEquationOfState Members

        /// <summary>
        /// Calculates the pressure $p$ via
        /// $p = \rho e (\kappa - 1) - \kappa p_\infty$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetPressure"/>
        /// </param>
        /// <returns>
        /// <latex>\rho e (\kappa - 1) - \kappa p_\infty</latex>
        /// where
        /// <latex mode="inline">\rho</latex> = <see name="StateVector.Density"/>,
        /// <latex mode="inline">e</latex> = <see name="StateVector.InnerEnergy"/>
        /// and <latex mode="inline">\kappa</latex> and
        /// <latex mode="inline">p_\infty</latex> are the heat capacity ratio
        /// and the reference pressure supplied to
        /// <see cref="StiffenedGasLaw.StiffenedGasLaw"/>, respectively.
        /// </returns>
        public double GetPressure(StateVector state) {
            return (kappa - 1.0) * state.Density * state.InnerEnergy - kappa * pi;
        }

        /// <summary>
        /// Calculates the temperature $T$ via
        /// $T = \frac{\kappa}{C_p} (e - \frac{p_\infty}{\rho}$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetTemperature"/>
        /// </param>
        /// <returns>
        /// <latex>\frac{\kappa}{C_p} (e - \frac{p_\infty}{\rho}</latex>
        /// where
        /// <latex mode="inline">\rho</latex> = <see name="StateVector.Density"/>,
        /// <latex mode="inline">e</latex> = <see name="StateVector.InnerEnergy"/>
        /// and <latex mode="inline">\kappa</latex>,
        /// <latex mode="inline">p_\infty</latex> and
        /// <latex mode="inline">C_p</latex> are the heat capacity ratio, the
        /// reference pressure and the heat capacity at constant pressure
        /// supplied to <see cref="StiffenedGasLaw.StiffenedGasLaw"/>, respectively.
        /// </returns>
        public double GetTemperature(StateVector state) {
            return kappa / C_p * (state.InnerEnergy - pi / state.Density);
        }

        /// <summary>
        /// Calculates the speed of sound $a$ via
        /// $a = \sqrt{\kappa \frac{p + \pi}{\rho}}$.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetSpeedOfSound"/>
        /// </param>
        /// <returns>
        /// <latex>\sqrt{\kappa \frac{p + \pi}{\rho}}</latex>
        /// where
        /// <latex mode="inline">\rho</latex> = <see name="StateVector.Density"/>,
        /// <latex mode="inline">p</latex> = <see name="StateVector.Pressure"/>
        /// and <latex mode="inline">\kappa</latex> and
        /// <latex mode="inline">\pi</latex> are the heat capacity ratio
        /// and the reference pressure supplied to
        /// <see cref="StiffenedGasLaw.StiffenedGasLaw"/>, respectively.
        /// </returns>
        public double GetSpeedOfSound(StateVector state) {
            return Math.Sqrt(kappa * (state.Pressure + pi) / state.Density);
        }

        /// <summary>
        /// Calculates the entropy $S$ via
        /// $S = \frac{p + \kappa \pi}{\rho^\kappa}$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetEntropy"/>
        /// </param>
        /// <returns>
        /// <latex>\frac{p + \kappa \pi}{\rho^\kappa}</latex>
        /// where
        /// <latex mode="inline">\rho</latex> = <see name="StateVector.Density"/>,
        /// <latex mode="inline">p</latex> = <see name="StateVector.Pressure"/>
        /// and <latex mode="inline">\kappa</latex> and
        /// <latex mode="inline">\pi</latex> are the heat capacity ratio
        /// and the reference pressure supplied to
        /// <see cref="StiffenedGasLaw.StiffenedGasLaw"/>, respectively.
        /// </returns>
        /// <remarks>
        /// Unverified.
        /// </remarks>
        public double GetEntropy(StateVector state) {
            return (state.Pressure + kappa * pi) / Math.Pow(state.Density, kappa);
        }

        #endregion
    }
}
