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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// [LowMach] Buoyant force.
    /// </summary>
    //public class Buoyancy : BoSSS.Foundation.IVolumeForm {
    public class Buoyancy : LinearSource {
        private Vector GravityDirection;
        private int SpatialComponent;
        private double Froude;
        private MaterialLaw EoS;
        private PhysicsMode physicsMode;
        private string[] m_ParameterOrdering;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="GravityDirection">Unit vector for spatial direction of gravity.</param>
        /// <param name="SpatialComponent">Spatial component of source.</param>
        /// <param name="Froude">Dimensionless Froude number.</param>
        /// <param name="physicsMode"></param>
        /// <param name="EoS">Equation of state for calculating density.</param>
        public Buoyancy(Vector GravityDirection, int SpatialComponent, double Froude, PhysicsMode physicsMode, MaterialLaw EoS) {
            // Check direction
            if ((GravityDirection.Abs() - 1.0) > 1.0e-13)
                throw new ArgumentException("Length of GravityDirection vector has to be 1.0");

            // Initialize
            this.GravityDirection = GravityDirection;
            this.SpatialComponent = SpatialComponent;
            this.Froude = Froude;
            this.EoS = EoS;
            this.physicsMode = physicsMode;

            switch (physicsMode) {
                case PhysicsMode.LowMach:
                    this.m_ParameterOrdering = new string[] { VariableNames.Temperature0 };
                    break;

                case PhysicsMode.Combustion:
                    this.m_ParameterOrdering = new string[] { VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 };
                    break;

                default:
                    throw new ApplicationException("wrong physicsmode");
            }
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="x"></param>
        /// <param name="parameters"></param>
        /// <param name="U"></param>
        /// <returns></returns>
        protected override double Source(double[] x, double[] parameters, double[] U) {
            double src = 0.0;

            //double rho = EoS.GetDensity(U[0]);
            double rho = EoS.GetDensity(parameters);

            src = 1.0 / (Froude * Froude) * rho * GravityDirection[SpatialComponent];

            return src;
        }

        /// <summary>
        /// Temperature
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return new string[] {
                    //VariableNames.Temperature
                };
            }
        }

        /// <summary>
        ///
        /// </summary>
        public override TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        ///
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }
    }

    /// <summary>
    /// [LowMach] Buoyant force.
    /// </summary>

    public class BuoyancyJacobi : IVolumeForm, ISupportsJacobianComponent {
        private Vector GravityDirection;
        private int SpatialComponent;
        private double Froude;
        private MaterialLaw EoS;
        private PhysicsMode physicsMode;
        private string[] m_ParameterOrdering;
        private string[] m_ArgumentOrdering;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="GravityDirection">Unit vector for spatial direction of gravity.</param>
        /// <param name="SpatialComponent">Spatial component of source.</param>
        /// <param name="Froude">Dimensionless Froude number.</param>
        /// <param name="physicsMode"></param>
        /// <param name="EoS">Equation of state for calculating density.</param>
        public BuoyancyJacobi(Vector GravityDirection, int SpatialComponent, double Froude, PhysicsMode physicsMode, MaterialLaw EoS, int noOfChemComponents) {
            // Check direction
            if ((GravityDirection.Abs() - 1.0) > 1.0e-13)
                throw new ArgumentException("Length of GravityDirection vector has to be 1.0");

            // Initialize
            this.GravityDirection = GravityDirection;
            this.SpatialComponent = SpatialComponent;
            this.Froude = Froude;
            this.EoS = EoS;
            this.physicsMode = physicsMode;

            switch (physicsMode) {
                case PhysicsMode.MixtureFraction:
                    this.m_ParameterOrdering = null;
                    this.m_ArgumentOrdering = new string[] { VariableNames.MixtureFraction };
                    break;

                case PhysicsMode.LowMach:
                    this.m_ParameterOrdering = null; // new string[] { VariableNames.Temperature0 };
                    this.m_ArgumentOrdering = new string[] { VariableNames.Temperature };
                    break;

                case PhysicsMode.Combustion:
                    this.m_ParameterOrdering = new string[] { VariableNames.ThermodynamicPressure, VariableNames.Rho, VariableNames.Mu, VariableNames.cp };// null;
                    string[] MFs = VariableNames.MassFractions(noOfChemComponents);
                    this.m_ArgumentOrdering = ArrayTools.Cat(new string[] { VariableNames.Temperature }, MFs);
                    break;

                default:
                    throw new ApplicationException("wrong physicsmode");
            }
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="x"></param>
        /// <param name="parameters"></param>
        /// <param name="U"></param>
        /// <returns></returns>
        protected double Source(double[] x, double[] parameters, double[] U) {
            double src = 0.0;

            double rho;
            switch (physicsMode) {
                case PhysicsMode.MixtureFraction:
                    rho = EoS.getDensityFromZ(U[0]);
                    break;

                case PhysicsMode.LowMach:
                case PhysicsMode.Combustion:
                    rho = EoS.GetDensity(U);
                    break;

                default:
                    throw new NotImplementedException("wrong PhysicsMode");
            }
            src = (1.0 / (Froude * Froude)) * rho * GravityDirection[SpatialComponent];

            return src;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return this.Source(cpv.Xglobal, cpv.Parameters, U) * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var SourceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { SourceDerivVol };

            //return new IEquationComponent[] { this};
        }

        /// <summary>
        /// Temperature
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        /// <summary>
        ///
        /// </summary>
        public virtual TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.AllOn;//TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        ///
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }
    }

    /// <summary>
    /// [LowMach] Buoyant force.
    /// </summary>

    public class DummyParameter : IVolumeForm, ISupportsJacobianComponent {
        private string[] m_ParameterOrdering;
        private string[] m_ArgumentOrdering;

        /// <summary>
        /// Ctor.
        /// </summary>
        public DummyParameter(double dummy) {
            this.m_ParameterOrdering = new string[] { VariableNames.ThermodynamicPressure, VariableNames.Rho};
            this.m_ArgumentOrdering = new string[] { };
        }
        /// <summary>
        /// Ctor.
        /// </summary>
        public DummyParameter() {
            this.m_ParameterOrdering = new string[] { VariableNames.ThermodynamicPressure, VariableNames.Rho, VariableNames.Mu, VariableNames.cp };
            this.m_ArgumentOrdering = new string[] { };
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="x"></param>
        /// <param name="parameters"></param>
        /// <param name="U"></param>
        /// <returns></returns>
        protected double Source(double[] x, double[] parameters, double[] U) {
            return 0.0;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return this.Source(cpv.Xglobal, cpv.Parameters, U) * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var SourceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { SourceDerivVol };
        }

        /// <summary>
        /// Temperature
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        /// <summary>
        ///
        /// </summary>
        public virtual TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.None;
            }
        }

        /// <summary>
        ///
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }
    }
}