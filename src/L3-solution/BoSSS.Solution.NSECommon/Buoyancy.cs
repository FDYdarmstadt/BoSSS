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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Utils;
using ilPSP;

namespace BoSSS.Solution.NSECommon {
    
    /// <summary>
    /// [LowMach] Buoyant force.
    /// </summary>
    //public class Buoyancy : BoSSS.Foundation.IVolumeForm {
    public class Buoyancy : LinearSource { 
        Vector GravityDirection;
        int SpatialComponent;
        double Froude;
        MaterialLaw EoS;
        PhysicsMode physicsMode;
        string[] m_ParameterOrdering;

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

            src =  1.0 / (Froude * Froude) * rho * GravityDirection[SpatialComponent]; 

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
        Vector GravityDirection;
        int SpatialComponent;
        double Froude;
        MaterialLaw EoS;
        PhysicsMode physicsMode;
        string[] m_ParameterOrdering;
        string[] m_ArgumentOrdering;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="GravityDirection">Unit vector for spatial direction of gravity.</param>
        /// <param name="SpatialComponent">Spatial component of source.</param>
        /// <param name="Froude">Dimensionless Froude number.</param>
        /// <param name="physicsMode"></param>
        /// <param name="EoS">Equation of state for calculating density.</param>
        public BuoyancyJacobi(Vector GravityDirection, int SpatialComponent, double Froude, PhysicsMode physicsMode, MaterialLaw EoS) {
            // Check direction
            if((GravityDirection.Abs() - 1.0) > 1.0e-13)
                throw new ArgumentException("Length of GravityDirection vector has to be 1.0");

            // Initialize
            this.GravityDirection = GravityDirection;
            this.SpatialComponent = SpatialComponent;
            this.Froude = Froude;
            this.EoS = EoS;
            this.physicsMode = physicsMode;

            switch(physicsMode) {
                case PhysicsMode.MixtureFraction:
                    this.m_ParameterOrdering = null;
                    this.m_ArgumentOrdering = new string[] { VariableNames.MixtureFraction };
                    break;
                case PhysicsMode.LowMach:
                    this.m_ParameterOrdering = null; // new string[] { VariableNames.Temperature0 };
                    this.m_ArgumentOrdering = new string[] { VariableNames.Temperature };
                    break;
                case PhysicsMode.Combustion:
                    this.m_ParameterOrdering = null;
                    this.m_ArgumentOrdering = new string[] { VariableNames.Temperature, VariableNames.MassFraction0, VariableNames.MassFraction1, VariableNames.MassFraction2, VariableNames.MassFraction3};
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
        public virtual  IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }


    }

}
