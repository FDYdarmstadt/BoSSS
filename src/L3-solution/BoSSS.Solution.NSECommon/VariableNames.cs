﻿/* =======================================================================
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

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Options for different physical applications.
    /// </summary>
    public enum PhysicsMode {

        /// <summary>
        /// Incompressible single-phase flows,
        /// i.e. solving incompressible continuity and momentum equation.
        /// </summary>
        Incompressible,

        /// <summary>
        /// Low Mach number single-phase flows,
        /// i.e. solving low Mach equations for
        /// mass, momentum and temperature.
        /// </summary>
        LowMach,

        /// <summary>
        /// Low Mach number single-phase flows with Combustion,
        /// i.e. solving low Mach equations for
        /// mass, momentum, temperature and mass fractions of 0 and 1.
        /// </summary>
        Combustion,

        /// <summary>
        /// Incompressible multiphase flows,
        /// i.e. solving equations for
        /// mass, momentum and level-set.
        /// </summary>
        Multiphase,


        Helical,

        /// <summary>
        /// Incompressible viscoelastic flows,
        /// i.e. solving equations for
        /// mass, momentum and constitutive 
        /// equations for viscoelastic material model.
        /// </summary>
        Viscoelastic,

        /// <summary>
        /// Reynolds-averaged Navier-Stokes (RANS)
        /// solver.
        /// Different Turbulence models will be supported in the future.
        /// </summary>
        RANS,


        /// <summary>
        ///  MixtureFraction solver.
        /// Used as a pre-step for calculating reactive flows.
        /// Equations for continuity, momentum and for a passive scalar Z are solved,
        /// while the Temperature and mass fraction fields are approximated by the Burke-Schumann limit
        /// </summary>
        MixtureFraction,


    }


    /// <summary>
    /// predefined equations names, that should be used by all Solvers
    /// </summary>
    public static class EquationNames {

        /// <summary>
        /// momentum equation component in x - direction; 
        /// </summary>
        public const string MomentumEquationX = "MomentumX";

        /// <summary>
        /// momentum equation component in y - direction; 
        /// </summary>
        public const string MomentumEquationY = "MomentumY";

        /// <summary>
        /// momentum equation component in z - direction; 
        /// </summary>
        public const string MomentumEquationZ = "MomentumZ";

        /// <summary>
        /// The name of the <paramref name="d"/>-th momentum equation component.
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D, {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string MomentumEquationComponent(int d) {
            switch (d) {
                case 0: return MomentumEquationX;
                case 1: return MomentumEquationY;
                case 2: return MomentumEquationZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }

        /// <summary>
        /// vector of momentum equation names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] MomentumEquations(int D) {
            if (D == 2)
                return new string[] { MomentumEquationX, MomentumEquationY };
            else if (D == 3)
                return new string[] { MomentumEquationX, MomentumEquationY, MomentumEquationZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// continuity equation
        /// </summary>
        public const string ContinuityEquation = "ContiEq";

        /// <summary>
        /// kinetic energy equation
        /// </summary>
        public const string KineticEnergyEquation = "KinEnergyEq";

        /// <summary>
        /// heat equation
        /// </summary>
        public const string HeatEquation = "HeatEq";


        /// <summary>
        /// heat equation
        /// </summary>
        public const string MixtureFractionEquation = "MixtureFractionEquation";



        /// <summary>
        /// auxiliary heat flux x - component
        /// </summary>
        public const string AuxHeatFluxX = "AuxHeatFluxX";

        /// <summary>
        /// auxiliary heat flux y - component
        /// </summary>
        public const string AuxHeatFluxY = "AuxHeatFluxY";

        /// <summary>
        /// auxiliary heat flux z - component
        /// </summary>
        public const string AuxHeatFluxZ = "AuxHeatFluxZ";

        /// <summary>
        /// 
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        static public string AuxHeatFluxComponent(int d) {
            switch (d) {
                case 0: return AuxHeatFluxX;
                case 1: return AuxHeatFluxY;
                case 2: return AuxHeatFluxZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }

        /// <summary>
        /// vector of the auxiliary heat flux
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] AuxHeatFlux(int D) {
            if (D == 2)
                return new string[] { AuxHeatFluxX, AuxHeatFluxY };
            else if (D == 3)
                return new string[] { AuxHeatFluxX, AuxHeatFluxY, AuxHeatFluxZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// Mass balance equation of chemical component 0 (usually fuel)
        /// </summary>
        public const string SpeciesMassBalance0 = "SpeciesMassBalance0";

        /// <summary>
        /// Mass balance equation of chemical component 1 (usually oxidizer)
        /// </summary>
        public const string SpeciesMassBalance1 = "SpeciesMassBalance1";

        /// <summary>
        /// Mass balance equation of chemical component 2
        /// </summary>
        public const string SpeciesMassBalance2 = "SpeciesMassBalance2";

        /// <summary>
        /// Mass balance equation of chemical component 3
        /// </summary>
        public const string SpeciesMassBalance3 = "SpeciesMassBalance3";

        /// <summary>
        /// Mass balance equation of chemical component  4
        /// </summary>
        public const string SpeciesMassBalance4 = "SpeciesMassBalance4";

        /// <summary>
        /// 
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        static public string SpeciesMassBalanceName(int component) {
            switch (component) {
                case 0: return SpeciesMassBalance0;
                case 1: return SpeciesMassBalance1;
                case 2: return SpeciesMassBalance2;
                case 3: return SpeciesMassBalance3;
                case 4: return SpeciesMassBalance4;
                default: throw new NotSupportedException("Name of equation for component " + component + " not supported yet");
            }
        }

        /// <summary>
        /// vector of the auxiliary heat flux
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] SpeciesMassBalanceNames(int NoOfComponents) {
            switch (NoOfComponents) {
                case 3:
                return new string[] { SpeciesMassBalance0, SpeciesMassBalance1, SpeciesMassBalance2 };
                case 4:
                return new string[] { SpeciesMassBalance0, SpeciesMassBalance1, SpeciesMassBalance2, SpeciesMassBalance3 };
                case 5:
                return new string[] { SpeciesMassBalance0, SpeciesMassBalance1, SpeciesMassBalance2, SpeciesMassBalance3, SpeciesMassBalance4 };
                default:
                throw new NotImplementedException("Solver for" + NoOfComponents + "components not supported!");
            }
        }



        /// <summary>
        /// xx - component constitutive equation for viscoelastic extra stress
        /// </summary>
        public const string ConstitutiveXX = "ConstitutiveXX";

        /// <summary>
        /// xy - component constitutive equation for viscoelastic extra stress
        /// </summary>
        public const string ConstitutiveXY = "ConstitutiveXY";

        /// <summary>
        /// yy - component constitutive equation for viscoelastic extra stress
        /// </summary>
        public const string ConstitutiveYY = "ConstitutiveYY";


        /// <summary>
        /// vector of the auxiliary heat flux
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] Constitutive(int D)
        {
            if (D == 2)
                return new string[] { ConstitutiveXX, ConstitutiveXY, ConstitutiveYY };
            else if (D == 3)
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ". 3D-case not implemented yet!");
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }
    }

    
    /// <summary>
    /// predefined variable names for velocity, pressure, etc. that should be used
    /// by all Navier-Stokes - Solvers
    /// </summary>
    /// <remarks>
    /// Follows the CGNS naming conventions, whenever possible (see http://www.grc.nasa.gov/WWW/cgns/CGNS_docs_rel31/sids/dataname.html);
    /// </remarks>
    public static class VariableNames {

        /// <summary>
        /// molecular viscosity, resp. dynamic viscosity
        /// </summary>
        public const string ViscosityMolecular = "ViscosityMolecular";

        /// <summary>
        /// velocity component in x - direction; see also <see cref="VelocityVector"/> and <see cref="Velocity_d"/>.
        /// </summary>
        public const string VelocityX = "VelocityX";

        /// <summary>
        /// velocity component in y - direction; see also <see cref="VelocityVector"/> and <see cref="Velocity_d"/>.
        /// </summary>
        public const string VelocityY = "VelocityY";

        /// <summary>
        /// velocity component in z - direction; see also <see cref="VelocityVector"/> and <see cref="Velocity_d"/>.
        /// </summary>
        public const string VelocityZ = "VelocityZ";

        /// <summary>
        /// Gravity/volume force component in x - direction; see also <see cref="GravityVector"/> and <see cref="Gravity_d"/>.
        /// </summary>
        public const string GravityX = "GravityX";

        /// <summary>
        /// Gravity/volume force component in y - direction; see also <see cref="GravityVector"/> and <see cref="Gravity_d"/>.
        /// </summary>
        public const string GravityY = "GravityY";

        /// <summary>
        /// Gravity/volume force component in z - direction; see also <see cref="GravityVector"/> and <see cref="Gravity_d"/>.
        /// </summary>
        public const string GravityZ = "GravityZ";

        /// <summary>
        /// Volume force component in x - direction; see also <see cref="VolumeForceVector"/>.
        /// </summary>
        public const string VolumeForceX = "VolumeForceX";

        /// <summary>
        /// Volume force component in y - direction; see also <see cref="VolumeForceVector"/>.
        /// </summary>
        public const string VolumeForceY = "VolumeForceY";

        /// <summary>
        /// Volume force component in z - direction; see also <see cref="VolumeForceVector"/>.
        /// </summary>
        public const string VolumeForceZ = "VolumeForceZ";

        /// <summary>
        /// vorticity component in x - direction; 
        /// </summary>
        public const string VorticityX = "VorticityX";

        /// <summary>
        /// vorticity component in y - direction; 
        /// </summary>
        public const string VorticityY = "VorticityY";

        /// <summary>
        /// vorticity component in z - direction; 
        /// </summary>
        public const string VorticityZ = "VorticityZ";

        /// <summary>
        /// vorticity in a 2D domain (only one component)
        /// </summary>
        public const string Vorticity2D = "Vorticity";


        /// <summary>
        /// convective part, x-component, i.e. \f$ \nabla \cdot ( u_1 \cdot \vec{u} ) \f$ 
        /// </summary>
        public const string ConvectiveX = "ConvectiveX";
        
        /// <summary>
        /// convective part, y-component, i.e. \f$ \nabla \cdot ( u_2 \cdot \vec{u} ) \f$ 
        /// </summary>
        public const string ConvectiveY = "ConvectiveY";

        /// <summary>
        /// convective part, z-component, i.e. \f$ \nabla \cdot ( u_3 \cdot \vec{u} ) \f$ 
        /// </summary>
        public const string ConvectiveZ = "ConvectiveZ";


        /// <summary>
        /// velocity magnitude
        /// </summary>
        public const string VelocityMagnitude = "VelocityMagnitude";

        /// <summary>
        /// vorticity magnitude
        /// </summary>
        public const string VorticityMagnitude = "VorticityMagnitude";

        /// <summary>
        /// intermediate velocity component of a splitting scheme in x - direction; 
        /// </summary>
        public const string VelocityIntermedX = "VelocityX*";

        /// <summary>
        /// intermediate velocity component of a splitting scheme in y - direction; 
        /// </summary>
        public const string VelocityIntermedY = "VelocityY*";

        /// <summary>
        /// intermediate velocity component of a splitting scheme in z - direction; 
        /// </summary>
        public const string VelocityIntermedZ = "VelocityZ*";

        /// <summary>
        /// correction velocity component of a splitting scheme in x - direction, e.g. in the SIMPLE algorithm; 
        /// </summary>
        public const string VelocityCorX = "VelocityX'";

        /// <summary>
        /// correction velocity component of a splitting scheme in y - direction, e.g. in the SIMPLE algorithm; 
        /// </summary>
        public const string VelocityCorY = "VelocityY'";

        /// <summary>
        /// correction velocity component of a splitting scheme in z - direction, e.g. in the SIMPLE algorithm; 
        /// </summary>
        public const string VelocityCorZ = "VelocityZ'";

        /// <summary>
        /// velocity (linearization point) component in x - direction;
        /// </summary>
        public const string Velocity0X = "Velocity0X";

        /// <summary>
        /// velocity (linearization point) component in y - direction; 
        /// </summary>
        public const string Velocity0Y = "Velocity0Y";

        /// <summary>
        /// velocity (linearization point) component in z - direction; 
        /// </summary>
        public const string Velocity0Z = "Velocity0Z";

        /// <summary>
        /// cell-wise mean value of velocity (linearization point) component in x - direction; see also <see cref="VelocityVector"/>.
        /// </summary>
        public const string Velocity0X_Mean = "Velocity0X_Mean";

        /// <summary>
        /// cell-wise mean value of velocity (linearization point) component in y - direction; see also <see cref="VelocityVector"/>.
        /// </summary>
        public const string Velocity0Y_Mean = "Velocity0Y_Mean";

        /// <summary>
        /// cell-wise mean value of velocity (linearization point) component in z - direction; see also <see cref="VelocityVector"/>.
        /// </summary>
        public const string Velocity0Z_Mean = "Velocity0Z_Mean";

        /// <summary>
        /// cell-wise mean value of momentum (linearization point) component in x - direction; used in LowMach SIMPLE.
        /// </summary>
        public const string Momentum0X_Mean = "Momentum0X_Mean";

        /// <summary>
        /// cell-wise mean value of momentum (linearization point) component in y - direction; used in LowMach SIMPLE.
        /// </summary>
        public const string Momentum0Y_Mean = "Momentum0Y_Mean";

        /// <summary>
        /// cell-wise mean value of momentum (linearization point) component in z - direction; used in LowMach SIMPLE.
        /// </summary>
        public const string Momentum0Z_Mean = "Momentum0Z_Mean";

        /// <summary>
        /// boundary value for velocity component in x - direction; see also <see cref="BoundaryVelocityVector"/> and <see cref="BoundaryVelocity_d"/>.
        /// </summary>
        public const string BoundaryVelocityX = "BoundaryVelocityX";

        /// <summary>
        /// boundary value for velocity component in y - direction; see also <see cref="BoundaryVelocityVector"/> and <see cref="BoundaryVelocity_d"/>.
        /// </summary>
        public const string BoundaryVelocityY = "BoundaryVelocityY";

        /// <summary>
        /// boundary value for velocity component in z - direction; see also <see cref="BoundaryVelocityVector"/> and <see cref="BoundaryVelocity_d"/>.
        /// </summary>
        public const string BoundaryVelocityZ = "BoundaryVelocityZ";

        /// <summary>
        /// gradient in x-direction of velocity component in x-direction
        /// </summary>
        public const string VelocityX_GradientX = "VelocityX_GradientX";

        /// <summary>
        /// gradient in y-direction of velocity component in x-direction
        /// </summary>
        public const string VelocityX_GradientY = "VelocityX_GradientY";

        /// <summary>
        /// gradient in z-direction of velocity component in x-direction
        /// </summary>
        public const string VelocityX_GradientZ = "VelocityX_GradientZ";

        /// <summary>
        /// gradient in x-direction of velocity component in y-direction
        /// </summary>
        public const string VelocityY_GradientX = "VelocityY_GradientX";

        /// <summary>
        /// gradient in y-direction of velocity component in y-direction
        /// </summary>
        public const string VelocityY_GradientY = "VelocityY_GradientY";

        /// <summary>
        /// gradient in z-direction of velocity component in y-direction
        /// </summary>
        public const string VelocityY_GradientZ = "VelocityY_GradientZ";

        /// <summary>
        /// gradient in x-direction of velocity component in z-direction
        /// </summary>
        public const string VelocityZ_GradientX = "VelocityZ_GradientX";

        /// <summary>
        /// gradient in y-direction of velocity component in z-direction
        /// </summary>
        public const string VelocityZ_GradientY = "VelocityZ_GradientY";

        /// <summary>
        /// gradient in z-direction of velocity component in z-direction
        /// </summary>
        public const string VelocityZ_GradientZ = "VelocityZ_GradientZ";


        ///// <summary>
        ///// variable name for a single component of the veolcityY gradient
        ///// </summary>
        //public static string VelocityYGradient(string Component) { return "VelocityYGradient" + Component; ; }

        //public static string[] VelocityXGradient() { return new string[] { VelocityXGradient("X"), VelocityXGradient("Y") }; }
        ///// <summary>
        ///// variable name for the velocity gradient of v
        ///// </summary>
        //static public string[] VelocityGradient() {
        //    return new string[] { VelocityGradient("X", "X"), VelocityGradient("X", "Y"), VelocityGradient("Y", "X"), VelocityGradient("Y", "Y") };
        //}
        ///// <summary>
        ///// variable name for a single component of the veolcity gradient of V
        ///// </summary>
        //public static string VelocityGradient(string VelocityComponent, string Derivative) { return "Velocity"+VelocityComponent+ "Gradient"+Derivative; }

        ///// <summary>
        ///// variable name for the velocity gradient of v
        ///// </summary>
        //static public string[] VelocityGradient() {
        //    return new string[] {  VelocityGradient("X", "X"), VelocityGradient("X", "Y"), VelocityGradient("Y", "X"), VelocityGradient("Y", "Y") };
        //}


        /// <summary>
        /// variable name for kinetic energy
        /// </summary>
        public const string KineticEnergy = "KineticEnergy";


        /// <summary>
        /// variable name for Pressure
        /// </summary>
        public const string Pressure = "Pressure";

        /// <summary>
        /// variable name for Pressure parameter
        /// </summary>
        public const string Pressure0 = "Pressure0";


        /// <summary>
        /// x - component of the pressure gradient
        /// </summary>
        public const string PressureGradX = "PressureGradX";

        /// <summary>
        /// y - component of the pressure gradient
        /// </summary>
        public const string PressureGradY = "PressureGradY";

        /// <summary>
        /// z - component of the pressure gradient
        /// </summary>
        public const string PressureGradZ = "PressureGradZ";

        /// <summary>
        /// variable name for the pressure gradient
        /// </summary>
        static public string[] PressureGradient(int D) {
            switch (D) {
                case 1: return new string[] { PressureGradX };
                case 2: return new string[] { PressureGradX, PressureGradY };
                case 3: return new string[] { PressureGradX, PressureGradY, PressureGradZ };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// Components of the pressure gradient
        /// </summary>
        static public string PressureGradientComponent(int i) {
            switch (i) {
                case 0: return PressureGradX;
                case 1: return PressureGradY;
                case 2: return PressureGradZ;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }


        /// <summary>
        /// variable name for Pressure Correction (e.g. in SIMPLE)
        /// </summary>
        public const string PresCor = "Pressure'";

        /// <summary>
        /// variable for thermodynamic pressure in Low-Mach flows, i.e. p0.
        /// </summary>
        public const string ThermodynamicPressure = "ThermodynamicPressure";

        /// <summary>
        /// variable name for density
        /// </summary>
        public const string Rho = "Density";


        /// <summary>
        /// variable name for density
        /// </summary>
        public const string Mu = "Viscosity";

        /// <summary>
        /// variable name for the mixture heat capacity
        /// </summary>
        public const string cp = "MixtureHeatCapacity";

        /// <summary>
        /// variable name for temperature
        /// </summary>
        public const string Temperature = "Temperature";


        /// <summary>
        /// Phoretic field/concentration
        /// </summary>
        public const string Phoretic = "PhoreticField";


        /// <summary>
        /// variable name for temperature (linearization point)
        /// </summary>
        public const string Temperature0 = "Temperature0";

        /// <summary>
        /// variable name for mean value of temperature (linearization point)
        /// </summary>
        public const string Temperature0Mean = "Temperature0Mean";

        /// <summary>
        /// x - component of the Temperature gradient
        /// </summary>
        public const string TemperatureGradient0 = "TemperatureGradient[0]";

        /// <summary>
        /// y - component of the Temperature gradient
        /// </summary>
        public const string TemperatureGradient1 = "TemperatureGradient[1]";

        /// <summary>
        /// z - component of the Temperature gradient
        /// </summary>
        public const string TemperatureGradient2 = "TemperatureGradient[2]";
         

        /// <summary>
        /// variable name for the Gradient of the Temperature
        /// </summary>
        static public string[] TemperatureGradient(int D) {
            switch (D) {
                case 1: return new string[] { TemperatureGradient0 };
                case 2: return new string[] { TemperatureGradient0, TemperatureGradient1 };
                case 3: return new string[] { TemperatureGradient0, TemperatureGradient1, TemperatureGradient2 };
                default: throw new NotSupportedException("unsupported spatial dimension.");
            }
        }

        /// <summary>
        /// Components of the Temperature gradient
        /// </summary>
        static public string TemperatureGradientComponent(int i) {
            switch (i) {
                case 0: return TemperatureGradient0;
                case 1: return TemperatureGradient1;
                case 2: return TemperatureGradient2;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// x - component of the heat flux
        /// </summary>
        public const string HeatFluxX = "HeatFluxX";

        /// <summary>
        /// y - component of the  heat flux
        /// </summary>
        public const string HeatFluxY = "HeatFluxY";

        /// <summary>
        /// z - component of the  heat flux
        /// </summary>
        public const string HeatFluxZ = "HeatFluxZ";

        /// <summary>
        /// variable name for the heat flux vector
        /// </summary>
        static public string[] HeatFluxVector(int D) {
            switch (D) {
                case 1: return new string[] { HeatFluxX };
                case 2: return new string[] { HeatFluxX, HeatFluxY };
                case 3: return new string[] { HeatFluxX, HeatFluxY, HeatFluxZ };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// Components of the heat flux vector
        /// </summary>
        static public string HeatFluxVectorComponent(int i) {
            switch (i) {
                case 0: return HeatFluxX;
                case 1: return HeatFluxY;
                case 2: return HeatFluxZ;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// x - component of the heat flux
        /// </summary>
        public const string ResidualAuxHeatFluxX = "ResAuxHeatFluxX";

        /// <summary>
        /// y - component of the  heat flux
        /// </summary>
        public const string ResidualAuxHeatFluxY = "ResAuxHeatFluxY";

        /// <summary>
        /// z - component of the  heat flux
        /// </summary>
        public const string ResidualAuxHeatFluxZ = "ResAuxHeatFluxZ";

        /// <summary>
        /// variable name for the heat flux vector
        /// </summary>
        static public string[] ResidualAuxHeatFluxVector(int D) {
            switch (D) {
                case 1: return new string[] { ResidualAuxHeatFluxX };
                case 2: return new string[] { ResidualAuxHeatFluxX, ResidualAuxHeatFluxY };
                case 3: return new string[] { ResidualAuxHeatFluxX, ResidualAuxHeatFluxY, ResidualAuxHeatFluxZ };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// Components of the heat flux vector
        /// </summary>
        static public string ResidualAuxHeatFluxVectorComponent(int i) {
            switch (i) {
                case 0: return ResidualAuxHeatFluxX;
                case 1: return ResidualAuxHeatFluxY;
                case 2: return ResidualAuxHeatFluxZ;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }


        /// <summary>
        /// x - component of the heat flux
        /// </summary>
        public const string HeatFlux0X = "HeatFlux0X";

        /// <summary>
        /// y - component of the  heat flux
        /// </summary>
        public const string HeatFlux0Y = "HeatFlux0Y";

        /// <summary>
        /// z - component of the  heat flux
        /// </summary>
        public const string HeatFlux0Z = "HeatFlux0Z";

        /// <summary>
        /// variable name for the volumetric heat source
        /// </summary>
        static public string HeatSource = "HeatSource";
        
        /// <summary>
        /// variable name for the heat flux vector
        /// </summary>
        static public string[] HeatFlux0Vector(int D) {
            switch (D) {
                case 1: return new string[] { HeatFlux0X };
                case 2: return new string[] { HeatFlux0X, HeatFlux0Y };
                case 3: return new string[] { HeatFlux0X, HeatFlux0Y, HeatFlux0Z };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// Components of the heat flux vector
        /// </summary>
        static public string HeatFlux0VectorComponent(int i) {
            switch (i) {
                case 0: return HeatFlux0X;
                case 1: return HeatFlux0Y;
                case 2: return HeatFlux0Z;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        static public string MassFluxExtension = "MassFluxExtension";
               

        /// <summary>
        /// for XNSE - discontinuous level set field (used for level set evolution)
        /// </summary>
        public const string LevelSetDG = "PhiDG";

        /// <summary>
        /// variable name for the phase dividing interface
        /// </summary>
        public const string LevelSetCG = "Phi";

        /// <summary>
        /// level-set names if more than one level-set is used.
        /// </summary>
        static public string LevelSetDGidx(int iLevSet) {
            if(iLevSet < 0 || iLevSet >= 4)
                throw new ArgumentOutOfRangeException("Invalid level-set index: " + iLevSet);
            if(iLevSet == 0)
                return LevelSetDG;
            else
                return $"Phi{iLevSet+1}DG";
        }

        /// <summary>
        /// level-set names if more than one level-set is used.
        /// </summary>
        static public string LevelSetCGidx(int iLevSet) {
            if(iLevSet < 0 || iLevSet >= 4)
                throw new ArgumentOutOfRangeException("Invalid level-set index: " + iLevSet);
            if(iLevSet == 0)
                return LevelSetCG;
            else
                return $"Phi{iLevSet+1}";
        }


        /// <summary>
        /// variable name for a single Level Set
        /// </summary>
        public const string LevelSet = "LevelSet";

        /// <summary>
        /// convention for variable names of properties defined at the interface 
        /// </summary>
        public static string AsLevelSetVariable(string levelSetName, string variable)
        {
            return variable + "@" + levelSetName;
        }

        /// <summary>
        /// convention for variable names of properties defined at the interface 
        /// </summary>
        public static IList<string> AsLevelSetVariable(string levelSetName, IList<string> variables)
        {
            for (int i = 0; i < variables.Count; ++i)
            {
                variables[i] += "@" + levelSetName;
            }
            return variables;
        }

        /// <summary>
        /// x - component of the Level-Set gradient
        /// </summary>
        public const string LevelSetGradient0 = "LevelSetGradient[0]";

        /// <summary>
        /// y - component of the Level-Set gradient
        /// </summary>
        public const string LevelSetGradient1 = "LevelSetGradient[1]";

        /// <summary>
        /// z - component of the Level-Set gradient
        /// </summary>
        public const string LevelSetGradient2 = "LevelSetGradient[2]";

        /// <summary>
        /// variable name for the Gradient of a Level-Set
        /// </summary>
        static public string[] LevelSetGradient(int D) {
            switch (D) {
                case 1: return new string[] { LevelSetGradient0 };
                case 2: return new string[] { LevelSetGradient0, LevelSetGradient1 };
                case 3: return new string[] { LevelSetGradient0, LevelSetGradient1, LevelSetGradient2 };
                default: throw new NotSupportedException("unsupported spatial dimension.");
            }
        }

        /// <summary>
        /// Components of the Level-Set gradient
        /// </summary>
        static public string LevelSetGradientComponent(int i) {
            switch (i) {
                case 0: return LevelSetGradient0;
                case 1: return LevelSetGradient1;
                case 2: return LevelSetGradient2;
                default: throw new NotSupportedException("unsupported spatial dimension.");
            }
        }


        /// <summary>
        /// Interface normal component in x - direction;
        /// </summary>
        public const string NormalX = "NormalX";

        /// <summary>
        /// Interface normal component in y - direction;
        /// </summary>
        public const string NormalY = "NormalY";

        /// <summary>
        /// Interface normal component in z - direction;
        /// </summary>
        public const string NormalZ = "NormalZ";

        /// <summary>
        /// mean curvature of the interface
        /// </summary>
        public const string Curvature = "Curvature";

        /// <summary>
        /// artificial surface force (usually only used in manufactured solutions) - x component
        /// </summary>
        public const string SurfaceForceX = "SurfaceForceX";

        /// <summary>
        /// artificial surface force (usually only used in manufactured solutions) - y component
        /// </summary>
        public const string SurfaceForceY = "SurfaceForceY";

        /// <summary>
        /// artificial surface force (usually only used in manufactured solutions) - z component
        /// </summary>
        public const string SurfaceForceZ = "SurfaceForceZ";

        /// <summary>
        /// additional pressure near walls due to intermolecular forces
        /// </summary>
        public const string DisjoiningPressure = "DisjoiningPressure";

        /// <summary>
        /// variable name for the interface normal
        /// </summary>
        static public string[] NormalVector(int D) {
            switch (D) {
                case 1: return new string[] { NormalX };
                case 2: return new string[] { NormalX, NormalY };
                case 3: return new string[] { NormalX, NormalY, NormalZ };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// Components of the interface normal
        /// </summary>
        static public string NormalVectorComponent(int i) {
            switch (i) {
                case 0: return NormalX;
                case 1: return NormalY;
                case 2: return NormalZ;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// variable name for the surface force
        /// </summary>
        static public string[] SurfaceForceVector(int D) {
            switch (D) {
                case 1: return new string[] { SurfaceForceX };
                case 2: return new string[] { SurfaceForceX, SurfaceForceY };
                case 3: return new string[] { SurfaceForceX, SurfaceForceY, SurfaceForceZ };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// Components of the surface force
        /// </summary>
        static public string SurfaceForceComponent(int i) {
            switch (i) {
                case 0: return SurfaceForceX;
                case 1: return SurfaceForceY;
                case 2: return SurfaceForceZ;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }


        /// <summary>
        ///  x - component of cell-wise mean value of the Level-Set gradient
        /// </summary>
        public static string MeanLevelSetGradient0 = "LevelSetGradient[0]Mean";

        /// <summary>
        ///  y - component of cell-wise mean value of the Level-Set gradient
        /// </summary>
        public static string MeanLevelSetGradient1 = "LevelSetGradient[1]Mean";
        
        /// <summary>
        ///  z - component of cell-wise mean value of the Level-Set gradient
        /// </summary>
        public static string MeanLevelSetGradient2 = "LevelSetGradient[2]Mean";


        /// <summary>
        /// variable name for the cell-wise mean value of the Level-Set gradient
        /// </summary>
        static public string[] MeanLevelSetGradient(int D) {
            switch (D) {
                case 1: return new string[] { MeanLevelSetGradient0 };
                case 2: return new string[] { MeanLevelSetGradient0, MeanLevelSetGradient1 };
                case 3: return new string[] { MeanLevelSetGradient0, MeanLevelSetGradient1, MeanLevelSetGradient2 };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }        

        /// <summary>
        /// components of cell-wise mean value of the Level-Set gradient
        /// </summary>
        static public string MeanLevelSetGradientComponent(int i) {
            switch (i) {
                case 0: return  MeanLevelSetGradient0 ;
                case 1: return  MeanLevelSetGradient1 ;
                case 2: return  MeanLevelSetGradient2 ;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }


        /// <summary>
        /// variable names for multiple Level Sets
        /// </summary>
        static public string Phi_n(int n = 0) {
            return "LevelSet_" + n;
        }

        /// <summary>
        /// variable name for scalar (linearization point)
        /// </summary>
        public const string Phi0 = "Phi0";

        /// <summary>
        /// variable name for mean value of scalar (linearization point)
        /// </summary>
        public const string Phi0Mean = "Phi0Mean";



        /// <summary>
        /// x - component of the extension velocity
        /// </summary>
        public const string ExtensionVelocityX = "ExtensionVelocityX";

        /// <summary>
        /// x - component of the extension velocity
        /// </summary>
        public const string ExtensionVelocityY = "ExtensionVelocityY";

        /// <summary>
        /// x - component of the extension velocity
        /// </summary>
        public const string ExtensionVelocityZ = "ExtensionVelocityZ";

        /// <summary>
        /// variable name of the extension velocity
        /// </summary>
        static public string[] ExtensionVelocity(int D) {
            switch (D) {
                case 1: return new string[] { ExtensionVelocityX };
                case 2: return new string[] { ExtensionVelocityX, ExtensionVelocityY };
                case 3: return new string[] { ExtensionVelocityX, ExtensionVelocityY, ExtensionVelocityZ };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// Components of the extension velocity 
        /// </summary>
        static public string ExtensionVelocityComponent(int d) {
            switch (d) {
                case 0: return ExtensionVelocityX;
                case 1: return ExtensionVelocityY;
                case 2: return ExtensionVelocityZ;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }


        /// <summary>
        /// variable name for the MassFraction of component 0
        /// </summary>
        public const string MassFraction0 = "MassFraction0";
        /// <summary>
        /// variable name for the mixture fraction (z) 
        /// </summary>
        public const string MixtureFraction = "MixtureFraction";

        /// <summary>
        /// variable name for the MassFraction of component 0 at linearization point
        /// </summary>
        public const string MassFraction0_0 = "MassFraction0_0";

        /// <summary>
        /// variable name for the MassFractionMean of component 0
        /// </summary>
        public const string MassFraction0Mean= "MassFraction0Mean";

        /// <summary>
        /// variable name for MassFraction of component 1
        /// </summary>
        public const string MassFraction1 = "MassFraction1";

        /// <summary>
        /// variable name for MassFraction of component 1 at linearization point
        /// </summary>
        public const string MassFraction1_0 = "MassFraction1_0";

        /// <summary>
        /// variable name for the MassFractionMean of component 1
        /// </summary>
        public const string MassFraction1Mean = "MassFraction1Mean";

        /// <summary>
        /// variable name for MassFraction of component  2 
        /// </summary>
        public const string MassFraction2 = "MassFraction2";

        /// <summary>
        /// variable name for MassFraction of component  2 at linearization point
        /// </summary>
        public const string MassFraction2_0 = "MassFraction2_0";

        /// <summary>
        /// variable name for the MassFractionMean of component 2
        /// </summary>
        public const string MassFraction2Mean = "MassFraction2Mean";

        /// <summary>
        /// variable name for MassFraction of component 3
        /// </summary>
        public const string MassFraction3 = "MassFraction3";

        /// <summary>
        /// variable name for MassFraction of component 3 at linearization point
        /// </summary>
        public const string MassFraction3_0 = "MassFraction3_0";

        /// <summary>
        /// variable name for the MassFractionMean of component 3
        /// </summary>
        public const string MassFraction3Mean = "MassFraction3Mean";

        /// <summary>
        /// variable name for MassFraction of component 4
        /// </summary>
        public const string MassFraction4 = "MassFraction4";

        /// <summary>
        /// variable name for MassFraction of component 3 at linearization point
        /// </summary>
        public const string MassFraction4_0 = "MassFraction4_0";

        /// <summary>
        /// variable name for the MassFractionMean of component 3
        /// </summary>
        public const string MassFraction4Mean = "MassFraction4Mean";

        /// <summary>
        /// variable name for the arrhenius reaction term 
        /// </summary>
        public const string ReactionRate = "ReactionRate";


        /// <summary>
        /// the names of all mass fractions listed in an array
        /// </summary>
        static public string[] MassFractions(int NumberOfSpecies) {
            switch (NumberOfSpecies) {
                case 1: return new string[] { MassFraction0};
                case 2: return new string[] { MassFraction0, MassFraction1 };
                case 3: return new string[] { MassFraction0, MassFraction1, MassFraction2 };
                case 4: return new string[] { MassFraction0, MassFraction1, MassFraction2, MassFraction3 };
                case 5: return new string[] { MassFraction0, MassFraction1, MassFraction2, MassFraction3, MassFraction4 };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// the names of all mass fractions from last timestep listed in an array
        /// </summary>
        static public string[] MassFractions_t0(int NumberOfSpecies) {
            switch (NumberOfSpecies) {
                case 1: return new string[] { MassFraction0 + "_t0" };
                case 2: return new string[] { MassFraction0 + "_t0", MassFraction1 + "_t0" };
                case 3: return new string[] { MassFraction0 + "_t0", MassFraction1 + "_t0", MassFraction2 + "_t0" };
                case 4: return new string[] { MassFraction0 + "_t0", MassFraction1 + "_t0", MassFraction2 + "_t0", MassFraction3 + "_t0" };
                case 5: return new string[] { MassFraction0 + "_t0", MassFraction1 + "_t0", MassFraction2 + "_t0", MassFraction3 + "_t0", MassFraction4 + "_t0" };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// the names of all mass fractions (at linearization point) listed in an array
        /// </summary>
        static public string[] MassFractions0(int NumberOfSpecies) {
            switch (NumberOfSpecies) {
                case 1: return new string[] { MassFraction0_0 };
                case 2: return new string[] { MassFraction0_0, MassFraction1_0 };
                case 3: return new string[] { MassFraction0_0, MassFraction1_0, MassFraction2_0 };
                case 4: return new string[] { MassFraction0_0, MassFraction1_0, MassFraction2_0, MassFraction3_0 };
                case 5: return new string[] { MassFraction0_0, MassFraction1_0, MassFraction2_0, MassFraction3_0, MassFraction4_0 };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// the names of all mean mass fractions listed in an array
        /// </summary>
        static public string[] MassFractionsMean(int NumberOfSpecies) {
            switch (NumberOfSpecies) {
                case 3: return new string[] { MassFraction0Mean, MassFraction1Mean, MassFraction2Mean };
                case 4: return new string[] { MassFraction0Mean, MassFraction1Mean, MassFraction2Mean, MassFraction3Mean };
                case 5: return new string[] { MassFraction0Mean, MassFraction1Mean, MassFraction2Mean, MassFraction3Mean , MassFraction4Mean };
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// the name of the n-th mass fraction
        /// </summary>
        static public string MassFraction_n(int n) {
            switch (n) {
                case 0: return MassFraction0;
                case 1: return MassFraction1;
                case 2: return MassFraction2;
                case 3: return MassFraction3;
                case 4: return MassFraction4;
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }

        /// <summary>
        /// the name of the n-th mass fraction
        /// </summary>
        static public string MassFraction_t0_n(int n) {
            switch (n) {
                case 0: return MassFraction0 + "_t0";
                case 1: return MassFraction1 + "_t0";
                case 2: return MassFraction2 + "_t0";
                case 3: return MassFraction3 + "_t0";
                case 4: return MassFraction4 + "_t0";
                default: throw new NotSupportedException("unsupported number of species.");
            }
        }


        /// <summary>
        /// The name of the <paramref name="d"/>-th velocity component.
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D, {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string Velocity_d(int d) {
            switch (d) {
                case 0: return VelocityX;
                case 1: return VelocityY;
                case 2: return VelocityZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }


        /// <summary>
        /// The name of the <paramref name="d"/>-th gravity (volume force) component.
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D, {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string Gravity_d(int d) {
            switch (d) {
                case 0: return GravityX;
                case 1: return GravityY;
                case 2: return GravityZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }
        
        /// <summary>
        /// the name of the <paramref name="d"/>-th convection component
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D, {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string Convective_d(int d) {
            switch (d) {
                case 0: return ConvectiveX;
                case 1: return ConvectiveY;
                case 2: return ConvectiveZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }
        
        /// <summary>
        /// the name of the <paramref name="d"/>-th vorticity component
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string Vorticity_d(int d) {
            switch (d) {
                case 0: return VorticityX;
                case 1: return VorticityY;
                case 2: return VorticityZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }
        
        /// <summary>
        /// the name of the <paramref name="d"/>-th intermediate velocity component (of a splitting scheme)
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D, {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string VelocityIntermed_d(int d) {
            switch (d) {
                case 0: return VelocityIntermedX;
                case 1: return VelocityIntermedY;
                case 2: return VelocityIntermedZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }
        
        /// <summary>
        /// the name of the <paramref name="d"/>-th velocity (linearization point) component
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D, {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string Velocity0_d(int d) {
            switch (d) {
                case 0: return Velocity0X;
                case 1: return Velocity0Y;
                case 2: return Velocity0Z;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }

        /// <summary>
        /// vector of velocity names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] VelocityVector(int D) {
            if (D == 2)
                return new string[] { VelocityX, VelocityY };
            else if (D == 3)
                return new string[] { VelocityX, VelocityY, VelocityZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }  

        /// <summary>
        /// vector of gravity/volume force names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] GravityVector(int D) {
            if (D == 2)
                return new string[] { GravityX, GravityY };
            else if (D == 3)
                return new string[] { GravityX, GravityY, GravityZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// component of gravity/volume force names
        /// </summary>
        public static string VolumeForce_d(int d) {
            return VolumeForceVector(3)[d];
        }


        /// <summary>
        /// vector of gravity/volume force names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] VolumeForceVector(int D) {
            if (D == 2)
                return new string[] { VolumeForceX, VolumeForceY };
            else if (D == 3)
                return new string[] { VolumeForceX, VolumeForceY, VolumeForceZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of convective names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] ConvectiveVector(int D) {
            if (D == 2)
                return new string[] { ConvectiveX, ConvectiveY };
            else if (D == 3)
                return new string[] { ConvectiveX, ConvectiveY, ConvectiveZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of velocity names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] VorticityVector(int D) {
            if (D == 2)
                return new string[] { Vorticity2D };
            else if (D == 3)
                return new string[] { VorticityX, VorticityY, VorticityZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of intermediate velocity names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] VelocityIntermedVector(int D) {
            if (D == 2)
                return new string[] { VelocityIntermedX, VelocityIntermedY };
            else if (D == 3)
                return new string[] { VelocityIntermedX, VelocityIntermedY, VelocityIntermedZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }
        
        /// <summary>
        /// vector of velocity names (linearization point)
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] Velocity0Vector(int D) {
            if (D == 2)
                return new string[] { Velocity0X, Velocity0Y };
            else if (D == 3)
                return new string[] { Velocity0X, Velocity0Y, Velocity0Z };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// the name of the <paramref name="d"/>-th boundary velocity component
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D, {0,1,2} in 3D;
        /// </param>
        /// <returns></returns>
        static public string BoundaryVelocity_d(int d) {
            switch (d) {
                case 0: return BoundaryVelocityX;
                case 1: return BoundaryVelocityY;
                case 2: return BoundaryVelocityZ;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }
        
        /// <summary>
        /// vector of boundary velocity names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] BoundaryVelocityVector(int D) {
            if (D == 2)
                return new string[] { BoundaryVelocityX, BoundaryVelocityY };
            else if (D == 3)
                return new string[] { BoundaryVelocityX, BoundaryVelocityY, BoundaryVelocityZ };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of mean velocity (linearization point) names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] Velocity0MeanVector(int D) {
            if (D == 2)
                return new string[] { Velocity0X_Mean, Velocity0Y_Mean };
            else if (D == 3)
                return new string[] { Velocity0X_Mean, Velocity0Y_Mean, Velocity0Z_Mean };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of mean momentum (linearization point) names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] Momentum0MeanVector(int D) {
            if (D == 2)
                return new string[] { Momentum0X_Mean, Momentum0Y_Mean };
            else if (D == 3)
                return new string[] { Momentum0X_Mean, Momentum0Y_Mean, Momentum0Z_Mean };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of velocity gradient names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[,] Velocity_GradientVector(int D) {
            if (D == 2)
                return new string[,] { { VelocityX_GradientX, VelocityX_GradientY }, { VelocityY_GradientX, VelocityY_GradientY } };
            else if (D == 3)
                return new string[,] { { VelocityX_GradientX, VelocityX_GradientY, VelocityX_GradientZ }, { VelocityY_GradientX, VelocityY_GradientY, VelocityY_GradientZ }, { VelocityZ_GradientX, VelocityZ_GradientY, VelocityZ_GradientZ } };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of velocity gradient names of velocity in x-direction
        /// </summary>
        public static string[] VelocityX_GradientVector() {
            return new string[] { VelocityX_GradientX, VelocityX_GradientY };
        }

        /// <summary>
        /// vector of velocity gradient names of velocity in y-direction
        /// </summary>
        public static string[] VelocityY_GradientVector() {
            return new string[] { VelocityY_GradientX, VelocityY_GradientY };
        }

        /// <summary>
        /// XX Tensor component of the extra stress tensor
        /// </summary>
        public const string StressXX = "StressXX";

        /// <summary>
        /// YY Tensor component of the extra stress tensor
        /// </summary>
        public const string StressYY = "StressYY";

        /// <summary>
        /// ZZ Tensor component of the extra stress tensor
        /// </summary>
        public const string StressZZ = "StressZZ";

        /// <summary>
        /// YZ Tensor component of the extra stress tensor
        /// </summary>
        public const string StressYZ = "StressYZ";

        /// <summary>
        /// XZ Tensor component of the extra stress tensor
        /// </summary>
        public const string StressXZ = "StressXZ";

        /// <summary>
        /// XY Tensor component of the extra stress tensor
        /// </summary>
        public const string StressXY = "StressXY";

        /// <summary>
        /// XX Tensor component of the extra stress tensor as parameter
        /// </summary>
        public const string StressXXP = "StressXXP";

        /// <summary>
        /// YY Tensor component of the extra stress tensor as parameter
        /// </summary>
        public const string StressYYP = "StressYYP";

        /// <summary>
        /// ZZ Tensor component of the extra stress tensor as parameter
        /// </summary>
        public const string StressZZP = "StressZZP";

        /// <summary>
        /// YZ Tensor component of the extra stress tensor as parameter
        /// </summary>
        public const string StressYZP = "StressYZP";

        /// <summary>
        /// XZ Tensor component of the extra stress tensor as parameter
        /// </summary>
        public const string StressXZP = "StressXZP";

        /// <summary>
        /// XY Tensor component of the extra stress tensor as parameter
        /// </summary>
        public const string StressXYP = "StressXYP";

        /// <summary>
        /// Extra stress tensor
        /// </summary>
        /// <param name="row">
        /// row index
        /// </param>
        /// /// <param name="column">
        /// column index
        /// </param>
        static public string ExtraStress(int row, int column)
        {
            switch (row)
            {
                case 0:
                    switch (column)
                    {
                        case 0: return StressXX;
                        case 1: return StressXY;
                        default: throw new NotSupportedException("unsupported spatial dimension");
                    }
                case 1:
                    switch (column)
                    {
                        case 0: return StressXY;
                        case 1: return StressYY;
                        default: throw new NotSupportedException("unsupported spatial dimension");
                    }
            }
            throw new NotSupportedException("unsupported spatial dimension");
        }

        /// <summary>
        /// Extra stress tensor in 3D
        /// </summary>
        /// <param name="row">
        /// row index
        /// </param>
        /// /// <param name="column">
        /// column index
        /// </param>
        static public string ExtraStress3D(int row, int column)
        {
            switch (row)
            {
                case 0:
                    switch (column)
                    {
                        case 0: return StressXX;
                        case 1: return StressXY;
                        case 2: return StressXZ;
                        default: throw new NotSupportedException("unsupported spatial dimension");
                    }
                case 1:
                    switch (column)
                    {
                        case 0: return StressXY;
                        case 1: return StressYY;
                        case 2: return StressYZ;
                        default: throw new NotSupportedException("unsupported spatial dimension");
                    }
                case 2:
                    switch (column)
                    {
                        case 0: return StressXZ;
                        case 1: return StressYZ;
                        case 2: return StressZZ;
                        default: throw new NotSupportedException("unsupported spatial dimension");
                    }
            }
            throw new NotSupportedException("unsupported spatial dimension");
        }
        /// <summary>
        /// Velocity gradient tensor XX component
        /// </summary>
        public const string VelocityXGradientX = "VelocityXGradientX";

        /// <summary>
        /// Velocity gradient tensor XY component
        /// </summary>
        public const string VelocityXGradientY = "VelocityXGradientY";

        /// <summary>
        /// Velocity gradient tensor YX component
        /// </summary>
        public const string VelocityYGradientX = "VelocityYGradientX";

        /// <summary>
        /// Velocity gradient tensor YY component
        /// </summary>
        public const string VelocityYGradientY = "VelocityYGradientY";


        /// <summary>
        /// Helical 
        /// </summary>
        public const string u = "ur";

        /// <summary>
        /// Helical
        /// </summary>
        public const string v = "uxi";

        /// <summary>
        /// Helical: call dierkes
        /// </summary>
        public const string w = "ueta";

        /// <summary>
        /// Max Sigma
        /// </summary>
        public const string MaxSigma = "MaximalSigma";

        /// <summary>
        /// vector of orientation of a rigid object, used in XNSERO. Currently only in 2D
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] OrientationVector(int D) {
            if (D == 2)
                return new string[] { OrientationVectorX, OrientationVectorY };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// the name of the <paramref name="d"/>-th orientationVector component
        /// </summary>
        /// <param name="d">
        /// spatial component {0,1} in 2D;
        /// </param>
        /// <returns></returns>
        static public string OrientationVectorComponent(int d) {
            switch (d) {
                case 0: return OrientationVectorX;
                case 1: return OrientationVectorY;
                default: throw new NotSupportedException("unsupported spatial dimension: d = " + d + ".");
            }
        }

        /// <summary>
        /// vector of orientation of a rigid object, used in XNSERO - x component
        /// </summary>
        public const string OrientationVectorX = "OrientationVectorX";

        /// <summary>
        /// vector of orientation of a rigid object, used in XNSERO - y component
        /// </summary>
        public const string OrientationVectorY = "OrientationVectorY";
    }
}
