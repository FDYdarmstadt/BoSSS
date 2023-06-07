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
using System.Diagnostics;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP.Utils;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Implementation of the time derivative as a linearized source term in the low-Mach combustion solver.
    /// Based on the implicit Euler scheme.
    /// </summary>
    public class MassMatrixComponent : BoSSS.Solution.Utils.LinearSource {
        string[] m_ArgumentOrdering;
        IList<string> m_ParameterOrdering;
        MaterialLaw EoS;
        int m_SpatialDimension = 2;
        int NumberOfReactants = 3;
        int j;
        PhysicsMode physicsMode;
        /// <summary>
        /// Ctor for variable density flows
        /// </summary> 
        /// <param name="EoS">The material law</param>
        /// <param name="energy">Set conti: true for the energy equation</param>
        /// <param name="ArgumentOrdering"></param>
        /// <param name="TimeStepSize"></param>
        public MassMatrixComponent(MaterialLaw EoS, String[] ArgumentOrdering, int j, PhysicsMode _physicsMode, int spatDim, int NumberOfReactants) {
            this.j = j;
            this.EoS = EoS;
            m_ParameterOrdering = EoS.ParameterOrdering;
            this.m_SpatialDimension = spatDim;
            this.NumberOfReactants = NumberOfReactants;

            int SpatDim = 2;
            int numberOfReactants = 3;
            this.physicsMode = _physicsMode;
            switch (_physicsMode) {
                case PhysicsMode.Multiphase:
                m_ArgumentOrdering = new string[] { VariableNames.LevelSet };
                m_ParameterOrdering = ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDim), VariableNames.Velocity0MeanVector(SpatDim));
                break;
                case PhysicsMode.LowMach:
                m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.Temperature);
                if (EoS == null)
                    throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                else
                    this.EoS = EoS;
                break;
                case PhysicsMode.MixtureFraction:
                m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.MixtureFraction);

                if (EoS == null)
                    throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                else
                    this.EoS = EoS;
                break;
                case PhysicsMode.Combustion:
                m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(SpatDim), VariableNames.Temperature, VariableNames.MassFractions(numberOfReactants - 1)); // u,v,w,T, Y0,Y1,Y2,Y3  as variables (Y4 is calculated as Y4 = 1- (Y0+Y1+Y2+Y3)
                if (EoS == null)
                    throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
                else
                    this.EoS = EoS;
                break;
                default:
                throw new NotImplementedException();
            }
        }
        /// <summary>
        /// ctor for constant density flows
        /// </summary>
        /// <param name="EoS"></param>
        /// <param name="TimeStepSize"></param>
        /// <param name="ArgumentOrdering"></param>
        /// <param name="energy"></param>
        /// <param name="conti"></param>
        public MassMatrixComponent(String[] ArgumentOrdering) {
            m_ArgumentOrdering = ArgumentOrdering;
            this.EoS = null;
            m_ParameterOrdering = null;
        }



        /// <summary>
        /// The argument
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }

        /// <summary>
        /// Paramaters used to compute the density
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
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
            double mult = 1.0;
            double rho = 1.0;

            switch (physicsMode) {
                case PhysicsMode.Incompressible:
                break;
                case PhysicsMode.MixtureFraction:
                rho = EoS.getDensityFromZ(U[m_SpatialDimension]);
                break;
                case PhysicsMode.LowMach:
                double[] DensityArgumentsIn = U.GetSubVector(m_SpatialDimension, 1);
                rho = EoS.GetDensity(DensityArgumentsIn);
                break;
                case PhysicsMode.Combustion:
                double[] DensityArgumentsIn2 = U.GetSubVector(m_SpatialDimension, NumberOfReactants);
                rho = EoS.GetDensity(DensityArgumentsIn2);
                break;
                default:
                throw new NotImplementedException("PhysicsMode not implemented");
            }

            return mult * rho * U[j];


        }
    }






    /// <summary>
    /// Implementation of the time derivative as a linearized source term in the low-Mach combustion solver.
    /// Based on the implicit Euler scheme.
    /// </summary>
    public class MassMatrixLowMachComponent : IVolumeForm, ISupportsJacobianComponent {
        IList<string> m_ArgumentOrdering;
        IList<string> m_ParameterOrdering;
        MaterialLaw EoS;
        int m_SpatialDimension;
        int NumberOfReactants;
        int j;



        /// <summary>
        /// Ctor for variable density flows based on the mixture fraction approach
        /// </summary> 
        /// <param name="EoS">The material law</param>
        /// <param name="energy">Set conti: true for the energy equation</param>
        /// <param name="varname"></param>
        /// <param name="TimeStepSize"></param>
        public MassMatrixLowMachComponent(MaterialLaw EoS, string varname, int spatDim) {
            this.EoS = EoS;
            this.m_SpatialDimension = spatDim;
            m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(spatDim), VariableNames.Pressure,VariableNames.MixtureFraction);
            this.j = m_ArgumentOrdering.IndexOf(varname);

            if (varname == VariableNames.Pressure && j != 2)
                throw new Exception("!!!");
     

            if (EoS == null)
                throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
            else
                this.EoS = EoS;
        }

        /// <summary>
        /// Ctor for variable density flows
        /// </summary> 
        /// <param name="EoS">The material law</param>
        /// <param name="energy">Set conti: true for the energy equation</param>
        /// <param name="varname"></param>
        /// <param name="TimeStepSize"></param>
        public MassMatrixLowMachComponent(MaterialLaw EoS, string varname, int spatDim, int NumberOfReactants) {
            this.EoS = EoS;
            this.m_SpatialDimension = spatDim;
            this.NumberOfReactants = NumberOfReactants;
            m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(spatDim), VariableNames.Pressure, VariableNames.Temperature, VariableNames.MassFractions(NumberOfReactants));
            this.j = m_ArgumentOrdering.IndexOf(varname);

            if (varname == VariableNames.Pressure && j != 2)
                throw new Exception("!!!");
            if (j < 0 || j > NumberOfReactants + 3)
                throw new ArgumentOutOfRangeException();

            if (EoS == null)
                throw new ApplicationException("EoS has to be given for Low-Mach flows to calculate density.");
            else
                this.EoS = EoS;
        }

        public IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }

        public IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DerivVol };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double rho;
            if(EoS is MaterialLawMixtureFractionNew) {
                //double[] DensityArgumentsIn = U.GetSubVector(m_SpatialDimension + 1,1);
                rho = EoS.getDensityFromZ(U[3]);              

            } else {
                double[] DensityArgumentsIn = U.GetSubVector(m_SpatialDimension + 1, NumberOfReactants + 1);
                rho = EoS.GetDensity(DensityArgumentsIn);
            }
            
            double ret = rho * U[j] * V;

            if (j == 2)
                ret = 0.0;

            //if (j == 2)
            //    ret = rho * V ;

            return ret;
        }

        /// <summary>
        /// Active terms are <see cref="TermActivationFlags.UxV"/> and
        /// <see cref="TermActivationFlags.V"/>
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.UxV;
            }
        }
    }
}
