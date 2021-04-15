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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation;
using ilPSP;

namespace BoSSS.Solution.NSECommon.Operator.Convection {
    public class ConvectionAtIB : ILevelSetForm {
        public ConvectionAtIB(int _d, int _D, LevelSetTracker LsTrk, double _LFFA, IncompressibleBoundaryCondMap _bcmap, double fluidDensity, bool UseMovingMesh, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter){
            m_LsTrk = LsTrk;
            m_D = _D;
            m_d = _d;
            LFFA = _LFFA;
            //varMode = _varMode;
            fDensity = fluidDensity;
            m_UseMovingMesh = UseMovingMesh;
            m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;

            NegFlux = new LinearizedConvection(_D, _bcmap, _d);
            //NegFlux = new ConvectionInBulk_LLF(_D, _bcmap, _d, fluidDensity, 0, _LFFA, double.NaN, LsTrk);
            //NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"), null);
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
        }
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        bool m_UseLevelSetVelocityParameter;

        LevelSetTracker m_LsTrk;
        int m_D;
        int m_d;

        double fDensity;
        double LFFA;
        bool m_UseMovingMesh;

        // Use Fluxes as in Bulk Convection
        LinearizedConvection NegFlux;


        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Velocity_d(m_d) }; }
        }

        public IList<string> ParameterOrdering {
            get {
                var ret = ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));
                if(m_UseLevelSetVelocityParameter)
                    ret = ret.Cat(VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.VelocityVector(m_D)));
                return ret;
            }
        }
     
        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        /// <summary>
        /// Species ID of the solid
        /// </summary>
        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId(m_SolidSpecies); }
        }

        /// <summary>
        /// Species ID of the fluid; 
        /// </summary>
        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId(m_FluidSpc); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

       

        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {




            if(m_UseMovingMesh == true) {
                return 0.0;
            } else {
                BoSSS.Foundation.CommonParams inp = cp;

                /*
                double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
                double wLevSet = parameters_P[2];
                pRadius = parameters_P[3];

                double[] uLevSet_temp = new double[1];
                if(m_d == 0) {
                    uLevSet_temp[0] = uLevSet[0] + pRadius * wLevSet * -cp.Normal[1];
                } else {
                    uLevSet_temp[0] = uLevSet[1] + pRadius * wLevSet * cp.Normal[0];
                }

                //Outer values for Velocity and VelocityMean
                inp.Parameters_OUT[0] = uLevSet[0] + pRadius * wLevSet * -cp.Normal[1];
                inp.Parameters_OUT[1] = uLevSet[1] + pRadius * wLevSet * cp.Normal[0];
                // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
                inp.Parameters_OUT[2] = 0;
                inp.Parameters_OUT[3] = 0;
                */

                double[] uLevSet;
                if(m_UseLevelSetVelocityParameter) {
                    uLevSet = inp.Parameters_IN.GetSubVector(m_D * 2, m_D);
                } else {
                    uLevSet = new double[m_D];
                }

                double[] uLevSet_temp = new double[1];
                uLevSet_temp[0] = uLevSet[m_d];
                
                inp.Parameters_OUT = new double[inp.D * 2];
                Array.Copy(uLevSet, 0, inp.Parameters_OUT, 0, m_D);
                Array.Copy(uLevSet, 0, inp.Parameters_OUT, m_D, m_D);

                double FlxNeg = this.NegFlux.InnerEdgeForm(ref inp, U_Neg, uLevSet_temp, null, null, v_Neg, 0, null, null);
                if(FlxNeg.IsNaNorInf())
                    throw new ArithmeticException("NaN/Inf in immersed boundary convection");

                return FlxNeg*fDensity;
            }
        }
    
    
        
    }

    public class ConvectionAtIB_Newton : ILevelSetForm, ISupportsJacobianComponent {
        public ConvectionAtIB_Newton(int _d, int _D, LevelSetTracker LsTrk, double _LFFA, IncompressibleBoundaryCondMap _bcmap, double fluidDensity, bool UseMovingMesh, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter) {
            m_LsTrk = LsTrk;
            m_D = _D;
            m_d = _d;
            LFFA = _LFFA;
            //varMode = _varMode;
            fDensity = fluidDensity;
            m_UseMovingMesh = UseMovingMesh;
            m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;

            NegFlux = new LinearizedConvectionJacobi(_D, _bcmap, _d);
            //NegFlux = new ConvectionInBulk_LLF(_D, _bcmap, _d, fluidDensity, 0, _LFFA, double.NaN, LsTrk);
            //NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"), null);
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
        }
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        bool m_UseLevelSetVelocityParameter;

        LevelSetTracker m_LsTrk;
        int m_D;
        int m_d;

        double fDensity;
        double LFFA;
        bool m_UseMovingMesh;

        // Use Fluxes as in Bulk Convection
        LinearizedConvectionJacobi NegFlux;


        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(m_D); }
        }

        public IList<string> ParameterOrdering {
            get {
                var ret = new string[] { };
                if (m_UseLevelSetVelocityParameter)
                    ret = ret.Cat(VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.VelocityVector(m_D)));
                return ret;
            }
        }

        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        /// <summary>
        /// Species ID of the solid
        /// </summary>
        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId(m_SolidSpecies); }
        }

        /// <summary>
        /// Species ID of the fluid; 
        /// </summary>
        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId(m_FluidSpc); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }



        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {




            if (m_UseMovingMesh == true) {
                return 0.0;
            } else {
                BoSSS.Foundation.CommonParams inp = cp;

                /*
                double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
                double wLevSet = parameters_P[2];
                pRadius = parameters_P[3];

                double[] uLevSet_temp = new double[1];
                if(m_d == 0) {
                    uLevSet_temp[0] = uLevSet[0] + pRadius * wLevSet * -cp.Normal[1];
                } else {
                    uLevSet_temp[0] = uLevSet[1] + pRadius * wLevSet * cp.Normal[0];
                }

                //Outer values for Velocity and VelocityMean
                inp.Parameters_OUT[0] = uLevSet[0] + pRadius * wLevSet * -cp.Normal[1];
                inp.Parameters_OUT[1] = uLevSet[1] + pRadius * wLevSet * cp.Normal[0];
                // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
                inp.Parameters_OUT[2] = 0;
                inp.Parameters_OUT[3] = 0;
                */

                double[] uLevSet;
                if (m_UseLevelSetVelocityParameter) {
                    uLevSet = inp.Parameters_IN.GetSubVector(0, m_D);
                } else {
                    uLevSet = new double[m_D];
                }

                double[] uLevSet_temp = uLevSet;

                //inp.Parameters_OUT = new double[inp.D * 2];
                //Array.Copy(uLevSet, 0, inp.Parameters_OUT, 0, m_D);
                //Array.Copy(uLevSet, 0, inp.Parameters_OUT, m_D, m_D);

                double FlxNeg = this.NegFlux.InnerEdgeForm(ref inp, U_Neg, uLevSet_temp, null, null, v_Neg, 0, null, null);
                if (FlxNeg.IsNaNorInf())
                    throw new ArithmeticException("NaN/Inf in immersed boundary convection");

                return FlxNeg * fDensity;
            }
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }
}
