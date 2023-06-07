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
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;

namespace BoSSS.Solution.NSECommon.Operator.Continuity {
    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class DivergenceAtIB : ILevelSetForm, ISupportsJacobianComponent {
        public DivergenceAtIB(int _D, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter, double _vorZeichen = 1.0) {
            this.D = _D;
            this.LevelSetIndex = iLevSet;
            this.PositiveSpecies = SolidSpecies;
            this.NegativeSpecies = FluidSpc;
            this.m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;
            this.scale = _vorZeichen;
        }

        int D;
        bool m_UseLevelSetVelocityParameter;
        // this scale is important! E.g. the XNSE multiplies the whole continuity equation by -1. this scale has to be set accordingly.
        double scale = 1.0;

        /// <summary>
        /// the penalty flux
        /// </summary>
        static double DirichletFlux(double UxN_in, double UxN_out) {
            return (UxN_in - UxN_out);
        }

        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double uAxN = GenericBlas.InnerProd(U_Neg, cp.Normal);

            double uBxN;
            if (m_UseLevelSetVelocityParameter) {
                uBxN = cp.Normal.InnerProd(cp.Parameters_IN);
            } else {
                uBxN = 0.0;
            }

            //if(cp.Parameters_IN.L2Norm() > 1.0e-10)
            //    throw new Exception("nonzero vel at ib");

            double FlxNeg = -DirichletFlux(uAxN, uBxN); // flux on A-side
            return scale * FlxNeg * v_Neg;
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(this.D);
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if (m_UseLevelSetVelocityParameter)
                    return VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(LevelSetIndex), VariableNames.VelocityVector(D));
                else
                    return null;
            }
        }

        public int LevelSetIndex {
            get;
            private set;
        }

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }

    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class DivergenceAtIB_LowMach : ILevelSetForm, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;

        public DivergenceAtIB_LowMach(int _D, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter, double _vorZeichen = 1.0) {
            this.D = _D;
            //this.m_LsTrk = lsTrk;
            this.LevelSetIndex = iLevSet;
            this.PositiveSpecies = SolidSpecies;
            this.NegativeSpecies = FluidSpc;
            this.m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;
            this.scale = _vorZeichen;
        }

        int D;
        bool m_UseLevelSetVelocityParameter;
        // this scale is important! E.g. the XNSE multiplies the whole continuity equation by -1. this scale has to be set accordingly.
        double scale = 1.0;

        /// <summary>
        /// the penalty flux
        /// </summary>
        static double DirichletFlux(double UxN_in, double UxN_out) {
            return (UxN_in - UxN_out);
        }

        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double uAxN = GenericBlas.InnerProd(U_Neg, cp.Normal);

            double uBxN;
            if (m_UseLevelSetVelocityParameter) {
                uBxN = cp.Normal.InnerProd(cp.Parameters_IN);
            } else {
                uBxN = 0.0;
            }

            //if(cp.Parameters_IN.L2Norm() > 1.0e-10)
            //    throw new Exception("nonzero vel at ib");

            double FlxNeg = -DirichletFlux(uAxN, uBxN); // flux on A-side
            return scale * FlxNeg * v_Neg;
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(this.D);
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if (m_UseLevelSetVelocityParameter)
                    return VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(LevelSetIndex), VariableNames.VelocityVector(D));
                else
                    return null;
            }
        }

        public int LevelSetIndex {
            get;
            private set;
        }

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }




}
