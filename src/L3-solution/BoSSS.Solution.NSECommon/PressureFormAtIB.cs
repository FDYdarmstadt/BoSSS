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
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Solution.NSECommon.Operator.Pressure { 
    
    
    public class PressureFormAtIB : ILevelSetForm, ISupportsJacobianComponent {
        public PressureFormAtIB(int _d, int _D, LevelSetTracker LsTrk, int iLevSet, string FluidSpc, string SolidSpecies) {
            m_d = _d;
            m_D = _D;
            m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentException();
        
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
        }
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;

        LevelSetTracker m_LsTrk;
        int m_d;
        int m_D;

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Pressure };
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
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] pA, double[] pB, double[,] Grad_pA, double[,] Grad_pB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            return vA * inp.Normal[m_d] * pA[0];
        }


        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }
}
