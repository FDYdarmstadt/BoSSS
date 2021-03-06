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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Solution.XNSECommon.Operator.Pressure {
    
    /// <summary>
    /// 
    /// </summary>
    public class PressureFormAtLevelSet : ILevelSetForm, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;

        public PressureFormAtLevelSet(int _d, int _D, 
            double _wA = 1.0, double _wB = 1.0, 
            bool _freeSurface = false, double _pFree = 0.0) {

            m_d = _d;
            m_D = _D;
            //m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentException();

            wA = _wA;
            wB = _wB;

            freeSurface = _freeSurface;
            pFree = _pFree;
        }

        int m_d;
        int m_D;

        double wA;
        double wB;

        bool freeSurface;
        double pFree;

        public double InnerEdgeForm(ref CommonParams inp, double[] pA, double[] pB, double[,] Grad_pA, double[,] Grad_pB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            if (!freeSurface) {
                return -(vB - vA) * inp.Normal[m_d] * 0.5 * (pB[0] + pA[0]);
            } else {
                //return -(pFree * vB - pA[0] * vA) * inp.Normal[m_d];
                return -(vB - vA) * inp.Normal[m_d] * 0.5 * (pB[0] + pA[0]);
            }
            
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Pressure };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }

     
    }


}
