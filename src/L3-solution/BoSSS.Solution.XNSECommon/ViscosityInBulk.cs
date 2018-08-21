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
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using System.Collections;
using ilPSP;
using System.Diagnostics;
using BoSSS.Foundation;

namespace BoSSS.Solution.XNSECommon.Operator.Viscosity {
    
    
    public class ViscosityInBulk_GradUTerm : BoSSS.Solution.NSECommon.swipViscosity_Term1, IEquationComponentSpeciesNotification {

        public ViscosityInBulk_GradUTerm(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, int d, int D, double _muA, double _muB,
            double _betaA = 0.0, double _betaB = 0.0)
            : base(penalty, d, D, bcMap, NSECommon.ViscosityOption.ConstantViscosity, constantViscosityValue: double.NegativeInfinity) {
            muA = _muA;
            muB = _muB;
            //betaA = _betaA;
            //betaB = _betaB;
            base.m_alpha = sw;
            this.m_bcMap = bcMap;
            base.velFunction = null;
            this.m_penalty = penalty;
        }

        double muA;
        double muB;

        double betaA;
        double betaB;

        double currentMu = double.NaN;
        double complementMu = double.NaN;

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

        /// <summary>
        /// multiplier for the penalty computation
        /// </summary>
        double m_penalty;

        ///// <summary>
        ///// length scales for cells in order to determine the penalty parameter.
        ///// </summary>
        //MultidimensionalArray m_LenScales;

        public void SetParameter(string speciesName, SpeciesId SpcId) {
            switch (speciesName) {
                case "A": currentMu = muA; complementMu = muB; SetBndfunction("A"); m_beta = betaA;  break;
                case "B": currentMu = muB; complementMu = muA; SetBndfunction("B"); m_beta = betaB;  break;
                default: throw new ArgumentException("Unknown species.");
            }

            //double muFactor; // the WTF factor
            //if (jCellOut >= 0)
            //    muFactor = 1.0;
            //else
            //    muFactor = Math.Max(currentMu, complementMu) / currentMu;
            double muFactor = Math.Max(currentMu, complementMu) / currentMu;
            base.m_penalty_base = this.m_penalty * muFactor;
        }

        void SetBndfunction(string S) {
            int D = base.m_D;
            base.velFunction = D.ForLoop(d => this.m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + S]);
        }

        
        //protected override double penalty(int jCellIn, int jCellOut) {
            
        //    double muFactor; // the WTF factor
        //    if(jCellOut >= 0)
        //        muFactor = 1.0; 
        //    else
        //        muFactor = Math.Max(currentMu, complementMu) / currentMu;
        //    double penaltySizeFactor_A = 1.0 / this.m_LenScales[jCellIn];
        //    double penaltySizeFactor_B = jCellOut >= 0 ?  1.0 / this.m_LenScales[jCellOut] : 0;
        //    Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
        //    Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
        //    Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
        //    Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
        //    double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);
        //    return this.m_penalty * penaltySizeFactor * muFactor;  // Bem.: die Viskosität wird in der swipViscosity_Term1.InnerEdgeForm(...) dazumultipliziert.
            
        //    //return (base.penalty(jCellIn, jCellOut, cj)/currentMu)*Math.Max(currentMu, complementMu);
        //}
         

        protected override double Viscosity(double[] Parameters) {
            return currentMu;
        }

        //public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {

        //    m_slipLengths = (MultidimensionalArray)cs.UserDefinedValues["SlipLengths"];
        //}

    }

    public class ViscosityInBulk_GradUtranspTerm : BoSSS.Solution.NSECommon.swipViscosity_Term2, IEquationComponentSpeciesNotification {

        public ViscosityInBulk_GradUtranspTerm(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, int d, int D, double _muA, double _muB, 
            double _betaA = 0.0, double _betaB = 0.0)
            : base(penalty, d, D, bcMap, NSECommon.ViscosityOption.ConstantViscosity, constantViscosityValue: double.NegativeInfinity) {
            muA = _muA;
            muB = _muB;
            //betaA = _betaA;
            //betaB = _betaB;
            base.m_alpha = sw;
            base.velFunction = null;
            this.m_bcMap = bcMap;
            this.m_penalty = penalty;
        }

        double muA;
        double muB;

        double betaA;
        double betaB;

        double currentMu = double.NaN;
        double complementMu = double.NaN;
        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

        /// <summary>
        /// multiplier for the penalty computation
        /// </summary>
        double m_penalty;

        ///// <summary>
        ///// length scales for cells in order to determine the penalty parameter.
        ///// </summary>
        //MultidimensionalArray m_LenScales;


        public void SetParameter(string speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": currentMu = muA; complementMu = muB; SetBndfunction("A"); m_beta = betaA;  break;
                case "B": currentMu = muB; complementMu = muA; SetBndfunction("B"); m_beta = betaB; break;
                default: throw new ArgumentException("Unknown species.");
            }

            

            //double muFactor; // the WTF factor
            //if (jCellOut >= 0)
            //    muFactor = 1.0;
            //else
            //    muFactor = Math.Max(currentMu, complementMu) / currentMu;
            double muFactor = Math.Max(currentMu, complementMu) / currentMu;
            base.m_penalty_base = this.m_penalty * muFactor;
        }

        void SetBndfunction(string S) {
            int D = base.m_D;
            base.velFunction = D.ForLoop(d => this.m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + S]);
        }


        //protected override double penalty(int jCellIn, int jCellOut) {
            
        //    double muFactor; // the WTF factor
        //    if(jCellOut >= 0)
        //        muFactor = 1.0;
        //    else
        //        muFactor = Math.Max(currentMu, complementMu) / currentMu;
        //    double penaltySizeFactor_A = 1.0 / this.m_LenScales[jCellIn];
        //    double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / this.m_LenScales[jCellOut] : 0;
        //    Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
        //    Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
        //    Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
        //    Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
        //    double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);
        //    return this.m_penalty * penaltySizeFactor * muFactor;  // Bem.: die Viskosität wird in der swipViscosity_Term2.InnerEdgeForm(...) dazumultipliziert.
            
        //    //return (base.penalty(jCellIn, jCellOut, cj)/currentMu)*Math.Max(currentMu, complementMu);
        //}
        

        protected override double Viscosity(double[] Parameters) {
            return currentMu;
        }

        //public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
        //    m_LenScales = cs.CellLengthScales;
        //}
    }

    public class ViscosityInBulk_divTerm : BoSSS.Solution.NSECommon.swipViscosity_Term3, IEquationComponentSpeciesNotification {

        public ViscosityInBulk_divTerm(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, int d, int D, double _muA, double _muB)
            : base(penalty, d, D, bcMap, NSECommon.ViscosityOption.ConstantViscosity, constantViscosityValue: double.NegativeInfinity) {
            muA = _muA;
            muB = _muB;
            base.m_alpha = sw;
            base.velFunction = null;
            this.m_bcMap = bcMap;
            this.m_penalty = penalty;
        }

        double muA;
        double muB;
        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

        double currentMu = double.NaN;
        double complementMu = double.NaN;

        /// <summary>
        /// multiplier for the penalty computation
        /// </summary>
        double m_penalty;

        ///// <summary>
        ///// length scales for cells in order to determine the penalty parameter.
        ///// </summary>
        //MultidimensionalArray m_LenScales;

        public void SetParameter(string speciesName, SpeciesId SpcIds) {
            switch (speciesName) {
                case "A": currentMu = muA; complementMu = muB; SetBndfunction("A");  break;
                case "B": currentMu = muB; complementMu = muA; SetBndfunction("B");  break;
                default: throw new ArgumentException("Unknown species.");
            }

            double muFactor = Math.Max(currentMu, complementMu) / currentMu;
            base.m_penalty_base = this.m_penalty * muFactor;
        }
               

        void SetBndfunction(string S) {
            int D = base.m_D;
            base.velFunction = D.ForLoop(d => this.m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + S]);
        }

        protected override double Viscosity(double[] Parameters) {
            return currentMu;
        }

        /*
        protected override double penalty(int jCellIn, int jCellOut) {
            double penaltySizeFactor_A = 1.0 / this.m_LenScales[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / this.m_LenScales[jCellOut] : 0;
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);
            return this.m_penalty * penaltySizeFactor; // Bem.: die Viskosität wird in der swipViscosity_Term3.InnerEdgeForm(...) dazumultipliziert.

            // return (base.penalty(jCellIn, jCellOut, cj)/currentMu)*Math.Max(currentMu, complementMu);
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            m_LenScales = cs.CellLengthScales;
        }
        */
    }
    
}
