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
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Foundation.XDG {
    
    /// <summary>
    /// method parameters (aka. arguments) for level set integrands
    /// </summary>
	public sealed class LevSetIntParams {
		
        /// <summary>
		/// integration nodes in physical coordinates
		/// </summary>
		/// <remarks>
		/// <list type="bullet">
		///   <item>1st index: cell</item>
		///   <item>2nd index: quadrature node</item>
		///   <item>3rd index: spatial dimension</item>
		/// </list>
		/// </remarks>			
		public MultidimensionalArray X;

		/// <summary>
		/// normals on level set at quadrature points
		/// </summary>
		/// <list type="bullet">
		///   <item>1st index: cell</item>
		///   <item>2nd index: quadrature node</item>
		///   <item>3rd index: spatial direction</item>
		/// </list>
		public MultidimensionalArray Normal;

		/// <summary>
        /// Values of parameter fields at quadrature nodes on positive side, i.e. where the level-set function is positive.
		/// </summary>
		public MultidimensionalArray[] ParamsPos;

        /// <summary>
        /// Values of parameter fields at quadrature nodes on negative side, i.e. where the level-set function is negative.
        /// </summary>
        public MultidimensionalArray[] ParamsNeg;
      

        /// <summary>
        /// 1st item to integrate, i.e. cell/edge offset;
        /// </summary>
        public int i0;

        /// <summary>
        /// number of items
        /// </summary>
        public int Len {
            get {
                int L = Normal.GetLength(0);
                return L;
            }
        }

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        /// <summary>
        /// Access to the level-set tracker.
        /// </summary>
        public LevelSetTracker LsTrk;
	}

    /// <summary>
    /// common input parameters for the abstract functions
    /// </summary>
    public struct CommonParamsLs {

        /// <summary>
        /// Position vector, in global/physical coordinates at which the integrand should be evaluated.
        /// Note: depending on the quadrature rule, that point is not necessarily located on the zero-level-set.
        /// </summary>
        public double[] x;

        /// <summary>
        /// Normal vector, parallel to the gradient of the level-set function.
        /// </summary>
        public double[] n;

        /// <summary>
        /// Guess what?
        /// </summary>
        public int SpatialDim {
            get {
                Debug.Assert(x.Length == n.Length);
                return n.Length;
            }
        }

        /// <summary>
        /// Values of parameter variables on negative side, i.e. where the level-set function is negative.
        /// </summary>
        public double[] ParamsNeg;

        /// <summary>
        /// Values of parameter variables on positive side, i.e. where the level-set function is positive.
        /// </summary>
        public double[] ParamsPos;

        /// <summary>
        /// cell index
        /// </summary>
        public int jCell;

        /// <summary>
        /// Physical time.
        /// </summary>
        public double time;

        ///// <summary>
        ///// A characteristic length scale for the negative part of the cut cell, i.e. the respective part of the cut cell where the level-set field is negative.
        ///// </summary>
        //public double NegCellLengthScale;

        ///// <summary>
        ///// A characteristic length scale for the positive part of the cut cell, i.e. the respective part of the cut cell where the level-set field is positive.
        ///// </summary>
        //public double PosCellLengthScale;
    }

    /// <summary>
    /// this interface should be implemented by bulk equation components which require to switch coefficients based on species.
    /// </summary>
    public interface IEquationComponentSpeciesNotification {
        void SetParameter(string speciesName, SpeciesId SpcId);
    }

    /// <summary>
    /// this interface should be implemented by bulk equation components which are only valid in one species
    /// </summary>
    public interface ISpeciesFilter {

        /// <summary>
        /// the species in which the bulk equation component is valid
        /// </summary>
        SpeciesId validSpecies { get; }
    }

    /// <summary>
    /// An integrand on the level set.
    /// </summary>
    public interface ILevelSetForm : IEquationComponent {


        /// <summary>
        /// index of the species-separating level set.
        /// </summary>
        int LevelSetIndex { get; }


        /// <summary>
        /// regarding this integrand, the species on the positive side of level set number <see cref="LevelSetIndex"/>
        /// </summary>
        SpeciesId PositiveSpecies { get; }

        /// <summary>
        /// guess what?
        /// </summary>
        SpeciesId NegativeSpecies { get; }


        TermActivationFlags LevelSetTerms { get; }

        double LevelSetForm(ref CommonParamsLs inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB);
    }

    /// <summary>
    /// Interface for equation components which require e.g. grid and/or problem-dependent coefficients,
    /// e.g. cell length scales
    /// </summary>
    /// <seealso cref="XSpatialOperator.XEvaluatorBase.OperatorCoefficients"/>
    public interface ILevelSetEquationComponentCoefficient : IEquationComponent {

        /// <summary>
        /// Passes various coefficients to the equation component.
        /// </summary>
        /// <param name="csA">
        /// Coefficient set related to negative species (<see cref="NegativeSpecies"/>)
        /// </param>
        /// <param name="csB">
        /// Coefficient set related to positive species (<see cref="PositiveSpecies"/>)
        /// </param>
        /// <param name="DomainDGdeg">
        /// actual polynomial degree of domain variables
        /// </param>
        /// <param name="TestDGdeg">
        /// actual polynomial degree of codomain variables
        /// </param>
        void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg);

    }

    /// <summary>
    /// The XDG-counterpart of <see cref="BoSSS.Foundation.IEdgeform_UxV"/>
    /// </summary>
    public interface ILevelSetForm_UxV : ILevelSetForm {
       
        void LevelSetForm_UxV(LevSetIntParams inp, MultidimensionalArray Koeff_UxV);
    }

    /// <summary>
    /// The XDG-counterpart of <see cref="BoSSS.Foundation.IEdgeform_UxV"/>
    /// </summary>
    public interface ILevelSetForm_GradUxV : ILevelSetForm {
        void LevelSetForm_GradUxV(LevSetIntParams inp, MultidimensionalArray Koeff_GradUxV);

    }

    public interface ILevelSetForm_UxGradV : ILevelSetForm {
        void LevelSetForm_UxGradV(LevSetIntParams inp, MultidimensionalArray Koeff_UxGradV);

    }

    public interface ILevelSetForm_GradUxGradV : ILevelSetForm {
        void LevelSetForm_GradUxGradV(LevSetIntParams inp, MultidimensionalArray Koeff_GradUxGradV);

    }

    public interface ILevelSetForm_V : ILevelSetForm {

        void LevelSetForm_V(LevSetIntParams inp, 
            MultidimensionalArray Koeff_V);

    }

    public interface ILevelSetForm_GradV : ILevelSetForm {

        void LevelSetForm_GradV(LevSetIntParams inp, 
            MultidimensionalArray Koeff_GradV);

    }


    public interface INonlinLevelSetForm_V : ILevelSetForm {

        void LevelSetForm_V(LevSetIntParams inp, 
            MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB,
            MultidimensionalArray Koeff_V);

    }

    public interface INonlinLevelSetForm_GradV : ILevelSetForm {

        void LevelSetForm_GradV(LevSetIntParams inp, 
            MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB,
            MultidimensionalArray Koeff_GradV);

    }







    
}
