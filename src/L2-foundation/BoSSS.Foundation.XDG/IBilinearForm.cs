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
        /// A characteristic length scale for the negative parts of the cut cells, i.e. the respective part of a cut cell where the level-set field is negative.
        /// </summary>
        public MultidimensionalArray NegCellLengthScale;


        /// <summary>
        /// A characteristic length scale for the positive parts of the cut cells, i.e. the respective part of a cut cell where the level-set field is positive.
        /// </summary>
        public MultidimensionalArray PosCellLengthScale;
        

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

        /// <summary>
        /// A characteristic length scale for the negative part of the cut cell, i.e. the respective part of the cut cell where the level-set field is negative.
        /// </summary>
        public double NegCellLengthScale;

        /// <summary>
        /// A characteristic length scale for the positive part of the cut cell, i.e. the respective part of the cut cell where the level-set field is positive.
        /// </summary>
        public double PosCellLengthScale;
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


    public interface ILevelSetForm_UxV : ILevelSetForm {
       
        void LevelSetForm_UxV(LevSetIntParams inp, MultidimensionalArray Koeff_UxV);
    }

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

        void LevelSetForm_V(LevSetIntParams inp, MultidimensionalArray Koeff_V);

    }

    public interface ILevelSetForm_GradV : ILevelSetForm {

        void LevelSetForm_GradV(LevSetIntParams inp, MultidimensionalArray Koeff_GradV);

    }







    /*

    /// <summary>
    /// a bilinear form on the level set;
    /// </summary>
    public interface IBilinearForm : ILevelSetIntegrand {



        /// <summary>
        /// 
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="Koeff_UxV">\f$ M_{u,v}\f$ </param>
        /// <param name="Koeff_NablaUxV"></param>
        /// <param name="Koeff_UxNablaV"></param>
        /// <param name="Koeff_NablaUxNablaV"></param>
        void EdgeForm(LevSetIntParams inp, MultidimensionalArray Koeff_UxV, MultidimensionalArray Koeff_NablaUxV, MultidimensionalArray Koeff_UxNablaV, MultidimensionalArray Koeff_NablaUxNablaV);

    }
    */

    /*
    public interface ILinear2ndDerivativeCouplingFlux : ILevelSetIntegrand {


        /// <summary>
        /// </summary>
        /// <param name="FunctionMatrix"></param>
        /// <param name="AffineOffset"></param>
        /// <param name="inparams">given as reference for performance reasons, 
        /// DO NOT WRITE to this structure;</param>
        void PrimalVar_LevelSetFlux(LevSetIntParams inparams,
                                    MultidimensionalArray FunctionMatrix,
                                    MultidimensionalArray AffineOffset);



        /// <summary>
        /// 
        /// </summary>
        /// <param name="PrimalVar_FunctionMatrix"></param>
        /// <param name="DerivVar_FunctionMatrix">
        /// 
        /// </param>
        /// <param name="AffineOffset"></param>
        /// <param name="inparams">given as reference for performance reasons, 
        /// DO NOT WRITE to this structure;</param>
        void DerivativVar_LevelSetFlux(LevSetIntParams inparams,
                                       MultidimensionalArray PrimalVar_FunctionMatrix, MultidimensionalArray DerivVar_FunctionMatrix,
                                       MultidimensionalArray AffineOffset);


    }
    */
}
