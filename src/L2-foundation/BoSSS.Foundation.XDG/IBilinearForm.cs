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
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.XDG {
    
   

    /// <summary>
    /// this interface should be implemented by bulk equation components which require to switch coefficients based on species.
    /// </summary>
    public interface IEquationComponentSpeciesNotification : IEquationComponent {

        /// <summary>
        /// called before the integration on respective species 
        /// </summary>
        void SetParameter(string speciesName, SpeciesId SpcId);
    }

    /// <summary>
    /// this interface should be implemented by bulk equation components which are only valid in one species
    /// </summary>
    public interface ISpeciesFilter : IEquationComponent {


        /// <summary>
        /// the species in which the bulk equation component is valid;
        /// Null deactivates the Filter, i.e. the component is integrated for all species.
        /// </summary>
        string ValidSpecies { get; }
    }

    /// <summary>
    /// Interface for equation components which require e.g. grid and/or problem-dependent coefficients,
    /// e.g. cell length scales
    /// </summary>
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
    /// An integrand on the level set; per definition a level-set is always an interior edge
    /// </summary>
    public interface ILevelSetForm : IInnerEdgeForm {

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

        /// <summary>
        /// Controls integration at the level-set.
        /// </summary>
        TermActivationFlags LevelSetTerms { get; }
    }

    /// <summary>
    /// provides access to the level set tracker
    /// </summary>
    public interface ILevelSetFormSetup { 
        /// <summary>
        /// Called before Integration
        /// </summary>
        void Setup(LevelSetTracker lsTrk);

    }

    /// <summary>
    /// The XDG-counterpart of <see cref="BoSSS.Foundation.IEdgeForm_UxV"/>
    /// </summary>
    public interface ILevelSetForm_UxV : ILevelSetForm, IInnerEdgeform_UxV {
       
        //void LevelSetForm_UxV(LevSetIntParams inp, MultidimensionalArray Koeff_UxV);
    }

    /// <summary>
    /// The XDG-counterpart of <see cref="BoSSS.Foundation.IEdgeForm_UxV"/>
    /// </summary>
    public interface ILevelSetForm_GradUxV : ILevelSetForm, IInnerEdgeform_GradUxV {
        //void LevelSetForm_GradUxV(LevSetIntParams inp, MultidimensionalArray Koeff_GradUxV);

    }

    public interface ILevelSetForm_UxGradV : ILevelSetForm, IInnerEdgeform_UxGradV {

    }

    public interface ILevelSetForm_GradUxGradV : ILevelSetForm, IInnerEdgeform_GradUxGradV {

    }
    public interface ILevelSetForm_V : ILevelSetForm, IInnerEdgeSource_V {

    }

    public interface ILevelSetForm_GradV : ILevelSetForm, IInnerEdgeSource_GradV {

    }


    public interface INonlinLevelSetForm_V : ILevelSetForm, INonlinInnerEdgeForm_V {

    }

    public interface INonlinLevelSetForm_GradV : ILevelSetForm, INonlinInnerEdgeForm_GradV {

    }







    
}
