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


namespace BoSSS.Application.XNSE_Solver {

    
    /// <summary>
    /// Approximation \f$ \mathcal A^\vartheta\f$ 
    /// of predictor matrix \f$ \mathbf{A}^\vartheta_{C} - \frac{1}{Re} \mathbf{A}_{D}\f$  in correction step of the SIMPLE-algorithm.
    /// </summary>        
    public enum ApproxPredictor {

        /// <summary>
        /// Have a guess...
        /// </summary>
        MassMatrix,

        /// <summary>
        /// No approximation of the predictor matrix is performed, the predictor matrix itself is used;
        /// thus, its inverse is a dense matrix.
        /// Therefore, this option may only be usefull for small testcases. 
        /// If used without under-relaxation the SIMPLE algorithm 
        /// must converge, in the Stokes case, in one step, since it is mathematically equivalent to a saddle-point solver.
        /// </summary>
        Exact, 
        

        /// <summary>
        /// Diagonal matrix of predictor.
        /// </summary>
        Diagonal,


        /// <summary>
        /// Block diagonal matrix of predictor.
        /// </summary>
        BlockDiagonal,

        /// <summary>
        /// "Row-Sum" of the Blocks.
        /// </summary>
        BlockSum,


        /// <summary>
        /// Neumann-Series of the Matrix
        /// If the Matrix has the property
        /// \lim\limits_{n \rightarrow \infty}{(I-A)^n} = 0
        /// the Inverse of the Matrix can be approximated as 
        /// A^{-1} = \sum_{n = 0}^\infty (I - A)^n
        /// </summary>
        Neumann,


        /// <summary>
        /// Calculates the Operator Without the Fluxes over Inner edges
        /// </summary>
        LocalizedOperator

        

    }
}
