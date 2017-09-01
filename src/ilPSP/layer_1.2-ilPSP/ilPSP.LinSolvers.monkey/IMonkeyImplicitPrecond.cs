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

namespace ilPSP.LinSolvers.monkey {

    /// <summary>
    /// basic definition of an implicit preconditioner in monkey;
    /// </summary>
    public interface IMonkeyImplicitPrecond : ilPSP.LinSolvers.IImplicitPrecond {

        /// <summary>
        /// computes, if neccessary, data for preconditioning (such as incomplete LU decomposition)
        /// </summary>
        /// <param name="OrigMatrix">
        /// original matrix of the linear system
        /// </param>
        /// <param name="mt">
        /// matrix type will be specified by thew solver that hosts the preconditioner;
        /// If the preconditioner needs to create his own matrices, it shouch use this storage format;
        /// However, this is only a 'recommendation'.
        /// </param>
        /// <param name="dev">
        /// factory for creating matrices and vector objects
        /// </param>
        void Initialize(IMutableMatrixEx OrigMatrix, Device dev, MatrixType mt);
        
        /// <summary>
        /// 1st pass of implicit preconditioning preconditioning: 
        /// </summary>
        /// <param name="pc_output">
        /// output of preconditioning, must be locked (see <see cref="LockAbleObject.IsLocked"/>);
        /// on exit of the <see cref="DoPrecond"/>-method, this vector contains the 
        /// result of the preconditioning operation.
        /// </param>
        /// <param name="pc_input">
        /// input for preconditioning, must be locked (see <see cref="LockAbleObject.IsLocked"/>);
        /// before the nesting solver calls <see cref="DoPrecond"/> to perform the precond., it 
        /// places the input vector in this object.
        /// </param>
        /// <param name="mtx">
        /// original matrix of the solver, in locked state
        /// </param>
        /// <param name="dev">
        /// factory  to create temporary vectors and matrices
        /// </param>
        void CreateTempObjects(VectorBase pc_output, VectorBase pc_input, MatrixBase mtx, Device dev);
        
        /// <summary>
        /// cleanup; 
        /// </summary>
        void ReleaseTempObjects();
        
        /// <summary>
        /// in every iteration of the solver, this method is called to perform the 
        /// preconditionig, i.e. to compute z = M*b, where z, M and b are Vatrices and 
        /// vectors specified in <see cref="CreateTempObjects"/>
        /// </summary>
        void DoPrecond();
    }
}
