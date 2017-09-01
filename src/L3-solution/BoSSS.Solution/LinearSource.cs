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

using System.Collections.Generic;
/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using BoSSS.Foundation;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Utility class which helps the user in creating the function matrices that are needed
    /// in an implementation of <see cref="IVolumeForm"/>; The user has to override
    /// <see cref="Source(double[],double[],double[])"/>, where he is able to implement the linear function as 
    /// an algebraic formula. All function matrices and offsets (aka. intercept) are constructed from the user-defined
    /// functions by this class.
    /// </summary>
    public abstract class LinearSource : IVolumeForm {

        /// <summary>
        /// not in use, returning null
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        /// <summary>
        /// to be implemented by user 
        /// </summary>
        abstract public IList<string> ArgumentOrdering {
            get;
        }


        /// <summary>
        /// override this method to implement the linear source
        /// </summary>
        abstract protected double Source(double[] x, double[] parameters, double[] U);



        /// <summary>
        /// arguments used for the linear source functions to extract the function matrix
        /// </summary>
        double[] m_Arguments;


        /// <summary>
        /// helper function to initialize <see cref="m_Arguments"/>
        /// </summary>
        private void AllocateArrays() {
            if (m_Arguments == null) {
                m_Arguments = new double[ArgumentOrdering.Count];
            }
        }

        /// <summary>
        /// Active terms are <see cref="TermActivationFlags.UxV"/> and
        /// <see cref="TermActivationFlags.V"/>
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        /// translates <see cref="Source"/> into <see cref="IVolumeForm.VolumeForm"/>
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return this.Source(cpv.Xglobal, cpv.Parameters, U) * V;
        }
    }
}
