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
using System.Reflection;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using BoSSS.Solution;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Base class for variable sets.
    /// </summary>
    public abstract class BaseVariableSet {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="GridDat"></param>
        /// <param name="Control"></param>
        /// <param name="IOFields"></param>
        /// <param name="RegisteredFields"></param>
        public BaseVariableSet(IGridData GridDat, SIMPLEControl Control, ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields) {
            InitFromAttributes.CreateFieldsAuto(this, GridDat, Control.FieldOptions, BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic, IOFields, RegisteredFields);
        }

        /// <summary>
        /// checks if any field contains NaN's or Inf's and throws an exception if this is the case.
        /// </summary>
        public bool CheckForNanOrInf() {
            int bool_FoundNanOrInf_local = 0;
            int bool_FoundNanOrInf_global = 0;

            Type myT = this.GetType();
            FieldInfo[] flds = myT.GetFields();

            string ErrMsg = "";
            bool Err = false;

            foreach (FieldInfo fi in flds) {

                if (fi.FieldType == typeof(SinglePhaseField)) {
                    DGField f = (DGField)fi.GetValue(this);
                    if (f != null) {
                        try {
                            if (f.CheckForNanOrInf(true, true, false) >= 0)
                                bool_FoundNanOrInf_local = 1;
                        } catch (ArithmeticException ae) {
                            ErrMsg += ae.Message + "\n";
                            Err = true;
                        }
                    }
                }

                if (fi.FieldType == typeof(VectorField<SinglePhaseField>)) {
                    VectorField<SinglePhaseField> vf = (VectorField<SinglePhaseField>)fi.GetValue(this);
                    if (vf != null) {
                        try {
                            if (vf.CheckForNanOrInf(true, true, false) >= 0)
                                bool_FoundNanOrInf_local = 1;
                        } catch (ArithmeticException ae) {
                            ErrMsg += ae.Message + "\n";
                            Err = true;
                        }
                    }
                }
            }

            if (Err)
                throw new ArithmeticException(ErrMsg);

            unsafe {
                csMPI.Raw.Allreduce((IntPtr)(&bool_FoundNanOrInf_local), (IntPtr)(&bool_FoundNanOrInf_global), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
            }

            bool FoundNanOrInf_global = (bool_FoundNanOrInf_global > 0) ? true : false;

            return FoundNanOrInf_global;
        }
    }
}
