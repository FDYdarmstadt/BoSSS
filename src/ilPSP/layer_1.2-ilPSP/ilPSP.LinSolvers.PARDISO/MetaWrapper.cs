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

namespace ilPSP.LinSolvers.PARDISO {
    
    
    /// <summary>
    /// Another wrapper layer that encapsulates 
    /// </summary>
    class MetaWrapper {

        /// <summary>
        /// ctor.
        /// </summary>
        public MetaWrapper(Version __V) {
            Init(__V);
        }

        private void Init(Version __V) {
            switch (__V) {
                case Version.MKL: if(mkl == null) mkl = new Wrapper_MKL(); break;
                case Version.v4: if(v4 == null) v4 = new Wrapper_v4(); break;
                case Version.v5: if(v5 == null) v5 = new Wrapper_v5(); break;
            }
        }


        static Wrapper_MKL mkl;

        static Wrapper_v4 v4;

        static Wrapper_v5 v5;


        /// <summary>
        /// PARDISO interface, see PARDISO documentation
        /// </summary>
        public unsafe int PARDISOINIT(void* pt, int* mtype, int* iparm, double* dparam) {
            if (mkl != null) {
                return mkl.PARDISOINIT(pt, mtype, iparm);
            } else if (v4 != null) {
                int error = 0;
                int solver = 0; // use sparse direct
                v4.PARDISOINIT(pt, mtype, &solver, iparm, dparam,&error);
                return error;
            } else if (v5 != null) {
                int error = 0;
                int solver = 0; // use sparse direct
                v5.PARDISOINIT(pt, mtype, &solver, iparm, dparam, &error);
                return error;
            }
            throw new NotImplementedException();
        }

        /// <summary>
        /// PARDISO interface, see PARDISO documentation
        /// </summary>
        public unsafe int PARDISO(void* pt, int* maxfct, int* mnum, int* mtype,
                                            int* phase, int* n,
                                            double* a, int* ia, int* ja,
                                            int* perm, int* nrhs, int* iparm,
                                            int* msglvl, double* b, double* x,
                                            int* error, double* dparam) {
            //
            if (mkl != null) {
                return mkl.PARDISO(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error);
            } else if (v4 != null) {
                v4.PARDISO(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error, dparam);
                return *error;
            } else if (v5 != null) {
                v5.PARDISO(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error, dparam);
                return *error;
            }
            throw new NotImplementedException();
        }


        public string PARDISOerror2string(int error) {
            if (mkl != null) {
                return mkl.PARDISOerror2string(error);
            } else if (v4 != null) {
                return v4.PARDISOerror2string(error);
            } else if (v5 != null) {
                return v5.PARDISOerror2string(error);
            }
            throw new NotImplementedException();
        }
    }
}
