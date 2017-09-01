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

using MPI.Wrappers;

namespace ilPSP.LinSolvers.HYPRE {

    /// <summary>
    /// object-oriented wrapper around HYPRE_IJVector objects
    /// </summary>
    public class IJVector : IDisposable {

        /// <summary>
        /// pointer/handle to the HYPRE_IJVector object
        /// </summary>
        internal Wrappers.T_IJVector m_IJVector;

        /// <summary>
        /// pointer/handle to the HYPRE_Vector object which is used by the HYPRE_IJVector object
        /// </summary>
        internal Wrappers.T_ParCRS_vector ParCRS_vector;

        /// <summary>
        /// constructs and initializes a new HYPRE_IJVector object
        /// </summary>
        /// <param name="partition">distribution of the vector over MPI processes</param>
        public IJVector(IPartitioning partition) {
            if (partition.IsMutable)
                throw new NotSupportedException();

            if (partition.TotalLength > (int.MaxValue - 2))
                throw new ApplicationException("unable to create HYPRE vector: no. of matrix rows is larger than HYPRE index type (32 bit signed int);");

            m_VectorPartition = partition;

            int jLower = (int) partition.i0;
            int jUpper = (int)partition.i0 + partition.LocalLength - 1;
            MPI_Comm comm = csMPI.Raw._COMM.WORLD;
            int Nupdate = partition.LocalLength;

            // create object
            HypreException.Check(Wrappers.IJVector.Create(comm, jLower, jUpper, out m_IJVector));
            HypreException.Check(Wrappers.IJVector.SetObjectType(m_IJVector, Wrappers.Constants.HYPRE_PARCSR));
            HypreException.Check( Wrappers.IJVector.Initialize(m_IJVector));

            // set values
            int nvalues = Math.Min(1024, Nupdate);
            int[] indices = new int[nvalues];
            double[] values = new double[nvalues];
            int i0 = (int)m_VectorPartition.i0;
            for (int i = 0; i < Nupdate; i += nvalues) {

                if (i + nvalues > Nupdate) {
                    nvalues = Nupdate - i;
                }

                for (int ii = 0; ii < nvalues; ii++) {
                    indices[ii] = i + ii + i0;
                    //if (vec == null)
                    //    values[ii] = mapping[i_ii];
                    //else
                    //    values[ii] = vec[i_ii];
                }

                HypreException.Check(Wrappers.IJVector.SetValues(m_IJVector, nvalues, indices, values));
            }

            // assable
            HypreException.Check(Wrappers.IJVector.Assemble(m_IJVector));
            HypreException.Check(Wrappers.IJVector.GetObject(m_IJVector, out ParCRS_vector));
        }

        /// <summary>
        /// calls <see cref="Dispose"/>;
        /// </summary>
        ~IJVector() {
            Dispose();
        }
        
        IPartitioning m_VectorPartition;

        /// <summary>
        /// partition of this vector (very informative, isn't it? 
        /// Name and Type of Property allready explain everything)
        /// </summary>
        public IPartitioning VectorPartition { get { return m_VectorPartition; } }

        /// <summary>
        /// set the vector values
        /// </summary>
        /// <typeparam name="vectype"></typeparam>
        /// <param name="vec"></param>
        public void SetValues<vectype>(vectype vec) 
            where vectype: IList<double> {

            int Nupdate = m_VectorPartition.LocalLength;
            int i0 = (int)m_VectorPartition.i0;

            int nvalues = Math.Min(1024, Nupdate);
            int[] indices = new int[nvalues];
            double[] values = new double[nvalues];
            for (int i = 0; i < Nupdate; i += nvalues) {

                if (i + nvalues > Nupdate) {
                    nvalues = Nupdate - i;
                }

                for (int ii = 0; ii < nvalues; ii++) {
                    indices[ii] = i + ii + i0;
                    values[ii] = vec[i + ii];
                }
                
                HypreException.Check(Wrappers.IJVector.SetValues(m_IJVector, nvalues, indices, values));
            }
        }

        /// <summary>
        /// copies 
        /// </summary>
        /// <typeparam name="vectype"></typeparam>
        /// <param name="vec"></param>
        public void GetValues<vectype>(vectype vec)
            where vectype : IList<double> {
            int Nupdate = m_VectorPartition.LocalLength;
            int i0 = (int)m_VectorPartition.i0;
            
            int nvalues = Math.Min(1024, Nupdate);
            int[] indices = new int[nvalues];
            double[] values = new double[nvalues];
            for (int i = 0; i < Nupdate; i += nvalues) {

                if (i + nvalues > Nupdate) {
                    nvalues = Nupdate - i;
                }

                for (int ii = 0; ii < nvalues; ii++) {
                    indices[ii] = i + ii + i0;
                }

                HypreException.Check(Wrappers.IJVector.GetValues(m_IJVector, nvalues, indices, values));

                for (int ii = 0; ii < nvalues; ii++) {
                    vec[i+ii] = values[ii];
                }
            }

        }

        /// <summary>
        /// destroy the hypre objects, if not allready done
        /// </summary>
        public void Dispose() {
            if (m_IJVector.p != IntPtr.Zero) {
                HypreException.Check(Wrappers.IJVector.Destroy(m_IJVector));
                m_IJVector.p = IntPtr.Zero;
                ParCRS_vector.p = IntPtr.Zero;
            }
        }
    }
}
