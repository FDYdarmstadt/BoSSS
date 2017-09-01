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
using System.Runtime.InteropServices;

namespace BoSSS.Platform {

#pragma warning disable 1591
    /// <summary>
    /// Some parts of the GNU scientific library which are used by certain
    /// parts of BoSSS
    /// </summary>
    public class GSL {

        private const string LIBNAME = "libgsl-0.dll";

        // Vectors
        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_vector_alloc")]
        public static extern IntPtr gsl_vector_alloc(int n);

        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_vector_set")]
        public static extern void gsl_vector_set(IntPtr v, int i, double x);

        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_vector_free")]
        public static extern void gsl_vector_free(IntPtr v);


        // Matrices
        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_matrix_alloc")]
        public static extern IntPtr gsl_matrix_alloc(int m, int n);

        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_matrix_set")]
        public static extern void gsl_matrix_set(IntPtr v, int i, int j, double x);

        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_matrix_free")]
        public static extern void gsl_matrix_free(IntPtr m);


        // Polynomial roots
        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_poly_complex_workspace_alloc")]
        public static extern IntPtr gsl_poly_complex_workspace_alloc(int noOfCoefficients);

        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_poly_complex_solve")]
        public static extern void gsl_poly_complex_solve(IntPtr coefficients, int noOfCoefficients, IntPtr workSpace, IntPtr result);

        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_poly_complex_workspace_free")]
        public static extern void gsl_poly_complex_workspace_free(IntPtr workSpace);


        // Linear least squares
        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_linalg_SV_decomp")]
        public static extern int gsl_linalg_SV_decompc(IntPtr A, IntPtr V, IntPtr S, IntPtr work);

        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_linalg_SV_solve")]
        public static extern int gsl_linalg_SV_solve(IntPtr U, IntPtr V, IntPtr S, IntPtr b, IntPtr x);

        // Dilagorithmic function
        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_sf_dilog")]
        public static extern double gsl_sf_dilog(double x);

        /// <summary>
        /// Jacobian elliptic functions sn(u|m), cn(u|m), dn(u|m) (see
        /// https://www.gnu.org/software/gsl/manual/html_node/Elliptic-Functions-_0028Jacobi_0029.html#Elliptic-Functions-_0028Jacobi_0029)
        /// </summary>
        /// <param name="u">
        /// Evaluation point
        /// </param>
        /// <param name="m">
        /// Modulus of the elliptic function. Note that the definition of this
        /// modulus is different from the Maple definition, for example, which
        /// uses sqrt(m) instead.
        /// </param>
        /// <param name="sn">
        /// Sinus of the Jacobi amplitude function
        /// JacobiAM(<paramref name="u"/>, <paramref name="m"/>)
        /// </param>
        /// <param name="cn">
        /// Co-sinus of the Jacobi amplitude function
        /// JacobiAM(<paramref name="u"/>, <paramref name="m"/>)
        /// </param>
        /// <param name="dn">
        /// Derivative of the Jacobi amplitude function
        /// JacobiAM(<paramref name="u"/>, <paramref name="m"/>)
        /// </param>
        /// <returns></returns>
        [DllImportAttribute(LIBNAME, EntryPoint = "gsl_sf_elljac_e")]
        public static extern int gsl_sf_elljac_e(double u, double m, double[] sn, double[] cn, double[] dn);
    }
#pragma warning restore 1591
}

