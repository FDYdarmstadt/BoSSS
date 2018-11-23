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
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP.LinSolvers.MUMPS
{
	public unsafe class MUMPS_csharp
	{
        [DllImport("dmumps")]
        extern static unsafe int* mumps_get_mapping();

        [DllImport("dmumps")]
		extern static unsafe int* mumps_get_pivnul_list();

		[DllImport("dmumps")]
		extern static unsafe int* mumps_get_sym_perm();

		[DllImport("dmumps")]
		extern static unsafe int* mumps_get_uns_perm();

		[DllImport("dmumps", EntryPoint = "dmumps_f77_")]
		extern static unsafe void MUMPS_F77(
		   int* job,
		   int* sym,
		   int* par,
		   int* comm_fortran,
		   int* n,
		   int* icntl,
		   double* cntl,
		   int* keep,
		   double* dkeep,
		   long* keep8,
		   int* nz,
		   int* irn,
		   int* irn_avail,
		   int* jcn,
		   int* jcn_avail,
		   double* a,
		   int* a_avail,
		   int* nz_loc,
		   int* irn_loc,
		   int* irn_loc_avail,
		   int* jcn_loc,
		   int* jcn_loc_avail,
		   double* a_loc,
		   int* a_loc_avail,
		   int* nelt,
		   int* eltptr,
		   int* eltptr_avail,
		   int* eltvar,
		   int* eltvar_avail,
		   double* a_elt,
		   int* a_elt_avail,
		   int* perm_in,
		   int* perm_in_avail,
		   double* rhs,
		   int* rhs_avail,
		   double* redrhs,
		   int* redrhs_avail,
		   int* info,
		   double* rinfo,
		   int* infog,
		   double* rinfog,
		   int* deficiency,
		   int* lwk_user,
		   int* size_schur,
		   int* listvar_schur,
		   int* listvar_schur_avail,
		   double* schur,
		   int* schur_avail,
		   double* wk_user,
		   int* wk_user_avail,
		   double* colsca,
		   int* colsca_avail,
		   double* rowsca,
		   int* rowsca_avail,
		   int* instance_number,
		   int* nrhs,
		   int* lrhs,
		   int* lredrhs,
		   double* rhs_sparse,
		   int* rhs_sparse_avail,
		   double* sol_loc,
		   int* sol_loc_avail,
		   int* irhs_sparse,
		   int* irhs_sparse_avail,
		   int* irhs_ptr,
		   int* irhs_ptr_avail,
		   int* isol_loc,
		   int* isol_loc_avail,
		   int* nz_rhs,
		   int* lsol_loc,
		   int* schur_mloc,
		   int* schur_nloc,
		   int* schur_lld,
		   int* schur_mblock,
		   int* schur_nblock,
		   int* schur_nprow,
		   int* schur_npcol,
		   int* ooc_tmpdir,
		   int* ooc_prefix,
		   int* write_problem,
		   int* ooc_tmpdirlen,
		   int* ooc_prefixlen,
		   int* write_problemlen
		   );

		public unsafe struct DMUMPS_STRUC_CS {

			public int sym, par, job;
			public int comm_fortran;    /* Fortran communicator */
			//public fixed int icntl[40];
			public int[] icntl;
			public int[] keep;
			public double[] cntl;
			public double[] dkeep;
			public long[] keep8;
			public int n; // order of the matrix

			public int nz_alloc; /* used in matlab interface to decide if we
								free + malloc when we have large variation */

			/* Assembled entry */
			public int nz; // nonzeros in matrix
			public int[] irn; // array of row indices
			public int[] jcn; // array of col indices
			public double[] a; // array of values

			/* Distributed entry */
			public int nz_loc;
			public int[] irn_loc;
			public int[] jcn_loc;
			public double[] a_loc;

			/* Element entry */
			public int nelt;
			public int[] eltptr;
			public int[] eltvar;
			public double[] a_elt;

			/* Ordering, if given by user */
			public int[] perm_in;

			/* Orderings returned to user */
			public int[] sym_perm;    /* symmetric permutation */
			public int[] uns_perm;    /* column permutation */

			/* Scaling (inout but complicated) */
			public double[] colsca;
			public double[] rowsca;
			public int colsca_from_mumps;
			public int rowsca_from_mumps;

			/* RHS, solution, ouptput data and statistics */
			public double[] rhs;
			public double[] redrhs, rhs_sparse;
			public double[] sol_loc;
			public int[] irhs_sparse, irhs_ptr, isol_loc;
			public int nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
			public int schur_mloc, schur_nloc, schur_lld;
			public int mblock, nblock, nprow, npcol;
			public int[] info, infog;
			public double[] rinfo, rinfog;

			/* Null space */
			public int deficiency;
			public int[] pivnul_list; // Row indices corresponding to the null pivots
			public int[] mapping;

			/* Schur */
			//public int size_schur;
			public int[] listvar_schur;
			public double[] schur; // On output of factorization contains the Schurmatrix

			/* Internal parameters */
			public int instance_number;
			public double[] wk_user;

			/* Version number: length=14 in FORTRAN + 1 for final \0 + 1 for alignment */
			public string version_number;
			/* For out-of-core */
			public string ooc_tmpdir;
			public string ooc_prefix;
			/* To save the matrix in matrix market format */
			public string write_problem;
			public int lwk_user;

		};

		unsafe static int Extract_Pointers_d(double* component, double* dummyPointer) {
			if (component == null) {
				component = dummyPointer;
				return 0;
			} else {
				return 1;
			}
		}

		unsafe static int Extract_Pointers_i(int* component, int* dummyPointer) {
			if (component == null) {
				component = dummyPointer;
				return 0;
			} else {
				return 1;
			}
		}

		public static void mumps_cs(ref DMUMPS_STRUC_CS mumps_par) {
			unsafe
			{
				fixed (double* rhs = mumps_par.rhs, sol_loc = mumps_par.sol_loc, redrhs = mumps_par.redrhs, a = mumps_par.a, rhs_sparse = mumps_par.rhs_sparse, a_loc = mumps_par.a_loc, a_elt = mumps_par.a_elt,
					colsca = mumps_par.colsca, rowsca = mumps_par.rowsca, schur = mumps_par.schur, wk_user = mumps_par.wk_user, cntl = mumps_par.cntl, dkeep = mumps_par.dkeep, rinfo = mumps_par.rinfo, rinfog = mumps_par.rinfog)
				{

					fixed (int* listvar_schur = mumps_par.listvar_schur, ooc_tmpdir = new int[255], ooc_prefix = new int[63], write_problem = new int[255], irn = mumps_par.irn, jcn = mumps_par.jcn,
						irn_loc = mumps_par.irn_loc, jcn_loc = mumps_par.jcn_loc, eltptr = mumps_par.eltptr, eltvar = mumps_par.eltvar, perm_in = mumps_par.perm_in, sym_perm = mumps_par.sym_perm,
						uns_perm = mumps_par.uns_perm, irhs_sparse = mumps_par.irhs_sparse, irhs_ptr = mumps_par.irhs_ptr, isol_loc = mumps_par.isol_loc, pivnul_list = mumps_par.pivnul_list, mapping = mumps_par.mapping
						, icntl = &(mumps_par.icntl[0]), keep = mumps_par.keep, info = mumps_par.info, infog = mumps_par.infog,

						//test
						job = &mumps_par.job, sym = &mumps_par.sym, par = &mumps_par.par, comm_fortran = &mumps_par.comm_fortran, n = &mumps_par.n, nz = &mumps_par.nz, nz_loc = &mumps_par.nz_loc, nelt = &mumps_par.nelt,
						deficiency = &mumps_par.deficiency, lwk_user = &mumps_par.lwk_user, instance_number = &mumps_par.instance_number, nrhs = &mumps_par.nrhs, lrhs = &mumps_par.lrhs, lredrhs = &mumps_par.lredrhs,
						nz_rhs = &mumps_par.nz_rhs, lsol_loc = &mumps_par.lsol_loc, schur_mloc = &mumps_par.schur_mloc, schur_nloc = &mumps_par.schur_nloc, schur_lld = &mumps_par.schur_lld, mblock = &mumps_par.mblock, nblock = &mumps_par.nblock,
						nprow = &mumps_par.nprow, npcol = &mumps_par.npcol
						)
					{
						fixed (long* keep8 = mumps_par.keep8)
						{
							{
								//int instance_number; mumps_par.instance_number = &instance_number; 
								//int*[] icntl = new int*[40];
								int size_schur; // = mumps_par.listvar_schur.Length;

								/*
								* The following local variables will 
								*  be passed to the F77 interface.
								*/
								int perm_in_avail;
								int listvar_schur_avail;
								int schur_avail;
								int wk_user_avail;
								int irn_avail, jcn_avail, a_avail, rhs_avail, redrhs_avail;
								/* These are actually used
								 * as booleans, but we stick
								 * to simple types for the
								 * C-F77 interface */
								int irn_loc_avail, jcn_loc_avail, a_loc_avail;
								int eltptr_avail, eltvar_avail, a_elt_avail;
								int colsca_avail, rowsca_avail;
								int irhs_ptr_avail, rhs_sparse_avail, sol_loc_avail;
								int irhs_sparse_avail, isol_loc_avail;
								//int* info; int* infog;
								// MUMPS_REAL* rinfo; MUMPS_REAL* rinfog;
								/* Other local variables */
								//MUMPS_REAL rdummy; MUMPS_REAL* rdummyp;
								//MUMPS_COMPLEX cdummy; MUMPS_COMPLEX* cdummyp;
								/* String lengths to be passed to Fortran by address */
								int ooc_tmpdirlen;
								int ooc_prefixlen;
								int write_problemlen;
								int version_number;
								int i;
								const int no = 0;
								const int yes = 1;
								int idummy; int* idummyp;
								double rdummy; double* rdummyp;
								idummyp = &idummy;
								rdummyp = &rdummy;
								/* [SDCZ]MUMPS_F77 always calls either
								 * MUMPS_NULLIFY_C_COLSCA or MUMPS_ASSIGN_C_COLSCA
								 * (and ROWSCA). The next two lines are thus not
								 * strictly necessary. */
								//MUMPS_COLSCA_STATIC = 0;
								//MUMPS_ROWSCA_STATIC = 0;

								if (mumps_par.job == -1) { /* job = -1: we just reset all pointers to 0 */ // Automatically done in csharp
																										   //mumps_par.irn = 0; mumps_par.jcn = 0; mumps_par.a = 0; mumps_par.rhs = 0; mumps_par.wk_user = 0;
																										   //mumps_par.redrhs = 0;
																										   //mumps_par.eltptr = 0; mumps_par.eltvar = 0; mumps_par.a_elt = 0; mumps_par.perm_in = 0; mumps_par.sym_perm = 0; mumps_par.uns_perm = 0; mumps_par.irn_loc = 0; mumps_par.jcn_loc = 0; mumps_par.a_loc = 0; mumps_par.listvar_schur = 0; mumps_par.schur = 0; mumps_par.mapping = 0; mumps_par.pivnul_list = 0; mumps_par.colsca = 0; mumps_par.colsca_from_mumps = 0;
																										   //mumps_par.rowsca = 0; mumps_par.colsca_from_mumps = 0; mumps_par.rhs_sparse = 0; mumps_par.irhs_sparse = 0; mumps_par.sol_loc = 0; mumps_par.irhs_ptr = 0; mumps_par.isol_loc = 0;

									// Configuring strings
									mumps_par.ooc_tmpdir = "NAME_NOT_INITIALIZED";
									mumps_par.ooc_prefix = "NAME_NOT_INITIALIZED";
									mumps_par.write_problem = "NAME_NOT_INITIALIZED";
									mumps_par.version_number = "5.0.2\0";
									//mumps_par->version_number[MUMPS_VERSION_MAX_LEN + 1] = '\0';
									/* Next line initializes scalars to arbitrary values.
									 * Some of those will anyway be overwritten during the
									 * call to Fortran routine [SDCZ]MUMPS_INIT_PHASE */
									mumps_par.n = 0; mumps_par.nz = 0; mumps_par.nz_loc = 0; mumps_par.nelt = 0; mumps_par.deficiency = 0; mumps_par.lwk_user = 0; size_schur = 0; mumps_par.lrhs = 0; mumps_par.lredrhs = 0; mumps_par.nrhs = 0; mumps_par.nz_rhs = 0; mumps_par.lsol_loc = 0;
									mumps_par.schur_mloc = 0; mumps_par.schur_nloc = 0; mumps_par.schur_lld = 0; mumps_par.mblock = 0; mumps_par.nblock = 0; mumps_par.nprow = 0; mumps_par.npcol = 0;
								}
								mumps_par.ooc_tmpdir = "NAME_NOT_INITIALIZED";
								mumps_par.ooc_prefix = "NAME_NOT_INITIALIZED";
								mumps_par.write_problem = "NAME_NOT_INITIALIZED";
								mumps_par.version_number = "5.0.2\0";
								ooc_tmpdirlen = mumps_par.ooc_tmpdir.Length;
								ooc_prefixlen = mumps_par.ooc_prefix.Length;
								write_problemlen = mumps_par.write_problem.Length;
								/* Avoid the use of strnlen which may not be
								 * available on all systems. Allow strings without
								 * \0 at the end, if the file is not found, the
								 * Fortran layer is responsible for raising an
								 * error.  */
								if (ooc_tmpdirlen > 255) {
									ooc_tmpdirlen = 255;
								}
								if (ooc_prefixlen > 63) {
									ooc_prefixlen = 63;
								}
								if (write_problemlen > 255) {
									write_problemlen = 255;
								}

								irn_avail = Extract_Pointers_i(irn, idummyp);
								jcn_avail = Extract_Pointers_i(jcn, idummyp);
								rhs_avail = Extract_Pointers_d(rhs, rdummyp);
								wk_user_avail = Extract_Pointers_d(wk_user, rdummyp);
								redrhs_avail = Extract_Pointers_d(redrhs, rdummyp);
								irn_loc_avail = Extract_Pointers_i(irn_loc, idummyp);
								jcn_loc_avail = Extract_Pointers_i(jcn_loc, idummyp);
								a_loc_avail = Extract_Pointers_d(a_loc, rdummyp);
								a_avail = Extract_Pointers_d(a, rdummyp);
								eltptr_avail = Extract_Pointers_i(eltptr, idummyp);
								eltvar_avail = Extract_Pointers_i(eltvar, idummyp);
								a_elt_avail = Extract_Pointers_d(a_elt, rdummyp);
								perm_in_avail = Extract_Pointers_i(perm_in, idummyp);
								listvar_schur_avail = Extract_Pointers_i(listvar_schur, idummyp);
								schur_avail = Extract_Pointers_d(schur, rdummyp);

								/* EXTRACT_POINTERS not adapted to rowsca and colsca */
								if (mumps_par.rowsca != null && mumps_par.rowsca_from_mumps == 0) {
									/* has been set by user and was not allocated in mumps */
									//rowsca = mumps_par->rowsca;
									rowsca_avail = yes;
								} else {
									/* FIXME: changing rowsca in C after an earlier call
									   where rowsca was computed by mumps is not possible. */
									//rowsca = rdummyp;
									rowsca_avail = no;
								}
								if (mumps_par.colsca != null && mumps_par.colsca_from_mumps == 0)
								  /* has been changed by user and was not allocated in mumps */
								  {
									//colsca = mumps_par->colsca;
									colsca_avail = yes;
								} else {
									/* FIXME: changing colsca in C after an earlier call
									   where colsca was computed by mumps is not possible */
									//colsca = rdummyp;
									colsca_avail = no;
								}

								rhs_sparse_avail = Extract_Pointers_d(rhs_sparse, rdummyp);
								sol_loc_avail = Extract_Pointers_d(sol_loc, rdummyp);
								irhs_sparse_avail = Extract_Pointers_i(irhs_sparse, idummyp);
								isol_loc_avail = Extract_Pointers_i(isol_loc, idummyp);
								irhs_ptr_avail = Extract_Pointers_i(irhs_ptr, idummyp);

								//EXTRACT_POINTERS(rhs_sparse, cdummyp);
								//EXTRACT_POINTERS(sol_loc, cdummyp);
								//EXTRACT_POINTERS(irhs_sparse, idummyp);
								//EXTRACT_POINTERS(isol_loc, idummyp);
								//EXTRACT_POINTERS(irhs_ptr, idummyp);
								/* printf("irn_avail,jcn_avail, rhs_avail, a_avail, eltptr_avail, eltvar_avail,a_elt_avail,perm_in_avail= %d %d %d %d %d %d %d \n", irn_avail,jcn_avail, rhs_avail, a_avail, eltptr_avail, eltvar_avail, a_elt_avail, perm_in_avail); */
								/*
								 * Extract integers (input) or pointers that are
								 * always allocated (such as ICNTL, INFO, ...)
								 */
								/* size_schur = mumps_par->size_schur; */
								//int instance_number = mumps_par->instance_number;
								//icntl = mumps_par->icntl;
								//cntl = mumps_par->cntl;
								//keep = mumps_par->keep;
								//dkeep = mumps_par->dkeep;
								//keep8 = mumps_par->keep8;
								//info = mumps_par->info;
								//infog = mumps_par->infog;
								//rinfo = mumps_par->rinfo;
								//rinfog = mumps_par->rinfog;

								for (i = 0; i < ooc_tmpdirlen; i++) {
									ooc_tmpdir[i] = (int)mumps_par.ooc_tmpdir[i];
								}
								for (i = 0; i < ooc_prefixlen; i++) {
									ooc_prefix[i] = (int)mumps_par.ooc_prefix[i];
								}
								for (i = 0; i < write_problemlen; i++) {
									write_problem[i] = (int)mumps_par.write_problem[i];
								}
								for (i = 0; i < mumps_par.version_number.Length; i++) {
									version_number = (int)mumps_par.version_number[i];
								}


								/* Call F77 interface */
								//   MUMPS_F77(&mumps_par.job, &mumps_par.sym, &mumps_par.par, &mumps_par.comm_fortran, &mumps_par.n, icntl, cntl, keep, dkeep, keep8, &mumps_par.nz
								//, irn, &irn_avail, jcn, &jcn_avail, a, &a_avail, &mumps_par.nz_loc, irn_loc, &irn_loc_avail, jcn_loc, &jcn_loc_avail
								//, a_loc, &a_loc_avail, &mumps_par.nelt, eltptr, &eltptr_avail, eltvar, &eltvar_avail, a_elt, &a_elt_avail, perm_in, &perm_in_avail,
								//rhs, &rhs_avail, redrhs, &redrhs_avail, info, rinfo, infog, rinfog, &mumps_par.deficiency, &mumps_par.lwk_user, &size_schur, listvar_schur,
								//&listvar_schur_avail, schur, &schur_avail, wk_user, &wk_user_avail, colsca, &colsca_avail, rowsca, &rowsca_avail, &mumps_par.instance_number,
								//&mumps_par.nrhs, &mumps_par.lrhs, &mumps_par.lredrhs, rhs_sparse, &rhs_sparse_avail, sol_loc, &sol_loc_avail, irhs_sparse, &irhs_sparse_avail, irhs_ptr, &irhs_ptr_avail,
								//isol_loc, &isol_loc_avail, &mumps_par.nz_rhs, &mumps_par.lsol_loc, &mumps_par.schur_mloc, &mumps_par.schur_nloc, &mumps_par.schur_lld, &mumps_par.mblock, &mumps_par.nblock,
								//&mumps_par.nprow, &mumps_par.npcol, ooc_tmpdir, ooc_prefix, write_problem, &ooc_tmpdirlen, &ooc_prefixlen, &write_problemlen);

								MUMPS_F77(job, sym, par, comm_fortran, n, icntl, cntl, keep, dkeep, keep8, nz
						   , irn, &irn_avail, jcn, &jcn_avail, a, &a_avail, nz_loc, irn_loc, &irn_loc_avail, jcn_loc, &jcn_loc_avail
						   , a_loc, &a_loc_avail, nelt, eltptr, &eltptr_avail, eltvar, &eltvar_avail, a_elt, &a_elt_avail, perm_in, &perm_in_avail,
						   rhs, &rhs_avail, redrhs, &redrhs_avail, info, rinfo, infog, rinfog, deficiency, lwk_user, &size_schur, listvar_schur,
						   &listvar_schur_avail, schur, &schur_avail, wk_user, &wk_user_avail, colsca, &colsca_avail, rowsca, &rowsca_avail, instance_number,
						   nrhs, lrhs, lredrhs, rhs_sparse, &rhs_sparse_avail, sol_loc, &sol_loc_avail, irhs_sparse, &irhs_sparse_avail, irhs_ptr, &irhs_ptr_avail,
						   isol_loc, &isol_loc_avail, nz_rhs, lsol_loc, schur_mloc, schur_nloc, schur_lld, mblock, nblock,
						   nprow, npcol, ooc_tmpdir, ooc_prefix, write_problem, &ooc_tmpdirlen, &ooc_prefixlen, &write_problemlen);






								//mapping and pivnul_list are usually 0 except if
								// *MUMPS_ASSIGN_MAPPING / MUMPS_ASSIGN_PIVNUL_LIST was called.
								// */


								//int * mumps_par_mapping = mumps_get_mapping();
								//if (mumps_par_mapping == null) {
								//    mumps_par.mapping = null;
								//    if (*deficiency != 0)
								//        throw new ApplicationException("I dont know");
								//} else {
								//    int LenBuf = *deficiency;
								//    mumps_par.mapping = new int[LenBuf];
								//    for (int ii = 0; ii < LenBuf; ii++)
								//        mumps_par.mapping[ii] = mumps_par_mapping[ii];
								//}

								int* mumps_par_mapping = mumps_get_mapping();
								mumps_par.sym_perm = new int[mumps_par.nz];
								if (mumps_par_mapping != null) {
									for (int ii = 0; ii < mumps_par.nz; ii++)
										mumps_par.mapping[ii] = mumps_par_mapping[ii];
								}

                                //int* mumps_par_pivnul_list = mumps_get_pivnul_list();
                                //if (mumps_par_pivnul_list == null) {
                                //    mumps_par.pivnul_list = null;
                                //    if (*deficiency != 0)
                                //        throw new ApplicationException("I dont know");
                                //} else {
                                //    //throw new NotImplementedException();
                                //    int LenBuf = *deficiency;
                                //    mumps_par.pivnul_list = new int[LenBuf];
                                //    for (int ii = 0; ii < LenBuf; ii++)
                                //        mumps_par.pivnul_list[ii] = mumps_par_pivnul_list[ii];

                                //}

								//int* mumps_par_pivnul_list = mumps_get_pivnul_list();
								//mumps_par.pivnul_list = new int[mumps_par.n];
								//if (mumps_par_pivnul_list != null) {
								//	for (int ii = 0; ii < mumps_par.n; ii++)
								//		mumps_par.pivnul_list[ii] = mumps_par_pivnul_list[ii];
								//}


								// to get permutations computed during analysis 
								int* mumps_par_sym_perm = mumps_get_sym_perm();
								mumps_par.sym_perm = new int[mumps_par.n];
								if (mumps_par_sym_perm != null) {
									for (int ii = 0; ii < mumps_par.n; ii++)
										mumps_par.sym_perm[ii] = mumps_par_sym_perm[ii];
								}

								int* mumps_par_uns_perm = mumps_get_uns_perm();
								mumps_par.uns_perm = new int[mumps_par.n];
								if (mumps_par_uns_perm != null) {
									for (int ii = 0; ii < mumps_par.n; ii++)
										mumps_par.uns_perm[ii] = mumps_par_uns_perm[ii];
								}

								// * colsca/rowsca can either be user data or have been modified
								// * within mumps by calls to MUMPS_ASSIGN_COLSCA and/or
								// * MUMPS_ASSIGN_ROWSCA. In all cases their address is containeds
								// * in MUMPS_COLSCA_STATIC and/or MUMPS_ROWSCA_STATIC.
								// *
								// * In case of a null pointer, we also reset mumps_par->rowsca/colsca
								// * to 0 (case of JOB=-2, the Fortran pointer will be NULL but the
								// * C pointer should also be null.
								// */

								//if (rowsca_avail == no) {
								//    mumps_par.rowsca = MUMPS_ROWSCA_STATIC;
								//    if (MUMPS_ROWSCA_STATIC) {
								//        /* remember that row Scaling was computed by MUMPS */
								//        mumps_par.rowsca_from_mumps = 1;
								//    }
								//}
								//if (colsca_avail == no) {
								//    mumps_par->colsca = MUMPS_COLSCA_STATIC;
								//    if (MUMPS_COLSCA_STATIC) {
								//        /* remember that column Scaling was computed by MUMPS */
								//        mumps_par.colsca_from_mumps = 1;
								//    }
								//}


							}
						}
					}
				}
			}

		}
	}
}
