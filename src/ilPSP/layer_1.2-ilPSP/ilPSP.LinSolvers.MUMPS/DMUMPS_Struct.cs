/* 
 * Copyright (C) 2010, Florian Kummer, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsmechanik
 *
 * Use, modification and distribution is subject to the Boost Software
 * License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *  
 * Authors: Florian Kummer
 * 
 */
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.MUMPS.Wrappers {
    //[StructLayout(LayoutKind.Sequential)]
    internal class DMUMPS_Struct {
        public int sym, par, job;
        public int _comm_fortran;    /* Fortran communicator */
        public int[] icntl = new int[40];
        public double[] cntl = new double[15];
        public int n;

        public int nz_alloc = 0; /* used in matlab interface to decide if we
                                    free + malloc when we have large variation */

        /* Assembled entry */
        public int nz;
        public int[] irn;
        public int[] jcn;
        public double[] a;

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

        /* Scaling (input only in this version) */
        public double[] colsca;
        public double[] rowsca;

        /* RHS, solution, ouptput data and statistics */
        public double[] rhs, redrhs, rhs_sparse, sol_loc;
        public int[] irhs_sparse, irhs_ptr, isol_loc;
        public int nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
        public int schur_mloc, schur_nloc, schur_lld;
        public int mblock, nblock, nprow, npcol;
        public int[] info = new int[40];
        public int[] infog = new int[40];
        public double[] rinfo = new double[20];
        public double[] rinfog = new double[20];

        /* Null space */
        public int deficiency;
        public int[] pivnul_list;
        public int[] mapping;

        /* Schur */
        public int size_schur;
        public int[] listvar_schur;
        public double[] schur;

        /* Internal parameters */
        public int instance_number;
        public double[] wk_user;

        /* Version number: length=14 in FORTRAN + 1 for final \0 + 1 for alignment */
        //public string version_number ;// = new version_number[MUMPS_VERSION_MAX_LEN + 1 + 1];
        /* For out-of-core */
        public string ooc_tmpdir;// = new char[256];
        public string ooc_prefix;// = new char[64];
        /* To save the matrix in matrix market format */
        public string write_problem;// = new char[256];
        public int lwk_user;
    }
}
