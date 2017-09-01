/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: quad_test.c,v 1.14 2011/04/15 08:04:21 zlb Exp $
 *
 * Numerical quadrature test. */

#include "phg.h"
#include <stdlib.h>
#include <math.h>

#if FT_PHG == FT___FLOAT128 && (!HAVE_LIBQUADMATH || !HAVE_QUADMATH_H)
#undef Pow
static FLOAT
Pow0(FLOAT base, int p)
{
    if (p == 0) {
	return 1.0Q;
    }
    else if (p == 1) {
	return base;
    }
    else if (p & 1) {
	return Pow0(base, p - 1) * base;
    }
    else {
	FLOAT a = Pow0(base, p / 2);
	return a * a; 
    }
}

static FLOAT
Pow(FLOAT base, FLOAT power)
{
    int p = (int)power;
    assert((FLOAT)p == power && p >= 0);
    return Pow0(base, p);
}
#endif

void Preamble(FILE* TriQrFile) {
fprintf(TriQrFile, "// rules from the PHG code, http://lsec.cc.ac.cn/phg/index_en.htm\n");
      fprintf(TriQrFile, "// @article{zhang2009set,\n");
      fprintf(TriQrFile, "// title={A set of symmetric quadrature rules on triangles and tetrahedra},\n");
      fprintf(TriQrFile, "// author={Zhang, Linbo and Cui, Tao and Liu, Hui and others},\n");
      fprintf(TriQrFile, "// journal={J. Comput. Math},\n");
      fprintf(TriQrFile, "// volume={27},\n");
      fprintf(TriQrFile, "// number={1},\n");
      fprintf(TriQrFile, "// pages={89--96},\n");
      fprintf(TriQrFile, "//   year={2009}\n");
      fprintf(TriQrFile, "// }\n");
      fprintf(TriQrFile, "// for the generation of thest rules, sizeof(FLOAT) = %d, EPSILON = %.32Qe\n", sizeof(FLOAT), FLOAT_EPSILON);
}

int
main(int argc, char *argv[])
{
    QUAD *quad;
    int i, p, np;
    FLOAT a, b, error;
    const FLOAT *pw, *pp;
    BOOLEAN dump_flag = TRUE;
    FILE* dumpfile;
    char dumpfileName[1000];

    FILE* TriQrFile, *TetraQrFile;

    FLOAT TriVol, TetraVol;
    FLOAT TriVtx[3][2] = { { 0, 4.0 / 3.0 }, { -2.0 * Sqrt(3.0) / 3.0, -2.0 / 3.0 }, { 2.0 * Sqrt(3.0) / 3.0, -2.0 / 3.0 } };
    FLOAT TetraVtx[4][3];
    FLOAT w;
    FLOAT X[3];
    int k,d;
    TetraVol = (16.0/27.0)*Sqrt(2.0)*Sqrt(3.0);
    TetraVtx[0][0] = 0.0;                     TetraVtx[0][1] = 0.0;           TetraVtx[0][2] =  Sqrt(2.0);
    TetraVtx[1][0] = 0.0;                     TetraVtx[1][1] = 4.0 / 3.0;     TetraVtx[1][2] = -Sqrt(2.0) / 3.0;
    TetraVtx[2][0] = -2.0 * Sqrt(3.0) / 3.0;  TetraVtx[2][1] = -2.0 / 3.0;    TetraVtx[2][2] = -Sqrt(2.0) /3.0;
    TetraVtx[3][0] = 2.0 * Sqrt(3.0) / 3.0;   TetraVtx[3][1] =  -2.0 / 3.0;   TetraVtx[3][2] = -Sqrt(2.0) / 3.0;
        
    TriVol = 4.0*Sqrt(3.0)/3.0;

    printf("Float size = %d\n",sizeof(FLOAT));

    if(dump_flag) {
      TriQrFile = fopen("TriQr.cs", "w");
      TetraQrFile = fopen("TetraQr.cs", "w");
      Preamble(TriQrFile);
      Preamble(TetraQrFile);
    }

    phgOptionsRegisterNoArg("-dump_rules", "Dump quadrature rules to stdout",
			    &dump_flag);

    phgInit(&argc, &argv);
    
    /* 1D quadrature lambda0^p + lambda1^p */
    fprintf(stderr, "\n========================== 1D quadrature rules.\n");
    for (p = 0; p <= 30; p++) {
	quad = phgQuadGetQuad1D(p);
	np = quad->npoints;
	pp = quad->points;
	pw = quad->weights;
	b = 0.;
	if (dump_flag) {
            //printf("dim = %d, order = %d, npoints = %d\n", 1, p, np);
	  sprintf(dumpfileName, "LineP%d.txt", p);
          dumpfile = fopen(dumpfileName, "w");
        }
	for (i = 0; i < np; i++, pp += 2) {
	    if (dump_flag)
	      fprintf(dumpfile,"%0.16le %0.16le %0.16le\n",
			(double)pp[0], (double)pp[1], (double)*pw);
	    b += (Pow(pp[0],p) + Pow(pp[1],p)) * *(pw++);
	}
	/* analytic: \int lambda^p = 1 / (p + 1) */
	a = 2. / (FLOAT)(p + 1) * 1.;
	fprintf(stderr,
		"    order %2d, points %4d, result = %le, rel. error = %le\n",
		p, np, (double)b, (double)Fabs(error=(a-b)/b));
	if (Fabs(error) > 10000.0 * FLOAT_EPSILON)
	    fprintf(stderr, "*** WARNING: possibly incorrect result\n");
        if (dump_flag) {
	  fclose(dumpfile);
        }
    }

    /* 2D quadrature lambda0^p + lambda1^p + lambda2^p */
    fprintf(stderr, "\n========================== 2D quadrature rules.\n");
    for (p = 0; p <= 30; p++) {
	quad = phgQuadGetQuad2D(p);
	np = quad->npoints;
	pp = quad->points;
	pw = quad->weights;
	b = 0.;
	if (dump_flag) {
	  //printf("dim = %d, order = %d, npoints = %d\n", 2, p, np);
           sprintf(dumpfileName, "TriP%d.txt", p);
          dumpfile = fopen(dumpfileName, "w");
 fprintf(TriQrFile,"{\n");
          fprintf(TriQrFile,"QuadRule q = QuadRule.CreateEmpty(this, %d, 2);\n",np);
        }
	for (i = 0; i < np; i++, pp += 3) {
	    if (dump_flag)
		fprintf(dumpfile,"%0.16le %0.16le %0.16le %0.16le\n", (double)pp[0],
			(double)pp[1], (double)pp[2], (double)*pw);
	    b += (Pow(pp[0],p) + Pow(pp[1],p) + Pow(pp[2],p)) * *(pw);

	    w = TriVol * (*pw);
          for(d = 0; d < 2; d++)
            X[d] = 0;
          for(k = 0; k < 3; k++)
            for(d = 0; d < 3; d++)
              X[d] += pp[k]*TriVtx[k][d];
          for(d = 0; d < 2; d++) {
            //quadmath_snprintf(buffy, sizeof(buffy), "%+-#*.32Qe", width, r);
	    fprintf(TriQrFile,"q.Nodes[%d, %d] = %0.32Qe; \n",i,d,X[d]);
	  }
           fprintf(TriQrFile,"q.Weights[%d] = %0.32Qe; \n",i,w);

	}
        fprintf(TriQrFile,"q.OrderOfPrecision = %d;\n",p );
        fprintf(TriQrFile,"m_QuadRules.Add(q);\n");

	/* analytic: \int lambda^p = 1 / ((p + 1) * (p + 2)) */
	a = 3. / (FLOAT)((p + 1) * (p + 2)) * 2.;
	fprintf(stderr,
		"    order %2d, points %4d, result = %le, rel. error = %le\n",
		p, np, (double)b, (double)Fabs(error=(a - b) / b));
	if (Fabs(error) > 10000.0 * FLOAT_EPSILON) {
	    fprintf(stderr, "*** WARNING: possibly incorrect result\n");
        if(dump_flag) {
	      fprintf(TriQrFile, "//*** WARNING: possibly incorrect result, error is %0.32Qe\n", Fabs(error));
            }
        }
        fprintf(TriQrFile,"}\n");
        if (dump_flag) {
	  fclose(dumpfile);
        }
    }

    /* 3D quadrature lambda0^p + lambda1^p + lambda2^p + lambda3^p */
    fprintf(stderr, "\n========================== 3D quadrature rules.\n");
    for (p = 0; p <= 30; p++) {
	quad = phgQuadGetQuad3D(p);
	np = quad->npoints;
	pp = quad->points;
	pw = quad->weights;
	b = 0.;
	if (dump_flag) {
	  //	    printf("dim = %d, order = %d, npoints = %d\n", 3, p, np);
          sprintf(dumpfileName, "TetraP%d.txt", p);
          dumpfile = fopen(dumpfileName, "w");
          fprintf(TetraQrFile,"{\n");
          fprintf(TetraQrFile,"QuadRule q = QuadRule.CreateEmpty(this, %d, 3);\n",np);
        }
	for (i = 0; i < np; i++, pp += 4) {
	  if (dump_flag) {
		fprintf(dumpfile,"%0.32le %0.32le %0.32le %0.32le %0.16le\n",
			(double)pp[0], (double)pp[1], (double)pp[2],
			(double)pp[3], (double)*pw);
	  }
          
          w = TetraVol * (*pw);
          for(d = 0; d < 3; d++)
            X[d] = 0;
          for(k = 0; k < 4; k++)
            for(d = 0; d < 3; d++)
              X[d] += pp[k]*TetraVtx[k][d];
          for(d = 0; d < 3; d++) {
            //quadmath_snprintf(buffy, sizeof(buffy), "%+-#*.32Qe", width, r);
	    fprintf(TetraQrFile,"q.Nodes[%d, %d] = %0.32Qe; \n",i,d,X[d]);
	  }
           fprintf(TetraQrFile,"q.Weights[%d] = %0.32Qe; \n",i,w);

	    b += (Pow(pp[0],p) + Pow(pp[1],p) + Pow(pp[2],p) + Pow(pp[3],p)) *
		 *(pw++);
	}
        fprintf(TetraQrFile,"q.OrderOfPrecision = %d;\n",p );
        fprintf(TetraQrFile,"m_QuadRules.Add(q);\n");
        
	/* analytic: \int lambda^p = 1 / ((p + 1) * (p + 2) * (p + 3)) */
	a = 4. / (FLOAT)((p + 1) * (p + 2) * (p + 3)) * 6.;
	fprintf(stderr,
		"    order %2d, points %4d, result = %le, rel. error = %le\n",
		p, np, (double)b, (double)Fabs(error=(a - b) / b));
	if (Fabs(error) > 10000.0 * FLOAT_EPSILON) {
	    fprintf(stderr, "*** WARNING: possibly incorrect result\n");
            if(dump_flag) {
	      fprintf(TetraQrFile, "//*** WARNING: possibly incorrect result, error is %0.32Qe\n", Fabs(error));
            }
        }
        fprintf(TetraQrFile,"}\n");
        if (dump_flag) {
	  fclose(dumpfile);
        }
    }

    phgFinalize();

    if(dump_flag) {
      fclose(TriQrFile);
      fclose(TetraQrFile);
    }
    return 0;
}
