using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ilPSP.Utils;
namespace ilPSP {
    partial class MultidimensionalArray {
        unsafe static private void Multiply_SumGOTO_GOTO(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            const int MAX_REC = 10;
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            if (RunDim > MAX_REC)
                throw new NotSupportedException("Recursion depth not supported");
            int* iS = stackalloc int[MAX_REC];
            int rec = 0;  // recursion counter
            double** pT_Old = stackalloc double*[MAX_REC];
            double** pA_Old = stackalloc double*[MAX_REC];
            double** pB_Old = stackalloc double*[MAX_REC];
            Beginning_Multiply_SumGOTO_GOTO:
            int Li = lenRun[rec];
            for (; iS[rec] < Li; iS[rec]++) {
                if (rec >= (RunDim - 1)) {
                    {
                        Multiply_InnerSum(SumDim, pT, pA, pB, lenSum, cycSumA, cycSumB, scl, Tscl);
                    }
                } else {
                    pT_Old[rec] = pT;
                    pA_Old[rec] = pA;
                    pB_Old[rec] = pB;
                    rec++;
                    goto Beginning_Multiply_SumGOTO_GOTO;
                }
                pT += cycRunT[rec];
                pA += cycRunA[rec];
                pB += cycRunB[rec];
            }
            iS[rec] = 0;
            rec--;
            if (rec >= 0) {
                pT = pT_Old[rec] + cycRunT[rec];
                pA = pA_Old[rec] + cycRunA[rec];
                pB = pB_Old[rec] + cycRunB[rec];
                iS[rec]++;
                goto Beginning_Multiply_SumGOTO_GOTO;
            }
        }
        unsafe static void Multiply_InnerSum(int SumDim, double* pT, double* pA, double* pB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            const int MAX_REC = 10;
            if (SumDim > MAX_REC)
                throw new NotSupportedException("Recursion depth not supported");
            int* iS = stackalloc int[MAX_REC];
            int rec = 0;  // recursion counter
            double** pA_Old = stackalloc double*[MAX_REC];
            double** pB_Old = stackalloc double*[MAX_REC];
            double acc = 0.0;
            Beginning_Multiply_InnerSum:
            int Li = lenSum[rec];
            for (; iS[rec] < Li; iS[rec]++) {
                if (rec >= (SumDim - 1)) {
                    acc += (*pB) * (*pA);
                } else {
                    pA_Old[rec] = pA;
                    pB_Old[rec] = pB;
                    rec++;
                    goto Beginning_Multiply_InnerSum;
                }
                pA += cycSumA[rec];
                pB += cycSumB[rec];
            }
            iS[rec] = 0;
            rec--;
            if (rec >= 0) {
                pA = pA_Old[rec] + cycSumA[rec];
                pB = pB_Old[rec] + cycSumB[rec];
                iS[rec]++;
                goto Beginning_Multiply_InnerSum;
            }
            *pT = acc * scl + Tscl * (*pT);
        }
        unsafe static private void Multiply_Sum0_GOTO(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            const int MAX_REC = 10;
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            if (RunDim > MAX_REC)
                throw new NotSupportedException("Recursion depth not supported");
            int* iS = stackalloc int[MAX_REC];
            int rec = 0;  // recursion counter
            double** pT_Old = stackalloc double*[MAX_REC];
            double** pA_Old = stackalloc double*[MAX_REC];
            double** pB_Old = stackalloc double*[MAX_REC];
            Beginning_Multiply_Sum0_GOTO:
            int Li = lenRun[rec];
            for (; iS[rec] < Li; iS[rec]++) {
                if (rec >= (RunDim - 1)) {
                    {
                        *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                    }
                } else {
                    pT_Old[rec] = pT;
                    pA_Old[rec] = pA;
                    pB_Old[rec] = pB;
                    rec++;
                    goto Beginning_Multiply_Sum0_GOTO;
                }
                pT += cycRunT[rec];
                pA += cycRunA[rec];
                pB += cycRunB[rec];
            }
            iS[rec] = 0;
            rec--;
            if (rec >= 0) {
                pT = pT_Old[rec] + cycRunT[rec];
                pA = pA_Old[rec] + cycRunA[rec];
                pB = pB_Old[rec] + cycRunB[rec];
                iS[rec]++;
                goto Beginning_Multiply_Sum0_GOTO;
            }
        }
        unsafe static private void Multiply_Sum1_GOTO(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            const int MAX_REC = 10;
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            if (RunDim > MAX_REC)
                throw new NotSupportedException("Recursion depth not supported");
            int* iS = stackalloc int[MAX_REC];
            int rec = 0;  // recursion counter
            double** pT_Old = stackalloc double*[MAX_REC];
            double** pA_Old = stackalloc double*[MAX_REC];
            double** pB_Old = stackalloc double*[MAX_REC];
            Beginning_Multiply_Sum1_GOTO:
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int Li = lenRun[rec];
            for (; iS[rec] < Li; iS[rec]++) {
                if (rec >= (RunDim - 1)) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += (*pB) * (*pA);
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                } else {
                    pT_Old[rec] = pT;
                    pA_Old[rec] = pA;
                    pB_Old[rec] = pB;
                    rec++;
                    goto Beginning_Multiply_Sum1_GOTO;
                }
                pT += cycRunT[rec];
                pA += cycRunA[rec];
                pB += cycRunB[rec];
            }
            iS[rec] = 0;
            rec--;
            if (rec >= 0) {
                pT = pT_Old[rec] + cycRunT[rec];
                pA = pA_Old[rec] + cycRunA[rec];
                pB = pB_Old[rec] + cycRunB[rec];
                iS[rec]++;
                goto Beginning_Multiply_Sum1_GOTO;
            }
        }
        unsafe static private void Multiply_Sum2_GOTO(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            const int MAX_REC = 10;
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            if (RunDim > MAX_REC)
                throw new NotSupportedException("Recursion depth not supported");
            int* iS = stackalloc int[MAX_REC];
            int rec = 0;  // recursion counter
            double** pT_Old = stackalloc double*[MAX_REC];
            double** pA_Old = stackalloc double*[MAX_REC];
            double** pB_Old = stackalloc double*[MAX_REC];
            Beginning_Multiply_Sum2_GOTO:
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            int Li = lenRun[rec];
            for (; iS[rec] < Li; iS[rec]++) {
                if (rec >= (RunDim - 1)) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            double* pA_Old_k1 = pA;
                            double* pB_Old_k1 = pB;
                            for (int k1 = K1; k1 > 0; k1--) {
                                acc += (*pB) * (*pA);
                                pA += cAk1;
                                pB += cBk1;
                            }
                            pA = pA_Old_k1 + cAk0;
                            pB = pB_Old_k1 + cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                } else {
                    pT_Old[rec] = pT;
                    pA_Old[rec] = pA;
                    pB_Old[rec] = pB;
                    rec++;
                    goto Beginning_Multiply_Sum2_GOTO;
                }
                pT += cycRunT[rec];
                pA += cycRunA[rec];
                pB += cycRunB[rec];
            }
            iS[rec] = 0;
            rec--;
            if (rec >= 0) {
                pT = pT_Old[rec] + cycRunT[rec];
                pA = pA_Old[rec] + cycRunA[rec];
                pB = pB_Old[rec] + cycRunB[rec];
                iS[rec]++;
                goto Beginning_Multiply_Sum2_GOTO;
            }
        }
        unsafe static private void Multiply_Sum0_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            {
                *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
            }
        }
        unsafe static private void Multiply_Sum1_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            {
                double acc = 0.0;
                double* pA_Old_k0 = pA;
                double* pB_Old_k0 = pB;
                for (int k0 = K0; k0 > 0; k0--) {
                    acc += (*pB) * (*pA);
                    pA += cAk0;
                    pB += cBk0;
                }
                *pT = acc * scl + (*pT) * Tscl;
                pA = pA_Old_k0;
                pB = pB_Old_k0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll1_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            {
                double acc = 0.0;
                acc += pA[0 * cAk0] * pB[0 * cBk0];
                *pT = acc * scl + (*pT) * Tscl;
            }
        }
        unsafe static private void Multiply_Sum1Unroll2_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            {
                double acc = 0.0;
                acc += pA[0 * cAk0] * pB[0 * cBk0];
                acc += pA[1 * cAk0] * pB[1 * cBk0];
                *pT = acc * scl + (*pT) * Tscl;
            }
        }
        unsafe static private void Multiply_Sum1Unroll3_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            {
                double acc = 0.0;
                acc += pA[0 * cAk0] * pB[0 * cBk0];
                acc += pA[1 * cAk0] * pB[1 * cBk0];
                acc += pA[2 * cAk0] * pB[2 * cBk0];
                *pT = acc * scl + (*pT) * Tscl;
            }
        }
        unsafe static private void Multiply_Sum2_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            {
                double acc = 0.0;
                double* pA_Old_k0 = pA;
                double* pB_Old_k0 = pB;
                for (int k0 = K0; k0 > 0; k0--) {
                    double* pA_Old_k1 = pA;
                    double* pB_Old_k1 = pB;
                    for (int k1 = K1; k1 > 0; k1--) {
                        acc += (*pB) * (*pA);
                        pA += cAk1;
                        pB += cBk1;
                    }
                    pA = pA_Old_k1 + cAk0;
                    pB = pB_Old_k1 + cBk0;
                }
                *pT = acc * scl + (*pT) * Tscl;
                pA = pA_Old_k0;
                pB = pB_Old_k0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll1_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            {
                double acc = 0.0;
                double* pA_Old_k0 = pA;
                double* pB_Old_k0 = pB;
                for (int k0 = K0; k0 > 0; k0--) {
                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                    pA += cAk0;
                    pB += cBk0;
                }
                *pT = acc * scl + (*pT) * Tscl;
                pA = pA_Old_k0;
                pB = pB_Old_k0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll2_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            {
                double acc = 0.0;
                double* pA_Old_k0 = pA;
                double* pB_Old_k0 = pB;
                for (int k0 = K0; k0 > 0; k0--) {
                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                    acc += pA[1 * cAk1] * pB[1 * cBk1];
                    pA += cAk0;
                    pB += cBk0;
                }
                *pT = acc * scl + (*pT) * Tscl;
                pA = pA_Old_k0;
                pB = pB_Old_k0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll3_FOR0(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 0);
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            {
                double acc = 0.0;
                double* pA_Old_k0 = pA;
                double* pB_Old_k0 = pB;
                for (int k0 = K0; k0 > 0; k0--) {
                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                    acc += pA[1 * cAk1] * pB[1 * cBk1];
                    acc += pA[2 * cAk1] * pB[2 * cBk1];
                    pA += cAk0;
                    pB += cBk0;
                }
                *pT = acc * scl + (*pT) * Tscl;
                pA = pA_Old_k0;
                pB = pB_Old_k0;
            }
        }
        unsafe static private void Multiply_Sum0_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum1_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += (*pB) * (*pA);
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll1_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    acc += pA[0 * cAk0] * pB[0 * cBk0];
                    *pT = acc * scl + (*pT) * Tscl;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll2_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    acc += pA[0 * cAk0] * pB[0 * cBk0];
                    acc += pA[1 * cAk0] * pB[1 * cBk0];
                    *pT = acc * scl + (*pT) * Tscl;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll3_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    acc += pA[0 * cAk0] * pB[0 * cBk0];
                    acc += pA[1 * cAk0] * pB[1 * cBk0];
                    acc += pA[2 * cAk0] * pB[2 * cBk0];
                    *pT = acc * scl + (*pT) * Tscl;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum2_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        double* pA_Old_k1 = pA;
                        double* pB_Old_k1 = pB;
                        for (int k1 = K1; k1 > 0; k1--) {
                            acc += (*pB) * (*pA);
                            pA += cAk1;
                            pB += cBk1;
                        }
                        pA = pA_Old_k1 + cAk0;
                        pB = pB_Old_k1 + cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll1_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += pA[0 * cAk1] * pB[0 * cBk1];
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll2_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += pA[0 * cAk1] * pB[0 * cBk1];
                        acc += pA[1 * cAk1] * pB[1 * cBk1];
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll3_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += pA[0 * cAk1] * pB[0 * cBk1];
                        acc += pA[1 * cAk1] * pB[1 * cBk1];
                        acc += pA[2 * cAk1] * pB[2 * cBk1];
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
                pT += cTi0;
                pA += cAi0;
                pB += cBi0;
            }
        }
        unsafe static private void MultiplyWTrafo_Sum0_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += (*pB) * (*pA);
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll1_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    acc += pA[0 * cAk0] * pB[0 * cBk0];
                    *pT = acc * scl + (*pT) * Tscl;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll2_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    acc += pA[0 * cAk0] * pB[0 * cBk0];
                    acc += pA[1 * cAk0] * pB[1 * cBk0];
                    *pT = acc * scl + (*pT) * Tscl;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll3_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    acc += pA[0 * cAk0] * pB[0 * cBk0];
                    acc += pA[1 * cAk0] * pB[1 * cBk0];
                    acc += pA[2 * cAk0] * pB[2 * cBk0];
                    *pT = acc * scl + (*pT) * Tscl;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        double* pA_Old_k1 = pA;
                        double* pB_Old_k1 = pB;
                        for (int k1 = K1; k1 > 0; k1--) {
                            acc += (*pB) * (*pA);
                            pA += cAk1;
                            pB += cBk1;
                        }
                        pA = pA_Old_k1 + cAk0;
                        pB = pB_Old_k1 + cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll1_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += pA[0 * cAk1] * pB[0 * cBk1];
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll2_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += pA[0 * cAk1] * pB[0 * cBk1];
                        acc += pA[1 * cAk1] * pB[1 * cBk1];
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll3_FOR1(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 1);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                {
                    double acc = 0.0;
                    double* pA_Old_k0 = pA;
                    double* pB_Old_k0 = pB;
                    for (int k0 = K0; k0 > 0; k0--) {
                        acc += pA[0 * cAk1] * pB[0 * cBk1];
                        acc += pA[1 * cAk1] * pB[1 * cBk1];
                        acc += pA[2 * cAk1] * pB[2 * cBk1];
                        pA += cAk0;
                        pB += cBk0;
                    }
                    *pT = acc * scl + (*pT) * Tscl;
                    pA = pA_Old_k0;
                    pB = pB_Old_k0;
                }
            }
        }
        unsafe static private void Multiply_Sum0_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += (*pB) * (*pA);
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll1_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        acc += pA[0 * cAk0] * pB[0 * cBk0];
                        *pT = acc * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll2_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        acc += pA[0 * cAk0] * pB[0 * cBk0];
                        acc += pA[1 * cAk0] * pB[1 * cBk0];
                        *pT = acc * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll3_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        acc += pA[0 * cAk0] * pB[0 * cBk0];
                        acc += pA[1 * cAk0] * pB[1 * cBk0];
                        acc += pA[2 * cAk0] * pB[2 * cBk0];
                        *pT = acc * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            double* pA_Old_k1 = pA;
                            double* pB_Old_k1 = pB;
                            for (int k1 = K1; k1 > 0; k1--) {
                                acc += (*pB) * (*pA);
                                pA += cAk1;
                                pB += cBk1;
                            }
                            pA = pA_Old_k1 + cAk0;
                            pB = pB_Old_k1 + cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll1_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += pA[0 * cAk1] * pB[0 * cBk1];
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll2_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += pA[0 * cAk1] * pB[0 * cBk1];
                            acc += pA[1 * cAk1] * pB[1 * cBk1];
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll3_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += pA[0 * cAk1] * pB[0 * cBk1];
                            acc += pA[1 * cAk1] * pB[1 * cBk1];
                            acc += pA[2 * cAk1] * pB[2 * cBk1];
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void MultiplyWTrafo_Sum0_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += (*pB) * (*pA);
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll1_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        acc += pA[0 * cAk0] * pB[0 * cBk0];
                        *pT = acc * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll2_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        acc += pA[0 * cAk0] * pB[0 * cBk0];
                        acc += pA[1 * cAk0] * pB[1 * cBk0];
                        *pT = acc * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll3_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        acc += pA[0 * cAk0] * pB[0 * cBk0];
                        acc += pA[1 * cAk0] * pB[1 * cBk0];
                        acc += pA[2 * cAk0] * pB[2 * cBk0];
                        *pT = acc * scl + (*pT) * Tscl;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            double* pA_Old_k1 = pA;
                            double* pB_Old_k1 = pB;
                            for (int k1 = K1; k1 > 0; k1--) {
                                acc += (*pB) * (*pA);
                                pA += cAk1;
                                pB += cBk1;
                            }
                            pA = pA_Old_k1 + cAk0;
                            pB = pB_Old_k1 + cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll1_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += pA[0 * cAk1] * pB[0 * cBk1];
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll2_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += pA[0 * cAk1] * pB[0 * cBk1];
                            acc += pA[1 * cAk1] * pB[1 * cBk1];
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll3_FOR2(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for (int k0 = K0; k0 > 0; k0--) {
                            acc += pA[0 * cAk1] * pB[0 * cBk1];
                            acc += pA[1 * cAk1] * pB[1 * cBk1];
                            acc += pA[2 * cAk1] * pB[2 * cBk1];
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
            }
        }
        unsafe static private void Multiply_Sum0_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += (*pB) * (*pA);
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll1_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            acc += pA[0 * cAk0] * pB[0 * cBk0];
                            *pT = acc * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll2_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            acc += pA[0 * cAk0] * pB[0 * cBk0];
                            acc += pA[1 * cAk0] * pB[1 * cBk0];
                            *pT = acc * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll3_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            acc += pA[0 * cAk0] * pB[0 * cBk0];
                            acc += pA[1 * cAk0] * pB[1 * cBk0];
                            acc += pA[2 * cAk0] * pB[2 * cBk0];
                            *pT = acc * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                double* pA_Old_k1 = pA;
                                double* pB_Old_k1 = pB;
                                for (int k1 = K1; k1 > 0; k1--) {
                                    acc += (*pB) * (*pA);
                                    pA += cAk1;
                                    pB += cBk1;
                                }
                                pA = pA_Old_k1 + cAk0;
                                pB = pB_Old_k1 + cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll1_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += pA[0 * cAk1] * pB[0 * cBk1];
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll2_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += pA[0 * cAk1] * pB[0 * cBk1];
                                acc += pA[1 * cAk1] * pB[1 * cBk1];
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll3_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += pA[0 * cAk1] * pB[0 * cBk1];
                                acc += pA[1 * cAk1] * pB[1 * cBk1];
                                acc += pA[2 * cAk1] * pB[2 * cBk1];
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void MultiplyWTrafo_Sum0_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += (*pB) * (*pA);
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll1_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            acc += pA[0 * cAk0] * pB[0 * cBk0];
                            *pT = acc * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll2_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            acc += pA[0 * cAk0] * pB[0 * cBk0];
                            acc += pA[1 * cAk0] * pB[1 * cBk0];
                            *pT = acc * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll3_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            acc += pA[0 * cAk0] * pB[0 * cBk0];
                            acc += pA[1 * cAk0] * pB[1 * cBk0];
                            acc += pA[2 * cAk0] * pB[2 * cBk0];
                            *pT = acc * scl + (*pT) * Tscl;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                double* pA_Old_k1 = pA;
                                double* pB_Old_k1 = pB;
                                for (int k1 = K1; k1 > 0; k1--) {
                                    acc += (*pB) * (*pA);
                                    pA += cAk1;
                                    pB += cBk1;
                                }
                                pA = pA_Old_k1 + cAk0;
                                pB = pB_Old_k1 + cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll1_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += pA[0 * cAk1] * pB[0 * cBk1];
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll2_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += pA[0 * cAk1] * pB[0 * cBk1];
                                acc += pA[1 * cAk1] * pB[1 * cBk1];
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll3_FOR3(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 3);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        {
                            double acc = 0.0;
                            double* pA_Old_k0 = pA;
                            double* pB_Old_k0 = pB;
                            for (int k0 = K0; k0 > 0; k0--) {
                                acc += pA[0 * cAk1] * pB[0 * cBk1];
                                acc += pA[1 * cAk1] * pB[1 * cBk1];
                                acc += pA[2 * cAk1] * pB[2 * cBk1];
                                pA += cAk0;
                                pB += cBk0;
                            }
                            *pT = acc * scl + (*pT) * Tscl;
                            pA = pA_Old_k0;
                            pB = pB_Old_k0;
                        }
                        pT += cTi2;
                        pA += cAi2;
                        pB += cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void Multiply_Sum0_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += (*pB) * (*pA);
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll1_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                acc += pA[0 * cAk0] * pB[0 * cBk0];
                                *pT = acc * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll2_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                acc += pA[0 * cAk0] * pB[0 * cBk0];
                                acc += pA[1 * cAk0] * pB[1 * cBk0];
                                *pT = acc * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum1Unroll3_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                acc += pA[0 * cAk0] * pB[0 * cBk0];
                                acc += pA[1 * cAk0] * pB[1 * cBk0];
                                acc += pA[2 * cAk0] * pB[2 * cBk0];
                                *pT = acc * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    double* pA_Old_k1 = pA;
                                    double* pB_Old_k1 = pB;
                                    for (int k1 = K1; k1 > 0; k1--) {
                                        acc += (*pB) * (*pA);
                                        pA += cAk1;
                                        pB += cBk1;
                                    }
                                    pA = pA_Old_k1 + cAk0;
                                    pB = pB_Old_k1 + cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll1_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll2_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                                    acc += pA[1 * cAk1] * pB[1 * cBk1];
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void Multiply_Sum2Unroll3_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            double* pT = pT_org, pA = pA_org, pB = pB_org;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = I0; i0 > 0; i0--) {
                double* pT_Old_i0 = pT;
                double* pA_Old_i0 = pA;
                double* pB_Old_i0 = pB;
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                                    acc += pA[1 * cAk1] * pB[1 * cBk1];
                                    acc += pA[2 * cAk1] * pB[2 * cBk1];
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
                pT = pT_Old_i0 + cTi0;
                pA = pA_Old_i0 + cAi0;
                pB = pB_Old_i0 + cBi0;
            }
        }
        unsafe static private void MultiplyWTrafo_Sum0_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                *pT = (*pB) * (*pA) * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += (*pB) * (*pA);
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll1_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 1);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                acc += pA[0 * cAk0] * pB[0 * cBk0];
                                *pT = acc * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll2_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 2);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                acc += pA[0 * cAk0] * pB[0 * cBk0];
                                acc += pA[1 * cAk0] * pB[1 * cBk0];
                                *pT = acc * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum1Unroll3_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 1);
            Debug.Assert(lenSum[0] == 3);
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                acc += pA[0 * cAk0] * pB[0 * cBk0];
                                acc += pA[1 * cAk0] * pB[1 * cBk0];
                                acc += pA[2 * cAk0] * pB[2 * cBk0];
                                *pT = acc * scl + (*pT) * Tscl;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            int K0 = lenSum[0];
            int K1 = lenSum[1];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    double* pA_Old_k1 = pA;
                                    double* pB_Old_k1 = pB;
                                    for (int k1 = K1; k1 > 0; k1--) {
                                        acc += (*pB) * (*pA);
                                        pA += cAk1;
                                        pB += cBk1;
                                    }
                                    pA = pA_Old_k1 + cAk0;
                                    pB = pB_Old_k1 + cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll1_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll2_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 2);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                                    acc += pA[1 * cAk1] * pB[1 * cBk1];
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        unsafe static private void MultiplyWTrafo_Sum2Unroll3_FOR4(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            double* pT, pA, pB;
            Debug.Assert(RunDim == 4);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            int I2 = lenRun[2];
            int cTi2 = cycRunT[2];
            int cAi2 = cycRunA[2];
            int cBi2 = cycRunB[2];
            int I3 = lenRun[3];
            int cTi3 = cycRunT[3];
            int cAi3 = cycRunA[3];
            int cBi3 = cycRunB[3];
            Debug.Assert(SumDim == 2);
            Debug.Assert(lenSum[1] == 3);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            int cAk1 = cycSumA[1];
            int cBk1 = cycSumB[1];
            for (int i0 = 0; i0 < I0; i0++) {
                int Trf_i0_T = (Trf_T[i0 * trfCycle_T + trfPreOffset_T] + trfPostOffset_T) * trfT0sw;
                int Trf_i0_A = (Trf_A[i0 * trfCycle_A + trfPreOffset_A] + trfPostOffset_A) * trfA0sw;
                int Trf_i0_B = (Trf_B[i0 * trfCycle_B + trfPreOffset_B] + trfPostOffset_B) * trfB0sw;
                if (Trf_i0_T < 0) continue;
                if (Trf_i0_A < 0) continue;
                if (Trf_i0_B < 0) continue;
                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);
                for (int i1 = I1; i1 > 0; i1--) {
                    double* pT_Old_i1 = pT;
                    double* pA_Old_i1 = pA;
                    double* pB_Old_i1 = pB;
                    for (int i2 = I2; i2 > 0; i2--) {
                        double* pT_Old_i2 = pT;
                        double* pA_Old_i2 = pA;
                        double* pB_Old_i2 = pB;
                        for (int i3 = I3; i3 > 0; i3--) {
                            {
                                double acc = 0.0;
                                double* pA_Old_k0 = pA;
                                double* pB_Old_k0 = pB;
                                for (int k0 = K0; k0 > 0; k0--) {
                                    acc += pA[0 * cAk1] * pB[0 * cBk1];
                                    acc += pA[1 * cAk1] * pB[1 * cBk1];
                                    acc += pA[2 * cAk1] * pB[2 * cBk1];
                                    pA += cAk0;
                                    pB += cBk0;
                                }
                                *pT = acc * scl + (*pT) * Tscl;
                                pA = pA_Old_k0;
                                pB = pB_Old_k0;
                            }
                            pT += cTi3;
                            pA += cAi3;
                            pB += cBi3;
                        }
                        pT = pT_Old_i2 + cTi2;
                        pA = pA_Old_i2 + cAi2;
                        pB = pB_Old_i2 + cBi2;
                    }
                    pT = pT_Old_i1 + cTi1;
                    pA = pA_Old_i1 + cAi1;
                    pB = pB_Old_i1 + cBi1;
                }
            }
        }
        /// <summary>
        /// Selects an optimized implementation for this particular case
        /// </summary>
        unsafe static private void Multiply_Dispatch(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {
            if (SumDim == 0) {
                switch (RunDim) {
                    case 0: throw new NotSupportedException();
                    case 1: Multiply_Sum0_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                    case 2: Multiply_Sum0_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                    case 3: Multiply_Sum0_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                    case 4: Multiply_Sum0_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                    default: Multiply_Sum0_GOTO(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                }
            } else if (SumDim == 1) {
                switch (RunDim) {
                    case 0: throw new NotSupportedException();
                    case 1: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: Multiply_Sum1Unroll1_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum1Unroll2_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum1Unroll3_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum1_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    case 2: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: Multiply_Sum1Unroll1_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum1Unroll2_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum1Unroll3_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum1_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    case 3: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: Multiply_Sum1Unroll1_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum1Unroll2_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum1Unroll3_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum1_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    case 4: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: Multiply_Sum1Unroll1_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum1Unroll2_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum1Unroll3_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum1_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    default: Multiply_Sum1_GOTO(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                }
            } else if (SumDim == 2) {
                switch (RunDim) {
                    case 0: throw new NotSupportedException();
                    case 1: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: Multiply_Sum2Unroll1_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum2Unroll2_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum2Unroll3_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum2_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    case 2: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: Multiply_Sum2Unroll1_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum2Unroll2_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum2Unroll3_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum2_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    case 3: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: Multiply_Sum2Unroll1_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum2Unroll2_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum2Unroll3_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum2_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    case 4: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: Multiply_Sum2Unroll1_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 2: Multiply_Sum2Unroll2_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            case 3: Multiply_Sum2Unroll3_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                            default: Multiply_Sum2_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                        }
                    }
                    break;
                    default: Multiply_Sum2_GOTO(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl); break;
                }
            } else {
                Multiply_SumGOTO_GOTO(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl);
            }
        }
        /// <summary>
        /// Selects an optimized implementation for this particular case
        /// </summary>
        unsafe static private void MultiplyWTrafo_Dispatch(int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf_T, int* Trf_A, int* Trf_B, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B) {
            if (SumDim == 0) {
                switch (RunDim) {
                    case 0: throw new NotSupportedException();
                    case 1: MultiplyWTrafo_Sum0_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                    case 2: MultiplyWTrafo_Sum0_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                    case 3: MultiplyWTrafo_Sum0_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                    case 4:
                    MultiplyWTrafo_Sum0_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                    throw new NotImplementedException();
                }
            } else if (SumDim == 1) {
                switch (RunDim) {
                    case 0: throw new NotSupportedException();
                    case 1: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum1Unroll1_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum1Unroll2_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum1Unroll3_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum1_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    case 2: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum1Unroll1_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum1Unroll2_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum1Unroll3_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum1_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    case 3: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum1Unroll1_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum1Unroll2_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum1Unroll3_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum1_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    case 4: {
                        switch (lenSum[0]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum1Unroll1_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum1Unroll2_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum1Unroll3_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum1_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    throw new NotImplementedException();
                }
            } else if (SumDim == 2) {
                switch (RunDim) {
                    case 0: throw new NotSupportedException();
                    case 1: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum2Unroll1_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum2Unroll2_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum2Unroll3_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum2_FOR1(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    case 2: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum2Unroll1_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum2Unroll2_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum2Unroll3_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum2_FOR2(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    case 3: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum2Unroll1_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum2Unroll2_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum2Unroll3_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum2_FOR3(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    case 4: {
                        switch (lenSum[1]) {
                            case 0: return;
                            case 1: MultiplyWTrafo_Sum2Unroll1_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 2: MultiplyWTrafo_Sum2Unroll2_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            case 3: MultiplyWTrafo_Sum2Unroll3_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                            default: MultiplyWTrafo_Sum2_FOR4(RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf_T, Trf_A, Trf_B, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B); break;
                        }
                    }
                    break;
                    throw new NotImplementedException();
                }
            } else {
                throw new NotImplementedException();
            }
        }
    }
}
