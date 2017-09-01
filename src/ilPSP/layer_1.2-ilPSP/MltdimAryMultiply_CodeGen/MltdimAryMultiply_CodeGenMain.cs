using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using System.IO;

namespace MltdimAryMultiply_CodeGen {
    class MltdimAryMultiply_CodeGenMain {

        delegate void D();
        
        static void Main(string[] args) {
            Out = new StreamWriter("..\\..\\..\\ilPSP\\MultidimensionalArray_CodeGenMultiply.cs");

            Prolog();


            // general version:
            // =================

            // 'arbitrary' number of indices, only limited by constant 'MAX_REC';
            // summation-loops is recursive-GOTO, running-loops are recursive GOTO
            RunCycle_GOTO(SumGoto, "SumGOTO", null); // running loops with GOTO
            SumGOTOimpl(); // summation with GOTO

            // now, the unrolling for the special cases....
            // ============================================

            // explicit for -- loops should be faster than recursive GOTO

            int RUN_UNROLL = 4; // max. number of running loops
            int SUM_UNROLL = 3; // max. number of summation loops (indices over which the sum is taken)

            RunCycle_GOTO(Sum0, "Sum0", null);
            RunCycle_GOTO(Sum1, "Sum1", Sum1Prep);
            RunCycle_GOTO(Sum2, "Sum2", Sum2Prep);
                       

            for (int l = 0; l <= RUN_UNROLL; l++) {

                // "normal" summation
                // ==================

                // without any summation index
                RunCycle(l, Sum0, "Sum0", null);
            
                // with one summation index
                RunCycle(l, Sum1, "Sum1", Sum1Prep);
                for (int k = 1; k <= SUM_UNROLL; k++) {
                    // if the range of the summation index is very small (1, 2, 3, or maybe a bit more) we unroll the inner loop too
                    RunCycle(l, delegate() { Sum1Unroll(k); }, string.Format("Sum1Unroll{0}", k), delegate() { Sum1UnrollPrep(k); });
                }

                // with two summation indices...
                RunCycle(l, Sum2, "Sum2", Sum2Prep);
                for(int k = 1; k <= SUM_UNROLL; k++) {
                    RunCycle(l, delegate() { Sum2Unroll(k); }, string.Format("Sum2Unroll{0}", k), delegate() { Sum2UnrollPrep(k); });
                }


                // with index transformation
                // =========================
                if(l > 0) { // index transformation is only availabel for running loops:
                    //         therefore it is meaningless if there is no running loop
                    RunCycleWithTrafo(l, Sum0, "Sum0", null);
                    RunCycleWithTrafo(l, Sum1, "Sum1", Sum1Prep);
                    for(int k = 1; k <= SUM_UNROLL; k++) {
                        RunCycleWithTrafo(l, delegate() { Sum1Unroll(k); }, string.Format("Sum1Unroll{0}", k), delegate() { Sum1UnrollPrep(k); });
                    }
                    RunCycleWithTrafo(l, Sum2, "Sum2", Sum2Prep);
                    for(int k = 1; k <= SUM_UNROLL; k++) {
                        RunCycleWithTrafo(l, delegate() { Sum2Unroll(k); }, string.Format("Sum2Unroll{0}", k), delegate() { Sum2UnrollPrep(k); });
                    }
                }
            }
            
            // dispatch to the specific implementations...
            // ===========================================

            Dispatcher(RUN_UNROLL, SUM_UNROLL, false);
            Dispatcher(RUN_UNROLL, SUM_UNROLL, true);

            // close the file
            // ==============

            Epilog();
            Out.Close();
        }

        static string PARAMS_TYPES = "int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl";
        static string PARAMS_ONLY = "RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl";
        static string PARAMS_TYPES_WTRF = "int RunDim, int SumDim, double* pT_org, double* pA_org, double* pB_org, int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl, int* Trf, int TrfEnd, int trfT0sw, int trfA0sw, int trfB0sw, int trfPreOffset_T, int trfCycle_T, int trfPostOffset_T, int trfPreOffset_A, int trfCycle_A, int trfPostOffset_A, int trfPreOffset_B, int trfCycle_B, int trfPostOffset_B";
        static string PARAMS_ONLY_WTRF = "RunDim, SumDim, pT_org, pA_org, pB_org, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scl, Tscl, Trf, TrfEnd, trfT0sw, trfA0sw, trfB0sw, trfPreOffset_T, trfCycle_T, trfPostOffset_T, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B";

        static void Dispatcher(int RUN_UNROLL, int SUM_UNROLL, bool withTrafo) {
            string mod = withTrafo ? "WTrafo" : "";
            string params_types = withTrafo ? PARAMS_TYPES_WTRF : PARAMS_TYPES;
            string params_only = withTrafo ? PARAMS_ONLY_WTRF : PARAMS_ONLY;


            Out.WriteLine("/// <summary>");
            Out.WriteLine("/// Selects an optimized implementation for this particular case");
            Out.WriteLine("/// </summary>");
            Out.WriteLine("unsafe static private void Multiply{1}_Dispatch({0}) {{", params_types, mod);

            // ###################################################################################################################
            Out.WriteLine("if (SumDim == 0) {");  // 0 SUMMATION CYCLE ###############################################
            Out.WriteLine("switch (RunDim) {");
            Out.WriteLine("case 0: throw new NotSupportedException();");
            for (int l = 1; l <= RUN_UNROLL; l++)
                Out.WriteLine("case {0}: Multiply{2}_Sum0_FOR{0}({1}); break;", l, params_only, mod);
            if(withTrafo)
                Out.WriteLine("throw new NotImplementedException();");
            else
                Out.WriteLine("default: Multiply{1}_Sum0_GOTO({0}); break;", params_only, mod);
            Out.WriteLine("}");
            // ###################################################################################################################
            Out.WriteLine("} else if (SumDim == 1) {");  // 1 SUMMATION CYCLE ###############################################
            
            Out.WriteLine("switch (RunDim) {");
            Out.WriteLine("case 0: throw new NotSupportedException();");
            for (int l = 1; l <= RUN_UNROLL; l++) {
                Out.WriteLine("case {0}: {{", l);

                Out.WriteLine("switch(lenSum[0]) {");
                Out.WriteLine("case 0: return;");

                for (int k = 1; k <= SUM_UNROLL; k++) {
                    Out.WriteLine("case {0}: Multiply{3}_Sum1Unroll{0}_FOR{1}({2}); break;", k, l, params_only, mod);
                }

                Out.WriteLine("default: Multiply{2}_Sum1_FOR{0}({1}); break;", l, params_only, mod);
                Out.WriteLine("}");
                Out.WriteLine("} break;");
            }
            if(withTrafo)
                Out.WriteLine("throw new NotImplementedException();");
            else
                Out.WriteLine("default: Multiply{1}_Sum1_GOTO({0}); break;", params_only, mod);
            Out.WriteLine("}");
            
            // ###################################################################################################################
            Out.WriteLine("} else if (SumDim == 2) {");  // >2 SUMMATION CYCLE (GOTO) #####################################

            Out.WriteLine("switch (RunDim) {");
            Out.WriteLine("case 0: throw new NotSupportedException();");
            for(int l = 1; l <= RUN_UNROLL; l++) {
                Out.WriteLine("case {0}: {{", l);

                Out.WriteLine("switch(lenSum[1]) {");
                Out.WriteLine("case 0: return;");
                
                for(int k = 1; k <= SUM_UNROLL; k++) {
                    Out.WriteLine("case {0}: Multiply{3}_Sum2Unroll{0}_FOR{1}({2}); break;", k, l, params_only, mod);
                }

                Out.WriteLine("default: Multiply{2}_Sum2_FOR{0}({1}); break;", l, params_only, mod);
                Out.WriteLine("}");
                Out.WriteLine("} break;");
            }
            if(withTrafo)
                Out.WriteLine("throw new NotImplementedException();");
            else
                Out.WriteLine("default: Multiply{1}_Sum2_GOTO({0}); break;", params_only, mod);
            Out.WriteLine("}");

            // ###################################################################################################################
            Out.WriteLine("} else {");  // >2 SUMMATION CYCLE (GOTO) ###############################################
            if(withTrafo)
                Out.WriteLine("throw new NotImplementedException();");
            else 
                Out.WriteLine("Multiply{1}_SumGOTO_GOTO({0});", params_only, mod);
            Out.WriteLine("}");


            Out.WriteLine("}");
        }


        static TextWriter Out;

        static void SumGoto() {
            Out.WriteLine("{");
            Out.WriteLine("Multiply_InnerSum(SumDim, pT, pA, pB, lenSum, cycSumA, cycSumB, scl, Tscl);");
            Out.WriteLine("}");
        }

        static void Sum0() {
            Out.WriteLine("{");
            Out.WriteLine("*pT = (*pB)*(*pA)*scl + (*pT)*Tscl;");
            Out.WriteLine("}");
        }

        static void Sum1() {
            Out.WriteLine("{");
            Out.WriteLine("double acc = 0.0;");
            Out.WriteLine("double* pA_Old_k0 = pA;");
            Out.WriteLine("double* pB_Old_k0 = pB;");
            Out.WriteLine("for (int k0 = K0; k0 > 0; k0--) {");
            Out.WriteLine("acc += (*pB)*(*pA);");
            Out.WriteLine("pA += cAk0;");
            Out.WriteLine("pB += cBk0;");
            Out.WriteLine("}");
            Out.WriteLine("*pT = acc*scl + (*pT)*Tscl;");
            Out.WriteLine("pA = pA_Old_k0;");
            Out.WriteLine("pB = pB_Old_k0;");
            Out.WriteLine("}");
        }

        static void Sum1Prep() {
            Out.WriteLine("Debug.Assert(SumDim == 1);");
            Out.WriteLine("int K0 = lenSum[0];");
            Out.WriteLine("int cAk0 = cycSumA[0];");
            Out.WriteLine("int cBk0 = cycSumB[0];");
        }

        static void Sum2() {
            Out.WriteLine("{");
            Out.WriteLine("    double acc = 0.0;");
            Out.WriteLine("    double* pA_Old_k0 = pA;");
            Out.WriteLine("    double* pB_Old_k0 = pB;");
            Out.WriteLine("    for (int k0 = K0; k0 > 0; k0--) {");
            Out.WriteLine("        double* pA_Old_k1 = pA;");
            Out.WriteLine("        double* pB_Old_k1 = pB;");
            Out.WriteLine("        for (int k1 = K1; k1 > 0; k1--) {");
            Out.WriteLine("            acc += (*pB)*(*pA);");
            Out.WriteLine("            pA += cAk1;");
            Out.WriteLine("            pB += cBk1;");
            Out.WriteLine("        }");
            Out.WriteLine("        pA = pA_Old_k1 + cAk0;");
            Out.WriteLine("        pB = pB_Old_k1 + cBk0;");
            Out.WriteLine("    }");
            Out.WriteLine("    *pT = acc*scl + (*pT)*Tscl;");
            Out.WriteLine("    pA = pA_Old_k0;");
            Out.WriteLine("    pB = pB_Old_k0;");
            Out.WriteLine("}");
        }



        static void Sum2Prep() {
            Out.WriteLine("Debug.Assert(SumDim == 2);");
            Out.WriteLine("int K0 = lenSum[0];");
            Out.WriteLine("int K1 = lenSum[1];");
            Out.WriteLine("int cAk0 = cycSumA[0];");
            Out.WriteLine("int cBk0 = cycSumB[0];");
            Out.WriteLine("int cAk1 = cycSumA[1];");
            Out.WriteLine("int cBk1 = cycSumB[1];");
        }

        static void Sum2Unroll(int k) {
            Out.WriteLine("{");
            Out.WriteLine("    double acc = 0.0;");
            Out.WriteLine("    double* pA_Old_k0 = pA;");
            Out.WriteLine("    double* pB_Old_k0 = pB;");
            Out.WriteLine("    for (int k0 = K0; k0 > 0; k0--) {");
            for(int i = 0; i < k; i++)
                Out.WriteLine("    acc += pA[{0}*cAk1]*pB[{0}*cBk1];", i);
            Out.WriteLine("        pA += cAk0;");
            Out.WriteLine("        pB += cBk0;");
            Out.WriteLine("    }");
            Out.WriteLine("    *pT = acc*scl + (*pT)*Tscl;");
            Out.WriteLine("    pA = pA_Old_k0;");
            Out.WriteLine("    pB = pB_Old_k0;");
            Out.WriteLine("}");
        }

        static void Sum2UnrollPrep(int k) {
            Out.WriteLine("Debug.Assert(SumDim == 2);");
            Out.WriteLine("Debug.Assert(lenSum[1] == {0});", k);
            Out.WriteLine("int K0 = lenSum[0];");
            Out.WriteLine("int cAk0 = cycSumA[0];");
            Out.WriteLine("int cBk0 = cycSumB[0];");
            Out.WriteLine("int cAk1 = cycSumA[1];");
            Out.WriteLine("int cBk1 = cycSumB[1];");
        }

        static void Sum1Unroll(int k) {
            Out.WriteLine("{");
            Out.WriteLine("double acc = 0.0;");
            for(int i = 0; i < k; i++)
                Out.WriteLine("acc += pA[{0}*cAk0]*pB[{0}*cBk0];", i);
            Out.WriteLine("*pT = acc*scl + (*pT)*Tscl;");
            Out.WriteLine("}");
        }

        static void Sum1UnrollPrep(int k) {
            Out.WriteLine("Debug.Assert(SumDim == 1);");
            Out.WriteLine("Debug.Assert(lenSum[0] == {0});", k);
            Out.WriteLine("int cAk0 = cycSumA[0];");
            Out.WriteLine("int cBk0 = cycSumB[0];");
        }


        static void RunCycle(int NoOfLoops, D inner, string innerName, D inner_prepare) {
            Out.WriteLine("unsafe static private void Multiply_{1}_FOR{0}({2}) {{", NoOfLoops, innerName, PARAMS_TYPES);
            Out.WriteLine("double* pT = pT_org, pA = pA_org, pB = pB_org;");
            Out.WriteLine("Debug.Assert(RunDim == {0});", NoOfLoops);
            for(int i = 0; i < NoOfLoops; i++) {
                Out.WriteLine("int I{0} = lenRun[{0}];", i);
                Out.WriteLine("int cTi{0} = cycRunT[{0}];", i);
                Out.WriteLine("int cAi{0} = cycRunA[{0}];", i);
                Out.WriteLine("int cBi{0} = cycRunB[{0}];", i);
            }

            if(inner_prepare != null)
                inner_prepare();

            for(int i = 0; i < NoOfLoops; i++) {
                Out.WriteLine("for (int i{0} = I{0}; i{0} > 0; i{0}--) {{", i);
                if(i < (NoOfLoops - 1)) {
                    Out.WriteLine("double* pT_Old_i{0} = pT;", i);
                    Out.WriteLine("double* pA_Old_i{0} = pA;", i);
                    Out.WriteLine("double* pB_Old_i{0} = pB;", i);
                }
            }

            inner();

            for(int i = NoOfLoops - 1; i >= 0; i--) {
                if(i < (NoOfLoops - 1)) {
                    Out.WriteLine("pT = pT_Old_i{0} + cTi{0};", i);
                    Out.WriteLine("pA = pA_Old_i{0} + cAi{0};", i);
                    Out.WriteLine("pB = pB_Old_i{0} + cBi{0};", i);
                } else {
                    Out.WriteLine("pT += cTi{0};", i);
                    Out.WriteLine("pA += cAi{0};", i);
                    Out.WriteLine("pB += cBi{0};", i);
                }
                Out.WriteLine("}");
            }

            Out.WriteLine("}");
        }

        static void RunCycleWithTrafo(int NoOfLoops, D inner, string innerName, D inner_prepare) {
            Out.WriteLine("unsafe static private void MultiplyWTrafo_{1}_FOR{0}({2}) {{", NoOfLoops, innerName, PARAMS_TYPES_WTRF);
            Out.WriteLine("double* pT, pA, pB;");
            Out.WriteLine("Debug.Assert(RunDim == {0});", NoOfLoops);
            for(int i = 0; i < NoOfLoops; i++) {
                Out.WriteLine("int I{0} = lenRun[{0}];", i);
                Out.WriteLine("int cTi{0} = cycRunT[{0}];", i);
                Out.WriteLine("int cAi{0} = cycRunA[{0}];", i);
                Out.WriteLine("int cBi{0} = cycRunB[{0}];", i);
            }

            if(inner_prepare != null)
                inner_prepare();

            for(int i = 0; i < NoOfLoops; i++) {
                if(i == 0) {
                    Out.WriteLine("for (int i{0} = 0; i{0} < I{0}; i{0}++) {{", i);
                } else {
                    Out.WriteLine("for (int i{0} = I{0}; i{0} > 0; i{0}--) {{", i);
                }
                if(i == 0) {
                    Out.WriteLine("Debug.Assert(i0 < TrfEnd);");
                    Out.WriteLine("int Trf_i0_T = (Trf[i0 * trfCycle_T + trfPreOffset_T]  + trfPostOffset_T) * trfT0sw;");
                    Out.WriteLine("int Trf_i0_A = (Trf[i0 * trfCycle_A + trfPreOffset_A]  + trfPostOffset_A) * trfA0sw;");
                    Out.WriteLine("int Trf_i0_B = (Trf[i0 * trfCycle_B + trfPreOffset_B]  + trfPostOffset_B) * trfB0sw;");
                    Out.WriteLine("if(Trf_i0_T < 0) continue;");
                    Out.WriteLine("if(Trf_i0_A < 0) continue;");
                    Out.WriteLine("if(Trf_i0_B < 0) continue;");
                    Out.WriteLine("pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0_T);");
                    Out.WriteLine("pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0_A);");
                    Out.WriteLine("pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0_B);");
                } else if(i < (NoOfLoops - 1)) {
                    Out.WriteLine("double* pT_Old_i{0} = pT;", i);
                    Out.WriteLine("double* pA_Old_i{0} = pA;", i);
                    Out.WriteLine("double* pB_Old_i{0} = pB;", i);
                }
            }

            inner();

            for(int i = NoOfLoops - 1; i >= 0; i--) {
                if(i == 0) {
                    
                } else if(i < (NoOfLoops - 1)) {
                    Out.WriteLine("pT = pT_Old_i{0} + cTi{0};", i);
                    Out.WriteLine("pA = pA_Old_i{0} + cAi{0};", i);
                    Out.WriteLine("pB = pB_Old_i{0} + cBi{0};", i);
                } else {
                    Out.WriteLine("pT += cTi{0};", i);
                    Out.WriteLine("pA += cAi{0};", i);
                    Out.WriteLine("pB += cBi{0};", i);
                }
                Out.WriteLine("}");
            }

            Out.WriteLine("}");
        }

        static void RunCycle_GOTO(D inner, string innerName, D inner_prepare) {
            string FuncName = string.Format("Multiply_{0}_GOTO", innerName);
            Out.WriteLine("unsafe static private void {0}({1}) {{", FuncName, PARAMS_TYPES);
            Out.WriteLine("const int MAX_REC = 10;");
            Out.WriteLine("double* pT = pT_org, pA = pA_org, pB = pB_org;");
            Out.WriteLine("if (RunDim > MAX_REC)");
            Out.WriteLine("throw new NotSupportedException(\"Recursion depth not supported\");");
            Out.WriteLine("int* iS = stackalloc int[MAX_REC];");
            Out.WriteLine("int rec = 0;  // recursion counter");
            Out.WriteLine("double** pT_Old = stackalloc double*[MAX_REC];");
            Out.WriteLine("double** pA_Old = stackalloc double*[MAX_REC];");
            Out.WriteLine("double** pB_Old = stackalloc double*[MAX_REC];");
            Out.WriteLine("Beginning_{0}:", FuncName);
            if (inner_prepare != null)
                inner_prepare();

            Out.WriteLine("int Li = lenRun[rec];");
            Out.WriteLine("for (; iS[rec] < Li; iS[rec]++) {");
            Out.WriteLine("if (rec >= (RunDim - 1)) {");
            inner();
            Out.WriteLine("} else {");
            Out.WriteLine("pT_Old[rec] = pT;");
            Out.WriteLine("pA_Old[rec] = pA;");
            Out.WriteLine("pB_Old[rec] = pB;");
            Out.WriteLine("rec++;");

            Out.WriteLine("goto Beginning_{0};", FuncName);

            Out.WriteLine("}");

            Out.WriteLine("pT += cycRunT[rec];");
            Out.WriteLine("pA += cycRunA[rec];");
            Out.WriteLine("pB += cycRunB[rec];");

            Out.WriteLine("}");
            Out.WriteLine("iS[rec] = 0;");
            Out.WriteLine("rec--;");
            Out.WriteLine("if (rec >= 0) {");
            Out.WriteLine("pT = pT_Old[rec] + cycRunT[rec];");
            Out.WriteLine("pA = pA_Old[rec] + cycRunA[rec];");
            Out.WriteLine("pB = pB_Old[rec] + cycRunB[rec];");
            Out.WriteLine("iS[rec]++;");
            Out.WriteLine("goto Beginning_{0};", FuncName);
            Out.WriteLine("}");


            Out.WriteLine("}");
        }

        static void SumGOTOimpl() {
            Out.WriteLine("unsafe static void Multiply_InnerSum(int SumDim, double* pT, double* pA, double* pB, int* lenSum, int* cycSumA, int* cycSumB, double scl, double Tscl) {");
            Out.WriteLine("const int MAX_REC = 10;");
            Out.WriteLine("if (SumDim > MAX_REC)");
            Out.WriteLine("    throw new NotSupportedException(\"Recursion depth not supported\");");
            Out.WriteLine("          int* iS = stackalloc int[MAX_REC];");
            Out.WriteLine("            int rec = 0;  // recursion counter");
            Out.WriteLine("double** pA_Old = stackalloc double*[MAX_REC];");
            Out.WriteLine("double** pB_Old = stackalloc double*[MAX_REC];");

            Out.WriteLine("double acc = 0.0;");

            Out.WriteLine("Beginning_Multiply_InnerSum:");

            Out.WriteLine("int Li = lenSum[rec];");
            Out.WriteLine("for (; iS[rec] < Li; iS[rec]++) {");
            Out.WriteLine("if (rec >= (SumDim - 1)) {");
            Out.WriteLine("acc += (*pB)*(*pA);");
            Out.WriteLine("} else {");
            Out.WriteLine("pA_Old[rec] = pA;");
            Out.WriteLine("pB_Old[rec] = pB;");
            Out.WriteLine("rec++;");

            Out.WriteLine("goto Beginning_Multiply_InnerSum;");

            Out.WriteLine("}");

            Out.WriteLine("pA += cycSumA[rec];");
            Out.WriteLine("pB += cycSumB[rec];");

            Out.WriteLine("}");
            Out.WriteLine("iS[rec] = 0;");
            Out.WriteLine("rec--;");
            Out.WriteLine("if (rec >= 0) {");
            Out.WriteLine("    pA = pA_Old[rec] + cycSumA[rec];");
            Out.WriteLine("    pB = pB_Old[rec] + cycSumB[rec];");
            Out.WriteLine("    iS[rec]++;");
            Out.WriteLine("    goto Beginning_Multiply_InnerSum;");
            Out.WriteLine("}");


            Out.WriteLine("*pT = acc*scl + Tscl*(*pT);");
            Out.WriteLine("}");

        }



        static void Prolog() {
            Out.WriteLine("using System;");
            Out.WriteLine("using System.Collections.Generic;");
            Out.WriteLine("using System.Diagnostics;");
            Out.WriteLine("using System.Linq;");
            Out.WriteLine("using ilPSP.Utils;");
            Out.WriteLine("namespace ilPSP {");
            Out.WriteLine("partial class MultidimensionalArray {");
        }


        static void Epilog() {
            Out.WriteLine("}");
            Out.WriteLine("}");
            Out.WriteLine();
        }
    }
}
