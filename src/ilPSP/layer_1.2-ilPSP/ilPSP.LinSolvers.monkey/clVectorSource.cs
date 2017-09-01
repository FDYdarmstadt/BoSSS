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

namespace ilPSP.LinSolvers.monkey.CL
{
    class clVectorSource
    {
        internal static string[] source = {
            "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n",

            "__kernel void scale(__global double* vector, double alpha, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   if(idx < size) {\n",
            "       vector[idx] *= alpha;\n",
            "   }\n",
            "}\n",
            
            "__kernel void acc(__global double* x, __global double* y, double alpha, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   if(idx < size) {\n",
            "       x[idx] += y[idx] * alpha;\n",
            "   }\n",
            "}\n",
            
            "__kernel void mew(__global double* x, __global double* y, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   if(idx < size) {\n",
            "       x[idx] *= y[idx];\n",
            "   }\n",
            "}\n",
            
            "__kernel void dnrm2(__global double* vector, __global double* result, __local double* sdata, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   int tid = get_local_id(0);\n",

            "   double value = vector[idx];\n",
            "   sdata[tid] = value * value;\n",
            "   if(idx + get_global_size(0) < size) {\n",
            "       value = vector[idx + get_global_size(0)];\n",
            "       sdata[tid] += value * value;\n",
            "   }\n",
            "   barrier(CLK_LOCAL_MEM_FENCE);\n",

            "   for(int s = get_local_size(0) / 2; s > 0; s >>= 1) {\n",
            "       if(tid < s) {\n",
            "           sdata[tid] += sdata[tid + s];\n",
            "       }\n",
            "       barrier(CLK_LOCAL_MEM_FENCE);\n",
            "   }\n",

            "   if(tid == 0) {\n",
            "       result[get_group_id(0)] = sdata[0];\n",
            "   }\n",
            "}\n",
            
            "__kernel void innerprod(__global double* x, __global double* y, __global double* result, __local double* sdata, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   int tid = get_local_id(0);\n",

            "   sdata[tid] = x[idx] * y[idx];\n",
            "   if(idx + get_global_size(0) < size) {\n",
            "       sdata[tid] += x[idx + get_global_size(0)] * y[idx + get_global_size(0)];\n",
            "   }\n",
            "   barrier(CLK_LOCAL_MEM_FENCE);\n",

            "   for(int s = get_local_size(0) / 2; s > 0; s >>= 1) {\n",
            "       if(tid < s) {\n",
            "           sdata[tid] += sdata[tid + s];\n",
            "       }\n",
            "       barrier(CLK_LOCAL_MEM_FENCE);\n",
            "   }\n",

            "   if(tid == 0) {\n",
            "       result[get_group_id(0)] = sdata[0];\n",
            "   }\n",
            "}\n",

            "__kernel void fillSendBuffer(__global double* sendBuffer, __global int* indices, __global double* data, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   if(idx < size) {\n",
            "       sendBuffer[idx] = data[indices[idx]];\n",
            "   }\n",
            "}" };
    }
}
