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
    class clMatrixSource
    {
        internal static string[] source = {
            "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n",

            "__kernel void csrMultiply(__global double* values, __global int* colIdx, __global int* rowStart, __global double* result, __global double* x, __local int* sharedRowStart, double alpha, double beta, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   int tid = get_local_id(0);\n",
            "   double rowacc = 0.0;\n",

            "   if(idx < size) {\n",
            "       sharedRowStart[tid] = rowStart[idx];\n",
            "   }\n",
            "   if(tid == 0) {\n",
            "       int loadIdx = min((int)((get_group_id(0) + 1) * get_local_size(0)), size);\n",
		    "       int storIdx = size % get_local_size(0) > 0 && idx + get_local_size(0) >= size ? size % get_local_size(0) : get_local_size(0);\n",
		    "       sharedRowStart[storIdx] = rowStart[loadIdx];\n",
	        "   }\n",
            "   barrier(CLK_LOCAL_MEM_FENCE);\n",

            "   if(idx < size) {\n",
            "       for(int i = sharedRowStart[tid]; i < sharedRowStart[tid + 1]; i++) {\n",
			"           rowacc += values[i] * x[colIdx[i]];\n",
            "       }\n",
		    "   }\n",
		    
		    "   result[idx] = result[idx] * beta + rowacc * alpha;\n",
	        "}\n",

            "__kernel void ccbcsrMultiply(__global double* cellData, __global double* xData, __global int* cellColIdx, __global double* result, __local double* sharedData, __local int* start, __local int* colIdx, double alpha, double beta, int cellsize, int cellrowsperblock, int cellsperrow, int stride, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   int tid = get_local_id(0);\n",
            "   int cellid = tid / cellsize;\n",
            "   int cid = tid % cellsize;\n",
            "   double rowacc = 0.0;\n",
            "   double value;\n",
            
            "   if(tid < cellrowsperblock) {\n",
            "       start[tid] = (cellrowsperblock * get_group_id(0) + tid) * cellsperrow;\n",
            "   }\n",
            "   barrier(CLK_LOCAL_MEM_FENCE);\n",

            "   for(int i = 0; i < cellsperrow; i++) {\n",
            "       if(tid < cellrowsperblock) {\n",
            "           colIdx[tid] = cellColIdx[start[tid] + i];\n",
            "       }\n",
            "       barrier(CLK_LOCAL_MEM_FENCE);\n",
            
            "       if(idx < size) {\n",
            "           sharedData[tid] = xData[colIdx[cellid] * cellsize + cid];\n",
            "       }\n",
            "       barrier(CLK_LOCAL_MEM_FENCE);\n",

            "       if(idx < size) {\n",
            "           for(int col = 0; col < cellsize; col++) {\n",
            "               value = cellData[(start[cellid] + i) * stride + col * cellsize + cid];\n",
            "               rowacc += value * sharedData[cellid * cellsize + col];\n",
            "           }\n",
            "       }\n",
            "       barrier(CLK_LOCAL_MEM_FENCE);\n",
            "   }\n",

            "   if(idx < size) {\n",
            "       result[idx] = result[idx] * beta + rowacc * alpha;\n",
            "   }\n",
            "}\n",

            "__kernel void ellMultiply(__global double* valData, __global int* colIdxData, __global double* xData, __global double* result, double alpha, double beta, int size, int colCount, int valStride, int colStride) {\n",
            "   int idx = get_global_id(0);\n",
            "   int tid = get_local_id(0);\n",
            "   valData    += get_group_id(0) * colCount * valStride;\n",
            "   colIdxData += get_group_id(0) * colCount * colStride;\n",
            "   int valIdx, colIdx;\n",

            "   if(idx < size) {\n",
            "       double rowacc = 0.0;\n",
            "       for(int col = 0; col < colCount; col++) {\n",
            "           valIdx = col * valStride + tid;\n",
            "           colIdx = col * colStride + tid;\n",
            "           rowacc += valData[valIdx] * xData[colIdxData[colIdx]];\n",
            "       }\n",
            "       result[idx] = result[idx] * beta + rowacc * alpha;\n",
            "   }\n",
            "}\n",

            "__kernel void mcellMultiply(__global double* valData, __global unsigned short* colIdxData, __global int* xSubStart, __global int* blockSubVector, __global double* xData, __global double* result, __local double* xSub, double alpha, double beta, int size, int colCount, int valStride, int colStride) {\n",
            "   int idx = get_global_id(0);\n",
            "   int tid = get_local_id(0);\n",
            "   valData    += get_group_id(0) * colCount * valStride;\n",
            "   colIdxData += get_group_id(0) * colCount * colStride;\n",

            "   __local int xStart;\n",
            "   __local int xLength;\n",
            "   int valIdx;\n",
            "   unsigned short colIdx;\n",

            "   if(tid == 0) {\n",
            "       xStart  = xSubStart[get_group_id(0)    ];\n",
            "       xLength = xSubStart[get_group_id(0) + 1] - xStart;\n",
            "   }\n",
            "   barrier(CLK_LOCAL_MEM_FENCE);\n",

            "   blockSubVector += xStart;\n",
            "   int ldIdx = tid;\n",
            "   while(ldIdx < xLength) {\n",
            "       xSub[ldIdx] = xData[blockSubVector[ldIdx]];\n",
            "       ldIdx += get_local_size(0);\n",
            "   }\n",
            "   barrier(CLK_LOCAL_MEM_FENCE);\n",

            "   if(idx < size) {\n",
            "       double rowacc = 0.0;\n",
            "       for(int col = 0; col < colCount; col++) {\n",
            "           valIdx = col * valStride + tid;\n",
            "           colIdx = col * colStride + tid;\n",
            "           rowacc += valData[valIdx] * xSub[colIdxData[colIdx]];\n",
            "       }\n",
            "       result[idx] = result[idx] * beta + rowacc * alpha;\n",
            "   }\n",
            "}\n",

            "__kernel void accumulateExternal(__global double* data, __global int* indices, __global double* rcvBuffer, double alpha, int size) {\n",
            "   int idx = get_global_id(0);\n",
            "   if(idx < size) {\n",
            "       data[indices[idx]] += rcvBuffer[idx] * alpha;\n",
            "   }\n",
            "}" };
    }
}
