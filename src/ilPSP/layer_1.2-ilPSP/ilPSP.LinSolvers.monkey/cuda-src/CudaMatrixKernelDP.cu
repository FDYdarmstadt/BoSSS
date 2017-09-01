/* 
 * Copyright (C) 2010, Florian Kummer, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsmechanik
 *
 * Use, modification and distribution is subject to the Boost Software
 * License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *  
 * Authors: Christoph Busold
 * 
 */

extern __shared__ char smem[];

extern "C" __global__ void sparseMultiply(double* values, int* colIdx, int* rowStart, double* result, double* x, double alpha, double beta, int size) {
	// Dynamically allocated shared memory, should be BlockDim.x + 1 ints (see cuFuncSetSharedSize host code)
	int* sharedRowStart = (int*)smem;
	
	// Indices
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	
	double rowacc = 0.0;
	
	// Each thread loads one element of rowStart
	if(idx < size) {
		sharedRowStart[tid] = rowStart[idx];
	}
	// The first thread loads additionally the next element, needed by the last thread
	if(tid == 0) {
		int loadIdx = min((blockIdx.x + 1) * blockDim.x, size);
		int storIdx = size % blockDim.x > 0 && idx + blockDim.x >= size ? size % blockDim.x : blockDim.x;
		sharedRowStart[storIdx] = rowStart[loadIdx];
	}
	__syncthreads();
	
	if(idx < size) {
		// Multiply and sum up data of this row
		for(int i = sharedRowStart[tid]; i < sharedRowStart[tid + 1]; i++) {
			rowacc += values[i] * x[colIdx[i]];
		}
		
		result[idx] = result[idx] * beta + rowacc * alpha;
	}
}

extern "C" __global__ void accumulateExternal(double* data, int* indices, double* rcvBuffer, double alpha, int size) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(idx < size) {
		data[indices[idx]] += rcvBuffer[idx] * alpha;
	}
}

// In this kernel each block computes multiple cell rows 
// IMPORTANT: All cell rows must have the same number of cells!
//            Otherwise sync in kernel might fail, causing crash or incorrect behaviour!
extern "C" __global__ void blockMultiply2(double* cellData, double* xData, int* cellColIdx, double* result, double alpha, double beta, int cellsize, int cellrowsperblock, int cellsperrow, int stride, int size) {
	// Dynamically allocated shared memory, should be blockDim.x doubles for xData
	double* sharedData = (double*)smem;
	// Start cell index of this thread
	int* start = (int*)&sharedData[blockDim.x];
	// Column of this thread's cell
	int* colIdx = (int*)&start[cellrowsperblock];
	
	// Global index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	// Number of the cell this thread is in
	int cellid = tid / cellsize;
	// Thread index inside this cell
	int cid = tid % cellsize;
	
	double rowacc = 0.0;
	double value;
	
	// Load start index for every cell row in this block
	if(tid < cellrowsperblock) {
		start[tid] = (cellrowsperblock * blockIdx.x + tid) * cellsperrow;
	}
	__syncthreads();
	
	// Loop over all cells, discard overlapping threads inside because of sync
	for(int i = 0; i < cellsperrow; i++) {
		// Load column index for every cell
		if(tid < cellrowsperblock) {
			colIdx[tid] = cellColIdx[start[tid] + i];
		}
		__syncthreads();
		
		// No overlapping threads
		if(idx < size) {
			// Load x at colIdx location into shared memory
			// colIdx * cellsize is the start index at this column
			// cid is the row index of this thread
			sharedData[tid] = xData[colIdx[cellid] * cellsize + cid];
		}
		__syncthreads();
		
		// No overlapping threads
		if(idx < size) {
			// Loop over all columns of this cell
			for(int col = 0; col < cellsize; col++) {
				// Load value of this column
				// cell * cellsize * cellsize is the start index of the current cell
				// col * cellsize is the start index of the current column
				// cid is the row index of this thread
				value = cellData[(start[cellid] + i) * stride + col * cellsize + cid];
				// Multiply value with x from sharedMemory
				// cellid * cellsize is the offset for the cell this thread is in
				// col is the column index of this loop cycle
				rowacc += value * sharedData[cellid * cellsize + col];
			}
		}
		
		__syncthreads();
	}
	
	// No overlapping threads
	if(idx < size) {
		// Write back result 
		result[idx] = result[idx] * beta + rowacc * alpha;
	}
}

// In this kernel each block computes one cell row (block size equals cell size)
extern "C" __global__ void blockMultiply(double* cellData, double* xData, int* cellColIdx, int* cellRowStart, double* result, double dia, int cellsize, int size) {
	double* sharedData = (double*)smem;
	__shared__ int colIdx;
	__shared__ int start;
	__shared__ int end;
	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	
	double rowacc = 0.0;
	double value;
	
	if(tid == 0) {
		start = cellRowStart[blockIdx.x    ];
		end   = cellRowStart[blockIdx.x + 1];
	}
	__syncthreads();
	
	for(int cell = start; cell < end; cell++) {
		if(tid == 0) {
			colIdx = cellColIdx[cell];
		}
		__syncthreads();
		
		if(idx < size) {
			sharedData[tid] = xData[colIdx * cellsize + tid];
		}
		__syncthreads();
		
		if(idx < size) {
			for(int col = 0; col < cellsize; col++) {
				value = cellData[cell * cellsize * cellsize + col * cellsize + tid];
				rowacc += value * sharedData[col];
			}
		}
	}
	
	if(idx < size) {
		rowacc += dia * xData[idx];
		result[idx] += rowacc;
	}
}

// ELLPACKmod format
extern "C" __global__ void ellMultiply(double* valData, int* colIdxData, double* xData, double* result, double alpha, double beta, int size, int colCount, int valStride, int colStride) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// Add offsets to the start of this block's value and column data pointers
	valData    += blockIdx.x * colCount * valStride;
	colIdxData += blockIdx.x * colCount * colStride;
	
	int valIdx;
	int colIdx;
	
	// No sync in this kernel, therefore overlapping threads are discarded here
	if(idx < size) {
		double rowacc = 0.0;
		
		// Loop over all columns
		for(int col = 0; col < colCount; col++) {
			// Index of the value and column index to load
			valIdx = col * valStride + threadIdx.x;
			colIdx = col * colStride + threadIdx.x;
			// Load value and multiply with x at column of this value
			rowacc += valData[valIdx] * xData[colIdxData[colIdx]];
		}
		
		// Write result back
		result[idx] = result[idx] * beta + rowacc * alpha;
	}
}

// ManualCacheELLPACK format
extern "C" __global__ void mcellMultiply(double* valData, unsigned short* colIdxData, int* xSubStart, int* blockSubVector, double* xData, double* result, double alpha, double beta, int size, int colCount, int valStride, int colStride) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	// Add offsets to the start of this block's value and column data pointers
	valData    += blockIdx.x * colCount * valStride;
	colIdxData += blockIdx.x * colCount * colStride;
	
	double* xSub = (double*)smem;
	__shared__ int xStart;
	__shared__ int xLength;
	
	int valIdx;
	unsigned short colIdx;
	
	if(tid == 0) {
		xStart  = xSubStart[blockIdx.x    ];
		xLength = xSubStart[blockIdx.x + 1] - xStart;
	}
	
	__syncthreads();
	
	blockSubVector += xStart;
	int ldIdx = tid;
	
	while(ldIdx < xLength) {
		xSub[ldIdx] = xData[blockSubVector[ldIdx]];
		ldIdx += blockDim.x;
	}
	
	__syncthreads();
	
	// No sync inside this loop, therefore overlapping threads are discarded here
	if(idx < size) {
		double rowacc = 0.0;
		
		// Loop over all columns
		for(int col = 0; col < colCount; col++) {
			// Index of the value and column index to load
			valIdx = col * valStride + tid;
			colIdx = col * colStride + tid;
			// Load value and multiply with x at column of this value
			rowacc += valData[valIdx] * xSub[colIdxData[colIdx]];
		}
		
		// Write result back
		result[idx] = result[idx] * beta + rowacc * alpha;
	}
}
