extern "C" __global__ void sparseMultiply(double* values, int* colIdx, int* rowStart, double* result, double* x, double alpha, double beta, int size) {
	// Dynamically allocated shared memory, should be BlockDim.x + 1 ints (see cuFuncSetSharedSize host code)
	extern __shared__ int sharedRowStart[];
	
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
		int lastIdx = min((blockIdx.x + 1) * blockDim.x, size);
		sharedRowStart[blockDim.x] = rowStart[lastIdx];
	}
	__syncthreads();
	
	if(idx < size) {
		// Multiply and sum up data of this row
		for(int i = sharedRowStart[tid]; i < sharedRowStart[tid + 1]; i++) {
			rowacc += values[i] * x[colIdx[i]];
		}

		result[idx] = rowacc * alpha;
	}
}
