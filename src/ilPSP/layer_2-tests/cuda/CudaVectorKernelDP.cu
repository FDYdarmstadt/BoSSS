extern "C" __global__ void scale(double* vector, double alpha, unsigned int size) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx < size) {
		vector[idx] *= alpha;
	}
}

extern "C" __global__ void acc(double* x, double* y, double alpha, int size) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx < size) {
		x[idx] += y[idx] * alpha;
	}
}

extern "C" __global__ void dnrm2(double* vector, double* result, int size) {
	__shared__ double sdata[256];

	int tid = threadIdx.x;
	int idx = blockIdx.x * blockDim.x * 2 + threadIdx.x;
	double value;

	sdata[tid] = 0.0;
	if(idx < size) {
		value = vector[idx];
		sdata[tid] += value * value;
	}
	if(idx + blockDim.x < size) {
		value = vector[idx + blockDim.x];
		sdata[tid] += value * value;
	}

	__syncthreads();

	for(int s = blockDim.x / 2; s > 0; s >>= 1) {
		if(tid < s) {
			sdata[tid] += sdata[tid + s];
		}

		__syncthreads();
	}

	if(tid == 0) {
		result[blockIdx.x] = sdata[0];
	}
}

extern "C" __global__ void innerprod(double* x, double* y, double* result, int size) {
	__shared__ double sdata[256];

	int tid = threadIdx.x;
	int idx = blockIdx.x * blockDim.x * 2 + threadIdx.x;

	sdata[tid] = 0.0;
	if(idx < size) {
		sdata[tid] += x[idx] * y[idx];
	}
	if(idx + blockDim.x < size) {
		sdata[tid] += x[idx + blockDim.x] * y[idx + blockDim.x];
	}

	__syncthreads();

	for(int s = blockDim.x / 2; s > 0; s >>= 1) {
		if(tid < s) {
			sdata[tid] += sdata[tid + s];
		}

		__syncthreads();
	}

	if(tid == 0) {
		result[blockIdx.x] = sdata[0];
	}
}
