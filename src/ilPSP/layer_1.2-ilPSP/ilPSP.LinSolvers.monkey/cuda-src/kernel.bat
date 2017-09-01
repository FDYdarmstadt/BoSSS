nvcc --ptx --compiler-bindir "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin" -I"C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\include" -o ./CudaMatrixKernelDP.ptx -arch compute_13 CudaMatrixKernelDP.cu
nvcc --ptx --compiler-bindir "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin" -I"C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\include" -o ./CudaVectorKernelDP.ptx -arch compute_13 CudaVectorKernelDP.cu
nvcc --cubin --compiler-bindir "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin" -I"C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\include" -o ./CudaMatrixKernelDP.cubin --ptxas-options=-v -arch sm_13 CudaMatrixKernelDP.cu
copy .\CudaMatrixKernelDP.ptx .\CudaMatrixKernelDP.txt
copy .\CudaVectorKernelDP.ptx .\CudaVectorKernelDP.txt