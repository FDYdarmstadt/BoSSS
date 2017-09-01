echo nvcc --ptx --compiler-bindir "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin" -I"C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\include" -o CudaMatrixKernelDP.ptx -arch compute_20 CudaMatrixKernelDP.cu
nvcc --ptx -o CudaMatrixKernelDP.ptx -arch compute_20 CudaMatrixKernelDP.cu
copy *.ptx .\vectorAddDrv\Win32\Debug
copy *.ptx .\vectorAddDrv\x64\Debug
copy *.ptx .\vectorAddDrvSharp\bin\Debug
