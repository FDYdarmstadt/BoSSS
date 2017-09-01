#gcc MonoClient.c `pkg-config --cflags --libs mono`
cp ../ilPSP.ExternalBinding.CodeGenC/bin/Debug/binding.* .
cp ../ilPSP.ExternalBinding/bin/Debug/*.dll .

rm *.o
gcc monkey.c -c `pkg-config --cflags mono`
gcc binding.c -c `pkg-config --cflags mono`
mpicc example.c -c 
mpif77 fexample.f -c

mpicc  monkey.o binding.o example.o  -o cExample.out `pkg-config --libs mono` 
mpif77 monkey.o binding.o fexample.o -o fExample.out `pkg-config --libs mono` 

