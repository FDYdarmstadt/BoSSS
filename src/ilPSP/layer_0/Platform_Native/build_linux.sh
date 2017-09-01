mpicc -fPIC -c BoSSS_MPI.c MPI_Exports2.c
mpif77 -shared  -o libPlatform_Native.so *.o
