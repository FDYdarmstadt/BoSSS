#!bin/bash

echo "---------- start of script ----------"

module load acml
module load gcc
module load openmpi/gcc/2.0.1

rm *.so

cd metis
rm *.o
mpicc -fPIC -c *.c
cd ..
mpicc -shared -o libmetis.so ./metis/*.o

cd parmetis
rm *.o
mpicc -fPIC -c *.c
cd ..
mpicc -shared -o libparmetis.so ./parmetis/*.o -L. -lmetis

echo "---------- end of script ----------"
