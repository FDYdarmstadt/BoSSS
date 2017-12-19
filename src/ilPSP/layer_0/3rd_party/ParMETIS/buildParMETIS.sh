#!bin/bash

echo "---------- start of script ----------"

module load acml
module load gcc
module load openmpi/gcc/1.6.6

rm *.so

cd METIS
rm *.o
mpicc -fPIC -c *.c
cd ..
mpicc -shared -o libmetis.so ./METIS/*.o

cd ParMETIS
rm *.o
mpicc -fPIC -c *.c
cd ..
mpicc -shared -o libparmetis.so ./ParMETIS/*.o -L. -lmetis

echo "---------- end of script ----------"