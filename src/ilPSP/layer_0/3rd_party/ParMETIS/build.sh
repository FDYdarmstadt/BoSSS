rm *.so 

cd METIS
rm *.o
mpicc -fPIC -c *.c
cd ..
mpicc -shared -o libmetis.so ./METIS/*.o

cd ParMETIS
rm *.o
mpicc -fPIC -c *.c
#mpicc -shared -o libparmetis.so *.o ../METIS/*.o
cd ..
mpicc -shared -o libparmetis.so ./ParMETIS/*.o -L. -lmetis
