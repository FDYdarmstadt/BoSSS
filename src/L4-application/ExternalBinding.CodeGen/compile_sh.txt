#ORIGIN_PATH=/mnt/c/Users/florian/Documents/BoSSS-master/public/src/L4-application/ExternalBinding.CodeGen/bin/Debug
#rm -r ./managed/*
#cp $ORIGIN_PATH/{*.dll,*.exe} ./managed/
rm *.o

SourceFiles=$(ls *.cpp)

for SourceFile in $SourceFiles; do
    g++ -c $SourceFile  `pkg-config --cflags --libs mono-2`
done
g++ *.o `pkg-config --cflags --libs mono-2` -o ExtBindingTest.out 