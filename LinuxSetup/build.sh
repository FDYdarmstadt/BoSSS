#!/bin/bash

if [ -z $BUILD_NUMBER ]
then
      echo "BUILD_NUMBER  is empty, setting to 0"
      BUILD_NUMBER=0
else
      echo "BUILD_NUMBER  is $BUILD_NUMBER"
fi

rm -rf ./BoSSS-Install
rm -rf ./bin
rm -rf ./doc
rm -f *.run
mkdir bin 
dotnet publish ../src/L4-application/PublicTestRunner/PublicTestRunner.csproj --configuration=Release --output ./bin/Release/net6.0
dotnet publish ../src/L4-application/PublicTestRunner/PublicTestRunner.csproj --configuration=Release --output ./bin/Release/net6.0

mkdir bin/native
cp -R $BOSSS_INSTALL/bin/native/* ./bin/native
# remove unnecessary directories present on build-machine
rm -fR ./bin/native/win
rm -fR ./bin/native/linux/amd64-openmpi-*
# patch newest native binaries
if [ -d "installer-tmp_linux-amd64-openmpi" ]; then
    # 
    #echo patching installer native binaries (windows)
    rm -fR ./bin/native/linux/amd64-openmpi/*
    cp -r ./installer-tmp_linux-amd64-openmpi/* ./bin/native/linux/amd64-openmpi/     
fi
mkdir ./doc
mkdir ./doc/ControlExamples
mkdir ./doc/ControlExamples/CNS
mkdir ./doc/ControlExamples/IBM
# mkdir ./doc/APIreference
cp ../doc/handbook/ControlExamples/IBM/* ./doc/ControlExamples/IBM/
cp ../doc/handbook/ControlExamples/CNS/* ./doc/ControlExamples/CNS/
# cp ../doc/notes/0030-BoSSSPad_Command_Overview/BoSSSPad_Command_Overview.pdf ./doc/
# cp ../doc/handbook/BoSSShandbook.pdf ./doc/
# cp ../doc/doxygen/html/* ./doc/APIreference/
# cp ../src/BoSSSpadGUI/InnoSetup/bin/* ./bin/
# Create selfextracting archive using makeself
mkdir BoSSS-Install
cp -R ./bin ./BoSSS-Install
cp -R ./doc ./BoSSS-Install
cp ./Setup.sh ./BoSSS-Install
makeself --notemp ./BoSSS-Install BoSSS-setup-$BUILD_NUMBER.run "BoSSS by Chair of Fluid Dynamics (FDY), TU Darmstadt" echo "BoSSS successfully extracted, please proceed by sourcing Setup.sh"
rm -r ./BoSSS-Install
