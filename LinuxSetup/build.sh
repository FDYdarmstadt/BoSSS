#!/bin/bash

if [ -z $BUILD_NUMBER ]
then
      echo "BUILD_NUMBER  is empty, setting to 0\n"
      BUILD_NUMBER=0
else
      echo "BUILD_NUMBER  is $BUILD_NUMBER\n"
fi

rm -rf ./BoSSS-Install
rm -rf ./bin
rm -rf ./doc
rm *.run
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/Utils/bcl/bin/Debug/bcl.exe                                            ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/Utils/bcl/bin/Release/bcl.exe                                           ./bin/Release  1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/Utils/AllSpark/bin/Debug/AllSpark.exe                                             ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/Utils/AllSpark/bin/Release/AllSpark.exe                                           ./bin/Release  1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/Utils/MatrixVisualizer/MatrixVisualizerVS17/bin/Debug/MatrixVisualizerVS17.dll    ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/Utils/MatrixVisualizer/MatrixVisualizerVS17/bin/Release/MatrixVisualizerVS17.dll  ./bin/Release  1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/BoSSSpad/bin/Debug/BoSSSpad.exe                   ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/BoSSSpad/bin/Release/BoSSSpad.exe                 ./bin/Release  1
# Include xml files for 'Describe' feature in BoSSSpad
cp ../src/L4-application/BoSSSpad/bin/Debug/*.xml           ./bin/Debug/
cp ../src/L4-application/BoSSSpad/bin/Release/*.xml         ./bin/Release/
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/PlotGenerator/bin/Debug/BoSSS.PlotGenerator.exe   ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/PlotGenerator/bin/Release/BoSSS.PlotGenerator.exe ./bin/Release  1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/NSE_SIMPLE/bin/Debug/NSE_SIMPLE.exe               ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/NSE_SIMPLE/bin/Release/NSE_SIMPLE.exe             ./bin/Release  1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/CNS/bin/Debug/CNS.exe                             ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/CNS/bin/Release/CNS.exe                           ./bin/Release  1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/IBM_Solver/bin/Debug/IBM_Solver.exe               ./bin/Debug    1
mono ../src/Utils/bcl/bin/Debug/bcl.exe deploy-at ../src/L4-application/IBM_Solver/bin/Release/IBM_Solver.exe             ./bin/Release  1
mkdir bin/native
cp -R $BOSSS_INSTALL/bin/native/* ./bin/native
# remove unnecessary directories present on build-machine
rm -fR ./bin/native/win
rm -fR ./bin/native/linux/amd64-openmpi-*
mkdir ./doc
mkdir ./doc/ControlExamples
mkdir ./doc/ControlExamples/CNS
mkdir ./doc/ControlExamples/IBM
mkdir ./doc/APIreference
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
rm -R BoSSS-Install
