rmdir .\Output /S /Q
rmdir .\bin /S /Q
rmdir .\doc /S /Q
bcl deploy-at ..\src\Utils\bcl\bin\Debug\bcl.exe                                             .\bin\Debug    1
bcl deploy-at ..\src\Utils\bcl\bin\Release\bcl.exe                                           .\bin\Release  1
bcl deploy-at ..\src\Utils\AllSpark\bin\Debug\AllSpark.exe                                             .\bin\Debug    1
bcl deploy-at ..\src\Utils\AllSpark\bin\Release\AllSpark.exe                                           .\bin\Release  1
bcl deploy-at ..\src\Utils\MatrixVisualizer\MatrixVisualizerVS15\bin\Debug\MatrixVisualizerVS15.dll    .\bin\Debug    1
bcl deploy-at ..\src\Utils\MatrixVisualizer\MatrixVisualizerVS15\bin\Release\MatrixVisualizerVS15.dll  .\bin\Release  1
bcl deploy-at ..\src\L4-application\BoSSSpad\bin\Debug\BoSSSpad.exe                   .\bin\Debug    1
bcl deploy-at ..\src\L4-application\BoSSSpad\bin\Release\BoSSSpad.exe                 .\bin\Release  1
::Include xml files for 'Describe' feature in BoSSSpad
xcopy ..\src\L4-application\BoSSSpad\bin\Debug\*.xml           .\bin\Debug /Y
xcopy ..\src\L4-application\BoSSSpad\bin\Release\*.xml         .\bin\Release /Y
bcl deploy-at ..\src\L4-application\PlotGenerator\bin\Debug\BoSSS.PlotGenerator.exe   .\bin\Debug    1
bcl deploy-at ..\src\L4-application\PlotGenerator\bin\Release\BoSSS.PlotGenerator.exe .\bin\Release  1
bcl deploy-at ..\src\L4-application\NSE_SIMPLE\bin\Debug\NSE_SIMPLE.exe               .\bin\Debug    1
bcl deploy-at ..\src\L4-application\NSE_SIMPLE\bin\Release\NSE_SIMPLE.exe             .\bin\Release  1
bcl deploy-at ..\src\L4-application\CNS\bin\Debug\CNS.exe                             .\bin\Debug    1
bcl deploy-at ..\src\L4-application\CNS\bin\Release\CNS.exe                           .\bin\Release  1
bcl deploy-at ..\src\L4-application\IBM_Solver\bin\Debug\IBM_Solver.exe               .\bin\Debug    1
bcl deploy-at ..\src\L4-application\IBM_Solver\bin\Release\IBM_Solver.exe             .\bin\Release  1
mkdir bin\native
xcopy "%BOSSS_INSTALL%\bin\native" .\bin\native /E /Y
mkdir doc
mkdir doc\ControlExamples
mkdir doc\ControlExamples\CNS
mkdir doc\ControlExamples\IBM
mkdir doc\APIreference
xcopy ..\doc\
xcopy ..\doc\handbook\ControlExamples\IBM\* .\doc\ControlExamples\IBM /Y
xcopy ..\doc\handbook\ControlExamples\CNS\* .\doc\ControlExamples\CNS /Y
xcopy ..\doc\notes\0030-BoSSSPad_Command_Overview\BoSSSPad_Command_Overview.pdf .\doc\ /Y
xcopy ..\doc\handbook\BoSSShandbook.pdf .\doc\ /Y
xcopy ..\doc\doxygen\html\* .\doc\APIreference /E /Q
xcopy ..\src\BoSSSpadGUI\InnoSetup\bin\* .\bin /e
ISCC.exe BoSSS-setup.iss