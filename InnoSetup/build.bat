rmdir .\Output /S /Q#
rmdir .\bin /S /Q
rmdir .\doc /S /Q
dotnet publish ..\src\L4-application\PublicTestRunner\PublicTestRunner.csproj --configuration=Release --output .\bin\Release\net8.0
mkdir bin\native
xcopy "%BOSSS_INSTALL%\bin\native" .\bin\native /E /Y
set "WIN_REDIST=%BOSSS_NATIVE_WINDOWS_REDIST%"
if not defined WIN_REDIST set "WIN_REDIST=C:\BoSSS-native-windows\WinRedist"
if not exist "%WIN_REDIST%\msmpisetup10.1.exe" (
    echo ERROR: Windows redistributables were not found at "%WIN_REDIST%".
    exit /b 1
)
mkdir ".\bin\native\win\redist"
xcopy "%WIN_REDIST%\*" ".\bin\native\win\redist\" /E /Y
if errorlevel 1 exit /b %ERRORLEVEL%
IF EXIST "installer-tmp_win-amd64" (
    echo patching installer native binaries (windows)
    rmdir .\bin\native\win\amd64\* /S /Q
    xcopy .\installer-tmp_win-amd64 .\bin\native\win\amd64\* /E /Y
) else (
    echo not using new native binaries (Windows)
)
IF EXIST "installer-tmp_linux-amd64-openmpi" (
    echo patching installer native binaries (Linux, amd64-openmpi)
    rmdir .\bin\native\linux\amd64-openmpi /S /Q
    xcopy .\installer-tmp_linux-amd64-openmpi .\bin\native\linux\amd64-openmpi\* /E /Y
) else (
    echo not using new native binaries (Linux, amd64-openmpi)
)
mkdir doc
mkdir doc\ControlExamples
mkdir doc\ControlExamples\CNS
xcopy ..\doc\
xcopy ..\doc\handbook\ControlExamples\CNS\* .\doc\ControlExamples\CNS /Y
xcopy ..\doc\notes\0030-BoSSSPad_Command_Overview\BoSSSPad_Command_Overview.pdf .\doc\ /Y
::xcopy ..\doc\handbook\BoSSShandbook.pdf .\doc\ /Y
::xcopy ..\doc\doxygen\html\* .\doc\APIreference /E /Q
::xcopy ..\src\BoSSSpadGUI\InnoSetup\bin\* .\bin /c /e
ISCC.exe BoSSS-setup.iss
