rmdir .\Output /S /Q#
rmdir .\bin /S /Q
rmdir .\doc /S /Q
dotnet publish ..\src\L4-application\PublicTestRunner\PublicTestRunner.csproj --configuration=Release --output .\bin\Release\net8.0
mkdir bin\native
xcopy "%BOSSS_INSTALL%\bin\native" .\bin\native /E /Y
rem Windows prerequisites are staged from the archived BoSSS-native build.
set "WIN_REDIST=installer-tmp_win-redist"
set "MPI_INSTALLER="
for /f "delims=" %%F in ('dir /b /a-d "%WIN_REDIST%\*mpi*" 2^>nul') do if not defined MPI_INSTALLER set "MPI_INSTALLER=%%F"
IF NOT DEFINED MPI_INSTALLER (
    echo ERROR: Missing MPI installer in "%WIN_REDIST%" ^(expected a filename containing "mpi"^).
    exit /b 1
)
set "REDIST_INSTALLER="
for /f "delims=" %%F in ('dir /b /a-d "%WIN_REDIST%\*redist*" 2^>nul') do if not defined REDIST_INSTALLER set "REDIST_INSTALLER=%%F"
IF NOT DEFINED REDIST_INSTALLER (
    echo ERROR: Missing Visual C++ redistributable in "%WIN_REDIST%" ^(expected a filename containing "redist"^).
    exit /b 1
)
echo replacing installer redistributables (Windows)
rmdir .\bin\native\win\redist /S /Q
mkdir .\bin\native\win\redist
xcopy "%WIN_REDIST%\*" .\bin\native\win\redist\ /E /Y
if errorlevel 1 exit /b %ERRORLEVEL%
copy /Y "%WIN_REDIST%\%MPI_INSTALLER%" ".\bin\native\win\redist\msmpisetup10.1.exe"
if errorlevel 1 exit /b %ERRORLEVEL%
copy /Y "%WIN_REDIST%\%REDIST_INSTALLER%" ".\bin\native\win\redist\vcredist_x64.exe"
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
