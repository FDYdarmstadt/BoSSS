@echo off
setlocal

set "NATIVE_DIR=%~1"
if not defined NATIVE_DIR set "NATIVE_DIR=%BOSSS_NATIVE_DIR%"
if not defined NATIVE_DIR (
    echo Usage: run-variant-tests.bat ^<BoSSS-native-BUILDS-directory^>
    exit /b 2
)
if not exist "%NATIVE_DIR%" (
    echo Native library directory does not exist: "%NATIVE_DIR%"
    exit /b 2
)
if not defined PETSC_TEST_CONFIGURATION set "PETSC_TEST_CONFIGURATION=Release"
if not defined PETSC_MPI_LAUNCHER set "PETSC_MPI_LAUNCHER=mpiexec"

set "PROJECT_DIR=%~dp0"
set "TEST_DLL=%PROJECT_DIR%bin\%PETSC_TEST_CONFIGURATION%\net8.0\PETScTest.dll"
dotnet build "%PROJECT_DIR%PETScTest.csproj" -c "%PETSC_TEST_CONFIGURATION%"
if errorlevel 1 exit /b %errorlevel%

call :run_serial "MPI, one rank and one thread" petsc_mpi 1
if errorlevel 1 exit /b %errorlevel%
call :run_mpi "MPI, two ranks" petsc_mpi 1 0
if errorlevel 1 exit /b %errorlevel%
call :run_serial "OpenMP, one thread" petsc_omp 1
if errorlevel 1 exit /b %errorlevel%
call :run_serial "OpenMP, two threads" petsc_omp 2
if errorlevel 1 exit /b %errorlevel%
call :run_mpi "MPI+OpenMP, two ranks and two threads per rank" petsc_mpi_omp 2 0
if errorlevel 1 exit /b %errorlevel%
call :run_mpi "OpenMP-only rejects two MPI ranks" petsc_omp 1 1
if errorlevel 1 exit /b %errorlevel%
call :run_serial "Direct LU through MPI variant, one rank" petsc_mpi 1 direct_lu
if errorlevel 1 exit /b %errorlevel%
call :run_serial "Direct Cholesky through OpenMP variant" petsc_omp 2 direct_cholesky
if errorlevel 1 exit /b %errorlevel%
call :run_serial "Direct LU through hybrid variant, one rank" petsc_mpi_omp 2 direct_lu
if errorlevel 1 exit /b %errorlevel%

echo All PETSc variant tests passed.
exit /b 0

:run_serial
echo == %~1 ==
set "BOSSS_NATIVE_DIR=%NATIVE_DIR%"
set "BOSSS_PETSC_VARIANT=%~2"
set "BOSSS_PETSC_SOLVER_MODE=%~4"
if not defined BOSSS_PETSC_SOLVER_MODE set "BOSSS_PETSC_SOLVER_MODE=krylov"
set "BOSSS_PETSC_EXPECT_OPENMP_MPI_REJECTION=0"
set "OMP_NUM_THREADS=%~3"
dotnet "%TEST_DLL%"
exit /b %errorlevel%

:run_mpi
echo == %~1 ==
set "BOSSS_NATIVE_DIR=%NATIVE_DIR%"
set "BOSSS_PETSC_VARIANT=%~2"
set "BOSSS_PETSC_SOLVER_MODE=krylov"
set "OMP_NUM_THREADS=%~3"
set "BOSSS_PETSC_EXPECT_OPENMP_MPI_REJECTION=%~4"
%PETSC_MPI_LAUNCHER% -n 2 dotnet "%TEST_DLL%"
exit /b %errorlevel%
