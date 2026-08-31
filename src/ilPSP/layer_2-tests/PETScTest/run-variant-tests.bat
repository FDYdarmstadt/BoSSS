@echo off
setlocal

if not defined PETSC_TEST_CONFIGURATION set "PETSC_TEST_CONFIGURATION=Release"
if not defined PETSC_MPI_LAUNCHER set "PETSC_MPI_LAUNCHER=mpiexec"

set "PROJECT_DIR=%~dp0"
set "TEST_DLL=%PROJECT_DIR%bin\%PETSC_TEST_CONFIGURATION%\net8.0\PETScTest.dll"
if not "%PETSC_TEST_SKIP_BUILD%"=="1" (
    dotnet build "%PROJECT_DIR%PETScTest.csproj" -c "%PETSC_TEST_CONFIGURATION%"
    if errorlevel 1 exit /b %errorlevel%
)
dotnet "%TEST_DLL%" --test-auto-selection
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
call :run_serial "Automatic selection, one rank and two threads" auto 2
if errorlevel 1 exit /b %errorlevel%
call :run_mpi "Automatic selection, two ranks and one thread" auto 1 0
if errorlevel 1 exit /b %errorlevel%
call :run_mpi "Automatic selection, two ranks and two threads per rank" auto 2 0
if errorlevel 1 exit /b %errorlevel%
call :run_serial "Direct LU through MPI variant, one rank" petsc_mpi 1 direct_lu
if errorlevel 1 exit /b %errorlevel%
call :run_serial "Direct Cholesky through OpenMP variant" petsc_omp 2 direct_cholesky
if errorlevel 1 exit /b %errorlevel%
call :run_serial "Direct LU through hybrid variant, one rank" petsc_mpi_omp 2 direct_lu
if errorlevel 1 exit /b %errorlevel%
for %%V in (petsc_mpi petsc_omp petsc_mpi_omp) do (
    for %%B in (mumps superlu_dist umfpack) do (
        if /I not "%%V/%%B"=="petsc_omp/superlu_dist" (
            call :run_serial "Optional PETSc %%B package (%%V, one rank)" %%V 1 direct_lu 4 %%B
            if errorlevel 1 exit /b 1
        )
    )
)
for %%B in (mumps superlu_dist) do (
    call :run_mpi "PETSc %%B direct backend (MPI, two ranks)" petsc_mpi 1 0 4 direct_lu %%B
    if errorlevel 1 exit /b 1
    call :run_mpi "PETSc %%B direct backend (MPI+OpenMP, two ranks)" petsc_mpi_omp 2 0 4 direct_lu %%B
    if errorlevel 1 exit /b 1
)
call :run_expected_direct_rejection "Default direct backend rejects two MPI ranks" 2 1 petsc_mpi direct_lu default
if errorlevel 1 exit /b %errorlevel%
call :run_expected_direct_rejection "UMFPACK direct backend rejects two MPI ranks" 2 1 petsc_mpi direct_lu umfpack
if errorlevel 1 exit /b %errorlevel%
call :run_expected_direct_rejection "UMFPACK rejects Cholesky" 1 1 petsc_omp direct_cholesky umfpack
if errorlevel 1 exit /b %errorlevel%
call :run_expected_direct_rejection "SuperLU_DIST rejects Cholesky" 1 1 petsc_omp direct_cholesky superlu_dist
if errorlevel 1 exit /b %errorlevel%
call :run_expected_direct_rejection "SuperLU_DIST requires an MPI-capable variant" 1 1 petsc_omp direct_lu superlu_dist
if errorlevel 1 exit /b %errorlevel%
call :run_mpi "Large MPI matrix, two ranks" petsc_mpi 1 0 2048
if errorlevel 1 exit /b %errorlevel%
call :run_serial "Large OpenMP matrix" petsc_omp 2 krylov 4096
if errorlevel 1 exit /b %errorlevel%
call :run_mpi "Large MPI+OpenMP matrix, two ranks" petsc_mpi_omp 2 0 2048
if errorlevel 1 exit /b %errorlevel%

echo All PETSc variant tests passed.
exit /b 0

:run_serial
echo == %~1 ==
set "SOLVER_MODE=%~4"
if not defined SOLVER_MODE set "SOLVER_MODE=krylov"
set "LOCAL_ROWS=%~5"
if not defined LOCAL_ROWS set "LOCAL_ROWS=4"
set "BACKEND=%~6"
set "BACKEND_ARG="
if defined BACKEND set "BACKEND_ARG=--backend %BACKEND% --allow-unavailable-direct-backend"
set "OMP_NUM_THREADS=%~3"
dotnet "%TEST_DLL%" --variant "%~2" --solver "%SOLVER_MODE%" --local-rows "%LOCAL_ROWS%" %BACKEND_ARG%
exit /b %errorlevel%

:run_mpi
echo == %~1 ==
set "LOCAL_ROWS=%~5"
if not defined LOCAL_ROWS set "LOCAL_ROWS=4"
set "SOLVER_MODE=%~6"
if not defined SOLVER_MODE set "SOLVER_MODE=krylov"
set "BACKEND=%~7"
set "BACKEND_ARG="
if defined BACKEND set "BACKEND_ARG=--backend %BACKEND% --allow-unavailable-direct-backend"
set "OMP_NUM_THREADS=%~3"
set "EXPECT_REJECTION_ARG="
if "%~4"=="1" set "EXPECT_REJECTION_ARG=--expect-openmp-mpi-rejection"
%PETSC_MPI_LAUNCHER% -n 2 dotnet "%TEST_DLL%" --variant "%~2" --solver "%SOLVER_MODE%" --local-rows "%LOCAL_ROWS%" %EXPECT_REJECTION_ARG% %BACKEND_ARG%
exit /b %errorlevel%

:run_expected_direct_rejection
echo == %~1 ==
set "OMP_NUM_THREADS=%~3"
%PETSC_MPI_LAUNCHER% -n %~2 dotnet "%TEST_DLL%" --variant "%~4" --solver "%~5" --backend "%~6" --expect-unsupported-direct-configuration
exit /b %errorlevel%
