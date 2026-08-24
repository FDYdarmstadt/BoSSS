#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
native_dir="${1:-${BOSSS_NATIVE_DIR:-}}"
configuration="${PETSC_TEST_CONFIGURATION:-Release}"

if [[ -z "${native_dir}" ]]; then
    echo "Usage: $0 <BoSSS-native-lib-directory>" >&2
    echo "Alternatively set BOSSS_NATIVE_DIR." >&2
    exit 2
fi
if [[ ! -d "${native_dir}" ]]; then
    echo "Native library directory does not exist: ${native_dir}" >&2
    exit 2
fi

read -r -a mpi_launcher <<< "${PETSC_MPI_LAUNCHER:-mpiexec}"
if [[ ${#mpi_launcher[@]} -eq 0 ]]; then
    echo "PETSC_MPI_LAUNCHER must name an MPI launcher." >&2
    exit 2
fi

project="${script_dir}/PETScTest.csproj"
test_dll="${script_dir}/bin/${configuration}/net8.0/PETScTest.dll"
dotnet build "${project}" -c "${configuration}"

run_serial() {
    local label="$1"
    local variant="$2"
    local threads="$3"
    local solver_mode="${4:-krylov}"
    echo "== ${label} =="
    env \
        BOSSS_NATIVE_DIR="${native_dir}" \
        BOSSS_PETSC_VARIANT="${variant}" \
        BOSSS_PETSC_SOLVER_MODE="${solver_mode}" \
        BOSSS_PETSC_EXPECT_OPENMP_MPI_REJECTION=0 \
        OMP_NUM_THREADS="${threads}" \
        dotnet "${test_dll}"
}

run_mpi() {
    local label="$1"
    local variant="$2"
    local threads="$3"
    local expect_rejection="${4:-0}"
    echo "== ${label} =="
    env \
        BOSSS_NATIVE_DIR="${native_dir}" \
        BOSSS_PETSC_VARIANT="${variant}" \
        BOSSS_PETSC_SOLVER_MODE=krylov \
        BOSSS_PETSC_EXPECT_OPENMP_MPI_REJECTION="${expect_rejection}" \
        OMP_NUM_THREADS="${threads}" \
        "${mpi_launcher[@]}" -n 2 dotnet "${test_dll}"
}

run_serial "MPI, one rank and one thread" petsc_mpi 1
run_mpi "MPI, two ranks" petsc_mpi 1
run_serial "OpenMP, one thread" petsc_omp 1
run_serial "OpenMP, two threads" petsc_omp 2
run_mpi "MPI+OpenMP, two ranks and two threads per rank" petsc_mpi_omp 2
run_mpi "OpenMP-only rejects two MPI ranks" petsc_omp 1 1
run_serial "Direct LU through MPI variant, one rank" petsc_mpi 1 direct_lu
run_serial "Direct Cholesky through OpenMP variant" petsc_omp 2 direct_cholesky
run_serial "Direct LU through hybrid variant, one rank" petsc_mpi_omp 2 direct_lu

echo "All PETSc variant tests passed."
