#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
configuration="${PETSC_TEST_CONFIGURATION:-Release}"

read -r -a mpi_launcher <<< "${PETSC_MPI_LAUNCHER:-mpiexec}"
if [[ ${#mpi_launcher[@]} -eq 0 ]]; then
    echo "PETSC_MPI_LAUNCHER must name an MPI launcher." >&2
    exit 2
fi

project="${script_dir}/PETScTest.csproj"
test_dll="${script_dir}/bin/${configuration}/net8.0/PETScTest.dll"
dotnet build "${project}" -c "${configuration}"
dotnet "${test_dll}" --test-auto-selection

run_serial() {
    local label="$1"
    local variant="$2"
    local threads="$3"
    local solver_mode="${4:-krylov}"
    local local_rows="${5:-4}"
    echo "== ${label} =="
    env \
        OMP_NUM_THREADS="${threads}" \
        dotnet "${test_dll}" \
            --variant "${variant}" \
            --solver "${solver_mode}" \
            --local-rows "${local_rows}"
}

run_mpi() {
    local label="$1"
    local variant="$2"
    local threads="$3"
    local expect_rejection="${4:-0}"
    local local_rows="${5:-4}"
    local extra_args=()
    if [[ "${expect_rejection}" == "1" ]]; then
        extra_args+=(--expect-openmp-mpi-rejection)
    fi
    echo "== ${label} =="
    env \
        OMP_NUM_THREADS="${threads}" \
        "${mpi_launcher[@]}" -n 2 dotnet "${test_dll}" \
            --variant "${variant}" \
            --solver krylov \
            --local-rows "${local_rows}" \
            "${extra_args[@]}"
}

run_serial "MPI, one rank and one thread" petsc_mpi 1
run_mpi "MPI, two ranks" petsc_mpi 1
run_serial "OpenMP, one thread" petsc_omp 1
run_serial "OpenMP, two threads" petsc_omp 2
run_mpi "MPI+OpenMP, two ranks and two threads per rank" petsc_mpi_omp 2
run_mpi "OpenMP-only rejects two MPI ranks" petsc_omp 1 1
run_serial "Automatic selection, one rank and two threads" auto 2
run_mpi "Automatic selection, two ranks and one thread" auto 1
run_mpi "Automatic selection, two ranks and two threads per rank" auto 2
run_serial "Direct LU through MPI variant, one rank" petsc_mpi 1 direct_lu
run_serial "Direct Cholesky through OpenMP variant" petsc_omp 2 direct_cholesky
run_serial "Direct LU through hybrid variant, one rank" petsc_mpi_omp 2 direct_lu
run_mpi "Large MPI matrix, two ranks" petsc_mpi 1 0 2048
run_serial "Large OpenMP matrix" petsc_omp 2 krylov 4096
run_mpi "Large MPI+OpenMP matrix, two ranks" petsc_mpi_omp 2 0 2048

echo "All PETSc variant tests passed."
