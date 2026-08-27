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
if [[ "${PETSC_TEST_SKIP_BUILD:-0}" != "1" ]]; then
    dotnet build "${project}" -c "${configuration}"
fi
dotnet "${test_dll}" --test-auto-selection

run_serial() {
    local label="$1"
    local variant="$2"
    local threads="$3"
    local solver_mode="${4:-krylov}"
    local local_rows="${5:-4}"
    local backend="${6:-}"
    local extra_args=()
    if [[ -n "${backend}" ]]; then
        extra_args+=(--backend "${backend}" --allow-unavailable-direct-backend)
    fi
    echo "== ${label} =="
    env \
        OMP_NUM_THREADS="${threads}" \
        dotnet "${test_dll}" \
            --variant "${variant}" \
            --solver "${solver_mode}" \
            --local-rows "${local_rows}" \
            "${extra_args[@]}"
}

run_mpi() {
    local label="$1"
    local variant="$2"
    local threads="$3"
    local expect_rejection="${4:-0}"
    local local_rows="${5:-4}"
    local solver_mode="${6:-krylov}"
    local backend="${7:-}"
    local extra_args=()
    if [[ "${expect_rejection}" == "1" ]]; then
        extra_args+=(--expect-openmp-mpi-rejection)
    fi
    if [[ -n "${backend}" ]]; then
        extra_args+=(--backend "${backend}" --allow-unavailable-direct-backend)
    fi
    echo "== ${label} =="
    env \
        OMP_NUM_THREADS="${threads}" \
        "${mpi_launcher[@]}" -n 2 dotnet "${test_dll}" \
            --variant "${variant}" \
            --solver "${solver_mode}" \
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
for variant in petsc_mpi petsc_omp petsc_mpi_omp; do
    for backend in mumps superlu_dist umfpack; do
        run_serial "Optional PETSc ${backend} package (${variant}, one rank)" "${variant}" 1 direct_lu 4 "${backend}"
    done
done
for backend in default mumps superlu_dist umfpack; do
    run_mpi "PETSc ${backend} direct backend (MPI, two ranks)" petsc_mpi 1 0 4 direct_lu "${backend}"
    run_mpi "PETSc ${backend} direct backend (MPI+OpenMP, two ranks)" petsc_mpi_omp 2 0 4 direct_lu "${backend}"
done
run_mpi "Large MPI matrix, two ranks" petsc_mpi 1 0 2048
run_serial "Large OpenMP matrix" petsc_omp 2 krylov 4096
run_mpi "Large MPI+OpenMP matrix, two ranks" petsc_mpi_omp 2 0 2048
echo "All PETSc variant tests passed."
