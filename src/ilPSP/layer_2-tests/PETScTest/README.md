# PETSc managed variant tests

The suite launches each PETSc configuration in a fresh process because PETSc
runtime configurations cannot be exchanged safely inside one process.

Covered scenarios:

- `petsc_mpi`: one rank/one thread and two MPI ranks
- `petsc_omp`: one rank with one and two OpenMP threads
- `petsc_mpi_omp`: two MPI ranks with two OpenMP threads per rank
- expected rejection of `petsc_omp` under two BoSSS MPI ranks
- typed direct LU and Cholesky through the regular API
- the optional MUMPS, SuperLU_DIST and UMFPACK packages on every native
  variant; unavailable packages are explicitly reported as skipped
- the distributed MUMPS and SuperLU_DIST backends under two MPI ranks for the
  MPI and hybrid variants

The executable also supports `--tutorial`, which runs a small SIP Poisson
application with the same PETSc Level-3 configuration. This complements the
matrix smoke tests with an actual BoSSS solver path.

On Linux:

```bash
bash ./run-variant-tests.sh /path/to/BoSSS-native/lib
```

For Slurm, the launcher can include site-specific arguments:

```bash
PETSC_MPI_LAUNCHER="srun --mpi=pmix" \
  bash ./run-variant-tests.sh /path/to/BoSSS-native/lib
```

On Windows, pass the parent `BUILDS` directory containing the three variant
subdirectories:

```bat
run-variant-tests.bat C:\src\BoSSS-native\BUILDS
```

Set `PETSC_TEST_CONFIGURATION` to choose Debug or Release. It defaults to
Release. `PETSC_MPI_LAUNCHER` defaults to `mpiexec` on both platforms.

## Direct-solver API

The options-driven iterative class is named `PETScKrylovSolver`; its Level-3
configuration is `PETScKrylovSolverConfig`. Direct and Krylov implementations
share `PETScSolverBase`, which owns matrix transfer, native lifecycle, solve
validation, convergence statistics, and disposal.

For an ilPSP matrix, use the typed API rather than constructing PETSc option
strings manually:

```csharp
using var solver = new PETScDirectSolver {
    NativeVariant = PETScNativeVariant.Mpi,
    Factorization = PETScDirectFactorization.LU,
    Backend = PETScDirectBackend.Default
};
solver.DefineMatrix(matrix);
SolverResult result = solver.Solve(solution, rightHandSide);
```

Level-3 configurations can use `PETScDirectSolverConfig`, or the convenience
codes `direct_petsc`, `direct_petsc_cholesky`, `direct_petsc_mumps`,
`direct_petsc_superlu_dist`, and `direct_petsc_umfpack`.

The test executable accepts `--backend default|mumps|superlu_dist|umfpack` for
direct solver modes. The test script runs every backend instead of assuming
which configurations a package supports. `--allow-unavailable-direct-backend`
turns only the explicit "backend was not compiled" setup error into a skip;
once a backend is initialized, an unsuccessful solve remains a test failure.
The ilPSP layer rejects the known sequential `default` and `umfpack` backends
before native setup when more than one MPI rank is active. It also rejects
Cholesky with the LU-only UMFPACK and SuperLU_DIST packages.
The variant scripts include negative tests for every one of these restrictions
and require the corresponding exception.

For native MUMPS lifecycle diagnosis, `--probe-mumps-lifecycle` runs MUMPS as
the first and only PETSc solver in a fresh process on a one-row-per-rank
identity matrix. Phase markers distinguish setup, solve, disposal and finalization.
`--probe-krylov-lifecycle` runs the identical matrix and lifecycle through CG
with Jacobi, providing a controlled comparison with the MUMPS path.
Adding `--pin-csr` keeps the managed CSR arrays pinned through native solver
destruction, allowing their lifetime to be tested independently.
The probe also prints its complete one-row CSR ownership data and the relevant
native libraries mapped into rank zero before `KSPSolve`.
