# PETSc managed variant tests

The suite launches each PETSc configuration in a fresh process because PETSc
runtime configurations cannot be exchanged safely inside one process.

Covered scenarios:

- `petsc_mpi`: one rank/one thread and two MPI ranks
- `petsc_omp`: one rank with one and two OpenMP threads
- `petsc_mpi_omp`: two MPI ranks with two OpenMP threads per rank
- expected rejection of `petsc_omp` under two BoSSS MPI ranks
- typed direct LU on the MPI and hybrid variants at one rank
- typed direct Cholesky on the OpenMP variant

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

The current portable PETSc libraries contain PETSc's default sequential LU and
Cholesky implementations. MUMPS, SuperLU_DIST, and UMFPACK are exposed by the
managed API but require rebuilding the native PETSc variants with the selected
package. Distributed direct solves should use a PETSc build containing MUMPS or
SuperLU_DIST.
