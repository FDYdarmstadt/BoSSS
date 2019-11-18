/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.SipPoisson
{

    /// <summary>
    /// Benchmark application, solves a Poisson problem using the symmetric interior penalty (SIP) method.
    /// </summary>
    public class SipPoissonMain : Application<SipControl> {

#pragma warning disable 649
        /// <summary>
        /// dependent variable
        /// </summary>
        [InstantiateFromControlFile("T", "T", IOListOption.ControlFileDetermined)]
        protected SinglePhaseField T;

        /// <summary>
        /// exact solution, to determine L2-Error, see also <see cref="SipControl.ExactSolution_provided"/>.
        /// </summary>
        [InstantiateFromControlFile("Tex", "Tex", IOListOption.Never)]
        protected SinglePhaseField Tex;

        /// <summary>
        /// dependent variable
        /// </summary>
        [InstantiateFromControlFile("RHS", "T", IOListOption.ControlFileDetermined)]
        protected SinglePhaseField RHS;
#pragma warning restore 649

        /// <summary>
        /// Solver residual, for a DG discretization which is one degree higher than the degree of the solution
        /// </summary>
        private SinglePhaseField ResiualKP1;

        /// <summary>
        /// error of the numerical solution
        /// </summary>
        private SinglePhaseField Error;

        /// <summary>
        /// DG field instantiation
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();

            ResiualKP1 = new SinglePhaseField(new Basis(this.GridData, T.Basis.Degree + 1), "ResidualKP1");
            base.IOFields.Add(ResiualKP1);

            Error = new SinglePhaseField(new Basis(this.GridData, Math.Max(T.Basis.Degree + 1, Tex.Basis.Degree)), "Error");
            base.m_IOFields.Add(Error);
        }

        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            _Main(args, false, delegate () {
                SipPoissonMain p = new SipPoissonMain();
                Console.WriteLine("ipPoisson: " + ilPSP.Environment.MPIEnv.MPI_Rank + " of " + ilPSP.Environment.MPIEnv.MPI_Size
                    + " on compute node '" + ilPSP.Environment.MPIEnv.Hostname + "';");
                return p;
            });
        }
 
        /// <summary>
        /// LHS of the equation <see cref="LaplaceMtx"/>*<see cref="T"/> + <see cref="LaplaceAffine"/> = <see cref="RHS"/>.
        /// </summary>
        BlockMsrMatrix LaplaceMtx;

        /// <summary>
        /// Part of the RHS which contains e.g. boundary conditions; still on LHS, must be subtracted from RHS of the equation.
        /// <see cref="LaplaceMtx"/>*<see cref="T"/> + <see cref="LaplaceAffine"/> = <see cref="RHS"/>
        /// </summary>
        double[] LaplaceAffine;

        /// <summary>
        /// Spatial operator to assemble <see cref="LaplaceMtx"/> and <see cref="LaplaceAffine"/>.
        /// </summary>
        SpatialOperator LapaceIp;

        /// <summary>
        /// Includes assembly of the matrix.
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            using (FuncTrace tr = new FuncTrace()) {

                // create operator
                // ===============
                {
                    double D = this.GridData.SpatialDimension;
                    double penalty_base = (T.Basis.Degree + 1) * (T.Basis.Degree + D) / D;
                    double penalty_factor = base.Control.penalty_poisson;

                    BoundaryCondMap<BoundaryType> PoissonBcMap = new BoundaryCondMap<BoundaryType>(this.GridData, this.Control.BoundaryValues, "T");

                    LapaceIp = new SpatialOperator(1, 1, QuadOrderFunc.SumOfMaxDegrees(), "T", "T");

                    MultidimensionalArray LengthScales;
                    if (this.GridData is GridData) {
                        LengthScales = ((GridData)GridData).Cells.cj;
                    } else if (this.GridData is AggregationGridData) {
                        LengthScales = ((AggregationGridData)GridData).AncestorGrid.Cells.cj;
                    } else {
                        throw new NotImplementedException();
                    }

                    var flux = new ipFlux(penalty_base * base.Control.penalty_poisson, LengthScales, PoissonBcMap);

                    LapaceIp.EquationComponents["T"].Add(flux);

                    LapaceIp.Commit();
                }
            }
        }

        /// <summary>
        /// computes <see cref="LaplaceMtx"/> and <see cref="LaplaceAffine"/>
        /// </summary>
        private void UpdateMatrices() {
            using (var tr = new FuncTrace()) {
                             
                // time measurement for matrix assembly
                Stopwatch stw = new Stopwatch();
                stw.Start();

                // Stats:
                {
                    int BlkSize = T.Mapping.MaxTotalNoOfCoordinatesPerCell;
                    int NoOfMtxBlocks = 0;
                    foreach (int[] Neigs in this.GridData.iLogicalCells.CellNeighbours) {
                        NoOfMtxBlocks++; //               diagonal block
                        NoOfMtxBlocks += Neigs.Length; // off-diagonal block
                    }
                    NoOfMtxBlocks = NoOfMtxBlocks.MPISum();

                    int MtxBlockSize = BlkSize * BlkSize;
                    int MtxSize = MtxBlockSize * NoOfMtxBlocks;

                    double MtxStorage = MtxSize * (8.0 + 4.0) / (1024 * 1024); // 12 bytes (double+int) per entry

                    Console.WriteLine("   System size:                 {0}", T.Mapping.TotalLength);
                    Console.WriteLine("   No of blocks:                {0}", T.Mapping.TotalNoOfBlocks);
                    Console.WriteLine("   No of blocks in matrix:      {0}", NoOfMtxBlocks);
                    Console.WriteLine("   DG coordinates per cell:     {0}", BlkSize);
                    Console.WriteLine("   Non-zeros per matrix block:  {0}", MtxBlockSize);
                    Console.WriteLine("   Total non-zeros in matrix:   {0}", MtxSize);
                    Console.WriteLine("   Approx. matrix storage (MB): {0}", MtxStorage);


                    base.QueryHandler.ValueQuery("MtxBlkSz", MtxBlockSize, true);
                    base.QueryHandler.ValueQuery("NNZMtx", MtxSize, true);
                    base.QueryHandler.ValueQuery("NNZblk", NoOfMtxBlocks, true);
                    base.QueryHandler.ValueQuery("MtxMB", MtxStorage, true);
                }


                Console.WriteLine("creating sparse system for {0} DOF's ...", T.Mapping.Ntotal);


                // quadrature domain
                var volQrSch = new CellQuadratureScheme(true, CellMask.GetFullMask(this.GridData, MaskType.Geometrical));
                var edgQrSch = new EdgeQuadratureScheme(true, EdgeMask.GetFullMask(this.GridData, MaskType.Geometrical));

                using (new BlockTrace("SipMatrixAssembly", tr)) {
                    LaplaceMtx = new BlockMsrMatrix(T.Mapping);
                    LaplaceAffine = new double[T.Mapping.LocalLength];

                    LapaceIp.ComputeMatrixEx(T.Mapping, null, T.Mapping,
                                             LaplaceMtx, LaplaceAffine,
                                             volQuadScheme: volQrSch, edgeQuadScheme: edgQrSch);
                }

                stw.Stop();

                //Print process information
                Console.WriteLine("done {0} sec.", stw.Elapsed.TotalSeconds);
                LaplaceMtx.GetMemoryInfo(out long AllocatedMem, out long UsedMem);
                Console.WriteLine("   Used   matrix storage (MB): {0}", UsedMem /(1024.0*1024));
                Console.WriteLine("   Alloc. matrix storage (MB): {0}", AllocatedMem/(1024.0*1024));
            }
        }

        /// <summary>
        /// control of mesh adaptation
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            if (this.Control.AdaptiveMeshRefinement && TimestepNo > 1) {

                int oldJ = this.GridData.CellPartitioning.TotalLength;

                double LocNormPow2 = this.ResiualKP1.CoordinateVector.L2NormPow2(); // norm of residual on this processor
                double TotNormPow2 = LocNormPow2.MPISum(); //                          norm of residual over all processors
                double MeanNormPow2PerCell = TotNormPow2 / oldJ; //                    mean norm per cell

                double maxSoFar = 0;
                int jMax = -1;
                for (int j = 0; j < oldJ; j++) {
                    double CellNorm = Error.Coordinates.GetRow(j).L2NormPow2();

                    if (CellNorm > maxSoFar) {
                        jMax = j;
                        maxSoFar = CellNorm;
                    }
                }


                int LevelIndicator(int j, int CurrentLevel) {
                    double CellNorm = this.ResiualKP1.Coordinates.GetRow(j).L2NormPow2();

                    if (j == jMax)
                        return CurrentLevel + 1;
                    else
                        return CurrentLevel;
                }

                bool AnyChange = GridRefinementController.ComputeGridChange((GridData)(this.GridData), null, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                int NoOfCellsToRefine = 0;
                int NoOfCellsToCoarsen = 0;
                if (AnyChange) {
                    int[] glb = (new int[] {
                        CellsToRefineList.Count,
                        Coarsening.Sum(L => L.Length),
                    }).MPISum();

                    NoOfCellsToRefine = glb[0];
                    NoOfCellsToCoarsen = glb[1];
                }

                // Update Grid
                // ===========

                if (AnyChange) {
                    Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                    newGrid = ((GridData)(this.GridData)).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                } else {
                    newGrid = null;
                    old2NewGrid = null;
                }
            } else {
                newGrid = null;
                old2NewGrid = null;
            }
        }

        /// <summary>
        /// Single run of the solver
        /// </summary>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {

                if (Control.ExactSolution_provided) {
                    Tex.Clear();
                    Tex.ProjectField(this.Control.InitialValues_Evaluators["Tex"]);

                    RHS.Clear();
                    RHS.ProjectField(this.Control.InitialValues_Evaluators["RHS"]);
                }

                if (Control.AdaptiveMeshRefinement == false) {
                    base.NoOfTimesteps = -1;
                    if (TimestepNo > 1)
                        throw new ApplicationException("steady-state-equation.");
                    base.TerminationKey = true;
                }

                // Update matrices
                // ---------------

                UpdateMatrices();


                // call solver
                // -----------
                double mintime, maxtime;
                bool converged;
                int NoOfIterations;

                LinearSolverCode solvercodes = this.Control.LinearSolver.SolverCode;
                switch (solvercodes) {

                    case LinearSolverCode.classic_cg:
                    case LinearSolverCode.classic_mumps:
                    case LinearSolverCode.classic_pardiso:
                    default:
                        ClassicSolve(out mintime, out maxtime, out converged, out NoOfIterations);
                        break;
                }

                Console.WriteLine("finished; " + NoOfIterations + " iterations.");
                Console.WriteLine("converged? " + converged);
                Console.WriteLine("Timespan: " + mintime + " to " + maxtime + " seconds");


                base.QueryHandler.ValueQuery("minSolRunT", mintime, true);
                base.QueryHandler.ValueQuery("maxSolRunT", maxtime, true);
                base.QueryHandler.ValueQuery("Conv", converged ? 1.0 : 0.0, true);
                base.QueryHandler.ValueQuery("NoIter", NoOfIterations, true);
                base.QueryHandler.ValueQuery("NoOfCells", this.GridData.CellPartitioning.TotalLength, true);
                base.QueryHandler.ValueQuery("DOFs", T.Mapping.TotalLength, true);
                base.QueryHandler.ValueQuery("BlockSize", T.Basis.Length, true);

                

                if (base.Control.ExactSolution_provided) {
                    Error.Clear();
                    Error.AccLaidBack(1.0, Tex);
                    Error.AccLaidBack(-1.0, T);

                    double L2_ERR = Error.L2Norm();
                    Console.WriteLine("\t\tL2 error on " + this.Grid.NumberOfCells + ": " + L2_ERR);
                    base.QueryHandler.ValueQuery("SolL2err", L2_ERR, true);

                }

                // evaluate residual
                // =================
                {
                    this.ResiualKP1.Clear();
                    if (this.Control.InitialValues_Evaluators.ContainsKey("RHS")) {
                        this.ResiualKP1.ProjectField(this.Control.InitialValues_Evaluators["RHS"]);
                    }

                    var ev = this.LapaceIp.GetEvaluator(T.Mapping, ResiualKP1.Mapping);
                    ev.Evaluate(-1.0, 1.0, ResiualKP1.CoordinateVector);
                }

                // return
                // ======

                return 0.0;
            }
        }

        /// <summary>
        /// Solution of the system
        /// <see cref="LaplaceMtx"/>*<see cref="T"/> + <see cref="LaplaceAffine"/> = <see cref="RHS"/>
        /// using a black-box solver
        /// </summary>
        private void ClassicSolve(out double mintime, out double maxtime, out bool Converged, out int NoOfIter) {

            mintime = double.MaxValue;
            maxtime = double.MinValue;
            Converged = false;
            NoOfIter = int.MaxValue;

            for (int i = 0; i < base.Control.NoOfSolverRuns; i++) {

                // create sparse solver
                // --------------------
                ISparseSolver ipSolver;
                LinearSolverCode solvercodes = this.Control.LinearSolver.SolverCode;

                switch (solvercodes) {
                    case LinearSolverCode.classic_pardiso:
                        ipSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver() {
                            CacheFactorization = true,
                            UseDoublePrecision = true
                        };
                        break;

                    case LinearSolverCode.classic_mumps:
                        ipSolver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
                        break;

                    case LinearSolverCode.classic_cg:
                        ipSolver = new ilPSP.LinSolvers.monkey.CG() {
                            MaxIterations = 1000000,
                            Tolerance = 1.0e-10,
                            DevType = ilPSP.LinSolvers.monkey.DeviceType.Cuda
                        };
                        break;

                    default:
                        throw new ArgumentException();
                }
                ipSolver.DefineMatrix(LaplaceMtx);

                // call solver
                // -----------
                T.Clear();

                Console.WriteLine("RUN " + i + ": solving system...");

                var RHSvec = RHS.CoordinateVector.ToArray();
                BLAS.daxpy(RHSvec.Length, -1.0, this.LaplaceAffine, 1, RHSvec, 1);

                T.Clear();
                SolverResult solRes = ipSolver.Solve(T.CoordinateVector, RHSvec);
                mintime = Math.Min(solRes.RunTime.TotalSeconds, mintime);
                maxtime = Math.Max(solRes.RunTime.TotalSeconds, maxtime);

                Converged = solRes.Converged;
                NoOfIter = solRes.NoOfIterations;

                Console.WriteLine("Pardiso phase 11: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_11.Elapsed.TotalSeconds);
                Console.WriteLine("Pardiso phase 22: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_22.Elapsed.TotalSeconds);
                Console.WriteLine("Pardiso phase 33: " + ilPSP.LinSolvers.PARDISO.PARDISOSolver.Phase_33.Elapsed.TotalSeconds);

                ipSolver.Dispose();
            }
        }

        /// <summary>
        /// Shutdown function
        /// </summary>
        protected override void Bye() {
            object SolL2err;
            if (this.QueryHandler.QueryResults.TryGetValue("SolL2err", out SolL2err)) {
                Console.WriteLine("Value of Query 'SolL2err' " + SolL2err.ToString());
            } else {
                Console.WriteLine("query 'SolL2err' not found.");
            }
        }

        /// <summary>
        /// default plotting
        /// </summary>
        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0) {
            DGField[] Fields = new DGField[] { T, Tex, RHS, ResiualKP1, Error };
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, "poisson_" + timestepNo, phystime, superSampling);
        }
    }
}
