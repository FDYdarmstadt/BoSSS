using MPI.Wrappers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Connectors;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Application.CahnHilliard;
using static BoSSS.Application.CahnHilliard.CahnHilliardMain;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.Grid.RefElements;

// TODO: test mit norm, jump norm zum laufen bringen
// erst dotnet, dann mono, dann openfoam

namespace BoSSS.Application.ExternalBinding {
    

    public class FixedOperators : IForeignLanguageProxy {

        /// <summary>
        /// 
        /// </summary>
        [CodeGenExport]
        public FixedOperators() {

        }


        IntPtr m_ForeignPtr;

        /// <summary>
        /// %
        /// </summary>
        public void _SetForeignPointer(IntPtr ptr) {
            if (ptr == IntPtr.Zero) {
                m_ForeignPtr = IntPtr.Zero;
            } else {

                if (m_ForeignPtr != IntPtr.Zero) {
                    throw new ApplicationException("already registered");
                }
                m_ForeignPtr = ptr;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public IntPtr _GetForeignPointer() {
            return m_ForeignPtr;
        }


        SinglePhaseField c;
        SinglePhaseField[] velocity;
        SinglePhaseField mu;
        OpenFoamMatrix _mtx;
        OpenFoamPatchField _ptch;
        double rho0 = 1;
        double rho1 = 2;
        ScalarFunction InitFunc()
        {
            // double rMin = 2.0e-3 / sqrt(noOfTotalCells) * 3.0 / sqrt(2);
            // double radius = 0.5e-3;
            double radius = 7;
            // double radius = rMin * 1.3;
            return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x, 2) + Math.Pow(z, 2)) + Math.Pow(radius, 1)) * Math.Sqrt(2)))).Vectorize();
        }

        /// <summary>
        /// Returns the auxiliary field appearing during the solution of the Cahn-Hilliard system
        /// </summary>
        [CodeGenExport]
        public OpenFoamDGField GetFlux() {
            int D = 3;
            OpenFoamDGField Flux = new OpenFoamDGField(_ptch.Grid, 0, D);
            int nCells = _ptch.Grid.NumberOfCells;

            for (int j = 0; j < nCells; j++) {
                double cMean = c.GetMeanValue(j);
                double C0 = (cMean + 1.0)/2.0;
                double C1 = 1.0 - C0;
                for (int d = 0; d < D; d++) {
                    double massFlux = (rho0 * C0 + rho1 * C1) / 2.0 * velocity[d].GetMeanValue(j);
                    if (d == 1 && velocity[d].GetMeanValue(j) < 1e-10)
                        Console.WriteLine("Hello from BoSSS " + velocity[d].GetMeanValue(j));
                    Flux.SetDGcoordinate(d, j, 0, massFlux);
                }
            }
            return Flux;
        }

        /// <summary>
        /// Returns the auxiliary field appearing during the solution of the Cahn-Hilliard system
        /// </summary>
        [CodeGenExport]
        public OpenFoamDGField GetMu() {
            int D = 3;
            OpenFoamDGField Mu0 = new OpenFoamDGField(_ptch.Grid, 0, D);
            int nCells = _ptch.Grid.NumberOfCells;

            for (int j = 0; j < nCells; j++) {
                Mu0.SetDGcoordinate(0, j, 0, mu.GetMeanValue(j));
            }
            return Mu0;
        }

        /// <summary>
        /// Return the norm of the field c
        /// </summary>
        public double Norm(OpenFoamMatrix mtx = null) {
            if (mtx != null) {
                OpenFoamDGField C = mtx.Fields[0];
                c = C.Fields[0] as SinglePhaseField;
                if (InitFunc() != null){
                    c.Clear();
                    c.ProjectField(InitFunc());
                }
            }
            return c.L2Norm();
        }
        /// <summary>
        /// Return the relative change of the norm of the field c
        /// </summary>
        public double NormRelChange() {
            double postNorm = c.L2Norm();
            SinglePhaseField c0 = new SinglePhaseField(c.Basis);
            c0.Clear();
            c0.ProjectField(InitFunc());
            double preNorm = c0.L2Norm();
            return (postNorm - preNorm)/postNorm;
        }

        /// <summary>
        /// Return the jump norm of the field c
        /// </summary>
        public double JumpNorm() {
            GridData grd = (GridData)c.GridDat;
            int D = grd.SpatialDimension;

            c.MPIExchange();

            double Unorm = 0;

            EdgeQuadrature.GetQuadrature(
                new int[] { D + 1 }, grd,
                (new EdgeQuadratureScheme()).Compile(grd, c.Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for (int j = 0; j < Length; j++) {
                        int iEdge = j + i0;
                        int jCell_IN = grd.Edges.CellIndices[iEdge, 0];
                        int jCell_OT = grd.Edges.CellIndices[iEdge, 1];
                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });

                        if (jCell_OT >= 0) {

                            int iTrafo_IN = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                            int iTrafo_OT = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];

                            MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                            MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);

                            NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN, false);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT, false);

                            c.Evaluate(jCell_IN, 1, NS_IN, uIN);
                            c.Evaluate(jCell_OT, 1, NS_OT, uOT);

                            uDiff.Acc(+1.0, uIN);
                            uDiff.Acc(-1.0, uOT);
                        } else {
                            uDiff.Clear();
                        }
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            Unorm = Unorm.MPISum();

            return Unorm.Sqrt();
        }

        /// <summary>
        /// Solves the Cahn-Hilliard equation
        /// This method only contains arguments that can be made available to OpenFOAM given the limitations of the mono-C-interface.
        /// For usage exclusively within BoSSS, the method <see cref="CahnHilliardInternal"/> is generally more convenient.
        /// </summary>
        [CodeGenExport]
        public void CahnHilliard(OpenFoamMatrix mtx, OpenFoamDGField U, OpenFoamPatchField ptch, OpenFoamPatchField ptchU) {

            // TODO sync from OpenFOAM
            // double epsilon = 1e-5; // capillary width
            //                        // double cahn = 0.005;
            //                        // double M = 1; // mobility parameter
            // double M = Math.Sqrt(epsilon); // mobility parameter
            // double sigma = 0.063;
            // double lam = 3 / (2 * Math.Sqrt(2)) * sigma * epsilon; // Holger's lambda
                                                              // double diff = M * lam;
            CahnHilliardParameters chParams = new CahnHilliardParameters(_stationary: true);
            CahnHilliardInternal(mtx, U, ptch, ptchU, null, chParams);
        }

        /// <summary>
        /// Solves the Cahn-Hilliard equation
        /// This method also contains arguments that cannot be made available to OpenFOAM due to limitations of the mono-C-interface.
        /// </summary>
        public void CahnHilliardInternal(OpenFoamMatrix mtx, OpenFoamDGField U, OpenFoamPatchField ptch, OpenFoamPatchField ptchU, ScalarFunction func = null, CahnHilliardParameters chParams = new CahnHilliardParameters()) {
            try {
                // {

                // System.Diagnostics.Debugger.Launch();
                // Console.WriteLine("Debugger is attached; press enter to continue");
                // Console.ReadLine();


                _mtx = mtx;
                _ptch = ptch;


                double tanh(double x)
                {
                    return Math.Tanh(x);
                }
                double pow(double x, int e)
                { // inefficient but who cares
                    // return Math.Tanh(x);
                    double acc = 1.0;
                    for (int i = 0; i < e; i++)
                        acc *= x;
                    return acc;
                }
                double pow2(double x)
                {
                    // return Math.Tanh(x);
                    return pow(x, 2);
                }
                double sqrt(double x)
                {
                    return Math.Sqrt(x);
                    // double Sqrt = x / 2.0;
                    // double temp = 0;
                    // while (Sqrt - temp > 1e-10)
                    // {
                    //     temp = Sqrt;
                    //     Sqrt = (x / temp + temp) / 2;
                    // }
                    // return Sqrt;
                }
                // grid, etc
                // =========

                GridData grd = mtx.ColMap.GridDat as GridData;
                // PlotGrid("grid.plt", grd);

                Basis b = mtx.ColMap.BasisS[0];
                Basis bMu = mtx.ColMap.BasisS[0];
                int noOfTotalCells = b.GridDat.Grid.NumberOfCells;

                OpenFoamDGField C = mtx.Fields[0];
                var u = U.Fields[0] as SinglePhaseField;
                var v = U.Fields[1] as SinglePhaseField;
                var w = U.Fields[2] as SinglePhaseField;
                velocity = new[] { u, v, w };
                // u.ProjectField(((_3D)((x, y, z) => 0.05)).Vectorize());
                // v.ProjectField(((_3D)((x, y, z) => 0.0)).Vectorize());
                // w.ProjectField(((_3D)((x, y, z) => 0.0)).Vectorize());
                // OpenFOAMGrid ofGrid = ptch.Grid;
                // OpenFoamDGField fields = new(ofGrid, b.Degree, 2);
                // OpenFoamMatrix fullMtx = new(ofGrid, fields);

                var map = new UnsetteledCoordinateMapping(b, bMu);

                int nParams = 7;
                // var op = new SpatialOperator(2, nParams, 2, QuadOrderFunc.Linear(), "c", "mu", "VelocityX", "VelocityY", "VelocityZ", "c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_mu");
                var op = new SpatialOperator(2, nParams, 2, QuadOrderFunc.NonLinearWithoutParameters(3), "c", "mu", "VelocityX", "VelocityY", "VelocityZ", "c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_mu");
                // var op = new XSpatialOperatorMk2(2, nParams, 2, QuadOrderFunc.Linear(), new List<string>{"a", "b"}, "c", "mu", "VelocityX", "VelocityY", "VelocityZ","c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_mu");
                // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "mu", "c0", "VelocityX", "VelocityY", "VelocityZ", "c_Res", "mu_Res");
                // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "mu", "c0","LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "c_Res", "mu_Res");


                // op.LinearizationHint = LinearizationHint.GetJacobiOperator;

                double lambda = 0.0;
                double penalty_const = 2.6;
                double diff = chParams.Diffusion;
                double cahn = chParams.Cahn;

                var RealLevSet = new LevelSet(b, "Levset");
                //var RealTracker = new LevelSetTracker((GridData)(b.GridDat), XQuadFactoryHelper.MomentFittingVariants.Saye, 2, new string[] { "a", "b" }, RealLevSet);
                //RealTracker.UpdateTracker(0.0);
                // RealLevSet.Acc(1.0, C);
                // RealLevSet

                int J = b.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                // for (int j = 0; j < J; j++)
                // {
                //     int N = b.GetLength(j);
                //     for (int n = 0; n < N; n++)
                //         RealLevSet.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                // }

                // TODO

                c = C.Fields[0] as SinglePhaseField;
                // ScalarFunction func() {
                //     double rMin = 2.0e-3 / sqrt(noOfTotalCells) * 3.0 / sqrt(2);
                //     double radius = 0.5e-3;
                //     // double radius = rMin * 1.3;
                //     return ((_3D)((x, y, z) => tanh((-sqrt(pow(x - 1.0e-3, 2) + pow(z - 0.0e-3, 2)) + pow(radius, 1)) * 50000))).Vectorize();
                //     // return ((_3D)((x, y, z) => tanh(((x - 0.0011) + 0.01 * z)*3750))).Vectorize();
                //     // return ((_3D)((x, y, z) => Math.Tanh(((x - 0.0011) + 0.01 * z) * 5500))).Vectorize();
                //     // return ((_3D)((x, y, z) => tanh(((x - 2.5) + 0.1 * (y - 2.5))/1))).Vectorize();
                // }
                if (c.L2Norm() < 1e-20){
                    if (func == null){
                        Console.WriteLine("No initialization function given - using droplet");
                        func = InitFunc();
                    }
                    Console.WriteLine("Zero order parameter field encountered - initializing with given function");
                    c.Clear();
                    u.Clear();
                    v.Clear();
                    w.Clear();
                    c.ProjectField(func);
                }
                mu = new SinglePhaseField(b);
                // for (int j = 0; j < J; j++)
                // {
                //     int N = b.GetLength(j);
                //     for (int n = 0; n < N; n++) {
                //         c.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                //         // mu.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                //     }
                // }

                // foreach (var val in new double[]{0, 0.5, 1, 1.5, 1.99, 2.01 ,10})
                //     Console.WriteLine("Correct Value: " + Math.Tanh(val) + "My value: " + tanh(val));

                // // var LSGrad = new MultidimensionalArray(noOfTotalCells);
                // var nsma = new MultidimensionalArray(2);
                // nsma.Allocate(noOfTotalCells, 3);
                // var ns = new NodeSet(Cube.Instance, nsma, false);
                // var LSGrad = new MultidimensionalArray(3);
                // LSGrad.Allocate(noOfTotalCells, ns.NoOfNodes, 3);
                // RealLevSet.EvaluateGradient(0, noOfTotalCells, ns, LSGrad); <- This causes a mono exception
                // var LSGradX = new SinglePhaseField(b);
                // var LSGradY = new SinglePhaseField(b);
                // var LSGradZ = new SinglePhaseField(b);
                // for (int i = 0; i < noOfTotalCells; i++)
                // {
                //     LSGradX.SetMeanValue(i, LSGrad[i, 0, 0]);
                //     LSGradY.SetMeanValue(i, LSGrad[i, 0, 1]);
                //     LSGradZ.SetMeanValue(i, LSGrad[i, 0, 2]);
                //     // LSGradX.SetMeanValue(i, LSGrad[0,i]);
                //     // LSGradY.SetMeanValue(i, LSGrad[1,i]);
                //     // LSGradZ.SetMeanValue(i, LSGrad[2,i]);
                // }
                // RealLevSet.EvaluateGradient
                var domfields = (IReadOnlyDictionary<string, DGField>)(new Dictionary<string, DGField>() { { "c", c }, { "mu", mu } });
                var paramfields = (IReadOnlyDictionary<string, DGField>)(new Dictionary<string, DGField>(){
                        {"VelocityX", u},
                        {"VelocityY", v},
                        {"VelocityZ", w},
                        {"c0", c},
                        // {"LevelSetGradient[0]", LSGradX},
                        // {"LevelSetGradient[1]", LSGradY},
                        // {"LevelSetGradient[2]", LSGradZ}
                        {"LevelSetGradient[0]", c},
                        {"LevelSetGradient[1]", c},
                        {"LevelSetGradient[2]", c}
                        });
                Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> GetNamedInputFields = delegate (DGField[] fields)
                {

                    return (domfields, paramfields);
                };
                RealLevSet.Clear();
                RealLevSet.Acc(1.0, c);
                LevelSetUpdater lsu = new LevelSetUpdater(grd, XQuadFactoryHelper.MomentFittingVariants.Classic,
                                                         6, new string[] { "a", "b" },
                                                         GetNamedInputFields,
                                                         RealLevSet, "c", ContinuityProjectionOption.None);

                var RealTracker = lsu.Tracker;
                // mu.Laplacian(-cahn, c);
                // mu.Acc(-1.0, c);
                // mu.ProjectPow(1.0, c, 3.0);
                RealLevSet.Clear();
                RealLevSet.Acc(1.0, c);
                RealTracker.UpdateTracker(0);

                SubGrid subgr = RealTracker.Regions.GetNearFieldSubgrid(6);
                SubGridBoundaryModes subgrbnd = 0;
                CellMask subgrMask = subgr.VolumeMask;
                CellMask fullMask = CellMask.GetFullMask(grd);
                SubGrid fullSubGrd = new SubGrid(fullMask);

                CellMask mask = null;
                SubGrid sgrid = null;
                // CellMask mask = fullMask;
                // SubGrid sgrid = fullSubGrd;
                // CellMask mask = subgrMask;
                // SubGrid sgrid = subgr;

                int noOfNarrowBandCells = 0;
                if (mask != null) {
                    foreach (bool cellInNarrowBand in mask.GetBitMask())
                    {
                        if (cellInNarrowBand)
                        {
                            noOfNarrowBandCells++;
                        }
                    }
                    if (noOfNarrowBandCells == 0){
                        Console.WriteLine("Solving only in a narrow band containing " + noOfNarrowBandCells + " of " + noOfTotalCells + " cells");
                        // throw new ApplicationException("No interface found");
                        // mask = fullMask;
                        // sgrid = fullSubGrd;
                        // Console.WriteLine("No narrow band cells detected, solving on the whole domain");
                    } else {
                        Console.WriteLine("Solving only in a narrow band containing " + noOfNarrowBandCells + " of " + noOfTotalCells + " cells");
                    }
                }

                System.Collections.BitArray subGridCellMask = mask?.GetBitMask();

                var CHCdiff = new CahnHilliardCDiff(ptch, penalty_const, diff, lambda, mask);
                var CHCconv = new CahnHilliardCConv(new[] { u, v, w }, mask);
                var CHMudiff = new CahnHilliardMuDiff(ptch, penalty_const, cahn, mask);
                var CHMusource = new CahnHilliardMuSource(mask);
                op.EquationComponents["Res_c"].Add(CHCdiff);
                op.EquationComponents["Res_c"].Add(CHCconv);
                op.EquationComponents["Res_mu"].Add(CHMusource);
                op.EquationComponents["Res_mu"].Add(CHMudiff);
                // op.LinearizationHint = LinearizationHint.GetJacobiOperator;

                double[] MassScales = new double[2];
                MassScales[0] = 1.0;
                // MassScales[1] = 1.0;
                op.TemporalOperator = new ConstantTemporalOperator(op, MassScales);

                op.LinearizationHint = LinearizationHint.GetJacobiOperator;

                op.Commit();

                SinglePhaseField Res_c = new SinglePhaseField(b);
                SinglePhaseField Res_mu = new SinglePhaseField(b);
                List<DGField> ParameterMap = new List<DGField>();
                for (int i = 0; i < nParams; i++)
                {
                    ParameterMap.Add(new SinglePhaseField(b));
                }
                for (int j = 0; j < J; j++)
                {
                    int N = b.GetLength(j);
                    for (int n = 0; n < N; n++)
                        ParameterMap[3].Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                }
                for (int j = 0; j < 3; j++)
                {
                    ParameterMap[j] = new[] { u, v, w }[j];
                }
                NonLinearSolverConfig nls = new NonLinearSolverConfig();
                var ls = new DirectSolver.Config();
                nls.SolverCode = NonLinearSolverCode.Newton;
                // nls.SolverCode = NonLinearSolverCode.Picard;
                // nls.ConvergenceCriterion = 1e-5;
                // nls.MaxSolverIterations = 100;
                nls.verbose = true;

                lsu.InitializeParameters(domfields, paramfields);

                // var tp = new Tecplot(grd.Grid.GridData, 3);
                Tecplot("plot.1", 0.0, 3, c, mu, RealLevSet, u, v, w);

                // TODO saye instead of hmf
                // XdgSubGridTimestepping TimeStepper = new XdgSubGridTimestepping(op,
                //                                          new SinglePhaseField[] { c, mu },
                //                                          new SinglePhaseField[] { Res_c, Res_mu },
                //                                          // TimeSteppingScheme.ExplicitEuler,
                //                                          TimeSteppingScheme.ImplicitEuler,
                //                                          sgrid,
                //                                          subgrbnd,
                //                                          LinearSolver: ls,
                //                                          NonLinearSolver: nls,
                //                                          _UpdateLevelset: (() => lsu),
                //                                          _LevelSetHandling: LevelSetHandling.LieSplitting,
                //                                          // _LevelSetHandling: LevelSetHandling.Coupled_Once,
                //                                          _AgglomerationThreshold: 0.0,
                //                                          _optTracker: RealTracker
                //                                          );

                // XdgTimestepping TimeStepperNoSG = new XdgTimestepping(opNoSG,
                //                                          new SinglePhaseField[]{cNoSG, muNoSG},
                //                                          new SinglePhaseField[]{Res_cNoSG, Res_muNoSG},
                //                                          // TimeSteppingScheme.ExplicitEuler,
                //                                          TimeSteppingScheme.ImplicitEuler,
                //                                          LinearSolver: ls,
                //                                          NonLinearSolver: nls,
                //                                          // _UpdateLevelset: (() => lsu),
                //                                          // _LevelSetHandling: LevelSetHandling.LieSplitting,
                //                                          // _LevelSetHandling: LevelSetHandling.Coupled_Once,
                //                                          _AgglomerationThreshold: 0.0,
                //                                          _optTracker: RealTrackerNoSG
                //                                          );

                XdgTimestepping TimeStepper = new XdgTimestepping(op,
                                                         new SinglePhaseField[]{c, mu},
                                                         new SinglePhaseField[]{Res_c, Res_mu},
                                                         // TimeSteppingScheme.CrankNicolson,
                                                         TimeSteppingScheme.ImplicitEuler,
                                                                  null,
                                                                  null,
                                                         LinearSolver: ls,
                                                         NonLinearSolver: nls,
                                                         // RealTracker,
                                                         lsu: (() => lsu)
                                                         // _LevelSetHandling: LevelSetHandling.LieSplitting,
                                                         // _LevelSetHandling: LevelSetHandling.Coupled_Once,
                                                         // _AgglomerationThreshold: 0.0
                                                         );

                double endTime = chParams.endT;
                double dt = chParams.dt;
                double time = 0.0;
                int t = 0;
                while (time < endTime)
                {
                    RealLevSet.Clear();
                    RealLevSet.Acc(1.0, c);
                    RealTracker.UpdateTracker(time);
                    TimeStepper.Solve(time, dt);

                    Tecplot("plot." + (t + 2), time, 3, c, mu, RealLevSet, u, v, w);

                    time += dt;
                    t++;
                    // if (time > 2e0)
                    //     dt *= 10;
                    // Tecplot("plot." + (t + 2), (t + 1) / timesteps, 3, c, mu, RealLevSet, u, v, w, cNoSG, muNoSG);

                }
            } catch (Exception e) {
                Console.WriteLine(e.GetType());
                Console.WriteLine(e.Message);
                Console.WriteLine(e.StackTrace);
                Console.WriteLine(e);
                throw e;
                }
            // }
        }

        /// <summary>
        /// 
        /// </summary>
        [CodeGenExport]
        // public void Laplacian(OpenFoamMatrix mtx) {
        public void Laplacian(OpenFoamMatrix mtx, OpenFoamPatchField ptch) {

            // grid, etc
            // =========

            GridData grd = mtx.ColMap.GridDat as GridData;
            // PlotGrid("grid.plt", grd);

            var b = mtx.ColMap.BasisS[0];
            var map = new UnsetteledCoordinateMapping(b);

            var L = new Laplace(1.3, ptch);
            var op = new SpatialOperator(1, 0, 1, QuadOrderFunc.Linear(), "T", "c1");

            op.EquationComponents["c1"].Add(L);
            op.Commit();

            // evaluate operator
            // =================

            var eval = op.GetMatrixBuilder(map, null, map);
            eval.ComputeMatrix(mtx, mtx.RHSbuffer);
            mtx.RHSbuffer.ScaleV(-1); // convert LHS affine vector to RHS

            Console.WriteLine("Computed Laplacian Matrix, norm is " + mtx.InfNorm());
            // Console.WriteLine("Computed Laplacian Matrix, RHS is ");
            // foreach (var elem in mtx.RHSbuffer)
            //     Console.Write(elem + " ");
            // Console.WriteLine();
        }

        static public void PlotGrid(string filename, IGridData grd) {


            string SanitizeName(string s) {
                char[] ot = s.ToCharArray();
                for(int k = 0; k < ot.Length; k++) {
                    if(char.IsWhiteSpace(ot[k])) {
                        ot[k] = '_';
                    }

                    if (ot[k] == '(')
                        ot[k] = 'L';
                    if (ot[k] == ')')
                        ot[k] = 'R';
                }
                return new string(ot);
            }



            var et2Name = grd.EdgeTagNames;
            Console.WriteLine($"Grid containing {et2Name.Count} EdgeTag names: ");
            int i = 0;
            foreach (var t in et2Name) { // loop over all different edge tag names...
                string name = t.Value;
                byte tag2color = t.Key;

                string sname = SanitizeName(name);

                if (name.Equals(sname)) {
                    Console.WriteLine($"   {i}: {name} -- tag = {tag2color}");
                } else {
                    Console.WriteLine($"   {i}: {name} -- tag = {tag2color}   (marked as '{sname}') in output file.");
                }
                i++;
            }


            var B0 = new Basis(grd, 0);
            SinglePhaseField[] bndyMarkers = new SinglePhaseField[et2Name.Count + 1];

            int[,] Edge2GeomCell = grd.iGeomEdges.CellIndices;
            int[] G2L = grd.iGeomCells.GeomCell2LogicalCell;
            byte[] EdgeTags = grd.iGeomEdges.EdgeTags;

            i = 0;
            foreach(var t in et2Name) { // loop over all different edge tag names...
                string name = t.Value;
                byte tag2color = t.Key;
                string sname = SanitizeName(name);


                var FI = new SinglePhaseField(B0, "Marker-" + sname);
                bndyMarkers[i] = FI;
                i++;

                for (int e = 0; e < EdgeTags.Length; e++) { // loop over edges...
                    byte tag_e = EdgeTags[e];

                    if(tag_e == tag2color) {
                        // mar cells next to edge e

                        foreach(int jG in Edge2GeomCell.GetRow(e)) {
                            if (jG < 0)
                                continue;

                            // convert geometrical cell index to logical cell index
                            int jL;
                            if (G2L == null)
                                jL = jG;
                            else
                                jL = G2L[jG];

                            // color respective cell
                            FI.SetMeanValue(jL, tag2color);
                        }

                    }
                }

            }

            var dummy = new SinglePhaseField(B0, "DummyData");
            bndyMarkers[bndyMarkers.Length - 1] = dummy;

            Tecplot(filename, 0.0, 0, bndyMarkers);
        }

        static public void Tecplot(string filename, double time, int supersampling, params BoSSS.Foundation.DGField[] flds) {
            if (flds == null || flds.Length <= 0) {
                Console.WriteLine("No DG fields specified - not writing any output files.");
                return;
            }

            if (supersampling > 3) {
                Console.WriteLine("Plotting with a supersampling greater than 3 is deactivated because it would very likely exceed this machines memory.");
                Console.WriteLine("Higher supersampling values are supported by external plot application.");
                supersampling = 3;
            }

            string directory = "~/";
            string FullPath;
            // if (directory == null || directory.Length <= 0) {
            //     directory = CurrentDocDir ?? "";
            //     FullPath = Path.Combine(directory, filename);
            // } else {
            FullPath = filename;
            // }

            Console.WriteLine("Writing output file {0}...", FullPath);


            BoSSS.Solution.Tecplot.Tecplot.PlotFields(flds, FullPath, time, supersampling);

        }

        /// <summary>
        /// SIP-form for the Laplacian
        /// </summary>
        class Laplace : BoSSS.Solution.NSECommon.SIPLaplace {

            OpenFoamPatchField _ptch;

            /// <summary>
            /// 
            /// </summary>
            public Laplace(double penalty_const, OpenFoamPatchField ptch)
                : base(penalty_const, "T")
            {
                _ptch = ptch;
            }

            /// <summary>
            /// always true
            /// </summary>
            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                return _ptch.IsDirichlet(inp.EdgeTag);
            }

            /// <summary>
            /// Dirichlet boundary value
            /// </summary>
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) {
                if (inp.EdgeTag == 0){
                    throw new ApplicationException("Edge Index of a boundary edge should not be zero");
                }
                return _ptch.Values[inp.EdgeTag - 1][0];
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
                return 0;
            }
        }

        class CahnHilliardMuSource : mu_Source {

            // OpenFoamPatchField _ptch;

            public CahnHilliardMuSource(CellMask Subgrid = null)
                : base("c"){
            }
            // public bool GetIsDiri(ref CommonParamsBnd inp){
            //     return this.IsDirichlet(ref inp);
            // }
        }

        class CahnHilliardMuDiff : mu_Diffusion {

            OpenFoamPatchField _ptch;

            public CahnHilliardMuDiff(OpenFoamPatchField ptch, double penalty_const, double __cahn, CellMask Subgrid = null)
                : base(3, penalty_const, __cahn, null, "c"){
                _ptch = ptch;
            }

            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                // return !_ptch.IsDirichlet(inp.EdgeTag);
                // return _ptch.IsDirichlet(inp.EdgeTag);
                return true;
                // return false;
            }

            /// <summary>
            /// Dirichlet boundary value
            /// </summary>
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) {
                if (inp.EdgeTag == 0){
                    throw new ApplicationException("Edge Index of a boundary edge should not be zero");
                }
                // Console.WriteLine("diriMu " + _ptch.Values[inp.EdgeTag - 1][0]);
                return 0.0;
                // return _ptch.Values[inp.EdgeTag - 1][0];
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
                // throw new Exception("g_neum should not be called");
                return 0;
            }
        }

        class CahnHilliardCDiff : __c_Diffusion {

            OpenFoamPatchField _ptch;

            public CahnHilliardCDiff(OpenFoamPatchField ptch, double penalty_const, double __diff, double __lambda, CellMask Subgrid = null)
                : base(3, penalty_const, __diff, __lambda, null){
                this._ptch = ptch;
            }
            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                return _ptch.IsDirichlet(inp.EdgeTag);
            }
            public bool GetIsDiri(int et){
                return _ptch.IsDirichlet(et);
            }

            /// <summary>
            /// Dirichlet boundary value
            /// </summary>
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) { // TODO generalize
                // throw new Exception("g_diri should not be called");
                if (inp.EdgeTag == 0){
                    throw new ApplicationException("Edge Index of a boundary edge should not be zero");
                }
                // Console.WriteLine("diriC " + _ptch.Values[inp.EdgeTag - 1][0]);
                return _ptch.Values[inp.EdgeTag - 1][0];
                // return 1.0;
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
                // throw new Exception("g_neum should not be called");
                return 0;
            }
        }
        class CahnHilliardCConv : __c_Flux {
            public CahnHilliardCConv(DGField[] Velocity, CellMask Subgrid = null)
                : base(3, null){} // TODO check if velocity is really communicated through parameter fields

        }

        // class CahnHilliardCSource : c_Source {

        //     OpenFoamPatchField _ptch;

        //     public CahnHilliardCSource(OpenFoamPatchField ptch)
        //         : base(1.0){
        //         this._ptch = ptch;
        //     }
        // }
    }
}
