using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
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
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;

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
        SinglePhaseField phi;
        OpenFoamMatrix _mtx;
        OpenFoamPatchField _ptch;
        double rho0 = 1;
        double rho1 = 2;

        /// <summary>
        /// Returns the auxiliary field appearing during the solution of the Cahn-Hilliard system
        /// </summary>
        [CodeGenExport]
        public OpenFoamDGField GetFlux() {
            int D = 3;
            OpenFoamDGField Flux = new(_ptch.Grid, 0, D);
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
        public OpenFoamDGField GetPhi() {
            int D = 3;
            OpenFoamDGField Phi0 = new(_ptch.Grid, 0, D);
            int nCells = _ptch.Grid.NumberOfCells;

            for (int j = 0; j < nCells; j++) {
                Phi0.SetDGcoordinate(0, j, 0, phi.GetMeanValue(j));
            }
            return Phi0;
        }

        /// <summary>
        /// Solves the Cahn-Hilliard equation
        /// </summary>
        [CodeGenExport]
        public void CahnHilliard(OpenFoamMatrix mtx, OpenFoamDGField U, OpenFoamPatchField ptch, OpenFoamPatchField ptchU) {
            try {
            // {

                // System.Diagnostics.Debugger.Launch();
                // Console.WriteLine("Debugger is attached; press enter to continue");
                // Console.ReadLine();

                _mtx = mtx;
                _ptch = ptch;
                double tanh(double x) {
                    // return x*x*x;
                    // return 0;
                    if (x > 0.25) {
                        return 1;
                    }
                    if (x < -0.25) {
                        return -1;
                    }
                    // if (x > 1.5){
                    //     return 1;
                    // }
                    // if (x < -1.5){
                    //     return -1;
                    // }
                    double fact(int n) {
                        double acc = 1;
                        for (int i = 0; i < n; i++) {
                            acc *= i + 1;
                        }
                        return acc;
                    }
                    // return Math.Tanh(x);
                    // double bernoulli(int m, int n){
                    double bernoulli(int n) {
                        if (n % 2 == 1 && n > 1) {
                            return 0;
                        }
                        double choose(int n, int k) {
                            if (k > n) {
                                throw new ArgumentException();
                            }
                            // double denom = (fact(k) * fact(n - k));
                            // if (denom <= 0){
                            //     Console.WriteLine(n);
                            //     Console.WriteLine(k);
                            //     Console.WriteLine(n-k);
                            //     Console.WriteLine(fact(k));
                            //     Console.WriteLine(fact(n-k));
                            //     throw new Exception();
                            // }
                            return fact(n) / (fact(k) * fact(n - k));
                        }
                        double acc = 0;
                        for (int k = 0; k <= n; k++) {
                            for (int r = 0; r <= k; r++) {
                                acc += (pow(-1, r) * choose(k, r) * pow(r, n) / (k + 1));
                            }
                        }
                        return acc;
                    }
                    double taylorTerm(int n) {
                        return pow(2.0, 2 * n) * (pow(2.0, 2 * n) - 1) * bernoulli(2 * n) * pow(x, 2 * n - 1) / fact(2 * n);
                    }
                    double acc = 0;
                    for (int i = 0; i < 8; i++) {
                        acc += taylorTerm(i);
                    }
                    if (acc > 1) {
                        return 1;
                    }
                    if (acc < -1) {
                        return -1;
                    }
                    return acc;
                }
                double pow(double x, int e) { // inefficient but who cares
                    // return Math.Tanh(x);
                    double acc = 1.0;
                    for (int i = 0; i < e; i++)
                        acc *= x;
                    return acc;
                }
                double pow2(double x) {
                    // return Math.Tanh(x);
                    return pow(x, 2);
                }
                double sqrt(double x) {
                    // return Math.Sqrt(x);
                    double sqrt = x / 2;
                    double temp = 0;
                    while (sqrt - temp > 1e-10) {
                        temp = sqrt;
                        sqrt = (x / temp + temp) / 2;
                    }
                    return sqrt;
                }
                // grid, etc
                // =========

                GridData grd = mtx.ColMap.GridDat as GridData;
                // PlotGrid("grid.plt", grd);

                Basis b = mtx.ColMap.BasisS[0];
                Basis bPhi = mtx.ColMap.BasisS[0];

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

                var map = new UnsetteledCoordinateMapping(b, bPhi);

                int nParams = 7;
                var op = new SpatialOperator(2, nParams, 2, QuadOrderFunc.Linear(), "c", "phi", "VelocityX", "VelocityY", "VelocityZ","c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_phi");
                // var op = new XSpatialOperatorMk2(2, nParams, 2, QuadOrderFunc.Linear(), new List<string>{"a", "b"}, "c", "phi", "VelocityX", "VelocityY", "VelocityZ","c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_phi");
                // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "phi", "c0", "VelocityX", "VelocityY", "VelocityZ", "c_Res", "phi_Res");
                // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "phi", "c0","LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "c_Res", "phi_Res");

                // TODO sync from OpenFOAM
                double lambda = 0.0;
                double penalty_const = 1.0;
                // double M = 5e-11*0.02; // mobility parameter
                // double epsilon = 0.5; // capillary width
                double epsilon = 1e-5; // capillary width
                // double cahn = 0.05;
                // double M = 1; // mobility parameter
                double M = sqrt(epsilon); // mobility parameter
                double sigma = 0.063;
                double lam = 3 / (2 * sqrt(2)) * sigma * epsilon; // Holger's lambda
                double diff = M * lam;

                double cahn = 1.0 / (epsilon * epsilon);

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
                // // ScalarFunction func() {
                // //     // double radius = 0.5e-3;
                // //     // return ((_3D)((x, y, z) => tanh((-sqrt(sqrt(sqrt((pow(x - 2.5e-3, 8) + pow(y - 2.5e-3, 8))))) + radius) / (sqrt(2.0) * cahn)))).Vectorize();
                // //     // return ((_3D)((x, y, z) => tanh(((x - 0.05) + 1 * (z))/0.01))).Vectorize();
                // //     return ((_3D)((x, y, z) => tanh(((x - 2.5) + 0.1 * (y - 2.5))/1))).Vectorize();
                // }
                // c.Clear();
                // Console.WriteLine("TODO: remove this later");
                // c.ProjectField(func());
                // SinglePhaseField c = new(b);
                phi = new(b);
                // for (int j = 0; j < J; j++)
                // {
                //     int N = b.GetLength(j);
                //     for (int n = 0; n < N; n++) {
                //         c.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                //         // phi.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                //     }
                // }

                // foreach (var val in new double[]{0, 0.5, 1, 1.5, 1.99, 2.01 ,10})
                //     Console.WriteLine("Correct Value: " + Math.Tanh(val) + "My value: " + tanh(val));

                // MultidimensionalArray LSGrad = RealLevSet.m_gradPhi;
                var domfields = (IReadOnlyDictionary<string, DGField>)(new Dictionary<string, DGField>() { { "c", c }, { "phi", phi } });
                var paramfields = (IReadOnlyDictionary<string, DGField>)(new Dictionary<string, DGField>(){
                        {"VelocityX", u},
                        {"VelocityY", v},
                        {"VelocityZ", w},
                        {"c0", c},
                        {"LevelSetGradient[0]", c}, // TODO
                        {"LevelSetGradient[1]", c},
                        {"LevelSetGradient[2]", c}
                        });
                Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> GetNamedInputFields = delegate (DGField[] fields) {

                    return (domfields, paramfields);
                };
                RealLevSet.Clear();
                RealLevSet.Acc(1.0, c);
                LevelSetUpdater lsu = new LevelSetUpdater(grd, XQuadFactoryHelper.MomentFittingVariants.Classic,
                                                         5, new string[] { "a", "b" },
                                                         GetNamedInputFields,
                                                         RealLevSet, "c", ContinuityProjectionOption.None);

                var RealTracker = lsu.Tracker;
                // phi.Laplacian(-cahn, c);
                // phi.Acc(-1.0, c);
                // phi.ProjectPow(1.0, c, 3.0);
                RealLevSet.Clear();
                RealLevSet.Acc(1.0, c);
                RealTracker.UpdateTracker(0.0);

                Console.WriteLine("TODO: go back to 5 later");
                SubGrid subgr = RealTracker.Regions.GetNearFieldSubgrid(1);
                SubGridBoundaryModes subgrbnd = 0;
                CellMask subgrMask = subgr.VolumeMask;

                var CHCdiff = new CahnHilliardCDiff(ptch, penalty_const, diff, lambda, subgrMask);
                var CHCconv = new CahnHilliardCConv(new[] { u, v, w }, subgrMask);
                // var CHPhidiff = new CahnHilliardPhiDiff(ptch, penalty_const, cahn);
                var CHPhidiff = new CahnHilliardPhiDiff(ptch, penalty_const, 1.0, subgrMask);
                var CHPhisource = new CahnHilliardPhiSource(cahn);
                op.EquationComponents["Res_c"].Add(CHCdiff);
                op.EquationComponents["Res_c"].Add(CHCconv);
                // op.EquationComponents["Res_c"].Add(CHCsource);
                op.EquationComponents["Res_phi"].Add(CHPhisource);
                op.EquationComponents["Res_phi"].Add(CHPhidiff);
                // op.LinearizationHint = LinearizationHint.GetJacobiOperator;

                double[] MassScales = new double[2];
                MassScales[0] = 1.0;
                // MassScales[1] = 1.0;
                // op.TemporalOperator = new ConstantTemporalOperator(op, MassScales);

                op.LinearizationHint = LinearizationHint.GetJacobiOperator;
                op.Commit();

                SinglePhaseField Res_c = new(b);
                SinglePhaseField Res_phi = new(b);
                List<DGField> ParameterMap = new();
                for (int i = 0; i < nParams; i++) {
                    ParameterMap.Add(new SinglePhaseField(b));
                }
                for (int j = 0; j < J; j++) {
                    int N = b.GetLength(j);
                    for (int n = 0; n < N; n++)
                        ParameterMap[3].Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                }
                for (int j = 0; j < 3; j++) {
                    ParameterMap[j] = new[] { u, v, w }[j];
                }
                NonLinearSolverConfig nls = new();
                var ls = new DirectSolver.Config();
                nls.SolverCode = NonLinearSolverCode.Newton;
                // nls.SolverCode = NonLinearSolverCode.Picard;
                // nls.ConvergenceCriterion = 1e-5;
                nls.verbose = true;

                // var ev = op.GetEvaluatorEx(
                //     new CoordinateMapping(grd), null, map);

                // ev.ActivateSubgridBoundary(subgr.VolumeMask, subgrbnd);

                // ev.Evaluate<CoordinateVector>(1.0, 1.0, );
                // ev.ComputeMatrix();
               
                lsu.InitializeParameters(domfields, paramfields);


                var tp = new Tecplot(grd.Grid.GridData, 3);
                Tecplot("plot.1", 0.0, 3, c, phi, RealLevSet, u, v, w);

                // TODO saye instead of hmf
                XdgSubGridTimestepping TimeStepper = new(op,
                                                         new SinglePhaseField[]{c, phi},
                                                         new SinglePhaseField[]{Res_c, Res_phi},
                                                         // TimeSteppingScheme.ExplicitEuler,
                                                         TimeSteppingScheme.ImplicitEuler,
                                                         subgr,
                                                         subgrbnd,
                                                         LinearSolver: ls,
                                                         NonLinearSolver: nls,
                                                         _UpdateLevelset: (() => lsu),
                                                         _LevelSetHandling: LevelSetHandling.LieSplitting,
                                                         _AgglomerationThreshold: 0.0,
                                                         _optTracker: RealTracker
                                                         );


                int timesteps = 3;
                double dt = 2e-3;
                for (int t = 0; t < timesteps; t++) {

                    // var eval = op.GetMatrixBuilder(map, ParameterMap, map);
                    // var eval = op.GetEvaluatorEx(new SinglePhaseField[]{c, phi}.ToArray(), ParameterMap, map);

                    // TimeStepper.ComputeOperatorMatrix(mtx, mtx.RHSbuffer, map, TimeStepper.CurrentState,null,t*dt,t);

                    // eval.ActivateSubgridBoundary(subgr.VolumeMask, subgrbnd);
                    // TimeStepper.ActivateSubgridBoundary(subgr.VolumeMask, subgrbnd);
                    // BlockMsrMatrix mat = new((IBlockPartitioning)mtx.ColPartition);
                    // CoordinateVector cv = new(c,phi);
                    // eval.ComputeMatrix(mat, cv);

                    // mat.RHSbuffer.ScaleV(-1); // convert LHS affine vector to RHS
                    // for (int i = 0; i < mat.NoOfRows; i++) {
                    //     if (mat.GetNoOfNonZerosPerRow(i) > 0) {
                    //         mat.SetDiagonalElement(i, 1.0);
                    //     }
                    // }
                    // eval.Evaluate(1.0,1.0, cv);
                    // eval.ComputeMatrix(mtx, mtx.RHSbuffer);
                    // eval.ComputeMatrix<OpenFoamMatrix,double[]>(mtx, null);

                    TimeStepper.Solve(dt * t / timesteps, dt * 1.0 / timesteps);
                    RealLevSet.Clear();
                    RealLevSet.Acc(1.0, c);
                    RealTracker.UpdateTracker(0.0);
                    Tecplot("plot." + (t + 2), (t + 1) / timesteps, 3, c, phi, RealLevSet, u, v, w);
                }
                // evaluate operator
                // =================

                // var eval = op.GetMatrixBuilder(map, ParameterMap, map);
                // eval.ComputeMatrix(fullMtx, fullMtx.RHSbuffer);
                // mtx.RHSbuffer.ScaleV(-1); // convert LHS affine vector to RHS

                // Console.WriteLine("Computed Cahn Hilliard Matrix, norm is " + fullMtx.InfNorm());
                // Console.WriteLine(fullMtx);
                // Console.WriteLine("Computed Laplacian Matrix, RHS is ");
                // foreach (var elem in mtx.RHSbuffer)
                //     Console.Write(elem + " ");
                // Console.WriteLine();


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

        class CahnHilliardPhiSource : phi_Source {

            // OpenFoamPatchField _ptch;

            public CahnHilliardPhiSource(double _cahn, CellMask Subgrid = null)
                : base(true, _cahn, Subgrid){
            }
            // public bool GetIsDiri(ref CommonParamsBnd inp){
            //     return this.IsDirichlet(ref inp);
            // }
        }

        class CahnHilliardPhiDiff : phi_Diffusion {

            OpenFoamPatchField _ptch;

            public CahnHilliardPhiDiff(OpenFoamPatchField ptch, double penalty_const, double __cahn, CellMask Subgrid = null)
                : base(3, penalty_const, __cahn, null, Subgrid){
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
                // Console.WriteLine("diriPhi " + _ptch.Values[inp.EdgeTag - 1][0]);
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

        class CahnHilliardCDiff : c_Diffusion {

            OpenFoamPatchField _ptch;

            public CahnHilliardCDiff(OpenFoamPatchField ptch, double penalty_const, double __diff, double __lambda, CellMask Subgrid = null)
                : base(3, penalty_const, __diff, __lambda, null, Subgrid: Subgrid){
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
        class CahnHilliardCConv : c_Flux {
            public CahnHilliardCConv(DGField[] Velocity, CellMask Subgrid = null)
                : base(3, Velocity, null, Subgrid: Subgrid){}

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
