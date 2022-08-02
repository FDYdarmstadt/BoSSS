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

        /// <summary>
        ///
        /// </summary>
        [CodeGenExport]
        // public void Laplacian(OpenFoamMatrix mtx) {
        // public void CahnHilliard(OpenFoamMatrix mtx, OpenFoamMatrix Umtx, OpenFoamPatchField ptch, OpenFoamPatchField ptchU) {
        public void CahnHilliard(OpenFoamMatrix mtx, OpenFoamPatchField ptch, OpenFoamPatchField ptchU) {
            try{

                // grid, etc
                // =========

                GridData grd = mtx.ColMap.GridDat as GridData;
                // PlotGrid("grid.plt", grd);

                Basis b = mtx.ColMap.BasisS[0];
                Basis bPhi = mtx.ColMap.BasisS[0];

                OpenFoamDGField C = mtx.Fields[0];
                // SinglePhaseField Phi = new(b);
                OpenFOAMGrid ofGrid = ptch.Grid;
                OpenFoamDGField fields = new(ofGrid, b.Degree, 2);
                OpenFoamMatrix fullMtx = new(ofGrid, fields);

                var map = new UnsetteledCoordinateMapping(b, bPhi);

                int nParams = 7;
                var op = new SpatialOperator(2, nParams, 2, QuadOrderFunc.Linear(), "c", "phi", "VelocityX", "VelocityY", "VelocityZ","c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_phi");
                // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "phi", "c0", "VelocityX", "VelocityY", "VelocityZ", "c_Res", "phi_Res");
                // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "phi", "c0","LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "c_Res", "phi_Res");


                var CHCdiff = new CahnHilliardCDiff(ptch);
                // var CHCsource = new CahnHilliardCSource(ptch);
                var CHPhidiff = new CahnHilliardPhiDiff(ptch);
                var CHPhisource = new CahnHilliardPhiSource();
                op.EquationComponents["Res_c"].Add(CHCdiff);
                // op.EquationComponents["Res_c"].Add(CHCsource);
                op.EquationComponents["Res_phi"].Add(CHPhisource);
                op.EquationComponents["Res_phi"].Add(CHPhidiff);
                // op.LinearizationHint = LinearizationHint.GetJacobiOperator;

                double[] MassScales = new double[2];
                MassScales[0] = 1.0;
                // MassScales[1] = 1.0;
                op.TemporalOperator = new ConstantTemporalOperator(op, MassScales);

                op.LinearizationHint = LinearizationHint.GetJacobiOperator;
                op.Commit();

                // var RealLevSet = new LevelSet(b, "Levset");
                // var RealTracker = new LevelSetTracker((GridData)(b.GridDat), XQuadFactoryHelper.MomentFittingVariants.Saye, 2, new string[] { "A", "B" }, RealLevSet);
                // RealTracker.UpdateTracker(0.0);
                // RealLevSet.Clear();
                // // RealLevSet.Acc(1.0, C);

                // int J = b.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                // for (int j = 0; j < J; j++)
                // {
                //     int N = b.GetLength(j);
                //     for (int n = 0; n < N; n++)
                //         RealLevSet.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                // }

                // TODO
                int J = b.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                SinglePhaseField c = new(b);
                SinglePhaseField phi = new(b);
                for (int j = 0; j < J; j++)
                {
                    int N = b.GetLength(j);
                    for (int n = 0; n < N; n++) {
                        c.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                        // phi.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                    }
                }
                double tanh(double x){
                    // return Math.Tanh(x);
                    return x;
                }
                double sqrt(double x){
                    // return Math.Sqrt(x);
                    double sqrt = x / 2;
                    double temp = 0;

                    // Iterate until sqrt is different of temp, that is updated on the loop
                    while(sqrt - temp > 1e-10){
                        // initially 0, is updated with the initial value of 128
                        // (on second iteration = 65)
                        // and so on
                        temp = sqrt;

                        // Then, replace values (256 / 128 + 128 ) / 2 = 65
                        // (on second iteration 34.46923076923077)
                        // and so on
                        sqrt = ( x/temp + temp) / 2;
                    }
                    return sqrt;
                }
                c.Clear();
                ScalarFunction func() {
                    double cahn = 1.0;
                    double radius = 0.1;
                    return ((_3D)((x, y, z) => tanh((-sqrt(((x - 0.3)*(x - 0.3)) + ((y - 0.3)*(y - 0.3))) + radius) / (sqrt(2.0) * cahn)))).Vectorize();
                }
                c.ProjectField(func());
                phi.Laplacian(-1.0, c);
                phi.Acc(-1.0, c);
                phi.ProjectPow(1.0, c, 3.0);
                SinglePhaseField Res_c = new(b);
                SinglePhaseField Res_phi = new(b);
                List<DGField> ParameterMap = new();
                for (int i = 0; i < nParams; i++)
                {
                    ParameterMap.Add(new SinglePhaseField(b));
                }
                // ParameterMap[3].Acc(1.0, C);
                for (int j = 0; j < J; j++)
                {
                    int N = b.GetLength(j);
                    for (int n = 0; n < N; n++)
                        ParameterMap[3].Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                }
                NonLinearSolverConfig nls = new();
                var ls = new DirectSolver.Config();
                nls.SolverCode = NonLinearSolverCode.Newton;
                XdgTimestepping TimeStepper = new(op,
                                                  new SinglePhaseField[]{c, phi},
                                                  new SinglePhaseField[]{Res_c, Res_phi},
                                                  // TimeSteppingScheme.ExplicitEuler,
                                                  TimeSteppingScheme.ImplicitEuler,
                                                  null,
                                                  null,
                                                  ls,
                                                  nls,
                                                  ParameterMap,
                                                  null);

                var tp = new Tecplot(grd.Grid.GridData, 3);
                for (int i = 0; i < 2; i++)
                    Tecplot("plot-pre" + i, 1.0, 3, new SinglePhaseField[2]{c, phi}[i]);

                for (int t = 0; t < 11; t++)
                    {
                        TimeStepper.Solve(t/10.0, 1.0/10.0);
                    for (int i = 0; i < 2; i++)
                        Tecplot("plot_time_" + t + "_field_" + i, 1.0, 3, new SinglePhaseField[2] { c, phi }[i]);
                }
                // evaluate operator
                // =================

                // var eval = op.GetMatrixBuilder(map, ParameterMap, map);
                // eval.ComputeMatrix(fullMtx, fullMtx.RHSbuffer);
                // mtx.RHSbuffer.ScaleV(-1); // convert LHS affine vector to RHS

                // var tp = new Tecplot(grd.Grid.GridData, 3);
                for (int i = 0; i < 2; i++)
                    Tecplot("plot" + i, 1.0, 3, new SinglePhaseField[2]{c, phi}[i]);
                // Console.WriteLine("Computed Cahn Hilliard Matrix, norm is " + fullMtx.InfNorm());
                // Console.WriteLine(fullMtx);
                // Console.WriteLine("Computed Laplacian Matrix, RHS is ");
                // foreach (var elem in mtx.RHSbuffer)
                //     Console.Write(elem + " ");
                // Console.WriteLine();
            } catch (Exception e) {
                // Console.WriteLine(e.GetType());
                // Console.WriteLine(e.Message);
                // Console.WriteLine(e.StackTrace);
                Console.WriteLine(e);
                throw e;
            }
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
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) { // TODO generalize
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

            public CahnHilliardPhiSource()
                : base(false, 1.0){
            }
            // public bool GetIsDiri(ref CommonParamsBnd inp){
            //     return this.IsDirichlet(ref inp);
            // }
        }

        // TODO change everything that requires boundary conditions
        class CahnHilliardPhiDiff : phi_Diffusion {

            OpenFoamPatchField _ptch;

            public CahnHilliardPhiDiff(OpenFoamPatchField ptch)
                : base(3, 1, 1, null){
                _ptch = ptch;
            }

            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                return !_ptch.IsDirichlet(inp.EdgeTag);
            }
            public bool GetIsDiri(int et){
                return _ptch.IsDirichlet(et);
            }

            /// <summary>
            /// Dirichlet boundary value
            /// </summary>
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) { // TODO generalize
                if (inp.EdgeTag == 0){
                    throw new ApplicationException("Edge Index of a boundary edge should not be zero");
                }
                // return 1.0;
                // Console.WriteLine("diriPhi " + _ptch.Values[inp.EdgeTag - 1][0]);
                return _ptch.Values[inp.EdgeTag - 1][0];
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

            public CahnHilliardCDiff(OpenFoamPatchField ptch)
                : base(3, 1, 1, 1, null){
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
                throw new Exception("g_neum should not be called");
                return 0;
            }
        }

        class CahnHilliardCSource : c_Source {

            OpenFoamPatchField _ptch;

            public CahnHilliardCSource(OpenFoamPatchField ptch)
                : base(1.0){
                this._ptch = ptch;
            }
        }
    }
}
