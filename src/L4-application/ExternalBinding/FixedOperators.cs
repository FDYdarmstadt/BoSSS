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

            // grid, etc
            // =========

            // Console.WriteLine("Test1");
            GridData grd = mtx.ColMap.GridDat as GridData;
            // PlotGrid("grid.plt", grd);

            Console.WriteLine("Test2");
            var b = mtx.ColMap.BasisS[0];
            var bPhi = mtx.ColMap.BasisS[0];

            var C = mtx.Fields[0];
            // SinglePhaseField Phi = new(b);
            OpenFOAMGrid ofGrid = ptch.Grid;
            OpenFoamDGField fields = new(ofGrid, b.Degree, 2);
            OpenFoamMatrix fullMtx = new(ofGrid, fields);

            var map = new UnsetteledCoordinateMapping(b, bPhi);
            // Console.WriteLine("Test3");

            // var L = new Laplace(1.3, ptch);
            // var CH = new Laplace(1.3, ptch);
            // var op = new SpatialOperator(3, 0, 1, QuadOrderFunc.Linear(), "c", "phi", "u", "c0");
            // var op = new SpatialOperator(2, 0, 1, QuadOrderFunc.Linear(), "c", "phi", "c0");
            // var op = new SpatialOperator(2, 7, 2, QuadOrderFunc.Linear(), "c", "phi", "c0", "VelocityX", "VelocityY", "VelocityZ", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "c_Res", "phi_Res");
            var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "phi", "c0", "VelocityX", "VelocityY", "VelocityZ", "c_Res", "phi_Res");
            // Console.WriteLine("Test4");

            // var CHconv = new c_Flux(3);

            var CHCdiff = new CahnHilliardCDiff(ptch);
            Console.WriteLine("Test4.1");
            var CHPhidiff = new CahnHilliardPhiDiff(ptch);
            Console.WriteLine("Test4.2");
            var CHPhisource = new CahnHilliardPhiSource();
            Console.WriteLine("Test5");
            op.EquationComponents["c_Res"].Add(CHCdiff);
            // op.EquationComponents["phi_Res"].Add(CHPhisource);
            op.EquationComponents["phi_Res"].Add(CHPhidiff);
            // op.LinearizationHint = LinearizationHint.GetJacobiOperator;
            Console.WriteLine("Test6");
            op.Commit();
            Console.WriteLine("Test6.1");

            // Console.WriteLine("dirichlet: " + ptch.IsDirichlet(1));
            // Console.WriteLine("dirichlet: " + ptch.IsDirichlet(2));
            // Console.WriteLine("dirichlet: " + ptch.IsDirichlet(3));
            // Console.WriteLine("dirichlet: " + CHCdiff.GetIsDiri(1));
            // Console.WriteLine("dirichlet: " + CHCdiff.GetIsDiri(2));
            // Console.WriteLine("dirichlet: " + CHCdiff.GetIsDiri(3));
            // Console.WriteLine("dirichlet: " + CHPhidiff.GetIsDiri(1));
            // Console.WriteLine("dirichlet: " + CHPhidiff.GetIsDiri(2));
            // Console.WriteLine("dirichlet: " + CHPhidiff.GetIsDiri(3));

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
            // evaluate operator
            // =================

            // Console.WriteLine(op.DomainVar.Count);
            // foreach (var elem in op.DomainVar){
            //     Console.WriteLine(elem);
            // }
            // Console.WriteLine(map.NoOfVariables);

            Console.WriteLine(op.CodomainVar.Count);
            foreach (var elem in op.CodomainVar){
                Console.WriteLine(elem);

            }
            Console.WriteLine(map.NoOfVariables);

            List<DGField> ParameterMap = new();
            for (int i = 0; i < 4; i++){
                ParameterMap.Add(new SinglePhaseField(b));
            }
            var eval = op.GetMatrixBuilder(map, ParameterMap, map);
            Console.WriteLine("Test7");
            eval.ComputeMatrix(fullMtx, fullMtx.RHSbuffer);
            Console.WriteLine("Test8");
            mtx.RHSbuffer.ScaleV(-1); // convert LHS affine vector to RHS

            Console.WriteLine("Computed Cahn Hilliard Matrix, norm is " + fullMtx.InfNorm());
            // Console.WriteLine("Computed Laplacian Matrix, RHS is ");
            // foreach (var elem in mtx.RHSbuffer)
            //     Console.Write(elem + " ");
            // Console.WriteLine();
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
                return _ptch.Values[inp.EdgeTag - 1][0];
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
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
                return _ptch.Values[inp.EdgeTag - 1][0];
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
                return 0;
            }
        }
    }
}
