using BoSSS.Application.DerivativeTest_XDG_Enriched;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;

namespace DerivativeTest_XDG_Enriched
{
    class Program
    {
        static void Main(string[] args)
        {
            BoSSS.Solution.Application.InitMPI();
            double[] xNodes = GenericBlas.Linspace(0, 1, 5);
            double[] tNodes = GenericBlas.Linspace(0, 1, 5);
            var grid = Grid2D.Cartesian2DGrid(xNodes, tNodes);
            int dgDegree = 1;
            int LsDegree = 1;
            int Quaddegree = 5; //we need a high Quadratur Degree as otherwise the BC on the Edges are integrated incorrectly
            double x0 = 0.25;
            double a = 0.5;
            Func<double, double> FlowFunc = (t => a);

            double cL = 1;
            double cR = 0;

            var LevelSet = new LevelSet(new Basis(grid, LsDegree), "LevelSet");
            var LsTrk = new LevelSetTracker(grid.GridData, XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "A", "B" }, LevelSet);
            XDGBasis basis = new XDGBasis(LsTrk, dgDegree);
            XDGField c = new XDGField(basis, "c");
            XDGField[] ConservativeVars = new XDGField[1];
            ConservativeVars[0] = c;
            var mapping = new CoordinateMapping(ConservativeVars);

            XDGBasis tmp_basis = new XDGBasis(LsTrk, c.Basis.Degree + 1);
            XDGField tmp_field = new XDGField(tmp_basis);
            var enriched_basis = new XDGBasis(LsTrk, c.Basis.Degree + 1);
            XDGField[] enriched_fields = new XDGField[ConservativeVars.Length];
            enriched_fields[0] = new XDGField(tmp_basis);
            enriched_fields[0].AccLaidBack(1, c);
            var enriched_mapping = new CoordinateMapping(enriched_fields);
            int length_R = enriched_fields[0].CoordinateVector.Length;


            var Op = new XSpatialOperatorMk2(new string[] { "c" }, new string[] { "codom1" }, (int[] a_1, int[] b_1, int[] c_1) => Quaddegree, new string[] { "A", "B" });
            //Bulk
            Op.EquationComponents["codom1"].Add(new ScalarTransportFlux("A", x0, cL, cR, FlowFunc));
            Op.EquationComponents["codom1"].Add(new ScalarTransportFlux("B", x0, cL, cR, FlowFunc));
            //Interface
            Op.EquationComponents["codom1"].Add(new UpwindFlux_XDG_Interface(LsTrk, FlowFunc));
            Op.Commit();


            LevelSet.ProjectField((x, t) => x - 0.3 - 0.45 * t);
            LsTrk.UpdateTracker(0.0);
            c.GetSpeciesShadowField("A").ProjectField((_2D)((x, t) => cL + 0.1));
            c.GetSpeciesShadowField("B").ProjectField((_2D)((x, t) => cR - 0.1));

            var FDJbuilder = Op.GetFDJacobianBuilder(mapping, null, enriched_mapping);
            var CheckMatrix = new BlockMsrMatrix(FDJbuilder.CodomainMapping, FDJbuilder.DomainMapping);
            var CheckAffine = new double[FDJbuilder.CodomainMapping.LocalLength];
            FDJbuilder.ComputeMatrix(CheckMatrix, CheckAffine);

            double RelTol = Math.Sqrt(BLAS.MachineEps); // be generous...
            if(RelTol <= 0.0)
                throw new ArithmeticException();
            var ScalarTransportMat = new BlockMsrMatrix(enriched_mapping, mapping);
            var ScalarTransportMat_Affine = new double[enriched_fields[0].CoordinateVector.Length];

            var MatBuilder = Op.GetMatrixBuilder(mapping, null, enriched_mapping);
            MatBuilder.ComputeMatrix(ScalarTransportMat, ScalarTransportMat_Affine);

            double MtxTol = Math.Max(CheckMatrix.InfNorm(), ScalarTransportMat.InfNorm());
            double AffTol = Math.Max(Math.Max(CheckAffine.MPI_L2Norm(), ScalarTransportMat_Affine.MPI_L2Norm()), MtxTol);
            var ErrMatrix = ScalarTransportMat.CloneAs();
            var ErrAffine = ScalarTransportMat_Affine.CloneAs();
            ErrMatrix.Acc(-1.0, CheckMatrix);
            ErrAffine.AccV(-1.0, CheckAffine);
            double LinfMtx = ErrMatrix.InfNorm();
            double L2Aff = ErrAffine.L2NormPow2().Sqrt();
            bool passed1 = (LinfMtx < MtxTol * RelTol);
            bool passed2 = (L2Aff < AffTol * RelTol);
            Console.WriteLine("Finite Difference Jacobian: Matrix/Affine delta norm {0} {1}, passed? {2} {3}", LinfMtx, L2Aff, passed1, passed2);
            // m_passed = m_passed && passed1;
            // m_passed = m_passed && passed2;

        }
    }
}
