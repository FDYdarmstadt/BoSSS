using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding {
    

    public static class FixedOperators_ {

        unsafe public static void Laplacian(ref int GridRef, 
            ref int DgDegree,

            out int ierr) {
            try {
                // grid, etc
                // =========

                GridData grd = (GridData)(Infrastructure.GetObject(GridRef));

                var b = new Basis(grd, DgDegree);
                var map = new UnsetteledCoordinateMapping(b);
                
                var L = new Laplace(1.3, grd.Cells.cj);
                var op = new SpatialOperator(1, 0, 1, QuadOrderFunc.Linear(), "T", "c1");
                op.EquationComponents["c1"].Add(L);
                op.Commit();

                // evaluate operator
                // =================

                var Mtx = new BlockMsrMatrix(map, map);
                double[] B = new double[map.LocalLength];

                var eval = op.GetMatrixBuilder(map, null, map);
                eval.ComputeMatrix(Mtx, B);

                // return data
                // ===========

                throw new NotImplementedException("todo");

            } catch (Exception e) {
                ierr = Infrastructure.ErrorHandler(e);
            }
            ierr = 0;
        }


        class Laplace : BoSSS.Solution.NSECommon.SIPLaplace {
            public Laplace(double penalty_const, MultidimensionalArray cj)
                : base(penalty_const, cj, "T") //
            {
                //m_boundaryCondMap = __boundaryCondMap;
                //m_bndFunc = m_boundaryCondMap.bndFunction["T"];
            }


            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                return true;
            }
        }


    }
}
