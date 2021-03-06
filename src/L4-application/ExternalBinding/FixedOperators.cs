﻿using BoSSS.Foundation;
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
        public void Laplacian(OpenFoamMatrix mtx) {

            // grid, etc
            // =========

            GridData grd = mtx.ColMap.GridDat as GridData;

            var b = mtx.ColMap.BasisS[0];
            var map = new UnsetteledCoordinateMapping(b);

            var L = new Laplace(1.3);
            var op = new SpatialOperator(1, 0, 1, QuadOrderFunc.Linear(), "T", "c1");
            op.EquationComponents["c1"].Add(L);
            op.Commit();

            // evaluate operator
            // =================

            var eval = op.GetMatrixBuilder(map, null, map);
            eval.ComputeMatrix(mtx, mtx.RHSbuffer);
            mtx.RHSbuffer.ScaleV(-1); // convert LHS affine vector to RHS

            Console.WriteLine("Computed Laplacian Matrix, norm is " + mtx.InfNorm());
        }

        /// <summary>
        /// SIP-form for the Laplacian
        /// </summary>
        class Laplace : BoSSS.Solution.NSECommon.SIPLaplace {

            /// <summary>
            /// 
            /// </summary>
            public Laplace(double penalty_const)
                : base(penalty_const, "T") //
            {
                //m_boundaryCondMap = __boundaryCondMap;
                //m_bndFunc = m_boundaryCondMap.bndFunction["T"];
            }

            /// <summary>
            /// always true
            /// </summary>
            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                return true;
            }
        }


    }
}
