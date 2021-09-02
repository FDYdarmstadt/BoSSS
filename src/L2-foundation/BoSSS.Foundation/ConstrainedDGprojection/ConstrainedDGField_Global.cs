using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Foundation.ConstrainedDGprojection {
    
    /// <summary>
    /// Projection of a DG field onto a continuous subspace of the DG space:
    /// Here, the projection onto a continuous sub-space by solving a global system;
    /// This results in the best continuous approximation in the L2-sense, but might be expensive to compute
    /// tan the patch-wise implementation <see cref="ConstrainedDgField_Patchwise"/>.
    /// </summary>
    public class ConstrainedDGField_Global : ConstrainedDGFieldMk3 {
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="b"><see cref="ConstrainedDGFieldMk3.Basis"/></param>
        /// <param name="__domainLimit">
        /// <see cref="ConstrainedDGFieldMk3.domainLimit"/>
        /// </param>
        public ConstrainedDGField_Global(Basis b, CellMask __domainLimit) : base(b, __domainLimit) {
            m_grd = (GridData) b.GridDat;
            mySolver = new ConstrainedProjectionInternal(this, base.domainLimit, base.domainLimit, false);
        }

        GridData m_grd;

        /// <summary>
        /// Release of internal solver
        /// </summary>
        public override void Dispose() {
            if(mySolver != null) {
                mySolver.Dispose();
                mySolver = null;
            }
        }

        ConstrainedDGFieldMk3.ConstrainedProjectionInternal mySolver;


        /// <summary>
        /// Projects some DG field <paramref name="orgDGField"/> onto the internal, continuous representation.
        /// </summary>
        /// <param name="orgDGField">
        /// input; unchanged on exit
        /// </param>
        override public void ProjectDGField(ConventionalDGField orgDGField) {
            using(new FuncTrace()) {
                MPICollectiveWatchDog.Watch();

                if(orgDGField.Basis.Degree > this.Basis.Degree)
                    throw new ArgumentException("continuous projection on a lower degree basis is not recommended");

                SetDGCoordinatesOnce(orgDGField, domainLimit);
                UpdateInternalProjection(csMPI.Raw._COMM.WORLD);

                this.mySolver.PerformProjection();

                UpdateInternalProjection(csMPI.Raw._COMM.WORLD);
            }
        }
    }
}
