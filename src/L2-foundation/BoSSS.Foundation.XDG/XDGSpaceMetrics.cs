using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Properties of the discrete XDG space. Note that the properties of discrete XDG space
    /// (e.g. measures of cut-cells and the mass matrix) depend on the 
    /// the chosen type (<see cref="XQuadFactoryHelper.MomentFittingVariants"/>) and order of the cut-cell quadrature, 
    /// therefore these are collected in this central object.
    /// Instances of this object are obtained via <see cref="LevelSetTracker.GetXDGSpaceMetrics(XQuadFactoryHelper.MomentFittingVariants, int, int)"/>.
    /// </summary>
    public sealed class XDGSpaceMetrics {

        XQuadFactoryHelper m_qfHelper;



        /// <summary>
        /// ctor.
        /// </summary>
        internal XDGSpaceMetrics(XQuadFactoryHelper qfHelper, int __quadorder) {
            CutCellQuadOrder = __quadorder;
            m_qfHelper = qfHelper;
        }

        /// <summary>
        /// Used as an input for the co0nstruction of <see cref="XQuadSchemeHelper"/>.
        /// </summary>
        public XQuadFactoryHelper XQuadFactoryHelper {
            get {
                return m_qfHelper;
            }
        }

        /// <summary>
        /// The owner object.
        /// </summary>
        public LevelSetTracker Tracker {
            get;
            private set;
        }


        /// <summary>
        /// The quadrature order used for cut cells volumes.
        /// </summary>
        public int CutCellQuadOrder {
            get;
            private set;
        }

        /// <summary>
        /// The type of quadrature which is be used for computing.
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType {
            get {
                return m_qfHelper.CutCellQuadratureType;
            }
        }

        SpeciesId[] m_SpeciesList;

        /// <summary>
        /// All species for which metrics are available, a subset of <see cref="LevelSetTracker.SpeciesIdS"/>.
        /// </summary>
        public IList<SpeciesId> SpeciesList {
            get {
                return m_SpeciesList.CloneAs();
            }
        }

        /// <summary>
        /// Cell measures and metrics for cut cells
        /// </summary>
        public CutCellMetrics CutCellMetrics {
            get {


            }
        }


    }
}
