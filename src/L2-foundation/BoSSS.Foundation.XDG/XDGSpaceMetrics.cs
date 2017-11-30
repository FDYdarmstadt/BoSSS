using ilPSP;
using ilPSP.Tracing;
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
        internal XDGSpaceMetrics(LevelSetTracker lsTrk, XQuadFactoryHelper qfHelper, int __quadorder, SpeciesId[] speciesIds) {
            using(new FuncTrace()) {
                // ----
                // init 
                // ----

                if(!speciesIds.IsSubsetOf(lsTrk.SpeciesIdS)) {
                    throw new ArgumentException();
                }
                CutCellQuadOrder = __quadorder;
                m_qfHelper = qfHelper;
                this.Tracker = lsTrk;
                this.m_SpeciesList = speciesIds.ToArray();

                // ---------------------
                // compute all the stuff
                // ---------------------

                m_CutCellMetrics = new CutCellMetrics(this);
                m_MassMatrixFactory = new MassMatrixFactory(this);
                m_XQuadSchemeHelper = new XQuadSchemeHelper(this);
            }
        }

        /// <summary>
        /// Provides access to quadrature factories; however, most of the time the user wants to use schemes, <see cref="XQuadSchemeHelper"/>.
        /// </summary>
        public XQuadFactoryHelper XQuadFactoryHelper {
            get {
                return m_qfHelper;
            }
        }

        XQuadSchemeHelper m_XQuadSchemeHelper;

        /// <summary>
        /// Provides access to quadrature schemes on cut-cell domains.
        /// </summary>
        public XQuadSchemeHelper XQuadSchemeHelper {
            get {
                return m_XQuadSchemeHelper;
            }
        }


        MassMatrixFactory m_MassMatrixFactory;
        
        /// <summary>
        /// Ye olde provider of the best mass matrices in town. 
        /// </summary>
        public MassMatrixFactory MassMatrixFactory {
            get {
                return m_MassMatrixFactory;
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

        CutCellMetrics m_CutCellMetrics;

        /// <summary>
        /// Cell measures and metrics for cut cells
        /// </summary>
        public CutCellMetrics CutCellMetrics {
            get {
                return m_CutCellMetrics;
            }
        }


    }
}
