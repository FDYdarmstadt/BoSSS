using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Properties of the discrete XDG space. Note that the properties of discrete XDG space
    /// (e.g. measures of cut-cells and the mass matrix) depend on the 
    /// the chosen type (<see cref="CutCellQuadratureMethod"/>) and order of the cut-cell quadrature, 
    /// therefore these are collected in this central object.
    /// Instances of this object are obtained via <see cref="LevelSetTracker.GetXDGSpaceMetrics(IEnumerable{SpeciesId}, int, int)"/>.
    /// </summary>
    public sealed class XDGSpaceMetrics {

        
        static XQuadFactoryHelperBase GetXQuadFactoryHelper(CutCellQuadratureMethod variant, LevelSetTracker.LevelSetData[] lsDatas) {
            XQuadFactoryHelperBase qfHelper;
            if(variant == CutCellQuadratureMethod.Algoim) {
                qfHelper = new XQuadFactoryHelperAlgoim(
                    lsDatas);
            } else {
                qfHelper = new XQuadFactoryHelper(
                    lsDatas,
                    variant);
            }
            return qfHelper;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        internal XDGSpaceMetrics(LevelSetTracker lsTrk, int __quadorder, SpeciesId[] speciesIds, int HistoryIndex) {
            using(new FuncTrace("XDGSpaceMetrics.ctor")) {
                // ----
                // init 
                // ----
#if TEST
                MPICollectiveWatchDog.WatchAtRelease();
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
#endif
                if (!speciesIds.IsSubsetOf(lsTrk.SpeciesIdS)) {
                    throw new ArgumentException();
                }
                CutCellQuadOrder = __quadorder;
                this.Tracker = lsTrk;
                this.m_SpeciesList = speciesIds.ToArray();

                m_LevelSetRegions = lsTrk.RegionsHistory[HistoryIndex];
                m_LevelSetData = lsTrk.DataHistories.Select(his => his[HistoryIndex]).ToList().AsReadOnly();
                foreach(var lsData in m_LevelSetData) {
                    Debug.Assert(lsData.HistoryIndex == HistoryIndex);
                }

                // ---------------------
                // compute all the stuff
                // ---------------------
                m_qfHelper = new XQuadFactoryHelperCached(
                    (int iThread) => GetXQuadFactoryHelper(lsTrk.CutCellQuadratureType, m_LevelSetData.ToArray()));
                m_XQuadSchemeHelper = new XQuadSchemeHelper(this);
                m_qfHelper.CreateRulesAndMPIExchgange(this.CutCellQuadOrder);

                m_CutCellMetrics = new CutCellMetrics(this);
                m_MassMatrixFactory = new MassMatrixFactory(this);
            }
        }

        XQuadFactoryHelperCached m_qfHelper;


        /// <summary>
        /// Provides access to quadrature factories; however, most of the time the user wants to use schemes, <see cref="XQuadSchemeHelper"/>.
        /// </summary>
        public IXQuadFactoryHelper XQuadFactoryHelper {
            get {
                return m_qfHelper;                
            }
        }
        

        readonly XQuadSchemeHelper m_XQuadSchemeHelper;

        /// <summary>
        /// Provides access to quadrature schemes on cut-cell domains.
        /// </summary>
        public XQuadSchemeHelper XQuadSchemeHelper {
            get {
                return m_XQuadSchemeHelper;
            }
        }
 

        readonly MassMatrixFactory m_MassMatrixFactory;
        
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
            set;
        }

        /// <summary>
        /// Underlying background mesh of the XDG space.
        /// </summary>
        public GridData GridDat {
            get {
                return Tracker.GridDat;
            }
        }

        readonly LevelSetTracker.LevelSetRegions m_LevelSetRegions;

        /// <summary>
        /// Constant during object lifetime.
        /// </summary>
        public LevelSetTracker.LevelSetRegions LevelSetRegions {
            get {
                return m_LevelSetRegions;
            }
        }

        ReadOnlyCollection<LevelSetTracker.LevelSetData> m_LevelSetData;

        /// <summary>
        /// Data, e.g. level set gradients, 
        /// constant during object lifetime. 
        /// - index: level set index
        /// </summary>
        public IList<LevelSetTracker.LevelSetData> LevelSetData {
            get {
                return m_LevelSetData;
            }
        }

        /// <summary>
        /// The number of level-sets used.
        /// </summary>
        public int NoOfLevelSets {
            get {
                Debug.Assert(m_LevelSetData.Count == Tracker.NoOfLevelSets);
                return m_LevelSetData.Count;
            }
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
        public CutCellQuadratureMethod CutCellQuadratureType {
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

        //SpeciesId[] m_TotalSpeciesList;

        /// <summary>
        /// All species <see cref="LevelSetTracker.SpeciesIdS"/>.
        /// </summary>
        public IList<SpeciesId> TotalSpeciesList {
            get {
                return Tracker.SpeciesIdS;
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

		/// <summary>
		/// Write all the quadrature rules (edge, volume and surface) as nodes + weights into vtp files
		/// </summary>
		public void WriteAllQuadratureRulesToVtp() {
			var cutCellMetrics = CutCellMetrics;
			cutCellMetrics.WriteSurfaceRulesToVtp();
			cutCellMetrics.WriteEdgeRulesToVtp();
			cutCellMetrics.WriteVolumeRulesToVtp();
		}

	}
}
