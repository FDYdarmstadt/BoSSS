﻿using BoSSS.Foundation.Grid.Classic;
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
    /// the chosen type (<see cref="XQuadFactoryHelperBase.MomentFittingVariants"/>) and order of the cut-cell quadrature, 
    /// therefore these are collected in this central object.
    /// Instances of this object are obtained via <see cref="LevelSetTracker.GetXDGSpaceMetrics(IEnumerable{SpeciesId}, int, int)"/>.
    /// </summary>
    public sealed class XDGSpaceMetrics {

        XQuadFactoryHelperBase m_qfHelper;

        /// <summary>
        /// ctor.
        /// </summary>
        internal XDGSpaceMetrics(LevelSetTracker lsTrk, XQuadFactoryHelperBase qfHelper, int __quadorder, SpeciesId[] speciesIds, int HistoyIndex) {
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
                m_qfHelper = qfHelper;
                this.Tracker = lsTrk;
                this.m_SpeciesList = speciesIds.ToArray();

                m_LevelSetRegions = lsTrk.RegionsHistory[HistoyIndex];
                m_LevelSetData = lsTrk.DataHistories.Select(his => his[HistoyIndex]).ToList().AsReadOnly();
                foreach(var lsData in m_LevelSetData) {
                    Debug.Assert(lsData.HistoryIndex == HistoyIndex);
                }

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
        public XQuadFactoryHelperBase XQuadFactoryHelper {
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

        LevelSetTracker.LevelSetRegions m_LevelSetRegions;

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
        /// - index: lecel set index
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
        public XQuadFactoryHelperBase.MomentFittingVariants CutCellQuadratureType {
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


    }
}
