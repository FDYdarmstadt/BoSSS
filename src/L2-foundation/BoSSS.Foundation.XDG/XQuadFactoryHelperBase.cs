using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature;
using ilPSP;
using System;

namespace BoSSS.Foundation.XDG {

    abstract public class XQuadFactoryHelperBase : IXQuadFactoryHelper {



        /// <summary>
        /// Used type of the HMF.
        /// </summary>
        public CutCellQuadratureMethod CutCellQuadratureType {
            get;
            protected set;
        }

        internal XQuadFactoryHelperBase(LevelSetTracker.LevelSetData[] lsDatas) {
            //var lsTrk = lsDatas[0].Tracker;
            int iHi = lsDatas[0].HistoryIndex;
            //if(lsDatas.Length != lsTrk.LevelSets.Count)
            //    throw new ArgumentException();
            for(int iLs = 0; iLs < lsDatas.Length; iLs++) {
                if(lsDatas[iLs].LevelSetIndex != iLs)
                    throw new ArgumentException();
                //if(!object.ReferenceEquals(lsDatas[iLs].Tracker, lsTrk))
                //    throw new ArgumentException();
                if(lsDatas[iLs].HistoryIndex != iHi)
                    throw new ArgumentException();
            }

            this.m_LevelSetDatas = lsDatas.CloneAs();
        }

        protected IGridData gdat {
            get {
                return this.m_LevelSetDatas[0].GridDat;
            }
        }


#if DEBUG
        public static bool CheckQuadRules = true;
#else
        public static bool CheckQuadRules = false;
#endif

        /// <summary>
        /// Triggers the creation of all quadrature rules, so that later, the can be retrieved from the cache.
        /// This is supposed to give a better runtime behavior, since the creation of quadrature rules is quite expensive.
        /// Furthermore, in the creation of quadrature rules there are some parts where MPI communication is needed or makes things more robust.
        /// One example are hanging node on MPI boundaries; There, a quadrature rule might be difficult to be created right locally.
        /// </summary>
        abstract public void CreateRulesAndMPIExchgange(int __quadorder);

        protected readonly LevelSetTracker.LevelSetData[] m_LevelSetDatas;

        /// <summary>
        /// Generates a quadrature rule factory for integrating over the zero-level-set surface.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref);
        // ref (1), on cell

        /// <summary>
        /// Generates a quadrature rule factory for integrating over a surface.
        /// The surface is defined by two conditions: levelset0 = 0 and on side jmp1 of levelset1
        /// <see cref="GetSurfaceFactory(int, RefElement)"/>
        /// </summary>
        /// <param name="levSetIndex0">Index of the primary level set defining the surface quadrature rule.</param>
        /// <param name="levSetIndex1">Index of the auxiliary level set for species differentiation.</param>
        /// <param name="jmp1">Specifies the side of <paramref name="levSetIndex1"/> that defines the phase of interest.</param>
        /// <param name="KrefVol">The reference element used for the quadrature rule.</param>
        /// <returns></returns>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);

        /// <summary>
        /// Generates a quadrature rule factory for the cut edge integrals.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefEdge);
        // ref (2), on edges

        /// <summary>
        /// Generates an edge quadrature rule factory for edges cut by two level sets. <see cref="GetEdgeRuleFactory(int, JumpTypes, RefElement)"/>
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);

        /// <summary>
        /// Generates a quadrature rule factory for the cut volume integrals.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref);
        // ref (3), on cells

        /// <summary>
        /// Generates a volume quadrature rule factory for cells cut by two level sets. <see cref="GetVolRuleFactory(int, JumpTypes, RefElement)"/>
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);

        /// <summary>
        /// Returns an edge rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e., for integrals 
        /// ```math 
        ///    \int_{ E \cap \mathfrak{I} } \ldots \textrm{dl} . 
        /// ```
        /// (elements on the zero-level-set surface), i.e. on $E  \cap \mathfrak{I} $ for each edge $ E $.
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefEdge);

        /// <summary>
        /// Returns a cell boundary rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e., for integrals 
        /// ```math 
        ///    \int_{ \partial K \cap \mathfrak{I} } \ldots \textrm{dl} . 
        /// ```
        /// (elements on the zero-level-set surface), i.e. on $  \partial K \cap \mathfrak{I} $ for each cell $ k $.
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        abstract public IQuadRuleFactory<CellBoundaryQuadRule> _GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol);


        /// <summary>
        /// Generates a volume quadrature rule factory for cells cut by two level sets. <see cref="GetSurfaceElement_BoundaryRuleFactory(int, RefElement)"/>
        /// The surface is defined by two conditions: levelset0 = 0 and on side jmp1 of levelset1
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);
        // ref (5a), on edge

        /// <summary>
        /// Generates a quadrature rule factory the intersection of levelset0 and levelset1 where levelset0 = levelset1 = 0
        /// This is a point in 2D, a line in 3D.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol);
        // ref (6), on cell

    }

}
