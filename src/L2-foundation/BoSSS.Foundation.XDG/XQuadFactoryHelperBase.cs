using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Foundation.XDG {
    abstract public class XQuadFactoryHelperBase {
        internal XQuadFactoryHelperBase(LevelSetTracker.LevelSetData[] lsDatas) {
            //var lsTrk = lsDatas[0].Tracker;
            int iHi = lsDatas[0].HistoryIndex;
            //if(lsDatas.Length != lsTrk.LevelSets.Count)
            //    throw new ArgumentException();
            for (int iLs = 0; iLs < lsDatas.Length; iLs++) {
                if (lsDatas[iLs].LevelSetIndex != iLs)
                    throw new ArgumentException();
                //if(!object.ReferenceEquals(lsDatas[iLs].Tracker, lsTrk))
                //    throw new ArgumentException();
                if (lsDatas[iLs].HistoryIndex != iHi)
                    throw new ArgumentException();
            }

            this.m_LevelSetDatas = lsDatas.CloneAs();
        }

#if DEBUG
        public static bool CheckQuadRules = true;
#else
        public static bool CheckQuadRules = false;
#endif

        protected readonly LevelSetTracker.LevelSetData[] m_LevelSetDatas;


        /// <summary>
        /// Generates a quadrature rule factory for the cut edge integrals.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefVol);
        // ref (2), on edges


        /// <summary>
        /// Generates a quadrature rule factory for the cut volume integrals.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref);
        // ref (3), on cells

        /// <summary>
        /// Generates a quadrature rule factory for integrating over the zero-level-set surface.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref);
        // ref (1), on cell


        /// <summary>
        /// Generates a quadrature rule factory the intersection of levelset0 and levelset1 where levelset0 = levelset1 = 0
        /// This is a point in 2D, a line in 3D.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory);
        // ref (6), on cell

        /// <summary>
        /// Returns a rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e. on \f$  K \cap \mathfrak{I}\f$ .
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="QuadRule"/>'s on edges
        /// </returns>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol);
        // ref (5a), on edge

        /// <summary>
        /// Returns a rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e. on \f$  K \cap \mathfrak{I}\f$ .
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="QuadRule"/>'s on edges
        /// </returns>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory);
        // ref (5a), on edge
    }

}
