using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using log4net.Core;
using System;
using System.Collections.Generic;
using System.Data;
using System.Text;

namespace BoSSS.Foundation.XDG {
    abstract public class XQuadFactoryHelperBase {

        /// <summary>
        /// Different variants of the moment-fitting procedure for the creation
        /// of the surface and volume quadrature rules.
        /// </summary>
        public enum MomentFittingVariants {

            /// <summary>
            /// The original method published in 2013 which uses a two-step
            /// procedure: The surface rules are created first and then used to
            /// create the volume rules
            /// </summary>
            Classic,

            /// <summary>
            /// One-step variant proposed by Florian (see XNSE paper, submitted
            /// 2015). Surface and volume rules are created using a single
            /// moment-fitting by additionally enforcing Gauss' theorem on the
            /// discrete level.
            /// </summary>
            OneStepGauss,

            /// <summary>
            /// Same as <see cref="OneStepGauss"/>, but additionally enforces
            /// Stokes' theorem on a discrete level.
            /// </summary>
            OneStepGaussAndStokes,

            /// <summary>
            /// Two step-procedure: using Stokes theorem to create surface rules, 
            /// and the Gauss theorem to create Volume rules.
            /// </summary>
            TwoStepStokesAndGauss,


            /// <summary>
            /// Only for debugging purpose, see <see cref="ExactCircleLevelSetIntegration"/>, <see cref="ExactCircleLevelSetIntegration.RADIUS"/>
            /// </summary>
            ExactCircle,

            /// <summary>
            /// Gaussian quadrature rules for <see cref="Square"/> and <see cref="Cube"/> elements,
            /// obtained throug recursive subdivision, as described in 
            /// (Saye 2015)
            /// </summary>
            /// <remarks>
            /// High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles,
            /// R. Saye, SIAM Journal on Scientific Computing, 2015
            /// </remarks>
            Saye,
            /// <summary>
            /// Gaussian quadrature rules for <see cref="Square"/> and <see cref="Cube"/> elements,
            /// obtained through recursive subdivision, as described in 
            /// (Saye 2022)
            /// </summary>
            Algoim,
        }

        /// <summary>
        /// Used type of the HMF.
        /// </summary>
        public MomentFittingVariants CutCellQuadratureType {
            get;
            protected set;
        }

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
        /// Generates a quadrature rule factory for integrating over the zero-level-set surface.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref);
        // ref (1), on cell

        /// <summary>
        /// Generates a quadrature rule factory for integrating over a surface.
        /// The surface is defined by two conditions: levelset0 = 0 and on side jmp1 of levelset1
        /// <see cref="GetSurfaceFactory(int, RefElement)"/>
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory);

        /// <summary>
        /// Generates a quadrature rule factory for the cut edge integrals.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefVol);
        // ref (2), on edges

        /// <summary>
        /// Generates an edge quadrature rule factory for edges cut by two level sets. <see cref="GetEdgeRuleFactory(int, JumpTypes, RefElement)"/>
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory);

        /// <summary>
        /// Generates a quadrature rule factory for the cut volume integrals.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref);
        // ref (3), on cells

        /// Generates a volume quadrature rule factory for cells cut by two level sets. <see cref="GetVolRuleFactory(int, JumpTypes, RefElement)"/>
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory);

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
        /// Generates a volume quadrature rule factory for cells cut by two level sets. <see cref="GetSurfaceElement_BoundaryRuleFactory(int, RefElement)"/>
        /// The surface is defined by two conditions: levelset0 = 0 and on side jmp1 of levelset1
        /// </summary>
        /// <param name="levSetIndex0"></param>
        /// <param name="levSetIndex1"></param>
        /// <param name="jmp1"></param>
        /// <param name="KrefVol"></param>
        /// <param name="backupFactory"></param>
        /// <returns></returns>
        abstract public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory);
        // ref (5a), on edge

        /// <summary>
        /// Generates a quadrature rule factory the intersection of levelset0 and levelset1 where levelset0 = levelset1 = 0
        /// This is a point in 2D, a line in 3D.
        /// </summary>
        abstract public IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory);
        // ref (6), on cell

    }

}
