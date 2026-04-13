using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature;
using ilPSP;
using System;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Different variants of the moment-fitting procedure for the creation
    /// of the surface and volume quadrature rules.
    /// </summary>
    public enum CutCellQuadratureMethod {

        /// <summary>
        /// The original hierarchical moment fitting method published in 2013 which uses a two-step 
        /// procedure: The surface rules are created first and then used to
        /// create the volume rules
        /// </summary>
        Classic = 0,

        /// <summary>
        /// One-step variant proposed by Florian (see XNSE paper, submitted
        /// 2015). Surface and volume rules are created using a single
        /// moment-fitting by additionally enforcing Gauss' theorem on the
        /// discrete level.
        /// </summary>
        OneStepGauss = 1,

        /// <summary>
        /// Same as <see cref="OneStepGauss"/>, but additionally enforces
        /// Stokes' theorem on a discrete level.
        /// </summary>
        OneStepGaussAndStokes = 2,

        /// <summary>
        /// Two step-procedure: using Stokes theorem to create surface rules, 
        /// and the Gauss theorem to create Volume rules.
        /// </summary>
        TwoStepStokesAndGauss = 3,

        /// <summary>
        /// Gaussian quadrature rules for <see cref="Square"/> and <see cref="Cube"/> elements,
        /// obtained throug recursive subdivision, as described in 
        /// (Saye 2015)
        /// </summary>
        /// <remarks>
        /// High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles,
        /// R. Saye, SIAM Journal on Scientific Computing, 2015
        /// </remarks>
        Saye = 5,

        /// <summary>
        /// Gaussian quadrature rules for <see cref="Square"/> and <see cref="Cube"/> elements,
        /// obtained through recursive subdivision, as described in 
        /// (Saye 2022)
        /// </summary>
        Algoim = 6,
    }

    public interface IXQuadFactoryHelper {
        CutCellQuadratureMethod CutCellQuadratureType { get; }

        /// <summary>
        /// Triggers the creation of all quadrature rules, so that later, the can be retrieved from the cache.
        /// This is supposed to give a better runtime behavior, since the creation of quadrature rules is quite expensive.
        /// Furthermore, in the creation of quadrature rules there are some parts where MPI communication is needed or makes things more robust.
        /// One example are hanging node on MPI boundaries; There, a quadrature rule might be difficult to be created right locally.
        /// </summary>
        void CreateRulesAndMPIExchgange(int __quadorder);


        /// <summary>
        /// Generates a quadrature rule factory for integrating over the zero-level-set surface.
        /// </summary>
        IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref);


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
        IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);

        /// <summary>
        /// Generates a quadrature rule factory for the cut edge integrals.
        /// </summary>
        IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefEdge);
        // ref (2), on edges

        /// <summary>
        /// Generates an edge quadrature rule factory for edges cut by two level sets. <see cref="GetEdgeRuleFactory(int, JumpTypes, RefElement)"/>
        /// </summary>
        IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);

        /// <summary>
        /// Generates a quadrature rule factory for the cut volume integrals.
        /// </summary>
        IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref);
        // ref (3), on cells

        /// <summary>
        /// Generates a volume quadrature rule factory for cells cut by two level sets. <see cref="GetVolRuleFactory(int, JumpTypes, RefElement)"/>
        /// </summary>
        IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);

        /// <summary>
        /// Returns an edge rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e., for integrals 
        /// ```math 
        ///    \int_{ E \cap \mathfrak{I} } \ldots \textrm{dl} . 
        /// ```
        /// (elements on the zero-level-set surface), i.e. on $E  \cap \mathfrak{I} $ for each edge $ E $.
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefEdge);

        /// <summary>
        /// Returns a cell boundary rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e., for integrals 
        /// ```math 
        ///    \int_{ \partial K \cap \mathfrak{I} } \ldots \textrm{dl} . 
        /// ```
        /// (elements on the zero-level-set surface), i.e. on $  \partial K \cap \mathfrak{I} $ for each cell $ k $.
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        IQuadRuleFactory<CellBoundaryQuadRule> _GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol);


        /// <summary>
        /// Generates a volume quadrature rule factory for cells cut by two level sets. <see cref="GetSurfaceElement_BoundaryRuleFactory(int, RefElement)"/>
        /// The surface is defined by two conditions: levelset0 = 0 and on side jmp1 of levelset1
        /// </summary>
        IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol);
        // ref (5a), on edge

        /// <summary>
        /// Generates a quadrature rule factory the intersection of levelset0 and levelset1 where levelset0 = levelset1 = 0
        /// This is a point in 2D, a line in 3D.
        /// </summary>
        IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol);
        // ref (6), on cell
    }

}
