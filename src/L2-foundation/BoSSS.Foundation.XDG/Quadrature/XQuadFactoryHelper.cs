/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.XDG.Quadrature;
using IntersectingQuadrature;
using static BoSSS.Foundation.XDG.XQuadFactoryHelper;

namespace BoSSS.Foundation.XDG {


    

    /// <summary>
    /// Auxiliary class that helps with the creation of XDG-quadrature schemes;
    /// instances can be obtained via <see cref="LevelSetTracker.GetXQuadFactoryHelper"/>.
    /// </summary>
    public class XQuadFactoryHelper : XQuadFactoryHelperBase {

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
            private set;
        }
        public MultiLevelSetBeckFactoryCreator zwoLSSayeFactories { get; private set; }

        /// <summary>
        /// ctor.
        /// </summary>
        internal XQuadFactoryHelper(LevelSetTracker.LevelSetData[] lsDatas, MomentFittingVariants momentFittingVariant) 
            : base(lsDatas) {


            this.CutCellQuadratureType = momentFittingVariant;
        }


        MultiLevelSetBruteForceQuadratureFactory zwoLSBruteForceFactories;


        //LevelSetTracker lsTrk;

        // -----------------------------------------------------
        // Factory creation

        Quadrature.HMF.LineAndPointQuadratureFactory[] LineAndPoint_in2D = null;
        IQuadRuleFactory<CellBoundaryQuadRule>[] CellFaceVolume_in3D = null;
        IQuadRuleFactory<CellBoundaryQuadRule>[] CellFaceSurface_in3D = null;

        /// <summary>
        /// Returns a rule for the edges of surface-elements (elements on the zero-level-set surface, 
        /// i.e. on \f$  K \cap \mathfrak{I}\f$ .
        /// (point integrals in 2D, Line integrals in 3D)
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="CellBoundaryQuadRule"/>'s
        /// </returns>
        IQuadRuleFactory<CellBoundaryQuadRule> _GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol) {
            int D = this.m_LevelSetDatas[0].GridDat.SpatialDimension;

            if (D == 2) {

                if (LineAndPoint_in2D == null)
                    LineAndPoint_in2D = new LineAndPointQuadratureFactory[this.m_LevelSetDatas.Length];

                if (LineAndPoint_in2D[levSetIndex] == null) {
                    LineAndPoint_in2D[levSetIndex] = new LineAndPointQuadratureFactory(
                        KrefVol,
                        this.m_LevelSetDatas[levSetIndex],
                        CutCellQuadratureType == MomentFittingVariants.OneStepGaussAndStokes);
                }

                return LineAndPoint_in2D[levSetIndex].GetPointFactory();
            } else {
                //throw new NotImplementedException("3d is not implemented yet");
                Debug.Assert(LineAndPoint_in2D == null);

                if (CellFaceSurface_in3D == null)
                    CellFaceSurface_in3D = new IQuadRuleFactory<CellBoundaryQuadRule>[this.m_LevelSetDatas.Length];
                if(CellFaceSurface_in3D[levSetIndex] == null) {
                    var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

                    switch (CutCellQuadratureType) {
                        case MomentFittingVariants.Saye:
                        CellFaceSurface_in3D[levSetIndex] = SayeFactories.SayeGaussRule_EdgeSurface3D(
                            this.m_LevelSetDatas[levSetIndex],
                            rootFindingAlgorithm);
                        break;
                        default:
                        var CoFaceQuadRuleFactory = new CutLineOnEdgeQuadRuleFactory(
                            this.m_LevelSetDatas[levSetIndex],
                            rootFindingAlgorithm,
                            JumpTypes.Heaviside);
                        CellFaceSurface_in3D[levSetIndex] = new LevelSetEdgeSurfaceQuadRuleFactory(
                            this.m_LevelSetDatas[levSetIndex],
                            CoFaceQuadRuleFactory,
                            JumpTypes.Heaviside);
                        //new LevelSetEdgeVolumeQuadRuleFactory(
                        //    lsTrk, levSetIndex, rootFindingAlgorithm, JumpTypes.Heaviside);
                        break;
                    }

                }
                return CellFaceSurface_in3D[levSetIndex];

            }
        }

        /// <summary>
        /// Returns a rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e. on \f$  K \cap \mathfrak{I}\f$ .
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="QuadRule"/>'s on edges
        /// </returns>
        override public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol) {
            var gdat = this.m_LevelSetDatas[levSetIndex].GridDat;
            int D = gdat.SpatialDimension;
            return new EdgeRuleFromCellBoundaryFactory(gdat,
                _GetSurfaceElement_BoundaryRuleFactory(levSetIndex, KrefVol),
                this.m_LevelSetDatas[levSetIndex].Region.GetCutCellMask4LevSet(levSetIndex));
        }

        /// <summary>
        /// Returns a rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e. on \f$  K \cap \mathfrak{I}\f$ .
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="QuadRule"/>'s on edges
        /// </returns>
        public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory)
        {
            //switch (CutCellQuadratureType)
            //{
                //case MomentFittingVariants.Saye:
                //    if (zwoLSSayeFactories == null)
                //    {
                //        zwoLSSayeFactories = new MultiLevelSetBeckFactoryCreator(m_LevelSetDatas);
                //    }
                //    return zwoLSSayeFactories.GetEdgePointRuleFactory(levSetIndex0, levSetIndex1, jmp1, backupFactory);
                //default:
                    //old stuff 
                    if (zwoLSBruteForceFactories == null)
                    {
                        zwoLSBruteForceFactories = new MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetEdgePointRuleFactory(levSetIndex0, levSetIndex1, jmp1, backupFactory);
            //}
        }



        /// <summary>
        /// Quadrature rule on cell boundaries
        /// </summary>
        IQuadRuleFactory<CellBoundaryQuadRule> GetCellFaceFactory(int levSetIndex, RefElement Kref, JumpTypes jumpType) {
            int D = this.m_LevelSetDatas[0].GridDat.SpatialDimension;

            if (D == 2) {
                if (jumpType != JumpTypes.Heaviside && jumpType != JumpTypes.OneMinusHeaviside)
                    throw new NotSupportedException();
                Debug.Assert(CellFaceVolume_in3D == null);

                if (LineAndPoint_in2D == null)
                    LineAndPoint_in2D = new LineAndPointQuadratureFactory[this.m_LevelSetDatas.Length];

                if (LineAndPoint_in2D[levSetIndex] == null) {
                    LineAndPoint_in2D[levSetIndex] = new LineAndPointQuadratureFactory(Kref, this.m_LevelSetDatas[levSetIndex], true);
                }

                return LineAndPoint_in2D[levSetIndex].GetLineFactory(jumpType == JumpTypes.Heaviside ? true : false);
            } else if (D == 3) {
                Debug.Assert(LineAndPoint_in2D == null);
                if (CellFaceVolume_in3D == null)
                    CellFaceVolume_in3D = new IQuadRuleFactory<CellBoundaryQuadRule>[this.m_LevelSetDatas.Length];
                if (jumpType != JumpTypes.Heaviside)
                    throw new NotSupportedException();
                
                if(CellFaceVolume_in3D[levSetIndex] == null) {
                    var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);
                    switch (CutCellQuadratureType) {
                        case MomentFittingVariants.Saye:
                        if (CellFaceVolume_in3D == null)
                            CellFaceVolume_in3D = new SayeGaussEdgeRuleFactory[this.m_LevelSetDatas.Length];
                        CellFaceVolume_in3D[levSetIndex] = SayeFactories.SayeGaussRule_EdgeVolume3D(
                            this.m_LevelSetDatas[levSetIndex], rootFindingAlgorithm);
                        break;

                        case MomentFittingVariants.Algoim:
                            if (CellFaceVolume_in3D == null)
                                CellFaceVolume_in3D = new IQuadRuleFactory<CellBoundaryQuadRule>[this.m_LevelSetDatas.Length];

                            var factory = new AlgoimFactories(this.m_LevelSetDatas[levSetIndex], Kref);
                            CellFaceVolume_in3D[levSetIndex] = factory.GetEdgeVolumeFactory();
                            break;

                        default:
                        if (CellFaceVolume_in3D == null)
                            CellFaceVolume_in3D = new LevelSetEdgeVolumeQuadRuleFactory[this.m_LevelSetDatas.Length];
                        CellFaceVolume_in3D[levSetIndex] = new LevelSetEdgeVolumeQuadRuleFactory(
                            this.m_LevelSetDatas[levSetIndex], rootFindingAlgorithm, JumpTypes.Heaviside);
                        break;
                    }
                }
                return CellFaceVolume_in3D[levSetIndex];
            } else {
                throw new NotSupportedException();
            }
        }

        private void CheckJmp(JumpTypes jmp) {
            if (jmp != JumpTypes.Heaviside && jmp != JumpTypes.OneMinusHeaviside)
                throw new NotSupportedException();
        }

        /// <summary>
        /// Generates a quadrature rule factory for the cut edge integrals.
        /// </summary>
        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefVol) {

            //void Phi(int x, NodeSet nodes, MultidimensionalArray inU, MultidimensionalArray outU)
            //{
            //    ((LevelSet)m_LevelSetDatas[levSetIndex].LevelSet).EvaluateEdge(x, 1, nodes, inU, outU);
            //    inU.Scale(-1);
            //};

            //return new BruteForceQuadratureFactory(new BruteForceEdgeScheme(Phi));

            var gdat = this.m_LevelSetDatas[levSetIndex].GridDat;
            int D = gdat.SpatialDimension;

            if (!gdat.Grid.RefElements.Contains(KrefVol, (a, b) => object.ReferenceEquals(a, b)))
                throw new ArgumentException();

            CheckJmp(jmp);

            if (D == 2) {
                var r = new EdgeRuleFromCellBoundaryFactory(gdat,
                    GetCellFaceFactory(levSetIndex, KrefVol, jmp),
                    m_LevelSetDatas[levSetIndex].Region.GetCutCellMask4LevSet(levSetIndex));
                return r;
            } else {
                if (jmp == JumpTypes.Heaviside) {
                    var r = new EdgeRuleFromCellBoundaryFactory(gdat,
                        GetCellFaceFactory(levSetIndex, KrefVol, JumpTypes.Heaviside),
                        m_LevelSetDatas[levSetIndex].Region.GetCutCellMask4LevSet(levSetIndex));
                    return r;
                } else if (jmp == JumpTypes.OneMinusHeaviside) {

                    return new ComplementaryRuleFactory(GetEdgeRuleFactory(levSetIndex, JumpTypes.Heaviside, KrefVol));
                } else
                    throw new ArgumentOutOfRangeException("unsupported jump type");
            }
        }

        /// <summary>
        /// Generates an edge quadrature rule factory for edges cut by two level sets.
        /// </summary>
        override public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory)
        {
            switch (CutCellQuadratureType)
            {
                case MomentFittingVariants.Saye:
                    if (zwoLSSayeFactories == null)
                    {
                        zwoLSSayeFactories = new MultiLevelSetBeckFactoryCreator(m_LevelSetDatas);
                    }
                    return zwoLSSayeFactories.GetEdgeRuleFactory(levSetIndex0, jmp0, levSetIndex1, jmp1, backupFactory);

                default:

                    if (zwoLSBruteForceFactories == null)
                    {
                        zwoLSBruteForceFactories = new MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetEdgeRuleFactory(levSetIndex0, jmp0, levSetIndex1, jmp1, backupFactory);
            }
        }

        /// <summary>
        /// Generates a quadrature rule factory for the cut volume integrals.
        /// </summary>
        override public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref)
        {
            CheckJmp(jmp);
            var ctx = this.m_LevelSetDatas[levSetIndex].GridDat;

            if (jmp == JumpTypes.Heaviside)
            {
                if (m_SurfaceFactory == null)
                    m_SurfaceFactory = new IQuadRuleFactory<QuadRule>[m_LevelSetDatas.Length];
                if (m_VolumeFactory == null)
                    m_VolumeFactory = new IQuadRuleFactory<QuadRule>[m_LevelSetDatas.Length];

                if (m_VolumeFactory[levSetIndex] == null)
                {
                    switch (CutCellQuadratureType)
                    {
                        case MomentFittingVariants.Classic:
                            m_VolumeFactory[levSetIndex] = new LevelSetVolumeQuadRuleFactory(
                                this.m_LevelSetDatas[levSetIndex],
                                GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                                GetSurfaceFactory(levSetIndex, Kref),
                                jumpType: jmp);
                            break;

                        case MomentFittingVariants.OneStepGauss:
                        case MomentFittingVariants.OneStepGaussAndStokes:
                            {
                                bool bStokes = CutCellQuadratureType == MomentFittingVariants.OneStepGaussAndStokes;
                                LevelSetComboRuleFactory2 ComboRuleFactroy = new LevelSetComboRuleFactory2(
                                        this.m_LevelSetDatas[levSetIndex],
                                        this.GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                                        bStokes ? this._GetSurfaceElement_BoundaryRuleFactory(levSetIndex, Kref) : null,
                                        _UseAlsoStokes: bStokes,
                                        _SurfaceNodesOnZeroLevset: false,
                                        _DoCheck: CheckQuadRules);
                                m_VolumeFactory[levSetIndex] = ComboRuleFactroy.GetVolumeFactory();
                                m_SurfaceFactory[levSetIndex] = ComboRuleFactroy.GetSurfaceFactory();
                                break;
                            }

                        case MomentFittingVariants.TwoStepStokesAndGauss:
                        case MomentFittingVariants.ExactCircle:
                            {
                                m_VolumeFactory[levSetIndex] = (new LevelSetVolumeQuadRuleFactory2b(Kref,
                                        this.m_LevelSetDatas[levSetIndex],
                                        GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                                        GetSurfaceFactory(levSetIndex, Kref),
                                        jmp));
                                break;
                            }
                        case MomentFittingVariants.Saye:
                            var comboFactory = Quadrature.SayeFactories.SayeGaussRule_Combo(
                                this.m_LevelSetDatas[levSetIndex],
                                new LineSegment.SafeGuardedNewtonMethod(1e-14));
                            m_VolumeFactory[levSetIndex] = comboFactory.GetVolumeFactory();
                            m_SurfaceFactory[levSetIndex] = comboFactory.GetSurfaceFactory();
                            break;
                        case MomentFittingVariants.Algoim:
                            var algoimComboFactory = new Quadrature.AlgoimFactories(
                                    this.m_LevelSetDatas[levSetIndex],
                                    Kref);
                            m_VolumeFactory[levSetIndex] = algoimComboFactory.GetVolumeFactory();
                            m_SurfaceFactory[levSetIndex] = algoimComboFactory.GetSurfaceFactory();
                            break;
                        default:
                            throw new NotSupportedException(String.Format(
                                "Variant {0} not implemented.", CutCellQuadratureType));
                    }
                }

                Debug.Assert(m_VolumeFactory[levSetIndex] != null);
                return m_VolumeFactory[levSetIndex];
            }
            else if (jmp == JumpTypes.OneMinusHeaviside)
            {
                IQuadRuleFactory<QuadRule> ret;
                switch (CutCellQuadratureType)
                {
                    case MomentFittingVariants.Saye:
                        ret = Quadrature.SayeFactories.SayeGaussRule_NegativeVolume(this.m_LevelSetDatas[levSetIndex],
                                new LineSegment.SafeGuardedNewtonMethod(1e-14));
                        break;
                    case MomentFittingVariants.Algoim:
                        var algoimComboFactory = new Quadrature.AlgoimFactories(
                                this.m_LevelSetDatas[levSetIndex],
                                Kref,
                                true);
                        ret = algoimComboFactory.GetVolumeFactory();
                        break;
                    default:
                        ret = new ComplementaryRuleFactory(GetVolRuleFactory(levSetIndex, JumpTypes.Heaviside, Kref));
                        break;
                }
                Debug.Assert(ret != null);
                return ret;
            }
            else
            {
                throw new ArgumentOutOfRangeException("unsupported jump type");
            }
        }


        /// <summary>
        /// Generates a volume quadrature rule factory for cells cut by two level sets.
        /// </summary>
        public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory)
        {
            switch (CutCellQuadratureType)
            {
                case MomentFittingVariants.Saye:
                    if (zwoLSSayeFactories == null)
                    {
                        zwoLSSayeFactories = new MultiLevelSetBeckFactoryCreator(m_LevelSetDatas);
                    }
                    return zwoLSSayeFactories.GetVolRuleFactory(levSetIndex0, jmp0, levSetIndex1, jmp1, backupFactory);
                default:
                    if (zwoLSBruteForceFactories == null)
                    {
                        zwoLSBruteForceFactories = new MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetVolRuleFactory(levSetIndex0, jmp0, levSetIndex1, jmp1, backupFactory);
            }
        }

        /// <summary>
        /// Integration orders of all quadrature rules for volume integrals that have been cached so far
        /// </summary>
        public int[] GetCachedVolumeOrders(int levSetIdx) {
            /*
            switch (momentFittingVariant) {
                case MomentFittingVariants.Classic:
                case MomentFittingVariants.ExactCircle:
                case MomentFittingVariants.TwoStepStokesAndGauss:
                if (m_VolumeFactory == null || m_VolumeFactory[levSetIdx] == null)
                    return new int[0];
                else
                    return m_VolumeFactory[levSetIdx].GetCachedRuleOrders();

                case MomentFittingVariants.OneStepGauss:
                case MomentFittingVariants.OneStepGaussAndStokes:
                if (m_ComboRuleFactroy == null || m_ComboRuleFactroy[levSetIdx] == null)
                    return new int[0];
                else
                    return m_ComboRuleFactroy[levSetIdx].GetVolumeFactory().GetCachedRuleOrders();


                default:
                throw new NotImplementedException();
            }
            */

            if (m_VolumeFactory == null || m_VolumeFactory[levSetIdx] == null)
                return new int[0];
            else
                return m_VolumeFactory[levSetIdx].GetCachedRuleOrders();
        }

        //SurfaceStokes_2D[] m_StokesSurface2D;
        //LevelSetVolumeQuadRuleFactory2b[] m_VolumeFactory2b;
        //LevelSetComboRuleFactory2[] m_ComboRuleFactroy;


        IQuadRuleFactory<QuadRule>[] m_SurfaceFactory = null;
        IQuadRuleFactory<QuadRule>[] m_VolumeFactory = null;

        /// <summary>
        /// Generates a quadrature rule factory for integrating over the zero-level-set surface.
        /// </summary>
        override public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref) {
            //if (m_ComboRuleFactroy == null)
            //    m_ComboRuleFactroy = new LevelSetComboRuleFactory2[this.lsTrk.LevelSets.Count];
            if (m_SurfaceFactory == null)
                m_SurfaceFactory = new IQuadRuleFactory<QuadRule>[this.m_LevelSetDatas.Length];
            if (m_VolumeFactory == null)
                m_VolumeFactory = new IQuadRuleFactory<QuadRule>[this.m_LevelSetDatas.Length];
            //if(m_StokesSurface2D == null)
            //    m_StokesSurface2D = new SurfaceStokes_2D[this.lsTrk.LevelSets.Count];

            if (m_SurfaceFactory[levSetIndex] == null) {
                switch (CutCellQuadratureType) {
                    case MomentFittingVariants.Classic:

                    m_SurfaceFactory[levSetIndex] = new LevelSetSurfaceQuadRuleFactory(
                         m_LevelSetDatas[levSetIndex],
                         GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside));
                    break;

                    case MomentFittingVariants.OneStepGauss:
                    case MomentFittingVariants.OneStepGaussAndStokes:
                    {
                        bool bStokes = CutCellQuadratureType == MomentFittingVariants.OneStepGaussAndStokes;
                        var ComboRuleFactroy = new LevelSetComboRuleFactory2(
                                m_LevelSetDatas[levSetIndex],
                                this.GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                                bStokes ? this._GetSurfaceElement_BoundaryRuleFactory(levSetIndex, Kref) : null,
                                _SurfaceNodesOnZeroLevset: false,
                                _DoCheck: CheckQuadRules,
                                _UseAlsoStokes: bStokes);

                        m_VolumeFactory[levSetIndex] = ComboRuleFactroy.GetVolumeFactory();
                        m_SurfaceFactory[levSetIndex] = ComboRuleFactroy.GetSurfaceFactory();
                        break;
                    }
                                        
                    case MomentFittingVariants.TwoStepStokesAndGauss:
                        m_SurfaceFactory[levSetIndex] = (new SurfaceStokes_2D(
                            m_LevelSetDatas[levSetIndex],
                            this.GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                            this._GetSurfaceElement_BoundaryRuleFactory(levSetIndex, Kref),
                            _SurfaceNodesOnZeroLevset: false,
                            _DoCheck: CheckQuadRules)).GetSurfaceFactory();
                        break;

                    case MomentFittingVariants.ExactCircle:
                        return new ExactCircleLevelSetIntegration(levSetIndex, this.m_LevelSetDatas[levSetIndex].GridDat, Kref);
                    case MomentFittingVariants.Saye:
                        var comboFactory = Quadrature.SayeFactories.SayeGaussRule_Combo(
                                this.m_LevelSetDatas[levSetIndex],
                                new LineSegment.SafeGuardedNewtonMethod(1e-14));
                        m_VolumeFactory[levSetIndex] = comboFactory.GetVolumeFactory();
                        m_SurfaceFactory[levSetIndex] = comboFactory.GetSurfaceFactory();
                        break;
                    case MomentFittingVariants.Algoim:
                        var algoimComboFactory = new Quadrature.AlgoimFactories(
                                this.m_LevelSetDatas[levSetIndex],
                                Kref);
                        m_VolumeFactory[levSetIndex] = algoimComboFactory.GetVolumeFactory();
                        m_SurfaceFactory[levSetIndex] = algoimComboFactory.GetSurfaceFactory();
                        break;
                    default:
                        throw new NotSupportedException(String.Format(
                            "Variant {0} not implemented.", CutCellQuadratureType));
                }
            }

            return m_SurfaceFactory[levSetIndex];
        }

        /// <summary>
        /// Generates a quadrature rule factory for integrating over a surface.
        /// The surface is defined by two conditions: levelset0 = 0 and on side jmp1 of levelset1
        /// </summary>
        public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory)
        {
            switch (CutCellQuadratureType)
            {
                case MomentFittingVariants.Saye:
                    if (zwoLSSayeFactories == null)
                    {
                        zwoLSSayeFactories = new MultiLevelSetBeckFactoryCreator(m_LevelSetDatas);
                    }
                    return zwoLSSayeFactories.GetSurfaceFactory(levSetIndex0, levSetIndex1, jmp1, backupFactory);
                default:
                    if (zwoLSBruteForceFactories == null)
                    {
                        zwoLSBruteForceFactories = new MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetSurfaceFactory(levSetIndex0,
                        levSetIndex1,
                        jmp1, backupFactory);
            }
        }


        /// <summary>
        /// Generates a quadrature rule factory the intersection of levelset0 and levelset1 where levelset0 = levelset1 = 0
        /// This is a point in 2D, a line in 3D.
        /// </summary>
        public IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            //switch (CutCellQuadratureType)
            //{
            //    case MomentFittingVariants.OneStepGauss:
            //        if (zwoLSSayeFactories == null)
            //        {
            //            zwoLSSayeFactories = new MultiLevelSetBeckFactoryCreator(m_LevelSetDatas);
            //        }
            //        return zwoLSSayeFactories.GetIntersectionFactory(levSetIndex0, levSetIndex1, backupFactory);
            //    default:
                    if (zwoLSBruteForceFactories == null) {
                zwoLSBruteForceFactories = new MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
            }
            return zwoLSBruteForceFactories.GetIntersectionFactory(levSetIndex0, levSetIndex1, backupFactory);
                    //}
        }

        /// <summary>
        /// Integration orders of all quadrature rules for volume integrals that have been cached so far
        /// </summary>
        public int[] GetCachedSurfaceOrders(int levSetIdx) {
            /*
            switch (momentFittingVariant) {
                case MomentFittingVariants.Classic:
                if (m_SurfaceFactory == null || m_SurfaceFactory[levSetIdx] == null)
                    return new int[0];
                else
                    return m_SurfaceFactory[levSetIdx].GetCachedRuleOrders();

                case MomentFittingVariants.OneStepGauss:
                case MomentFittingVariants.OneStepGaussAndStokes:
                if (m_ComboRuleFactroy == null || m_ComboRuleFactroy[levSetIdx] == null)
                    return new int[0];
                else
                    return m_ComboRuleFactroy[levSetIdx].GetSurfaceFactory().GetCachedRuleOrders();

                default:
                throw new NotImplementedException();
            }
            */
            if (m_SurfaceFactory == null || m_SurfaceFactory[levSetIdx] == null)
                return new int[0];
            else
                return m_SurfaceFactory[levSetIdx].GetCachedRuleOrders();
        }


        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

        /// <summary>
        /// Creates, from a rule for the positive domain (<see cref="JumpTypes.Heaviside"/>)
        /// the rule for the negative domain and vice-versa.
        /// </summary>
        class ComplementaryRuleFactory : IQuadRuleFactory<QuadRule> {

            public ComplementaryRuleFactory(IQuadRuleFactory<QuadRule> orgRule) {
                m_orgrule = orgRule;
            }

            IQuadRuleFactory<QuadRule> m_orgrule;

            /// <summary>
            /// If there are any cached rules, this method returns their order.
            /// </summary>
            public int[] GetCachedRuleOrders() {
                return m_orgrule.GetCachedRuleOrders();
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
                if (mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("Expecting a geometrical mask.");

                QuadRule fullRule = RefElement.GetQuadratureRule(order);
                int L1 = fullRule.NoOfNodes;
                int D = fullRule.SpatialDim;

                var otherRule = m_orgrule.GetQuadRuleSet(mask, order);
                var ret = new List<IChunkRulePair<QuadRule>>(otherRule.Count());
                foreach (var x in otherRule) {

                    Chunk chk = x.Chunk;
                    QuadRule qr = x.Rule;
                    int L2 = qr.NoOfNodes;

                    Debug.Assert(qr.SpatialDim == fullRule.SpatialDim);

                    QuadRule compQr = new QuadRule();
                    compQr.OrderOfPrecision = qr.OrderOfPrecision;

                    compQr.Nodes = new NodeSet(this.RefElement, L1 + L2, D, true);
                    compQr.Weights = MultidimensionalArray.Create(L1 + L2);


                    compQr.Nodes.SetSubArray(fullRule.Nodes, new int[] { 0, 0 }, new int[] { L1 - 1, D - 1 });
                    compQr.Weights.SetSubArray(fullRule.Weights, new int[] { 0 }, new int[] { L1 - 1 });
                    compQr.Nodes.SetSubArray(qr.Nodes, new int[] { L1, 0 }, new int[] { L1 + L2 - 1, D - 1 });
                    compQr.Weights.AccSubArray(-1, qr.Weights, new int[] { L1 }, new int[] { L1 + L2 - 1 });

                    compQr.Nodes.LockForever();

                    ret.Add(new ChunkRulePair<QuadRule>(chk, compQr));
                }

                return ret;
            }

            public RefElement RefElement {
                get {
                    return m_orgrule.RefElement;
                }
            }
        }
    }
}
