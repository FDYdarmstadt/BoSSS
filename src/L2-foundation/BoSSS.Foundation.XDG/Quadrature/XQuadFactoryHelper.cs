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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Foundation.XDG.Quadrature.Saye;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;


namespace BoSSS.Foundation.XDG {


    

    /// <summary>
    /// Auxiliary class that helps with the creation of XDG-quadrature schemes;
    /// instances can be obtained via <see cref="LevelSetTracker.GetXQuadFactoryHelper"/>.
    /// </summary>
    public class XQuadFactoryHelper : XQuadFactoryHelperBase {

        


        /// <summary>
        /// ctor.
        /// </summary>
        internal XQuadFactoryHelper(LevelSetTracker.LevelSetData[] lsDatas, CutCellQuadratureMethod momentFittingVariant) 
            : base(lsDatas) {

            if(momentFittingVariant == CutCellQuadratureMethod.Algoim)
                throw new ArgumentException("for Algoim, a dedicated factory (XQuadFactoryHelperAlgoim) must be used", nameof(momentFittingVariant));

            if(momentFittingVariant == CutCellQuadratureMethod.OneStepGauss 
                || momentFittingVariant == CutCellQuadratureMethod.OneStepGaussAndStokes 
                || momentFittingVariant == CutCellQuadratureMethod.TwoStepStokesAndGauss
                || momentFittingVariant == CutCellQuadratureMethod.ExactCircle) {
                //if(lsDatas.Length > 1)
                //    throw new ArgumentException($"cut cell quadrature method '{momentFittingVariant}' does not support more than one level-set.", nameof(momentFittingVariant));
                if(lsDatas[0].GridDat.SpatialDimension != 2)
                    throw new ArgumentException($"cut cell quadrature method '{momentFittingVariant}' only supports 2D.", nameof(momentFittingVariant));
            }

            if(momentFittingVariant == CutCellQuadratureMethod.Classic) {
                //if(lsDatas.Length > 1)
                //    throw new ArgumentException($"cut cell quadrature method '{momentFittingVariant}' does not support more than one level-set.", nameof(momentFittingVariant));
            }


            this.CutCellQuadratureType = momentFittingVariant;
            
            this.m_CutEdges4LevelSet = new EdgeMask[lsDatas.Length];
            for(int iLs = 0; iLs < lsDatas.Length; iLs++) {
                var cutCellSubGrid = lsDatas[iLs].Region.GetCutCellSubgrid4LevSet(iLs);
                var innerCut = cutCellSubGrid.InnerEdgesMask;
                var bndyCut = cutCellSubGrid.AllEdgesMask.Intersect((gdat as GridData).BoundaryEdges);

                m_CutEdges4LevelSet[iLs] = innerCut.Union(bndyCut).ToGeometicalMask();
            }
    }


        /// <summary>
        /// - index: level-set index
        /// </summary>
        EdgeMask[] m_CutEdges4LevelSet;


        /// <summary>
        /// Triggers the creation of all quadrature rules, so that later, the can be retrieved from the cache.
        /// This is supposed to give a better runtime behavior, since the creation of quadrature rules is quite expensive.
        /// Furthermore, in the creation of quadrature rules there are some parts where MPI communication is needed or makes things more robust.
        /// One example are hanging node on MPI boundaries; There, a quadrature rule might be difficult to be created right locally.
        /// </summary>
        override public void CreateRulesAndMPIExchgange(int __quadorder) {
            using (var tr = new FuncTrace()) {
                MPICollectiveWatchDog.WatchAtRelease(csMPI.Raw._COMM.WORLD);


                // populate the caches...
                foreach(var KrefEdge in this.gdat.iGeomEdges.EdgeRefElements) {
                    foreach(var levSetIndex in m_LevelSetDatas.Select(ls => ls.LevelSetIndex)) {
                        this.GetSurfaceElement_BoundaryRuleFactory(levSetIndex, KrefEdge);
                        this.GetEdgeRuleFactory(levSetIndex, JumpTypes.Heaviside, KrefEdge);
                        this.GetEdgeRuleFactory(levSetIndex, JumpTypes.OneMinusHeaviside, KrefEdge);
                    }
                }
                
                // perform the MPI exchange
                var allFactories = m_SurfaceElement_BoundaryRuleFactory.Values.ToArray();
                allFactories = allFactories.Cat(m_EdgeRuleFactory1.Values);
                foreach (EdgeRuleFromCellBoundaryFactory f in  allFactories) {
                    f.CreateRulesAndMPIExchgange(__quadorder);
                }
            }
        }



        Quadrature.BruteForce.MultiLevelSetBruteForceQuadratureFactory zwoLSBruteForceFactories;
        //public Quadrature.Beck.MultiLevelSetBeckFactoryCreator zwoLSSayeFactories { get; private set; }


        // -----------------------------------------------------
        // Factory creation

        /// <summary>
        /// - 1st index: level set index
        /// - 2nd index: reference element index (index into <see cref="IGeometricalCellsData.RefElements"/>)
        /// </summary>
        Quadrature.HMF.LineAndPointQuadratureFactory[,] LineAndPoint_in2D = null;
        
        
        IQuadRuleFactory<CellBoundaryQuadRule>[,] CellFaceVolume_in3D = null;
        IQuadRuleFactory<CellBoundaryQuadRule>[,] CellFaceSurface_in3D = null;

        /// <summary>
        /// Returns a rule for the edges of surface-elements (elements on the zero-level-set surface, 
        /// i.e., on
        /// ```math 
        ///     \partial K \cap \mathfrak{I} .
        /// ```
        /// (point integrals in 2D, Line integrals in 3D)
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="CellBoundaryQuadRule"/>'s
        /// </returns>
        public override IQuadRuleFactory<CellBoundaryQuadRule> _GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol) {
            int D = gdat.SpatialDimension;
            int iKref = Array.IndexOf(gdat.iGeomCells.RefElements, KrefVol);
            int NoOfRefElements = gdat.iGeomCells.RefElements.Length;

            if (D == 2) {

                if (LineAndPoint_in2D == null)
                    LineAndPoint_in2D = new LineAndPointQuadratureFactory[this.m_LevelSetDatas.Length, NoOfRefElements];

                if (LineAndPoint_in2D[levSetIndex, iKref] == null) {
                    LineAndPoint_in2D[levSetIndex, iKref] = new LineAndPointQuadratureFactory(
                        KrefVol,
                        this.m_LevelSetDatas[levSetIndex],
                        true);// CutCellQuadratureType == MomentFittingVariants.OneStepGaussAndStokes);

                }

                return LineAndPoint_in2D[levSetIndex, iKref].GetPointFactory();
            } else {
                //throw new NotImplementedException("3d is not implemented yet");
                Debug.Assert(LineAndPoint_in2D == null);

                if (CellFaceSurface_in3D == null)
                    CellFaceSurface_in3D = new IQuadRuleFactory<CellBoundaryQuadRule>[this.m_LevelSetDatas.Length, NoOfRefElements];
                if (CellFaceSurface_in3D[levSetIndex, iKref] == null) {
                    var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

                    switch (CutCellQuadratureType) {
                        case CutCellQuadratureMethod.Saye: {
                            CellFaceSurface_in3D[levSetIndex, iKref] = SayeFactories.SayeGaussRule_EdgeSurface3D(
                                this.m_LevelSetDatas[levSetIndex],
                                rootFindingAlgorithm);
                            //CellFaceSurface_in3D[levSetIndex] = SayeFactories.SayeGaussRule_EdgeSurface3D(
                            //    this.m_LevelSetDatas[levSetIndex],
                            //    rootFindingAlgorithm);
                            break;
                        }
                        default: {
                            var CoFaceQuadRuleFactory = new CutLineOnEdgeQuadRuleFactory(
                                this.m_LevelSetDatas[levSetIndex],
                                rootFindingAlgorithm,
                                JumpTypes.Heaviside);
                            CellFaceSurface_in3D[levSetIndex, iKref] = new LevelSetEdgeSurfaceQuadRuleFactory(
                                this.m_LevelSetDatas[levSetIndex],
                                CoFaceQuadRuleFactory,
                                JumpTypes.Heaviside);
                            //new LevelSetEdgeVolumeQuadRuleFactory(
                            //    lsTrk, levSetIndex, rootFindingAlgorithm, JumpTypes.Heaviside);
                            break;
                        }
                    }

                }
                return CellFaceSurface_in3D[levSetIndex, iKref];

            }
        }

        Dictionary<(int iLevSet, int iKref), EdgeRuleFromCellBoundaryFactory> m_SurfaceElement_BoundaryRuleFactory = new Dictionary<(int iLevSet, int iKref), EdgeRuleFromCellBoundaryFactory>();



        /// <summary>
        /// Returns a rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e., for integrals 
        /// ```math 
        ///    \int_{ E \cap \mathfrak{I} \textrm{dl} . 
        /// ```
        /// This are point integrals in 2D and line integrals in 3D.
        /// 
        /// Internally, this method uses the <see cref="EdgeRuleFromCellBoundaryFactory"/>
        /// to split up a cell boundary quadrature for $`     \int_{ \partial K \cap \mathfrak{I} \textrm{dl}  `$ 
        /// into the respective edges.
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="QuadRule"/>'s on edges
        /// </returns>
        override public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefEdge) {
            int iKref = Array.IndexOf(gdat.iGeomEdges.EdgeRefElements, KrefEdge);
            if(iKref < 0)
                throw new ArgumentException("Expecting an **Edge** reference element.");

            var key = (levSetIndex, iKref);
            if(!m_SurfaceElement_BoundaryRuleFactory.ContainsKey(key)) {

                var gdat = this.m_LevelSetDatas[levSetIndex].GridDat;
                int D = gdat.SpatialDimension;

                var maxDom = m_CutEdges4LevelSet[levSetIndex];
                if(gdat.iGeomCells.RefElements.Length > 1)
                    maxDom = maxDom.Intersect(gdat.Edges.GetEdges4RefElement(KrefEdge));


                m_SurfaceElement_BoundaryRuleFactory.Add((levSetIndex, iKref),
                    new EdgeRuleFromCellBoundaryFactory(gdat, D > 2,
                                                        KrefEdge,
                                                        gdat.iGeomCells.RefElements.Select(KrefVol => _GetSurfaceElement_BoundaryRuleFactory(levSetIndex, KrefVol)),
                                                        maxDom));
            }

            return m_SurfaceElement_BoundaryRuleFactory[key];
        }

        /// <summary>
        /// Returns a rule factory for the boundary of surface-elements 
        /// (elements on the zero-level-set surface), i.e. on \f$  K \cap \mathfrak{I}\f$ .
        /// This are point integrals in 2D and line integrals in 3D.
        /// </summary>
        /// <returns>
        /// the returned factory produces <see cref="QuadRule"/>'s on edges
        /// </returns>
        public override IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            switch (CutCellQuadratureType) {
                case CutCellQuadratureMethod.Saye:
                    //var r = GetSurfaceElement_BoundaryRuleFactory_algoim(levSetIndex0, levSetIndex1, jmp1, KrefVol, backupFactory);
                    //return r;
                    return Quadrature.Intersecting.IntersectingQuadratureFactories.EdgePoint(m_LevelSetDatas[levSetIndex0], m_LevelSetDatas[levSetIndex1], jmp1);

                // rem: for the surface element:
                //return IntersectingQuadratureFactories.Surface(m_LevelSetDatas[levSetIndex0], m_LevelSetDatas[levSetIndex1], jmp1);
                default:
                    if (zwoLSBruteForceFactories == null) {
                        zwoLSBruteForceFactories = new Quadrature.BruteForce.MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetEdgePointRuleFactory(levSetIndex0, levSetIndex1, jmp1, backupFactory);
            }
        }

        /// <summary>
        /// Quadrature rule on cell boundaries
        /// </summary>
        public IQuadRuleFactory<CellBoundaryQuadRule> GetCellFaceFactory(int levSetIndex, RefElement Kref, JumpTypes jumpType) {
            int D = gdat.SpatialDimension;
            //int D = this.m_LevelSetDatas[0].GridDat.SpatialDimension;
            int iKref = Array.IndexOf(gdat.iGeomCells.RefElements, Kref);
            int NoOfKref = gdat.iGeomCells.RefElements.Length;

            if (D == 2) {
                if (jumpType != JumpTypes.Heaviside && jumpType != JumpTypes.OneMinusHeaviside)
                    throw new NotSupportedException();
                Debug.Assert(CellFaceVolume_in3D == null);

                if (LineAndPoint_in2D == null)
                    LineAndPoint_in2D = new LineAndPointQuadratureFactory[this.m_LevelSetDatas.Length, gdat.iGeomCells.RefElements.Length];

                if (LineAndPoint_in2D[levSetIndex, iKref] == null) {
                    LineAndPoint_in2D[levSetIndex, iKref] = new LineAndPointQuadratureFactory(Kref, this.m_LevelSetDatas[levSetIndex], true);
                }

                return LineAndPoint_in2D[levSetIndex, iKref].GetLineFactory(jumpType == JumpTypes.Heaviside ? true : false);
            } else if (D == 3) {
                Debug.Assert(LineAndPoint_in2D == null);
                if (CellFaceVolume_in3D == null)
                    CellFaceVolume_in3D = new IQuadRuleFactory<CellBoundaryQuadRule>[this.m_LevelSetDatas.Length, NoOfKref];
                if (jumpType != JumpTypes.Heaviside)
                    throw new NotSupportedException();

                if (CellFaceVolume_in3D[levSetIndex, iKref] == null) {
                    var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);
                    switch (CutCellQuadratureType) {
                        case CutCellQuadratureMethod.Saye: {
                            if(CellFaceVolume_in3D == null)
                                CellFaceVolume_in3D = new SayeGaussEdgeRuleFactory[this.m_LevelSetDatas.Length, NoOfKref];
                            CellFaceVolume_in3D[levSetIndex, iKref] = SayeFactories.SayeGaussRule_EdgeVolume3D(
                                this.m_LevelSetDatas[levSetIndex], rootFindingAlgorithm);
                            //if (CellFaceVolume_in3D == null)
                            //    CellFaceVolume_in3D = new SayeGaussEdgeRuleFactory[this.m_LevelSetDatas.Length];
                            //CellFaceVolume_in3D[levSetIndex] = SayeFactories.SayeGaussRule_EdgeVolume3D(
                            //    this.m_LevelSetDatas[levSetIndex], rootFindingAlgorithm);
                            break;
                        }

                        default:
                            if (CellFaceVolume_in3D == null)
                                CellFaceVolume_in3D = new LevelSetEdgeVolumeQuadRuleFactory[this.m_LevelSetDatas.Length, NoOfKref];
                            CellFaceVolume_in3D[levSetIndex, iKref] = new LevelSetEdgeVolumeQuadRuleFactory(
                                this.m_LevelSetDatas[levSetIndex], rootFindingAlgorithm, JumpTypes.Heaviside);
                            break;
                    }
                }
                return CellFaceVolume_in3D[levSetIndex, iKref];
            } else {
                throw new NotSupportedException();
            }
        }

        private void CheckJmp(JumpTypes jmp) {
            if (jmp != JumpTypes.Heaviside && jmp != JumpTypes.OneMinusHeaviside)
                throw new NotSupportedException();
        }


        Dictionary<(int iLevelSet, JumpTypes jmp, int iKref),EdgeRuleFromCellBoundaryFactory> m_EdgeRuleFactory1 = new Dictionary<(int iLevelSet,JumpTypes, int), EdgeRuleFromCellBoundaryFactory>();

        /// <summary>
        /// Generates a quadrature rule factory for the cut edge integrals.
        /// </summary>
        public override IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefEdge) {


            var gdat = this.m_LevelSetDatas[levSetIndex].GridDat;
            int D = gdat.SpatialDimension;

            int iKref = Array.IndexOf(gdat.iGeomEdges.EdgeRefElements, KrefEdge);
            if(iKref < 0)
                throw new ArgumentException($"Expecting an **Edge** reference element, but got {KrefEdge}.");

            CheckJmp(jmp);

            if (D == 2) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // In 2D, we use the `EdgeRuleFromCellBoundaryFactory` for both `JumpTypes.Heaviside` and `JumpTypes.OneMinusHeaviside`, i.e., for both sub-domains
                // (this means, for both sub-domains, we have nodes which are within the sub-domain and therefore only positive weights, at least for Saye-Rules)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                var key = (levSetIndex, jmp, iKref);

                if (!m_EdgeRuleFactory1.ContainsKey(key)) {
                    var maxDom = this.m_CutEdges4LevelSet[levSetIndex];
                    if (gdat.iGeomCells.RefElements.Length > 1)
                        maxDom = maxDom.Intersect(gdat.Edges.GetEdges4RefElement(KrefEdge));

                
                    m_EdgeRuleFactory1.Add(key, new EdgeRuleFromCellBoundaryFactory(gdat, true,
                                                                                    KrefEdge,
                                                                                    gdat.iGeomCells.RefElements.Select( KrefVol => GetCellFaceFactory(levSetIndex, KrefVol, jmp)),
                                                                                    maxDom));
                    //var r = new EdgeRuleFromCellBoundaryFactory(gdat, true,
                    //    GetCellFaceFactory(levSetIndex, KrefVol, jmp),
                    //    m_LevelSetDatas[levSetIndex].Region.GetCutCellMask4LevSet(levSetIndex).ToGeometicalMask());
                }
                return m_EdgeRuleFactory1[key];
            } else {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // In 3D, we use the `ComplementaryRuleFactory` for `JumpTypes.OneMinusHeaviside`
                // (this means, we have negative weights and nodes outside of the sub-domain;
                // this can be a problem, e.g. for compressible methods.
                // maybe, such cases have not been tried yet)
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                if (jmp == JumpTypes.Heaviside) {
                    var key = (levSetIndex, JumpTypes.Heaviside, iKref);

                    if (!m_EdgeRuleFactory1.ContainsKey(key)) {
                        var maxDom = this.m_CutEdges4LevelSet[levSetIndex];
                        if (gdat.iGeomCells.RefElements.Length > 1)
                            maxDom = maxDom.Intersect(gdat.Edges.GetEdges4RefElement(KrefEdge));

                        m_EdgeRuleFactory1.Add(key, new EdgeRuleFromCellBoundaryFactory(gdat, true,
                                                                                        KrefEdge,
                                                                                        gdat.iGeomCells.RefElements.Select(KrefVol => GetCellFaceFactory(levSetIndex, KrefVol, jmp)),
                                                                                        maxDom));
                    }
                    return m_EdgeRuleFactory1[key];
                    //var r = new EdgeRuleFromCellBoundaryFactory(gdat, true,
                    //    GetCellFaceFactory(levSetIndex, KrefVol, JumpTypes.Heaviside),
                    //    m_LevelSetDatas[levSetIndex].Region.GetCutCellMask4LevSet(levSetIndex).ToGeometicalMask());
                    //return r;

                } else if (jmp == JumpTypes.OneMinusHeaviside) {

                    return new ComplementaryRuleFactory(GetEdgeRuleFactory(levSetIndex, JumpTypes.Heaviside, KrefEdge));
                } else
                    throw new ArgumentOutOfRangeException("unsupported jump type");
            }
        }

        /// <summary>
        /// Generates an edge quadrature rule factory for edges cut by two level sets.
        /// </summary>
        public override IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            switch (CutCellQuadratureType) {
                case CutCellQuadratureMethod.Saye: {
                    return Quadrature.Intersecting.IntersectingQuadratureFactories.Edge(m_LevelSetDatas[levSetIndex0], jmp0, m_LevelSetDatas[levSetIndex1], jmp1);
                }
                default: {
                    if(zwoLSBruteForceFactories == null) {
                        zwoLSBruteForceFactories = new Quadrature.BruteForce.MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetEdgeRuleFactory(levSetIndex0, jmp0, levSetIndex1, jmp1, backupFactory);
                }
            }
        }

        /// <summary>
        /// Generates a quadrature rule factory for the cut volume integrals.
        /// </summary>
        override public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref) {
            CheckJmp(jmp);
            var ctx = this.m_LevelSetDatas[levSetIndex].GridDat;

            if (jmp == JumpTypes.Heaviside) {
                if (m_SurfaceFactory == null)
                    m_SurfaceFactory = new IQuadRuleFactory<QuadRule>[m_LevelSetDatas.Length];
                if (m_VolumeFactory == null)
                    m_VolumeFactory = new IQuadRuleFactory<QuadRule>[m_LevelSetDatas.Length];

                if (m_VolumeFactory[levSetIndex] == null) {
                    switch (CutCellQuadratureType) {
                        case CutCellQuadratureMethod.Classic: {
                            m_VolumeFactory[levSetIndex] = new LevelSetVolumeQuadRuleFactory(
                                this.m_LevelSetDatas[levSetIndex],
                                GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                                GetSurfaceFactory(levSetIndex, Kref),
                                jumpType: jmp);
                            break;
                        }

                        case CutCellQuadratureMethod.OneStepGauss:
                        case CutCellQuadratureMethod.OneStepGaussAndStokes: {
                                bool bStokes = CutCellQuadratureType == CutCellQuadratureMethod.OneStepGaussAndStokes;

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

                        case CutCellQuadratureMethod.TwoStepStokesAndGauss:
                        case CutCellQuadratureMethod.ExactCircle: {

                                m_VolumeFactory[levSetIndex] = (new LevelSetVolumeQuadRuleFactory2b(Kref,
                                        this.m_LevelSetDatas[levSetIndex],
                                        GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                                        GetSurfaceFactory(levSetIndex, Kref),
                                        jmp));
                                break;
                            }
                        case CutCellQuadratureMethod.Saye:
                            var comboFactory = Quadrature.Saye.SayeFactories.SayeGaussRule_Combo(
                                this.m_LevelSetDatas[levSetIndex],
                                new LineSegment.SafeGuardedNewtonMethod(1e-14));
                            m_VolumeFactory[levSetIndex] = comboFactory.GetVolumeFactory();
                            m_SurfaceFactory[levSetIndex] = comboFactory.GetSurfaceFactory();
                            break;
                        //case CutCellQuadratureMethod.Algoim:
                        //    var algoimComboFactory = new Quadrature.Algoim.AlgoimFactories(
                        //            this.m_LevelSetDatas[levSetIndex],
                        //            Kref);
                        //    m_VolumeFactory[levSetIndex] = algoimComboFactory.GetVolumeFactory();
                        //    m_SurfaceFactory[levSetIndex] = algoimComboFactory.GetSurfaceFactory();
                        //    break;
                        default:
                            throw new NotSupportedException(String.Format(
                                "Variant {0} not implemented.", CutCellQuadratureType));
                    }
                }

                Debug.Assert(m_VolumeFactory[levSetIndex] != null);
                return m_VolumeFactory[levSetIndex];
            } else if (jmp == JumpTypes.OneMinusHeaviside) {
                IQuadRuleFactory<QuadRule> ret;
                switch(CutCellQuadratureType) {
                    case CutCellQuadratureMethod.Saye: {

                        ret = Quadrature.Saye.SayeFactories.SayeGaussRule_NegativeVolume(this.m_LevelSetDatas[levSetIndex],
                                new LineSegment.SafeGuardedNewtonMethod(1e-14));
                        break;
                    }
                    //case CutCellQuadratureMethod.Algoim: {
                    //    var algoimComboFactory = new Quadrature.Algoim.AlgoimFactories(
                    //            this.m_LevelSetDatas[levSetIndex],
                    //            Kref,
                    //            true);
                    //    ret = algoimComboFactory.GetVolumeFactory();
                    //    break;
                    //}
                    default: {
                        ret = new ComplementaryRuleFactory(GetVolRuleFactory(levSetIndex, JumpTypes.Heaviside, Kref));
                        break;
                    }
                }
                Debug.Assert(ret != null);
                return ret;
            } else {
                throw new ArgumentOutOfRangeException("unsupported jump type");
            }
        }


        /// <summary>
        /// Generates a volume quadrature rule factory for cells cut by two level sets.
        /// </summary>

        public override IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            switch(CutCellQuadratureType) {
                case CutCellQuadratureMethod.Saye: {
                    return Quadrature.Intersecting.IntersectingQuadratureFactories.Volume(m_LevelSetDatas[levSetIndex0], jmp0, m_LevelSetDatas[levSetIndex1], jmp1);
                    //if(zwoLSSayeFactories == null) {
                    //    zwoLSSayeFactories = new Quadrature.Beck.MultiLevelSetBeckFactoryCreator(m_LevelSetDatas);
                    //}
                    //return zwoLSSayeFactories.GetVolRuleFactory(levSetIndex0, jmp0, levSetIndex1, jmp1, backupFactory);
                }
                default: {
                    if(zwoLSBruteForceFactories == null) {
                        zwoLSBruteForceFactories = new Quadrature.BruteForce.MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetVolRuleFactory(levSetIndex0, jmp0, levSetIndex1, jmp1, backupFactory);
                }
            }
        }

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
                    case CutCellQuadratureMethod.Classic: {

                        m_SurfaceFactory[levSetIndex] = new LevelSetSurfaceQuadRuleFactory(
                             m_LevelSetDatas[levSetIndex],
                             GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside));
                        break;
                    }
                    case CutCellQuadratureMethod.OneStepGauss:
                    case CutCellQuadratureMethod.OneStepGaussAndStokes: {
                            bool bStokes = CutCellQuadratureType == CutCellQuadratureMethod.OneStepGaussAndStokes;
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
                                        
                    case CutCellQuadratureMethod.TwoStepStokesAndGauss: {
                        m_SurfaceFactory[levSetIndex] = (new SurfaceStokes_2D(
                            m_LevelSetDatas[levSetIndex],
                            this.GetCellFaceFactory(levSetIndex, Kref, JumpTypes.Heaviside),
                            this._GetSurfaceElement_BoundaryRuleFactory(levSetIndex, Kref),
                            _SurfaceNodesOnZeroLevset: false,
                            _DoCheck: CheckQuadRules)).GetSurfaceFactory();
                        break;
                    }
                    case CutCellQuadratureMethod.ExactCircle: {
                        return new ExactCircleLevelSetIntegration(levSetIndex, this.m_LevelSetDatas[levSetIndex].GridDat, Kref);
                    }
                    case CutCellQuadratureMethod.Saye: {
                        var comboFactory = Quadrature.Saye.SayeFactories.SayeGaussRule_Combo(
                                this.m_LevelSetDatas[levSetIndex],
                                new LineSegment.SafeGuardedNewtonMethod(1e-14));
                        m_VolumeFactory[levSetIndex] = comboFactory.GetVolumeFactory();
                        m_SurfaceFactory[levSetIndex] = comboFactory.GetSurfaceFactory();
                        break;
                    }
                    //case CutCellQuadratureMethod.Algoim: {
                    //    var algoimComboFactory = new Quadrature.Algoim.AlgoimFactories(
                    //            this.m_LevelSetDatas[levSetIndex],
                    //            Kref);
                    //    m_VolumeFactory[levSetIndex] = algoimComboFactory.GetVolumeFactory();
                    //    m_SurfaceFactory[levSetIndex] = algoimComboFactory.GetSurfaceFactory();
                    //    break;
                    //}
                    default: {
                        throw new NotSupportedException(String.Format(
                            "Variant {0} not implemented.", CutCellQuadratureType));
                    }
                }
            }

            return m_SurfaceFactory[levSetIndex];
        }

        /// <summary>
        /// Generates a quadrature rule factory for integrating over a surface.
        /// The surface is defined by two conditions: <paramref name="levSetIndex0"/> = 0 and on side <paramref name="jmp1"/> of <paramref name="levSetIndex1"/>
        /// </summary>
        public override IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            switch(CutCellQuadratureType) {
                case CutCellQuadratureMethod.Saye: {
                    return Quadrature.Intersecting.IntersectingQuadratureFactories.Surface(m_LevelSetDatas[levSetIndex0], m_LevelSetDatas[levSetIndex1], jmp1);
                }
                default: {
                    if(zwoLSBruteForceFactories == null) {
                        zwoLSBruteForceFactories = new Quadrature.BruteForce.MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetSurfaceFactory(levSetIndex0,
                        levSetIndex1,
                        jmp1, backupFactory);
                }
            }
        }


        /// <summary>
        /// Generates a quadrature rule factory the intersection of <paramref name="levSetIndex0"/>-th and <paramref name="levSetIndex1"/>-th level-set.
        /// This is a point in 2D, a line in 3D.
        /// </summary>
        public override IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            switch (CutCellQuadratureType) {
                case CutCellQuadratureMethod.Saye:
                    return Quadrature.Intersecting.IntersectingQuadratureFactories.Intersection(m_LevelSetDatas[levSetIndex0], m_LevelSetDatas[levSetIndex1]);
                default:
                    if (zwoLSBruteForceFactories == null) {
                        zwoLSBruteForceFactories = new Quadrature.BruteForce.MultiLevelSetBruteForceQuadratureFactory(m_LevelSetDatas);
                    }
                    return zwoLSBruteForceFactories.GetIntersectionFactory(levSetIndex0, levSetIndex1, backupFactory);
            }
        }



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

                    compQr.Nodes.SetSubArray(fullRule.Nodes, [0, 0], [L1 - 1, D - 1]);
                    compQr.Weights.SetSubArray(fullRule.Weights, [0], [L1 - 1]);
                    compQr.Nodes.SetSubArray(qr.Nodes, [L1, 0], [L1 + L2 - 1, D - 1]);
                    compQr.Weights.AccSubArray(-1, qr.Weights, [L1], [L1 + L2 - 1]);

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
