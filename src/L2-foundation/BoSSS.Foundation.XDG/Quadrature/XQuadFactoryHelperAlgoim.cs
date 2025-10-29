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
using BoSSS.Foundation.XDG.Quadrature.Algoim;
using IntersectingQuadrature;
using static BoSSS.Foundation.XDG.XQuadFactoryHelperBase;

namespace BoSSS.Foundation.XDG {




    /// <summary>
    /// Auxiliary class that helps with the creation of XDG-quadrature schemes for Algoim <see cref="AlgoimFactories"/>;
    /// instances can be obtained via <see cref="LevelSetTracker.GetXQuadFactoryHelper"/>.
    /// </summary>
    public class XQuadFactoryHelperAlgoim : XQuadFactoryHelperBase {

        Dictionary<RefElement, AlgoimDoubleCutFactories> DoubleCutFactories = new Dictionary<RefElement, AlgoimDoubleCutFactories>();

        public XQuadFactoryHelperAlgoim(LevelSetTracker.LevelSetData[] lsDatas) : base(lsDatas) {

            //there are some methods explicitly relying on this property
            this.CutCellQuadratureType = CutCellQuadratureMethod.Algoim;

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




        public override IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefVol) {
            CheckIfVolElement(KrefVol);
            bool negativeLevelSet = CheckJmp(jmp);
            var gdat = this.m_LevelSetDatas[levSetIndex].GridDat;
            var KrefEdge = gdat.iGeomEdges.EdgeRefElements.Single();

            var maxDom = this.m_CutEdges4LevelSet[levSetIndex];
            if(gdat.iGeomCells.RefElements.Length > 1)
                maxDom = maxDom.Intersect(gdat.Edges.GetEdges4RefElement(KrefEdge));

            var algoimFactory = new AlgoimFactories(m_LevelSetDatas[levSetIndex], KrefVol, negativeLevelSet);

            var r = new EdgeRuleFromCellBoundaryFactory(gdat, true,
                    KrefEdge,
                    [ algoimFactory.GetCellBoundaryVolumeFactory() ],
                    maxDom);

            return r;
        }

		public override IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
			CheckIfVolElement(KrefVol);
			var jumps = new JumpTypes[] { jmp0, jmp1 };

            if (!DoubleCutFactories.ContainsKey(KrefVol))
                DoubleCutFactories.Add(KrefVol, new AlgoimDoubleCutFactories(m_LevelSetDatas, KrefVol));

			var algoimFactory = DoubleCutFactories[KrefVol];
			var cellBoundaryFac = algoimFactory.GetCellBoundaryVolumeFactory(jumps);

            //var cutDom1 = m_LevelSetDatas[levSetIndex0].Region.GetCutCellMask4LevSet(levSetIndex0);
            //var cutDom2 = m_LevelSetDatas[levSetIndex1].Region.GetCutCellMask4LevSet(levSetIndex1);
            //var cutDom = cutDom1.Intersect(cutDom2).ToGeometicalMask();
            //var _cutDom = cutDom.Intersect(m_LevelSetDatas[levSetIndex0].GridDat.Cells.GetCells4Refelement(KrefVol));

            var maxdom = m_CutEdges4LevelSet[levSetIndex0].Intersect(m_CutEdges4LevelSet[levSetIndex1]);

            var KrefEdge = gdat.iGeomEdges.EdgeRefElements.Single();

			var r = new EdgeRuleFromCellBoundaryFactory(gdat, true,
		                                                KrefEdge,
                                                        [ cellBoundaryFac ],
		                                                maxdom);

			return r;
		}

        public override void CreateRulesAndMPIExchgange(int __quadorder) {
            throw new NotImplementedException();
        }


        public override IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            throw new NotImplementedException();
        }

        Quadrature.HMF.LineAndPointQuadratureFactory[] LineAndPoint_in2D = null;

        public override IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol) {
            var gdat = this.m_LevelSetDatas[levSetIndex].GridDat;
            int D = gdat.SpatialDimension;
            var KrefEdge = gdat.iGeomEdges.EdgeRefElements.Single();

            var maxDom = m_CutEdges4LevelSet[levSetIndex];
            if(gdat.iGeomCells.RefElements.Length > 1)
                maxDom = maxDom.Intersect(gdat.iGeomEdges.GetEdges4RefElement(KrefEdge));


            if(D == 2) {
                if(LineAndPoint_in2D == null)
                    LineAndPoint_in2D = new LineAndPointQuadratureFactory[this.m_LevelSetDatas.Length];

                if(LineAndPoint_in2D[levSetIndex] == null) {
                    LineAndPoint_in2D[levSetIndex] = new LineAndPointQuadratureFactory(
                        KrefVol,
                        this.m_LevelSetDatas[levSetIndex],
                        true);
                }

                //return LineAndPoint_in2D[levSetIndex].GetPointFactory();

                return new EdgeRuleFromCellBoundaryFactory(gdat, D > 2, KrefEdge,
                        [ LineAndPoint_in2D[levSetIndex].GetPointFactory() ],
                        maxDom);
            } else {
                CheckIfVolElement(KrefVol);
                

                var algoimFactory = new AlgoimFactories(m_LevelSetDatas[levSetIndex], KrefVol);

                var r = new EdgeRuleFromCellBoundaryFactory(gdat, D > 2, KrefEdge,
                        [ algoimFactory.GetCellBoundarySurfaceFactory() ],
                        maxDom);

                return r;
            }
        }

        public override IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
			CheckIfVolElement(KrefVol);
			JumpTypes[] jumps = new JumpTypes[] { jmp1, jmp1 };
			jumps[levSetIndex0] = JumpTypes.Implicit;

			if (!DoubleCutFactories.ContainsKey(KrefVol))
				DoubleCutFactories.Add(KrefVol, new AlgoimDoubleCutFactories(m_LevelSetDatas, KrefVol));

			var algoimFactory = DoubleCutFactories[KrefVol];
			var cellBoundaryFac = algoimFactory.GetCellBoundarySurfaceFactory(jumps);

			var cutDom1 = m_LevelSetDatas[levSetIndex0].Region.GetCutCellMask4LevSet(levSetIndex0);
			var cutDom2 = m_LevelSetDatas[levSetIndex1].Region.GetCutCellMask4LevSet(levSetIndex1);
			var cutDom = cutDom1.Intersect(cutDom2).ToGeometicalMask();
			var _cutDom = cutDom.Intersect(m_LevelSetDatas[levSetIndex0].GridDat.Cells.GetCells4Refelement(KrefVol));


            var KrefEdge = gdat.iGeomEdges.EdgeRefElements.Single();
			
            var maxDom = m_CutEdges4LevelSet[levSetIndex0].Intersect(m_CutEdges4LevelSet[levSetIndex1]);
            if(gdat.iGeomCells.RefElements.Length > 1)
                maxDom = maxDom.Intersect(gdat.iGeomEdges.GetEdges4RefElement(KrefEdge));


            var r = new EdgeRuleFromCellBoundaryFactory(gdat, gdat.SpatialDimension > 2,
                                                        KrefEdge,
                                                        [ cellBoundaryFac ],
                                                        maxDom);

			return r;
		}

        public override IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref) {
            CheckIfVolElement(Kref);
            var algoimFactory = new AlgoimFactories(m_LevelSetDatas[levSetIndex], Kref);
            return algoimFactory.GetSurfaceFactory(); 
        }

        public override IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol, IQuadRuleFactory<QuadRule> backupFactory) {
            CheckIfVolElement(KrefVol);
            JumpTypes[] jumps = new JumpTypes[] { jmp1, jmp1 };
            jumps[levSetIndex0] = JumpTypes.Implicit;

			if (!DoubleCutFactories.ContainsKey(KrefVol))
				DoubleCutFactories.Add(KrefVol, new AlgoimDoubleCutFactories(m_LevelSetDatas, KrefVol));

			var algoimFactory = DoubleCutFactories[KrefVol];
			return algoimFactory.GetSurfaceFactory(jumps);
        }

        public override IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref) {
            CheckIfVolElement(Kref);
            bool negativeLevelSet = CheckJmp(jmp);
            var algoimFactory = new AlgoimFactories(m_LevelSetDatas[levSetIndex], Kref, negativeLevelSet);
            return algoimFactory.GetVolumeFactory();
        }

        public override IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement Kref, IQuadRuleFactory<QuadRule> backupFactory) {
			CheckIfVolElement(Kref);
			JumpTypes[] jumps = new JumpTypes[] { jmp0, jmp1};

			if (!DoubleCutFactories.ContainsKey(Kref))
				DoubleCutFactories.Add(Kref, new AlgoimDoubleCutFactories(m_LevelSetDatas, Kref));

			var algoimFactory = DoubleCutFactories[Kref];
			return algoimFactory.GetVolumeFactory(jumps);
        }

        /// <summary>
        /// Checks if the jump type is implemented and determines if negativeVolume is desired
        /// </summary>
        /// <param name="jmp"></param>
        /// <exception cref="NotImplementedException"></exception>
        private bool CheckJmp(JumpTypes jmp) {
            if (jmp == JumpTypes.Heaviside)
                return false; //we are looking for positive level set values, i.e., ls(x)>0
            else if (jmp == JumpTypes.OneMinusHeaviside)
                return true; //we are looking for negative level set values, i.e., ls(x)<0
			else
                throw new NotImplementedException();
        }

        /// <summary>
        /// Checks if reference element exists
        /// </summary>
        /// <exception cref="ArgumentException"></exception>
        private void CheckIfVolElement(RefElement KrefVol) {
            if (!gdat.iGeomCells.RefElements.Contains(KrefVol, (a, b) => object.ReferenceEquals(a, b)))
                throw new ArgumentException($"Expecting a volume/cell reference element, but got {KrefVol.GetType().Name}");
        }

    }




}
