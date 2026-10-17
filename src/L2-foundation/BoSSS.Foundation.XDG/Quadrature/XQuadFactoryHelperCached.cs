using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature;
using ilPSP;
using IntersectingQuadrature.Tensor;
using NUnit.Framework.Internal;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.XDG {


    /// <summary>
    /// Auxiliary class that which provides caching and thread-parallel creation of XDG-quadrature schemes 
    /// created by some other quadrature factory helper class (<see cref="XQuadFactoryHelperBase"/>)
    /// 
    /// The caching is actually achieved by classes 
    /// - <see cref="CachingQuadRuleFactory{T}"/> and
    /// - <see cref="EdgeRuleFromCellBoundaryFactory"/>
    /// 
    /// Note: in order to work properly, 
    /// 1. the function <see cref="CreateRulesAndMPIExchgange(int)"/> has to be called once and 
    /// 2. and all quadrature rule factories have to be executed once; this is achieved by initializing <see cref="CutCellMetrics"/>.
    /// </summary>
    public class XQuadFactoryHelperCached : IXQuadFactoryHelper {

        internal XQuadFactoryHelperCached(Func<int,XQuadFactoryHelperBase> originalFactoryFactory) {
            int NumThreads = ilPSP.Environment.NumThreads;
            m_OriginalFactories = new XQuadFactoryHelperBase[NumThreads];
            for(int i = 0; i < NumThreads; i++)
                m_OriginalFactories[i] = originalFactoryFactory(i);
            this.CutCellQuadratureType = m_OriginalFactories[0].CutCellQuadratureType;

            m_IntersectionRuleFactory = new CacheGeneric<(int levSetIndex0, int levSetIndex1, RefElement KrefVol), QuadRule>(
                (iThread, t3) => m_OriginalFactories[iThread].GetIntersectionRuleFactory(t3.levSetIndex0, t3.levSetIndex1, t3.KrefVol)
            );


            m_SurfaceFactory = new CacheGeneric<(int levSetIndex, RefElement KrefEdge), QuadRule>(
                (iThread, tt) => m_OriginalFactories[iThread].GetSurfaceFactory(tt.levSetIndex, tt.KrefEdge)
            );
            m_SurfaceFactory1 = new CacheGeneric<(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule>( 
                (iThread, t4) => m_OriginalFactories[iThread].GetSurfaceFactory(t4.levSetIndex0, t4.levSetIndex1, t4.jmp1, t4.KrefVol)
            );
            m_VolRuleFactory = new CacheGeneric<(int levSetIndex, JumpTypes jmp, RefElement Kref), QuadRule>(
                (iThread, t3) => m_OriginalFactories[iThread].GetVolRuleFactory(t3.levSetIndex, t3.jmp, t3.Kref)
            );
            m_VolRuleFactory1 = new CacheGeneric<(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule>(
                (iThread, t5) => m_OriginalFactories[iThread].GetVolRuleFactory(t5.levSetIndex0, t5.jmp0, t5.levSetIndex1, t5.jmp1, t5.KrefVol)
            );
            m_SurfaceElement_BoundaryRuleFactory1b = new CacheGeneric<(int levSetIndex, RefElement KrefVol), CellBoundaryQuadRule>(
                (iThread, tt) => m_OriginalFactories[iThread]._GetSurfaceElement_BoundaryRuleFactory(tt.levSetIndex, tt.KrefVol)
            );

            m_EdgeRuleFactory = new CacheGeneric<(int levSetIndex, JumpTypes jmp, RefElement KrefEdge), QuadRule>(
                (iThread, t3) => m_OriginalFactories[iThread].GetEdgeRuleFactory(t3.levSetIndex, t3.jmp, t3.KrefEdge)
            );
            m_EdgeRuleFactory1 = new CacheGeneric<(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule>(
                (iThread, t5) => m_OriginalFactories[iThread].GetEdgeRuleFactory(t5.levSetIndex0, t5.jmp0, t5.levSetIndex1, t5.jmp1, t5.KrefVol)
            );
            m_SurfaceElement_BoundaryRuleFactory = new CacheGeneric<(int levSetIndex, RefElement KrefEdge), QuadRule>(
                (iThread, tt) => m_OriginalFactories[iThread].GetSurfaceElement_BoundaryRuleFactory(tt.levSetIndex, tt.KrefEdge)
            );
            m_SurfaceElement_BoundaryRuleFactory1 = new CacheGeneric<(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule>(
                (iThread, t4) => m_OriginalFactories[iThread].GetSurfaceElement_BoundaryRuleFactory(t4.levSetIndex0, t4.levSetIndex1, t4.jmp1, t4.KrefVol)
            );
        }

        readonly XQuadFactoryHelperBase[] m_OriginalFactories;


        public CutCellQuadratureMethod CutCellQuadratureType {
            get;
            private set;
        }


        class CacheGeneric<Tkey,TQuadRule> where TQuadRule : QuadRule {

            public CacheGeneric(Func<int,Tkey,IQuadRuleFactory<TQuadRule>> factoryFactory) {
                m_FactoryFactory = factoryFactory;
            }

            readonly Func<int,Tkey,IQuadRuleFactory<TQuadRule>> m_FactoryFactory;

            Dictionary<Tkey, IQuadRuleFactory<TQuadRule>> m_QuadRuleFactory = new Dictionary<Tkey, IQuadRuleFactory<TQuadRule>>();

            public IQuadRuleFactory<TQuadRule> GetFactory(Tkey key) {
                if(!m_QuadRuleFactory.ContainsKey(key)) {
                    var fac0 = m_FactoryFactory(0, key);
                    if(fac0 is IQuadRuleFactory_ext<TQuadRule> fac0_ext) { 
                        if(fac0_ext.RuleIsCached)
                            return fac0; // bypass the caching, because `fac0` performs its own caching
                    }

                    m_QuadRuleFactory.Add(
                        key,
                        new CachingQuadRuleFactory<TQuadRule>(
                            (int iThread) => m_FactoryFactory(iThread, key)
                        )
                    );
                }
                return m_QuadRuleFactory[key];
            }
        }


        public void CreateRulesAndMPIExchgange(int __quadorder) {
            m_OriginalFactories[0].CreateRulesAndMPIExchgange(__quadorder);
        }



        readonly CacheGeneric<(int levSetIndex0, int levSetIndex1, RefElement KrefVol), QuadRule> m_IntersectionRuleFactory;
        public IQuadRuleFactory<QuadRule> GetIntersectionRuleFactory(int levSetIndex0, int levSetIndex1, RefElement KrefVol) {
            return m_IntersectionRuleFactory.GetFactory((levSetIndex0, levSetIndex1, KrefVol));
        }

        readonly CacheGeneric<(int levSetIndex, RefElement Kref), QuadRule> m_SurfaceFactory;
        public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex, RefElement Kref) {
            return m_SurfaceFactory.GetFactory((levSetIndex, Kref));
        }

        readonly CacheGeneric<(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule> m_SurfaceFactory1;
        public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol) {
            return m_SurfaceFactory1.GetFactory((levSetIndex0, levSetIndex1, jmp1, KrefVol));
        }

        readonly CacheGeneric<(int levSetIndex, JumpTypes jmp, RefElement Kref), QuadRule> m_VolRuleFactory;
        public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex, JumpTypes jmp, RefElement Kref) {
            return m_VolRuleFactory.GetFactory((levSetIndex, jmp, Kref));
        }

        readonly CacheGeneric<(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule> m_VolRuleFactory1;
        public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol) {
            return m_VolRuleFactory1.GetFactory((levSetIndex0, jmp0, levSetIndex1, jmp1, KrefVol));
        }

        readonly CacheGeneric<(int levSetIndex, RefElement KrefVol), CellBoundaryQuadRule> m_SurfaceElement_BoundaryRuleFactory1b;      
        public IQuadRuleFactory<CellBoundaryQuadRule> _GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefVol) {
            return m_SurfaceElement_BoundaryRuleFactory1b.GetFactory((levSetIndex, KrefVol));
        }

        readonly CacheGeneric<(int levSetIndex, JumpTypes jmp, RefElement KrefEdge), QuadRule> m_EdgeRuleFactory;
        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex, JumpTypes jmp, RefElement KrefEdge) {
            return m_EdgeRuleFactory.GetFactory((levSetIndex, jmp, KrefEdge));
        }

        readonly CacheGeneric<(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule> m_EdgeRuleFactory1;
        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol) {
            return m_EdgeRuleFactory1.GetFactory((levSetIndex0, jmp0, levSetIndex1, jmp1, KrefVol));
        }

        readonly CacheGeneric<(int levSetIndex, RefElement KrefEdge), QuadRule> m_SurfaceElement_BoundaryRuleFactory;
        public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex, RefElement KrefEdge) {
            return m_SurfaceElement_BoundaryRuleFactory.GetFactory((levSetIndex, KrefEdge));
        }

        readonly CacheGeneric<(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol), QuadRule> m_SurfaceElement_BoundaryRuleFactory1;
        public IQuadRuleFactory<QuadRule> GetSurfaceElement_BoundaryRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol) {
            return m_SurfaceElement_BoundaryRuleFactory1.GetFactory((levSetIndex0, levSetIndex1, jmp1, KrefVol));
        }
    }


}