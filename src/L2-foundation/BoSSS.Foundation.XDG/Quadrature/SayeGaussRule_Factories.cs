using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;

namespace BoSSS.Foundation.XDG.Quadrature
{
    interface ISayeGaussComboRule<S, T>
        where S : IPsi
        where T : SayeArgument<S>
    {
        RefElement RefElement { get; }
        int order { set; }
        QuadRule[] ComboEvaluate(int cell, T arg);
        T CreateStartSetup();
    }

    class SayeGaussComboFactory<S, T>
        where S : IPsi
        where T : SayeArgument<S>
    {
        ISayeGaussComboRule<S, T> comboRule;

        IQuadRuleFactory<QuadRule> surfaceRuleFactory;
        IQuadRuleFactory<QuadRule> volumeRuleFactory;

        int ID = int.MinValue;

        IEnumerable<IChunkRulePair<QuadRule>>[] comboQuadRuleSet = new IEnumerable<IChunkRulePair<QuadRule>>[2];  

        SayeGaussComboFactory(ISayeGaussComboRule<S,T> ComboRule)
        {
            comboRule = ComboRule;
            surfaceRuleFactory = new ComboFactory(this, ComboFactory.mode.surface);
            volumeRuleFactory = new ComboFactory(this, ComboFactory.mode.volume);

        }
        
        IQuadRuleFactory<QuadRule> GetSurfaceRule()
        {
            return surfaceRuleFactory;
        }

        IQuadRuleFactory<QuadRule> GetVolumeRule()
        {
            return volumeRuleFactory;
        }

        void CalculateComboQuadRuleSet(ExecutionMask mask, int order)
        {
            comboRule.order = 2 * order;
            List<ChunkRulePair<QuadRule>>[] result = new [] {
                new List<ChunkRulePair<QuadRule>>(),
                new List<ChunkRulePair<QuadRule>>()
            };

            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    T arg = comboRule.CreateStartSetup();
                    QuadRule[] sayeRule = comboRule.ComboEvaluate(cell, arg);
                    ChunkRulePair<QuadRule> sayePair_surface = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule[0]);
                    ChunkRulePair<QuadRule> sayePair_volume = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule[1]);
                    result[0].Add(sayePair_surface);
                    result[1].Add(sayePair_volume);
                }
            }
            comboQuadRuleSet = result;
        }

        class ComboFactory :
            IQuadRuleFactory<QuadRule>
        {
            public enum mode { surface, volume };
            SayeGaussComboFactory<S, T> factory;
            mode integrationMode;

            public ComboFactory(SayeGaussComboFactory<S, T> Factory, mode IntegrationMode)
            {
                factory = Factory;
                integrationMode = IntegrationMode;
            }

            public RefElement RefElement => factory.comboRule.RefElement;

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
            {
                int ID = GetID();
                if(ID != factory.ID)
                {
                    factory.CalculateComboQuadRuleSet(mask, order);
                    factory.ID = ID;
                }
                switch (integrationMode)
                {
                    case mode.surface:
                        return factory.comboQuadRuleSet[0];
                    case mode.volume:
                        return factory.comboQuadRuleSet[1];
                    default:
                        throw new NotSupportedException();
                }
            }

            int GetID()
            {
                return 0;
            }

            public int[] GetCachedRuleOrders()
            {
                throw new NotImplementedException();
            }
        }

    }

    interface ISayeGaussRule<S, T>
        where S : IPsi
        where T : SayeArgument<S>
    {
        RefElement RefElement { get; }
        int order { set; }
        QuadRule Evaluate(int cell, T arg);
        T CreateStartSetup();
    }

    class SayeGaussRuleFactory<S, T> :
        IQuadRuleFactory<QuadRule>
        where S : IPsi
        where T : SayeArgument<S>
    {
        ISayeGaussRule<S, T> rule;

        public SayeGaussRuleFactory(ISayeGaussRule<S, T> Rule)
        {
            rule = Rule;
        }

        public RefElement RefElement => rule.RefElement;

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            rule.order = 2 * order;
            var result = new List<ChunkRulePair<QuadRule>>();

            int number = 0;
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    T arg = rule.CreateStartSetup();
                    QuadRule sayeRule = rule.Evaluate(cell, arg);
                    ChunkRulePair<QuadRule> sayePair = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule);
                    result.Add(sayePair);
                    ++number;
                }
            }
            stopWatch.Stop();
            long ts = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Number Of Cells " + number);
            Console.WriteLine("RunTime " + ts + "ms");
            return result;
        }
    }

    public static class SayeFactories
    {
        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 3D case
        /// </summary>
        /// 
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Volume3D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder)
        {
            ISayeGaussRule<LinearPSI<Cube>, LinearSayeSpace<Cube>> rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Volume
                );
            return new SayeGaussRuleFactory<LinearPSI<Cube>, LinearSayeSpace<Cube>>(rule);
        }

        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 3D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_LevelSet3D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussRule<LinearPSI<Cube>, LinearSayeSpace<Cube>> rule = new SayeFactory_Cube(_lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Surface
                );
            return new SayeGaussRuleFactory<LinearPSI<Cube>, LinearSayeSpace<Cube>>(rule);
        }

        /// <summary>
        /// Gauss rules for \f$ \int_{\frakA \cap K_j } \ldots \dV \f$ in the 2D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Volume2D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder,
            SayeFactory_Square.QuadratureMode mode
            )
        {
            ISayeGaussRule<LinearPSI<Square>, SayeSquare> rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                mode
                );
            return new SayeGaussRuleFactory<LinearPSI<Square>, SayeSquare>(rule);
        }

        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 2D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_LevelSet2D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder)
        {
            ISayeGaussRule<LinearPSI<Square>, SayeSquare> rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.Surface
                );
            return new SayeGaussRuleFactory<LinearPSI<Square>, SayeSquare>(rule);
        }

    }
}
