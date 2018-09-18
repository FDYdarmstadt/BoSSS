using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;

namespace BoSSS.Foundation.XDG.Quadrature
{
    public interface ISayeGaussComboRule<S, T>
        : ISayeGaussRule<S, T>
        where S : IPsi
        where T : SayeArgument<S>
    {
        QuadRule[] ComboEvaluate(int cell, T arg);
    }

    public class SayeGaussComboRuleFactory<S, T>
        where S : IPsi
        where T : SayeArgument<S>
    {
        ISayeGaussComboRule<S, T> comboRule;
        List<ChunkRulePair<QuadRule>>[] rulez;

        //Holds the factories instantiated by CalculateComboQuadRuleSet(...) 
        IQuadRuleFactory<QuadRule> surfaceRuleFactory;
        IQuadRuleFactory<QuadRule> volumeRuleFactory;

        class Status
        {
            public bool initialized;
            public int order;
        }

        Status ComboStatus;

        /// <summary>
        /// Calculates surface and volume quadrature rules in one step, which is faster when both rules 
        /// are needed. Before calling GetSurfaceRule() or GetVolumeRule() to receive the respective factories, call 
        /// CalculateComboQuadRuleSet(...).
        /// </summary>
        /// <param name="ComboRule"></param>
        public SayeGaussComboRuleFactory(ISayeGaussComboRule<S,T> ComboRule)
        {
            comboRule = ComboRule;
            rulez = new[] {
                new List<ChunkRulePair<QuadRule>>(),
                new List<ChunkRulePair<QuadRule>>()
            };
            ComboStatus = new Status
            {
                initialized = false,
                order = 0
            };
            
            volumeRuleFactory = new ComboFactoryWrapper(CalculateComboQuadRuleSet, rulez[0], comboRule.RefElement, ComboStatus);
            surfaceRuleFactory = new ComboFactoryWrapper(CalculateComboQuadRuleSet, rulez[1], comboRule.RefElement, ComboStatus);

        }

        /// <summary>
        /// Returns factory for surface rules
        /// </summary>
        /// <returns></returns>
        public IQuadRuleFactory<QuadRule> GetSurfaceRule()
        {
            return surfaceRuleFactory;
        }

        /// <summary>
        /// Returns factory for volume rules
        /// </summary>
        /// <returns></returns>
        public IQuadRuleFactory<QuadRule> GetVolumeRule()
        {
            return volumeRuleFactory;
        }

        /// <summary>
        /// Run this 
        /// </summary>
        /// <param name="mask"></param>
        /// <param name="order"></param>
        void CalculateComboQuadRuleSet(ExecutionMask mask, int order)
        {
            comboRule.order = 2 * order;
            rulez[0].Clear();
            rulez[1].Clear();
            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    T arg = comboRule.CreateStartSetup();
                    QuadRule[] sayeRule = comboRule.ComboEvaluate(cell, arg);
                    ChunkRulePair<QuadRule> sayePair_volume = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule[0]);
                    ChunkRulePair<QuadRule> sayePair_surface = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule[1]);
                    rulez[0].Add(sayePair_volume);
                    rulez[1].Add(sayePair_surface);
                }
            }
        }

        //
        class ComboFactoryWrapper :
            IQuadRuleFactory<QuadRule>
        {
            Action<ExecutionMask, int> ComboRule;
            private IEnumerable<IChunkRulePair<QuadRule>> rule;
            private RefElement refElem;
            ExecutionMask initialMask;
            Status RuleStatus;

            public ComboFactoryWrapper(
                Action<ExecutionMask, int> comboRule, 
                IEnumerable<IChunkRulePair<QuadRule>> Rule, 
                RefElement RefElem, 
                Status ruleStatus)
            {
                rule = Rule;
                refElem = RefElem;
                ComboRule = comboRule;
                RuleStatus = ruleStatus;
            }

            public RefElement RefElement => refElem;

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask subMask, int Order)
            {
                if(!RuleStatus.initialized || Order != RuleStatus.order)
                {
                    ComboRule(subMask, Order);
                    RuleStatus.order = Order;
                    initialMask = subMask;
                    RuleStatus.initialized = true;
                    return rule;
                }
                else
                {
                    //Check if subMask is initialMask. If so, return all rules
                    if (subMask.Equals(initialMask))
                    {
                        return rule;
                    }
                    //If not, filter rules to fit subMask
                    else
                    {
                        Debug.Assert(subMask.IsSubMaskOf(initialMask));
                        List<IChunkRulePair<QuadRule>> subRulez = new List<IChunkRulePair<QuadRule>>(subMask.Count());

                        IEnumerator<int> initialMask_enum = initialMask.GetItemEnumerator();
                        int i = 0;
                        BitArray subMask_bitmask = subMask.GetBitMask();
                        while (initialMask_enum.MoveNext())
                        {
                            if (subMask_bitmask[initialMask_enum.Current])
                            {
                                subRulez.Add(rule.ElementAt(i));
                            }
                            ++i;
                        }
                        return subRulez;
                        
                    }
                }
            }

            public int[] GetCachedRuleOrders()
            {
                throw new NotImplementedException();
            }
        }

    }

    public interface ISayeGaussRule<S, T>
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

    /// <summary>
    /// Holds available factories for Saye quadrature. Idea: Implement as Singleton factory. 
    /// </summary>
    public static class SayeFactories
    {
        #region Single QuadRules
        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 3D case
        /// </summary>
        /// 
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Volume3D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder
            )
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
            ISayeGaussRule<LinearPSI<Cube>, LinearSayeSpace<Cube>> rule = new SayeFactory_Cube(
                _lsData,
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
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussRule<LinearPSI<Square>, SayeSquare> rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.Surface
                );
            return new SayeGaussRuleFactory<LinearPSI<Square>, SayeSquare>(rule);
        }

        #endregion

        #region Combo QuadRules

        public static SayeGaussComboRuleFactory<LinearPSI<Square>, SayeSquare> SayeGaussRule_Combo2D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussComboRule<LinearPSI<Square>, SayeSquare> rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.Combo
                );
            return new SayeGaussComboRuleFactory<LinearPSI<Square>, SayeSquare>(rule);
        }

        public static SayeGaussComboRuleFactory<LinearPSI<Cube>, LinearSayeSpace<Cube>> SayeGaussRule_Combo3D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussComboRule<LinearPSI<Cube>, LinearSayeSpace<Cube>> rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Combo
                );
            return new SayeGaussComboRuleFactory<LinearPSI<Cube>, LinearSayeSpace<Cube>>(rule);
        }

        #endregion
    }
}
