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
    public interface ISayeGaussRule
    {
        RefElement RefElement { get; }
        int order { set; }
        QuadRule Evaluate(int cell);
    }

    public interface ISayeGaussComboRule
        : ISayeGaussRule
    {
        QuadRule[] ComboEvaluate(int cell);
    }

    public class SayeGaussComboRuleFactory
    {
        ISayeGaussComboRule comboRule;
        List<ChunkRulePair<QuadRule>>[] rulez;

        //Holds the factories instantiated by CalculateComboQuadRuleSet(...) 
        IQuadRuleFactory<QuadRule> surfaceRuleFactory;
        IQuadRuleFactory<QuadRule> volumeRuleFactory;

        class Status
        {
            public bool initialized;
            public int order;
            public ExecutionMask initialMask;
            public Mode RecalculationMode;
        }

        Status ComboStatus;

        public enum Mode { RecalculateOnEveryVolumeCall, RecalculateWhenInstanced};


        /// <summary>
        /// Calculates surface and volume quadrature rules in one step, which is faster when both rules 
        /// are needed. Before calling GetSurfaceRule() or GetVolumeRule() to receive the respective factories, call 
        /// CalculateComboQuadRuleSet(...).
        /// </summary>
        /// <param name="ComboRule"></param>
        public SayeGaussComboRuleFactory(ISayeGaussComboRule ComboRule, Mode recalculationMode)
        {
            comboRule = ComboRule;

            rulez = new[] {
                new List<ChunkRulePair<QuadRule>>(),
                new List<ChunkRulePair<QuadRule>>()
            };

            ComboStatus = new Status
            {
                initialized = false,
                order = 0,
                RecalculationMode = recalculationMode
            };
            
            volumeRuleFactory = new ComboFactoryWrapper(
                CalculateComboQuadRuleSet, 
                rulez[0], 
                comboRule.RefElement, 
                ComboStatus,
                ComboFactoryWrapper.QuadratureType.Volume);
            surfaceRuleFactory = new ComboFactoryWrapper(
                CalculateComboQuadRuleSet, 
                rulez[1], 
                comboRule.RefElement, 
                ComboStatus,
                ComboFactoryWrapper.QuadratureType.Surface);
        }

        /// <summary>
        /// Returns factory for surface rules
        /// </summary>
        /// <returns></returns>
        public IQuadRuleFactory<QuadRule> GetSurfaceFactory()
        {
            Console.WriteLine("Calling Surface Factory \n");
            return surfaceRuleFactory;
        }

        /// <summary>
        /// Returns factory for volume rules
        /// </summary>
        /// <returns></returns>
        public IQuadRuleFactory<QuadRule> GetVolumeFactory()
        {
            Console.WriteLine("Calling Volume Factory \n");
            return volumeRuleFactory;
        }

        /// <summary>
        /// Run this 
        /// </summary>
        /// <param name="mask"></param>
        /// <param name="order"></param>
        void CalculateComboQuadRuleSet(ExecutionMask mask, int order)
        {
            comboRule.order = order;
            rulez[0].Clear();
            rulez[1].Clear();
            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    QuadRule[] sayeRule = comboRule.ComboEvaluate(cell);
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
            public enum QuadratureType {Surface, Volume};

            QuadratureType quadMode; 

            Action<ExecutionMask, int> ComboRuleEvaluator;

            private IEnumerable<IChunkRulePair<QuadRule>> rule;

            private RefElement refElem;

            Status RuleStatus;

            public ComboFactoryWrapper(
                Action<ExecutionMask, int> comboRuleEvaluator,
                IEnumerable<IChunkRulePair<QuadRule>> Rule,
                RefElement RefElem,
                Status ruleStatus,
                QuadratureType quadType)
            {
                rule = Rule;
                refElem = RefElem;
                ComboRuleEvaluator = comboRuleEvaluator;
                RuleStatus = ruleStatus;
                quadMode = quadType;
            }

            public RefElement RefElement => refElem;

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
            {
                if(!RuleStatus.initialized || order != RuleStatus.order)
                {
                    Console.WriteLine("Initialize {0} \n", quadMode);
                    return InitializeRule(mask, order);
                }
                else
                {
                    switch (RuleStatus.RecalculationMode)
                    {
                        case Mode.RecalculateWhenInstanced:
                            Console.WriteLine("Reusing {0} \n", quadMode);
                            return UseExistingRule(mask);
                            //break;
                        case Mode.RecalculateOnEveryVolumeCall:
                            if (quadMode == QuadratureType.Volume)
                            {
                                Console.WriteLine("Recalcing {0} \n", quadMode);
                                return InitializeRule(mask, order);
                            }
                            else
                            {
                                Console.WriteLine("Reusing {0} \n", quadMode);
                                return UseExistingRule(mask);
                            }
                            //break;
                        default:
                            throw new NotImplementedException();
                            //break;
                    }
                }
            }

            IEnumerable<IChunkRulePair<QuadRule>> InitializeRule(ExecutionMask subMask, int Order)
            {
                ComboRuleEvaluator(subMask, Order);
                RuleStatus.order = Order;
                RuleStatus.initialMask = subMask;
                RuleStatus.initialized = true;
                return rule;
            }

            IEnumerable<IChunkRulePair<QuadRule>> UseExistingRule(ExecutionMask subMask)
            {
                if (subMask.Equals(RuleStatus.initialMask))
                {
                    return rule;
                }
                //If not, filter rules to fit subMask
                else
                {
                    Debug.Assert(subMask.IsSubMaskOf(RuleStatus.initialMask));
                    List<IChunkRulePair<QuadRule>> subRulez = new List<IChunkRulePair<QuadRule>>(subMask.Count());

                    IEnumerator<int> initialMask_enum = RuleStatus.initialMask.GetItemEnumerator();
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

            public int[] GetCachedRuleOrders()
            {
                throw new NotImplementedException();
            }

        }
    }

    class SayeGaussRuleFactory :
        IQuadRuleFactory<QuadRule>
    {
        ISayeGaussRule rule;

        public SayeGaussRuleFactory(ISayeGaussRule Rule)
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
            rule.order = order;
            var result = new List<ChunkRulePair<QuadRule>>();

            int number = 0;
            //Stopwatch stopWatch = new Stopwatch();
            //stopWatch.Start();
            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    QuadRule sayeRule = rule.Evaluate(cell);
                    ChunkRulePair<QuadRule> sayePair = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule);
                    result.Add(sayePair);
                    ++number;
                }
            }
            //stopWatch.Stop();
            //long ts = stopWatch.ElapsedMilliseconds;
            //Console.WriteLine("Number Of Cells " + number);
            //Console.WriteLine("RunTime " + ts + "ms");
            return result;
        }
    }

    /// <summary>
    /// Holds available factories for Saye quadrature.
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
            ISayeGaussRule rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Volume
                );
            return new SayeGaussRuleFactory(rule);
        }

        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 3D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_LevelSet3D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussRule rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Surface
                );
            return new SayeGaussRuleFactory(rule);
        }

        /// <summary>
        /// Gauss rules for \f$ \int_{\frakA \cap K_j } \ldots \dV \f$ in the 2D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Volume2D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussRule rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.PositiveVolume
                );
            return new SayeGaussRuleFactory(rule);
        }

        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 2D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_LevelSet2D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussRule rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.Surface
                );
            return new SayeGaussRuleFactory(rule);
        }

        #endregion

        #region Combo QuadRules

        public static SayeGaussComboRuleFactory SayeGaussRule_Combo2D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussComboRule rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.Combo
                );
            return new SayeGaussComboRuleFactory(
                rule, 
                SayeGaussComboRuleFactory.Mode.RecalculateWhenInstanced);
        }

        public static SayeGaussComboRuleFactory SayeGaussRule_Combo3D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussComboRule rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Combo
                );
            return new SayeGaussComboRuleFactory(
                rule, 
                SayeGaussComboRuleFactory.Mode.RecalculateOnEveryVolumeCall);
        }

        #endregion
    }
}
