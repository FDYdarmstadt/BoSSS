#undef LOG_TIME

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;


namespace BoSSS.Foundation.XDG.Quadrature {
    public class SayeGaussComboRuleFactory {
        ISayeGaussComboRule comboRule;
        List<ChunkRulePair<QuadRule>>[] rulez;

        //Holds the factories instantiated by CalculateComboQuadRuleSet(...) 
        IQuadRuleFactory<QuadRule> surfaceRuleFactory;
        IQuadRuleFactory<QuadRule> volumeRuleFactory;

        class Status {
            public bool Initialized;
            public int Order;
            public ExecutionMask MaxGrid;
            public RefElement ReferenceElement;
        }

        Status ComboStatus;


        /// <summary>
        /// Calculates surface and volume quadrature rules in one step, which is faster when both rules 
        /// are needed. Before calling GetSurfaceRule() or GetVolumeRule() to receive the respective factories, call 
        /// CalculateComboQuadRuleSet(...).
        /// </summary>
        public SayeGaussComboRuleFactory(ISayeGaussComboRule ComboRule, CellMask maxGrid) {
            comboRule = ComboRule;

            rulez = new[] {
                new List<ChunkRulePair<QuadRule>>(),
                new List<ChunkRulePair<QuadRule>>()
            };

            ComboStatus = new Status {
                Initialized = false,
                Order = -1,
                MaxGrid = maxGrid,
                ReferenceElement = comboRule.RefElement
            };

            volumeRuleFactory = new SayeFactoryWrapper(
                CalculateComboQuadRuleSet,
                rulez[0],
                ComboStatus);
            surfaceRuleFactory = new SayeFactoryWrapper(
                CalculateComboQuadRuleSet,
                rulez[1],
                ComboStatus);
        }

        /// <summary>
        /// Returns factory for surface rules
        /// </summary>
        /// <returns></returns>
        public IQuadRuleFactory<QuadRule> GetSurfaceFactory() {
            return surfaceRuleFactory;
        }

        /// <summary>
        /// Returns factory for volume rules
        /// </summary>
        /// <returns></returns>
        public IQuadRuleFactory<QuadRule> GetVolumeFactory() {
            return volumeRuleFactory;
        }

        /// <summary>
        /// Run this 
        /// </summary>
        /// <param name="mask"></param>
        /// <param name="order"></param>
        void CalculateComboQuadRuleSet(ExecutionMask mask, int order) {
            comboRule.order = order;
            rulez[0].Clear();
            rulez[1].Clear();
            //Find quadrature nodes and weights in each cell/chunk
#if LOG_TIME
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
#endif
            foreach (Chunk chunk in mask) {
                foreach (int cell in chunk.Elements) {
                    QuadRule[] sayeRule = comboRule.ComboEvaluate(cell);
                    ChunkRulePair<QuadRule> sayePair_volume = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule[0]);
                    ChunkRulePair<QuadRule> sayePair_surface = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule[1]);
                    rulez[0].Add(sayePair_volume);
                    rulez[1].Add(sayePair_surface);
                }
            }
#if LOG_TIME
            stopWatch.Stop();
            long ts = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Calculated combo cutcell rule : {0}ms", ts);
#endif
        }

        class SayeFactoryWrapper : IQuadRuleFactory<QuadRule> {
            Action<ExecutionMask, int> comboRuleEvaluator;

            private IEnumerable<IChunkRulePair<QuadRule>> rule;

            Status ruleStatus;

            public SayeFactoryWrapper(
                Action<ExecutionMask, int> comboRuleEvaluator,
                IEnumerable<IChunkRulePair<QuadRule>> Rule,
                Status ruleStatus) {
                rule = Rule;
                this.comboRuleEvaluator = comboRuleEvaluator;
                this.ruleStatus = ruleStatus;
            }

            public RefElement RefElement => ruleStatus.ReferenceElement;

            //int counter = 0;

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
                if (!ruleStatus.Initialized || order != ruleStatus.Order) {
                    InitializeRule(order);
                    return UseExistingRule(mask);
                } else {
                    return UseExistingRule(mask);
                }
            }

            void InitializeRule(int Order) {
                comboRuleEvaluator(ruleStatus.MaxGrid, Order);
                ruleStatus.Order = Order;
                ruleStatus.Initialized = true;
            }

            IEnumerable<IChunkRulePair<QuadRule>> UseExistingRule(ExecutionMask subMask) {
                if (subMask.Equals(ruleStatus.MaxGrid)) {
                    return rule;
                }
                //If not, filter rules to fit subMask
                else {
                    //Debug.Assert(subMask.IsSubMaskOf(RuleStatus.initialMask));
                    if (!subMask.IsSubMaskOf(ruleStatus.MaxGrid)) {
                        throw new Exception("subMask is probably empty.");
                    }
                    List<IChunkRulePair<QuadRule>> subRulez = new List<IChunkRulePair<QuadRule>>(subMask.Count());

                    IEnumerator<int> initialMask_enum = ruleStatus.MaxGrid.GetItemEnumerator();
                    int i = 0;
                    BitArray subMask_bitmask = subMask.GetBitMask();
                    while (initialMask_enum.MoveNext()) {
                        if (subMask_bitmask[initialMask_enum.Current]) {
                            subRulez.Add(rule.ElementAt(i));
                        }
                        ++i;
                    }
                    return subRulez;
                }
            }

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

        }
    }
}
