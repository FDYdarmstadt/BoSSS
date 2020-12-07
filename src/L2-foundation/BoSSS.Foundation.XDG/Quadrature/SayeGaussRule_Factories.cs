#undef LOG_ACTIONS
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

    class LevelSetOnEdgeIntercepter : ISayeGaussComboRule
    {
        ISayeGaussComboRule comboRule;

        LevelSetOnEdgeRule specialRule;

        int intOrder;

        public LevelSetOnEdgeIntercepter(ISayeGaussComboRule rule, LevelSetTracker.LevelSetData data)
        {
            specialRule = new LevelSetOnEdgeRule(data);
            comboRule = rule;
        }

        public RefElement RefElement => comboRule.RefElement;

        public int order
        {
            set
            {
                intOrder = value;
                comboRule.order = value;
            }
        }

        public QuadRule[] ComboEvaluate(int cell)
        {
            if (specialRule.IsSpecialCell(cell))
            {
                var rule = specialRule.ComboQuadRule( intOrder, cell);
                return new QuadRule[]
                {
                    rule.surfaceRule.Rule,
                    rule.volumeRule.Rule
                };
            }
            else
            {
                return comboRule.ComboEvaluate(cell);
            }
        }

        public QuadRule Evaluate(int cell)
        {
            throw new NotImplementedException();
        }
    }

    class VolumeOnEdgeIntercepter : ISayeGaussRule
    {
        ISayeGaussRule sayeRule;

        LevelSetOnEdgeRule specialRule;

        int intOrder;

        public VolumeOnEdgeIntercepter(ISayeGaussRule rule, LevelSetTracker.LevelSetData data)
        {
            specialRule = new LevelSetOnEdgeRule(data);
            sayeRule = rule;
        }

        public RefElement RefElement => sayeRule.RefElement;
        
        public int order
        {
            set
            {
                intOrder = value;
                sayeRule.order = value;
            }
        }

        public QuadRule Evaluate(int cell)
        {
            if (specialRule.IsSpecialCell(cell))
            {
                var rule = specialRule.VolumeQuadRule(intOrder, cell);
                return rule.Rule;
            }
            else
            {
                return sayeRule.Evaluate(cell);
            }
        }
    }

    class SurfaceOnEdgeIntercepter : ISayeGaussRule
    {
        ISayeGaussRule sayeRule;

        LevelSetOnEdgeRule specialRule;

        int intOrder;

        public SurfaceOnEdgeIntercepter(ISayeGaussRule rule, LevelSetTracker.LevelSetData data)
        {
            specialRule = new LevelSetOnEdgeRule(data);
            sayeRule = rule;
        }

        public RefElement RefElement => sayeRule.RefElement;
        
        public int order
        {
            set
            {
                intOrder = value;
                sayeRule.order = value;
            }
        }

        public QuadRule Evaluate(int cell)
        {
            if (specialRule.IsSpecialCell(cell))
            {
                var rule = specialRule.SurfaceQuadRule(intOrder, cell);
                return rule.Rule;
            }
            else
            {
                return sayeRule.Evaluate(cell);
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


#if LOG_TIME
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
#endif
            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    QuadRule sayeRule = rule.Evaluate(cell);
                    ChunkRulePair<QuadRule> sayePair = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule);
                    result.Add(sayePair);
                }
            }
#if LOG_TIME
            stopWatch.Stop();
            long ts = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Calculated cutcell rule : {0}ms", ts);
#endif
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
            IRootFindingAlgorithm RootFinder)
        {
            ISayeGaussRule rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.PositiveVolume);
            ISayeGaussRule wrappedRule = new VolumeOnEdgeIntercepter(rule, _lsData);
            return new SayeGaussRuleFactory(wrappedRule);
        }

        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 3D case
        /// </summary>
        /// 
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_NegativeVolume3D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder) {
            ISayeGaussRule rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.NegativeVolume);
            ISayeGaussRule wrappedRule = new VolumeOnEdgeIntercepter(rule, _lsData);
            return new SayeGaussRuleFactory(wrappedRule);
        }

        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 3D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Surface3D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder)
        {
            ISayeGaussRule rule = new SayeFactory_Cube(
                _lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Surface
                );
            ISayeGaussRule wrappedRule = new SurfaceOnEdgeIntercepter(rule, _lsData);
            return new SayeGaussRuleFactory(wrappedRule);
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
            ISayeGaussRule wrappedRule = new VolumeOnEdgeIntercepter(rule, _lsData);
            return new SayeGaussRuleFactory(wrappedRule);
        }

        /// <summary>
        /// Gauss rules for \f$ \int_{\frakA \cap K_j } \ldots \dV \f$ in the 2D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_NegativeVolume2D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder
            ) {
            ISayeGaussRule rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.NegativeVolume
                );
            ISayeGaussRule wrappedRule = new VolumeOnEdgeIntercepter(rule, _lsData);
            return new SayeGaussRuleFactory(wrappedRule);
        }

        /// <summary>
        /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 2D case
        /// </summary>
        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Surface2D( 
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder)
        {
            ISayeGaussRule rule = new SayeFactory_Square(
                _lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.Surface);
            ISayeGaussRule wrappedRule = new SurfaceOnEdgeIntercepter(rule, _lsData);
            return new SayeGaussRuleFactory(wrappedRule);
        }


        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Surface(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder)
        {
            RefElement refElem = _lsData.GridDat.Grid.GetRefElement(0);
            if (refElem is Grid.RefElements.Square)
            {
                return SayeGaussRule_Surface2D(_lsData, RootFinder);
            }
            else if (refElem is Grid.RefElements.Cube)
            {
                return SayeGaussRule_Surface3D(_lsData, RootFinder);
            }
            else
            {
                throw new NotImplementedException("Saye quadrature not available for this RefElement");
            }
        }

        public static IQuadRuleFactory<QuadRule> SayeGaussRule_Volume(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder)
        {
            RefElement refElem = _lsData.GridDat.Grid.GetRefElement(0);
            if (refElem is Grid.RefElements.Square)
            {
                return SayeGaussRule_Volume2D(_lsData, RootFinder);
            }
            else if (refElem is Grid.RefElements.Cube)
            {
                return SayeGaussRule_Volume3D(_lsData, RootFinder);
            }
            else
            {
                throw new NotImplementedException("Saye quadrature not available for this RefElement");
            }
        }

        public static IQuadRuleFactory<QuadRule> SayeGaussRule_NegativeVolume(
           LevelSetTracker.LevelSetData _lsData,
           IRootFindingAlgorithm RootFinder) {
            RefElement refElem = _lsData.GridDat.Grid.GetRefElement(0);
            if (refElem is Grid.RefElements.Square) {
                return SayeGaussRule_NegativeVolume2D(_lsData, RootFinder);
            } else if (refElem is Grid.RefElements.Cube) {
                return SayeGaussRule_NegativeVolume3D(_lsData, RootFinder);
            } else {
                throw new NotImplementedException("Saye quadrature not available for this RefElement");
            }
        }

        #endregion

        #region Combo QuadRules

        public static SayeGaussComboRuleFactory SayeGaussRule_Combo2D(
            LevelSetTracker.LevelSetData lsData,
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussComboRule rule = new SayeFactory_Square(
                lsData,
                RootFinder,
                SayeFactory_Square.QuadratureMode.Combo
                );
            ISayeGaussComboRule wrappedRule = new LevelSetOnEdgeIntercepter(rule, lsData);

            CellMask maxGrid = lsData.GridDat.Cells.GetCells4Refelement(rule.RefElement).Intersect(
                lsData.Region.GetCutCellMask().ToGeometicalMask());
            return new SayeGaussComboRuleFactory(
                wrappedRule,
                maxGrid);
        }

        public static SayeGaussComboRuleFactory SayeGaussRule_Combo3D(
            LevelSetTracker.LevelSetData lsData,
            IRootFindingAlgorithm RootFinder
            )
        {
            ISayeGaussComboRule rule = new SayeFactory_Cube(
                lsData,
                RootFinder,
                SayeFactory_Cube.QuadratureMode.Combo
                );
            ISayeGaussComboRule wrappedRule = new LevelSetOnEdgeIntercepter(rule, lsData);


            CellMask maxGrid = lsData.GridDat.Cells.GetCells4Refelement(rule.RefElement).Intersect(
                lsData.Region.GetCutCellMask().ToGeometicalMask());
            return new SayeGaussComboRuleFactory(
                wrappedRule, 
                maxGrid);
        }

        public static SayeGaussComboRuleFactory SayeGaussRule_Combo(
           LevelSetTracker.LevelSetData _lsData,
           IRootFindingAlgorithm RootFinder)
        {
            RefElement refElem = _lsData.GridDat.Grid.GetRefElement(0);
            if (refElem is Grid.RefElements.Square)
            {
                return SayeGaussRule_Combo2D(_lsData, RootFinder);
            }
            else if (refElem is Grid.RefElements.Cube)
            {
                return SayeGaussRule_Combo3D(_lsData, RootFinder);
            }
            else
            {
                throw new NotImplementedException("Saye quadrature not available for this RefElement");
            }
        }

    #endregion
    }
}
