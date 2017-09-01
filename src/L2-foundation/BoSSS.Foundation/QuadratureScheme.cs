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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// A pairing of a quadrature rule factory and its assigned domain.
    /// </summary>
    /// <typeparam name="TQuadRule">
    /// The type of the quadrature rules created by the factory.
    /// </typeparam>
    /// <typeparam name="TDomain">
    /// The type of the domain to be considered (e.g. cells or edges).
    /// </typeparam>
    public interface IFactoryDomainPair<out TQuadRule, out TDomain>
        where TQuadRule : QuadRule
        where TDomain : ExecutionMask {

        /// <summary>
        /// A factory for quadrature rules
        /// </summary>
        IQuadRuleFactory<TQuadRule> RuleFactory {
            get;
        }

        /// <summary>
        /// The domain of integration
        /// </summary>
        TDomain Domain {
            get;
        }

        /// <summary>
        /// An optional order that overrides the order passed to
        /// <see cref="RuleFactory"/>.
        /// </summary>
        int? Order {
            get;
        }
    }

    /// <summary>
    /// Defines a mapping between quadrature rule factories (see
    /// <see cref="IQuadRuleFactory{T}"/>) and the domains in which they should
    /// be used. The main purpose of implementations is to provide a measure to
    /// represent different domains and their different requirements on
    /// quadrature rules by a single object that can optimize the use of
    /// resources (see <see cref="Compile"/>)
    /// </summary>
    /// <typeparam name="TQuadRule">
    /// The type of the individual quadrature rules
    /// </typeparam>
    /// <typeparam name="TDomain">
    /// The type of the individual domains.
    /// </typeparam>
    public interface IQuadratureScheme<out TQuadRule, out TDomain>
        where TQuadRule : QuadRule
        where TDomain : ExecutionMask {

        /// Background integration domain. This domain is dominant over all
        /// optional domains specified in the factory chain (see
        /// <see cref="FactoryChain"/>) in the sense that all operations will
        /// be restricted to this domain.
        TDomain Domain {
            get;
        }

        /// <summary>
        /// Quadrature rule factories which should be applied. The higher
        /// quadrature rules (those at higher indices) dominate, i.e. they
        /// overwrite the lower ones.
        /// </summary>
        IEnumerable<IFactoryDomainPair<TQuadRule, TDomain>> FactoryChain {
            get;
        }

        /// <summary>
        /// Compiles the chain of factories (see <see cref="FactoryChain"/>)
        /// into a composite quadrature rules, i.e. a fixed mapping of chunks
        /// and quadrature rules.
        /// </summary>
        /// <param name="gridData">
        /// Information about the grid.
        /// </param>
        /// <param name="order">
        /// The minimum desired order of accuracy (if applicable)
        /// </param>
        /// <returns>
        /// A quadrature
        /// </returns>
        ICompositeQuadRule<TQuadRule> Compile(IGridData gridData, int order);
    }

    /// <summary>
    /// Base class for all quadrature schemes. Inherit from this class in order
    /// to fix a specific combination of the type of quadrature rule and the
    /// type of the occurring domains.
    /// </summary>
    /// <typeparam name="TQuadRule">
    /// <see cref="IQuadratureScheme{S, T}"/>
    /// </typeparam>
    /// <typeparam name="TDomain">
    /// <see cref="IQuadratureScheme{S, T}"/>
    /// </typeparam>
    /// <seealso cref="IQuadratureScheme{S, T}"/>
    public abstract class QuadratureScheme<TQuadRule, TDomain> : IQuadratureScheme<TQuadRule, TDomain>
        where TQuadRule : QuadRule
        where TDomain : ExecutionMask {

        /// <summary>
        /// <see cref="IQuadratureScheme{S, T}.FactoryChain"/>
        /// </summary>
        private readonly List<IFactoryDomainPair<TQuadRule, TDomain>> factoryChain
            = new List<IFactoryDomainPair<TQuadRule, TDomain>>();


        /// <summary>
        /// Creates a scheme with a given background domain.
        /// </summary>
        /// <param name="domain">
        /// Optional: A background domain (see
        /// <see cref="IQuadratureScheme{S, T}.Domain"/>). If null, the domain
        /// created by <see cref="GetDefaultDomain"/> will be used.
        /// </param>
        /// <param name="UseDefaultFactories">
        /// if true, quadrature rule factories for default (Gaussian) rules will be added for the domain <paramref name="domain"/>;
        /// if false, the user must add factories for all items in the domain.
        /// </param>
        public QuadratureScheme(bool UseDefaultFactories, TDomain domain = null) {
            this.Domain = domain;
            this.m_UseDefaultFactories = UseDefaultFactories;
        }

        bool m_UseDefaultFactories;

        #region IQuadratureInstruction<TQuadRule,TDomain> Members

        /// <summary>
        /// The background domain, see
        /// <see cref="IQuadratureScheme{S, T}.Domain"/>
        /// </summary>
        public TDomain Domain {
            get;
            private set;
        }

        /// <summary>
        /// <see cref="IQuadratureScheme{S, T}.FactoryChain"/>
        /// </summary>
        public IEnumerable<IFactoryDomainPair<TQuadRule, TDomain>> FactoryChain {
            get {
                return factoryChain;
            }
        }

        /// <summary>
        /// Allows to add factories to <see cref="FactoryChain"/>. Factories
        /// added <b>later</b> will <b>dominate</b> over the ones that have
        /// been added earlier. That is, in case of conflicts (more than one
        /// factory assigned to the same (sub-)domain), the one added later
        /// will be chosen.
        /// </summary>
        /// <param name="factory">
        /// The factory to be added.
        /// </param>
        /// <param name="domain">
        /// The domain on which <paramref name="factory"/> should act.
        /// </param>
        /// <param name="order">
        /// An optional order that overrides the order passed to
        /// <see cref="Compile"/>.
        /// </param>
        /// <returns>
        /// This object (to allow for fluent addition of multiple factories).
        /// </returns>
        public QuadratureScheme<TQuadRule, TDomain> AddFactoryDomainPair(
            IQuadRuleFactory<TQuadRule> factory, TDomain domain = null, int? order = null) {
            //
            if (factory == null)
                throw new ArgumentNullException();
            factoryChain.Add(new FactoryDomainPair(factory, domain, order));
            return this;
        }

        /// <summary>
        /// Constructs the <see cref="ExecutionMask"/> of the appropriate type
        /// containing all cells/edges with reference element <paramref name="E"/>
        /// </summary>
        /// <param name="E"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        TDomain GetDomainForRefElement(Grid.RefElements.RefElement E, IGridData g) {
            ExecutionMask em = null;

            if (typeof(TDomain) == typeof(EdgeMask) || typeof(TDomain).IsSubclassOf(typeof(EdgeMask)))
                em = g.iLogicalEdges.GetEdges4RefElement(E);
            else if (typeof(TDomain) == typeof(CellMask) || typeof(TDomain).IsSubclassOf(typeof(CellMask)))
                em = g.iLogicalCells.GetCells4Refelement(E);
            else
                throw new NotImplementedException();


            return ((TDomain)(em));
        }

        /// <summary>
        /// <see cref="IQuadratureScheme{S, T}.Compile"/>
        /// </summary>
        public ICompositeQuadRule<TQuadRule> Compile(IGridData gridData, int order) {

            // set domain
            TDomain baseDomain = Domain ?? GetDefaultDomain(gridData);

            // identify the reference elements
            // ===============================

            RefElement[] RefElements;
            //ID
            RefElements = GetRefElements(gridData);

            // define standard factories, if necessary
            // =======================================
            IEnumerable<IFactoryDomainPair<TQuadRule, TDomain>> factoryDomainPairs;
            if (m_UseDefaultFactories) {
                var defaultFactories = RefElements.Select(
                        RefElm => (new FactoryDomainPair(
                            GetDefaultRuleFactory(gridData, RefElm),
                            GetDomainForRefElement(RefElm, gridData).Intersect(baseDomain)))
                            );

                var _factoryDomainPairs = new List<IFactoryDomainPair<TQuadRule, TDomain>>();
                _factoryDomainPairs.AddRange(defaultFactories);
                _factoryDomainPairs.AddRange(FactoryChain);
                factoryDomainPairs = _factoryDomainPairs;
            } else {
                factoryDomainPairs = FactoryChain;
            }

            // apply factories
            // ---------------
            CompositeQuadRule<TQuadRule> fullRule = new CompositeQuadRule<TQuadRule>();
            int i = -1;
            foreach (var factoryDomainPair in factoryDomainPairs) {
                i++;
                TDomain currentDomain = baseDomain;
                var RefElm = factoryDomainPair.RuleFactory.RefElement;

                if (factoryDomainPair.Domain == null) {
                    currentDomain = GetDomainForRefElement(RefElm, gridData).Intersect(baseDomain);
                } else {
                    currentDomain = baseDomain.Intersect(factoryDomainPair.Domain);
                    Debug.Assert(currentDomain.Except(GetDomainForRefElement(RefElm, gridData)).NoOfItemsLocally <= 0);
                }

                // check the type of factory
                if (!RefElements.Contains(factoryDomainPair.RuleFactory.RefElement)) {
                    string simplexNames = "";
                    for (int _i = 0; _i < RefElements.Length; _i++) {
                        simplexNames += RefElements[_i].GetType().Name;
                        if (_i < RefElements.Length - 1)
                            simplexNames += ",";
                    }

                    throw new NotSupportedException("Found an illegal factory in the chain during quadrature scheme compilation: factories are required to work on '"
                        + simplexNames
                        + "'-elements, but this one works on '"
                        + factoryDomainPair.RuleFactory.RefElement.GetType().ToString() + "'-elements.");
                }

                //// Causes parallel issues for classic HMF
                //// -> deactivated by Björn until classic HMF is fixed
                //if(currentDomain != null && currentDomain.NoOfItemsLocally <= 0) {
                //    continue;
                //}

                CompositeQuadRule<TQuadRule> currentRule = CompositeQuadRule<TQuadRule>.Create(
                    factoryDomainPair.RuleFactory, factoryDomainPair.Order ?? order, currentDomain);

                int prevIE = -1;
                foreach (var CRP in currentRule) {
                    if (prevIE >= 0) {
                        if (CRP.Chunk.i0 < prevIE) {
                            throw new ApplicationException("Quadrature rule factory " + factoryDomainPair.RuleFactory.GetType().Name + " produced an un-sorted quadrature rule!");
                        }
                    }
                    prevIE = CRP.Chunk.JE;
                }

#if DEBUG
                //// Check removed since there are situations where it is valid _not_
                //// to return quadrature rule for a given cell/edge
                Debug.Assert(currentRule.NumberOfItems == currentDomain.NoOfItemsLocally);

                Debug.Assert(
                     currentRule.Where(w => (w.Rule.Nodes.IsLocked == false)).Count() <= 0,
                     "error in quadrature rule creation: factory delivered some rule non-locked node set.");


                Debug.Assert(
                    currentRule.Where(w => (w.Rule.NoOfNodes <= 0)).Count() <= 0,
                    "error in quadrature rule creation: factory delivered some rule with zero nodes.");

                foreach(IChunkRulePair<TQuadRule> cr in currentRule ) {
                    cr.Rule.Weights.CheckForNanOrInf();
                    cr.Rule.Nodes.CheckForNanOrInf();
                }

#endif

                if (i == 0) {
                    fullRule = currentRule;
                } else {
                    fullRule = CompositeQuadRule<TQuadRule>.Merge(fullRule, currentRule);
                }

                if (fullRule != null && fullRule.Count() > 0) {
                    int J = fullRule.Max(crp => crp.Chunk.JE);
                    System.Collections.BitArray ChunkTest = new System.Collections.BitArray(J);
                    foreach (var chunk in fullRule) {
                        int IE = chunk.Chunk.JE;
                        for (int ii = chunk.Chunk.i0; ii < IE; ii++) {
                            if (ChunkTest[ii])
                                throw new ArgumentException("More than one quadrature rule defined for integration item " + ii + ".");
                            ChunkTest[ii] = true;

                        }
                    }
                }

               
            }

            //// Check removed since there are situations where it is valid _not_
            //// to return quadrature rule for a given cell/edge
            //if (fullRule.NumberOfItems < baseDomain.NoOfItemsLocally) {
            //    throw new InvalidOperationException(
            //        "Insufficient quadrature rule factories for desired domain:"
            //        + " (domain has " + baseDomain.NoOfItemsLocally
            //        + " items, but quadrature rules were only found for "
            //        + fullRule.NumberOfItems + " items.");
            //}

            return fullRule;
        }

        private static Grid.RefElements.RefElement[] GetRefElements(IGridData gridData) {
            Grid.RefElements.RefElement[] RefElements;
            if (typeof(TDomain) == typeof(EdgeMask) || typeof(TDomain).IsSubclassOf(typeof(EdgeMask)))
                RefElements = gridData.iGeomEdges.EdgeRefElements;
            else if (typeof(TDomain) == typeof(CellMask) || typeof(TDomain).IsSubclassOf(typeof(CellMask)))
                RefElements = gridData.iGeomCells.RefElements;
            else
                throw new NotImplementedException();
            return RefElements;
        }

#endregion

        /// <summary>
        /// Implement this method by returning a default domain, e.g. all cells
        /// or edges of the grid.
        /// </summary>
        /// <param name="gridData">
        /// Information about the grid.
        /// </param>
        /// <returns>
        /// A default domain of integration.
        /// </returns>
        protected abstract TDomain GetDefaultDomain(IGridData gridData);

        /// <summary>
        /// Implement this method by returning a default factory that should
        /// be valid in the standard case (e.g., that returns Gauss rules for
        /// the reference element <paramref name="elem"/>)
        /// </summary>
        /// <param name="gridData">
        /// Information about the grid.
        /// </param>
        /// <param name="elem">
        /// reference element
        /// </param>
        /// <returns>
        /// A default quadrature rules factory.
        /// </returns>
        protected abstract IQuadRuleFactory<TQuadRule> GetDefaultRuleFactory(IGridData gridData, RefElement elem);

        /// <summary>
        /// Plain implementation of <see cref="IFactoryDomainPair{S, T}"/>
        /// </summary>
        private class FactoryDomainPair : IFactoryDomainPair<TQuadRule, TDomain> {

            /// <summary>
            /// Just stores the given values.
            /// </summary>
            /// <param name="ruleFactory">
            /// <see cref="RuleFactory"/>
            /// </param>
            /// <param name="domain">
            /// <see cref="Domain"/>
            /// </param>
            /// <param name="order">
            /// The desired order of the quadrature rule. Note that this is
            /// order not necessarily achieved by all
            /// <paramref name="ruleFactory"/>s
            /// </param>
            public FactoryDomainPair(IQuadRuleFactory<TQuadRule> ruleFactory, TDomain domain, int? order = null) {
                this.RuleFactory = ruleFactory;
                this.Domain = domain;
                this.Order = order;
            }

#region IFactoryDomainPair<TQuadRule,TDomain> Members

            /// <summary>
            /// <see cref="IFactoryDomainPair{S, T}.RuleFactory"/>
            /// </summary>
            public IQuadRuleFactory<TQuadRule> RuleFactory {
                get;
                private set;
            }

            /// <summary>
            /// <see cref="IFactoryDomainPair{S, T}.Domain"/>
            /// </summary>
            public TDomain Domain {
                get;
                private set;
            }

            /// <summary>
            /// An optional order that overrides the order passed to
            /// <see cref="RuleFactory"/>.
            /// </summary>
            public int? Order {
                get;
                private set;
            }

#endregion
        }
    }

    /// <summary>
    /// Extension methods for <see cref="IQuadratureScheme{A,B}"/> and various
    /// classes which implement this interface.
    /// </summary>
    public static class QuadratureScheme_Ext {

        /// <summary>
        /// adds a factory to a quadrature scheme
        /// </summary>
        public static R AddFactory<R, TQuadRule, TDomain>(this R scheme, IQuadRuleFactory<TQuadRule> factory, TDomain domain = null)
            where R : QuadratureScheme<TQuadRule, TDomain>
            where TQuadRule : QuadRule
            where TDomain : ExecutionMask {

            scheme.AddFactoryDomainPair(factory, domain);
            return scheme;
        }

        /// <summary>
        /// adds a factory to a cell-boundary quadrature scheme
        /// </summary>
        public static CellBoundaryQuadratureScheme AddFactory<TQuadRule>(this CellBoundaryQuadratureScheme scheme, IQuadRuleFactory<TQuadRule> factory)
            where TQuadRule : CellBoundaryQuadRule {
            scheme.AddFactoryDomainPair(factory, default(CellMask));
            return scheme;
        }

        /// <summary>
        /// adds a factory to a cell quadrature scheme
        /// </summary>
        public static CellQuadratureScheme AddFactory<TQuadRule>(this CellQuadratureScheme scheme, IQuadRuleFactory<TQuadRule> factory)
            where TQuadRule : QuadRule {
            scheme.AddFactoryDomainPair(factory, default(CellMask));
            return scheme;
        }

        /// <summary>
        /// adds a factory to an edge quadrature scheme
        /// </summary>
        public static EdgeQuadratureScheme AddFactory<TQuadRule>(this EdgeQuadratureScheme scheme, IQuadRuleFactory<TQuadRule> factory)
            where TQuadRule : QuadRule {
            scheme.AddFactoryDomainPair(factory, default(EdgeMask));
            return scheme;
        }

        /// <summary>
        /// adds a fixed rule to the quadrature scheme <paramref name="scheme"/>;
        /// this can be used to explicitly specify a quadrature rule, the choice of the order during the compilation of the scheme
        /// (see <see cref="QuadratureScheme{A,B}.Compile"/>) will have no effect.
        /// </summary>
        public static R AddFixedRule<R, TQuadRule, TDomain>(this R scheme, TQuadRule fixedRule, TDomain domain = null)
            where R : QuadratureScheme<TQuadRule, TDomain>
            where TQuadRule : QuadRule
            where TDomain : ExecutionMask {

            scheme.AddFactoryDomainPair(new FixedRuleFactory<TQuadRule>(fixedRule), domain);
            return scheme;
        }

        /// <summary>
        /// adds multiple fixed rules to the quadrature scheme <paramref name="scheme"/>;
        /// this can be used to explicitly specify a quadrature rule, the choice of the order during the compilation of the scheme
        /// (see <see cref="QuadratureScheme{A,B}.Compile"/>) will have no effect.
        /// </summary>
        public static R AddFixedRuleS<R, TQuadRule, TDomain>(this R scheme, IEnumerable<TQuadRule> fixedRules, TDomain[] domain = null)
            where R : QuadratureScheme<TQuadRule, TDomain>
            where TQuadRule : QuadRule
            where TDomain : ExecutionMask {

            if (domain != null)
                if (domain.Count() != fixedRules.Count())
                    throw new ArgumentException();
            for (int i = 0; i < fixedRules.Count(); i++) {
                scheme.AddFactoryDomainPair(new FixedRuleFactory<TQuadRule>(fixedRules.ElementAt(i)), (domain != null) ? domain[i] : null);
            }
            return scheme;
        }

        /// <summary>
        /// adds multiple fixed rules to the quadrature scheme <paramref name="scheme"/>;
        /// this can be used to explicitly specify a quadrature rule, the choice of the order during the compilation of the scheme
        /// (see <see cref="QuadratureScheme{A,B}.Compile"/>) will have no effect.
        /// </summary>
        public static CellQuadratureScheme AddFixedRuleS<TQuadRule>(this CellQuadratureScheme scheme, IEnumerable<TQuadRule> fixedRules, CellMask[] domain = null)
            where TQuadRule : QuadRule {

            if (domain != null)
                if (domain.Count() != fixedRules.Count())
                    throw new ArgumentException();
            for (int i = 0; i < fixedRules.Count(); i++) {
                scheme.AddFactoryDomainPair(new FixedRuleFactory<TQuadRule>(fixedRules.ElementAt(i)), (domain != null) ? domain[i] : null);
            }
            return scheme;
        }

        /// <summary>
        /// adds rules with fixed quadrature order to the quadrature scheme <paramref name="scheme"/>;
        /// this can be used to explicitly specify the quadrature order, the choice of the order during the compilation of the scheme
        /// (see <see cref="QuadratureScheme{A,B}.Compile"/>) will have no effect.
        /// </summary>
        public static CellQuadratureScheme AddFixedOrderRules(this CellQuadratureScheme scheme, IGridData GridDat, int order) {
            var QRs = GridDat.iGeomCells.RefElements.Select(Kref => Kref.GetQuadratureRule(order));
            return scheme.AddFixedRuleS<QuadRule>(QRs);
        }

        /// <summary>
        /// adds rules with fixed quadrature order to the quadrature scheme <paramref name="scheme"/>;
        /// this can be used to explicitly specify the quadrature order, the choice of the order during the compilation of the scheme
        /// (see <see cref="QuadratureScheme{A,B}.Compile"/>) will have no effect.
        /// </summary>
        public static EdgeQuadratureScheme AddFixedOrderRules(this EdgeQuadratureScheme scheme, IGridData GridDat, int order) {
            var QRs = GridDat.iGeomEdges.EdgeRefElements.Select(Kref => Kref.GetQuadratureRule(order));
            return scheme.AddFixedRuleS<EdgeQuadratureScheme, QuadRule, EdgeMask>(QRs);
        }

        /// <summary>
        /// adds a standard (Gaussian) quadrature rule to the quadrature scheme <paramref name="scheme"/>;
        /// the actual order of the quadrature rule will be determined by the order which is passed 
        /// at the scheme compilation (see <see cref="QuadratureScheme{A,B}.Compile"/>);
        /// </summary>
        public static R AddStandardRule<R, TDomain>(this R scheme, Grid.RefElements.RefElement s, TDomain domain)
            where R : QuadratureScheme<QuadRule, TDomain>
            where TDomain : ExecutionMask {

            scheme.AddFactoryDomainPair(new StandardQuadRuleFactory(s), domain);
            //    throw new NotImplementedException();
            return scheme;
        }

        /// <summary>
        /// adds a standard(Gaussian) quadrature rule to the quadrature scheme <paramref name="scheme"/>.
        /// </summary>
        public static CellQuadratureScheme AddStandardRule(this CellQuadratureScheme scheme, Grid.RefElements.RefElement s) {
            scheme.AddFactoryDomainPair(new StandardQuadRuleFactory(s), null);
            return scheme;
        }

        /// <summary>
        /// adds a standard(Gaussian) quadrature rule to the quadrature scheme <paramref name="scheme"/>.
        /// </summary>
        public static EdgeQuadratureScheme AddStandardRule(this EdgeQuadratureScheme scheme, Grid.RefElements.RefElement s) {
            scheme.AddFactoryDomainPair(new StandardQuadRuleFactory(s), null);
            return scheme;
        }

        /// <summary>
        /// Compilation of a cell volume scheme, even if it is null!
        /// </summary>
        public static ICompositeQuadRule<QuadRule> SaveCompile(this CellQuadratureScheme scheme, IGridData g, int order) {
            // sometimes, a bit of repetition seems easier ...
            scheme = (scheme ?? new CellQuadratureScheme(true));
            return scheme.Compile(g, order);
        }

        /// <summary>
        /// Compilation of an edge scheme, even if it is null!
        /// </summary>
        public static ICompositeQuadRule<QuadRule> SaveCompile(this EdgeQuadratureScheme scheme, IGridData g, int order) {
            // sometimes, a bit of repetition seems easier ...
            scheme = (scheme ?? new EdgeQuadratureScheme(true));
            return scheme.Compile(g, order);
        }
    }

    /// <summary>
    /// Quadrature scheme for standard volume integration over cells.
    /// </summary>
    public sealed class CellQuadratureScheme : QuadratureScheme<QuadRule, CellMask> {

        /// <summary>
        /// Constructs an empty quadrature scheme.
        /// </summary>
        /// <param name="domain">
        /// Background domain.
        /// </param>
        /// <param name="UseDefaultFactories">
        /// if true, quadrature rule factories for default (Gaussian) rules will be added for the domain <paramref name="domain"/>;
        /// if false, the user must add factories for all items in the domain.
        /// </param>
        public CellQuadratureScheme(bool UseDefaultFactories = true, CellMask domain = null)
            : base(UseDefaultFactories, domain) {
        }

        /// <summary>
        /// Convenience constructor that allows for the construction of a
        /// scheme with a predefined factory. Equivalent to creating an empty
        /// scheme and calling
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>(<paramref name="factory"/>, <paramref name="domain"/>).
        /// </summary>
        /// <param name="factory">
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>
        /// </param>
        /// <param name="domain">
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>
        /// </param>
        public CellQuadratureScheme(IQuadRuleFactory<QuadRule> factory, CellMask domain = null)
            : base(false, domain) {
            AddFactoryDomainPair(factory, domain);
        }

        /// <summary>
        /// <see cref="CellMask.GetFullMask"/>
        /// </summary>
        /// <param name="gridData">
        /// <see cref="QuadratureScheme{S, T}.GetDefaultDomain"/>
        /// </param>
        /// <returns>
        /// A mask containing all cells of the grid.
        /// </returns>
        protected override CellMask GetDefaultDomain(IGridData gridData) {
            var ret = CellMask.GetFullMask(gridData);
            return ret;
        }

        /// <summary>
        /// <see cref="StandardQuadRuleFactory"/>
        /// </summary>
        /// <param name="gridData">
        /// <see cref="QuadratureScheme{S, T}.GetDefaultRuleFactory"/>
        /// </param>
        /// <param name="elem">
        /// reference element, must be one of <see cref="BoSSS.Foundation.Grid.GridCommons.RefElements"/>
        /// </param>
        /// <returns>
        /// An instance of <see cref="StandardQuadRuleFactory"/> for the volume
        /// simplex.
        /// </returns>
        protected override IQuadRuleFactory<QuadRule> GetDefaultRuleFactory(IGridData gridData, RefElement elem) {
            return new StandardQuadRuleFactory(elem);
        }

        //new public  CellQuadratureScheme AddFactory(IQuadRuleFactory<QuadRule> f, CellMask dom = null) {
        //    base.AddFactory(f, dom);
        //    return this;
        //}


    }

    /// <summary>
    /// Quadrature scheme for standard integration over edges.
    /// </summary>
    public sealed class EdgeQuadratureScheme : QuadratureScheme<QuadRule, EdgeMask> {

        /// <summary>
        /// Constructs an empty quadrature scheme.
        /// </summary>
        /// <param name="domain">
        /// Background domain.
        /// </param>
        /// <param name="UseDefaultFactories">
        /// if true, quadrature rule factories for default (Gaussian) rules will be added for the domain <paramref name="domain"/>;
        /// if false, the user must add factories for all items in the domain.
        /// </param>
        public EdgeQuadratureScheme(bool UseDefaultFactories = true, EdgeMask domain = null)
            : base(UseDefaultFactories, domain) {
        }

        /// <summary>
        /// Convenience constructor that allows for the construction of a
        /// scheme with a predefined factory. Equivalent to creating an empty
        /// scheme and calling
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>(<paramref name="factory"/>, <paramref name="domain"/>).
        /// </summary>
        /// <param name="factory">
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>
        /// </param>
        /// <param name="domain">
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>
        /// </param>
        public EdgeQuadratureScheme(IQuadRuleFactory<QuadRule> factory, EdgeMask domain = null)
            : base(false, domain) {
            AddFactoryDomainPair(factory, domain);
        }

        /// <summary>
        /// <see cref="EdgeMask.GetFullMask"/>
        /// </summary>
        /// <param name="gridData">
        /// <see cref="QuadratureScheme{S, T}.GetDefaultDomain"/>
        /// </param>
        /// <returns>
        /// A mask containing all edges of the grid.
        /// </returns>
        protected override EdgeMask GetDefaultDomain(IGridData gridData) {
            return EdgeMask.GetFullMask(gridData);
        }

        /// <summary>
        /// <see cref="StandardQuadRuleFactory"/>
        /// </summary>
        /// <param name="gridData">
        /// <see cref="QuadratureScheme{S, T}.GetDefaultRuleFactory"/>
        /// </param>
        /// <returns>
        /// An instance of <see cref="StandardQuadRuleFactory"/> for the edge
        /// simplex.
        /// </returns>
        /// <param name="Kref">
        /// an edge reference element, must be one of <see cref="BoSSS.Foundation.Grid.GridData.EdgeData.EdgeRefElements"/>.
        /// </param>
        protected override IQuadRuleFactory<QuadRule> GetDefaultRuleFactory(IGridData gridData, RefElement Kref) {
            return new StandardQuadRuleFactory(Kref);
        }
    }

    /// <summary>
    /// Quadrature scheme for standard integration over the boundaries of
    /// cells. That is, for a given <b>cell</b>, integration is performed over
    /// all edges adjacent to this cell.
    /// </summary>
    public sealed class CellBoundaryQuadratureScheme : QuadratureScheme<CellBoundaryQuadRule, CellMask> {

        /// <summary>
        /// Constructs an empty quadrature scheme.
        /// </summary>
        /// <param name="domain">
        /// Background domain.
        /// </param>
        /// <param name="UseDefaultFactories">
        /// if true, quadrature rule factories for default (Gaussian) rules will be added for the domain <paramref name="domain"/>;
        /// if false, the user must add factories for all items in the domain.
        /// </param>
        public CellBoundaryQuadratureScheme(bool UseDefaultFactories, CellMask domain = null)
            : base(UseDefaultFactories, domain) {
        }

        /// <summary>
        /// Convenience constructor that allows for the construction of a
        /// scheme with a predefined factory. Equivalent to creating an empty
        /// scheme and calling
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>(<paramref name="factory"/>, <paramref name="domain"/>).
        /// </summary>
        /// <param name="factory">
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>
        /// </param>
        /// <param name="domain">
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>
        /// </param>
        public CellBoundaryQuadratureScheme(IQuadRuleFactory<CellBoundaryQuadRule> factory, CellMask domain = null)
            : base(false, domain) {
            AddFactoryDomainPair(factory, domain);
        }

        /// <summary>
        /// <see cref="CellMask.GetFullMask"/>
        /// </summary>
        /// <param name="gridData">
        /// <see cref="QuadratureScheme{S, T}.GetDefaultDomain"/>
        /// </param>
        /// <returns>
        /// A mask containing all cells of the grid.
        /// </returns>
        protected override CellMask GetDefaultDomain(IGridData gridData) {
            return CellMask.GetFullMask(gridData);
        }

        /// <summary>
        /// <see cref="CellBoundaryFromEdgeRuleFactory{T}"/>
        /// </summary>
        /// <param name="gridData">
        /// <see cref="QuadratureScheme{S, T}.GetDefaultRuleFactory"/>
        /// </param>
        /// <param name="elem">
        /// reference element, must be one of <see cref="BoSSS.Foundation.Grid.GridCommons.RefElements"/>
        /// </param>
        /// <returns>
        /// An instance of <see cref="CellBoundaryFromEdgeRuleFactory{T}"/>
        /// based on a <see cref="StandardQuadRuleFactory"/> for the edge
        /// simplex.
        /// </returns>
        protected override IQuadRuleFactory<CellBoundaryQuadRule> GetDefaultRuleFactory(IGridData gridData, RefElement elem) {
            return new StandardCellBoundaryQuadRuleFactory(elem);
        }
    }


    /// <summary>
    /// quadrature scheme for the <see cref="DoubleEdgeQuadrature"/>-class.
    /// </summary>
    public sealed class DoubleEdgeQuadratureScheme : QuadratureScheme<DoubleEdgeQuadRule, EdgeMask> {


        /// <summary>
        /// Constructs an empty quadrature scheme.
        /// </summary>
        /// <param name="domain">
        /// Background domain.
        /// </param>
        /// <param name="UseDefaultFactories">
        /// if true, quadrature rule factories for default (Gaussian) rules will be added for the domain <paramref name="domain"/>;
        /// if false, the user must add factories for all items in the domain.
        /// </param>
        public DoubleEdgeQuadratureScheme(bool UseDefaultFactories, EdgeMask domain = null)
            : base(UseDefaultFactories, domain) {
        }

        /// <summary>
        /// Convenience constructor that allows for the construction of a
        /// scheme with a predefined factory. Equivalent to creating an empty
        /// scheme and calling
        /// <see cref="QuadratureScheme{S, T}.AddFactoryDomainPair"/>(<paramref name="factory"/>, <paramref name="domain"/>).
        /// </summary>
        public DoubleEdgeQuadratureScheme(IQuadRuleFactory<DoubleEdgeQuadRule> factory, EdgeMask domain = null)
            : base(false, domain) {
            AddFactoryDomainPair(factory, domain);
        }

        /// <summary>
        /// all edges of the grid.
        /// </summary>
        protected override EdgeMask GetDefaultDomain(IGridData gridData) {
            return EdgeMask.GetFullMask(gridData);
        }

        /// <summary>
        /// returns a <see cref="StandardDoubleEdgeRuleFactory"/>-object.
        /// </summary>
        protected override IQuadRuleFactory<DoubleEdgeQuadRule> GetDefaultRuleFactory(IGridData gridData, RefElement elem) {
            return new StandardDoubleEdgeRuleFactory(gridData, elem);
        }
    }
}
