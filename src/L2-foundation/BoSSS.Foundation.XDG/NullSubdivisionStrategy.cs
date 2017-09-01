/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System.Collections.Generic;
using BoSSS.Foundation.Grid;

namespace BoSSS.Foundation.XDG.Quadrature.Subdivision {

    /// <summary>
    /// A subdivisions strategy that does not subdivide anything. Mainly exists
    /// in order to avoid error-prone null-checks (Null object pattern).
    /// </summary>
    public class NullSubdivisionStrategy : ISubdivisionStrategy {

        /// <summary>
        /// The simplex (not) to be subdivided.
        /// </summary>
        private readonly Simplex simplex;

        /// <summary>
        /// Constructs the thing.
        /// </summary>
        /// <param name="simplex">
        /// The simplex (not) to be subdivided.
        /// </param>
        public NullSubdivisionStrategy(Simplex simplex) {
            this.simplex = simplex;
        }

        #region ISubdivisionStrategy Members

        /// <summary>
        /// <see cref="ISubdivisionStrategy.Simplex"/>
        /// </summary>
        public Simplex Simplex {
            get;
            private set;
        }

        /// <summary>
        /// Uses <see cref="SubdivisionNode.NullSubdivisionNode"/> to create an
        /// identity transformation for the given simplex. The chunks of
        /// <paramref name="mask"/> are not altered.
        /// </summary>
        /// <param name="mask">
        /// <see cref="ISubdivisionStrategy.GetSubdivisionNodes"/>
        /// </param>
        /// <returns>
        /// <see cref="ISubdivisionStrategy.GetSubdivisionNodes"/>
        /// </returns>
        public IEnumerable<KeyValuePair<Chunk, IEnumerable<SubdivisionNode>>> GetSubdivisionNodes(ExecutionMask mask) {
            foreach (Chunk chunk in mask) {
                yield return new KeyValuePair<Chunk, IEnumerable<SubdivisionNode>>(
                    chunk,
                    new SubdivisionNode[] { SubdivisionNode.NullSubdivisionNode(simplex.SpatialDimension) });
            }
        }

        #endregion
    }
}
