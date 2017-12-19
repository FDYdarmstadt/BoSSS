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
using System.Runtime.InteropServices;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// Represents a (one-dimensional) line segment in two- or
    /// three-dimensional space. The line segment generally lives in the
    /// reference coordinate system and is parametrized via a parameter t
    /// ranging from -1 to 1. Here, t=-1 marks the start point and t=1 marks
    /// the end point of the segment.
    /// </summary>
    public class LineSegment : IEquatable<LineSegment>, IObserver<LevelSetTracker.LevelSetRegions> {

        /// <summary>
        /// Minimal distance between two points.
        /// </summary>
        private const double EPSILON = 1e-14;

        /// <summary>
        /// Return an instance of <see cref="SafeGuardedNewtonMethod"/> where
        /// the tolerance is set to 1e-10
        /// </summary>
        public static readonly IRootFindingAlgorithm DefaultRootFindingAlgorithm
            = new SafeGuardedNewtonMethod(1e-14);

        /// <summary>
        /// The Cartesian length of this segment
        /// </summary>
        private double? length = null;

        /// <summary>
        /// Cache for the roots found on this line segment. The keys are given
        /// by a level set - cell pair.
        /// </summary>
        private Dictionary<Tuple<ILevelSet, int>, double[]> rootsCache = new Dictionary<Tuple<ILevelSet, int>, double[]>();

        RefElement m_Kref;

        /// <summary>
        /// Dingsbums
        /// </summary>
        public RefElement RefElement {
            get {
                return m_Kref;
            }
        }

        /// <summary>
        /// Constructs a new segment with the given start and end points
        /// </summary>
        /// <param name="spatialDimension">
        /// The spatial dimension of start and end points
        /// </param>
        /// <param name="start">
        /// The start point of the segment
        /// </param>
        /// <param name="end">
        /// The end point of the segment
        /// </param>
        /// <param name="iVertexStart">
        /// Optional vertex index of the start point in some list of vertices
        /// containing <paramref name="start"/>. Used by some callers to evade
        /// equality comparisons of double values.
        /// </param>
        /// <param name="iVertexEnd">
        /// Optional vertex index of the end point in some list of vertices
        /// containing <paramref name="end"/>. Used by some callers to evade
        /// equality comparisons of double values.
        /// </param>
        /// <param name="Kräf">
        /// Reference element.
        /// </param>
        /// <param name="rootFindingAlgorithm">
        /// The desired root finding algorithm. By default,
        /// <see cref="DefaultRootFindingAlgorithm"/> will be used.
        /// </param>
        public LineSegment(
            int spatialDimension,
            RefElement Kräf,
            double[] start,
            double[] end,
            int iVertexStart = -1,
            int iVertexEnd = -1,
            IRootFindingAlgorithm rootFindingAlgorithm = null) {

            if (start.Length != spatialDimension)
                throw new ArgumentException();
            if (end.Length != spatialDimension)
                throw new ArgumentException();

            this.m_Kref = Kräf;
            this.iVertexStart = iVertexStart;
            this.iVertexEnd = iVertexEnd;
            this.Start = start;
            this.End = end;
            RootFindingAlgorithm = rootFindingAlgorithm ?? DefaultRootFindingAlgorithm;
        }

        /// <summary>
        /// Spatial dimension of the computational domain.
        /// </summary>
        public int SpatialDimension {
            get {
                return Start.Length;
            }
        }

        /// <summary>
        /// Start coordinates
        /// </summary>
        public double[] Start {
            get;
            private set;
        }

        /// <summary>
        /// Optional vertex index of the start point in some list of vertices
        /// containing <see cref="Start"/>. Used by some callers to evade
        /// equality comparisons of double values.
        /// </summary>
        public int iVertexStart {
            get;
            private set;
        }

        /// <summary>
        /// End coordinates
        /// </summary>
        public double[] End {
            get;
            private set;
        }

        /// <summary>
        /// Optional vertex index of the end point in some list of vertices
        /// containing <see cref="End"/>. Used by some callers to evade
        /// equality comparisons of double values.
        /// </summary>
        public int iVertexEnd {
            get;
            private set;
        }

        /// <summary>
        /// The L2 distance between <see cref="Start"/> and <see cref="End"/>
        /// </summary>
        public double Length {
            get {
                if (!length.HasValue) {
                    length = Start.L2Distance(End);
                }
                return length.Value;
            }
        }

        /// <summary>
        /// The coefficients of the polynomials that represents the projection
        /// of a set of <see cref="SpatialDimension"/>-dimensional polynomials
        /// onto this line segment (see <see cref="ProjectBasisPolynomials"/>)
        /// <list type="bullet">
        ///     <item>1st index: reference element index</item>
        ///     <item>2nd index: polynomial index </item>
        ///     <item>3rd index: coefficient index</item>
        /// </list>
        /// </summary>
        public MultidimensionalArray ProjectedPolynomialCoefficients {
            get;
            private set;
        }

        /// <summary>
        /// Splits this segment at each point defined by
        /// <paramref name="splitPoints"/>.
        /// </summary>
        /// <param name="splitPoints">
        /// The split points given as a parameter
        /// \f$ t \in ]-1, 1[\f$  of a parametrization that varies
        /// linearly between <see cref="Start"/> and <see cref="End"/>. In
        /// particular, t=-1 is the start point and t=1 the end point.
        /// </param>
        /// <returns>
        /// A list of (<paramref name="splitPoints"/>.Length + 1) consecutive
        /// line segments.
        /// </returns>
        public LineSegment[] Split(double[] splitPoints) {
            LineSegment[] subSegments = new LineSegment[splitPoints.Length + 1];
            double[] currentStart;
            double[] currentEnd = Start;

            for (int i = 0; i < splitPoints.Length; i++) {
                currentStart = currentEnd;
                currentEnd = GetPointOnSegment(splitPoints[i]);

                subSegments[i] = new LineSegment(SpatialDimension, this.m_Kref, currentStart, currentEnd);
                subSegments[i].StartCoord = i > 0 ? splitPoints[i - 1] : -1.0;
                subSegments[i].EndCoord = splitPoints[i];
            }

            subSegments[splitPoints.Length] = new LineSegment(SpatialDimension, this.m_Kref, currentEnd, this.End);
            if (splitPoints.Length > 0) {
                subSegments[splitPoints.Length].StartCoord = splitPoints[splitPoints.Length - 1];
            } else {
                subSegments[splitPoints.Length].StartCoord = -1.0;
            }
            subSegments[splitPoints.Length].EndCoord = 1.0;

            return subSegments;
        }

        /// <summary>
        /// The start coordinate expressed in terms of the parametrization
        /// described in <see cref="Split"/> of the <b>parent</b> element 
        /// </summary>
        internal double StartCoord;

        /// <summary>
        /// The end coordinate expressed in terms of the parametrization
        /// described in <see cref="Split"/> of the <b>parent</b> element
        /// </summary>
        internal double EndCoord;

        /// <summary>
        /// Evaluates the <see cref="SpatialDimension"/>-dimensional
        /// coordinate of the point with parameter value t.
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public double[] GetPointOnSegment(double t) {
            Debug.Assert(
                t.Abs() <= 1.0,
                "Point out of range");

            double[] result = new double[SpatialDimension];
            for (int d = 0; d < SpatialDimension; d++) {
                result[d] = 0.5 * ((End[d] - Start[d]) * t + Start[d] + End[d]);
            }
            return result;
        }

        /// <summary>
        /// Determines the parameter value t for a given <paramref name="pt"/>.
        /// </summary>
        /// <param name="pt"></param>
        /// <returns></returns>
        public double GetSegmentCoordinateForPoint(double[] pt) {
            if (pt.Length != this.SpatialDimension)
                throw new ArgumentException();
            var St = this.Start;
            var En = this.End;
            double nom = 0, denom = 0;
            for (int d = this.SpatialDimension - 1; d >= 0; d--) {
                double a = En[d] - St[d];
                nom += (pt[d] - St[d]) * a;
                denom += a * a;
            }
            return 2.0 * (nom / denom) - 1.0;
        }

        /// <summary>
        /// Returns a reversed line segment
        /// </summary>
        public LineSegment Inverse {
            get {
                LineSegment inverse = new LineSegment(
                    this.Start.Length,
                    this.m_Kref,
                    this.End.CloneAs(),
                    this.Start.CloneAs(),
                    this.iVertexEnd,
                    this.iVertexStart,
                    this.RootFindingAlgorithm);
                if (this.ProjectedPolynomialCoefficients != null) {
                    inverse.ProjectedPolynomialCoefficients =
                        this.ProjectedPolynomialCoefficients.CloneAs();
                }

                return inverse;
            }
        }

        /// <summary>
        /// Computes the projection of all
        /// <see cref="SpatialDimension"/>-dimensional basis polynomials in
        /// <paramref name="basis"/> onto this line segment and stores the
        /// result in <see cref="ProjectedPolynomialCoefficients"/>
        /// </summary>
        /// <param name="basis">
        /// A basis containing the polynomials to be projected
        /// </param>
        /// <remarks>
        /// Don't ask me (B. Müller) how it works, so please don't touch it. I
        /// remember that I coded it but even at that time, I didn't fully
        /// understand <b>why</b> it works.
        /// </remarks>
        public void ProjectBasisPolynomials(Basis basis) {
            int NoOfRefElm = 1;
            int iKref = Array.IndexOf(basis.GridDat.iGeomCells.RefElements, this.m_Kref);
            int noOfPolynomials = basis.Polynomials[iKref].Count;
            int noOfCoefficientsPerDimension = basis.Degree + 1;

            MultidimensionalArray[] T = new MultidimensionalArray[SpatialDimension];
            ProjectedPolynomialCoefficients = MultidimensionalArray.Create(NoOfRefElm, noOfPolynomials, noOfCoefficientsPerDimension);
            
            // Construct transformations
            for (int d = 0; d < SpatialDimension; d++) {
                double a = 0.5 * (End[d] - Start[d]);
                double b = 0.5 * (Start[d] + End[d]);

                T[d] = MultidimensionalArray.Create(noOfCoefficientsPerDimension, noOfCoefficientsPerDimension);
                for (int i = 0; i < noOfCoefficientsPerDimension; i++) {
                    for (int j = 0; j < noOfCoefficientsPerDimension; j++) {
                        if (i > j) {
                            continue;
                        }

                        T[d][i, j] = j.Choose(i) * Math.Pow(a, i) * Math.Pow(b, j - i);
                    }
                }
            }

            {
                for (int p = 0; p < noOfPolynomials; p++) {
                    Polynomial currentPolynomial = basis.Polynomials[iKref][p];
                    double[] coefficients = new double[noOfCoefficientsPerDimension];

                    // Transform coefficients to D-dimensional "matrix"
                    MultidimensionalArray originalCoefficients = MultidimensionalArray.Create(
                        Enumerable.Repeat(noOfCoefficientsPerDimension, SpatialDimension).ToArray());
                    for (int j = 0; j < currentPolynomial.Coeff.Length; j++) {
                        int[] exponents = ArrayTools.GetRow(currentPolynomial.Exponents, j);
                        originalCoefficients[exponents] += currentPolynomial.Coeff[j];
                    }

                    // Do projection
                    switch (SpatialDimension) {
                        case 1:
                            T[0].gemv(1.0, originalCoefficients.Storage, 0.0, coefficients);
                            break;

                        case 2:
                            MultidimensionalArray coefficientMatrix =
                                MultidimensionalArray.Create(noOfCoefficientsPerDimension, noOfCoefficientsPerDimension);
                            coefficientMatrix.Set(originalCoefficients);
                            coefficientMatrix = T[1] * (T[0] * coefficientMatrix).TransposeInPlace();

                            // Only left upper triangle can possibly be populated
                            for (int i = 0; i < noOfCoefficientsPerDimension; i++) {
                                for (int j = 0; j < i + 1; j++) {
                                    coefficients[i] += coefficientMatrix[i - j, j];
                                }
                            }
                            break;

                        case 3:
                            MultidimensionalArray matrix = MultidimensionalArray.Create(
                                noOfCoefficientsPerDimension, noOfCoefficientsPerDimension);

                            MultidimensionalArray[] Ty = new MultidimensionalArray[noOfCoefficientsPerDimension];
                            for (int i = 0; i < noOfCoefficientsPerDimension; i++) {
                                matrix.Set(originalCoefficients.ExtractSubArrayShallow(-1, -1, i));
                                Ty[i] = T[1] * (T[0] * matrix).TransposeInPlace();
                            }

                            // Only left upper triangle can possibly be populated
                            MultidimensionalArray tempCoefficients = MultidimensionalArray.Create(noOfCoefficientsPerDimension, Ty[0].NoOfRows);
                            for (int i = 0; i < Ty[0].NoOfRows; i++) {
                                for (int j = 0; j < i + 1; j++) {
                                    int index = i - j;

                                    for (int k = 0; k < noOfCoefficientsPerDimension; k++) {
                                        tempCoefficients[k, i] += Ty[k][index, j];
                                    }
                                }
                            }

                            MultidimensionalArray transformedTempCoefficients = T[2] * tempCoefficients;
                            for (int i = 0; i < noOfCoefficientsPerDimension; i++) {
                                for (int j = 0; j < i + 1; j++) {
                                    coefficients[i] += transformedTempCoefficients[i - j, j];
                                }
                            }
                            break;

                        default:
                            throw new ApplicationException("Invalid spatial dimension");
                    }

                    ProjectedPolynomialCoefficients.ExtractSubArrayShallow(0, p, -1).AccVector(1.0, coefficients);
                }
            }
        }

        /// <summary>
        /// Determines the length of the segment in cell
        /// <paramref name="cell"/> of the physical domain.
        /// </summary>
        /// <param name="context"></param>
        /// <param name="cell"></param>
        /// <returns></returns>
        public double GetPhysicalLength(GridData context, int cell) {
            if (!context.Cells.IsCellAffineLinear(cell))
                throw new NotImplementedException("todo: nonlinear cell implementation");

            MultidimensionalArray startGlobal = MultidimensionalArray.Create(1, 1, SpatialDimension);
            context.TransformLocal2Global(
                //MultidimensionalArray.CreateWrapper(Start, 1, SpatialDimension),
                new NodeSet(this.m_Kref, this.Start),
                cell,
                1,
                startGlobal,
                0);

            MultidimensionalArray endGlobal = MultidimensionalArray.Create(1, 1, SpatialDimension);
            context.TransformLocal2Global(
                //MultidimensionalArray.CreateWrapper(End, 1, SpatialDimension),
                new NodeSet(this.m_Kref, this.End),
                cell,
                1,
                endGlobal,
                0);

            return startGlobal.Storage.L2Distance(endGlobal.Storage);
        }

        /// <summary>
        /// The root-finding algorithm to determine roots on this segment.
        /// </summary>
        public IRootFindingAlgorithm RootFindingAlgorithm {
            get;
            private set;
        }

        /// <summary>
        /// Determines all roots of <paramref name="levelSet"/> on this segment
        /// in cell <paramref name="cell"/>.
        /// </summary>
        /// <param name="levelSet">
        /// The level set whose roots are of interest
        /// </param>
        /// <param name="cell">
        /// The cell where the evaluation takes place
        /// </param>
        /// <param name="iKref">
        /// Reference element index.
        /// </param>
        /// <returns>
        /// The parameter t of each root of <paramref name="levelSet"/> on this
        /// segment in cell <paramref name="cell"/>
        /// </returns>
        public double[] GetRoots(ILevelSet levelSet, int cell, int iKref) {
            return RootFindingAlgorithm.GetRoots(this, levelSet, cell, iKref);

            // Caching currently deactivated since it didn't really help...
            //Tuple<ILevelSet, int> key = Tuple.Create(levelSet, cell);
            //double[] roots;
            //if (!rootsCache.TryGetValue(key, out roots)) {
            //    roots = RootFindingAlgorithm.GetRoots(this, levelSet, cell);
            //    rootsCache.Add(key, roots);
            //}
            //return roots;
        }

        #region IEquatable<LineSegment> Members

        /// <summary>
        /// Segments are considered equal if the start and end coordinates are
        /// equal (w.r.t. <see cref="EPSILON"/>). Here, the direction of the
        /// segment is ignored which means that two segments are considered
        /// equal even if start and points are exchanged.
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(LineSegment other) {
            if (other == null) {
                return false;
            }

            return (EqualsImp(other) || EqualsImp(other.Inverse));
        }

        /// <summary>
        /// Implementation of <see cref="Equals"/>
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        private bool EqualsImp(LineSegment other) {
            if (Start.RoughlyEquals(other.Start, EPSILON)
                && End.RoughlyEquals(other.End, EPSILON)) {
                return true;
            } else {
                return false;
            }
        }

        #endregion

        /// <summary>
        /// Interface for approximate one-dimensional root-finding algorithms.
        /// </summary>
        public interface IRootFindingAlgorithm {

            /// <summary>
            /// The desired tolerance of the root-finding algorithm
            /// </summary>
            double Tolerance {
                get;
            }

            /// <summary>
            /// Determines the roots of <paramref name="levelSet"/> on the
            /// given <paramref name="segment"/> in cell <paramref name="cell"/>
            /// </summary>
            /// <param name="segment"></param>
            /// <param name="levelSet"></param>
            /// <param name="cell"></param>
            /// <param name="iKref"></param>
            /// <returns>
            /// Determines the parameter t for each root of
            /// <paramref name="levelSet"/> on <paramref name="segment"/> in
            /// the given <paramref name="cell"/>.
            /// </returns>
            double[] GetRoots(LineSegment segment, ILevelSet levelSet, int cell, int iKref);
        }

        /// <summary>
        /// A root-finding algorithm making use of the GNU Scientific Library
        /// (GSL; cf. <see cref="GSL"/>).
        /// </summary>
        public class GSLRootFindingAlgorithm : IRootFindingAlgorithm {

            /// <summary>
            /// Constructs a root-finder with the given
            /// <paramref name="tolerance"/>
            /// </summary>
            /// <param name="tolerance"></param>
            public GSLRootFindingAlgorithm(double tolerance) {
                this.Tolerance = tolerance;
            }

            #region IRootFindingAlgorithm Members

            /// <summary>
            /// Tolerance of the algorithm
            /// </summary>
            public double Tolerance {
                get;
                private set;
            }

            /// <summary>
            /// Finds the roots using <see cref="GSL.gsl_poly_complex_solve"/>
            /// </summary>
            /// <param name="segment"></param>
            /// <param name="levelSet"></param>
            /// <param name="cell"></param>
            /// <param name="iKref"></param>
            /// <returns></returns>
            public double[] GetRoots(LineSegment segment, ILevelSet levelSet, int cell, int iKref) {
                LevelSet levelSetField = levelSet as LevelSet;
                if (levelSetField == null) {
                    throw new NotImplementedException("Method currently only works for polynomial level sets");
                }

                int maxNoOfCoefficientsPerDimension = levelSetField.Basis.Degree + 1;
                int noOfPolynomials = segment.ProjectedPolynomialCoefficients.GetLength(1);

                double[] coefficients = new double[maxNoOfCoefficientsPerDimension];
                for (int i = 0; i < noOfPolynomials; i++) {
                    double dgCoefficient = levelSetField.Coordinates[cell, i];
                    for (int j = 0; j < maxNoOfCoefficientsPerDimension; j++) {
                        coefficients[j] += dgCoefficient * segment.ProjectedPolynomialCoefficients[iKref, i, j];
                    }
                }

                // Make sure "leading" coefficient (i.e., the last element of the
                // list of coefficients) is not too small since this will make gsl
                // crash (or lead to bogus results)
                int newLength = coefficients.Length;
                while (newLength > 0 && Math.Abs(coefficients[newLength - 1]) < EPSILON) {
                    newLength--;
                }
                if (newLength != coefficients.Length) {
                    coefficients = coefficients.Take(newLength).ToArray();
                }

                // Make sure polynomial is not constant since this will make gsl crash
                if (newLength < 2) {
                    return new double[0];
                }

                int maxNoOfRoots = coefficients.Length - 1;
                double[] roots = new double[2 * maxNoOfRoots];

                GCHandle coefficientHandle = GCHandle.Alloc(coefficients, GCHandleType.Pinned);
                GCHandle rootsHandle = GCHandle.Alloc(roots, GCHandleType.Pinned);

                IntPtr workSpace = GSL.gsl_poly_complex_workspace_alloc(coefficients.Length);
                GSL.gsl_poly_complex_solve(
                    coefficientHandle.AddrOfPinnedObject(),
                    coefficients.Length,
                    workSpace,
                    rootsHandle.AddrOfPinnedObject());
                GSL.gsl_poly_complex_workspace_free(workSpace);

                coefficientHandle.Free();
                rootsHandle.Free();

                List<double> realRoots = new List<double>(levelSetField.Basis.Degree);
                for (int i = 0; i < maxNoOfRoots; i++) {
                    double realPart = roots[2 * i];
                    double imaginaryPart = roots[2 * i + 1];

                    // Exclude |$realPart| == 1.0 since it doesn't matter
                    if (imaginaryPart == 0.0 && Math.Abs(realPart) <= 1.0 - EPSILON) {
                        realRoots.Add(realPart);
                    }
                }

                return realRoots.OrderBy(d => d).ToArray();
            }

            #endregion
        }

        /// <summary>
        /// Implementation of a SIMPLE safe-guarded Newton algorithm following
        /// Press et al.: Numerical recipes (Chapter 9.4)
        /// </summary>
        public class SafeGuardedNewtonMethod : IRootFindingAlgorithm {

            /// <summary>
            /// Constructs a new root-finder
            /// </summary>
            /// <param name="tolerance"></param>
            public SafeGuardedNewtonMethod(double tolerance) {
                this.Tolerance = tolerance;
            }

            #region IRootFindingAlgorithm Members

            /// <summary>
            /// The selected tolerance
            /// </summary>
            public double Tolerance {
                get;
                private set;
            }

            /// <summary>
            /// Uses a safe-guarded Newton method to find all roots on
            /// <paramref name="segment"/>
            /// </summary>
            /// <param name="segment"></param>
            /// <param name="levelSet"></param>
            /// <param name="cell"></param>
            /// <param name="iKref"></param>
            /// <returns></returns>
            public double[] GetRoots(LineSegment segment, ILevelSet levelSet, int cell, int iKref) {
                LevelSet levelSetField = levelSet as LevelSet;
                if (levelSetField == null) {
                    throw new NotImplementedException("Method currently only works for polynomial level sets");
                }

                int maxNoOfCoefficientsPerDimension = levelSetField.Basis.Degree + 1;
                int noOfPolynomials = segment.ProjectedPolynomialCoefficients.GetLength(1);

                double[] coefficients = new double[maxNoOfCoefficientsPerDimension];
                for (int i = 0; i < noOfPolynomials; i++) {
                    double dgCoefficient = levelSetField.Coordinates[cell, i];
                    for (int j = 0; j < maxNoOfCoefficientsPerDimension; j++) {
                        coefficients[j] += dgCoefficient * segment.ProjectedPolynomialCoefficients[iKref, i, j];
                    }
                }
                //if(coefficients.L2NormPow2() < this.Tolerance) {
                //    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                //    // special case :
                //    // the zero-level-set is probably parallel to this line segment
                //    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                //    return new double[] { -1.0, 1.0 };
                //}

                double[] roots;
                unsafe {
                    fixed (double* pCoeff = &coefficients[coefficients.Length - 1]) {
                        double xLow = -1.0;
                        double xHigh = 1.0;

                        int NO_OF_BRACKETS = 8;
                        List<double> lowerBounds = new List<double>();
                        List<double> upperBounds = new List<double>();
                        double dx2 = 2.0 / NO_OF_BRACKETS;
                        double x = xLow;
                        double fOld = Eval(x, pCoeff, coefficients.Length);
                        for (int i = 0; i < NO_OF_BRACKETS; i++) {
                            x += dx2;
                            double fNew = Eval(x, pCoeff, coefficients.Length);

                            if (i == 0 && fOld.Abs() < Tolerance) {
                                lowerBounds.Add(x - dx2);
                                upperBounds.Add(x);
                                fOld = fNew;
                                continue;
                            } else {
                                if (fNew.Abs() < Tolerance) {
                                    lowerBounds.Add(x - dx2);
                                    upperBounds.Add(x);

                                    x += dx2;
                                    fOld = Eval(x, pCoeff, coefficients.Length);
                                    continue;
                                }
                            }

                            if (fNew * fOld <= 0.0) {
                                lowerBounds.Add(x - dx2);
                                upperBounds.Add(x);
                            }
                            fOld = fNew;
                        }

                        // Actual Newton-Raphson
                        int MAX_ITERATIONS = 50;

                        roots = new double[lowerBounds.Count];
                        for (int j = 0; j < lowerBounds.Count; j++) {
                            xLow = lowerBounds[j];
                            xHigh = upperBounds[j];

                            double fLeft = Eval(xLow, pCoeff, coefficients.Length);
                            double fRight = Eval(xHigh, pCoeff, coefficients.Length);

                            if (fLeft.Abs() < Tolerance) {
                                roots[j] = xLow;
                                break;
                            }

                            if (fRight.Abs() < Tolerance) {
                                roots[j] = xHigh;
                                break;
                            }

                            if (fLeft.Sign() == fRight.Sign()) {
                                throw new Exception();
                            }

                            if (fLeft > 0.0) {
                                xLow = upperBounds[j];
                                xHigh = lowerBounds[j];
                            }

                            double root = 0.5 * (xLow + xHigh);
                            double f;
                            double df;
                            Eval(root, pCoeff, coefficients.Length, out f, out df);

                            double dxOld = (xHigh - xLow).Abs();
                            double dx = dxOld;

                            int i = 0;
                            while (true) {
                                if (i > MAX_ITERATIONS) {
                                    throw new Exception("Max iterations exceeded");
                                }

                                double a = ((root - xHigh) * df - f) * ((root - xLow) * df - f);
                                if (a > 0.0 || 2.0 * Math.Abs(f) > Math.Abs(dxOld * df)) {
                                    // Newton out of range or too slow -> Bisect
                                    dxOld = dx;
                                    dx = 0.5 * (xHigh - xLow);
                                    root = xLow + dx;
                                } else {
                                    // Take Newton step
                                    dxOld = dx;
                                    dx = -f / df;

                                    //// Convergence acceleration according to Yao2014
                                    //double fDelta = Eval(root + dx, pCoeff, coefficients.Length);
                                    //dx = -(f + fDelta) / df;

                                    root += dx;
                                }

                                Eval(root, pCoeff, coefficients.Length, out f, out df);

                                if (Math.Abs(f) <= Tolerance) {
                                    roots[j] = root;
                                    break;
                                }

                                if (f < 0.0) {
                                    xLow = root;
                                } else {
                                    xHigh = root;
                                }

                                i++;
                            }
                        }
                    }
                }

                return roots;
            }

            #endregion

            /// <summary>
            /// Evaluates the polynomial p with coefficients stored in an array
            /// <b>ending</b> at <paramref name="pCoeff"/>
            /// </summary>
            /// <param name="x">One-dimensional coordinate</param>
            /// <param name="pCoeff">
            /// Pointer to the last coefficient of a polynomial (i.e., the one
            /// with the highest exponent)
            /// </param>
            /// <param name="numberOfCoefficients">
            /// The number of coefficients
            /// </param>
            /// <returns>
            /// p(x)
            /// </returns>
            private unsafe double Eval(double x, double* pCoeff, int numberOfCoefficients) {
                double p = *pCoeff;
                for (int i = 1; i < numberOfCoefficients; i++) {
                    p = p * x + *(--pCoeff);
                }
                return p;
            }

            /// <summary>
            /// Evaluates the polynomial p with coefficients stored in an array
            /// <b>ending</b> at <paramref name="pCoeff"/>, as well as its
            /// derivative p' at <paramref name="x"/>
            /// </summary>
            /// <param name="x">One-dimensional coordinate</param>
            /// <param name="pCoeff">
            /// Pointer to the last coefficient of a polynomial (i.e., the one
            /// with the highest exponent)
            /// </param>
            /// <param name="numberOfCoefficients">
            /// The number of coefficients
            /// </param>
            /// <param name="p">p(x)</param>
            /// <param name="dp">p'(x)</param>
            private unsafe void Eval(double x, double* pCoeff, int numberOfCoefficients, out double p, out double dp) {
                p = *pCoeff;
                dp = 0.0;
                for (int i = 1; i < numberOfCoefficients; i++) {
                    dp = dp * x + p;
                    p = p * x + *(--pCoeff);
                }
            }
        }


        #region IObserver<LevelSetStatus> Members

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void OnCompleted() {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void OnError(Exception error) {
        }

        /// <summary>
        /// Clear cache
        /// </summary>
        /// <param name="status">
        /// Not used.
        /// </param>
        public void OnNext(LevelSetTracker.LevelSetRegions status) {
            rootsCache.Clear();
        }

        #endregion
    }
}
