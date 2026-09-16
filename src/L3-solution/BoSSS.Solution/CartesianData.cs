using ilPSP;
using ilPSP.Utils;
using System;
using System.Runtime.Serialization;

namespace BoSSS.Solution.Control {

    /// <summary>
    /// This class allows to interpolate boundarey/intitial data interpolated from a cartesian mesh
    /// </summary>
    [DataContract]
    [Serializable]
    public class CartesianData : IBoundaryAndInitialData {

        [DataMember]
        private MultidimensionalArray m_CartesianFieldData;

        [DataMember]
        private double[][] m_coordNodes;

        [DataMember]
        private double[] m_spacing;

        [DataMember]
        private bool[] m_isEquidistant;

        /// <summary>
        /// Private ctor for 
        /// </summary>
        private CartesianData() {

        }

        /// <summary>
        /// Special constructor, load 2D data from a textfile.
        /// </summary>
        public CartesianData(string FieldData_file, string xNodes_file, string yNodes_file) //
            : this(IMatrixExtensions.LoadFromTextFile(FieldData_file),
                  new double[][] {
                      VectorIO.LoadFromTextFile(xNodes_file),
                      VectorIO.LoadFromTextFile(yNodes_file)
                  }) {

        }

        /// <summary>
        /// Constructor
        /// </summary>
        public CartesianData(MultidimensionalArray CartesianFieldData, double[][] coordNodes) {
            if(CartesianFieldData == null)
                throw new ArgumentNullException(nameof(CartesianFieldData));
            if(coordNodes == null)
                throw new ArgumentNullException(nameof(coordNodes));
            if(CartesianFieldData.Dimension < 1 || CartesianFieldData.Dimension > 3)
                throw new ArgumentException("Cartesian data must have one, two, or three dimensions.", nameof(CartesianFieldData));
            if(CartesianFieldData.Dimension != coordNodes.Length)
                throw new ArgumentException($"expecting equal number of node vectors (got {coordNodes.Length}) and dimension of field data (is {CartesianFieldData.Dimension})");

            m_CartesianFieldData = CartesianFieldData.CloneAs();
            m_coordNodes = new double[coordNodes.Length][];
            m_spacing = new double[coordNodes.Length];
            m_isEquidistant = new bool[coordNodes.Length];

            for(int iDim = 0; iDim < CartesianFieldData.Dimension; iDim++) {
                if(coordNodes[iDim] == null)
                    throw new ArgumentException($"node vector for dimension {iDim} must not be null", nameof(coordNodes));
                if(coordNodes[iDim].Length == 0)
                    throw new ArgumentException($"node vector for dimension {iDim} must not be empty", nameof(coordNodes));
                if(CartesianFieldData.GetLength(iDim) != coordNodes[iDim].Length)
                    throw new ArgumentException($"dimension {iDim}, mismatch between node vector length ({coordNodes[iDim].Length}) and data dimension ({CartesianFieldData.GetLength(iDim)})");

                m_coordNodes[iDim] = coordNodes[iDim].CloneAs();
                for(int iNode = 0; iNode < m_coordNodes[iDim].Length; iNode++) {
                    if(double.IsNaN(m_coordNodes[iDim][iNode]) || double.IsInfinity(m_coordNodes[iDim][iNode]))
                        throw new ArgumentException($"node vector for dimension {iDim} contains a non-finite value", nameof(coordNodes));
                    if(iNode > 0 && m_coordNodes[iDim][iNode] <= m_coordNodes[iDim][iNode - 1])
                        throw new ArgumentException($"node vector for dimension {iDim} must be strictly increasing", nameof(coordNodes));
                }

                if(m_coordNodes[iDim].Length > 1) {
                    double firstSpacing = m_coordNodes[iDim][1] - m_coordNodes[iDim][0];
                    double tolerance = Math.Max(1.0, Math.Max(Math.Abs(m_coordNodes[iDim][0]), Math.Abs(m_coordNodes[iDim][m_coordNodes[iDim].Length - 1])))
                        * GenericBlas.MachineEps * 1000.0;
                    bool equidistant = true;
                    for(int iNode = 2; iNode < m_coordNodes[iDim].Length; iNode++) {
                        if(Math.Abs((m_coordNodes[iDim][iNode] - m_coordNodes[iDim][iNode - 1]) - firstSpacing) > tolerance) {
                            equidistant = false;
                            break;
                        }
                    }
                    m_spacing[iDim] = firstSpacing;
                    m_isEquidistant[iDim] = equidistant;
                }
            }
        }

        private static void FindInterval(double coordinate, double[] nodes, bool equidistant, double spacing, out int lower, out double fraction) {
            if(double.IsNaN(coordinate) || double.IsInfinity(coordinate))
                throw new ArgumentException("Evaluation coordinates must be finite.", nameof(coordinate));

            if(nodes.Length == 1 || coordinate <= nodes[0]) {
                lower = 0;
                fraction = 0.0;
                return;
            }
            if(coordinate >= nodes[nodes.Length - 1]) {
                lower = nodes.Length - 1;
                fraction = 0.0;
                return;
            }

            int found;
            if(equidistant) {
                lower = (int)Math.Floor((coordinate - nodes[0]) / spacing);
                if(lower < 0)
                    lower = 0;
                if(lower >= nodes.Length - 1)
                    lower = nodes.Length - 2;
                found = lower;
            } else {
                found = Array.BinarySearch(nodes, coordinate);
                if(found >= 0) {
                    lower = found;
                    fraction = 0.0;
                    return;
                }
                lower = ~found - 1;
            }

            fraction = (coordinate - nodes[lower]) / (nodes[lower + 1] - nodes[lower]);
        }

        private double Interpolate(double[] coordinates, int[] lower, double[] fraction) {
            for(int iDim = 0; iDim < m_coordNodes.Length; iDim++)
                FindInterval(coordinates[iDim], m_coordNodes[iDim], m_isEquidistant[iDim], m_spacing[iDim], out lower[iDim], out fraction[iDim]);

            int next0 = lower[0] + (m_coordNodes[0].Length > 1 ? 1 : 0);
            int next1 = lower[1] + (m_coordNodes.Length > 1 && m_coordNodes[1].Length > 1 ? 1 : 0);
            int next2 = lower[2] + (m_coordNodes.Length > 2 && m_coordNodes[2].Length > 1 ? 1 : 0);
            switch(m_coordNodes.Length) {
                case 1:
                    if(m_coordNodes[0].Length == 1)
                        return m_CartesianFieldData[0];
                    return m_CartesianFieldData[lower[0]] * (1.0 - fraction[0])
                        + m_CartesianFieldData[lower[0] + 1] * fraction[0];
                case 2:
                    double v00 = m_CartesianFieldData[lower[0], lower[1]];
                    double v10 = m_CartesianFieldData[next0, lower[1]];
                    double v01 = m_CartesianFieldData[lower[0], next1];
                    double v11 = m_CartesianFieldData[next0, next1];
                    double x0 = v00 * (1.0 - fraction[0]) + v10 * fraction[0];
                    double x1 = v01 * (1.0 - fraction[0]) + v11 * fraction[0];
                    return x0 * (1.0 - fraction[1]) + x1 * fraction[1];
                default:
                    double value = 0.0;
                    for(int corner = 0; corner < 8; corner++) {
                        int i0 = lower[0] + ((corner & 1) == 0 ? 0 : next0 - lower[0]);
                        int i1 = lower[1] + ((corner & 2) == 0 ? 0 : next1 - lower[1]);
                        int i2 = lower[2] + ((corner & 4) == 0 ? 0 : next2 - lower[2]);
                        double weight = ((corner & 1) == 0 ? 1.0 - fraction[0] : fraction[0])
                            * ((corner & 2) == 0 ? 1.0 - fraction[1] : fraction[1])
                            * ((corner & 4) == 0 ? 1.0 - fraction[2] : fraction[2]);
                        value += m_CartesianFieldData[i0, i1, i2] * weight;
                    }
                    return value;
            }
        }

        public double Evaluate(double[] X, double t) {
            if(X == null)
                throw new ArgumentNullException(nameof(X));
            if(X.Length != m_coordNodes.Length)
                throw new ArgumentException($"Expected {m_coordNodes.Length} coordinates, got {X.Length}.", nameof(X));
            return Interpolate(X, new int[3], new double[3]);
        }

        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            if(input == null)
                throw new ArgumentNullException(nameof(input));
            if(output == null)
                throw new ArgumentNullException(nameof(output));
            if(input.Dimension != 2 || input.GetLength(1) != m_coordNodes.Length)
                throw new ArgumentException($"Input must have dimensions (number of points, {m_coordNodes.Length}).", nameof(input));
            if(output.Dimension != 1 || output.GetLength(0) != input.GetLength(0))
                throw new ArgumentException("Output must have one entry for each input point.", nameof(output));

            int[] lower = new int[3];
            double[] fraction = new double[3];
            double[] coordinates = new double[m_coordNodes.Length];
            for(int iPoint = 0; iPoint < input.GetLength(0); iPoint++) {
                for(int iDim = 0; iDim < coordinates.Length; iDim++)
                    coordinates[iDim] = input[iPoint, iDim];
                output[iPoint] = Interpolate(coordinates, lower, fraction);
            }
        }



    }
}
