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
using System.Text;
using System.Xml;
using System.Globalization;

namespace BoSSS.Solution.Utils.Spline {

    /// <summary>
    /// A cubic B-Spline, which can be used for encoding
    /// boundary values.
    /// </summary>
    /// <remarks>
    /// The B-Spline depends only on one spatial coordinate and is 
    /// constant in the other (two) directions, i.e.
    /// it is a mapping <br/>
    /// (x,y,z)  ->  f(x) <br/>
    /// or <br/>
    /// (x,y,z)  ->  f(y) <br/>
    /// or <br/>
    /// (x,y,z)  ->  f(z). <br/>
    /// </remarks>
    [Serializable]
    public class Spline {

        double[] M;
        private double[] c;
        private double[] d;
        private int dim;
        private double[] nodes;
        private double[] values;

        /// <summary>
        /// defines the behavior of the spline when its argument is out of bounds,
        /// i.e. when x (or y or z) is lower than the 0-th or greater than the last x-Node (or y- or z-Node).
        /// </summary>
        public enum OutOfBoundsBehave {
            /// <summary>
            /// first and last polynomial will be extrapolated
            /// </summary>
            Extrapolate,

            /// <summary>
            /// the argument is clamped, i.e. the evaluation of the spline will give
            /// the value at the 0-th or last node, respectively. 
            /// </summary>
            Clamp
        }

        private OutOfBoundsBehave m_oobb;
             

        /// <summary>
        /// Constructs a B-Spline based on a set of <paramref name="nodes"/> and <paramref name="values"/>;
        /// </summary>
        /// <param name="nodes">
        /// length must be at least 2, and equal to the length of <paramref name="values"/>;
        /// </param>
        /// <param name="values">
        /// length must be at least 2, and equal to the length of <paramref name="nodes"/>;
        /// </param>
        /// <param name="d">
        /// Coordinate axis, on which the <paramref name="nodes"/>
        /// <list type="bullet">
        /// <item>0: <paramref name="nodes"/> are values on the x - Axis</item>
        /// <item>1: <paramref name="nodes"/> are values on the y - Axis</item>
        /// <item>2: <paramref name="nodes"/> are values on the z - Axis</item>
        /// </list>
        /// </param>
        /// <param name="oobb">
        /// </param>
        public Spline(double[] nodes, double[] values, int d, OutOfBoundsBehave oobb) {
            M = new double[nodes.Length];
            this.nodes = nodes;
            this.values = values;
            dim = d;
            m_oobb = oobb;
            double[] h = new double[nodes.Length - 1];
            double[,] A = new double[nodes.Length, 3];         //stores the coefficients of the equation system with unknown M (A*M = b)
            double[] b = new double[nodes.Length];             //stores the solution vector of the equation system
            ConstructorCommon(ref h, ref A, ref b);
        }

        /// <summary>
        /// Constructs a B-Spline from an XML element. The XML element is parsed in order to extract the needed information, namely
        /// the interpolation points and their values respectively.
        /// </summary>
        /// <param name="xmlelm">
        /// Input is the given XmlElement, normally representing a spline
        /// </param>
        public Spline(XmlElement xmlelm) {
            parse(xmlelm);                                      //parse the given element to extract the information needed
            M = new double[nodes.Length];
            double[] h = new double[this.nodes.Length - 1];
            double[,] A = new double[this.nodes.Length, 3];         //stores the coefficients of the equation system with unknown M (A*M = b)
            double[] b = new double[this.nodes.Length];             //stores the solution vector of the equation system
            ConstructorCommon(ref h, ref A, ref b);
        }

        /// <summary>
        /// This method is used as a private constructor by the other constructors. 
        /// It calculates the needed initial parameters and calls the function that finds the moments M,
        /// solving a linear equation system (LES)
        /// </summary>
        /// <param name="h">
        /// contains the distances between the interpolation nodes
        /// </param>
        /// <param name="A">
        /// This matrix contains the coefficients, needed to solve the LES</param>
        /// <param name="b">
        /// This is the vector from the right side in the linear equation system</param>
        private void ConstructorCommon(ref double[] h, ref double[,] A, ref double[] b) {
            //calculates the distance between each two neighboring nodes and stores it in the array h
            h = this.calculateDistanceH();

            //initialising the solution vector b. According to the natural border conditions b[n] = b[0] = 0
            b = this.calculateB(h);

            //initialising the matrix A. 
            //According to the natural border conditions A[0,2] = A[n,0] = 0 && A[0,1] = A[n,1] = 1
            A = this.initA(h);
            M = solveEquation(A, b);            //calculates the moment M using the Gauss-Elimination Method
            computeCoefficients(h);              //should compute the coefficients in the arrays d and c that are needed for computing the spline
        }

        /// <summary>
        /// Evaluates the B-spline at point <paramref name="X"/>
        /// </summary>
        /// <param name="X">
        /// The coordinate vector, either (x,y) in 2D, (x,y,z) in 3D or
        /// just (x) in 1D.
        /// </param>
        /// <returns>
        /// the calculated value for the spline at the given point</returns>
        /// <remarks>
        /// The B-Spline depends only on one spatial coordinate and is 
        /// constant in the other (two) directions, i.e.
        /// it is a mapping <br/>
        /// (x,y,z)  ->  f(x) <br/>
        /// or <br/>
        /// (x,y,z)  ->  f(y) <br/>
        /// or <br/>
        /// (x,y,z)  ->  f(z). <br/>
        /// </remarks>
        public double Evaluate(params double[] X) {

            double x = X[dim];              //dim is the global variable that is set to 0 for x, 1 for y, 2 for z 
            
            if( m_oobb == OutOfBoundsBehave.Clamp) {
                if(x <= this.nodes[0])
                    return this.values[0];

                if (x >= this.nodes[this.nodes.Length - 1])
                    return this.values[this.nodes.Length - 1];
            }

            return Evaluate1D(x);
        }

        /// <summary>
        /// Spline evaluation core
        /// </summary>
        /// <param name="x">argument</param>
        /// <returns>spline value at <paramref name="x"/></returns>
        public double Evaluate1D(double x) {
            double spline;
            int interval = binarySearch(x);   //the interval where the given point is, is calculated with the aid of the binary search
            if (interval >= this.nodes.Length - 1) {
                interval = this.nodes.Length - 2;
            }
            double xI = this.nodes[interval];
            double xI1 = this.nodes[interval + 1];
            double MI = M[interval];
            double MI1 = M[interval + 1];

            spline = ((double)1 / 6) * ((((xI1 - x) * (xI1 - x) * (xI1 - x)) / (xI1 - xI)) * MI +
                    (((x - xI) * (x - xI) * (x - xI)) / (xI1 - xI)) * MI1) + c[interval] * (x - xI) + d[interval];

            return spline;
        }

        /// <summary>
        /// this function finds the interval which x belongs to
        /// </summary>
        /// <param name="x">
        /// The node that is to be placed between the interpolation points
        /// </param>
        /// <returns>
        /// To which interval this node belongs taking into account the interpolation points
        /// </returns>
        private int binarySearch(double x) {
            int min = 0;
            int max = this.nodes.Length - 1;
            int mid = 0;
            do {
                mid = (min + max) / 2;
                if (x > this.nodes[mid]) {
                    min = mid + 1;
                } else {
                    max = mid - 1;
                }
            }
            while (mid > 0 && min <= max && (!(this.nodes[mid] >= x && this.nodes[mid - 1] < x)) && (!(this.nodes[mid] < x && this.nodes[mid + 1] >= x)));
            if (mid > 0 && this.nodes[mid] >= x && this.nodes[mid - 1] < x) {
                return (mid - 1);
            }
            return mid;
        }


        ///<summary>
        ///Solves the equation Ax = b, where A is a diagonal matrix and b is the solution vector, using backwards substitution
        ///</summary>
        ///<param name="A">
        ///This is the matrix A from the equation Ax=b that is to be solved
        ///</param>
        ///<param name="b">
        ///b is the solution vector
        ///</param>
        ///
        private double[] solveEquation(double[,] A, double[] b) {

            //manipulating the matrix in order to get only zeros under the diagonal. A[1,] should be treated separately
            //A[1,0] = 0;
            A[b.Length - 1, 0] = 0;
            double factor;
            for (int i = 1; i < b.Length; i++) {
                factor = A[i, 0] / A[i - 1, 1];
                A[i, 0] = 0;
                A[i, 1] = A[i, 1] - A[i - 1, 2] * factor;
                b[i] = b[i] - b[i - 1] * factor;
            }
            //manipulating the matrix in order to get only ones on the diagonal
            for (int i = 1; i < b.Length; i++) {
                double factor1 = A[i, 1];
                A[i, 1] = A[i, 1] / factor1;
                A[i, 2] = A[i, 2] / factor1;
                b[i] = b[i] / factor1;
            }
            //manipulating the matrix so that we get only zeros above the diagonal
            for (int i = b.Length - 2; i > 0; i--) {
                double factor2 = A[i, 2] / A[i + 1, 1];
                A[i, 2] = A[i, 2] - A[i + 1, 1] * factor2;
                b[i] = b[i] - b[i + 1] * factor2;
            }
            return b;
        }

        /// <summary>
        /// This function calculates the distance between each two neighbouring nodes and stores it in the array h
        /// </summary>
        /// <returns>
        /// the calculated array containing the distances
        /// </returns>
        private double[] calculateDistanceH() {
            double[] h = new double[this.nodes.Length - 1];
            h[0] = this.nodes[1] - this.nodes[0];
            for (int i = 1; i < nodes.Length - 1; i++) {
                h[i] = this.nodes[i + 1] - this.nodes[i];
            }
            return h;
        }

        /// <summary>
        /// This function initializes the solution vector b. According to the natural border conditions b[n] = b[0] = 0
        /// </summary>
        /// <param name="h">
        /// the distances between the input nodes are needed to compute b
        /// </param>
        /// <returns>
        /// the computed array b
        /// </returns>
        private double[] calculateB(double[] h) {
            double[] b = new double[this.nodes.Length];
            b[0] = 0;
            b[nodes.Length - 1] = 0;
            for (int i = 1; i < nodes.Length - 1; i++) {
                b[i] = (this.values[i + 1] - this.values[i]) / h[i] - (this.values[i] - this.values[i - 1]) / h[i - 1];
            }
            return b;
        }

        /// <summary>
        /// This function initializes the matrix A. According to the natural border conditions A[0,2] = A[n,0] = 0 &amp; &amp; A[0,1] = A[n,1] = 1
        /// </summary>
        /// <param name="h">
        /// </param>
        /// <returns></returns>
        private double[,] initA(double[] h) {
            double[,] A = new double[this.nodes.Length, 3];         //stores the coefficients of the equation system with unknown M (A*M = b)
            A[0, 0] = 0;
            A[nodes.Length - 1, 2] = 0;
            A[0, 1] = 1;
            A[nodes.Length - 1, 1] = 1;
            A[0, 2] = 0;
            A[nodes.Length - 1, 0] = 0;
            for (int i = 1; i < nodes.Length - 1; i++) {
                A[i, 0] = h[i - 1] / 6;
                A[i, 1] = (h[i - 1] + h[i]) / 3;
                A[i, 2] = h[i] / 6;
            }
            return A;
        }

        /// <summary>
        /// this function computes the coefficients stored in the arrays c and d, needed to evaluate the spline
        /// </summary>
        /// <param name="h">
        /// This array contains the distances between the interpolation nodes
        /// </param>
        private void computeCoefficients(double[] h) {
            d = new double[h.Length];
            c = new double[h.Length];
            for (int i = 0; i < d.Length; i++) {
                d[i] = this.values[i] - ((h[i] * h[i]) / 6) * M[i];
                c[i] = (this.values[i + 1] - this.values[i]) / h[i] - (h[i] / 6) * (M[i + 1] - M[i]);
            }
        }

        /// <summary>
        /// This function parses an XmlElement in order to extract the nodes and values for the spline
        /// </summary>
        /// <param name="xmlelem">
        /// the given XmlElement, normally the "spline" element
        /// </param>
        private void parse(XmlElement xmlelem) {
            XmlNodeList children = xmlelem.ChildNodes;
            int actualCount = 0;
            int countDepaxisComments = 0;
            double[] xTmp = new double[children.Count * 2];
            double[] yTmp = new double[children.Count * 2];
            char[] separator = new char[] { ' ', '\t', '\n', '\r' };
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            foreach (XmlNode child in children) {
                if (child.NodeType == XmlNodeType.Element) {
                    if (child.Name == "depaxis") {
                        XmlNodeList depaxisChildren = child.ChildNodes;

                        foreach (XmlNode depaxisChild in depaxisChildren) {
                            //the element depaxis cannot have any children, otherwise an exception is thrown
                            if (depaxisChild.NodeType == XmlNodeType.Comment) {
                                countDepaxisComments++;
                            } else {
                                if (depaxisChild.NodeType != XmlNodeType.Text) {
                                    throw new ApplicationException("The element depaxis should contain 'x', 'y' or 'z'!");
                                } else {
                                    string[] depChild = depaxisChild.InnerText.Split(separator, StringSplitOptions.RemoveEmptyEntries);
                                    
                                    string depAxis = depChild[0].ToLower();
                                    switch (depAxis) {
                                        case "x":
                                            this.dim = 0;
                                            break;

                                        case "y":
                                            this.dim = 1;
                                            break;

                                        case "z":
                                            this.dim = 2;
                                            break;

                                        default:
                                            throw new ApplicationException("The element depaxis should contain 'x', 'y' or 'z'!");
                                    }
                                }
                            }
                        }
                        

                        if (depaxisChildren.Count - countDepaxisComments != 1) {
                            throw new ApplicationException("The element depaxis can contain only 'x', 'y', 'z'!");
                        }
                    } else if (child.Name == "pt") {
                        // we extract the values of a pt element, storing them into xTmp and yTmp
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        string tmp = child.InnerText;
                        string[] tmpValues = tmp.Split(separator, StringSplitOptions.RemoveEmptyEntries);
                        xTmp[actualCount] = Convert.ToDouble(tmpValues[0], provider);
                        yTmp[actualCount] = Convert.ToDouble(tmpValues[1], provider);
                        actualCount = actualCount + 1;
                    } else if (child.Name == "OutOfBoundsBehave") {
                        //
                        // 

                        XmlNodeList depaxisChildren = child.ChildNodes;
                        if (depaxisChildren.Count - countDepaxisComments != 1) 
                            throw new ApplicationException("The no Childs allowed for element 'OutOfBoundsBehave';");

                        m_oobb = (OutOfBoundsBehave) Enum.Parse(typeof(OutOfBoundsBehave),child.InnerText);

                    } else {
                        throw new ApplicationException("unknown child element '" + child.Name + "' for spline;");
                    }
                }
            }
            this.nodes = new double[actualCount];
            this.values = new double[actualCount];
            double isAscending = 0;
            //copy the values from xTmp and yTmp to this.nodes and this.values 
            this.nodes[0] = xTmp[0];
            this.values[0] = yTmp[0];
            for (int i = 1; i < actualCount; i++) {
                this.nodes[i] = xTmp[i];
                this.values[i] = yTmp[i];
                isAscending = this.nodes[i] - this.nodes[i - 1];
                //if the interpolating nodes are not in ascending order, an exception is thrown
                if (isAscending <= 0) {
                    throw new ApplicationException("The nodes of the spline must be given in a strictly ascending order!");
                }
                isAscending = 0;
            }
        }
    }
}
