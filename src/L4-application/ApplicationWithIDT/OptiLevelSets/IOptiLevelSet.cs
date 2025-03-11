using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using BoSSS.Solution.Utils;
using System.Collections;
using BoSSS.Solution.Statistic;
using ilPSP.LinSolvers;

namespace ApplicationWithIDT.OptiLevelSets {
    /// <summary>
    /// Defines an interface for the representation of the Level Set targeted by the optimizer.
    /// - In order to spare work, each Level Set representation is projected onto a DG Level Set linked to all solver routines, this needs to be done every-time before a boSSS routine is used (e.g. residual evaluation is needed).
    /// </summary>
    public interface IOptiLevelSet : ICloneable {
        /// <summary>
        /// Obtains the Grid the LevelSet is based on
        /// </summary>
        /// <returns></returns>
        IGridData GetGrid();
        /// <summary>
        /// returns Number of Components
        /// </summary>
        /// <returns>Number of Components N</returns>
        int GetLength();
        /// <summary>
        /// sets a parameter $a_i$
        /// </summary>
        /// <param name="index">i param index</param>
        /// <param name="val">new value</param>
        void SetParam(int index, double val);
        /// gets a parameter $a_i$
        /// </summary>
        /// <param name="index">i param index</param>
        /// <returns>param value $a_i$</returns>
        double GetParam(int index);
        /// <summary>
        /// Determines whether the DOF of the given index influences the NearBand of the LevelSet. If it does not the Jacobian entires do not need to be computed.
        /// </summary>
        /// <param name="index">i param index</param>
        /// <returns></returns>
        bool IsInNearBand(int index);
        /// gets a parameter $a_i$
        /// </summary>
        /// <param name="index">i param index</param>
        /// <returns>param value $a_i$</returns>
        double[] GetParamsAsArray();
        /// <summary>
        /// copies params from one OptiLevelSet to another
        /// </summary>
        /// <param name="source">LevelSet to be copied from</param>
        void CopyParamsFrom(IOptiLevelSet source);

        /// <summary>
        /// accumulator for param
        /// 
        /// $a_i \rightarrow a_i+\text{acc}$
        /// </summary>
        /// <param name="index">i param index</param>
        /// <param name="acc">amount to be added</param>
        void AccToParam(int index, double acc);

        /// <summary>
        /// gets the Name of one of the parameters, mostly unneeded.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        string GetParamName(int index);

        /// <summary>
        /// removes oscillations in the zero iso-contour (unused so far)
        /// </summary>
        void Reinitialize(double L, double kappa_S);

        /// <summary>
        /// Transforms this into a SinglePhaseField of some degree
        /// </summary>
        /// <param name="degree">of SInglePhaseField</param>
        /// <returns></returns>
        SinglePhaseField ToSinglePhaseField(int degree);

        /// <summary>
        /// This method assembles a transformation Matrix that maps from the IOptiLevelSetCooridinates into the DGLevelSet coordinates
        /// </summary>
        /// <param name="targetLS"></param>
        void AssembleTransMat(LevelSet targetLS);
        /// <summary>
        /// Method used to project this object onto the DG-LevelSet which is used in the rest of the code
        /// </summary>
        /// <param name="targetLS"></param>
        void ProjectOntoLevelSet(LevelSet targetLS);
        /// <summary>
        /// Prints the Coordinates of this IOptiLevelSet
        /// </summary>
        void Print();
        /// <summary>
        /// This function takes an arbitrary LevelSet and Projects it onto this object using a LevelSet object defined on one cell.
        /// </summary>
        /// <param name="sourceLS">source Level Set</param>
        void ProjectFromForeignLevelSet(SinglePhaseField sourceLS);

        /// <summary>
        /// Projects the member LevelSet onto this object
        /// </summary>
        /// <param name="sourceLS">member LevelSet </param>
        /// <exception cref="ArgumentException"> TransMat was not assembled or assembled with a different LevelSet</exception>
        /// <exception cref="NotImplementedException">only works if optiLevelSet has an orthonormal Basis</exception>
        void ProjectFromLevelSet(ConventionalDGField sourceLS);
        /// <summary>
        /// helper Function to obtain the value of a Function described by the polynomial object 
        /// </summary>
        /// <param name="x">point to be evaluated</param>
        /// <param name="a">coefficient</param>
        /// <param name="p">Polynomial</param>
        /// <returns> function value a * p(x) </returns>
        public static double FuncFromPolynomial(double[] x, double a, Polynomial p) {
            double ret = 0;
            for(int i = 0; i < p.Coeff.Length; i++) {
                ret += p.Coeff[i] * Math.Pow(x[0], p.Exponents[i, 0]) * Math.Pow(x[1], p.Exponents[i, 1]);
            }
            return a * ret;
        }
        /// <summary>
        /// Transforms points living in a variable rectangular domain [xmin,xmax] x [ymin,ymax] into the reference domain [-1,1] x[-1,1]
        /// </summary>
        /// <param name="x"> point to be transformed</param>
        /// <param name="xMin"></param>
        /// <param name="xMax"></param>
        /// <param name="yMin"></param>
        /// <param name="yMax"></param>
        /// <returns>transformed point </returns>
        public static double[] TransformToReferenceElement(double[] x, double xMin, double xMax, double yMin, double yMax) {
            double[] ret = new double[x.Length];
            var scaleX = xMax - xMin;
            var scaleY = yMax - yMin;
            ret[0] = (2.0 * (x[0] - xMin) - scaleX) / scaleX;
            ret[1] = (2.0 * (x[1] - yMin) - scaleY) / scaleY;
            return ret;
        }

        /// <summary>
        /// A test that checks if the LevelSet is Orthonormal by computing the L2-Scalar Product of each pair of components
        /// </summary>
        /// <returns> a bool that is true if orthonormality is given</returns>
        bool TestOrthonormality();
        public static string NameFromPolynomial(Polynomial p) {
            string ret = "";
            for(int i = 0; i < p.Coeff.Length; i++) {
                if(i < p.Coeff.Length - 1) {
                    if(p.Exponents[i, 0] > 0) {
                        if(p.Exponents[i, 1] > 0) {
                            ret += Math.Round(p.Coeff[i], 2) + "x^" + p.Exponents[i, 0] + "y^" + p.Exponents[i, 1] + " + ";
                        } else {
                            ret += Math.Round(p.Coeff[i], 2) + "x^" + p.Exponents[i, 0] + " + ";
                        }

                    } else {
                        if(p.Exponents[i, 1] > 0) {
                            ret += Math.Round(p.Coeff[i], 2) + "y^" + p.Exponents[i, 1] + " + ";
                        } else {
                            ret += Math.Round(p.Coeff[i], 2) + " + ";
                        }
                    }
                } else {
                    if(p.Exponents[i, 0] > 0) {
                        if(p.Exponents[i, 1] > 0) {
                            ret += Math.Round(p.Coeff[i], 2) + "x^" + p.Exponents[i, 0] + "y^" + p.Exponents[i, 1];
                        } else {
                            ret += Math.Round(p.Coeff[i], 2) + "x^" + p.Exponents[i, 0];
                        }

                    } else {
                        if(p.Exponents[i, 1] > 0) {
                            ret += Math.Round(p.Coeff[i], 2) + "y^" + p.Exponents[i, 1];
                        } else {
                            ret += Math.Round(p.Coeff[i], 2);
                        }
                    }
                }


            }


            return ret;
        }
        /// <summary>
        /// This method projects a scalar function onto the OptiLevelSet. It is only usable if the OptiLevelSet has an Orthonormal basis (constructed by the <see cref="CreateONBLevelSet"/> method)
        /// 
        /// The projection is done:
        /// 1. initializing a SinglePhasephield using the member grid 
        /// 2. Then the projection matrix is assembled (projection from the OptiLevelSet field onto the SInglephasefield)
        /// 3. The Transpose of that matrix is multiplied with the coordinate Vector of the Singlephasefield so the coordinateVector for the optiLevelSet is Obtained
        /// </summary>
        /// <param name="initialShockPostion"></param>
        void ProjectFromFunction(Func<double[], double> initialShockPostion);
        /// <summary>
        /// Get the polynomial degree of the level set
        /// </summary>
        /// <returns></returns>
        int GetDegree();
        /// <summary>
        /// get the mean value in a cell
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        double GetMeanValue(int i);
        /// <summary>
        /// Transform to LevelSet object
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        LevelSet ToLevelSet(int v);
        /// <summary>
        /// Assemble a level set Tracker
        /// </summary>
        void AssembleTracker();
        /// <summary>
        /// Calc level-set specific norm
        /// </summary>
        /// <param name="levelSetStepCoordinates"></param>
        /// <returns></returns>
        double Norm(double[] levelSetStepCoordinates);
        /// <summary>
        /// a method that returns the level-set specific regularization Matrix
        /// </summary>
        /// <returns></returns>
        MsrMatrix GetRegMatrix();
    }
}
