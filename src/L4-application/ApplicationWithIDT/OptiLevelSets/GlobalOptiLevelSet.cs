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
using System.Linq;
using ilPSP.LinSolvers;

namespace ApplicationWithIDT.OptiLevelSets {
    /// <summary>
    /// Global Level Set
    /// - defined by a global function phi(x,y)= a_1 *f_1(x,y) + a_2*f_2(x,y) + .... with some DOFS (a_1,a_2,...) and custom basis functions (f_1,f_2,...)
    /// </summary>
    public class GlobalOptiLevelSet : IOptiLevelSet {

        /// <summary>
        /// custom names of the parameters
        /// </summary>
        public List<string> m_ParamNames {
            get;
            private set;
        }
        /// <summary>
        /// Matrix that transforms this OptiLevelSet to the LevelSet used by BoSSS
        /// </summary>
        public MultidimensionalArray m_TransMat {
            get;
            private set;
        }
        /// <summary>
        /// current coordinates
        /// </summary>
        public List<double> m_ParamValues {
            get;
            private set;
        }
        /// <summary>
        /// Number of basis functions
        /// </summary>
        public int Length {
            get;
            private set;
        }
        /// <summary>
        /// Spatial dimension of the underlying domain
        /// </summary>
        public int dim {
            get;
            private set;
        }
        /// <summary>
        /// the basis functions
        /// </summary>
        public List<Func<double[], double, double>> m_phi {
            get;
            private set;
        }
        /// <summary>
        /// the Grid used in the Method
        /// </summary>
        public IGrid m_grid {
            get;
            private set;
        }

        /// <summary>
        /// the degree of the LevelSet (if polynomial)
        /// </summary>
        private int degree;



        /// <summary>
        ///  a grid with one cell covering the computational domain
        /// </summary>
        private Grid2D globalGrid;

        /// <summary>
        /// Some mapping (used for Agglomeration only)
        /// </summary>
        public CoordinateMapping m_mapping {
            get;
            private set;
        }

        /// <summary>
        /// a hack to compute agglomeration...
        /// </summary>
        public SinglePhaseField[] m_fields {
            get;
            private set;
        }

        /// <summary>
        /// indicating if the basis is Orthonormal
        /// </summary>
        public bool IsOrthonormal = false;

        public GlobalOptiLevelSet(IGrid grid, int dim = 2) {
            this.dim = dim;
            m_ParamNames = new List<string>();
            m_ParamValues = new List<double>();
            m_phi = new List<Func<double[], double, double>>();
            Length = 0;
            m_grid = grid;
            degree = 0;
            var nodes = m_grid.iGridData.iVertices.Coordinates;
            var coord = nodes.ExtractSubArrayShallow(0, -1).To1DArray();
            double yMin = coord[1];
            double yMax = coord[1];
            double xMax = coord[0];
            double xMin = coord[0];
            for(int i = 1; i < nodes.GetLength(0); i++) {
                coord = nodes.ExtractSubArrayShallow(i, -1).To1DArray();
                if(coord[1] < yMin) {
                    yMin = coord[1];
                }
                if(coord[1] > yMax) {
                    yMax = coord[1];
                }
                if(coord[0] < xMin) {
                    xMin = coord[0];
                }
                if(coord[0] > xMax) {
                    xMax = coord[0];
                }
            }
            double[] xNodes = new double[] { xMin, xMax };
            double[] yNodes = new double[] { yMin, yMax };
            globalGrid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
            m_fields = new SinglePhaseField[1];
        }
        public GlobalOptiLevelSet(IGrid grid, List<string> m_ParamNames, List<double> m_ParamValues, List<Func<double[], double, double>> m_phi, int degree, bool isOrthormal = false) {
            if(m_ParamNames.Count != m_ParamValues.Count) {
                throw new ArgumentException("OptiLevelSet: number of Params must equal Number of Param Names");
            } else if(m_ParamNames.Count != m_phi.Count) {
                throw new ArgumentException("OptiLevelSet: number of phiComponents must equal Number of Param Names");
            }
            this.degree = degree;
            this.m_phi = m_phi;
            this.m_ParamNames = m_ParamNames;
            this.m_ParamValues = m_ParamValues;
            IsOrthonormal = isOrthormal;
            m_grid = grid;
            m_fields = new SinglePhaseField[1];
            var nodes = m_grid.iGridData.iVertices.Coordinates;
            var coord = nodes.ExtractSubArrayShallow(0, -1).To1DArray();
            double yMin = coord[1];
            double yMax = coord[1];
            double xMax = coord[0];
            double xMin = coord[0];
            for(int i = 1; i < nodes.GetLength(0); i++) {
                coord = nodes.ExtractSubArrayShallow(i, -1).To1DArray();
                if(coord[1] < yMin) {
                    yMin = coord[1];
                }
                if(coord[1] > yMax) {
                    yMax = coord[1];
                }
                if(coord[0] < xMin) {
                    xMin = coord[0];
                }
                if(coord[0] > xMax) {
                    xMax = coord[0];
                }
            }
            double[] xNodes = new double[] { xMin, xMax };
            double[] yNodes = new double[] { yMin, yMax };
            globalGrid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

            Length = m_ParamNames.Count;
        }

        /// <summary>
        /// as the LevelSet is global this is always true
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public bool IsInNearBand(int index) {
            return true;
        }


        /// <summary>
        /// Creates an Global ONB on a square domain
        /// </summary>
        /// <param name="xMin">smallest x value of Square </param>
        /// <param name="xMax">biggest x value of Square</param> 
        /// <param name="yMin">smallest y value of Square</param>
        /// <param name="yMax">biggest y value of Square</param>
        /// <param name="degree"> polynomial degree of ONB </param>
        public void CreateGlobalONBOnSquareDomain(double xMin, double xMax, double yMin, double yMax, int degree) {
            var grid = Grid2D.Cartesian2DGrid(new double[] { xMin, xMax }, new double[] { yMin, yMax });
            var PolyList = grid.GridData.ChefBasis.GetOrthonormalPolynomials(degree);
            m_ParamNames.Clear();
            m_ParamValues.Clear();
            m_phi.Clear();
            this.degree = degree;
            foreach(Polynomial p in PolyList[0]) {
                m_phi.Add(ConvertPolynomialToFunc(p));
                m_ParamValues.Add(1);
                m_ParamNames.Add("p.ToString()");
            }

        }

        private Func<double[], double, double> ConvertPolynomialToFunc(Polynomial p) {


            throw new NotImplementedException();
        }


        /// <summary>
        /// standard clone method
        /// </summary>
        /// <returns>cloned OptiLevelSet </returns>
        public object Clone() {
            GlobalOptiLevelSet r = new GlobalOptiLevelSet(m_grid, dim);

            // done via for loop to ensure we don't get a shallow copy
            for(int i = 0; i < Length; i++) {
                r.m_ParamNames.Add(m_ParamNames[i]);
                r.m_ParamValues.Add(m_ParamValues[i]);
                r.m_phi.Add(m_phi[i]);
            }
            r.Length = Length;
            r.degree = degree;
            r.IsOrthonormal = IsOrthonormal;
            r.globalGrid = globalGrid;//here
            r.m_TransMat = m_TransMat; // and here we want a shallow copy
            r.m_fields = m_fields;
            return r;
        }

        /// <summary>
        /// guess what?
        /// </summary>
        public GlobalOptiLevelSet CloneAs() {
            return (GlobalOptiLevelSet)Clone();
        }
        /// <summary>
        /// Evaluates the OptilevelSet at a point x:
        /// : $\sum_i a_i \varphi_i(x) $
        /// </summary>
        /// <param name="x"> point in space</param>
        /// <returns> value of LevelSet at point x</returns>
        public double Evaluator(double[] x) {
            double ret = 0;
            for(int i = 0; i < Length; i++) {
                ret += m_phi[i](x, m_ParamValues[i]);
            }
            return ret;
        }
        /// <summary>
        /// Converts Evaluator to a Func object
        /// </summary>
        /// <returns>Evaluator of LevelSet</returns>
        public Func<double[], double> GetEvaluator() {
            Func<double[], double> f = Evaluator;
            return f;
        }
        /// <summary>
        /// Converts Evaluator into a vectorized Version 
        /// </summary>
        /// <returns></returns>
        public ScalarFunction GetVectorizedEvaluator() {
            Func<double[], double> f = Evaluator;
            return delegate (MultidimensionalArray inp, MultidimensionalArray res) {
                int D = inp.GetLength(1);
                double[] X = new double[D];

                for(int i = 0; i < inp.GetLength(0); i++) {
                    for(int d = 0; d < D; d++)
                        X[d] = inp[i, d];

                    res[i] = f(X);
                }
            };
        }
        /// <summary>
        /// adds a new basis function/component $\varphi_{N+1}$ to LevelSet
        /// $$\varphi = \sum^{N}_i a_i \varphi_i \rightarrow \varphi = \sum^{N+1}_i a_i \varphi_i $$
        /// </summary>
        /// <param name="phi_component"> basis function \varphi_{N+1}</param>
        /// <param name="ParamName"> name of the parameter </param>
        /// <param name="ParamValue"> initial value a_{N+1} </param>
        public void AddComponent(Func<double[], double, double> phi_component, string ParamName, double ParamValue = 1) {
            // 2D test
            double test = 0;
            double[] test_x = new double[dim];
            try {
                test = phi_component(test_x, ParamValue);
            } catch {
                throw new IndexOutOfRangeException("OptiLevelSet: The new Component you added operates on dimensions higher then " + dim);
            }
            m_phi.Add(phi_component);
            m_ParamNames.Add(ParamName);
            m_ParamValues.Add(ParamValue);
            Length = m_phi.Count;
        }
        /// <summary>
        /// Erases a Component
        /// </summary>
        /// <param name="index">index of component to be erased</param>
        public void EraseComponent(int index) {
            if(index < 0) {
                throw new ArgumentException("OptiLevelSet: negative Index cannot be removed");
            } else if(index >= Length) {
                throw new ArgumentException("OptiLevelSet: Index larger then number of Parameters");
            }
            m_phi.RemoveAt(index);
            m_ParamNames.RemoveAt(index);
            m_ParamValues.RemoveAt(index);
            Length = m_phi.Count;
        }
        /// <summary>
        /// returns Number of Components
        /// </summary>
        /// <returns>Number of Components N</returns>
        public int GetLength() {
            return Length;
        }
        /// <summary>
        /// sets a parameter $a_i$
        /// </summary>
        /// <param name="index">i param index</param>
        /// <param name="val">new value</param>
        public void SetParam(int index, double val) {
            m_ParamValues[index] = val;
        }
        /// <summary>
        /// gets a parameter $a_i$
        /// </summary>
        /// <param name="index">i param index</param>
        /// <returns>param value $a_i$</returns>
        public double GetParam(int index) {
            return m_ParamValues[index];
        }
        /// <summary>
        /// copies params from one IOptiLevelSet to this
        /// </summary>
        /// <param name="source">LevelSet to be copied from</param>
        public void CopyParamsFrom(IOptiLevelSet source) {
            if(GetLength() == source.GetLength()) {
                for(int i = 0; i < Length; i++) {
                    m_ParamNames[i] = source.GetParamName(i);
                    m_ParamValues[i] = source.GetParam(i);
                }
            } else {
                throw new ArgumentException("OptiLevelSet.CopyParamsFrom: source length, " + source.GetLength() + "  differs from target length: " + Length);
            }
        }

        /// <summary>
        /// accumulator for param
        /// 
        /// $a_i \rightarrow a_i+\text{acc}$
        /// </summary>
        /// <param name="index">i param index</param>
        /// <param name="acc">amount to be added</param>
        public void AccToParam(int index, double acc) {
            m_ParamValues[index] += acc;
        }

        public string GetParamName(int index) {
            return m_ParamNames[index];
        }
        public SinglePhaseField ToSinglePhaseField(int degree) {
            var basis = new Basis(globalGrid, degree);
            var ret = new SinglePhaseField(basis);
            var f = GetVectorizedEvaluator();
            ret.ProjectField(f);
            return ret;
        }
        public LevelSet ToLevelSet(int degree) {
            var basis = new Basis(globalGrid, degree);
            var ret = new LevelSet(basis, "OptiLevelSet");
            var f = GetVectorizedEvaluator();
            ret.ProjectField(f);
            return ret;
        }

        public void UpdateFieldAndMapping() {
            m_fields = new SinglePhaseField[1];
            m_fields[0] = ToSinglePhaseField(degree);
            m_mapping = new CoordinateMapping(m_fields);
        }
        public void AssembleTransMat(LevelSet targetLS) {
            m_TransMat = MultidimensionalArray.Create(targetLS.CoordinateVector.Length, GetLength());
            GlobalOptiLevelSet tmp = CloneAs();
            LevelSet LSbackup = targetLS.CloneAs();
            //Set all params of mp to zero
            for(int i = 0; i < GetLength(); i++) {
                tmp.SetParam(i, 0);
            }
            for(int i = 0; i < GetLength(); i++) {
                tmp.SetParam(i, 1);
                var f = tmp.GetVectorizedEvaluator();
                LSbackup.ProjectField(f);
                for(int j = 0; j < LSbackup.CoordinateVector.Length; j++) {
                    m_TransMat[j, i] = LSbackup.CoordinateVector[j];
                }
                tmp.SetParam(i, 0);
            }
        }
        public void ProjectOntoLevelSet(LevelSet targetLS) {
            if(targetLS.CoordinateVector.Length == m_TransMat.NoOfRows) {
                double[] result = new double[targetLS.CoordinateVector.Length];
                m_TransMat.MatVecMul(1.0, m_ParamValues, 1.0, result);
                targetLS.CoordinateVector.Clear();
                targetLS.CoordinateVector.Acc(1.0, result);
            } else {
                throw new ArgumentException("TransMat was not assembled or assembled with a different LevelSet");
            }

        }
        public void Print() {
            for(int i = 0; i < m_ParamNames.Count - 1; i++) {
                Console.Write(m_ParamNames[i] + "=" + m_ParamValues[i] + " ,");
            }
            Console.WriteLine(m_ParamNames[m_ParamNames.Count] + "=" + m_ParamValues[m_ParamNames.Count]);
        }
        /// <summary>
        /// This function creates an global OptiLevelSet with orthonormal Basis functions on a rectangular domain with variable degree.
        /// This is achieved by extracting polynomial basis functions from Polynomial lists.
        /// </summary>
        /// <param name="xMin"></param>
        /// <param name="xMax"></param>
        /// <param name="yMin"></param>
        /// <param name="yMax"></param>
        /// <param name="degree"></param>
        /// <returns></returns>
        public static GlobalOptiLevelSet CreateONBLevelSet(double xMin, double xMax, double yMin, double yMax, int degree) {
            var LSO = new GlobalOptiLevelSet(Grid2D.Cartesian2DGrid(new double[] { xMin, xMax }, new double[] { yMin, yMax }));
            var grid = Grid2D.Cartesian2DGrid(new double[] { xMin, xMax }, new double[] { yMin, yMax });
            var pList = grid.GridData.ChefBasis.GetOrthonormalPolynomials(degree);
            var tmp_fld = new SinglePhaseField(new Basis(grid.GridData, degree));
            for(int i = 0; i < tmp_fld.Basis.Polynomials[0].Count; i++) {
                var p = tmp_fld.Basis.Polynomials[0][i];
                Func<double[], double> func = x => FuncFromPolynomial(TransformToReferenceElement(x, xMin, xMax, yMin, yMax), 1.0, p);
                //tmp_fld.ProjectFunction(1.0,(x) => FuncFromPolynomial(TransformToReferenceElement(x, xMin, xMax, yMin, yMax), 1.0, p), new CellQuadratureScheme(),null);
                tmp_fld.ProjectField(x => FuncFromPolynomial(TransformToReferenceElement(x, xMin, xMax, yMin, yMax), 1.0, p));
                var OneOverVol = 1 / tmp_fld.L2Norm();
                //Console.WriteLine(OneOverVol);
                Func<double[], double, double> polyFunc = (x, a) => OneOverVol * FuncFromPolynomial(TransformToReferenceElement(x, xMin, xMax, yMin, yMax), a, p);

                LSO.AddComponent(polyFunc, NameFromPolynomial(p), 1.0);
            }
            LSO.IsOrthonormal = true;
            LSO.m_grid = grid;
            LSO.degree = degree;
            return LSO;
        }
        /// <summary>
        /// This function takes an arbitrary LevelSet and Projects it onto this object using a LevelSet object defined on one cell.
        /// </summary>
        /// <param name="sourceLS">source Level Set</param>
        public void ProjectFromForeignLevelSet(SinglePhaseField sourceLS) {

            var tmpLS = new LevelSet(new Basis(globalGrid, degree), "tmpLS"); //m_grid should be made out of one cell
            if(globalGrid.NumberOfCells > 1) {
                throw new NotSupportedException("wrong grid used - should have one cell only");
            }
            var TransMat = MultidimensionalArray.Create(tmpLS.CoordinateVector.Length, GetLength());
            GlobalOptiLevelSet tmp = CloneAs();
            LevelSet LSbackup = tmpLS.CloneAs();
            //Set all params of mp to zero
            for(int i = 0; i < GetLength(); i++) {
                tmp.SetParam(i, 0);
            }
            for(int i = 0; i < GetLength(); i++) {
                tmp.SetParam(i, 1);
                var f = tmp.GetVectorizedEvaluator();
                LSbackup.ProjectField(f);
                for(int j = 0; j < LSbackup.CoordinateVector.Length; j++) {
                    TransMat[j, i] = LSbackup.CoordinateVector[j];
                }
                tmp.SetParam(i, 0);
            }
            var TransMatTransposed = TransMat.TransposeTo();

            tmpLS.Clear();
            //here we project the sourceLS onto the LevelSet with only one cell
            tmpLS.ProjectFromForeignGrid(1.0, sourceLS);
            //Then we project onto our OptiLevelSet
            double[] result = new double[GetLength()];
            TransMatTransposed.MatVecMul(1.0, tmpLS.CoordinateVector, 1.0, result);
            m_ParamValues.Clear();
            m_ParamValues.AddRange(result);
        }
        /// <summary>
        /// Projects the member LevelSet onto this object
        /// </summary>
        /// <param name="sourceLS">member LevelSet </param>
        /// <exception cref="ArgumentException"> TransMat was not assembled or assembled with a different LevelSet</exception>
        /// <exception cref="NotImplementedException">only works if optiLevelSet has an orthonormal Basis</exception>
        public void ProjectFromLevelSet(ConventionalDGField sourceLS) {
            if(IsOrthonormal == true) {
                if(sourceLS.CoordinateVector.Length == m_TransMat.NoOfRows) {
                    double[] result = new double[GetLength()];
                    var TransMatTransposed = m_TransMat.TransposeTo();
                    TransMatTransposed.MatVecMul(1.0, sourceLS.CoordinateVector, 1.0, result);
                    m_ParamValues.Clear();
                    m_ParamValues.AddRange(result);
                } else {
                    throw new ArgumentException("TransMat was not assembled or assembled with a different LevelSet");
                }
            } else {
                throw new NotImplementedException("Reverse Transformation Only Implemented for Orthonormal LevelSet");
            }
        }
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
        /// this helper function gives an evaluator where only one component is active (this is only used to check if the basis is truly orthonormal)
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public Func<double[], double> GetOneComponent(int index) {
            double[] copy = new double[m_ParamValues.Count];
            m_ParamValues.CopyTo(copy);
            for(int i = 0; i < GetLength(); i++) {
                SetParam(i, 0);
            }
            SetParam(index, 1.0);
            var ev = GetEvaluator().CloneAs();
            return ev;
        }

        /// <summary>
        /// A test that checks if the LevelSet is Orthonormal by computing the L2-Scalar Product of each pair of components
        /// </summary>
        /// <returns> a bool that is true if orthonormality is given</returns>
        public bool TestOrthonormality() {
            bool isOrtho = false;
            var LSO_Clone = CloneAs();
            var field = new LevelSet(new Basis(m_grid, degree), "testing_field");
            double tol = 1e-14;
            for(int i = 0; i < GetLength(); i++) {
                var ev = GetOneComponent(i);
                for(int j = 0; j < GetLength(); j++) {
                    var ev2 = LSO_Clone.GetOneComponent(j);
                    field.ProjectField(x => ev2(x) * ev(x));
                    double res = field.IntegralOver(new CellMask(m_grid.iGridData, new BitArray(new bool[] { true })));
                    if(i == j && Math.Abs(res - 1) < tol) {
                        isOrtho = true;
                    } else if(i == j && Math.Abs(res - 1) > tol) {
                        isOrtho = false;
                        Console.WriteLine("OptiLevelSet.TestOrthonormality(): Component " + i + " is not normalized!");
                        break;
                    } else if(i != j && Math.Abs(res - 0) < tol) {
                        isOrtho = true;
                    } else {
                        isOrtho = false;
                        Console.WriteLine("OptiLevelSet.TestOrthonormality(): Components " + i + ", " + j + " are not normalized!");
                        break;
                    }
                }
            }
            return isOrtho;
        }
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
        public void ProjectFromFunction(Func<double[], double> initialShockPostion) {
            var tmp_fld = new SinglePhaseField(new Basis(m_grid, degree));
            if(IsOrthonormal == true) {
                var TransMat = MultidimensionalArray.Create(tmp_fld.CoordinateVector.Length, GetLength());
                GlobalOptiLevelSet tmp = CloneAs();
                var LSbackup = tmp_fld.CloneAs();
                //Set all params of mp to zero
                for(int i = 0; i < GetLength(); i++) {
                    tmp.SetParam(i, 0);
                }
                for(int i = 0; i < GetLength(); i++) {
                    tmp.SetParam(i, 1);
                    var f = tmp.GetVectorizedEvaluator();
                    LSbackup.ProjectField(f);
                    for(int j = 0; j < LSbackup.CoordinateVector.Length; j++) {
                        TransMat[j, i] = LSbackup.CoordinateVector[j];
                    }
                    tmp.SetParam(i, 0);
                }
                tmp_fld.ProjectField(initialShockPostion);
                double[] result = new double[GetLength()];
                var TransMatTransposed = TransMat.TransposeTo();
                TransMatTransposed.MatVecMul(1.0, tmp_fld.CoordinateVector, 1.0, result);
                m_ParamValues.Clear();
                m_ParamValues.AddRange(result);
            } else {
                throw new ArgumentException("LevelSet not Orthonormal");
            }

        }

        public double[] GetParamsAsArray() {
            return m_ParamValues.ToArray();
        }

        public IGridData GetGrid() {
            return m_grid.iGridData;
        }

        public int GetDegree() {
            return degree;
        }

        public double GetMeanValue(int i) {
            return ToSinglePhaseField(degree).GetMeanValue(i);
        }

        /// <summary>
        /// does Nothing as a Tracker is not needed here
        /// </summary>
        public void AssembleTracker() {
        }

        public double Norm(double[] levelSetStepCoordinates) {
            //Old version of calculation of norm (works better for ECCOMAS Test Case Scalar Advection)
            //double norm = 0;
            //for(int i = 0; i < this.GetLength(); i++) {
            //    norm += levelSetStepCoordinates[i];
            //}
            return levelSetStepCoordinates.Sum().Sqrt();
        }

        /// <summary>
        /// here nothing happens, as it is not clear how so far
        /// </summary>
        public void Reinitialize(double L, double kappa_S) {

        }
        public MsrMatrix GetRegMatrix() {
            var length_l = this.GetLength();
            var reg = new MsrMatrix(length_l, length_l, 1, 1);
            IMutuableMatrixEx_Extensions.AccEyeSp(reg, 1.0);
            return reg;
        }
    }
}
