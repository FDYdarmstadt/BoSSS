using BoSSS.Foundation;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.SpecificSolutions {


    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class RotatingDiskBoundaryLayerConditionsHAM_VelocityX : IBoundaryAndInitialData {

        /// <summary>
        /// rotation rate
        /// </summary>
        [JsonProperty]
        double m_OmegaStar;

        [JsonProperty]
        double m_kinViscosity;

        [JsonProperty]
        double m_density;


        /// <summary>
        /// HAM coefficeint for radial velocity U
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableU;

        /// <summary>
        ///  HAM coefficeint for azimuthal velocity V
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableV;


        public void SetData(MultidimensionalArray _CoeffVariableU, MultidimensionalArray _CoeffVariableV,
            double OmegaStar, double kinViscosity, double density) {

            m_OmegaStar = OmegaStar;
            m_kinViscosity = kinViscosity;
            m_density = density;

            if (_CoeffVariableU.Dimension != 2)
                throw new ArgumentException("_CoeffVariableU: Expecting Dimension 2");
            if (_CoeffVariableU.Lengths[0] != _CoeffVariableU.Lengths[1])
                throw new ArgumentException("_CoeffVariableU: Expecting same lengths for both dimensions");

            if (_CoeffVariableV.Dimension != 2)
                throw new ArgumentException("_CoeffVariableV: Expecting Dimension 2");
            if (_CoeffVariableV.Lengths[0] != _CoeffVariableV.Lengths[1])
                throw new ArgumentException("_CoeffVariableV: Expecting same lengths for both dimensions");

            if (_CoeffVariableU.Lengths[0] != _CoeffVariableV.Lengths[0])
                throw new ArgumentException("Expecting same lengths for both coeff inputs");


            this.CoeffVariableU = _CoeffVariableU;
            this.CoeffVariableV = _CoeffVariableV;

        }


        /// <summary>
        /// Evaluation
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public double Evaluate(double[] X, double t) {

            if (X.Length != 3)
                throw new NotSupportedException("only supported for 3D");

            double m_Lstar = Math.Sqrt(m_kinViscosity / m_OmegaStar);

            double z = X[2];
            double eta = z / m_Lstar;

            int coeffDim = CoeffVariableU.Lengths[0];

            double velU = 0.0; double velV = 0.0;
            for (int n = 0; n < coeffDim; n++) {
                for (int i = 0; i < coeffDim; i++) {
                    velU += CoeffVariableU[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                    velV += CoeffVariableV[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                }
            }

            double rStarOnDisk = Math.Sqrt(X[0].Pow2() + X[1].Pow2());
            double theta = Math.Atan2(X[1], X[0]);

            double returnVal = velU * Math.Cos(theta) - velV * Math.Sin(theta);
            returnVal *= rStarOnDisk * m_OmegaStar;

            return returnVal;
        }


        /// <summary>
        /// Vectorized evaluation
        /// </summary>
        /// <param name="input"></param>
        /// <param name="time"></param>
        /// <param name="output"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }



        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as RotatingDiskBoundaryLayerConditionsHAM_VelocityX;
            if (other == null)
                return false;


            if (m_density != other.m_density)
                return false;
            if (m_kinViscosity != other.m_kinViscosity)
                return false;
            if (m_OmegaStar != other.m_OmegaStar)
                return false;

            if (!CoeffVariableU.Equals(other.CoeffVariableU))
                return false;
            if (!CoeffVariableV.Equals(other.CoeffVariableV))
                return false;


            return true;
        }
    }


    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class RotatingDiskBoundaryLayerConditionsHAM_VelocityU : IBoundaryAndInitialData {

        /// <summary>
        /// rotation rate
        /// </summary>
        [JsonProperty]
        double m_OmegaStar;

        [JsonProperty]
        double m_kinViscosity;

        [JsonProperty]
        double m_density;


        /// <summary>
        /// HAM coefficeint for radial velocity U
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableU;


        public void SetData(MultidimensionalArray _CoeffVariableU, 
            double OmegaStar, double kinViscosity, double density) {

            m_OmegaStar = OmegaStar;
            m_kinViscosity = kinViscosity;
            m_density = density;

            if (_CoeffVariableU.Dimension != 2)
                throw new ArgumentException("_CoeffVariableU: Expecting Dimension 2");
            if (_CoeffVariableU.Lengths[0] != _CoeffVariableU.Lengths[1])
                throw new ArgumentException("_CoeffVariableU: Expecting same lengths for both dimensions");


            this.CoeffVariableU = _CoeffVariableU;

        }


        /// <summary>
        /// Evaluation
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public double Evaluate(double[] X, double t) {

            if (X.Length != 3)
                throw new NotSupportedException("only supported for 3D");

            double m_Lstar = Math.Sqrt(m_kinViscosity / m_OmegaStar);

            double z = X[2];
            double eta = z / m_Lstar;

            int coeffDim = CoeffVariableU.Lengths[0];

            double velU = 0.0; 
            for (int n = 0; n < coeffDim; n++) {
                for (int i = 0; i < coeffDim; i++) {
                    velU += CoeffVariableU[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                }
            }

            double rStarOnDisk = Math.Sqrt(X[0].Pow2() + X[1].Pow2());

            double returnVal = rStarOnDisk * m_OmegaStar * velU;

            return returnVal;
        }


        /// <summary>
        /// Vectorized evaluation
        /// </summary>
        /// <param name="input"></param>
        /// <param name="time"></param>
        /// <param name="output"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }



        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as RotatingDiskBoundaryLayerConditionsHAM_VelocityU;
            if (other == null)
                return false;


            if (m_density != other.m_density)
                return false;
            if (m_kinViscosity != other.m_kinViscosity)
                return false;
            if (m_OmegaStar != other.m_OmegaStar)
                return false;

            if (!CoeffVariableU.Equals(other.CoeffVariableU))
                return false;


            return true;
        }
    }



    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class RotatingDiskBoundaryLayerConditionsHAM_VelocityY : IBoundaryAndInitialData {

        /// <summary>
        /// rotation rate
        /// </summary>
        [JsonProperty]
        double m_OmegaStar;

        [JsonProperty]
        double m_kinViscosity;

        [JsonProperty]
        double m_density;


        /// <summary>
        /// HAM coefficeint for radial velocity U
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableU;

        /// <summary>
        ///  HAM coefficeint for azimuthal velocity V
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableV;


        public void SetData(MultidimensionalArray _CoeffVariableU, MultidimensionalArray _CoeffVariableV,
            double OmegaStar, double kinViscosity, double density) {

            m_OmegaStar = OmegaStar;
            m_kinViscosity = kinViscosity;
            m_density = density;

            if (_CoeffVariableU.Dimension != 2)
                throw new ArgumentException("_CoeffVariableU: Expecting Dimension 2");
            if (_CoeffVariableU.Lengths[0] != _CoeffVariableU.Lengths[1])
                throw new ArgumentException("_CoeffVariableU: Expecting same lengths for both dimensions");

            if (_CoeffVariableV.Dimension != 2)
                throw new ArgumentException("_CoeffVariableV: Expecting Dimension 2");
            if (_CoeffVariableV.Lengths[0] != _CoeffVariableV.Lengths[1])
                throw new ArgumentException("_CoeffVariableV: Expecting same lengths for both dimensions");

            if (_CoeffVariableU.Lengths[0] != _CoeffVariableV.Lengths[0])
                throw new ArgumentException("Expecting same lengths for both coeff inputs");


            this.CoeffVariableU = _CoeffVariableU;
            this.CoeffVariableV = _CoeffVariableV;

        }


        /// <summary>
        /// Evaluation
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public double Evaluate(double[] X, double t) {

            if (X.Length != 3)
                throw new NotSupportedException("only supported for 3D");

            double m_Lstar = Math.Sqrt(m_kinViscosity / m_OmegaStar);

            double z = X[2];
            double eta = z / m_Lstar;

            int coeffDim = CoeffVariableU.Lengths[0];

            double velU = 0.0; double velV = 0.0;
            for (int n = 0; n < coeffDim; n++) {
                for (int i = 0; i < coeffDim; i++) {
                    velU += CoeffVariableU[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                    velV += CoeffVariableV[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                }
            }

            double rStarOnDisk = Math.Sqrt(X[0].Pow2() + X[1].Pow2());
            double theta = Math.Atan2(X[1], X[0]);

            double returnVal = velV * Math.Cos(theta) + velU * Math.Sin(theta);
            returnVal *= rStarOnDisk * m_OmegaStar;

            return returnVal;
        }


        /// <summary>
        /// Vectorized evaluation
        /// </summary>
        /// <param name="input"></param>
        /// <param name="time"></param>
        /// <param name="output"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }



        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as RotatingDiskBoundaryLayerConditionsHAM_VelocityY;
            if (other == null)
                return false;


            if (m_density != other.m_density)
                return false;
            if (m_kinViscosity != other.m_kinViscosity)
                return false;
            if (m_OmegaStar != other.m_OmegaStar)
                return false;

            if (!CoeffVariableU.Equals(other.CoeffVariableU))
                return false;
            if (!CoeffVariableV.Equals(other.CoeffVariableV))
                return false;


            return true;
        }
    }


    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class RotatingDiskBoundaryLayerConditionsHAM_VelocityV : IBoundaryAndInitialData {

        /// <summary>
        /// rotation rate
        /// </summary>
        [JsonProperty]
        double m_OmegaStar;

        [JsonProperty]
        double m_kinViscosity;

        [JsonProperty]
        double m_density;


        /// <summary>
        ///  HAM coefficeint for azimuthal velocity V
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableV;


        public void SetData(MultidimensionalArray _CoeffVariableV,
            double OmegaStar, double kinViscosity, double density) {

            m_OmegaStar = OmegaStar;
            m_kinViscosity = kinViscosity;
            m_density = density;


            if (_CoeffVariableV.Dimension != 2)
                throw new ArgumentException("_CoeffVariableV: Expecting Dimension 2");
            if (_CoeffVariableV.Lengths[0] != _CoeffVariableV.Lengths[1])
                throw new ArgumentException("_CoeffVariableV: Expecting same lengths for both dimensions");


            this.CoeffVariableV = _CoeffVariableV;

        }


        /// <summary>
        /// Evaluation
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public double Evaluate(double[] X, double t) {

            if (X.Length != 3)
                throw new NotSupportedException("only supported for 3D");

            double m_Lstar = Math.Sqrt(m_kinViscosity / m_OmegaStar);

            double z = X[2];
            double eta = z / m_Lstar;

            int coeffDim = CoeffVariableV.Lengths[0];

            double velV = 0.0;
            for (int n = 0; n < coeffDim; n++) {
                for (int i = 0; i < coeffDim; i++) {
                    velV += CoeffVariableV[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                }
            }

            double rStarOnDisk = Math.Sqrt(X[0].Pow2() + X[1].Pow2());

            double returnVal = rStarOnDisk * m_OmegaStar * velV;

            return returnVal;
        }


        /// <summary>
        /// Vectorized evaluation
        /// </summary>
        /// <param name="input"></param>
        /// <param name="time"></param>
        /// <param name="output"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }



        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as RotatingDiskBoundaryLayerConditionsHAM_VelocityV;
            if (other == null)
                return false;


            if (m_density != other.m_density)
                return false;
            if (m_kinViscosity != other.m_kinViscosity)
                return false;
            if (m_OmegaStar != other.m_OmegaStar)
                return false;

            if (!CoeffVariableV.Equals(other.CoeffVariableV))
                return false;


            return true;
        }
    }



    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class RotatingDiskBoundaryLayerConditionsHAM_VelocityZ : IBoundaryAndInitialData {

        /// <summary>
        /// rotation rate
        /// </summary>
        [JsonProperty]
        double m_OmegaStar;

        [JsonProperty]
        double m_kinViscosity;

        [JsonProperty]
        double m_density;


        /// <summary>
        /// HAM coefficeint for radial velocity W
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableW;



        public void SetData(MultidimensionalArray _CoeffVariableW, 
            double OmegaStar, double kinViscosity, double density) {

            m_OmegaStar = OmegaStar;
            m_kinViscosity = kinViscosity;
            m_density = density;

            if (_CoeffVariableW.Dimension != 2)
                throw new ArgumentException("_CoeffVariableW: Expecting Dimension 2");
            if (_CoeffVariableW.Lengths[0] != _CoeffVariableW.Lengths[1])
                throw new ArgumentException("_CoeffVariableW: Expecting same lengths for both dimensions");


            this.CoeffVariableW = _CoeffVariableW;

        }


        /// <summary>
        /// Evaluation
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public double Evaluate(double[] X, double t) {

            if (X.Length != 3)
                throw new NotSupportedException("only supported for 3D");

            double m_Lstar = Math.Sqrt(m_kinViscosity / m_OmegaStar);

            double z = X[2];
            double eta = z / m_Lstar;

            int coeffDim = CoeffVariableW.Lengths[0];

            double velW = 0.0;
            for (int n = 0; n < coeffDim; n++) {
                for (int i = 0; i < coeffDim; i++) {
                    velW += CoeffVariableW[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                }
            }

            double returnVal = velW * Math.Sqrt(m_kinViscosity * m_OmegaStar);

            return returnVal;
        }


        /// <summary>
        /// Vectorized evaluation
        /// </summary>
        /// <param name="input"></param>
        /// <param name="time"></param>
        /// <param name="output"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }



        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as RotatingDiskBoundaryLayerConditionsHAM_VelocityZ;
            if (other == null)
                return false;


            if (m_density != other.m_density)
                return false;
            if (m_kinViscosity != other.m_kinViscosity)
                return false;
            if (m_OmegaStar != other.m_OmegaStar)
                return false;

            if (!CoeffVariableW.Equals(other.CoeffVariableW))
                return false;


            return true;
        }
    }

    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class RotatingDiskBoundaryLayerConditionsHAM_PressureP : IBoundaryAndInitialData {

        /// <summary>
        /// rotation rate
        /// </summary>
        [JsonProperty]
        double m_OmegaStar;

        [JsonProperty]
        double m_kinViscosity;

        [JsonProperty]
        double m_density;


        /// <summary>
        /// HAM coefficeint for radial pressure P
        /// </summary>
        [JsonProperty]
        MultidimensionalArray CoeffVariableP;



        public void SetData(MultidimensionalArray _CoeffVariableP,
            double OmegaStar, double kinViscosity, double density) {

            m_OmegaStar = OmegaStar;
            m_kinViscosity = kinViscosity;
            m_density = density;

            if (_CoeffVariableP.Dimension != 2)
                throw new ArgumentException("_CoeffVariableP: Expecting Dimension 2");
            if (_CoeffVariableP.Lengths[0] != _CoeffVariableP.Lengths[1])
                throw new ArgumentException("_CoeffVariableP: Expecting same lengths for both dimensions");


            this.CoeffVariableP = _CoeffVariableP;


        }


        /// <summary>
        /// Evaluation
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public double Evaluate(double[] X, double t) {

            if (X.Length != 3)
                throw new NotSupportedException("only supported for 3D");

            double m_Lstar = Math.Sqrt(m_kinViscosity / m_OmegaStar);

            double z = X[2];
            double eta = z / m_Lstar;

            int coeffDim = CoeffVariableP.Lengths[0];

            double valP = 0.0; 
            for (int n = 0; n < coeffDim; n++) {
                for (int i = 0; i < coeffDim; i++) {
                    valP += CoeffVariableP[n, i] * Math.Exp(-n * eta) * Math.Pow(eta, i);
                }
            }

            double returnVal = valP * m_density * m_kinViscosity * m_OmegaStar;

            return returnVal;
        }


        /// <summary>
        /// Vectorized evaluation
        /// </summary>
        /// <param name="input"></param>
        /// <param name="time"></param>
        /// <param name="output"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }



        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as RotatingDiskBoundaryLayerConditionsHAM_PressureP;
            if (other == null)
                return false;


            if (m_density != other.m_density)
                return false;
            if (m_kinViscosity != other.m_kinViscosity)
                return false;
            if (m_OmegaStar != other.m_OmegaStar)
                return false;

            if (!CoeffVariableP.Equals(other.CoeffVariableP))
                return false;


            return true;
        }
    }


    /// <summary>
    /// evaluated flow field according to cartesian coordiante system
    /// </summary>
    [Serializable]
    public enum evalFlowFields {

        velocityX = 0,

        velocityY = 1,

        velocityZ = 2,

        pressureP = 3
    }


    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class RotatingDiskBoundaryLayerConditions : IBoundaryAndInitialData {

        /// <summary>
        /// rotation rate
        /// </summary>
        [JsonProperty]
        double m_OmegaStar;

        [JsonProperty]
        double m_kinViscosity;

        [JsonProperty]
        double m_density;

        /// <summary>
        /// Lstar = \sqrt(\nu / \OmegaStar)
        /// </summary>
        [JsonProperty]
        double m_Lstar;


        /// <summary>
        /// dimensionless similarity coordinate
        /// </summary>
        [JsonProperty]
        double[] zValues;


        /// <summary>
        /// similarity variable
        /// </summary>
        [JsonIgnore]
        public evalFlowFields selectedEvalFlowField {
            get {
                return m_evalFlowField;
            }
            set {
                m_evalFlowField = value;
            }
        }

        [JsonProperty]
        evalFlowFields m_evalFlowField;


        /// <summary>
        /// Similarity variable U
        /// </summary>
        [JsonProperty]
        double[] SimilarityVariableU;

        /// <summary>
        /// Similarity variable U
        /// </summary>
        [JsonProperty]
        double[] SimilarityVariableV;

        /// <summary>
        /// Similarity variable U
        /// </summary>
        [JsonProperty]
        double[] SimilarityVariableW;

        /// <summary>
        /// Similarity variable U
        /// </summary>
        [JsonProperty]
        double[] SimilarityVariableP;


        //[JsonProperty]
        //int Udirection;
        //[JsonProperty]
        //int Vdirection;
        //[JsonProperty]
        //int Wdirection;


        /// <summary>
        /// operating point, i.e. center of local coordinate system 
        /// </summary>
        //[JsonProperty]
        //double m_radiusOP;

        /// <summary>
        /// switches between a rotating disk and plain moving wall 
        /// </summary>
        //[JsonProperty]
        //bool m_isRotating;


        public void SetData(double[] z, double[] simVarU, double[] simVarV, double[] simVarW, double[] simVarP, 
            double OmegaStar, double kinViscosity, double density) {

            //m_radiusOP = radiusOP;
            //m_isRotating = rotating;

            //Udirection = Udir;
            //Vdirection = Vdir;
            //Wdirection = Wdir;

            m_OmegaStar = OmegaStar;
            m_kinViscosity = kinViscosity;
            m_density = density;
            m_Lstar = Math.Sqrt(m_kinViscosity / m_OmegaStar);

            int L = z.Length;
            for (int i = 1; i < L; i++) {
                if (z[i] < z[i - 1])
                    throw new ArgumentException("Expecting z input array values to be strictly increasing");
            }

            zValues = z;

            if (simVarU.Length != L)
                throw new ArgumentException("Expecting all input arrays to be of equal length");
            if (simVarV.Length != L)
                throw new ArgumentException("Expecting all input arrays to be of equal length");
            if (simVarW.Length != L)
                throw new ArgumentException("Expecting all input arrays to be of equal length");
            if (simVarP.Length != L)
                throw new ArgumentException("Expecting all input arrays to be of equal length");


            SimilarityVariableU = simVarU;
            SimilarityVariableV = simVarV;
            SimilarityVariableW = simVarW;
            SimilarityVariableP = simVarP;

        }


        /// <summary>
        /// Evaluation
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public double Evaluate(double[] Xstar, double t) {

            if (Xstar.Length != 3)
                throw new NotSupportedException();

            double zStar = Xstar[2];
            double zVal = zStar / m_Lstar;

            (int i1, int i2) = binarySearch(zVal, zValues);

            //double svVal = Interpolate(zVal, zValues[i1], zValues[i2], evalSimilarityVariable[i1], evalSimilarityVariable[i2]);

            double returnVal = 0.0;
            switch (selectedEvalFlowField) {
                case evalFlowFields.velocityX: {
                        double svValU = Interpolate(zVal, zValues[i1], zValues[i2], SimilarityVariableU[i1], SimilarityVariableU[i2]);
                        double svValV = Interpolate(zVal, zValues[i1], zValues[i2], SimilarityVariableV[i1], SimilarityVariableV[i2]);

                        //double[] XonDisk = new double[] { Xstar[0] + m_radiusOP, Xstar[1] };
                        double rStarOnDisk = Math.Sqrt(Xstar[0].Pow2() + Xstar[1].Pow2()); 
                        double theta = Math.Atan2(Xstar[1], Xstar[0]);

                        returnVal = svValU * Math.Cos(theta) - svValV * Math.Sin(theta);
                        returnVal *= rStarOnDisk * m_OmegaStar;
                        break;
                    }
                case evalFlowFields.velocityY: {
                        double svValU = Interpolate(zVal, zValues[i1], zValues[i2], SimilarityVariableU[i1], SimilarityVariableU[i2]);
                        double svValV = Interpolate(zVal, zValues[i1], zValues[i2], SimilarityVariableV[i1], SimilarityVariableV[i2]);

                        //double[] XonDisk = new double[] { Xstar[0] + m_radiusOP, Xstar[1] };
                        double rStarOnDisk = Math.Sqrt(Xstar[0].Pow2() + Xstar[1].Pow2());
                        double theta = Math.Atan2(Xstar[1], Xstar[0]);

                        returnVal = svValV * Math.Cos(theta) + svValU * Math.Sin(theta);
                        returnVal *= rStarOnDisk * m_OmegaStar;
                        break;
                    }
                case evalFlowFields.velocityZ: {
                        double svValW = Interpolate(zVal, zValues[i1], zValues[i2], SimilarityVariableW[i1], SimilarityVariableW[i2]);
                        returnVal = svValW * Math.Sqrt(m_kinViscosity * m_OmegaStar);
                        break;
                    }
                case evalFlowFields.pressureP: {
                        double svValP = Interpolate(zVal, zValues[i1], zValues[i2], SimilarityVariableP[i1], SimilarityVariableP[i2]);
                        returnVal = svValP * m_density * m_kinViscosity * m_OmegaStar;
                        break;
                    }
            }

            return returnVal;
        }

        /// <summary>
        /// Vectorized evaluation
        /// </summary>
        /// <param name="input"></param>
        /// <param name="time"></param>
        /// <param name="output"></param>
        /// <exception cref="NotImplementedException"></exception>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }


        double Interpolate(double zVal, double z1, double z2, double v1, double v2) {
            if (Math.Abs(z2 - z1) <= Math.Abs(v2 - v1) * 1e-10)
                return 0.5 * (v2 + v1);

            double k = (v2 - v1) / (z2 - z1);
            double y = k * (zVal - z1) + v1;

            return y;
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
        static private (int i1, int i2) binarySearch(double x, double[] nodes) {
            if (x <= nodes[0]) {
                return (0, 0);
            } else if (x >= nodes[nodes.Length - 1]) {
                return (nodes.Length - 1, nodes.Length - 1);
            } else {

                int idx = binarySearch_impl(x, nodes);
                if (x > nodes[idx + 1])
                    throw new ApplicationException("overshoot");
                if (x < nodes[idx])
                    throw new ApplicationException("undershoot");
                return (idx, idx + 1);
            }
        }

        static private int binarySearch_impl(double x, double[] nodes) {
            //if(binS_mem >= 0 && this.nodes[binS_mem + 1] >= x && this.nodes[binS_mem] < x)
            //    return binS_mem;

            int min = 0;
            int max = nodes.Length - 1;
            int mid = 0;
            do {
                mid = (min + max) / 2;
                if (x > nodes[mid]) {
                    min = mid + 1;
                } else {
                    max = mid - 1;
                }
            }
            while (mid > 0 && min <= max && (!(nodes[mid] >= x && nodes[mid - 1] < x)) && (!(nodes[mid] < x && nodes[mid + 1] >= x)));
            if (mid > 0 && nodes[mid] >= x && nodes[mid - 1] < x) {
                //binS_mem = mid - 1;
                return (mid - 1);
            }
            //binS_mem = mid;
            return mid;
        }


        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as RotatingDiskBoundaryLayerConditions;
            if (other == null)
                return false;

            if (!SimilarityVariableU.ListEquals(other.SimilarityVariableU))
                return false;
            if (!SimilarityVariableV.ListEquals(other.SimilarityVariableV))
                return false;
            if (!SimilarityVariableW.ListEquals(other.SimilarityVariableW))
                return false;
            if (!SimilarityVariableP.ListEquals(other.SimilarityVariableP))
                return false;


            if (m_density != other.m_density)
                return false;
            if (m_kinViscosity != other.m_kinViscosity)
                return false;
            if (m_Lstar != other.m_Lstar)
                return false;
            if (m_OmegaStar != other.m_OmegaStar)
                return false;
            if (m_evalFlowField != other.m_evalFlowField)
                return false;

            //if (Udirection != other.Udirection)
            //    return false;
            //if (Vdirection != other.Vdirection)
            //    return false;
            //if (Wdirection != other.Wdirection)
            //    return false;

            if (!zValues.ListEquals(other.zValues))
                return false;


            return true;
        }
    }
}
