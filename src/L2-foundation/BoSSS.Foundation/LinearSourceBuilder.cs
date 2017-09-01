using System;
using System.Collections.Generic;
using System.Text;
using BoSSS.Foundation.LinAlg;

namespace BoSSS.Foundation.Utils {
    
    /*
    /// <summary>
    /// Implementing a <see cref="ILinearSource"/>-object maybe a little bit
    /// cumbersome, so
    /// this class can be used to define a <see cref="ILinearSource"/> from a
    /// <see cref="INonlinearSource"/>-object; 
    /// This class automatically evaluates
    /// a function matrix and an offset vector. 
    /// Of course, the functions defined by the <see cref="INonlinearSource"/>-object 
    /// must affine linear to produce a correct result. This class doesn't tests this issue;
    /// </summary>
    public sealed class LinearSourceBuilder : ILinearSource {


        
        
        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="src">a implementation of a source function that must be 
        /// affin-linear in a mathematical sence;</param>
        /// <param name="a_SourceTimeInDependent"><see cref="ILinearSource.SourceTimeIndependent"/></param>
        /// <param name="a_SourceSpaceInDependent"><see cref="ILinearSource.SourceSpaceInDependent"/></param>
        public LinearSourceBuilder(INonlinearSource src, 
            bool a_SourceTimeInDependent, bool a_SourceSpaceInDependent) {
            p_SourceSpaceInDependent = a_SourceSpaceInDependent;
            p_SourceTimeIndependent = a_SourceTimeInDependent;
            m_Source = src;

            m_Arguments = new double[m_Source.ArgumentOrdering.Length];
        }


        /// <summary>
        /// field containing value of coresponding Property
        /// </summary>
        private bool p_SourceSpaceInDependent = false;

        /// <summary>
        /// field containing value of coresponding Property
        /// </summary>
        private bool p_SourceTimeIndependent;


        /// <summary>
        /// linear martices are constructed from this object.
        /// </summary>
        internal INonlinearSource m_Source;

        #region ILinearSource Member

 


        /// <summary>
        /// <see cref="ILinearFlux.BorderFluxTimeInDependent"/>
        /// can be specifyed in the constructor of this object;
        /// </summary>
        public bool SourceSpaceInDependent {
            get { return p_SourceSpaceInDependent; }
        }


        /// <summary>
        /// <see cref="ILinearFlux.BorderFluxSpaceInDependent"/>
        /// can be specifyed in the constructor of this object;
        /// </summary>
        public bool SourceTimeIndependent {
            get { return p_SourceTimeIndependent; }
        }


        /// <summary>
        /// arguments used for the nonlinear flux functions to extract the function matrix
        /// </summary>
        double[] m_Arguments;


        /// <summary>
        /// This method can be used to test if the construction of the linear affine source
        /// was correct
        /// </summary>
        /// <returns>if the given nonlinear source is representable as an affine linear function,
        /// this should return 0 or a very small value</returns>
        public double Test() {
            Vector2D X;
            Vector2D n;

            int L = m_Source.ArgumentOrdering.Length;
            int I = 20*L;

            Random rand = new Random();

            double[] U = new double[L];
            
            double[] mat = new double[L];
            double o;

            double errsum = 0.0;
            for (int i = 0; i < I; i++) {

                // populate Arguments with random values
                // -------------------------------------
                X.x = rand.NextDouble()*10.0 - 5.0;
                X.y = rand.NextDouble()*10.0 - 5.0;

                for (int l = L-1; l>= 0; l--) {
                    U[l] = rand.NextDouble()*2.0 - 1.0;
                }

                n.x = rand.NextDouble()*2.0 - 1.0;
                n.y = rand.NextDouble()*2.0 - 1.0;
                n.Normalize();

                
                // "nonlinear" version
                double nls = m_Source.Source(0.0,X,U);

                // linear version
                this.Source(0.0, X, mat, out o);
                double lins = o;
                for (int l = L-1; l>= 0; l--) lins += mat[l]*U[l];

                errsum += Math.Abs(nls-lins);
            }

            return errsum;
        }


        /// <summary>
        /// <see cref="ILinearSource.Source"/>
        /// </summary>
        public void Source(double time, Vector2D x, double[] FunctionMatrix, out double AffineOffset) {
            int l = m_Arguments.Length;
            Array.Clear(m_Arguments, 0, l);

            AffineOffset = m_Source.Source(time, x, m_Arguments);
            for (int i = 0; i < l; i++) {
                m_Arguments[i] = 1.0;
                FunctionMatrix[i] = m_Source.Source(time, x, m_Arguments) - AffineOffset;
                m_Arguments[i] = 0.0;
            }
        }

        #endregion

        #region IEquationCommon Member


        /// <summary>
        /// returns the same argument ordering like the given <see cref="INonlinearSource"/>
        /// </summary>
        public Field[] ArgumentOrdering {
            get {
                return m_Source.ArgumentOrdering;
            }
        }

        #endregion
    }
     */
}
