using System;
using System.Collections.Generic;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Utils {


    
    /// <summary>
    /// Implementing a <see cref="ILinearFlux"/>-object maybe a little bit
    /// cumbersome, so
    /// this class can be used to define a <see cref="ILinearFlux"/> from a
    /// <see cref="INonlinearFlux"/>-object; 
    /// This class automatically evaluates
    /// a function matrix and an offset vector. 
    /// Of course, the functions defined by the <see cref="INonlinearFlux"/>-object 
    /// must affine linear to produce a correct result. This class doesn't tests this issue;
    /// </summary>
    public sealed class LinearFluxBuilder : ILinearFlux {

        /// <summary>
        /// constructor;
        /// </summary>
        /// <param name="a_Flux">
        /// affine-linear flux functions, implemented as
        /// <see cref="INonlinearFlux"/>
        /// </param>
        /// <param name="a_FluxTimeInDependent"><see cref="ILinearFlux.FluxTimeInDependent"/></param>
        /// <param name="a_FluxSpaceInDependent"><see cref="ILinearFlux.FluxSpaceInDependent"/></param>
        /// <param name="a_BorderFluxSpaceInDependent"><see cref="ILinearFlux.BorderFluxSpaceInDependent"/></param>
        /// <param name="a_BorderFluxTimeInDependent"><see cref="ILinearFlux.FluxTimeInDependent"/></param>
        public LinearFluxBuilder(INonlinearFlux a_Flux, 
            bool a_FluxTimeInDependent,  bool a_FluxSpaceInDependent,
            bool a_BorderFluxTimeInDependent, bool a_BorderFluxSpaceInDependent) {
            m_Flux = a_Flux;
            p_ArgumentOrdering = m_Flux.ArgumentOrdering;
            m_Arguments = new double[p_ArgumentOrdering.Length];
            m_Arguments2 = new double[p_ArgumentOrdering.Length];
            p_FluxTimeInDependent = a_FluxTimeInDependent;
            p_FluxSpaceInDependent = a_FluxSpaceInDependent;
            p_BorderFluxSpaceInDependent = a_BorderFluxSpaceInDependent;
            p_BorderFluxTimeInDependent = a_BorderFluxTimeInDependent;
        }


//        bool Own = true;
//        bool Neighbour = true;
        
        /// <summary>
        /// 
        /// </summary>
        public bool Local = false;
        
        /// <summary>
        /// 
        /// </summary>
        public bool _Flux = true;






        /// <summary>
        /// linear martices are constructed from this object.
        /// </summary>
        internal INonlinearFlux m_Flux;


        /// <summary>
        /// field containing value of coresponding Property
        /// </summary>
        private bool p_BorderFluxTimeInDependent = false;

        /// <summary>
        /// <see cref="ILinearFlux.BorderFluxTimeInDependent"/>
        /// can be specifyed in the constructor of this object;
        /// </summary>
        public bool BorderFluxTimeInDependent {
            get { return p_BorderFluxTimeInDependent; }
        }

        /// <summary>
        /// field containing value of coresponding Property
        /// </summary>
        private bool p_BorderFluxSpaceInDependent;

        /// <summary>
        /// <see cref="ILinearFlux.BorderFluxSpaceInDependent"/>
        /// can be specifyed in the constructor of this object;
        /// </summary>
        public bool BorderFluxSpaceInDependent {
            get { return p_BorderFluxSpaceInDependent; }
        }


        /// <summary>
        /// arguments used for the nonlinear flux functions to extract the function matrix
        /// </summary>
        double[] m_Arguments;
        
        /// <summary>
        /// arguments used for the nonlinear flux functions to extract the function matrix
        /// </summary>
        double[] m_Arguments2;


        /// <summary>
        /// <see cref="INonlinearFlux.BorderEdgeFlux"/>
        /// </summary>
        public void BorderEdgeFlux(double time, Vector2D x, Vector2D n, double[] FunctionMatrix, out double AffineOffset) {
            Array.Clear(m_Arguments, 0, m_Arguments.Length);
            AffineOffset = 0.0;

            AffineOffset = m_Flux.BorderEdgeFlux(time, x, n, m_Arguments);
            for (int i = m_Arguments.Length-1; i >= 0; i--) {
                m_Arguments[i] = 1;
                FunctionMatrix[i] = m_Flux.BorderEdgeFlux(time, x, n, m_Arguments);
                FunctionMatrix[i] -= AffineOffset;
                m_Arguments[i] = 0;
            }


            if ((!_Flux)) {
                Array.Clear(FunctionMatrix, 0, FunctionMatrix.Length);
                AffineOffset = 0.0;
            }
        }

        /// <summary>
        /// field containing value of coresponding Property
        /// </summary>
        bool p_FluxTimeInDependent;

        /// <summary>
        /// <see cref="ILinearFlux.FluxTimeInDependent"/>
        /// can be specifyed in the constructor of this object;
        /// </summary>
        public bool FluxTimeInDependent {
            get { return p_FluxTimeInDependent; }
        }

        /// <summary>
        /// field containing value of coresponding Property
        /// </summary>
        bool p_FluxSpaceInDependent;

        /// <summary>
        /// <see cref="ILinearFlux.FluxSpaceInDependent"/>
        /// can be specifyed in the constructor of this object;
        /// </summary>
        public bool FluxSpaceInDependent {
            get { return p_FluxSpaceInDependent; }
        }

        /// <summary>
        /// <see cref="INonlinearFlux.InnerEdgeFlux"/>
        /// </summary>
        public void InnerEdgeFlux(double time, Vector2D x, Vector2D n, double[] FunctionMatrixIn, double[] FunctionMatrixOut, out double AffineOffset) {
            Array.Clear(m_Arguments, 0, m_Arguments.Length);
            Array.Clear(m_Arguments2, 0, m_Arguments.Length);

            AffineOffset = m_Flux.InnerEdgeFlux(time, x, n, m_Arguments, m_Arguments2);
            for (int i = m_Arguments.Length-1; i >= 0; i--) {
                m_Arguments[i] = 1;
                FunctionMatrixIn[i] = m_Flux.InnerEdgeFlux(time, x, n, m_Arguments, m_Arguments2);
                FunctionMatrixIn[i] -= AffineOffset;
                m_Arguments[i] = 0;

                m_Arguments2[i] = 1;
                FunctionMatrixOut[i] = m_Flux.InnerEdgeFlux(time, x, n, m_Arguments, m_Arguments2);
                FunctionMatrixOut[i] -= AffineOffset;
                m_Arguments2[i] = 0;
            }

            if (!_Flux) {
                Array.Clear(FunctionMatrixIn, 0, FunctionMatrixIn.Length);
                Array.Clear(FunctionMatrixOut, 0, FunctionMatrixOut.Length);
                AffineOffset = 0.0;
            }
            //if (!Neighbour) {
            //    Array.Clear(FunctionMatrixOut, 0, FunctionMatrixOut.Length);
            //    AffineOffset = 0.0;
            //}
            //if (!Local) {
            //    Array.Clear(FunctionMatrixIn, 0, FunctionMatrixIn.Length);
            //    AffineOffset = 0.0;
            //}

        }

        /// <summary>
        /// <see cref="INonlinearFlux.Flux"/>
        /// </summary>
        public void Flux(double time, Vector2D x, Vector2D[] FunctionMatrix, out Vector2D AffineOffset) {
            Array.Clear(m_Arguments, 0, m_Arguments.Length);

            AffineOffset = m_Flux.Flux(time, x, m_Arguments);
            for (int i = m_Arguments.Length-1; i >= 0; i--) {
                m_Arguments[i] = 1;
                FunctionMatrix[i] = m_Flux.Flux(time, x, m_Arguments);
                FunctionMatrix[i].Sub(AffineOffset);
                m_Arguments[i] = 0;
            }


            if (!Local) {
                Array.Clear(FunctionMatrix, 0, FunctionMatrix.Length);
                AffineOffset = new Vector2D(0.0, 0.0);
            }

        }


        /// <summary>
        /// field containing value of coresponding Property
        /// </summary>
        Field[] p_ArgumentOrdering;

        /// <summary>
        /// <see cref="IEquationCommon.ArgumentOrdering"/>
        /// Simply taken from the <see cref="INonlinearFlux"/>-object givven as an argument
        /// for the construction of this object.
        /// </summary>
        public Field[] ArgumentOrdering {
            get { return p_ArgumentOrdering; }
        }

        /// <summary>
        /// This method can be used to test if the construction of the linear affine flux
        /// was correct
        /// </summary>
        /// <returns>if the given nonlinear flux is representable as an affine linear function,
        /// this should return 0 or a very small value</returns>
        public double Test() {
            Vector2D X;
            Vector2D n;

            int L = m_Flux.ArgumentOrdering.Length;
            int I = 20*L;

            Random rand = new Random();

            double[] Uin = new double[L];
            double[] Uout = new double[L];

            Vector2D o;
            Vector2D[] m1 = new Vector2D[L];
            Vector2D[] m2 = new Vector2D[L];

            double oSc;
            double[] m1Sc = new double[L];
            double[] m2Sc = new double[L];

            double errsum = 0.0;
            for (int i = 0; i < I; i++) {

                // populate Arguments with random values
                // -------------------------------------
                X.x = rand.NextDouble()*10.0 - 5.0;
                X.y = rand.NextDouble()*10.0 - 5.0;

                for (int l = L-1; l>= 0; l--) {
                    Uin[l] = rand.NextDouble()*2.0 - 1.0;
                    Uout[l] = rand.NextDouble()*2.0 - 1.0;
                }

                n.x = rand.NextDouble()*2.0 - 1.0;
                n.y = rand.NextDouble()*2.0 - 1.0;
                n.Normalize();





                // test "BorderEdgeFlux"
                // ---------------------
                {
                    // "nonlinear" version
                    double nbe = m_Flux.BorderEdgeFlux(0.0, X, n, Uin);

                    // linear version
                    this.BorderEdgeFlux(0.0, X, n, m1Sc, out oSc);

                    double lbe = 0.0;
                    lbe += oSc;

                    for (int l = L-1; l>= 0; l--) lbe += m1Sc[l]*Uin[l];

                    nbe -= lbe; // now, "nbe" should be equal to "lbe"
                    errsum += Math.Abs(nbe);
                }

                // test "Flux"
                // -----------
                {
                    // "nonlinear" version
                    Vector2D nf = m_Flux.Flux(0.0, X, Uin);

                    // linear version
                    this.Flux(0.0, X, m1, out o);

                    Vector2D lf = new Vector2D();
                    lf.Acc(o);

                    for (int l = L-1; l>= 0; l--) lf.Acc(m1[l], Uin[l]);

                    nf.Sub(lf); // now, "nf" should be equal to "lf"
                    errsum += nf.Abs();
                }

                // test "Flux"
                // -----------
                {
                    // "nonlinear" version
                    double nif = m_Flux.InnerEdgeFlux(0.0, X, n, Uin, Uout);

                    // linear version
                    this.InnerEdgeFlux(0.0, X, n, m1Sc, m2Sc, out oSc);

                    double lif = 0.0;
                    lif += oSc;

                    for (int l = L-1; l>= 0; l--) lif += m1Sc[l]*Uin[l];
                    for (int l = L-1; l>= 0; l--) lif += m2Sc[l]*Uout[l];

                    nif -= lif; // now, "nif" should be equal to "lif"
                    errsum += Math.Abs(nif);
                }
            }
            
            return errsum;
        }

    }
}
