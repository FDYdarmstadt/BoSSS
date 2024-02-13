using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;

namespace BoSSS.Solution.NSECommon{
    
    public abstract class LinearizedScalarConvectionJacobiBase : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {

        protected int m_SpatialDimension;

        protected int argumentIndex;

        protected virtual double GetScalar(params double[] Arguments){
            return Arguments[m_SpatialDimension];
        }

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="SpatDim">Spatial dimension (either 2 or 3)</param>
        /// <param name="BcMap"></param>
        /// <param name="EoS">Null for multiphase. Has to be given for Low-Mach and combustion to calculate density.</param>
        /// <param name="Argument">Variable name of the argument (e.g. "Temperature" or "MassFraction0")</param>
        public LinearizedScalarConvectionJacobiBase(int SpatDim) {
            m_SpatialDimension = SpatDim;
        }

        /// <summary>
        /// flux at inner edges
        /// </summary>        
        protected double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0;

            double ScalarIn = GetScalar(Uin);
            double ScalarOut = GetScalar(Uout);

            for (int i = 0; i < m_SpatialDimension; ++i) {
                r += 0.5 * ScalarIn * Uin[i] * inp.Normal[i];
                r += 0.5 * ScalarOut * Uout[i] * inp.Normal[i];
            }


            // Calculate dissipative part
            // ==========================
            double[] VelocityMeanIn = Uin.GetSubVector(0, m_SpatialDimension); ////////////////////////////TODO CHECK!!!!!!!!!!!!!!!!!!!!
            double[] VelocityMeanOut = Uout.GetSubVector(0, m_SpatialDimension);

            double LambdaIn;
            double LambdaOut;
            double TemperatureMeanIn = Uin[m_SpatialDimension];
            double TemperatureMeanOut = Uout[m_SpatialDimension];

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal, false, ScalarIn);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal, false, ScalarOut);

            double Lambda = Math.Max(LambdaIn, LambdaOut);

            r += 0.5 * Lambda * (ScalarIn - ScalarOut);
            if(double.IsNaN(r))
                throw new NotFiniteNumberException();

            return r;
        }

        /// <summary>
        /// flux at the boundary
        /// </summary>
        protected abstract double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin); 

        /// <summary>
        /// returns
        /// \f[ 
        ///   \vec{v} \cdot \phi,
        /// \f]
        /// where \f$ \vec{v}\f$  is the linearization point.
        /// </summary>
        protected void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {

            for (int i = 0; i < m_SpatialDimension; ++i) {
                output[i] = U[i] * GetScalar(U);
                if(double.IsNaN(output[i]))
                    throw new NotFiniteNumberException();
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return this.InnerEdgeFlux(ref inp, _uIN, _uOUT) * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return this.BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            int D = GradV.Length;
            double acc = 0;
            var buf = new double[D];
            this.Flux(ref cpv, U, buf);
            for(int d = 0; d < D; d++)
                acc += buf[d] * GradV[d];
            return -acc;
        }



        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DerivEdg, DerivVol };
        }

        /// <summary>
        /// Level-Set for multiphase.
        /// Temperature for convection in the temperature equation
        /// Name of the computed MassFraction in MassFraction balances (e.g. MassFraction0)
        /// </summary>
        public abstract IList<string> ArgumentOrdering{ get; }

        public abstract IList<string> ParameterOrdering{ get; } 

        /// <summary>
        /// <see cref="IEdgeForm.BoundaryEdgeTerms"/>
        /// </summary>
        virtual public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        /// <see cref="IEdgeForm.InnerEdgeTerms"/>
        /// </summary>
        virtual public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /// <summary>
        /// <see cref="IVolumeForm.VolTerms"/>
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxGradV | TermActivationFlags.GradV;
            }
        }
    }
}
