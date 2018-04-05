using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;


namespace BoSSS.Solution.LevelSetTools.FourierLevelSet {

    /// <summary>
    /// Coefficients for AdamsBashforth schemes (AB)
    /// </summary>
    public class ABSchemeCoeff {

        /// <summary>
        /// Ab-coeffs for f[0], f[-1], f[-2], ...
        /// </summary>
        public double[] coeffs;

        /// <summary>
        /// number of AB-steps
        /// </summary>
        public int steps
        {
            get
            {
                return coeffs.Length;
            }
        }

        /// <summary>
        /// empty constructor
        /// </summary>
        public ABSchemeCoeff() {
        }

        /// <summary>
        /// returns the coefficients for a AdamsBashforth scheme of order <paramref name="order"/>.
        /// </summary>
        /// <param name="order"></param>
        static public ABSchemeCoeff AB(int order) {
            ABSchemeCoeff ABsc = new ABSchemeCoeff();

            switch (order) {
                case 1: {
                        ABsc.coeffs = new double[] { 1.0 };
                        break;
                    }
                case 2: {
                        ABsc.coeffs = new double[] { 3.0 / 2.0, -1.0 / 2.0 };
                        break;
                    }
                case 3: {
                        ABsc.coeffs = new double[] { 23.0 / 12.0 , - 4.0 / 3.0 , 5.0 / 12.0 };
                        break;
                    }
                case 4: {
                        ABsc.coeffs = new double[] { 55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -3.0 / 8.0 };
                        break;
                    }
                default:
                    throw new NotImplementedException("AB scheme of order " + order + " is not supported.");
            }

            if (Math.Abs(ABsc.coeffs.Sum() - 1.0) > 1.0e-10)
                throw new ApplicationException();

            return ABsc;

        }

    }


    /// <summary>
    /// A Fourier Levelset timestepper for the explicit multistep AdamsBashforth (AB) schemes
    /// </summary>
    public class FLS_Timestepper_AB : FourierLevSetTimestepper {


        ABSchemeCoeff[] ABsc_chain;

        public FLS_Timestepper_AB(FourierLevSetControl Control, double[] currentState, DelComputeChangerate _DelCompChange, DelEvolveFourier _DelEvolveFourier, int ABorder)
            : base(Control, currentState, _DelCompChange, _DelEvolveFourier) {

            ABsc_chain = new ABSchemeCoeff[ABorder];
            for (int i = 0; i < ABorder; i++) {
                ABsc_chain[i] = ABSchemeCoeff.AB(i + 1);
            }

        }

        public override void moveLevelSet(double dt, ConventionalDGField[] velocity) {

            //// compute the current change rate
            //double[] current_changerate = DelCompChange(dt, velocity);
            //changerates.setCurrent(current_changerate);

            //for (int s = 0; s < ABsc_chain[changerates.m_capacity].steps; s++) {
            //    FLSproperty[0].AccV(dt * ABsc_chain[changerates.m_capacity].coeffs[s], changerates[-s]);
            //}

            //DelEvolveFourier(FLSproperty[0].CloneAs());

        }


        public override void updateFourierLevSet() {

            //if (changerates[0] != null) {
            //    if (changerates.m_capacity < ABsc_chain.Length - 1)
            //        changerates.increaseCapacity();
            //    changerates.push();
            //}
        }

    }
}
