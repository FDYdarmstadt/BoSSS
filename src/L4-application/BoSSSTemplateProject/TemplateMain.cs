using System;
using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.IO;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Solution.Template {

    /// <summary>
    /// A very basic implementation of upwinding.
    /// </summary>
    class SimpleUpwindFlux : LinearFlux {

        /// <summary>
        /// Hard-coded argument list for the flux.
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" };
            }
        }

        internal Func<double[], double, double>[] m_ConvectionVelocityField;


        /// <summary>
        /// Numerical flux at boundary edges.
        /// </summary>
        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            double uA = Uin[0];
            double uB = 0.0; // hard-coded boundary value

            int D = inp.D;

            double[] Vel = new double[D];
            Debug.Assert(m_ConvectionVelocityField.Length == D);
            double VelxN = 0;
            for(int d = 0; d < D; d++) {
                Vel[d] = m_ConvectionVelocityField[d](inp.X, inp.time);
                VelxN += Vel[d] * inp.Normale[d];
            }

            double uUpwnd;
            if(VelxN >= 0)
                uUpwnd = uA;
            else
                uUpwnd = uB;

            double Flx = uUpwnd * VelxN;

            return Flx;
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            double u = U[0];
            int D = inp.D;

            Debug.Assert(m_ConvectionVelocityField.Length == D);
            for(int d = 0; d < D; d++) {
                output[d] = m_ConvectionVelocityField[d](inp.Xglobal, inp.time)*u;
            }
            
        }

        /// <summary>
        /// Numerical flux at interior edges.
        /// </summary>
        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            double uA = Uin[0];
            double uB = Uout[0];

            int D = inp.D;

            double[] Vel = new double[D];
            Debug.Assert(m_ConvectionVelocityField.Length == D);
            double VelxN = 0;
            for(int d = 0; d < D; d++) {
                Vel[d] = m_ConvectionVelocityField[d](inp.X, inp.time);
                VelxN += Vel[d] * inp.Normale[d];
            }

            double uUpwnd;
            if(VelxN >= 0)
                uUpwnd = uA;
            else
                uUpwnd = uB;

            double Flx = uUpwnd*VelxN;

            return Flx;
        }
    }

    /// <summary>
    /// The control class for the application, where the simulation parameters are specified.
    /// </summary>
    public class TemplateControl : Control.AppControl {


        /// <summary>
        /// The only application-specific parameter for this application: a space and time dependent convection velocity.
        /// </summary>
        public Func<double[], double, double>[] ConvectionVelocityField;

    }



    /// <summary>
    /// Main application class.
    /// </summary>
    public class TemplateMain : Application<TemplateControl> {

        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {
            Application<TemplateControl>._Main(args, true, null, delegate () {
                return new TemplateMain();
            });
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            throw new NotImplementedException();
        }

        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
            throw new NotImplementedException();
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            throw new NotImplementedException();
        }

    }
}
