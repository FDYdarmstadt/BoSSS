using BoSSS.Foundation;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using NSEVariableNames = BoSSS.Solution.NSECommon.VariableNames;

namespace ZwoLevelSetSolver.ContactLine {
    abstract class ContactLineForm : IVolumeForm, ISupportsJacobianComponent {

        protected int D;

        public ContactLineForm(int D) {
            this.D = D;
        }

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public IList<string> ArgumentOrdering {
            get {
                return NSEVariableNames.VelocityVector(D);
            }
        }

        public virtual IList<string> ParameterOrdering {
            get {
                string[] parameters = NSEVariableNames.NormalVector(D);
                parameters = parameters.Cat(NSEVariableNames.AsLevelSetVariable(
                    VariableNames.SolidLevelSetCG, NSEVariableNames.NormalVector(D)));
                parameters = parameters.Cat(NSEVariableNames.MaxSigma);
                return parameters;
            }
        }

        public abstract double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV);

        protected Vector FluidSurfaceNormal(ref CommonParamsVol cpv) {
            return SurfaceNormal(cpv.Parameters, D, 0);
        }

        protected Vector SolidSurfaceNormal(ref CommonParamsVol cpv) {
            return SurfaceNormal(cpv.Parameters, D, D);
        }

        protected static double Sigma(ref CommonParamsVol cpv) {
            return cpv.Parameters[cpv.Parameters.Length - 1];
        }

        static Vector SurfaceNormal(double[] parameters, int D, int offSet) {

            Vector N = new Vector(D);

            for(int d = 0; d < D; d++) {
                N[d] = parameters[d + offSet];
            }
            N = N.Normalize();
            return N;
        }

        protected static Vector Tangent(Vector Nsurf, Vector Nedge) {
            Debug.Assert(Nsurf.Dim == Nedge.Dim);

            int D = Nsurf.Dim;

            Vector tau = new Vector(D);
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if(d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            tau = tau.Normalize();
            return tau;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            // only parameter dependent, leave this empty
            return new IEquationComponent[] { };
        }
    }
}
