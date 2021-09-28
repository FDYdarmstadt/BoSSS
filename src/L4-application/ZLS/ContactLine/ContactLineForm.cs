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

        protected static Vector SurfaceNormal(double[] parameters, int D, int offSet) {

            Vector N = new Vector(D);

            for (int d = 0; d < D; d++) {
                N[d] = parameters[d + offSet];
            }
            N.NormalizeInPlace();
            return N;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            // only parameter dependent, leave this empty
            return new IEquationComponent[] { };
        }
    }
}
