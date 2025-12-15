/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XNSECommon {

    public class InertiaForceTermsInRotSystem : BulkEquation {

        string speciesName;

        string codomainName;

        int d;

        int D;

        double rho;

        public InertiaForceTermsInRotSystem(
            string spcName,
            int d,
            int D,
            double[] angVelocity,
            INSE_Configuration config) {

            speciesName = spcName;
            codomainName = EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(2));

            this.d = d;
            this.D = D;
            if (D != 3)
                throw new ArgumentOutOfRangeException("only supported for 3D");
            if (d < 0 || d >= 2)
                throw new ArgumentOutOfRangeException("only supported for 0 and 1 component");

            PhysicalParameters physParams = config.getPhysParams;

            // set species arguments
            double rhoSpc;
            switch (spcName) {
                case "A": { rhoSpc = physParams.rho_A; break; }
                case "B": { rhoSpc = physParams.rho_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            rho = rhoSpc;

            var centrifugalForce = new CentrifugalForceTerm_RotatingDisk(d, rhoSpc, angVelocity[2]);
            AddComponent(centrifugalForce);

            var coriolisForce = new CoriolisForceTerm_RotatingDisk(d, rhoSpc, angVelocity[2]);
            AddComponent(coriolisForce);

        }

        public override string SpeciesName => speciesName;

        public override double MassScale => rho;

        public override string CodomainName => codomainName;

    }



    public class CentrifugalForceTerm_RotatingDisk : InertiaForceTermsInRotSystemComponent {


        public CentrifugalForceTerm_RotatingDisk(int comp, double density, double angularVel) : base(comp, density, angularVel) { }


        public override TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }


        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            (double rOnDsik, double theta) = GetCylindricCoords(cpv.Xglobal);
            double dirComp = (m_comp == 0) ? Math.Cos(theta) : Math.Sin(theta);
            return m_density * m_angularVelocity.Pow2() * rOnDsik * dirComp * V;
        }

    }


    public class CoriolisForceTerm_RotatingDisk : InertiaForceTermsInRotSystemComponent {


        public CoriolisForceTerm_RotatingDisk(int comp, double density, double angularVel) : base(comp, density, angularVel) { }


        public override TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.VelocityX, VariableNames.VelocityY };
            }
        }


        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            (double rOnDsik, double theta) = GetCylindricCoords(cpv.Xglobal);
            double radialComp = (m_comp == 0) ? - U[1] * Math.Cos(theta) : - U[1] * Math.Sin(theta);
            double azimuthComp = (m_comp == 0) ? - U[0] * Math.Sin(theta) : U[0] * Math.Cos(theta);
            return m_density * 2.0 * m_angularVelocity * (radialComp + azimuthComp) * V;
        }

    }



    public abstract class InertiaForceTermsInRotSystemComponent : IVolumeForm, ISupportsJacobianComponent {

        protected int m_comp;

        protected double m_density;

        protected double m_angularVelocity;

        public InertiaForceTermsInRotSystemComponent(int comp, double density, double angularVel) {

            if (comp < 0 || comp > 1)
                throw new ArgumentOutOfRangeException("only supported for components 0 and 1");

            this.m_comp = comp; 
            this.m_density = density;
            this.m_angularVelocity = angularVel;
        }


        public virtual TermActivationFlags VolTerms { get; }

        public virtual IList<string> ArgumentOrdering { get; }

        public virtual IList<string> ParameterOrdering {
            get {
                return new string[0];
            }
        }

        public virtual double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            throw new NotImplementedException();
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        protected (double radiusOnDisk, double theta) GetCylindricCoords(double[] Xglobal) {
            double rOnDisk = Math.Sqrt(Xglobal[0].Pow2() + Xglobal[1].Pow2());
            double theta = Math.Atan2(Xglobal[1], Xglobal[0]);

            return (rOnDisk, theta);
        }

    }

}
