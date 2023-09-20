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

    public class CylindricalCoordTerms_NSEmomentum : BulkEquation {

        string speciesName;

        string codomainName;

        int d;

        int D;

        double rho;

        public CylindricalCoordTerms_NSEmomentum(
            string spcName,
            int d,
            int D,
            INSE_Configuration config) {

            speciesName = spcName;
            codomainName = EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));

            this.d = d;
            this.D = D;
            if (D != 3)
                throw new ArgumentOutOfRangeException("only supported for 3D");
            if (d < 0 || d >= D)
                throw new ArgumentOutOfRangeException();

            PhysicalParameters physParams = config.getPhysParams;

            // set species arguments
            double rhoSpc, muSpc;
            switch (spcName) {
                case "A": { rhoSpc = physParams.rho_A; muSpc = physParams.mu_A; break; }
                case "B": { rhoSpc = physParams.rho_B; muSpc = physParams.mu_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            rho = rhoSpc;

            //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(2));
            if (d == 0 || d == 1) { 
                var convectiveTerm = new CylindricalCoordTerm_convectiveTerm(d, rhoSpc);
                AddComponent(convectiveTerm);
                var viscousTerm = new CylindricalCoordTerm_viscousTerm(d, muSpc);
                AddComponent(viscousTerm);
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => rho;

        public override string CodomainName => codomainName;
    }


    public class CylindricalCoordTerm_convectiveTerm : IVolumeForm, ISupportsJacobianComponent {

        int m_d;

        double m_density;


        public CylindricalCoordTerm_convectiveTerm(int d, double density) {
            this.m_d = d;
            this.m_density = density;
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public virtual IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(2);
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return VariableNames.Velocity0Vector(2);
            }
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double VelComp = 0.0;

            if (m_d == 0) {
                VelComp = - cpv.Parameters[1] * U[1];
            }

            if (m_d == 1) {
                VelComp = 0.5 * (cpv.Parameters[0] * U[1]) + 0.5 * (cpv.Parameters[1] * U[0]);
            }

            return m_density * (VelComp / cpv.Xglobal[0]) * V;
        }

    }

    public class CylindricalCoordTerm_viscousTerm : IVolumeForm, ISupportsJacobianComponent {

        int m_d;

        double m_viscosity;


        public CylindricalCoordTerm_viscousTerm(int d, double viscosity) {
            this.m_d = d;
            this.m_viscosity = viscosity;
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public virtual IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(2);
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[0];
            }
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            return m_viscosity * (U[m_d] / cpv.Xglobal[0].Pow2()) * V;
        }

    }



    public class CylindricalCoordTerms_NSEconti : BulkEquation {

        string speciesName;

        string codomainName;


        public CylindricalCoordTerms_NSEconti(
            string spcName) {

            speciesName = spcName;
            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityX);


            var contiTerm = new CylindricalCoordTerm_contiTerm();
            AddComponent(contiTerm);

        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 0.0;

        public override string CodomainName => codomainName;
    }

    public class CylindricalCoordTerm_contiTerm : IVolumeForm, ISupportsJacobianComponent {


        public CylindricalCoordTerm_contiTerm() {
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.VelocityX };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[0];
            }
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            return (U[0] / cpv.Xglobal[0]) * V;
        }

    }
}
