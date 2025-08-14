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

using System;
using System.Collections.Generic;
using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Control;


namespace BoSSS.Application.XNSE_Solver.DongBoundaryConditionTests {


    public static class KovasznayFlowSolutions {

        public static readonly Formula KovasznayFlow_u = new Formula(
            "VelX",
            false,
            "double VelX(double[] X) { " +
            "double lambda = 20.0 - 2.0 * Math.Sqrt(Math.PI.Pow2() + 100.0); " +
            "double velX = 1.0 - (Math.Exp(lambda * X[0]) * Math.Cos(2.0 * Math.PI * X[1]));" +
            "return velX; } "
        );

        public static readonly Formula KovasznayFlow_v = new Formula(
            "VelY",
            false,
            "double VelY(double[] X) { " +
            "double lambda = 20.0 - 2.0 * Math.Sqrt(Math.PI.Pow2() + 100.0); " +
            "double velY = (lambda/(2.0 * Math.PI)) * (Math.Exp(lambda * X[0]) * Math.Sin(2.0 * Math.PI * X[1]));" +
            "return velY; } "
        );

        public static readonly Formula KovasznayFlow_p = new Formula(
            "Pres",
            false,
            "double Pres(double[] X) { " +
            "double lambda = 20.0 - 2.0 * Math.Sqrt(Math.PI.Pow2() + 100.0); " +
            "double Pres = 0.5 * (1.0 - Math.Exp(2.0 * lambda * X[0]));" +
            "return Pres; } "
        );

    }


    public class DongBoundaryCondition_KovasznayFlow : IEdgeForm, ISupportsJacobianComponent, ISpeciesFilter {


        public DongBoundaryCondition_KovasznayFlow(string spcName, int d, double rho, double mu, IncompressibleMultiphaseBoundaryCondMap BcMap) {
            ValidSpecies = spcName;
            m_d = d;
            m_rho = rho;
            m_mu = mu;

            m_BcMap = BcMap;
        }

        int m_d;
        double m_rho;
        double m_mu;

        IncompressibleMultiphaseBoundaryCondMap m_BcMap;


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Flx_InCell = 0.0;

            double nu = m_mu / m_rho;
            double lambda = (1.0 / (2.0 * nu)) - Math.Sqrt(1.0 / (4.0 * nu.Pow2()) + 4.0 * Math.PI.Pow2());

            double x = inp.X[0];
            double y = inp.X[1];

            double p = 0.5 * (1.0 - Math.Exp(2 * lambda * x));

            double u = 1.0 - Math.Exp(lambda * x) * Math.Cos(2.0 * Math.PI * y);
            double v = (lambda / (2.0 * Math.PI)) * Math.Exp(lambda * x) * Math.Sin(2.0 * Math.PI * y);

            //Matrix(2, 2, [[-lambda*exp(lambda*x)*cos(2*pi*y), 2*exp(lambda*x)*pi*sin(2*pi*y)], [1/2*lambda^2*exp(lambda*x)*sin(2*pi*y)/pi, lambda*exp(lambda*x)*cos(2*pi*y)]])

            double u_x = -lambda * Math.Exp(lambda * x) * Math.Cos(2.0 * Math.PI * y);
            double u_y = 2.0 * Math.Exp(lambda * x) * Math.PI * Math.Sin(2.0 * Math.PI * y);

            double v_x = 0.5 * lambda.Pow2() * Math.Exp(lambda * x) * Math.Sin(2.0 * Math.PI * y) / Math.PI;
            double v_y = lambda * Math.Exp(lambda * x) * Math.Cos(2.0 * Math.PI * y);

            double[,] GradU = new double[2, 2] { { u_x, u_y }, { v_x, v_y } };

            if (m_BcMap.EdgeTag2Type[inp.EdgeTag] == IncompressibleBcType.Pressure_Outlet) {
                // pressure gradient
                Flx_InCell += -p * inp.Normal[m_d];
                // viscous terms 
                for (int d = 0; d < 2; d++) {
                    Flx_InCell += m_mu * (GradU[m_d, d]) * inp.Normal[d]; // + GradU[d, m_d]) * inp.Normal[d];
                }
            }

            if (m_BcMap.EdgeTag2Type[inp.EdgeTag] == IncompressibleBcType.Dong_OutFlow) {

                // pressure gradient
                Flx_InCell += -p * inp.Normal[m_d];

                for (int d = 0; d < 2; d++) {
                    // viscous terms 
                    Flx_InCell += m_mu * (GradU[m_d, d]) * inp.Normal[d]; // + GradU[d, m_d]) * inp.Normal[d];
                }

                // Dong term
                double U0 = 1.0;
                double delta = 1.0 / 20.0;
                double ndotu = (u * inp.Normal[0] + v * inp.Normal[1]);
                double Sout = 0.5 * (1.0 - Math.Tanh(ndotu / (U0 * delta)));
                double uAbs2 = u.Pow2() + v.Pow2();

                Flx_InCell -= 0.5 * uAbs2 * Sout * inp.Normal[m_d];

            }

            return -Flx_InCell * _vA;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            return 0.0;
        }

        public IList<string> ArgumentOrdering { get { return new string[0]; } }

        public IList<string> ParameterOrdering { get { return new string[0]; } }


        virtual public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.V; }
        }

        virtual public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.None; }
        }

        public string ValidSpecies {
            get;
            private set;
        }

        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { new EdgeFormDifferentiator(this, SpatialDimension) };
        }


    }


    public class ManufacturedSolutionComp : BulkEquation {

        public ManufacturedSolutionComp(int d, IncompressibleMultiphaseBoundaryCondMap BcMap, XNSE_OperatorConfiguration config, string spcName) {
            m_CodomainName = EquationNames.MomentumEquationComponent(d);
            m_SeciesName = spcName;

            PhysicalParameters physParams = config.getPhysParams;

            double rhoSpc, muSpc;
            switch (spcName) {
                case "A": { rhoSpc = physParams.rho_A; muSpc = physParams.mu_A; break; }
                case "B": { rhoSpc = physParams.rho_B; muSpc = physParams.mu_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }
            m_rho = rhoSpc;

            AddComponent(new DongBoundaryCondition_KovasznayFlow(spcName, d, rhoSpc, muSpc, BcMap));
        }


        string m_CodomainName;

        string m_SeciesName;

        double m_rho;

        public override string CodomainName { get { return m_CodomainName; } }

        public override string SpeciesName { get { return m_SeciesName; } }

        public override double MassScale { get { return m_rho; } }


    }

}
