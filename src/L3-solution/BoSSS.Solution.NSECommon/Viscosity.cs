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
using System.Linq;
using System.Text;
using BoSSS.Solution.Utils;
using BoSSS.Foundation;

namespace BoSSS.Solution.NSECommon {
    /*

    /// <summary>
    /// viscous operators for incompressible fluids with constant density and viscosity, i.e.
    /// discretization of the operator \f$  u \mapsto \frac{-1}{\textrm{Re}} \cdot \Delta u \f$  with
    /// boundary conditions
    /// \f$ 
    /// \left\{
    /// \begin{array}{rcll}
    /// u_d|_{\Gamma_\textrm{Inl} \cup \Gamma_\textrm{Wall}}  &amp; = &amp;
    ///             \left[ \vec{u}_\textrm{Inl,Wall} \right]_d  
    ///      &amp;  \textrm{ on } \Gamma_\textrm{Inl} \cup \Gamma_\textrm{Wall} \\
    /// (\vec{n} \cdot \nabla u_d)|_{\Gamma_\textrm{Out} \cup \Gamma_\textrm{POlt}}  &amp; = &amp;
    ///             0                          
    ///      &amp;  \textrm{ on } \Gamma_\textrm{Out} \cup \Gamma_\textrm{POlt} \\
    /// \end{array}
    /// \right. .
    /// \f$ 
    /// </summary>
    public class Viscosity : BoSSS.Solution.ZooOfOperators._2ndOrder.ipLaplace {

        double m_neg1overRe;
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        int m_d;
        Func<double[],double,double>[] velFunction;

        bool m_UseBoundaryVelocityParameter;
        bool m_VariableDensity;

        /// <summary>
        /// Ctor for common part of incompressible and low Mach number flows.
        /// </summary>
        /// <param name="penaltyMultiplyer"></param>
        /// <param name="bcmap"></param>
        /// <param name="d">velocity component index;
        /// affects only the boundary condition (and therefore, only the affine offset of the operator).
        /// </param>
        /// <param name="Rey">
        /// Reynolds number
        /// </param>
        /// <param name="UseBoundaryVelocityParameter">
        /// true, if (an offset to) the boundary velocity is supplied in form of parameter variables
        /// </param>
        /// <param name="VariableDensity">
        /// True for low Mach number flows.
        /// </param>
        private Viscosity(double penaltyMultiplyer, double Rey, IncompressibleBoundaryCondMap bcmap, int d, bool UseBoundaryVelocityParameter, bool VariableDensity)
            : base(penaltyMultiplyer, VariableNames.Velocity_d(d)) //
        {

            if (double.IsNaN(Rey) || double.IsInfinity(Rey) || Rey <= 0)
                throw new ArgumentException("Reynolds number out of range: " + Rey);

            m_neg1overRe = -1.0 / Rey;
            m_BcMap = bcmap;
            m_d = d;
            velFunction = bcmap.bndFunction[VariableNames.Velocity_d(d)];

            m_UseBoundaryVelocityParameter = UseBoundaryVelocityParameter;
            m_VariableDensity = VariableDensity;
        }

        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>
        /// <param name="penaltyMultiplyer"></param>
        /// <param name="Rey"></param>
        /// <param name="bcmap"></param>
        /// <param name="d"></param>
        /// <param name="UseBoundaryVelocityParameter"></param>
        public Viscosity(double penaltyMultiplyer, double Rey, IncompressibleBoundaryCondMap bcmap, int d, bool UseBoundaryVelocityParameter)
            : this(penaltyMultiplyer, Rey, bcmap, d, UseBoundaryVelocityParameter, false) {

            if (m_UseBoundaryVelocityParameter)
                parameterOrdering = new string[] { VariableNames.BoundaryVelocity_d(d) };
            else
                parameterOrdering = null;
        }

        MaterialLaw m_EoS = null;

        /// <summary>
        /// Ctor for low Mach number flows.
        /// </summary>
        /// <param name="penaltyMultiplyer"></param>
        /// <param name="Rey"></param>
        /// <param name="bcmap"></param>
        /// <param name="d"></param>
        /// <param name="EoS"></param>
        public Viscosity(double penaltyMultiplyer, double Rey, IncompressibleBoundaryCondMap bcmap, int d, MaterialLaw EoS)
            : this(penaltyMultiplyer, Rey, bcmap, d, false, true) {

            m_EoS = EoS;
            parameterOrdering = new string[] { VariableNames.Phi };
        }

        /// <summary>
        /// the value for inhomogeneous Dirichlet boundary conditions are taken from the parameter "g"
        /// </summary>
        override public IList<string> ParameterOrdering {
            get {
                return parameterOrdering;
            }
        }

        private string[] parameterOrdering;

        /// <summary>
        /// Dirichlet boundary value: the given velocity at the boundary.
        /// </summary>
        protected override double g_Diri(ref CommonParamsBnd inp) {
            Func<double[],double,double> boundVel = velFunction[inp.EdgeTag];
            double ret = boundVel(inp.X, 0);

            if (m_UseBoundaryVelocityParameter)
                ret += inp.Parameters_IN[0];

            return ret;
        }

        /// <summary>
        /// false for outflow and pressure outlet regions
        /// </summary>
        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            IncompressibleBcType edgType = m_BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                    // Atmospheric outlet/pressure outflow: hom. Neumann
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    return false;

                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++
                    return true;

                default:
                    throw new NotImplementedException("unsupported/unknown b.c. - missing implementation;");
            }
        }

        /// <summary>
        /// \f$ \frac{-1}/\frac{\textrm{Re}}\f$ 
        /// resp. \f$ \frac{-1}/\frac{\textrm{Re}} \eta\f$  for low Mach number flows.
        /// </summary>
        override public double Nu(double[] x, double[] parameters) {
            double nu;

            if (m_VariableDensity) {
                //dynamic viscosity
                double eta = m_EoS.GetViscosity(parameters[0]);
                nu = m_neg1overRe * eta;
            } else {
                nu = m_neg1overRe;
            }

            return nu;
        }
    }
    */
}
