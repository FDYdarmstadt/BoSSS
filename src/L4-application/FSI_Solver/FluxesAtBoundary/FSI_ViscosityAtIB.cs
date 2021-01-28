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
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using BoSSS.Foundation;
using ilPSP.Utils;
using BoSSS.Platform.LinAlg;
using FSI_Solver;
using ilPSP;
using System.Linq;

namespace BoSSS.Solution.NSECommon.Operator.Viscosity {

    public class FSI_ViscosityAtIB : ILevelSetForm, ILevelSetEquationComponentCoefficient {
        public FSI_ViscosityAtIB(int currentDim, int spatialDim, LevelSetTracker levelSetTracker, double penalty, double fluidViscosity, Func<Vector, FSI_ParameterAtIB> getParticleParams) {
            m_penalty = penalty;
            m_LsTrk = levelSetTracker;
            muA = fluidViscosity;
            component = currentDim;
            m_GetParticleParams = getParticleParams;
            m_D = spatialDim;
        }

        private readonly int component;
        private readonly int m_D;
        private readonly LevelSetTracker m_LsTrk;

        /// <summary>
        /// Describes: 0: velX, 1: velY, 2: rotVel, 3: particleradius, 4: active_stress, 5: first scaling parameter, 6: second scaling parameter, 7: current angle
        /// </summary>
        private readonly Func<Vector, FSI_ParameterAtIB> m_GetParticleParams;

        /// <summary>
        /// Viscosity in species A
        /// </summary>
        private readonly double muA;
        private readonly double m_penalty; // safty factor
        double m_penalty_degree; // DG degree scaling
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate"/>;
        /// basic adjustment of factors for penalty computation
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = csA.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_degree = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
        }
        
        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
        /// </summary>
        double Penalty(int jCellIn) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];

            double penaltySizeFactor = penaltySizeFactor_A;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(m_penalty_degree));

            double scaledPenalty = penaltySizeFactor * m_penalty_degree * m_penalty;
            if(scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

        }




        private enum BoundaryConditionType {
            passive = 0,
            active = 1
        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            int dim = inp.Normal.Dim;
            double _penalty = Penalty(inp.jCellIn);

            // Particle parameters
            // =====================
            if (inp.X.IsNullOrEmpty())
                throw new Exception("X is null or empty");
            if (m_GetParticleParams == null)
                throw new Exception("m_GetParticleParams is null or empty");
            if (inp.X.Abs() < 0)
                throw new ArithmeticException("invalid length of position vector");
            FSI_ParameterAtIB coupling = m_GetParticleParams(inp.X);
            if (coupling == null)
                throw new Exception("coupling is null or empty");
            Vector orientation = new Vector(Math.Cos(coupling.Angle()), Math.Sin(coupling.Angle()));
            Vector orientationNormal = new Vector(-Math.Sin(coupling.Angle()), Math.Cos(coupling.Angle()));
            Vector activeStressVector = orientationNormal * inp.Normal > 0 ? new Vector(-coupling.ActiveStress() * inp.Normal[1], coupling.ActiveStress() * inp.Normal[0]) 
                                                                  : new Vector(coupling.ActiveStress() * inp.Normal[1], -coupling.ActiveStress() * inp.Normal[0]);
            BoundaryConditionType bcType = orientation * inp.Normal <= 0 || coupling.ActiveStress() == 0 ? BoundaryConditionType.passive : BoundaryConditionType.active;

            Debug.Assert(ArgumentOrdering.Count == dim);
            Debug.Assert(Grad_uA.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == dim);
            Debug.Assert(Grad_uB.GetLength(1) == dim);

            // Gradient of u and v 
            // =====================
            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < dim; d++) {
                Grad_uA_xN += Grad_uA[component, d] * inp.Normal[d];
                Grad_vA_xN += Grad_vA[d] * inp.Normal[d];
            }

            double returnValue = 0.0;

            // 3D for IBM_Solver
            // =====================
            if (dim == 3) {
                returnValue -= Grad_uA_xN * (vA);                                                    // consistency term
                returnValue -= Grad_vA_xN * (uA[component] - 0);                                     // symmetry term
                returnValue += _penalty * (uA[component] - 0) * (vA);                                // penalty term
                Debug.Assert(!(double.IsInfinity(returnValue) || double.IsNaN(returnValue)));
                return returnValue * muA;
            }

            // 2D
            // =====================
            Vector uAFict = coupling.VelocityAtPointOnLevelSet(); 

            switch (bcType) {
                case BoundaryConditionType.passive: {
                        for (int d = 0; d < dim; d++) {
                            returnValue -= muA * Grad_uA[component, d] * vA * inp.Normal[d];
                            returnValue -= muA * Grad_vA[d] * (uA[component] - uAFict[component]) * inp.Normal[d];
                        }
                        returnValue += muA * (uA[component] - uAFict[component]) * vA * _penalty;
                        break;
                    }

                case BoundaryConditionType.active: {
                        // normal direction, solid wall
                        for (int dN = 0; dN < dim; dN++) {
                            for (int dD = 0; dD < dim; dD++) {
                                // consistency term
                                returnValue -= muA * (inp.Normal[dN] * Grad_uA[dN, dD] * inp.Normal[dD]) * vA * inp.Normal[component];
                                // symmetry term
                                returnValue -= muA * (Grad_vA[dD] * inp.Normal[dD]) * (inp.Normal[dN] * uA[dN] - inp.Normal[dN] * uAFict[dN]) * inp.Normal[component];      
                            }
                            // penalty term
                            returnValue += muA * inp.Normal[dN] * (uA[dN] - uAFict[dN]) * inp.Normal[component] * vA * _penalty;                  
                        }
                        // tangential direction, active part
                        double[,] P = new double[dim, dim];
                        for (int d1 = 0; d1 < dim; d1++) {
                            for (int d2 = 0; d2 < dim; d2++) {
                                if (d1 == d2) {
                                    P[d1, d2] = 1 - inp.Normal[d1] * inp.Normal[d2];
                                }
                                else {
                                    P[d1, d2] = inp.Normal[d1] * inp.Normal[d2];
                                }
                            }
                        }
                        for (int d1 = 0; d1 < dim; d1++) {
                            for (int d2 = 0; d2 < dim; d2++) {
                                returnValue -= P[d1, d2] * activeStressVector[d2] * (P[d1, component] * vA); 
                            }
                        }
                        break;
                    }
            }
            Debug.Assert(!(double.IsInfinity(returnValue) || double.IsNaN(returnValue)));
            return returnValue;
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(this.m_D); }
        }


        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
    }
}