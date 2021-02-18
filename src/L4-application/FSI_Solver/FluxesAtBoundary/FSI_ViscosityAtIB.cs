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
using BoSSS.Application.FSI_Solver;
using System.Linq;
using ilPSP.Tracing;

namespace BoSSS.Solution.NSECommon.Operator.Viscosity {

    public class FSI_ViscosityAtIB : ILevelSetForm, ILevelSetEquationComponentCoefficient {
        public FSI_ViscosityAtIB(int currentDim, int spatialDim, LevelSetTracker levelSetTracker, double penalty, double fluidViscosity, Particle[] allParticles, double minGridLength) {
            m_penalty = penalty;
            LsTrk = levelSetTracker;
            FluidViscosity = fluidViscosity;
            Component = currentDim;
            MinGridLength = minGridLength;
            AllParticles = allParticles;
            SpatialDim = spatialDim;
        }

        private readonly int Component;
        private readonly int SpatialDim;
        private readonly double MinGridLength;
        private readonly LevelSetTracker LsTrk;
        private readonly Particle[] AllParticles;
        private readonly double FluidViscosity;
        private readonly double m_penalty; // safety factor
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
            Vector normalVector = inp.Normal;
            if (inp.X.IsNullOrEmpty())
                throw new Exception("X is null or empty");
            if (inp.X.Abs() < 0)
                throw new ArithmeticException("invalid length of position vector");

            // Particle parameters
            // =====================
            Particle currentParticle = AllParticles[0];
            for (int p = 0; p < AllParticles.Length; p++) {
                if (AllParticles[p].Contains(inp.X, 1.5 * MinGridLength)) {
                    currentParticle = AllParticles[p];
                    break;
                }
            }
            double activeStress = currentParticle.ActiveStress;
            double angle = currentParticle.Motion.GetAngle(0);
            Vector orientation = new Vector(Math.Cos(angle), Math.Sin(angle));
            Vector orientationNormal = new Vector(-orientation[1], orientation[0]);
            //Vector activeStressVector = orientationNormal * normalVector > 0 ? new Vector(-activeStress * normalVector[1], activeStress * normalVector[0]) : new Vector(activeStress * normalVector[1], -activeStress * normalVector[0]);
            Vector activeStressVector = -activeStress * orientation;
            BoundaryConditionType bcType = orientation * normalVector <= 0 || activeStress == 0 ? BoundaryConditionType.passive : BoundaryConditionType.active;

            Debug.Assert(ArgumentOrdering.Count == dim);
            Debug.Assert(Grad_uA.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == dim);
            Debug.Assert(Grad_uB.GetLength(1) == dim);

            // Gradient of u and v 
            // =====================
            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < dim; d++) {
                Grad_uA_xN += Grad_uA[Component, d] * inp.Normal[d];
                Grad_vA_xN += Grad_vA[d] * inp.Normal[d];
            }

            double returnValue = 0.0;

            // 3D for IBM_Solver
            // =====================
            if (dim == 3) {
                returnValue -= Grad_uA_xN * (vA);                                                    // consistency term
                returnValue -= Grad_vA_xN * (uA[Component] - 0);                                     // symmetry term
                returnValue += _penalty * (uA[Component] - 0) * (vA);                           // penalty term
                Debug.Assert(!(double.IsInfinity(returnValue) || double.IsNaN(returnValue)));
                return returnValue * FluidViscosity;
            }

            // 2D
            // =====================
            Vector radialVector = currentParticle.CalculateRadialVector(inp.X);
            Vector uAFict = new Vector(currentParticle.Motion.GetTranslationalVelocity(0)[0] - currentParticle.Motion.GetRotationalVelocity(0) * radialVector[1],
                              currentParticle.Motion.GetTranslationalVelocity(0)[1] + currentParticle.Motion.GetRotationalVelocity(0) * radialVector[0]);

            switch (bcType) {
                case BoundaryConditionType.passive: {
                    for (int d = 0; d < dim; d++) {
                        returnValue -= FluidViscosity * Grad_uA[Component, d] * vA * normalVector[d];
                        returnValue -= FluidViscosity * Grad_vA[d] * (uA[Component] - uAFict[Component]) * normalVector[d];
                    }
                    returnValue += FluidViscosity * (uA[Component] - uAFict[Component]) * vA * _penalty;
                    break;
                }

                case BoundaryConditionType.active: {
                    // normal direction, solid wall
                    for (int dN = 0; dN < dim; dN++) {
                        for (int dD = 0; dD < dim; dD++) {
                            // consistency term
                            returnValue -= FluidViscosity * (normalVector[dN] * Grad_uA[dN, dD] * normalVector[dD]) * vA * normalVector[Component];
                            // symmetry term
                            returnValue -= FluidViscosity * (Grad_vA[dD] * normalVector[dD]) * (normalVector[dN] * uA[dN] - normalVector[dN] * uAFict[dN]) * normalVector[Component];
                        }
                        // penalty term
                        returnValue += FluidViscosity * inp.Normal[dN] * (uA[dN] - uAFict[dN]) * normalVector[Component] * vA * _penalty;
                    }
                    // tangential direction, active part
                    double[,] P = new double[dim, dim];
                    for (int d1 = 0; d1 < dim; d1++) {
                        for (int d2 = 0; d2 < dim; d2++) {
                            if (d1 == d2) {
                                P[d1, d2] = 1 - normalVector[d1] * normalVector[d2];
                            } else {
                                P[d1, d2] = normalVector[d1] * normalVector[d2];
                            }
                        }
                    }
                    for (int d1 = 0; d1 < dim; d1++) {
                        for (int d2 = 0; d2 < dim; d2++) {
                            returnValue -= P[d1, d2] * activeStressVector[d2] * (P[d1, Component] * vA);
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
            get { return VariableNames.VelocityVector(this.SpatialDim); }
        }


        public SpeciesId PositiveSpecies {
            get { return LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return LsTrk.GetSpeciesId("A"); }
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