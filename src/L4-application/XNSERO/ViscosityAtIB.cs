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
using ilPSP;
using System.Linq;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Application.XNSERO_Solver {
    public class ViscosityAtIB : ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient, IMultitreadSafety {

        //LevelSetTracker m_LsTrk;

        public ViscosityAtIB(int _d, int _D, Particle[] Particles, double penalty_base, double _muA, int iLevSet, string FluidSpc, string SolidSpecies, bool UsePhoretic) {
            m_penalty_base = penalty_base;
            this.FluidViscosity = _muA;
            this.Component = _d;
            this.m_D = _D;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            this.m_UsePhoretic = UsePhoretic;
            this.Particles = Particles;
            this.ActiveStress = Particles[0].ActiveStress;
        }
        private readonly int m_iLevSet;
        private readonly string m_FluidSpc;
        private readonly string m_SolidSpecies;
        private readonly int Component;
        private readonly int m_D;
        private readonly bool m_UsePhoretic;
        private readonly double ActiveStress;
        private readonly Particle[] Particles;

        /// <summary>
        /// Viskosity in species A
        /// </summary>
        private readonly double FluidViscosity;

        /// <summary>
        /// safety factor
        /// </summary>
        private readonly double m_penalty_base;

        /// <summary>
        /// degree and spatial dimension
        /// </summary>
        private double m_penalty_degree;


        /// <summary>
        /// length scale for negative species
        /// </summary>
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate"/>
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = csA.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_degree = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;

            IdentifyParticleCache_P = null;
        }


        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
        /// </summary>
        protected double Penalty(int jCellIn) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];

            double penaltySizeFactor = penaltySizeFactor_A;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(m_penalty_degree));

            double scaledPenalty = penaltySizeFactor * m_penalty_degree * m_penalty_base;
            if(scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

        }

        private enum BoundaryConditionType {
            passive = 0,
            active = 1
        }

        static public bool write = true;

        // caching/acceleration 
        Vector IdentifyParticleCache_X;
        int IdentifyParticleCache_jCell;
        Particle IdentifyParticleCache_P;


        /// <summary>
        /// For a given position <paramref name="X"/>, 
        /// on the surface of some particle, this method returns the respective particle.
        /// </summary>
        /// <remarks>
        /// Note:
        /// - this is a brute-force approach which may be costly if a large number of particles is used
        /// - some primitive caching is used to accelerate the method
        /// </remarks>
        Particle IdentifyParticle(Vector X, int jCell) {
            if(Particles == null || Particles.Length < 0)
                return null;

            double h = NegLengthScaleS[jCell] * 2;
            double eps = h*h*1e-10;

            if(IdentifyParticleCache_P != null) {
                if((jCell == IdentifyParticleCache_jCell) &&  (X - IdentifyParticleCache_X).AbsSquare() < eps) {
                    return IdentifyParticleCache_P;
                }
            }

            // first, sort the particles according to distance;
            // it is most likely that the particle with closest distance is a "hit"
            // however, I have doubts that this actually accelerates this method
            var PartDist = new (Particle P, double dist)[Particles.Length];
            for(int i = 0; i < Particles.Length; i++) {
                Particle P = Particles[i];
                var Pos = P.Motion.GetPosition();
                double dist = (X - Pos).L2Norm();
                PartDist[i] = (P, dist);
            }
            Array.Sort(PartDist, (T1, T2) => Math.Sign(T1.dist - T2.dist));

            // test if the particle contains the point 
            for(int i = 0; i < Particles.Length; i++) {
                if(PartDist[i].P.Contains(X, h)) {
                    IdentifyParticleCache_jCell = jCell;
                    IdentifyParticleCache_X = X;
                    IdentifyParticleCache_P = PartDist[i].P;
                    return IdentifyParticleCache_P;
                }
            }

            // Fallback; something is weird anyway.
            return PartDist[0].P;
        }

        double GetAngle(Vector X, int jCell, out Particle P) {
            P = IdentifyParticle(X, jCell);
            if(P == null)
                return double.NaN;

            Vector R = X - P.Motion.GetPosition(0);
            double alpha = R.Angle2D();

            return alpha - P.Motion.GetAngle(0);
        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            Vector normalVector = inp.Normal;
            double _penalty = Penalty(inp.jCellIn);
            int dim = normalVector.Dim;

            //uA[0] = u
            //uA[1] = v
            //uA[2] = PhoreticField

            // Grad_u[2,0] = d PhoreticField /dx 
            // Grad_u[2,1] = d PhoreticField /dy
            

            Debug.Assert(uA.Length == ArgumentOrdering.Count);
            Debug.Assert(uB.Length == ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == dim);
            Debug.Assert(Grad_uB.GetLength(1) == dim);
            if(inp.X.IsNullOrEmpty())
                throw new ArgumentNullException("X is null or empty");
            if(inp.X.Abs() < 0)
                throw new ArithmeticException("invalid length of position vector");

            Vector uAFict = new(inp.Parameters_IN[0], inp.Parameters_IN[1]);
            Vector orientationVector = new(inp.Parameters_IN[2], inp.Parameters_IN[3]);
            Vector orientationNormal = new(-orientationVector[1], orientationVector[0]);
            Vector activeStressVector = new(orientationNormal * normalVector > 0 ? -ActiveStress * normalVector[1] : ActiveStress * normalVector[1], orientationNormal * normalVector > 0 ? (ActiveStress * normalVector[0]) : -ActiveStress * normalVector[0]);

            BoundaryConditionType bcType = (orientationVector * normalVector <= 0) ? BoundaryConditionType.passive : BoundaryConditionType.active;
            //BoundaryConditionType bcType = BoundaryConditionType.active;

            
            if(m_UsePhoretic) {

                Vector GradPhoretic = new Vector(Grad_uA[2, 0], Grad_uA[2, 1]);
                Vector GradPhoreticNormal = (inp.Normal * GradPhoretic) * inp.Normal;
                Vector GradPhoreticTangential = GradPhoretic - GradPhoreticNormal;

                /*
                double alpha = GetAngle(inp.X, inp.jCellIn, out var Particle);
                */

                Debug.Assert(ArgumentOrdering.Count == dim + 1);
                double phoreticVal = uA[dim];
                //Vector tangential = normalVector.Rotate2D(Math.PI * 0.5);

                //uAFict = uAFict + GradPhoreticTangential * 0.7;

                // todo: add computation of slip velocity.
                // Note: if the relation is non-linear, special treatment is required!
                if(write) {
                    Console.WriteLine($"todo: add computation of slip velocity; phoretic value is {phoreticVal}, gradient absolute value is {GradPhoretic.Abs()}");
                    write = false; // prevent end-less output; 
                }
            } else {
                Debug.Assert(ArgumentOrdering.Count == dim);
            }


            // Gradient of u and v 
            // =====================
            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for(int d = 0; d < dim; d++) {
                Grad_uA_xN += Grad_uA[Component, d] * inp.Normal[d];
                Grad_vA_xN += Grad_vA[d] * inp.Normal[d];
            }

            double returnValue = 0.0;

            // 3D for IBM_Solver
            // =====================
            if(dim == 3) {
                if (activeStressVector.Abs() != 0)
                    throw new Exception("active stress only defined in 2D");
                returnValue -= Grad_uA_xN * (vA);                                                    // consistency term
                returnValue -= Grad_vA_xN * (uA[Component] - 0);                                     // symmetry term
                returnValue += _penalty * (uA[Component] - 0) * (vA);                                // penalty term
                Debug.Assert(!(double.IsInfinity(returnValue) || double.IsNaN(returnValue)));
                return returnValue * FluidViscosity;
            }

            // 2D
            // =====================
            switch (bcType) {
                case BoundaryConditionType.passive: {

                    // Grad_aU[0,0] = du_dx; Grad_aU[0,1] = du_dy;
                    // Grad_aU[1,0] = dv_dx; Grad_aU[1,1] = dv_dy;

                    for(int d = 0; d < dim; d++) {
                        returnValue -= FluidViscosity * Grad_uA[Component, d] * vA * normalVector[d]; // consistency term: \/u*n*viscosity 
                        returnValue -= FluidViscosity * Grad_vA[d] * (uA[Component] - uAFict[Component]) * normalVector[d]; // symmetry term 
                    }
                    returnValue += FluidViscosity * (uA[Component] - uAFict[Component]) * vA * _penalty; // penalty tern
                    break;
                }

                case BoundaryConditionType.active: {
                    // normal direction, solid wall
                    for(int dN = 0; dN < dim; dN++) {
                        for(int dD = 0; dD < dim; dD++) {
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
                    for(int d1 = 0; d1 < dim; d1++) {
                        for(int d2 = 0; d2 < dim; d2++) {
                            if(d1 == d2) {
                                P[d1, d2] = 1 - normalVector[d1] * normalVector[d2];
                            } else {
                                P[d1, d2] = normalVector[d1] * normalVector[d2];
                            }
                        }
                    }
                    for(int d1 = 0; d1 < dim; d1++) {
                        for(int d2 = 0; d2 < dim; d2++) {
                            returnValue -= P[d1, d2] * activeStressVector[d2] * (P[d1, Component] * vA);
                        }
                    }
                    break;
                }
                default:
                throw new Exception("Unknown boundary condition at particle surface");
            }
            Debug.Assert(!(double.IsInfinity(returnValue) || double.IsNaN(returnValue)));
            return returnValue;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }

        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        public IList<string> ArgumentOrdering {
            get { 
                var ret = VariableNames.VelocityVector(this.m_D); 
                if(this.m_UsePhoretic) {
                    ret = ret.Append(VariableNames.Phoretic).ToArray();
                }
                return ret;
            }
        }

        /// <summary>
        /// Species ID of the solid
        /// </summary>
        public string PositiveSpecies {
            get { return m_SolidSpecies; }
        }

        /// <summary>
        /// Species ID of the fluid; 
        /// </summary>
        public string NegativeSpecies {
            get { return m_FluidSpc; }
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new[] { VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.Velocity_d(0)),
                               VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.Velocity_d(1)),
                               VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.OrientationVectorComponent(0)),
                               VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.OrientationVectorComponent(1))
                };
            }
        }

        public bool IsMultithreadSafe => false;

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public IEquationComponent CloneForThread() {
            return null;
        }

        public object GetPadlock() {
            return this.Particles[0];
        }
    }
}